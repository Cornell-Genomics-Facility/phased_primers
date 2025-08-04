#! /usr/bin/env python3

"""
This app generates phased PCR primers with balanced nucleotide complexity—that is, sequences with even representation of A, T, C, and G across positions. 
This helps reduce issues like base-calling noise and phasing problems in sequencing platforms like Illumina.
The app also includes a plot of the nucleotide balance by position of sequencing.

Author: Paul Munn, Genomics Innovation Hub, Cornell University
(based on a SHiny app written by Franziska Bonath: https://github.com/FranBonath/phased_primers_shiny/tree/main)

Version history:
- 07/28/2025: Original version
"""

import gradio as gr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random, tempfile, os, uuid
import string

# ───────── global variables ─────────
VERSION = "1.0"          # bump when you release updates
# ────────────────────────────────────


# Defines IUPAC degenerate bases:
# Used to interpret ambiguous bases during complexity evaluation
def compile_library_of_special_nucs():
    return {
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
        "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"],
        "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "C", "G", "T"]
    }

# This step normalizes special bases into standard nucleotides so we can assess actual base balance
def adjust_special_nucleotides(nuc_table):
    """Convert degenerate bases into fractional A/T/C/G counts."""
    special = compile_library_of_special_nucs()
    if nuc_table.shape[1] > 4:
        # Convert to float to avoid dtype clashes when adding fractions
        std = nuc_table.iloc[0, :4].astype(float).copy()
        for col in nuc_table.columns[4:]:
            if col in special:
                cnt = nuc_table.at[0, col]
                for nuc in special[col]:
                    std[nuc] += cnt / len(special[col])
        return pd.DataFrame([std])
    # If only A/T/C/G in table, still make sure they're floats
    return nuc_table.iloc[:, :4].astype(float)

# This is the core algorithm for creating phased primers with maximum nucleotide diversity at early positions
# Primer complexity is engineered and evaluated through:
# 1) Base balance: Calculated per-position relative abundance of A/T/C/G.
# 2) Degenerate bases: Interpreted to contribute partially to standard nucleotides.
# 3) Selection strategy: At each phase step, the least-represented base is chosen to improve uniformity.
# This ensures that when multiple phased primers are pooled, the sequencer sees a mix of nucleotides at early cycles, 
# reducing issues with signal calibration and increasing sequencing quality.
# 
# Detailed steps:
# 1) Split the primer into a list of characters:
#    → CCTAHGGGRBGCAGCAG → ['C','C','T','A','H',...]
# 2) Create amount_of_phasing + 1 identical copies of this primer.
# 3) Iteratively prepend nucleotides to the copies to maximize base diversity at early sequencing cycles:
#    For each new position (phase_count from amount_of_phasing down to 1):
#    - Take a subset of original primer of length phase_count.
#    - Combine with previously added phasing bases.
#    - Count the frequency of all nucleotides (including special).
#    - Normalize using adjust_special_nucleotides().
#    - Find the least-represented base at this position.
#    - If there's a tie, prefer one that’s in the current primer context.
#    - Add it to the beginning of appropriate phased copies.
# 4) Result: a list of staggered primers (with added nucleotides up front) that differ only by a few bases, 
#    but introduce complexity early in the read.
def phase_primer(primer, phasing):
    if phasing == 0:
        return [list(primer)]
    split = list(primer)
    phased = [split.copy() for _ in range(phasing + 1)]
    added = []
    special = compile_library_of_special_nucs()

    for phase_cnt in range(phasing, 0, -1):
        subset = split[:phase_cnt]
        combined = added + subset
        counts = pd.Series(combined).value_counts()
        all_nucs = ["A", "T", "C", "G"] + list(special.keys())
        data = {n: counts.get(n, 0) for n in all_nucs}

        df = adjust_special_nucleotides(pd.DataFrame([data]))
        min_val = df.min(axis=1).values[0]
        choices = df.columns[df.iloc[0] == min_val].tolist()
        choices_in_subset = [c for c in choices if c in subset]
        chosen = random.choice(choices_in_subset or choices)

        added = [chosen] + added
        for i in range(phasing - phase_cnt + 1, phasing + 1):
            phased[i] = [chosen] + phased[i]
    return phased

# Combines each phased primer with the adapter, and names them 
def design_primers(phased, adapter, raw):
    rows = []
    for i, p in enumerate(phased):
        phase_len = len(p) - len(raw)
        phasing = ''.join(p[:phase_len])
        gene = ''.join(p[phase_len:])
        full = adapter + phasing + gene
        rows.append(
            dict(primer_name=f"primer_{i+1}",
                 adapter=adapter,
                 phasing_bases=phasing,
                 gene_specific_primer=gene,
                 full_sequence=full)
        )
    return pd.DataFrame(rows)

def plot_nuc(primer_list):
    max_len = 12
    mat = np.array([p[:max_len] for p in primer_list if len(p) >= max_len])
    nucs = ["A", "T", "C", "G"]
    pos = mat.shape[1]
    counts = {n: [(col == n).sum() for col in mat.T] for n in nucs}

    fig, ax = plt.subplots(figsize=(10, 5))
    width, x = 0.2, np.arange(pos)
    for i, n in enumerate(nucs):
        vals = np.array(counts[n]) / len(primer_list) * 100
        ax.bar(x + i * width, vals, width, label=n)
    ax.set(
        xlabel="Base position from 5′ end (after phasing)",
        ylabel="Nucleotide frequency (%)",
        title="Nucleotide Composition by Position Across Phased Primers",
        ylim=(0, 100)
    )
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels([str(i+1) for i in range(pos)])
    ax.legend(title="Nucleotide", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    return fig

def get_color_map(chem):
    maps = {
        "Four-channels (HiSeq & MiSeq)" : {"red":["A","T"], "green":["G","T"]},
        "Two-channels (original SBS)"   : {"orange":["A"], "red":["C"], "green":["T"], "none":["G"]},
        "Two-channels (XLEAP-SBS)"      : {"blue":["A"], "cyan":["C"], "green":["T"], "none":["G"]},
        "One-channel (iSeq 100)"        : {"green":["A"], "red":["C"], "orange":["T"], "none":["G"]}
    }
    return maps[chem]

def plot_colors(primer_list, chemistry):
    max_len = 12
    mat = np.array([p[:max_len] for p in primer_list if len(p) >= max_len])
    pos = mat.shape[1]

    cmap = get_color_map(chemistry)

    # guarantee a bucket for "none" even if the chemistry dict doesn't include it
    if "none" not in cmap:
        cmap["none"] = []

    # base→color lookup
    base_to_color = {b: c for c, bases in cmap.items() for b in bases}

    colors = list(cmap.keys())
    counts = {c: [0] * pos for c in colors}

    for i in range(pos):
        for b in mat[:, i]:
            counts[base_to_color.get(b, "none")][i] += 1

    fig, ax = plt.subplots(figsize=(10, 5))
    width, x = 0.8 / len(colors), np.arange(pos)

    for idx, c in enumerate(colors):
        vals = np.array(counts[c]) / len(primer_list) * 100
        ax.bar(
            x + idx * width,
            vals,
            width=width,
            label=c,
            color=c if c != "none" else "gray"
        )

    ax.set(
        xlabel="Base position from 5′ end (after phasing)",
        ylabel="Fluorescent signal representation (%)",
        title=f"Color Channel Composition by Position\n({chemistry})",
        ylim=(0, 100),
    )
    ax.set_xticks(x + width * (len(colors) / 2))
    ax.set_xticklabels([str(i + 1) for i in range(pos)])
    ax.legend(title="Color", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    return fig

# ----------  util: save fig / df to temp files ----------
# This is a utility function to save a figure to a temporary file
# and return the path to the file.
# It uses the tempfile module to get a temporary directory and a unique file name.
# It then saves the figure to the file and closes it.
# Finally, it returns the path to the file.
def _tmp_path(filename: str) -> str:
    tmp_dir = tempfile.mkdtemp()           # unique dir per run
    return os.path.join(tmp_dir, filename)

def save_fig_fixed(fig, filename: str) -> str:
    path = _tmp_path(filename)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path

def save_csv_fixed(df, filename: str) -> str:
    path = _tmp_path(filename)
    df.to_csv(path, index=False)
    return path

# ----------  main gradio function ----------
def run_tool(phasing, primer, adapter, chemistry):
    primers   = phase_primer(primer, phasing)
    df        = design_primers(primers, adapter, primer)

    fig_nuc   = plot_nuc(primers)
    fig_color = plot_colors(primers, chemistry)

    # fixed names for downloads
    nuc_png   = save_fig_fixed(fig_nuc,   "nucleotide_percentages_plot.png")
    clr_png   = save_fig_fixed(fig_color, "color_percentages_plot.png")
    csv_path  = save_csv_fixed(df,        "primers.csv")

    return fig_nuc, fig_color, df, nuc_png, clr_png, csv_path

chemistries = [
    "Four-channels (HiSeq & MiSeq)",
    "Two-channels (original SBS)",
    "Two-channels (XLEAP-SBS)",
    "One-channel (iSeq 100)",
]

with gr.Blocks(css="""
  .input-panel {background:#f3f3f3;padding:16px;border-radius:8px;}
""") as demo:
    
    gr.Markdown(
        f"""
        <h2 style="margin-bottom:0.2em;">
            Phased PCR1 Amplicon Primer Designer
            <span style="font-weight: normal; font-size: 0.8em;">
            (version {VERSION})
            </span>
        </h2>
        <p style="font-size: 0.9em; margin-top: 0.5em;">
            This app generates phased PCR primers with balanced nucleotide complexity—that is, sequences with even representation of A, T, C, and G across positions. 
            This helps reduce issues like base-calling noise and phasing problems in sequencing platforms like Illumina.
            The app also includes a plot of the nucleotide balance by position of sequencing.
        </p>
        """
    )

    # ───── inputs on light‑gray panel ─────
    with gr.Column(elem_classes="input-panel"):
        phasing_slider = gr.Slider(
            minimum=0, maximum=21, value=0, step=1,
            label="Amount of primer phasing (integer)"
        )
        ticks = ''.join(f'<span>{i}</span>' for i in range(22))
        gr.HTML(
            f'<div style="display:flex;justify-content:space-between;'
            f'font-size:11px;padding:0 4px 8px 4px;">{ticks}</div>'
        )
        primer_in  = gr.Textbox("CCTAHGGGRBGCAGCAG",
                                label="Gene‑specific primer (5′→3′)")
        adapter_in = gr.Textbox("ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                                label="Adapter sequence (5′→3′)")
        chem_in    = gr.Dropdown(
            ["Four-channels (HiSeq & MiSeq)",
             "Two-channels (original SBS)",
             "Two-channels (XLEAP-SBS)",
             "One-channel (iSeq 100)"],
            value="Four-channels (HiSeq & MiSeq)",
            label="Sequencing chemistry"
        )
        run_btn = gr.Button("Generate Primers")

    # ───── outputs: plots side‑by‑side, then table & downloads ─────
    with gr.Row():
        nuc_plot  = gr.Plot(label="Nucleotide Composition by Position")
        clr_plot  = gr.Plot(label="Color Channel Composition by Position")
    table_out = gr.Dataframe(label="Generated Primers")
    file_nuc  = gr.File(label="Download nucleotide_percentages_plot.png")
    file_clr  = gr.File(label="Download color_percentages_plot.png")
    file_csv  = gr.File(label="Download primers.csv")

    run_btn.click(
        run_tool,
        inputs=[phasing_slider, primer_in, adapter_in, chem_in],
        outputs=[nuc_plot, clr_plot, table_out, file_nuc, file_clr, file_csv]
    )

if __name__ == "__main__":
    demo.launch()
    # demo.launch(server_name="0.0.0.0", server_port=7860, share=True)

