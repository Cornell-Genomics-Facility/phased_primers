#! /usr/bin/env python
"""
Phased-primer designer with optional user-defined phased primer input and algorithmic color balancing

Usage example:
    python phased_primers_with_algorithmic_color_balancing_v1.py 

Author: Paul Munn, Genomics Innovation Hub, Cornell University

Version history:
    - 07/28/2025: Original version (1.0)
    - 08/04/2025: Modified to add color balancing plot (1.1)
    - 08/10/2025: Modified to enable user entered phased primers (1.2)
    - 08/18/2025: Modified so that algorithm to choose nucleotieds favors color balancing over nucleotide diversity (1.3)
    - 08/20/2025: Added rules to prevent 4-in-a-row bases (e.g. CCCC) and 5-in-a-row C/G mix (e.g. no CGCGC or CCCGG) (1.3.1)
    - 08/22/2025: Added plot below nucleotide plot to show bases contributing to each bar (1.3.2)
    - 08/25/2025: Bug fixes to enforce rules to prevent 4-in-a-row bases and 5-in-a-row C/G mix (1.3.3)
"""

import gradio as gr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random, tempfile, os
import string

# ───────── global variables ─────────
VERSION = "1.3.3"
# ────────────────────────────────────

# Map degenerate bases
# R = A,G 
# Y = C,T 
# M = A,C 
# K = G,T 
# S = C, G
# W = A,T 
# H = A,C,T 
# B = C,G,T 
# V = A,C,G 
# D = A,G,T 
# N = A,C,G,T 
def compile_library_of_special_nucs():
    return {
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
        "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"],
        "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "C", "G", "T"]
    }

def adjust_special_nucleotides(nuc_table):
    special = compile_library_of_special_nucs()
    if nuc_table.shape[1] > 4:
        std = nuc_table.iloc[0, :4].astype(float).copy()
        for col in nuc_table.columns[4:]:
            if col in special:
                cnt = nuc_table.at[0, col]
                for nuc in special[col]:
                    std[nuc] += cnt / len(special[col])
        return pd.DataFrame([std])
    return nuc_table.iloc[:, :4].astype(float)

def build_custom_phased_list(primer: str, custom_seq: str):
    s = (custom_seq or "").strip().upper().replace(" ", "")
    L = len(s)
    # ["", s[-1:], s[-2:], ..., s[-L:]]
    return [list(s[L-k:] + primer) for k in range(0, L+1)]

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
def phase_primer(primer, phasing, chemistry):
    """
    All chemistries enforce:
      • No 4 identical bases in a row at the front.
      • No 5-long run consisting only of C/G at the front.

    Random base preferences when needed:
      • Two-channels (original SBS): use A/C/T only for random picks (occasionally pick G at random).
      • Two-channels (XLEAP-SBS)   : prefer C/T, allow A less often; occasionally pick G at random.
      • Others (Four-channel, One-channel): uniform among legal candidates.
    """
    if phasing == 0:
        return [list(primer)]

    def violates_constraints(current_added, cand):
        # No 4 identical bases in a row at the front
        if len(current_added) >= 3 and all(x == cand for x in current_added[:3]):
            return True
        # No 5-long run consisting only of C/G at the front
        if cand in {"C", "G"} and len(current_added) >= 4 and all(x in {"C", "G"} for x in current_added[:4]):
            return True
        return False

    def pick_with_weights(pool, weights_map):
        try:
            w = [weights_map.get(b, 1.0) for b in pool]
            if all(v <= 0 for v in w):
                return random.choice(pool)
            return random.choices(pool, weights=w, k=1)[0]
        except Exception:
            return random.choice(pool)

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
        df = adjust_special_nucleotides(pd.DataFrame([data]))  # columns A/T/C/G (floats)

        # minimal count tie-set among A/T/C/G
        # compute the tie set among bases with the current minimum count
        two_channel = chemistry in {"Two-channels (original SBS)", "Two-channels (XLEAP-SBS)"}
        cols_to_consider = ["A", "C", "T"] if two_channel else ["A", "C", "T", "G"]

        df_min = df[cols_to_consider]
        min_val = df_min.min(axis=1).iat[0]
        choices = df_min.columns[df_min.iloc[0] == min_val].tolist()

        # enforce constraints first
        legal_choices = [b for b in choices if not violates_constraints(added, b)]

        # prefer choices present in the current subset when possible
        if legal_choices:
            base_pool = [c for c in legal_choices if c in subset] or legal_choices
        else:
            all_legal = [b for b in ["A", "C", "T", "G"] if not violates_constraints(added, b)]
            base_pool = [c for c in all_legal if c in subset] or all_legal or choices

        # chemistry-specific random selection rules (on constrained pool)
        # print('base pool: ', base_pool)
        if len(base_pool) == 1:
            chosen = base_pool[0]
        elif chemistry == "Two-channels (original SBS)":
            # occasionally pick G at random
            weights = {"A": 1.0, "C": 1.0, "T": 1.0, "G": 0.5}  # G occasionally selected when random
            chosen = pick_with_weights(base_pool, weights)
        elif chemistry == "Two-channels (XLEAP-SBS)":
            # Prefer C/T; A less often; occasionally pick G at random
            weights = {"C": 1.0, "T": 1.0, "A": 0.5, "G": 0.5}  # G occasionally selected when random
            chosen = pick_with_weights(base_pool, weights)
        else:
            # Four-channel & One-channel: uniform among constrained pool
            chosen = random.choice(base_pool)

        # prepend chosen base into all downstream phased primers
        added = [chosen] + added
        for i in range(phasing - phase_cnt + 1, phasing + 1):
            phased[i] = [chosen] + phased[i]

    return phased

def phase_primer_old(primer, phasing, chemistry):
    """
    When a final random choice is needed among equally good bases:
      1) Four-channels (HiSeq & MiSeq): prefer the base that minimizes |(A+C) - (G+T)|
      2) Two-channels (original SBS):   prefer A/C/T; occasionally pick G.
      3) Two-channels (XLEAP-SBS):      prefer C/T; pick A less often; occasionally G.
        Additionally, forbid:
          • 'C' appearing 4 times in a row in the phasing prefix (no "CCCC")
          • any 5-long run consisting only of C/G (e.g., "CGCGC", "CCCGG", etc.)
      4) One-channel (iSeq 100):        prefer A/C/T
    """
    if phasing == 0:
        return [list(primer)]

    def pick_with_weights(pool, weights_map):
        """Weighted random from pool using weights_map (fallback uniform if missing)."""
        try:
            w = [weights_map.get(b, 1.0) for b in pool]
            # guard against all zeros
            if all(v <= 0 for v in w):
                return random.choice(pool)
            return random.choices(pool, weights=w, k=1)[0]
        except Exception:
            return random.choice(pool)

    def violates_xleap_constraints(current_added, cand):
        """
        We add the candidate to the *front* of `added`, which forms the phasing prefix.
        Only new runs that begin at the front can be created; deeper runs were validated earlier.
        - No 'CCCC' at the front.
        - No 5-long run of only C/G at the front.
        """
        # If adding C would create 'CCCC' at the front
        if cand == "C" and len(current_added) >= 3 and all(x == "C" for x in current_added[:3]):
            return True
        # If adding C or G would create a 5-long run of only C/G at the front
        if cand in {"C", "G"} and len(current_added) >= 4 and all(x in {"C", "G"} for x in current_added[:4]):
            return True
        return False

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

        df = adjust_special_nucleotides(pd.DataFrame([data]))  # columns A/T/C/G (floats)
        min_val = df.min(axis=1).values[0]
        choices = df.columns[df.iloc[0] == min_val].tolist()   # subset of {"A","T","C","G"}

        # Prefer nucleotides that already appear in the current primer subset (original behavior)
        sub_cands = [c for c in choices if c in subset]
        base_pool = sub_cands if sub_cands else choices

        # Chemistry-specific tie-breaking
        if len(base_pool) == 1:
            chosen = base_pool[0]
        elif chemistry == "Two-channels (original SBS)":
            # Prefer A/C/T, occasionally G
            weights = {"A": 1.0, "C": 1.0, "T": 1.0, "G": 0.05}
            chosen = pick_with_weights(base_pool, weights)

        elif chemistry == "Two-channels (XLEAP-SBS)":
            # Prefer C/T strongly; A less often; occasionally G
            weights = {"C": 1.0, "T": 1.0, "A": 0.3, "G": 0.05}
            # Enforce constraints by filtering illegal candidates first
            legal_pool = [b for b in base_pool if not violates_xleap_constraints(added, b)]
            pool = legal_pool if legal_pool else base_pool  # fallback if all illegal
            chosen = pick_with_weights(pool, weights)

        else:
            # Other chemistries: unchanged behavior (uniform random among base_pool)
            chosen = random.choice(base_pool)

        # Record and prepend into phased copies (same as original)
        added = [chosen] + added
        for i in range(phasing - phase_cnt + 1, phasing + 1):
            phased[i] = [chosen] + phased[i]

    return phased

def phase_primer_old_old(primer, phasing, chemistry):
    """
    Same algorithm as before, but when a final random choice is needed among equally good bases:
      1) Four-channels (HiSeq & MiSeq): prefer the base that minimizes |(A+C) - (G+T)|
      2) Two-channels (original SBS):   prefer A/C/T
      3) Two-channels (XLEAP-SBS):      prefer C/T
      4) One-channel (iSeq 100):        prefer A/C/T
    """
    if phasing == 0:
        return [list(primer)]

    def prioritized_choice(cands, subset, counts_series, chem):
        # First keep the existing behavior: prefer candidates present in the current subset
        sub_cands = [c for c in cands if c in subset]
        base_pool = sub_cands if sub_cands else cands

        if chem == "Four-channels (HiSeq & MiSeq)":
            # Choose the candidate that best balances A+C vs G+T if added now
            A = float(counts_series.get("A", 0.0))
            C = float(counts_series.get("C", 0.0))
            G = float(counts_series.get("G", 0.0))
            T = float(counts_series.get("T", 0.0))
            def imbalance(b):
                A2, C2, G2, T2 = A, C, G, T
                if b == "A": A2 += 1.0
                elif b == "C": C2 += 1.0
                elif b == "G": G2 += 1.0
                elif b == "T": T2 += 1.0
                return abs((A2 + C2) - (G2 + T2))
            imbalances = [(b, imbalance(b)) for b in base_pool]
            min_val = min(v for _, v in imbalances)
            best = [b for b, v in imbalances if v == min_val]
            return random.choice(best)

        elif chem == "Two-channels (original SBS)":
            pref = {"A", "C", "T"}  # prefer signal channels; G is no-color
            preferred = [b for b in base_pool if b in pref]
            return random.choice(preferred if preferred else base_pool)

        elif chem == "Two-channels (XLEAP-SBS)":
            pref = {"C", "T"}       # C/T are colored; A is blue, G dark but spec asks C/T priority
            preferred = [b for b in base_pool if b in pref]
            return random.choice(preferred if preferred else base_pool)

        elif chem == "One-channel (iSeq 100)":
            pref = {"A", "C", "T"}  # signal vs G (dark)
            preferred = [b for b in base_pool if b in pref]
            return random.choice(preferred if preferred else base_pool)

        # Fallback (unknown chemistry)
        return random.choice(base_pool)

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
        df = adjust_special_nucleotides(pd.DataFrame([data]))  # → columns A/T/C/G with fractional adds
        # choose among least represented standard bases
        min_val = df.min(axis=1).values[0]
        choices = df.columns[df.iloc[0] == min_val].tolist()   # subset of {"A","T","C","G"}

        # If multiple candidates remain, use chemistry-aware priority instead of pure random
        if len(choices) == 1:
            chosen = choices[0]
        else:
            chosen = prioritized_choice(choices, subset, df.iloc[0], chemistry)

        added = [chosen] + added
        for i in range(phasing - phase_cnt + 1, phasing + 1):
            phased[i] = [chosen] + phased[i]

    return phased

def design_primers(phased, adapter, raw):
    rows = []
    for i, p in enumerate(phased):
        phase_len = max(0, len(p) - len(raw))
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
    # Use the shortest primer length so every position has data across all primers
    if not primer_list:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "No primers to plot.", ha="center", va="center")
        ax.axis("off")
        return fig
    L = min(len(p) for p in primer_list)
    if L == 0:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "Primers have zero length.", ha="center", va="center")
        ax.axis("off")
        return fig

    rows = [[ch.upper() for ch in p[:L]] for p in primer_list]
    mat = np.array(rows)
    total = mat.shape[0]

    # Targets and colors
    nucs   = ["A", "C", "T", "G"]
    colors = {"A": "green", "C": "blue", "T": "red", "G": "orange"}

    # Degenerate-to-standard translation
    special = compile_library_of_special_nucs()

    # Floating counts for fractional contributions from degenerate bases
    counts = {n: np.zeros(L, dtype=float) for n in nucs}
    for j in range(L):
        col = mat[:, j]
        for b in col:
            if b in nucs:
                counts[b][j] += 1.0
            elif b in special:  # split among mapped bases
                mapped = special[b]
                frac = 1.0 / len(mapped)
                for m in mapped:
                    counts[m][j] += frac
            else:
                pass  # unknown character → ignore

    # ── Figure with a second subplot showing raw nucleotides per primer under the bars ──
    # ---- Dynamic sizing for the figure and subplots ----
    TOP_IN         = 4.2               # inches reserved for the bar plot
    PER_ROW_IN     = 0.3               # inches per primer row
    BOTTOM_MIN_IN  = 1.2               # clamp min height of the bottom panel
    BOTTOM_MAX_IN  = 8.0               # clamp max height of the bottom panel

    bottom_in = float(np.clip(total * PER_ROW_IN, BOTTOM_MIN_IN, BOTTOM_MAX_IN))
    fig = plt.figure(figsize=(10, TOP_IN + bottom_in))
    gs = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[TOP_IN, bottom_in],
        hspace=0.32  # ↑ more space between top and bottom subplots (was ~0.25–0.30)
    )
    ax  = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)

    # ---- (your existing bar plotting code for ax) ----
    x = np.arange(L)
    nucs   = ["A", "C", "T", "G"]
    colors = {"A": "green", "C": "blue", "T": "red", "G": "orange"}
    group_width = 0.55
    width = group_width / len(nucs)
    offsets = (np.arange(len(nucs)) - (len(nucs)-1)/2) * width

    for i, n in enumerate(nucs):
        vals = (counts[n] / total) * 100.0
        ax.bar(x + offsets[i], vals, width=width, label=n, color=colors[n])

    ax.set(
        xlabel="Base position from 5′ end (after phasing)",
        ylabel="Nucleotide frequency (%)",
        title="Nucleotide Composition by Position Across Phased Primers",
        ylim=(0, 100),
        xlim=(-0.5, L - 0.5),
    )
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(L)])
    ax.axhline(25, color="black", linestyle="--", linewidth=1)
    ax.legend(title="Nucleotide", bbox_to_anchor=(1.02, 1), loc="upper left")

    # ---- Bottom subplot: adaptive row spacing + font size ----
    # Make spacing slightly tighter as rows grow, looser when few rows
    row_spacing = float(np.interp(total, [5, 40], [1.80, 1.35]))  # ↑ larger at high totals
    fontsize    = float(np.interp(total, [5, 40], [10.0, 7.8]))   # optional: tiny shrink with many rows

    ax2.set_xlim(-0.5, L - 0.5)
    ax2.set_ylim(row_spacing * (total - 1) + 0.5, -0.5)  # accommodate spacing
    ax2.set_yticks([])
    ax2.set_xticks(x)
    ax2.set_xticklabels([str(i + 1) for i in range(L)])
    for spine in ["top", "right", "left"]:
        ax2.spines[spine].set_visible(False)

    # Draw letters (bold)
    for j in range(L):
        col = mat[:, j]
        for r, base in enumerate(col):
            c = colors.get(base, "#666666")
            y = r * row_spacing
            ax2.text(
                j, y, base,
                ha="center", va="center",
                fontsize=fontsize,
                fontweight="bold",
                fontfamily="DejaVu Sans Mono",
                color=c
            )

    ax2.set_xlabel("Base position (per-primer bases shown below each position)")
    plt.subplots_adjust(left=0.075, right=0.87)
    plt.tight_layout()
    return fig

def plot_nuc_old(primer_list):
    # Use the shortest primer length so every position has data across all primers
    if not primer_list:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "No primers to plot.", ha="center", va="center")
        ax.axis("off")
        return fig
    L = min(len(p) for p in primer_list)
    if L == 0:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "Primers have zero length.", ha="center", va="center")
        ax.axis("off")
        return fig

    rows = [[ch.upper() for ch in p[:L]] for p in primer_list]
    mat = np.array(rows)
    total = mat.shape[0]

    # Targets and colors
    nucs   = ["A", "C", "T", "G"]
    colors = {"A": "green", "C": "blue", "T": "red", "G": "orange"}

    # Degenerate-to-standard translation
    special = compile_library_of_special_nucs()

    # Floating counts for fractional contributions from degenerate bases
    counts = {n: np.zeros(L, dtype=float) for n in nucs}

    for j in range(L):
        col = mat[:, j]
        for b in col:
            if b in nucs:
                counts[b][j] += 1.0
            elif b in special:  # split among mapped bases
                mapped = special[b]
                frac = 1.0 / len(mapped)
                for m in mapped:
                    counts[m][j] += frac
            else:
                pass  # unknown character → ignore

    # ── Wider spacing between groups: shrink group width and center bars at integer x ──
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(L)
    group_width = 0.55                     # fraction of available space per group (shrink to widen gaps)
    width = group_width / len(nucs)        # bar width within a group
    offsets = (np.arange(len(nucs)) - (len(nucs)-1)/2) * width  # center bars around each x

    for i, n in enumerate(nucs):
        vals = (counts[n] / total) * 100.0
        ax.bar(x + offsets[i], vals, width=width, label=n, color=colors[n])

    ax.set(
        xlabel="Base position from 5′ end (after phasing)",
        ylabel="Nucleotide frequency (%)",
        title="Nucleotide Composition by Position Across Phased Primers",
        ylim=(0, 100),
        xlim=(-0.5, L - 0.5)               # keep first/last groups fully visible
    )
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(L)])

    # --- add a black dotted horizontal line at 25% ---
    ax.axhline(25, color="black", linestyle="--", linewidth=1)

    ax.legend(title="Nucleotide", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    return fig

# def get_color_map(chem):
#     maps = {
#         "Four-channels (HiSeq & MiSeq)" : {"red":["A","C"], "green":["G","T"], "black":["H","R","B","S","W","D","N","K","Y","M","V"]},
#         "Two-channels (original SBS)"   : {"orange":["A"], "red":["C"], "green":["T"], "none":["G"]},
#         "Two-channels (XLEAP-SBS)"      : {"blue":["A"], "cyan":["C"], "green":["T"], "none":["G"]},
#         "One-channel (iSeq 100)"        : {"green":["A"], "red":["C"], "orange":["T"], "none":["G"]},
#     }
#     return maps[chem]

def get_color_spec(chem):
    """
    Returns a list of groups: (legend_label, bar_color, bases_set)
    Note: No 'Deg' group. Degenerate bases are split fractionally among A/C/G/T
    and then aggregated into these groups.
    """
    if chem == "Four-channels (HiSeq & MiSeq)":
        # Two bars: A/C (red) and G/T (green)
        return [
            ("A / C", "red",   {"A", "C"}),
            ("G / T", "green", {"G", "T"}),
        ]
    elif chem == "Two-channels (original SBS)":
        return [
            ("A", "orange",   {"A"}),
            ("C", "red",      {"C"}),
            ("T", "green",    {"T"}),
            ("G", "gray",     {"G"}),
        ]
    elif chem == "Two-channels (XLEAP-SBS)":
        return [
            ("A", "blue",     {"A"}),
            ("C", "cyan",     {"C"}),
            ("T", "green",    {"T"}),
            ("G", "gray",     {"G"}),
        ]
    else:  # "One-channel (iSeq 100)"
        return [
            ("A", "green",    {"A"}),
            ("C", "blue",     {"C"}),
            ("T", "red",      {"T"}),
            ("G", "gray",     {"G"}),
        ]

def plot_colors(primer_list, chemistry):
    # Use the shortest length so every position has data across all primers
    if not primer_list:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "No primers to plot.", ha="center", va="center")
        ax.axis("off")
        return fig
    L = min(len(p) for p in primer_list)
    if L == 0:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "Primers have zero length.", ha="center", va="center")
        ax.axis("off")
        return fig

    rows = [[ch.upper() for ch in p[:L]] for p in primer_list]
    mat = np.array(rows)
    total = mat.shape[0]

    # Fractional translation of degenerate bases → A/C/G/T
    special = compile_library_of_special_nucs()
    bases = ["A", "C", "G", "T"]
    base_counts = {b: np.zeros(L, dtype=float) for b in bases}

    for j in range(L):
        col = mat[:, j]
        for b in col:
            if b in base_counts:
                base_counts[b][j] += 1.0
            elif b in special:
                mapped = special[b]
                frac = 1.0 / len(mapped)
                for m in mapped:
                    if m in base_counts:
                        base_counts[m][j] += frac
            else:
                pass  # unknown → ignore

    # Build plotting groups per chemistry
    groups = []
    if chemistry == "Two-channels (XLEAP-SBS)":
        # C contributes to BOTH A and T channels; remove separate C bar
        groups = [
            ("A / C", "blue",  base_counts["A"] + base_counts["C"]),
            ("T / C", "green", base_counts["T"] + base_counts["C"]),
            ("G",     "gray",  base_counts["G"]),
        ]
    elif chemistry == "Two-channels (original SBS)":
        # A contributes to BOTH C and T channels; remove separate A bar
        groups = [
            ("A / C", "red",   base_counts["A"] + base_counts["C"]),
            ("A / T", "green", base_counts["A"] + base_counts["T"]),
            ("G",     "gray",  base_counts["G"]),
        ]
    else:
        # Other chemistries from spec (no duplication)
        spec = get_color_spec(chemistry)  # list of (label, color, set_of_bases)
        for lab, color, bset in spec:
            arr = np.zeros(L, dtype=float)
            for b in bset:
                arr += base_counts[b]
            groups.append((lab, color, arr))

    # Plot with wider spacing between groups
    fig, ax = plt.subplots(figsize=(10, 5))
    G = len(groups)
    x = np.arange(L)
    group_width = 0.55                   # shrink per-position cluster to widen gaps
    width = group_width / G
    offsets = (np.arange(G) - (G - 1) / 2) * width

    # Allow y-axis to exceed 100% if a channel duplicates contributions (XLEAP/original SBS)
    max_pct = 0.0
    for i, (lab, color, arr) in enumerate(groups):
        vals = (arr / total) * 100.0
        max_pct = max(max_pct, float(np.max(vals)) if vals.size else 0.0)
        ax.bar(x + offsets[i], vals, width=width, label=lab, color=color)

    ax.set(
        xlabel="Base position from 5′ end (after phasing)",
        ylabel="Fluorescent signal representation (%)",
        title=f"Color Channel Composition by Position\n({chemistry})",
        xlim=(-0.5, L - 0.5),
        ylim=(0, max(100, max_pct * 1.05) if max_pct > 0 else 100),
    )
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(L)])

    # --- add a black dotted horizontal line at 25% ---
    ax.axhline(25, color="black", linestyle="--", linewidth=1)

    ax.legend(title="Nucleotides", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    return fig

def _tmp_path(filename: str) -> str:
    tmp_dir = tempfile.mkdtemp()
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

# ---------- custom-primer helpers ----------
ALLOWED = set("ACTG")

# def validate_custom_seq(seq: str) -> str:
#    seq = (seq or "").strip().upper().replace(" ", "")
#    if not seq:
#        raise gr.Error("Please enter at least one base (A/C/T/G).")
#    if any(ch not in ALLOWED for ch in seq):
#        raise gr.Error("Invalid characters: only A, C, T, G are allowed (case-insensitive).")
#    return seq

def _status_html(status: str) -> str:
    if not status:
        return ""
    color = "#b00020" if status.startswith("❌") else "#1b7e1b"
    return f"<div style='color:{color};font-size:0.9em;white-space:pre-wrap;'>{status}</div>"

def build_custom_phased_list(primer: str, custom_seq: str):
    s = (custom_seq or "").strip().upper().replace(" ", "")
    L = len(s)
    # ["", s[-1:], s[-2:], ..., s[-L:]]
    return [list(s[L-k:] + primer) for k in range(0, L+1)]

def _cg_violation_html(seq: str) -> str:
    """Return an HTML error if 'CCCC' or any 5-long C/G-only run exists; else empty string."""
    s = (seq or "").upper()
    if ("AAAA" in s) or ("CCCC" in s) or ("TTTT" in s) or ("GGGG" in s):
        return "<div style='color:#b00020'>❌ Bases cannot appear four times in a row (e.g., no “CCCC”).</div>"
    for i in range(len(s) - 4):
        if all(ch in {"C", "G"} for ch in s[i:i+5]):
            return "<div style='color:#b00020'>❌ No 5-base runs composed only of C/G (e.g., “CGCGC”, “CCCGG”).</div>"
    return ""

def _resolve_adapter(adapter_choice: str, custom_value: str):
    """
    Validate adapter input:
    - If 'Other (user entry)' is selected, allow ONLY A/C/T/G (case-insensitive).
    - Enforce no 'CCCC' and no 5-long C/G-only runs.
    - Convert to UPPERCASE on success.
    """
    if adapter_choice.startswith("Other"):
        seq = ''.join((custom_value or "").split()).upper()
        if not seq:
            return None, "<div style='color:#b00020'>❌ Enter a custom adapter sequence.</div>"
        invalid = sorted({ch for ch in seq if ch not in {"A", "C", "T", "G"}})
        if invalid:
            bad = ", ".join(invalid)
            return None, f"<div style='color:#b00020'>❌ Invalid character(s): {bad}. Use only A, C, T, G.</div>"
        v = _cg_violation_html(seq)
        if v:
            return None, v
        return seq, ""
    # preset: take everything before first space; still enforce constraints
    seq = adapter_choice.split(" ")[0].strip().upper()
    v = _cg_violation_html(seq)
    if v:
        return None, v
    return seq, ""

def run_tool_with_adapter(phasing, primer, adapter_choice, custom_adapter, chemistry, custom_mode, custom_seq):
    # Adapter validation
    adapter_seq, adapter_msg = _resolve_adapter(adapter_choice, custom_adapter)
    if adapter_seq is None:
        return None, None, pd.DataFrame(), None, None, None, adapter_msg, "", ""  # adapter_status, primer_status, phasing_status

    # Primer validation (primer can include degenerate bases; only enforce C/G run rules)
    primer_msg = _cg_violation_html(primer)
    if primer_msg:
        return None, None, pd.DataFrame(), None, None, None, "", primer_msg, ""

    # Custom phasing validation (only when On; A/C/T/G only is handled in run_tool, but add C/G run rules here)
    phasing_msg = ""
    if custom_mode == "On":
        s = (custom_seq or "").upper().replace(" ", "")
        # only enforce run rules if the characters are valid A/C/T/G (run_tool will also validate chars/length)
        if set(s) <= {"A","C","T","G"}:
            phasing_msg = _cg_violation_html(s)
            if phasing_msg:
                return None, None, pd.DataFrame(), None, None, None, "", "", phasing_msg

    # Delegate to the existing runner
    out = run_tool(phasing, primer, adapter_seq, chemistry, custom_mode, custom_seq)
    # Append empty status placeholders (adapter_status, primer_status, phasing_status) if all good
    if isinstance(out, tuple):
        return (*out[:-1], "", "", out[-1])  # keep previous phasing/status_inline message in the last slot
    return out

# ---------- main function ----------
def run_tool(phasing, primer, adapter, chemistry, custom_mode, custom_seq):
    phasing = int(phasing)
    status = ""
    primer = primer.upper()

    if custom_mode == "On":
        seq = (custom_seq or "").strip().upper().replace(" ", "")
        invalid = sorted({ch for ch in seq if ch not in {"A","C","T","G"}})
        if invalid:
            status = f"❌ Invalid character(s): {', '.join(invalid)}. Use only A, C, T, G."
            return None, None, pd.DataFrame(), None, None, None, _status_html(status)
        if len(seq) != phasing:
            status = f"❌ Length mismatch: entered {len(seq)} base(s) but slider is {phasing}."
            return None, None, pd.DataFrame(), None, None, None, _status_html(status)

        primers = build_custom_phased_list(primer, seq)
        status = f"✅ Using custom phased bases “{seq}”. Generated {len(primers)} primers."
    else:
        primers = phase_primer(primer, phasing, chemistry)
        status = f"✅ Using algorithmic phasing (phasing = {phasing}). Generated {len(primers)} primers."

    df        = design_primers(primers, adapter, primer)
    fig_nuc   = plot_nuc(primers)
    fig_color = plot_colors(primers, chemistry)

    return (
        fig_nuc,
        fig_color,
        df,
        save_fig_fixed(fig_nuc,   "nucleotide_percentages_plot.png"),
        save_fig_fixed(fig_color, "color_percentages_plot.png"),
        save_csv_fixed(df,        "phased_primers.csv"),
        _status_html(status),
    )

# ---------- UI ----------
chemistries = [
    "Four-channels (HiSeq & MiSeq)",
    "Two-channels (original SBS)",
    "Two-channels (XLEAP-SBS)",
    "One-channel (iSeq 100)",
]

# def _resolve_adapter(adapter_choice: str, custom_value: str):
#     """
#     Validate adapter input.
#     - If 'Other (user entry)' is selected, allow ONLY A/C/T/G (case-insensitive).
#     - Convert to UPPERCASE when passing through.
#     - On invalid input, return (None, <error_html>) so the UI shows an inline error.
#     """
#     if adapter_choice.startswith("Other"):
#         # remove all whitespace and uppercase
#         seq = ''.join((custom_value or "").split()).upper()
#         if not seq:
#             return None, "<div style='color:#b00020'>❌ Enter a custom adapter sequence.</div>"
#         invalid = sorted({ch for ch in seq if ch not in {"A", "C", "T", "G"}})
#         if invalid:
#             bad = ", ".join(invalid)
#             return None, (
#                 f"<div style='color:#b00020'>❌ Invalid character(s) in custom adapter: {bad}. "
#                 f"Use only A, C, T, G (case-insensitive).</div>"
#             )
#         return seq, ""
# 
#     # preset: take everything before the first space as the adapter sequence
#     seq = adapter_choice.split(" ")[0].strip().upper()
#     return seq, ""

# def run_tool_with_adapter(phasing, primer, adapter_choice, custom_adapter, chemistry, custom_mode, custom_seq):
#     adapter_seq, msg = _resolve_adapter(adapter_choice, custom_adapter)
#     if adapter_seq is None:  # show error inline and clear outputs
#         return None, None, pd.DataFrame(), None, None, None, msg
#     # delegate to your existing run_tool (expects the raw adapter string)
#     return run_tool(phasing, primer, adapter_seq, chemistry, custom_mode, custom_seq)

with gr.Blocks(css="""
  /* Inputs panel (light mode) */
  .input-panel { background:#f3f3f3; padding:16px; border-radius:8px; }
  /* Inputs panel (dark mode overrides) */
  html.dark .input-panel, body.dark .input-panel, [data-theme*="dark"] .input-panel {
    background:#2b2b2b !important; color:inherit; border:1px solid rgba(255,255,255,0.08);
  }
  .inline-row { align-items:center; gap:12px; }

  /* Header: title (left) and Cornell link (right) on the SAME TOP ROW, same <h2> size */
  .header-bar { display:flex; align-items:flex-start; justify-content:space-between; gap:12px; margin-bottom:8px; }
  .header-left  { margin:0; }
  .header-right { margin:0; }              /* keep default <h2> size so it matches the left */
  .header-right a { color:inherit; text-decoration:none; }
  .header-right a:hover { text-decoration:underline; }

  /* Ticks that adapt to theme */
  .tick-grid {
    display:grid;
    grid-template-columns:repeat(22,1fr);
    column-gap:35px;
    font-size:11px;
    padding:0 0 8px 0;
    margin-left:20px;
    justify-items:start;
    user-select:none;
    color:inherit;
  }
  .tick-grid > span { justify-self:start; }
  @media (prefers-color-scheme: dark)  { .tick-grid { color: rgba(255,255,255,.9); } }
  @media (prefers-color-scheme: light) { .tick-grid { color: #444; } }

  /* give the radio extra width */
  #custom-mode { min-width: 160px; }
               
  /* Download buttons styling */
  .downloads-row .gr-button { 
    padding: 10px 16px; 
    border-radius: 8px; 
    font-weight: 600; 
  }
  .downloads-row { gap: 12px; }
               
  /* Make adapter dropdown a bit wider, custom adapter narrower */
  #adapter-dd { min-width: 400px; }         /* nudge wider; adjust as needed */
  #adapter-custom { min-width: 370px; }     /* nudge wider; adjust as needed */
""") as demo:

    demo.queue()

    gr.HTML(
        f"""
        <div class="header-bar">
          <h2 class="header-left">
            Phased PCR1 Amplicon Primer Designer
            <span style="font-weight: normal; font-size: 0.8em;">(version {VERSION})</span>
          </h2>
          <h2 class="header-right">
            <a href="https://github.com/Cornell-Genomics-Facility" target="_blank" rel="noopener">
              Cornell University, Genomics Facility
            </a>
          </h2>
        </div>
        <p style="font-size: 0.9em; margin:6px 0 0 0;">
          Generate phased PCR primers with balanced early-cycle nucleotide diversity.
          &nbsp;<strong>Read the docs:</strong>
          <a href="https://github.com/Cornell-Genomics-Facility/phased_primers" target="_blank" rel="noopener">
            github.com/Cornell-Genomics-Facility/phased_primers
          </a>
        </p>
        """
    )

    with gr.Column(elem_classes="input-panel"):
        phasing_slider = gr.Slider(minimum=0, maximum=21, value=0, step=1,
                                   label="Amount of primer phasing (integer)")
        ticks = ''.join(f'<span>{i}</span>' for i in range(22))
        gr.HTML(f"<div class='tick-grid'>{ticks}</div>")

        # ROW 1: adapter inputs (scale=3) + adapter error (scale=1)
        with gr.Row(elem_classes="inline-row"):
            with gr.Column(scale=3):
                adapter_dd = gr.Dropdown(
                    choices=[
                        "ACACTCTTTCCCTACACGACGCTCTTCCGATCT (TruSeq-P5 Forward)",
                        "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (TruSeq-P7 Reverse)",
                        "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG (Nextera-P5 Forward)",
                        "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG (Nextera-P7 Reverse)",
                        "Other (user entry)",
                    ],
                    value="ACACTCTTTCCCTACACGACGCTCTTCCGATCT (TruSeq-P5 Forward)",
                    label="Adapter",
                    elem_id="adapter-dd"
                )
                adapter_custom = gr.Textbox(
                    value="",
                    label="Custom adapter",
                    placeholder="Select 'Other' to enter custom adapter sequence",
                    interactive=False,
                    elem_id="adapter-custom"
                )
            with gr.Column(scale=1):
                adapter_status = gr.HTML("")

        # ROW 2: primer input (scale=3) + primer error (scale=1)
        with gr.Row(elem_classes="inline-row"):
            with gr.Column(scale=3):
                primer_in  = gr.Textbox(
                    "CCTAHGGGRBGCAGCAG",
                    label="Gene-specific primer (5′→3′)",
                    placeholder="Enter custom primer"
                )
            with gr.Column(scale=1):
                primer_status = gr.HTML("")

        # Enable/disable custom adapter; clear status on change
        def _toggle_adapter(choice):
            if choice.startswith("Other"):
                return gr.update(interactive=True), gr.update(value="")
            else:
                return gr.update(value="", interactive=False), gr.update(value="")
        adapter_dd.change(_toggle_adapter, inputs=adapter_dd, outputs=[adapter_custom, adapter_status])

        # Live validate adapter custom entry
        def _validate_adapter_live(val):
            seq = ''.join((val or "").split()).upper()
            # Accept only ACTG; then enforce C/G run rules
            if not seq:
                return ""
            invalid = sorted({ch for ch in seq if ch not in {"A", "C", "T", "G"}})
            if invalid:
                bad = ", ".join(invalid)
                return f"<div style='color:#b00020'>❌ Invalid character(s): {bad}. Use only A, C, T, G.</div>"
            return _cg_violation_html(seq)
        adapter_custom.change(_validate_adapter_live, inputs=adapter_custom, outputs=adapter_status)

        # Live validate primer input 
        def _validate_primer_live(val):
            seq = ''.join((val or "").split()).upper()
            if not seq:
                return ""
            allowed = set("ACGTRYMKSWHBVDN")  # IUPAC DNA codes
            invalid = sorted({ch for ch in seq if ch not in allowed})
            if invalid:
                bad = ", ".join(invalid)
                return (
                    "<div style='color:#b00020'>❌ Invalid character(s): "
                    f"{bad}. Use only IUPAC DNA codes: A C G T R Y M K S W H B V D N.</div>"
                )
            return _cg_violation_html(seq)  # still enforce literal CCCC and 5×C/G runs
        primer_in.change(_validate_primer_live, inputs=primer_in, outputs=primer_status)

        # Row 2: Chemistry + custom phasing + inline phasing status
        with gr.Row(elem_classes="inline-row"):
            chem_in = gr.Dropdown(
                ["Four-channels (HiSeq & MiSeq)",
                 "Two-channels (original SBS)",
                 "Two-channels (XLEAP-SBS)",
                 "One-channel (iSeq 100)"],
                value="Four-channels (HiSeq & MiSeq)",
                label="Sequencing chemistry",
                scale=1
            )
            custom_mode = gr.Radio(
                choices=["Off", "On"], value="Off",
                label="Enter your own phased primer",
                scale=1, elem_id="custom-mode"
            )
            custom_seq = gr.Textbox(
                value="",
                label="Custom phasing bases",
                placeholder="Click 'On' to enter custom phasing bases",
                interactive=False,
                scale=1
            )
            with gr.Column(scale=1):
                status_inline = gr.HTML("")  # phasing status / general status

        def _toggle(mode):
            if mode == "On":
                return gr.update(interactive=True), gr.update(value="")
            else:
                return gr.update(value="", interactive=False), gr.update(value="")
        custom_mode.change(_toggle, inputs=custom_mode, outputs=[custom_seq, status_inline])

        # Optional live validation for custom phasing bases
        def _validate_custom_seq_live(val):
            s = (val or "").upper().replace(" ", "")
            if not s:
                return ""
            invalid = sorted({ch for ch in s if ch not in {"A","C","T","G"}})
            if invalid:
                bad = ", ".join(invalid)
                return f"<div style='color:#b00020'>❌ Invalid character(s): {bad}. Use only A, C, T, G.</div>"
            return _cg_violation_html(s)
        custom_seq.change(_validate_custom_seq_live, inputs=custom_seq, outputs=status_inline)

        run_btn = gr.Button("Generate Phased Primers", variant="primary")

    # Outputs
    nuc_plot  = gr.Plot(label="Nucleotide Composition by Position")
    clr_plot  = gr.Plot(label="Color Channel Composition by Position")
    gr.Markdown("<hr style='opacity:.2;'>")
    table_out = gr.Dataframe(label="Generated Primers")

    with gr.Row(elem_classes="downloads-row"):
        dl_nuc = gr.DownloadButton(label="Download nucleotide_percentages_plot.png", variant="primary")
        dl_clr = gr.DownloadButton(label="Download color_percentages_plot.png", variant="primary")
        dl_csv = gr.DownloadButton(label="Download primers.csv", variant="primary")

    # NOTE: now we output three status fields: adapter_status, primer_status, status_inline (phasing/general)
    run_btn.click(
        run_tool_with_adapter,
        inputs=[phasing_slider, primer_in, adapter_dd, adapter_custom, chem_in, custom_mode, custom_seq],
        outputs=[nuc_plot, clr_plot, table_out, dl_nuc, dl_clr, dl_csv, adapter_status, primer_status, status_inline]
    )

if __name__ == "__main__":
    demo.launch()
