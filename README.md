# Phased Primer Designer 📊🧬

Generate **phased PCR primers** with balanced early–cycle nucleotide diversity — an approach that minimizes base-calling noise and phasing/pre-phasing artifacts on Illumina instruments.  

[**▶ Launch the app on Hugging Face Spaces**](https://programmeratlarge-phased-primers.hf.space/) &nbsp;*(public, multi-user URL)*

---

## Why Phased Primers?

| Challenge on Illumina SBS | How phasing helps |
|---------------------------|-------------------|
| Low diversity in cycle 1–12 saturates color channels → poor color-matrix estimation | Pre-pending a few staggered bases lets A/T/C/G appear evenly in the first cycles |
| Color over-representation (e.g. all “A”), especially in amplicon sequencing | Multiple phased primers are pooled so sequencer sees a balanced signal |
| Four-, two- and one-channel chemistries need different base/fluor balancing | App shows per-chemistry **color plots** to verify balance |

---

## App Features

| Feature | Description |
|---------|-------------|
| **Balanced phasing algorithm** | Iteratively chooses the least-represented base for each added position, respecting degenerate IUPAC codes. |
| **Chemistry-aware color plots** | Visualise percent contribution to each color channel for: 1) Four-channel (HiSeq/MiSeq), 2) Two-channel original SBS, 3) Two-channel **XLEAP-SBS**, 4) One-channel (iSeq 100). |
| **Nucleotide composition plot** | Side-by-side bar chart of %A, %T, %C, %G for positions 1–12. |
| **Downloadables** | • `nucleotide_percentages_plot.png`  • `color_percentages_plot.png`  • `primers.csv` (adapter + phasing + gene-specific) |
| **Versioned UI** | Version string in the header (`VERSION` constant). |
| **Multi-user public link** | Hosted on HF Spaces; each browser session is isolated. |

---

## Quick Start

1. Open the app: <https://programmeratlarge-phased-primers.hf.space/>
2. **Set phasing amount** (0 – 21).
3. Paste your **gene-specific primer** and **adapter**.
4. Choose sequencing **chemistry**.
5. Click **Generate Primers**.
6. Review the two plots ➜ download PNG + CSV files.

---

## Repository Structure

```text
.
├─ app.py               # main Gradio application
├─ requirements.txt     # python deps (gradio, pandas, matplotlib …)
└─ README.md            # (this file)
```

## Credits 🙏
Original concept & R/Shiny implementation:
Franziska Bonath – https://github.com/FranBonath

Python/Gradio port & enhancements:
Paul Munn, Genomics Innovation Hub, Cornell University

## Updates

* Version 1.0.0: Initial commit (after conversion to Python / Gradio), and addition of color balancing plot

