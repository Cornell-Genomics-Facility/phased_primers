# Phased Primer Designer 📊🧬

Generate **phased PCR primers** with balanced early–cycle nucleotide diversity — an approach that minimizes base-calling noise and phasing/pre-phasing artifacts on Illumina instruments.  

[**▶ Launch server version of the app from Hugging Face Spaces**](https://cornell-genomics-facility-phased-primers.hf.space) &nbsp;*(faster, more interactive, server version)*

[**▶ Launch local version of the app from GitHub**](https://cornell-genomics-facility.github.io/phased_primers/) &nbsp;*(slower, less interactive, runs localy in your browser)*

**After designing phased primers:** confirm primer design by checking each primer for self-complementarity (especially 3′) that can lead to self-dimers, hairpins, and cross-dimerization with other primers involving the new phased nucleotides.

---

## Why Phased Primers?

| Challenge on Illumina SBS | How phasing helps |
|---------------------------|-------------------|
| Low diversity in cycle 1–12 saturates color channels → poor color-matrix estimation | Pre-pending a few staggered bases lets A/T/C/G appear evenly in the first cycles |
| Color over-representation (e.g. all “A”), especially in amplicon sequencing | Multiple phased primers are pooled so sequencer sees a balanced signal |
| Four-, two- and one-channel chemistries need different base/fluor balancing | App shows per-chemistry **color plots** to verify balance |

For an in-depth guide to amplicon libraries using phased primers, see the tech note by [Franzika Bonath](https://github.com/FranBonath): [Tech Note: Increased Complexity of Amplicon Libraries using Phased Primers](https://ngisweden.scilifelab.se/2021/01/tech-note-increase-complexity-of-amplicon-libraries-using-phased-primers/).

---

## App Features

| Feature | Description |
|---------|-------------|
| **Chemistry-aware balanced phasing algorithm** | Iteratively chooses the least-represented base for each added position, respecting degenerate IUPAC codes. |
| **Chemistry-aware color plots** | Visualise percent contribution to each color channel for: 1) Four-channel (HiSeq/MiSeq), 2) Two-channel original SBS, 3) Two-channel **XLEAP-SBS**, 4) One-channel (iSeq 100). |
| **Nucleotide composition plot** | Side-by-side bar chart of %A, %T, %C, %G. |
| **Downloadables** | • `nucleotide_percentages_plot.png`  • `color_percentages_plot.png`  • `phased_primers.csv` (adapter + phasing + gene-specific) |
| **Multi-user public link** | Hosted on HF Spaces; each browser session is isolated. |

---

## Illumina Chemistry

Illumina currently ships **four distinct sequencing-by-synthesis (SBS) chemistries**:

| Chemistry | Instruments | Color logic | Compatibility rule |
|-----------|-------------|-------------|--------------------|
| **Four-channel** | HiSeq, MiSeq | Red laser detects **A/C**, green laser detects **G/T** | Each cycle must contain **≥1 red** **and** **≥1 green** signal |
| **Two-channel (original SBS)** | NovaSeq 6000 (S1–S4 flow cells), NextSeq 550, MiniSeq | **A/C = red**, **A/T = green**, **G = no color** | Each cycle must show **any color signal** (red/green) |
| **Two-channel (XLEAP-SBS)** | NovaSeq X/X-Plus | **A/C = blue**, **T/C = green**, **G = no color** | Require “≥1 green pixel” rule; **indices starting with GG are incompatible** |
| **Two-channel (XLEAP-SBS)** | MiSeq i100 | **A/C = blue**, **T/A = green**, **G = no color** | Require “≥1 green pixel” rule; **indices starting with GG are incompatible** |
| **One-channel** | iSeq 100 | Fluor mixing encoded in image intensity, not color | No color check; require **≥1 A or C or T** per cycle |

> ℹ️  Always consult the latest Illumina documentation when designing index or primer pools for a specific platform.

The app renders a **color-balance plot** for each chemistry so you can verify compatibility before ordering primers.

---

## Requirements

Only needed if you clone the repository to run the app locally

| Component | Tested version |
|-----------|----------------|
| **Python** | 3.9 – 3.12 |
| **gradio** | 5.x |
| **pandas** | 2.x |
| **numpy** | 1.26+ |
| **matplotlib** | 3.8+ |

Install with:

```bash
pip install -r requirements.txt
```

---

## Quick Start

1. Open the app: <https://cornell-genomics-facility-phased-primers.hf.space/>
2. **Set phasing amount** (0 – 21. See "**How many phased primers?**" below).
3. Paste your **gene-specific primer** and **adapter**.
4. Choose sequencing **chemistry**.
5. Click **Generate Primers**.
6. Review the two plots ➜ download PNG + CSV files.

**How many phased primers?** Increase the phasing amount until, in the **Color channel by position** plot, every cycle’s required signal(s) sits **above the dashed threshold** across **all positions**.
*The dashed threshold is a general guideline, not a strict pass/fail — optimal levels can vary with instrument, library complexity, and base composition. The threshold value itself differs by chemistry (and is set accordingly in the app).*

* **XLEAP-SBS:** only the **green** channel must exceed the dashed threshold.
* **Original SBS:** both **red** and **green** channels must exceed the dashed threshold.
* **Four-channel:** both **red** and **green** channels must exceed the dashed threshold.
* **One-channel (iSeq 100):** the combined **A+C+T** signal must exceed the dashed threshold.

---

## Repository Structure

```text
.
├─ app.py               # main Gradio application
├─ requirements.txt     # python deps (gradio, pandas, matplotlib …)
└─ README.md            # (this file)
```

---

## Credits 🙏
Original concept & R/Shiny implementation:
Franziska Bonath – https://github.com/FranBonath

Python/Gradio port & enhancements:
Genomics Innovation Hub, Cornell University

---

## Update history:

* Version 1.0.0: Initial commit (after conversion to Python / Gradio)
* Version 1.1.0: Addition of color balancing plot
* Version 1.2.0: Enable user to enter custom phased primers
* Version 1.3.0: Modifications to phasing algorithm to make it chemistry-aware
* Version 1.3.1: Added rules to prevent 4-in-a-row bases (e.g. CCCC) and 5-in-a-row C/G mix (e.g. no CGCGC or CCCGG)
* Version 1.3.2: Added plot below nucleotide plot to show bases contributing to each bar
* Version 1.3.3: Bug fixes to enforce rules to prevent 4-in-a-row bases and 5-in-a-row C/G mix
* Version 1.3.4: Updated rules to prevent 4-in-a-row bases and 5-in-a-row C/G mix
* Version 1.3.5: Added chemistry-aware threshold to plots
* Version 1.3.6: Added XLEAP-SBS MiSeq i100 chemistry
* Version 1.4.0: Added "Multiple primers" custom-phasing mode (one phasing-bases textbox per primer)
