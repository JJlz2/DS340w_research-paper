# Eccentric Phase Duration: Reproducing a Meta-Analysis and Testing Its Predictions
DS340W research project · Li Zhu

The parent paper, Amdi & King (2025), is a meta-analysis of how the duration of the eccentric (lowering) phase of resistance training affects maximal strength, muscle hypertrophy and countermovement jump (CMJ) height. This project reproduces it, then treats it as a prediction model: it checks whether studies the meta-analysis never saw land inside its 90% prediction intervals, and compares its intervals with a Bayesian model and a random forest.

## Running the code on a new computer
You need **R 4.4 or newer** (https://cran.r-project.org) and, ideally, **RStudio** (https://posit.co/download/rstudio-desktop/). An internet connection is needed the first time, to install the R packages.

1. On GitHub, click **Code → Download ZIP**, then unzip it (on Windows: right-click → **Extract All**).
2. Open the unzipped folder and double-click **`DS340w_research-paper.Rproj`**. RStudio opens with the working directory set to the repository folder.
3. In the RStudio console, run:
   ```r
   source("code/novel/run_all.R")
   ```
   This installs any missing R packages, then runs every script in order: first the reproduction of the parent paper, then the novel approach. It ends by printing Tables I–IV of the paper and a list of every file it wrote. It takes about a minute and a half once the packages are installed.

Without RStudio, open a terminal in the repository folder and run `Rscript code/novel/run_all.R`.

## Where the results go
Every output is saved in **`results/`**, with tables and figures kept apart and each split into the reproduction and the novel approach:

```
results/
  tables/
    reproduction/   A_reproduction_check.csv
    novel/          B_external_test_hao2026.csv, C_loso_*.csv, D_update_before_after.csv,
                    E_model_comparison_*.csv
  figures/
    reproduction/   Figure3-5 forest plots, the orchard plot, parent_diagnostic_plots.pdf
    novel/          fig1-fig3
  logs/             a copy of each script's console output, e.g. 04_compare_models.txt
```

After `run_all.R` finishes in RStudio, every table is also in the **Environment** pane under its file name (click one to view it), and the **Files** pane opens at `results/`.

Every table is rebuilt the same each time: the random forest uses a fixed seed (340) and one thread, and everything else is deterministic. A fresh download, run with newly installed packages, rebuilt every table byte for byte. The figures look the same but may differ slightly in font rendering between computers.

## Papers
| Role | Paper | File |
|---|---|---|
| Parent | Amdi, C. H., & King, A. (2025). The effect of eccentric phase duration on maximal strength, muscle hypertrophy and countermovement jump height: A systematic review and meta-analysis. *Journal of Sports Sciences*, 43(20), 2447–2464. doi:10.1080/02640414.2025.2535198 | `parent_paper/Amdi_King_2025_eccentric_phase_duration.pdf` |
| Reference 1 | Hao, M., Qi, H., Zhao, L., & Han, W. (2026). The effect of squat training with different eccentric contraction tempos on lower limb muscle strength. *PLOS One*, 21(7), e0354725. doi:10.1371/journal.pone.0354725 | `reference_paper_1/journal.pone.0354725.pdf` |
| Reference 2 | Boykin, J. R., et al. (2025). Comparative analysis of acute eccentric contraction-induced changes to the skeletal muscle transcriptome in young and aged mice and humans. *AJP-Regulatory, Integrative and Comparative Physiology*, 328(1), R45–R58. doi:10.1152/ajpregu.00224.2024 | `reference_paper_2/boykin-et-al-2024-….pdf` |

## Repository layout
| Path | Contents |
|---|---|
| `code/run_analysis.R` | The authors' OSF analysis script, made runnable with two fixes (see `code/README.md`) |
| `code/generate_forest_plots.R` | Rebuilds the paper's forest plots (Figures 3–5), which the authors' shared code doesn't include |
| `code/generate_orchard_plot.R` | An orchard plot of the strength results by training status, with prediction intervals |
| `code/novel/` | The novel approach, steps 0–5, and `run_all.R`, which runs everything (see below) |
| `data/` | Inputs: the parent paper's dataset, the raw participant data from Hao et al. (2026), train / test / validation splits of the parent dataset, and the parent paper's published Tables 1–3; see `data/README.md` |
| `data/processed/` | Hao et al. (2026) converted to the parent datasheet's format, and the saved PubMed update search |
| `results/` | Every output: tables, figures and logs (see above) |
| `parent_paper/`, `reference_paper_1/`, `reference_paper_2/` | The three papers |
| `s6uqn-osfstorage-archive/` | The parent paper's original OSF materials, unmodified: datasheet, R script, pre-registration and supplementary files |

## Part 1: reproduction of the parent paper
| Script | What it does | Writes |
|---|---|---|
| `code/run_analysis.R` | The authors' own analysis: every meta-analysis in the paper, sensitivity analyses and diagnostics | its console output; in `run_all.R`, `results/logs/run_analysis.txt` and `results/figures/reproduction/parent_diagnostic_plots.pdf` |
| `code/generate_forest_plots.R` | Figures 3–5 | `results/figures/reproduction/Figure*.png` |
| `code/generate_orchard_plot.R` | Strength results by training status, with prediction intervals | `results/figures/reproduction/Orchard_strength_training_status.png` |

The authors' script ships with `rho = 0.4`, their sensitivity-analysis value, so its output and the figures above use 0.4 (see `code/README.md`). Step 3 of the novel approach refits the analyses with their main value, `rho = 0.6`, and checks them against the published numbers: all nine match within 0.01 (`results/tables/reproduction/A_reproduction_check.csv`).

## Part 2: novel approach
**Question:** where the meta-analysis rests on only a few studies, are its 90% prediction intervals right? Studies it never saw should land inside them about 9 times in 10.

**Tests**, all scored on the interval for an *observed* effect (the published prediction interval widened by the new study's own sampling error):
- an external test on Hao et al. (2026), a trial published after the parent paper's search that meets all of its inclusion criteria;
- leave-one-study-out on the original nine studies (29 held-out effect sizes);
- a refit of the meta-analysis with Hao et al. added.

**Models compared on the same test cases:** the parent's multi-level meta-analysis, a Bayesian multi-level meta-analysis with half-normal(0.5) priors on the between-study SDs, and a quantile regression forest (`ranger`). The plain average of the other studies is a baseline. Scores: coverage, mean width, interval score (width plus a penalty for each miss; lower is better) and mean absolute error.

**Main results**
- All nine of the parent's analyses were reproduced within 0.01 of the published values.
- Hao's strength result (g = −0.85) fell outside the parent's 90% prediction interval for trained lifters (0.07 to 0.60) but inside the widened interval. Hao's CMJ result (g = −3.41) fell outside both.
- In leave-one-study-out, the parent's intervals held 97% of results (wider than needed, mostly because CMJ has only three studies), the Bayesian model's 90% with the lowest interval score (3.05 vs 3.75), and the random forest's 72%. No model predicted the effect itself better than the plain average of the other studies.
- With Hao added, the 90% confidence intervals for both of the parent's moderate-certainty conclusions (trained strength and CMJ) include zero.

### The steps (in `code/novel/`)
| Step | File | What it does | Writes |
|---|---|---|---|
| 0 | `00_install_packages.R` | Installs any missing R packages into your personal library | — |
| 1 | `01_pubmed_update_search.py` | Reruns the parent's PubMed search for newer records (optional; Python 3 and internet) | `data/processed/pubmed_update_search_rerun.csv` |
| 2 | `02_process_hao2026.R` | Turns Hao et al.'s raw participant data into effect-size rows | `data/processed/hao2026_rows.csv` |
| 3 | `03_test_predictions.R` | Reproduces the parent's results, then runs the external test, leave-one-study-out and the refit | `results/tables/reproduction/A_*`, `results/tables/novel/B`–`D_*`, `results/figures/novel/` |
| 4 | `04_compare_models.R` | Compares the parent, Bayesian and random forest models (about a minute) | `results/tables/novel/E_model_comparison_*.csv` |
| 5 | `05_paper_tables.R` | Prints Tables I–IV and the other numbers in the paper from the result tables | — |

`run_all.R` runs step 0, Part 1, then steps 2–5, saving each script's console output to `results/logs/`. Each script can also be run on its own from the repository folder, e.g. `Rscript code/novel/03_test_predictions.R`. Shared code (the authors' effect-size equations and model, unchanged) is in `code/novel/meta_helpers.R`.

**R packages:** `readxl`, `dplyr`, `tidyr`, `metafor`, `clubSandwich`, `ggplot2`, `ranger`, and `tidyverse` for the authors' script. Step 0 installs them. Tested with R 4.4.2, metafor 4.8-0, clubSandwich 0.6.2, ranger 0.18.0, dplyr 1.2.1, tidyr 1.3.2 and ggplot2 4.0.3.

**The PubMed search (step 1)** is the one result that can change over time. The search used in the paper was run on 12 November 2026 and is saved as `data/processed/pubmed_update_search.csv` (144 records, one eligible: Hao et al., 2026). A rerun searches the same dates but writes a separate file, then lists any records that differ from the saved search. PubMed adds subject terms to records weeks after they first appear, so a rerun can find a few older records that did not match the search on 12 November.

### Where each number in the paper comes from
| In the paper | File in `results/tables/` |
|---|---|
| Table I: leave-one-study-out comparison of the models | `novel/E_model_comparison_summary.csv` (leave-one-study-out, all outcomes) |
| Table II: the parent's intervals vs. Hao et al. | `novel/B_external_test_hao2026.csv` |
| Table III: leave-one-study-out by outcome | `novel/E_model_comparison_summary.csv` (by outcome) |
| Table IV: before and after adding Hao et al. | `novel/D_update_before_after.csv` |
| Reproduction of the parent's nine analyses | `reproduction/A_reproduction_check.csv` |
| Interval scores on Hao et al. | `novel/E_model_comparison_summary.csv` (external test) |
| Coverage of the published intervals, case by case | `novel/C_loso_predictions.csv`, `novel/C_loso_coverage_summary.csv` |
