# Eccentric Phase Duration: Reproducing a Meta-Analysis and Testing Its Predictions
DS340W research project · Li Zhu

The parent paper, Amdi & King (2025), is a meta-analysis of how the duration of the eccentric (lowering) phase of resistance training affects maximal strength, muscle hypertrophy and countermovement jump height. This repository holds the papers, datasets, figures and the parent paper's analysis code, made runnable.

## Papers
| Role | Paper | File |
|---|---|---|
| Parent | Amdi, C. H., & King, A. (2025). The effect of eccentric phase duration on maximal strength, muscle hypertrophy and countermovement jump height: A systematic review and meta-analysis. *Journal of Sports Sciences*, 43(20), 2447–2464. doi:10.1080/02640414.2025.2535198 | `The effect of eccentric phase duration ….pdf` |
| Reference 1 | Hao, M., Qi, H., Zhao, L., & Han, W. (2026). The effect of squat training with different eccentric contraction tempos on lower limb muscle strength. *PLOS One*, 21(7), e0354725. doi:10.1371/journal.pone.0354725 | `reference_paper_1/journal.pone.0354725.pdf` |
| Reference 2 | Boykin, J. R., et al. (2025). Comparative analysis of acute eccentric contraction-induced changes to the skeletal muscle transcriptome in young and aged mice and humans. *AJP-Regulatory, Integrative and Comparative Physiology*, 328(1), R45–R58. doi:10.1152/ajpregu.00224.2024 | `reference_paper_2/boykin-et-al-2024-….pdf` |

## Repository layout
| Path | Contents |
|---|---|
| `data/` | The parent paper's dataset, the raw participant data from Hao et al. (2026), and train / test / validation splits of the parent dataset; see `data/README.md` |
| `code/run_analysis.R` | The authors' OSF analysis script, made runnable with two fixes (see `code/README.md`) |
| `code/generate_forest_plots.R` | Rebuilds the paper's forest plots (Figures 3–5), which the authors' shared code doesn't include |
| `code/generate_orchard_plot.R` | An orchard plot of the strength results by training status, with prediction intervals |
| `figures/` | The reproduced Figures 3–5 |
| `figures and data/` | The orchard plot, and the parent paper's Tables 1–3 as CSV |
| `s6uqn-osfstorage-archive/` | The parent paper's original OSF materials, unmodified: datasheet, R script, pre-registration and supplementary files |
| `code/novel/` | The novel approach, in four steps (see below) |
| `data/processed/` | Hao et al. (2026) converted to the parent datasheet's format, and the updated PubMed search |
| `results/` | Tables from the novel approach: reproduction check, external test, leave-one-study-out and model comparison |
| `figures/novel/` | Figures from the novel approach |

## Reproducing the parent paper
1. Open `DS340w_research-paper.Rproj` in RStudio. This sets the working directory to the repository root.
2. Source `code/run_analysis.R`. The R packages it needs are `readxl`, `metafor`, `clubSandwich`, `tidyverse` and `dplyr`.
3. For the figures, source `code/generate_forest_plots.R` and `code/generate_orchard_plot.R`. The orchard plot also needs `ggplot2`.

## Novel approach (in progress)
- **Hypothesis:** where the meta-analysis rests on only 3–4 studies, its 90% prediction intervals are too narrow. Studies it never saw will land outside them more often than 1 in 10, and adding those studies will change its conclusions.
- **Method:** treat the meta-analysis as a prediction model and score its 90% prediction intervals on unseen data:
  - an external test on Hao et al. (2026), a trial published after the parent paper's search that meets all of its inclusion criteria;
  - leave-one-study-out on the original nine studies;
  - a refit with the new study added.
- **Models compared on the same test cases:** the parent's multi-level meta-analysis, a Bayesian multi-level meta-analysis with half-normal priors, and a quantile regression forest (`ranger`). The plain average of the other studies is a baseline.

### Running the novel approach
Open `DS340w_research-paper.Rproj` (or work from the repository root) and run, in order:
1. `python code/novel/01_pubmed_update_search.py`: reruns the parent paper's PubMed search for newer records (optional; needs internet). Writes `data/processed/pubmed_update_search.csv`.
2. `code/novel/02_process_hao2026.R`: turns Hao et al.'s raw data into effect-size rows. Writes `data/processed/hao2026_rows.csv`.
3. `code/novel/03_test_predictions.R`: reproduces the parent's results, then runs the external test, leave-one-study-out and the refit. Writes `results/A`–`D` and `figures/novel/`.
4. `code/novel/04_compare_models.R`: compares the three models. Writes `results/E_model_comparison_*.csv`. Takes a few minutes.

R packages: `readxl`, `metafor`, `clubSandwich`, `dplyr`, `ggplot2`, `ranger`.
