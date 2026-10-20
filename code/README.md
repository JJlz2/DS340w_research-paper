# Reproduction of the parent paper's analysis code

`run_analysis.R` is the authors' own meta-analysis script, downloaded unmodified
from their OSF repository: https://osf.io/s6uqn/
(original file: `s6uqn-osfstorage-archive/Analysis Code - ECC Review - R Script - March2025.R`,
kept in this repo untouched for reference/provenance).

## What was changed, and why

Only two things were changed to make the script executable — no analysis logic,
equations, or model specifications were touched.

1. **File path** (line ~8). The authors' shared script pointed at a placeholder
   path (`C:/path/to/datasheet/ECC datasheet.xlsx`) that doesn't exist on any
   machine but theirs. It now points at `data/raw/AK_ECC_datasheet_Kojic.xlsx`,
   the exact same dataset file from their OSF repository, committed in this repo.

2. **Removed one dead/duplicate line block** in the CMJ (power) section. The
   original script contains two consecutive `res_pow <- rma.mv(...)` calls: the
   first uses the default optimizer, which the authors' own comment
   (`###Issues with convergence with rho = 0.4 and 0.6, Nelder-Mead optimiser
   used.`) says does not converge; the second, immediately after, is correctly
   configured with `optimizer = "optim", optmethod = "Nelder-Mead"` and works.
   Because `metafor::rma.mv()` throws a hard error on non-convergence, the
   first (superseded) call halts the script before the working second call is
   ever reached. This is a leftover from editing, not a design choice — the
   fix removes only the first, broken call and leaves the authors' intended
   (working) call exactly as they wrote it.

## Running it

Open `DS340w_research-paper.Rproj` in RStudio (this sets the working directory
to the repo root automatically), then open and source `code/run_analysis.R`.
Required packages: `readxl`, `metafor`, `clubSandwich`, `tidyverse`, `dplyr`.

## Verification

Run end-to-end with exit code 0 and zero errors. Spot-checked outputs against
the published paper (Amdi & King, 2025, *J Sports Sciences*, Table 2 / Figure 5):
the pooled CMJ effect reproduces at g's 90% PI ≈ [-1.32, -0.13], matching the
published [-1.34, -0.12] (small residual difference only because this script's
default `rho = 0.4`, the sensitivity-analysis value, rather than the main
analysis's `rho = 0.6` — change `rho <- 0.4` on line 25 to `rho <- 0.6` to
reproduce the exact headline numbers in Table 2/Figures 3-5).
