### ============================================================
### Runs everything from a fresh download, in order, and puts every
### output in one place, results/:
###
###   results/tables/reproduction/   the reproduction check (CSV)
###   results/tables/novel/          the novel approach's tables (CSV)
###   results/figures/reproduction/  the parent paper's figures, rebuilt
###   results/figures/novel/         the novel approach's figures
###   results/logs/                  a copy of each script's console output
###
### Order:
###   step 0   install any missing R packages
###   Part 1, reproduction of the parent paper:
###     code/run_analysis.R            the authors' own analysis script
###     code/generate_forest_plots.R   Figures 3-5 (forest plots)
###     code/generate_orchard_plot.R   an orchard plot of the strength results
###   Part 2, novel approach:
###     step 2   turn Hao et al. (2026)'s raw data into effect-size rows
###     step 3   reproduction check, external test, leave-one-study-out, update
###     step 4   compare the parent, Bayesian and random forest models (~1 min)
###     step 5   print the numbers reported in the paper, to check against it
### Step 1, the PubMed search, is Python and needs internet: run it on its own.
###
### In RStudio: open DS340w_research-paper.Rproj, then run
###   source("code/novel/run_all.R")
### Afterwards every results table is in the Environment pane (click one to
### view it) and the Files pane shows results/.
### From a terminal in the repo folder:  Rscript code/novel/run_all.R
### ============================================================

if (!file.exists("data/raw/AK_ECC_datasheet_Kojic.xlsx"))
  stop("Run this from the repository folder: open DS340w_research-paper.Rproj in RStudio first.",
       call. = FALSE)

folders <- file.path("results", c("tables/reproduction", "tables/novel",
                                  "figures/reproduction", "figures/novel", "logs"))
for (d in folders) dir.create(d, showWarnings = FALSE, recursive = TRUE)

### Runs one script in its own environment. Its console output is shown
### (unless show = FALSE) and copied to results/logs/<script>.txt; warnings
### and messages are written out as they happen, so they land in the log too.
run_step <- function(file, env = new.env(), show = TRUE) {
  log <- file.path("results", "logs", sub("\\.R$", ".txt", basename(file)))
  cat("\n==================== ", file, " ====================\n", sep = "")
  sink(log, split = show)
  on.exit(sink(), add = TRUE)
  withCallingHandlers(source(file, local = env),
    warning = function(w) { cat("Warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
    message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") })
  invisible(env)
}

run_step("code/novel/00_install_packages.R")

### ---- Part 1: reproduction of the parent paper ----
### The authors' script prints every model in full and draws Cook's distance,
### leverage and funnel plots; both go to files instead of the screen. The two
### figure scripts then reuse its fitted models from the same environment.
parent <- new.env()
pdf("results/figures/reproduction/parent_diagnostic_plots.pdf", width = 8, height = 6)
tryCatch(run_step("code/run_analysis.R", env = parent, show = FALSE), finally = dev.off())
cat("The authors' full output is in results/logs/run_analysis.txt, and their diagnostic plots in\n",
    "results/figures/reproduction/parent_diagnostic_plots.pdf\n", sep = "")
run_step("code/generate_forest_plots.R", env = parent)
run_step("code/generate_orchard_plot.R", env = parent)

### ---- Part 2: novel approach ----
for (step in c("02_process_hao2026.R", "03_test_predictions.R",
               "04_compare_models.R", "05_paper_tables.R"))
  run_step(file.path("code", "novel", step))

### ---- Where everything is ----
cat("\n==================== Done ====================\nAll results are in results/:\n")
for (d in folders) {
  cat("  ", sub("^results/", "", d), "/\n", sep = "")
  for (f in list.files(d)) cat("      ", f, "\n", sep = "")
}

### Load every table into the session, so in RStudio each one appears in the
### Environment pane under its file name; click it to view
for (f in list.files("results/tables", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE))
  assign(sub("\\.csv$", "", basename(f)), read.csv(f, check.names = FALSE), envir = globalenv())
if (interactive()) {
  cat("\nEvery table is also in the Environment pane (top right): click one to view it.\n")
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    try(rstudioapi::filesPaneNavigate(file.path(getwd(), "results")), silent = TRUE)
}
