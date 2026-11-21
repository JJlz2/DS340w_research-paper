### ============================================================
### Novel approach, step 5 (Presentation)
###
### Prints the numbers reported in the paper, in the paper's own layout
### and rounding, straight from the result files that steps 3 and 4 just
### wrote. Use it to check the code's output against the paper.
###
### Inputs:  the A, B, D and E tables in results/tables/ (from steps 3 and 4),
###          data/processed/pubmed_update_search.csv (from step 1)
### Run from the repo root:  Rscript code/novel/05_paper_tables.R
### ============================================================

if (!file.exists("results/tables/novel/E_model_comparison_summary.csv"))
  stop("Run this from the repository folder, after steps 3 and 4 ",
       "(or source code/novel/run_all.R).", call. = FALSE)

read <- function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
num  <- function(x) sprintf("%.2f", x)
pct  <- function(x) ifelse(is.na(x), "/", sprintf("%.0f%%", 100 * x))
yn   <- function(x) ifelse(x, "Yes", "No")
show <- function(title, d) {
  cat("\n", title, "\n", sep = "")
  print(d, row.names = FALSE, right = FALSE)
}

summ <- read("results/tables/novel/E_model_comparison_summary.csv")
loso <- subset(summ, Test == "Leave-one-study-out (9 original studies)")
ext  <- subset(summ, Test == "External: Hao et al. (2026)")
pick <- function(d, outcome, model) d[d$Outcome == outcome & d$Model == model, ]

### Table I: leave-one-study-out, all 29 held-out results
models <- c("Parent meta-analysis"                    = "Parent model",
            "Bayesian meta-analysis"                  = "Bayesian (prior scale 0.5)",
            "Bayesian, prior scale 1.0 (sensitivity)" = "Bayesian (prior scale 1.0)",
            "Random forest (ML)"                      = "Random forest",
            "Average of other studies"                = "Average of other trials")
t1 <- do.call(rbind, lapply(names(models), function(m) pick(loso, "All outcomes", m)))
show("Table I. Leave-one-study-out comparison (29 test cases)",
     data.frame(Model    = unname(models),
                MAE      = num(t1$`Mean absolute error`),
                Coverage = pct(t1$`Coverage (target 0.90)`),
                Width    = ifelse(is.na(t1$`Mean interval width`), "/", num(t1$`Mean interval width`)),
                IS       = ifelse(is.na(t1$`Mean interval score (lower is better)`), "/",
                                  num(t1$`Mean interval score (lower is better)`))))

### Table II: the parent paper's intervals vs. Hao et al. (2026)
b <- read("results/tables/novel/B_external_test_hao2026.csv")
b <- b[match(c("Strength: all studies", "Strength: trained", "CMJ: all studies"), b$Prediction), ]
show("Table II. Parent paper's 90% intervals vs. Hao et al. (2026)",
     data.frame(Analysis = b$Prediction,
                `Pooled g` = num(b$`Parent pooled g`),
                `Published PI` = paste0(num(b$`Parent 90% PI lower`), ", ", num(b$`Parent 90% PI upper`)),
                `Widened interval` = paste0(num(b$`90% interval for observed effect, lower`), ", ",
                                            num(b$`90% interval for observed effect, upper`)),
                `Hao g` = num(b$`Hao g`),
                `Inside (PI / widened)` = paste(yn(b$`Inside published PI`), "/",
                                                yn(b$`Inside interval for observed effect`)),
                check.names = FALSE))

### Table III: leave-one-study-out by outcome, as "coverage; interval score"
outs <- c("Strength", "Hypertrophy", "CMJ", "All outcomes")
cell <- function(o, m) {
  r <- pick(loso, o, m)
  paste0(pct(r$`Coverage (target 0.90)`), "; ", num(r$`Mean interval score (lower is better)`))
}
cases <- vapply(outs, function(o) pick(loso, o, "Parent meta-analysis")$`Test cases`, numeric(1))
show("Table III. Leave-one-study-out results by outcome (coverage; interval score)",
     data.frame(`Outcome (cases)` = paste0(outs, " (", cases, ")"),
                `Parent model` = vapply(outs, cell, "", m = "Parent meta-analysis"),
                `Bayesian (prior scale 0.5)` = vapply(outs, cell, "", m = "Bayesian meta-analysis"),
                `Random forest` = vapply(outs, cell, "", m = "Random forest (ML)"),
                check.names = FALSE))

### Table IV: the parent's analyses before and after adding Hao et al. (2026)
d <- read("results/tables/novel/D_update_before_after.csv")
rows <- c("Strength: all studies", "Strength: trained", "Strength: sets to failure",
          "Strength: multi-joint", "CMJ: all studies")
version <- function(v) { x <- subset(d, Version == v); x[match(rows, x$Model), ] }
before  <- version("Parent (9 studies)")
after   <- version("Updated (+ Hao 2026)")
g_ci <- function(r) sprintf("%.2f [%.2f, %.2f]", r$`Pooled g`, r$`90% CI lower`, r$`90% CI upper`)
show("Table IV. Analyses before and after adding Hao et al. (2026)",
     data.frame(Analysis = rows,
                `Before: g [90% CI]` = g_ci(before),
                `After: g [90% CI]` = g_ci(after),
                `tau^2 before -> after` = paste(num(before$`Between-trial variance`), "->",
                                                num(after$`Between-trial variance`)),
                check.names = FALSE))

### Numbers quoted in the text
a <- read("results/tables/reproduction/A_reproduction_check.csv")
diffs <- a$`Largest difference from published`
cat("\nIn the text:\n")
cat(sprintf("  Reproduction: %d of %d parent analyses within 0.01 of the published values (largest difference %.4f).\n",
            sum(diffs < 0.01), length(diffs), max(diffs)))
cat(sprintf("  External test on Hao et al., mean interval score: parent %s, Bayesian %s, random forest %s.\n",
            num(pick(ext, "All outcomes", "Parent meta-analysis")$`Mean interval score (lower is better)`),
            num(pick(ext, "All outcomes", "Bayesian meta-analysis")$`Mean interval score (lower is better)`),
            num(pick(ext, "All outcomes", "Random forest (ML)")$`Mean interval score (lower is better)`)))
pubmed <- read("data/processed/pubmed_update_search.csv")
cat(sprintf("  PubMed update search (saved, run 12 Nov 2026): %d records, %d eligible (Hao et al., 2026).\n",
            nrow(pubmed), sum(grepl("^eligible", pubmed$`Screening status`))))
