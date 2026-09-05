### ============================================================
### Generates Figures 3-5 (forest plots) from the parent paper.
###
### These plots are NOT part of the authors' shared OSF code - only
### the underlying statistics are (see code/README.md for the two
### changes needed to make that script runnable). This script builds
### the actual figures from the exact same fitted model objects and
### raw effect-size data that run_analysis.R computes, via
### metafor::forest()/addpoly(), matching the published layout:
### individual study rows (with SEG/LEG Mean/SD columns), subgroup
### summary diamonds, an overall summary diamond, and the paper's
### region of practical equivalence (ROPE) as red dashed lines.
###
### Run from the repo root (or open via the .Rproj).
### ============================================================

library(metafor)

### Reuse every object already computed by the reproduction script -
### source it once to populate this session. (This also reproduces
### its own funnel/diagnostic plots to a stray Rplots.pdf, which is
### gitignored and harmless.)
source("code/run_analysis.R")

dir.create("figures", showWarnings = FALSE)

### 90% CI for a clubSandwich-robust rma.mv model, computed exactly
### as run_analysis.R does throughout (t-distribution, robust df).
robust_ci90 <- function(robust_obj) {
  est <- as.numeric(robust_obj$beta)
  se  <- as.numeric(robust_obj$se)
  df  <- robust_obj$df
  tcrit <- qt(0.95, df = df)
  c(estimate = est, lower = est - tcrit * se, upper = est + tcrit * se)
}

### ================================================================
### Figure 3: Maximal strength, subgrouped by training status
### ================================================================

str_dat <- data_combined_str
### The original script's long-format reshape (mike_rows / mike_rows2)
### renames g2/var_g2 to g/var_g for Mike et al.'s second (6 s) arm,
### but never swaps its raw display columns - slow_delta_mean/sd stay
### the 4 s arm's values on both of Mike's rows. Effect sizes are
### unaffected (computed correctly upstream); only these display
### columns and labels need disambiguating here.
str_dat$leg_mean <- ifelse(!is.na(str_dat$Comparison) & str_dat$Comparison == "slow2 v norm",
                           str_dat$slow2_delta_mean, str_dat$slow_delta_mean)
str_dat$leg_sd   <- ifelse(!is.na(str_dat$Comparison) & str_dat$Comparison == "slow2 v norm",
                           str_dat$slow2_delta_sd, str_dat$slow_delta_sd)
str_dat$label <- ifelse(str_dat$descriptor == "Mike et al. (2017)",
                        paste0(str_dat$descriptor, " (", ifelse(str_dat$Comparison == "slow2 v norm", "6 secs", "4 secs"), ") ", str_dat$measure),
                        paste(str_dat$descriptor, str_dat$measure))

ut_dat <- str_dat[str_dat$training_status == "ut", ]
t_dat  <- str_dat[str_dat$training_status == "t", ]
n_ut <- nrow(ut_dat)
n_t  <- nrow(t_dat)

### Row layout, top to bottom: untrained studies -> untrained diamond
### -> gap -> trained studies -> trained diamond -> gap -> overall diamond
rows_ut     <- rev(seq_len(n_ut)) + n_t + 5
row_ut_diam <- min(rows_ut) - 1
rows_t      <- rev(seq_len(n_t)) + 2
row_t_diam  <- min(rows_t) - 1
row_overall <- row_t_diam - 2

all_dat  <- rbind(ut_dat, t_dat)
all_rows <- c(rows_ut, rows_t)

ut_ci      <- robust_ci90(res_ut_robust)
t_ci       <- robust_ci90(res_t_robust)
overall_ci <- robust_ci90(res_pooled_str_robust)

png("figures/Figure3_strength_forest.png", width = 1400, height = 1100, res = 130)
par(mar = c(4, 4, 1, 2))

forest(x = all_dat$g, vi = all_dat$var_g,
       slab = all_dat$label,
       ilab = cbind(sprintf("%.2f", all_dat$norm_delta_mean), sprintf("%.2f", all_dat$norm_delta_sd),
                    sprintf("%.2f", all_dat$leg_mean), sprintf("%.2f", all_dat$leg_sd)),
       ilab.xpos = c(-5.4, -4.7, -3.6, -2.9),
       rows = all_rows,
       ylim = c(row_overall - 1.5, max(rows_ut) + 3.5),
       xlim = c(-9, 5),
       alim = c(-3, 3),
       level = 90,
       refline = 0,
       digits = 2,
       cex = 0.75,
       xlab = "Estimate (Hedges' g) [90% CI]",
       header = "Author(s) and Measure")

addpoly(ut_ci["estimate"], ci.lb = ut_ci["lower"], ci.ub = ut_ci["upper"],
        row = row_ut_diam, mlab = "Subgroup RE Model (Untrained)", cex = 0.75)
addpoly(t_ci["estimate"], ci.lb = t_ci["lower"], ci.ub = t_ci["upper"],
        row = row_t_diam, mlab = "Subgroup RE Model (Trained)", cex = 0.75)
addpoly(overall_ci["estimate"], ci.lb = overall_ci["lower"], ci.ub = overall_ci["upper"],
        row = row_overall, mlab = "RE Model for All Studies", cex = 0.8)

text(-9, max(rows_ut) + 1.7, "Untrained Participants", pos = 4, font = 2, cex = 0.85)
text(-9, max(rows_t) + 1.7,  "Trained Participants",   pos = 4, font = 2, cex = 0.85)
text(c(-5.05, -3.25), max(rows_ut) + 2.6, c("SEG(Δ)", "LEG(Δ)"), font = 2, cex = 0.8)
text(c(-5.4, -4.7, -3.6, -2.9), max(rows_ut) + 2.1, c("Mean", "SD", "Mean", "SD"), cex = 0.7)

abline(v = c(-0.25, 0.25), col = "red", lty = 2)

dev.off()

cat("Figure 3 (strength) written.\n")

### ================================================================
### Figure 4: Muscle hypertrophy, subgrouped by study design
### ================================================================

hyp_dat <- data_combined_hyp
hyp_dat$label <- paste(hyp_dat$descriptor, hyp_dat$measure)

within_dat  <- hyp_dat[hyp_dat$study_design == "within", ]
between_dat <- hyp_dat[hyp_dat$study_design == "between", ]
n_within  <- nrow(within_dat)
n_between <- nrow(between_dat)

rows_within      <- rev(seq_len(n_within)) + n_between + 5
row_within_diam  <- min(rows_within) - 1
rows_between     <- rev(seq_len(n_between)) + 2
row_between_diam <- min(rows_between) - 1
row_overall_hyp  <- row_between_diam - 2

all_hyp      <- rbind(within_dat, between_dat)
all_hyp_rows <- c(rows_within, rows_between)

within_ci  <- robust_ci90(res_within_hyp_robust)
between_ci <- robust_ci90(res_between_hyp_robust)
overall_hyp_ci <- robust_ci90(res_pooled_hyp_robust)

png("figures/Figure4_hypertrophy_forest.png", width = 1400, height = 1300, res = 130)
par(mar = c(4, 4, 1, 2))

forest(x = all_hyp$g, vi = all_hyp$var_g,
       slab = all_hyp$label,
       ilab = cbind(sprintf("%.2f", all_hyp$norm_delta_mean), sprintf("%.2f", all_hyp$norm_delta_sd),
                    sprintf("%.2f", all_hyp$slow_delta_mean), sprintf("%.2f", all_hyp$slow_delta_sd)),
       ilab.xpos = c(-6.2, -5.3, -3.9, -3.0),
       rows = all_hyp_rows,
       ylim = c(row_overall_hyp - 1.5, max(rows_within) + 3.5),
       xlim = c(-10, 5),
       alim = c(-3, 3),
       level = 90,
       refline = 0,
       digits = 2,
       cex = 0.7,
       xlab = "Estimate (Hedges' g) [90% CI]",
       header = "Author(s) and Measure")

addpoly(within_ci["estimate"], ci.lb = within_ci["lower"], ci.ub = within_ci["upper"],
        row = row_within_diam, mlab = "Subgroup RE Model (Within-participants)", cex = 0.7)
addpoly(between_ci["estimate"], ci.lb = between_ci["lower"], ci.ub = between_ci["upper"],
        row = row_between_diam, mlab = "Subgroup RE Model (Between-participants)", cex = 0.7)
addpoly(overall_hyp_ci["estimate"], ci.lb = overall_hyp_ci["lower"], ci.ub = overall_hyp_ci["upper"],
        row = row_overall_hyp, mlab = "RE Model for All Studies", cex = 0.75)

text(-10, max(rows_within) + 1.7, "Within-participant", pos = 4, font = 2, cex = 0.8)
text(-10, max(rows_between) + 1.7, "Between-participant", pos = 4, font = 2, cex = 0.8)
text(c(-5.75, -3.45), max(rows_within) + 2.6, c("SEG(Δ)", "LEG(Δ)"), font = 2, cex = 0.75)
text(c(-6.2, -5.3, -3.9, -3.0), max(rows_within) + 2.1, c("Mean", "SD", "Mean", "SD"), cex = 0.65)

abline(v = c(-0.28, 0.28), col = "red", lty = 2)

dev.off()

cat("Figure 4 (hypertrophy) written.\n")

### ================================================================
### Figure 5: Countermovement jump (CMJ) height - no subgroups
### ================================================================

pow_dat <- df_long_pow
### Same Mike et al. (2017) display-column quirk as strength (see above).
pow_dat$leg_mean <- ifelse(!is.na(pow_dat$Comparison) & pow_dat$Comparison == "slow2 v norm",
                           pow_dat$slow2_delta_mean, pow_dat$slow_delta_mean)
pow_dat$leg_sd   <- ifelse(!is.na(pow_dat$Comparison) & pow_dat$Comparison == "slow2 v norm",
                           pow_dat$slow2_delta_sd, pow_dat$slow_delta_sd)
pow_dat$label <- ifelse(pow_dat$descriptor == "Mike et al. (2017)",
                        paste0(pow_dat$descriptor, " (", ifelse(pow_dat$Comparison == "slow2 v norm", "6 secs", "4 secs"), ") ", pow_dat$measure),
                        paste(pow_dat$descriptor, pow_dat$measure))

n_pow <- nrow(pow_dat)
rows_pow    <- rev(seq_len(n_pow)) + 2
row_overall_pow <- min(rows_pow) - 2

overall_pow_ci <- robust_ci90(res_pow_robust)

png("figures/Figure5_CMJ_forest.png", width = 1300, height = 650, res = 130)
par(mar = c(4, 4, 1, 2))

forest(x = pow_dat$g, vi = pow_dat$var_g,
       slab = pow_dat$label,
       ilab = cbind(sprintf("%.2f", pow_dat$norm_delta_mean), sprintf("%.2f", pow_dat$norm_delta_sd),
                    sprintf("%.2f", pow_dat$leg_mean), sprintf("%.2f", pow_dat$leg_sd)),
       ilab.xpos = c(-4.3, -3.7, -2.7, -2.1),
       rows = rows_pow,
       ylim = c(row_overall_pow - 1.5, max(rows_pow) + 2.5),
       xlim = c(-7, 3),
       alim = c(-2, 2),
       level = 90,
       refline = 0,
       digits = 2,
       cex = 0.8,
       xlab = "Estimate (Hedges' g) [90% CI]",
       header = "Author(s) and Year")

addpoly(overall_pow_ci["estimate"], ci.lb = overall_pow_ci["lower"], ci.ub = overall_pow_ci["upper"],
        row = row_overall_pow, mlab = "RE Model for All Studies", cex = 0.85)

text(c(-4.0, -2.4), max(rows_pow) + 1.7, c("SEG(Δ)", "LEG(Δ)"), font = 2, cex = 0.85)
text(c(-4.3, -3.7, -2.7, -2.1), max(rows_pow) + 1.3, c("Mean", "SD", "Mean", "SD"), cex = 0.75)

abline(v = c(-0.11, 0.11), col = "red", lty = 2)

dev.off()

cat("Figure 5 (CMJ) written.\n")
cat("All three forest plots written to figures/\n")
