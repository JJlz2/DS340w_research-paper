### ============================================================
### Novel approach, step 3 of 4 (Validation, Modeling, Presentation)
###
### The parent meta-analysis (Amdi & King, 2025) reports a 90%
### prediction interval (PI) with every pooled effect: the range where
### a new study's effect is expected to fall. This script treats the
### meta-analysis as a prediction model and scores those predictions
### on data it never saw:
###
###   A. Reproduce the parent's published estimates (rho = 0.6).
###   B. External test: does Hao et al. (2026), a trial published after
###      the parent's search, land inside the parent's 90% PIs?
###   C. Leave-one-study-out (LOSO): refit without each original study
###      and check whether that study lands inside the predicted range.
###   D. Update: refit with Hao et al. added and compare before/after.
###
### Each check is scored two ways:
###   - against the PI exactly as the parent reports it, which predicts
###     a new study's TRUE effect; and
###   - against the interval for an OBSERVED effect, which also adds the
###     new study's own sampling error (its var_g). An observed g is
###     always noisier than a true effect, so this is the fair test.
###
### The effect-size equations and model specification are copied from
### the authors' own script (code/run_analysis.R); nothing is changed.
### They live in code/novel/meta_helpers.R, shared with step 4.
###
### Inputs:  data/raw/AK_ECC_datasheet_Kojic.xlsx   (parent dataset)
###          data/processed/hao2026_rows.csv        (from step 2)
### Outputs: results/tables/reproduction/A_*.csv, results/tables/novel/B-D_*.csv,
###          results/figures/novel/*.png
### Run from the repo root:  Rscript code/novel/03_test_predictions.R
### ============================================================

if (!file.exists("data/processed/hao2026_rows.csv"))
  stop("Run this from the repository folder (open DS340w_research-paper.Rproj in RStudio first), ",
       "after step 2.", call. = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(metafor)
  library(clubSandwich)
  library(dplyr)
  library(ggplot2)
})
source("code/novel/meta_helpers.R")

parent <- read_excel("data/raw/AK_ECC_datasheet_Kojic.xlsx")
hao    <- read.csv("data/processed/hao2026_rows.csv", stringsAsFactors = FALSE)

for (d in c("results/tables/reproduction", "results/tables/novel", "results/figures/novel"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

### ------------------------------------------------------------
### Data
### ------------------------------------------------------------
str_parent <- effect_sizes(parent, "strength")
hyp_parent <- effect_sizes(parent, "hypertrophy")
cmj_parent <- effect_sizes(parent, "power")
hao_es     <- es_between(hao)

### Strength subgroups used by the parent (Table 2)
subgroups <- list(
  "Strength: all studies"     = function(d) d,
  "Strength: trained"         = function(d) filter(d, training_status == "t"),
  "Strength: untrained"       = function(d) filter(d, training_status == "ut"),
  "Strength: multi-joint"     = function(d) filter(d, multi_joint == "mj"),
  "Strength: single-joint"    = function(d) filter(d, multi_joint == "sj"),
  "Strength: sets to failure" = function(d) filter(d, match_vl_stf == "stf"),
  "Strength: volume-load"     = function(d) filter(d, match_vl_stf == "vl")
)

### ------------------------------------------------------------
### A. Reproduce the published estimates (parent Table 2)
### ------------------------------------------------------------
published <- data.frame(
  model = c(names(subgroups), "Hypertrophy: all studies", "CMJ: all studies"),
  pub_est   = c(0.25, 0.33, 0.15, 0.08, 0.40, 0.26, 0.25, 0.05, -0.73),
  pub_ci_lo = c(-0.17, 0.07, -0.80, -0.73, -0.22, -0.73, 0.04, -0.22, -1.35),
  pub_ci_hi = c(0.67, 0.60, 1.10, 0.90, 1.03, 1.25, 0.45, 0.33, -0.12),
  pub_pi_lo = c(-0.86, 0.07, -2.00, -1.69, -0.68, -1.96, 0.04, -0.54, -1.35),
  pub_pi_hi = c(1.36, 0.60, 2.31, 1.86, 1.48, 2.49, 0.45, 0.65, -0.12))

parent_models <- list()
for (nm in names(subgroups)) parent_models[[nm]] <- fit_meta(subgroups[[nm]](str_parent))
parent_models[["Hypertrophy: all studies"]] <- fit_meta(hyp_parent)
parent_models[["CMJ: all studies"]]         <- fit_meta(cmj_parent)

reproduction <- bind_rows(lapply(names(parent_models), function(nm)
  cbind(model = nm, parent_models[[nm]]))) %>%
  left_join(published, by = "model") %>%
  mutate(max_abs_diff = pmax(abs(est - pub_est), abs(ci_lo - pub_ci_lo), abs(ci_hi - pub_ci_hi),
                             abs(pi_lo - pub_pi_lo), abs(pi_hi - pub_pi_hi)))
write.csv(plain_names(reproduction), "results/tables/reproduction/A_reproduction_check.csv", row.names = FALSE)

### ------------------------------------------------------------
### B. External test on Hao et al. (2026)
### Only predictions that apply to Hao's design are tested: Hao's
### participants were trained, did sets to failure, and squatted
### (multi-joint), so the untrained / volume-load / single-joint
### subgroups make no prediction about it.
### ------------------------------------------------------------
roped <- c(strength = 0.25, power = 0.11)   # the parent's practically-trivial bounds
tests <- data.frame(
  prediction = c("Strength: all studies", "Strength: trained",
                 "Strength: sets to failure", "Strength: multi-joint", "CMJ: all studies"),
  outcome = c("strength", "strength", "strength", "strength", "power"))

external <- bind_rows(lapply(seq_len(nrow(tests)), function(i) {
  m   <- parent_models[[tests$prediction[i]]]
  new <- hao_es %>% filter(outcome == tests$outcome[i])
  oi  <- obs_interval(m, new$var_g)
  data.frame(prediction = tests$prediction[i],
             parent_est = m$est, parent_ci_lo = m$ci_lo, parent_ci_hi = m$ci_hi,
             parent_pi_lo = m$pi_lo, parent_pi_hi = m$pi_hi,
             obs_lo = oi[["lo"]], obs_hi = oi[["hi"]],
             hao_g = new$g, hao_se = sqrt(new$var_g),
             hao_ci_lo = new$g - qnorm(0.95) * sqrt(new$var_g),
             hao_ci_hi = new$g + qnorm(0.95) * sqrt(new$var_g),
             inside_reported_pi = new$g >= m$pi_lo & new$g <= m$pi_hi,
             inside_observed_interval = new$g >= oi[["lo"]] & new$g <= oi[["hi"]],
             rope = unname(roped[tests$outcome[i]]))
}))
write.csv(plain_names(external), "results/tables/novel/B_external_test_hao2026.csv", row.names = FALSE)

### ------------------------------------------------------------
### C. Leave-one-study-out on the original nine studies
### ------------------------------------------------------------
loso_outcome <- function(dat, label) {
  bind_rows(lapply(unique(dat$descriptor), function(s) {
    train <- dat %>% filter(descriptor != s)
    test  <- dat %>% filter(descriptor == s)
    m <- tryCatch(fit_meta(train), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    bind_rows(lapply(seq_len(nrow(test)), function(i) {
      oi <- obs_interval(m, test$var_g[i])
      data.frame(outcome = label, held_out = s, measure = test$measure[i],
                 comparison = test$comparison[i], observed_g = test$g[i],
                 predicted = m$est, pi_lo = m$pi_lo, pi_hi = m$pi_hi,
                 obs_lo = oi[["lo"]], obs_hi = oi[["hi"]],
                 train_studies = m$studies, train_df = m$df)
    }))
  }))
}

loso <- bind_rows(loso_outcome(str_parent, "Strength"),
                  loso_outcome(hyp_parent, "Hypertrophy"),
                  loso_outcome(cmj_parent, "CMJ")) %>%
  mutate(error = observed_g - predicted,
         inside_reported_pi = observed_g >= pi_lo & observed_g <= pi_hi,
         inside_observed_interval = observed_g >= obs_lo & observed_g <= obs_hi)
write.csv(plain_names(loso), "results/tables/novel/C_loso_predictions.csv", row.names = FALSE)

coverage <- loso %>%
  group_by(outcome) %>%
  summarise(test_cases = n(), held_out_studies = n_distinct(held_out),
            coverage_reported_pi = mean(inside_reported_pi),
            coverage_observed_interval = mean(inside_observed_interval),
            mean_abs_error = mean(abs(error)), .groups = "drop") %>%
  bind_rows(loso %>% summarise(outcome = "All outcomes", test_cases = n(),
                               held_out_studies = n_distinct(paste(outcome, held_out)),
                               coverage_reported_pi = mean(inside_reported_pi),
                               coverage_observed_interval = mean(inside_observed_interval),
                               mean_abs_error = mean(abs(error))))
write.csv(plain_names(coverage), "results/tables/novel/C_loso_coverage_summary.csv", row.names = FALSE)

### ------------------------------------------------------------
### D. Update: add Hao et al. (2026) and refit
### ------------------------------------------------------------
str_updated <- bind_rows(str_parent, hao_es %>% filter(outcome == "strength"))
cmj_updated <- bind_rows(cmj_parent, hao_es %>% filter(outcome == "power"))

update_models <- c("Strength: all studies", "Strength: trained",
                   "Strength: sets to failure", "Strength: multi-joint")
updated <- bind_rows(
  bind_rows(lapply(update_models, function(nm)
    cbind(model = nm, version = "Parent (9 studies)", parent_models[[nm]]))),
  bind_rows(lapply(update_models, function(nm)
    cbind(model = nm, version = "Updated (+ Hao 2026)", fit_meta(subgroups[[nm]](str_updated))))),
  cbind(model = "CMJ: all studies", version = "Parent (9 studies)", parent_models[["CMJ: all studies"]]),
  cbind(model = "CMJ: all studies", version = "Updated (+ Hao 2026)", fit_meta(cmj_updated)))
write.csv(plain_names(updated), "results/tables/novel/D_update_before_after.csv", row.names = FALSE)

### ------------------------------------------------------------
### Figures (shared look: deck colours, recessive grid, one g axis)
### ------------------------------------------------------------
col_parent <- "#0077A8"   # parent paper's predictions
col_new    <- "#C4562A"   # new evidence (Hao et al., 2026)
ink        <- "#171B23"
ink_soft   <- "#4B5563"
rope_fill  <- "#ECEFEA"

theme_deck <- function() {
  theme_minimal(base_size = 15, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.4),
          axis.text = element_text(colour = ink_soft),
          axis.title = element_text(colour = ink_soft, size = 13),
          strip.text = element_text(colour = ink, face = "bold", hjust = 0, size = 14),
          legend.position = "top", legend.justification = "left",
          legend.text = element_text(colour = ink_soft, size = 13),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(10, 16, 10, 10))
}

### Figure 1: the external test
f1_rows <- c("Strength: all studies", "Strength: trained",
             "Strength: sets to failure", "CMJ: all studies")
f1 <- external %>%
  filter(prediction %in% f1_rows) %>%
  mutate(label = factor(prediction, levels = rev(f1_rows)), y = as.numeric(label))
lab_pi   <- "Parent paper: 90% prediction interval (dot = pooled g)"
lab_obs  <- "Same, widened by the new study's own sampling error"
lab_new  <- "Hao et al. (2026): observed g with 90% CI"
lab_rope <- "Parent's 'practically trivial' range"
f1_points <- bind_rows(
  f1 %>% transmute(y = y + 0.14, x = parent_est, series = lab_pi),
  f1 %>% transmute(y = y - 0.14, x = hao_g, series = lab_new))
p1 <- ggplot() +
  geom_rect(data = f1, aes(xmin = -rope, xmax = rope, ymin = y - 0.42, ymax = y + 0.42,
                           fill = lab_rope), colour = NA) +
  geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = f1, aes(x = obs_lo, xend = obs_hi, y = y + 0.14, yend = y + 0.14,
                              colour = lab_obs), linewidth = 4.2, alpha = 0.22, lineend = "round") +
  geom_segment(data = f1, aes(x = parent_pi_lo, xend = parent_pi_hi, y = y + 0.14, yend = y + 0.14,
                              colour = lab_pi), linewidth = 1.1, lineend = "round") +
  geom_segment(data = f1, aes(x = hao_ci_lo, xend = hao_ci_hi, y = y - 0.14, yend = y - 0.14,
                              colour = lab_new), linewidth = 1.1, lineend = "round") +
  geom_point(data = f1_points, aes(x = x, y = y, shape = series), size = 5.2, colour = "white",
             show.legend = FALSE) +
  geom_point(data = f1_points, aes(x = x, y = y, colour = series, shape = series), size = 3.6) +
  scale_colour_manual(values = setNames(c(col_parent, col_parent, col_new), c(lab_pi, lab_obs, lab_new)),
                      breaks = c(lab_pi, lab_obs, lab_new),
                      guide = guide_legend(nrow = 3, order = 1,
                                           override.aes = list(shape = c(16, NA, 18),
                                                               linewidth = c(1.1, 4.2, 1.1),
                                                               alpha = c(1, 0.22, 1)))) +
  scale_shape_manual(values = setNames(c(16, 18), c(lab_pi, lab_new)), guide = "none") +
  scale_fill_manual(values = setNames(rope_fill, lab_rope), guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = seq(-4, 3, 1)) +
  scale_y_continuous(breaks = f1$y, labels = as.character(f1$label), expand = expansion(add = 0.5)) +
  labs(x = "Hedges' g   (negative favours the shorter eccentric phase, positive the longer)", y = NULL) +
  theme_deck() +
  theme(legend.box = "vertical", legend.spacing.y = unit(2, "pt"))
ggsave("results/figures/novel/fig1_external_test_hao2026.png", p1, width = 9.5, height = 5.6, dpi = 200, bg = "white")

### Figure 2: LOSO, each held-out effect against its predicted range
f2 <- loso %>%
  mutate(case = paste0(sub(" et al\\.", "", held_out), " - ", gsub("_", " ", measure),
                       ifelse(comparison == "slow2 v norm", " (2nd arm)", "")),
         result = ifelse(inside_observed_interval, "Inside predicted range", "Outside predicted range"),
         outcome = factor(outcome, levels = c("Strength", "Hypertrophy", "CMJ"))) %>%
  group_by(outcome) %>% arrange(observed_g, .by_group = TRUE) %>%
  mutate(case = factor(case, levels = unique(case))) %>% ungroup()
loso_plot <- function(d) {
  ggplot(d, aes(y = case)) +
    geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.4, linetype = "dashed") +
    geom_segment(aes(x = obs_lo, xend = obs_hi, yend = case),
                 colour = col_parent, linewidth = 0.9, alpha = 0.55, lineend = "round") +
    geom_point(aes(x = observed_g, shape = result, colour = result), size = 2.8, stroke = 1.1) +
    scale_shape_manual(values = c("Inside predicted range" = 16, "Outside predicted range" = 4)) +
    scale_colour_manual(values = c("Inside predicted range" = ink, "Outside predicted range" = col_new)) +
    scale_x_continuous(breaks = seq(-5, 3, 1))
}
p2 <- loso_plot(f2) +
  facet_grid(outcome ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Hedges' g of the held-out study (dot) vs. the 90% range predicted from the other studies (bar)",
       y = NULL) +
  theme_deck() +
  theme(axis.text.y = element_text(size = 11),
        strip.text.y = element_text(angle = 0, hjust = 0))
ggsave("results/figures/novel/fig2_loso_coverage.png", p2, width = 10, height = 8, dpi = 200, bg = "white")

### Figure 2b: the strength panel on its own, sized for a slide
p2b <- loso_plot(filter(f2, outcome == "Strength")) +
  labs(x = "Hedges' g: held-out study (dot) vs. 90% range predicted without it (bar)", y = NULL) +
  theme_deck() +
  theme(axis.text.y = element_text(size = 13))
ggsave("results/figures/novel/fig2b_loso_strength.png", p2b, width = 9, height = 5.4, dpi = 200, bg = "white")

### Figure 3: before vs after adding the new study
f3 <- updated %>%
  filter(model %in% c("Strength: all studies", "Strength: trained",
                      "Strength: sets to failure", "CMJ: all studies")) %>%
  mutate(label = factor(model, levels = rev(c("Strength: all studies", "Strength: trained",
                                              "Strength: sets to failure", "CMJ: all studies"))),
         y = as.numeric(label) + ifelse(version == "Parent (9 studies)", 0.13, -0.13))
p3 <- ggplot(f3, aes(colour = version)) +
  geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.4, linetype = "dashed") +
  geom_segment(aes(x = pi_lo, xend = pi_hi, y = y, yend = y), linewidth = 0.8, alpha = 0.6,
               lineend = "round") +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = y, yend = y), linewidth = 2.6, lineend = "round") +
  geom_point(aes(x = est, y = y), size = 3.2, colour = "white") +
  geom_point(aes(x = est, y = y), size = 2.2) +
  scale_colour_manual(values = c("Parent (9 studies)" = col_parent, "Updated (+ Hao 2026)" = col_new)) +
  scale_x_continuous(breaks = seq(-4, 2, 1)) +
  scale_y_continuous(breaks = 1:4, labels = levels(f3$label), expand = expansion(add = 0.5)) +
  labs(x = "Pooled Hedges' g   (thick bar: 90% CI, thin bar: 90% PI)", y = NULL) +
  theme_deck()
ggsave("results/figures/novel/fig3_update_before_after.png", p3, width = 9.5, height = 4.8, dpi = 200, bg = "white")

### ------------------------------------------------------------
### Console summary
### ------------------------------------------------------------
options(width = 160)
cat("\n=== A. Reproduction (rho = 0.6) ===\n")
print(plain_names(reproduction %>% select(model, studies, k, est, ci_lo, ci_hi, pi_lo, pi_hi, tau2, max_abs_diff)),
      digits = 3, row.names = FALSE)
cat("\n=== Hao et al. (2026) effect sizes ===\n")
print(plain_names(hao_es %>% mutate(se = sqrt(var_g)) %>% select(outcome, measure, g, se)),
      digits = 3, row.names = FALSE)
cat("\n=== B. External test ===\n")
print(plain_names(external %>% select(prediction, parent_pi_lo, parent_pi_hi, obs_lo, obs_hi, hao_g,
                                      inside_reported_pi, inside_observed_interval)), digits = 3, row.names = FALSE)
cat("\n=== C. LOSO coverage (target ~0.90) ===\n")
print(plain_names(coverage), digits = 3, row.names = FALSE)
cat("\n=== D. Before vs after adding Hao et al. (2026) ===\n")
print(plain_names(updated %>% select(model, version, studies, k, est, ci_lo, ci_hi, pi_lo, pi_hi, tau2)),
      digits = 3, row.names = FALSE)

### The result tables left in the R session get the same plain-English
### column names as the saved files, for View() in RStudio.
reproduction <- plain_names(reproduction)
external     <- plain_names(external)
loso         <- plain_names(loso)
coverage     <- plain_names(coverage)
updated      <- plain_names(updated)
