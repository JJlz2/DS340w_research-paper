### ============================================================
### Novel approach, step 4 of 4 (Modeling, Validation)
###
### Step 3 scored the parent's own model on studies it never saw. This
### step asks whether a different model would have predicted those
### studies better. Every model is scored on exactly the same cases:
###   - leave-one-study-out (LOSO) on the original nine studies, one
###     outcome at a time, as in step 3 (29 held-out effects); and
###   - the external test on Hao et al. (2026), fitted on all nine.
###
### Models
###   1. Parent: the authors' multi-level meta-analysis (as in step 3).
###   2. Bayesian: the same multi-level structure and sampling
###      covariance, but the two between-study SDs get weakly
###      informative half-normal priors instead of REML point
###      estimates, so the variation between studies can't collapse
###      to zero when there are only a few trials.
###   3. Random forest (ML): a quantile regression forest predicting
###      each effect from the study features in the datasheet.
###   Reference: the plain average of the other studies' effects
###   (point prediction only), to show whether any model beats it.
###
### Scores, all for the OBSERVED effect (the fair test from step 3):
###   - mean absolute error of the point prediction;
###   - coverage of the 90% interval (target 0.90) and its mean width;
###   - interval score (Gneiting & Raftery, 2007): the width plus a
###     penalty for each miss, so a model can't win just by being
###     vague. Lower is better.
###
### Like-for-like: both meta-analyses predict Hao et al. from the
### all-studies models used in LOSO. The parent's subgroup predictions
### (e.g. trained only) are in step 3's output B.
###
### Inputs:  as step 3
### Outputs: results/E_model_comparison_*.csv
### Run from the repo root:  Rscript code/novel/04_compare_models.R
### ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(metafor)
  library(clubSandwich)
  library(dplyr)
  library(ranger)
})
source("code/novel/meta_helpers.R")

parent <- read_excel("data/raw/AK_ECC_datasheet_Kojic.xlsx")
hao    <- read.csv("data/processed/hao2026_rows.csv", stringsAsFactors = FALSE)

outcomes  <- c(Strength = "strength", Hypertrophy = "hypertrophy", CMJ = "power")
parent_es <- bind_rows(lapply(outcomes, function(o) effect_sizes(parent, o)))
hao_es    <- es_between(hao)

models <- c(parent   = "Parent meta-analysis",
            bayes    = "Bayesian meta-analysis",
            rf       = "Random forest (ML)",
            average  = "Average of other studies",
            bayes_s1 = "Bayesian, prior scale 1.0 (sensitivity)")

### ------------------------------------------------------------
### Model 2: Bayesian multi-level meta-analysis
### Same structure as fit_meta(): effects nested in studies, V from
### vcalc(). The overall mean has a flat prior and is integrated out
### exactly; the study-level and effect-level SDs get independent
### half-normal(0, prior_scale) priors and are integrated over a fine
### grid. With prior_scale = 0.5, 95% of the prior mass is below about
### 1 SD unit: sizeable variation between studies is allowed, absurd
### values are not. Scale 1.0 is run as a sensitivity check.
### A new study's observed effect is predicted by a mixture of normals,
### one per grid point, weighted by that point's posterior probability.
### ------------------------------------------------------------
fit_bayes <- function(dat, prior_scale = 0.5, steps = 150) {
  dat <- dat %>% mutate(es_id = row_number())
  V  <- as.matrix(vcalc(vi = var_g, cluster = descriptor, obs = es_id, rho = rho, data = dat))
  y  <- dat$g
  k  <- length(y)
  ZZ <- outer(dat$descriptor, dat$descriptor, "==") * 1      # same-study indicator
  grid <- expand.grid(tau_s = seq(0, 6 * prior_scale, length.out = steps + 1),
                      tau_e = seq(0, 6 * prior_scale, length.out = steps + 1))
  post <- t(vapply(seq_len(nrow(grid)), function(i) {
    S   <- V + grid$tau_s[i]^2 * ZZ + diag(grid$tau_e[i]^2, k)
    R   <- chol(S)
    Si1 <- backsolve(R, forwardsolve(t(R), rep(1, k)))       # S^-1 1
    Siy <- backsolve(R, forwardsolve(t(R), y))               # S^-1 y
    a   <- sum(Si1)
    mu  <- sum(Siy) / a                                      # GLS mean
    q   <- sum(y * Siy) - mu^2 * a                           # residual quadratic form
    c(mu = mu, var_mu = 1 / a,
      loglik = -sum(log(diag(R))) - 0.5 * log(a) - 0.5 * q)  # restricted (REML) likelihood
  }, numeric(3)))
  logpost <- post[, "loglik"] + dnorm(grid$tau_s, 0, prior_scale, log = TRUE) +
                                dnorm(grid$tau_e, 0, prior_scale, log = TRUE)
  w    <- exp(logpost - max(logpost))
  keep <- w > 1e-10
  comps <- data.frame(w = w[keep] / sum(w[keep]), mu = post[keep, "mu"],
                      var_mu = post[keep, "var_mu"],
                      tau2 = grid$tau_s[keep]^2 + grid$tau_e[keep]^2)
  ### For checking against REML: the grid point with the highest likelihood
  top <- which.max(post[, "loglik"])
  list(comps = comps, est = sum(comps$w * comps$mu),
       reml_est = post[top, "mu"], reml_tau2 = grid$tau_s[top]^2 + grid$tau_e[top]^2)
}

### 90% interval for an observed effect with sampling variance v
bayes_interval <- function(b, v) {
  sd   <- sqrt(b$comps$var_mu + b$comps$tau2 + v)
  cdf  <- function(x) sum(b$comps$w * pnorm(x, b$comps$mu, sd))
  lims <- b$est + c(-10, 10) * max(sd)
  q    <- function(p) uniroot(function(x) cdf(x) - p, lims, tol = 1e-8)$root
  c(lo = q(0.05), hi = q(0.95))
}

### ------------------------------------------------------------
### Model 3: random forest (quantile regression forest, Meinshausen 2006)
### One forest across all outcomes, with the outcome as a feature, so
### it learns from every row. Features are the study characteristics
### in the parent's datasheet; muscle_group is left out because it is
### blank for 13 of the 27 rows. Rows are weighted by 1 / var_g, as in
### the meta-analysis. Settings were fixed before seeing any results:
### ranger defaults, plus more trees and proper unordered factors.
### Its 5th-95th percentiles are an interval for an observed effect,
### because the forest is trained on observed effects.
### ------------------------------------------------------------
rf_data <- function(d) {
  data.frame(g               = d$g,
             outcome         = factor(d$outcome, levels = c("strength", "hypertrophy", "power")),
             training_status = factor(d$training_status, levels = c("t", "ut")),
             study_duration  = d$study_duration,
             match_vl_stf    = factor(d$match_vl_stf, levels = c("stf", "vl")),
             multi_joint     = factor(d$multi_joint, levels = c("mj", "sj")))
}

fit_rf <- function(train) {
  ### ranger draws the leaf values behind its quantiles with R's own RNG,
  ### not its seed argument, so set both to make the intervals repeatable
  set.seed(340)
  ranger(g ~ ., data = rf_data(train), case.weights = 1 / train$var_g,
         num.trees = 2000, quantreg = TRUE, respect.unordered.factors = "order",
         seed = 340, num.threads = 1)
}

rf_predict <- function(rf, new) {
  nd <- rf_data(new)
  q  <- predict(rf, nd, type = "quantiles", quantiles = c(0.05, 0.95))$predictions
  data.frame(predicted = predict(rf, nd)$predictions, lo = q[, 1], hi = q[, 2])
}

### ------------------------------------------------------------
### Score every model on one set of unseen cases
### ------------------------------------------------------------
predict_all <- function(train, test, rf) {
  m    <- fit_meta(train)
  b    <- fit_bayes(train)
  b1   <- fit_bayes(train, prior_scale = 1)
  rf_p <- rf_predict(rf, test)
  bind_rows(lapply(seq_len(nrow(test)), function(i) {
    v  <- test$var_g[i]
    oi <- obs_interval(m, v)
    bi <- bayes_interval(b, v)
    b1i <- bayes_interval(b1, v)
    data.frame(outcome = test$outcome[i], study = test$descriptor[i],
               measure = test$measure[i], comparison = test$comparison[i],
               observed_g = test$g[i], model = unname(models),
               predicted = c(m$est, b$est, rf_p$predicted[i], mean(train$g), b1$est),
               lo = c(oi[["lo"]], bi[["lo"]], rf_p$lo[i], NA, b1i[["lo"]]),
               hi = c(oi[["hi"]], bi[["hi"]], rf_p$hi[i], NA, b1i[["hi"]]))
  }))
}

### LOSO: the forest is trained once per held-out study, without any of
### that study's rows (all outcomes), so nothing about it leaks in.
studies    <- unique(parent_es$descriptor)
rf_without <- lapply(setNames(nm = studies), function(s) fit_rf(filter(parent_es, descriptor != s)))

loso <- bind_rows(lapply(outcomes, function(o) {
  dat <- filter(parent_es, outcome == o)
  bind_rows(lapply(unique(dat$descriptor), function(s)
    predict_all(filter(dat, descriptor != s), filter(dat, descriptor == s), rf_without[[s]])))
})) %>% mutate(test = "Leave-one-study-out (9 original studies)")

### External test: everything fitted on all nine studies
rf_all   <- fit_rf(parent_es)
external <- bind_rows(lapply(seq_len(nrow(hao_es)), function(i)
  predict_all(filter(parent_es, outcome == hao_es$outcome[i]), hao_es[i, ], rf_all))) %>%
  mutate(test = "External: Hao et al. (2026)")

interval_score <- function(y, lo, hi, alpha = 0.10) {
  (hi - lo) + (2 / alpha) * pmax(lo - y, 0) + (2 / alpha) * pmax(y - hi, 0)
}

predictions <- bind_rows(loso, external) %>%
  mutate(outcome = names(outcomes)[match(outcome, outcomes)],
         model = factor(model, levels = unname(models)),
         error = observed_g - predicted,
         inside = observed_g >= lo & observed_g <= hi,
         width = hi - lo,
         interval_score = interval_score(observed_g, lo, hi))
write.csv(predictions, "results/E_model_comparison_predictions.csv", row.names = FALSE)

score <- function(d) {
  d %>% summarise(cases = n(),
                  mean_abs_error = mean(abs(error)),
                  coverage_90 = mean(inside),
                  mean_width = mean(width),
                  mean_interval_score = mean(interval_score), .groups = "drop")
}
summary_tbl <- bind_rows(
  predictions %>% group_by(test, model) %>% score() %>% mutate(outcome = "All outcomes"),
  predictions %>% filter(test == first(loso$test)) %>% group_by(test, outcome, model) %>% score()) %>%
  select(test, outcome, everything()) %>%
  arrange(desc(test), factor(outcome, levels = c("All outcomes", names(outcomes))), model)
write.csv(summary_tbl, "results/E_model_comparison_summary.csv", row.names = FALSE)

### ------------------------------------------------------------
### Check: with the priors removed, the Bayesian code's likelihood
### should peak where the authors' REML fit lands (to grid resolution)
### ------------------------------------------------------------
reml_check <- bind_rows(lapply(names(outcomes), function(lbl) {
  dat <- filter(parent_es, outcome == outcomes[[lbl]])
  m <- fit_meta(dat)
  b <- fit_bayes(dat)
  data.frame(outcome = lbl, rma_est = m$est, grid_est = b$reml_est,
             rma_tau2 = m$tau2, grid_tau2 = b$reml_tau2)
}))

### ------------------------------------------------------------
### Console summary
### ------------------------------------------------------------
options(width = 160)
cat("\n=== Check: Bayesian likelihood vs the authors' REML fit ===\n")
print(reml_check, digits = 3, row.names = FALSE)
cat("\n=== Model comparison (90% intervals for the observed effect; lower interval score is better) ===\n")
print(as.data.frame(summary_tbl), digits = 3, row.names = FALSE)
cat("\n=== External test, case by case ===\n")
print(predictions %>% filter(test == first(external$test)) %>%
        select(outcome, model, observed_g, predicted, lo, hi, inside),
      digits = 3, row.names = FALSE)
