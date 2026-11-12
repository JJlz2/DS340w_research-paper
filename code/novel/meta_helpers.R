### ============================================================
### Shared pieces of the novel approach, sourced by steps 3 and 4:
### the authors' effect-size equations and their meta-analysis model,
### copied from code/run_analysis.R with nothing changed.
###
### Requires metafor, clubSandwich and dplyr to be loaded first.
### ============================================================

### rho = 0.6 is the authors' main-analysis value. (Their shared script
### ships with the 0.4 sensitivity value; see code/README.md.)
rho <- 0.6

meta_cols <- c("descriptor", "outcome", "measure", "study_design",
               "training_status", "study_duration", "match_vl_stf", "multi_joint")

### ------------------------------------------------------------
### Effect sizes: the authors' equations (Borenstein et al. numbering)
### ------------------------------------------------------------
es_between <- function(d) {
  d <- d %>% mutate(
    SDpooled  = sqrt(((slow_n - 1) * slow_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow_n + norm_n - 2)),   # 4.19
    SDpooled2 = sqrt(((slow2_n - 1) * slow2_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow2_n + norm_n - 2)),
    d1 = (slow_delta_mean - norm_delta_mean) / SDpooled,                                                          # 4.18
    d2 = (slow2_delta_mean - norm_delta_mean) / SDpooled2,
    j1 = 1 - 3 / (4 * (slow_n + norm_n - 2) - 1),                                                                 # 4.22
    j2 = 1 - 3 / (4 * (slow2_n + norm_n - 2) - 1),
    var_d1 = (slow_n + norm_n) / (slow_n * norm_n) + d1^2 / (2 * (slow_n + norm_n)),                              # 4.20
    var_d2 = (slow2_n + norm_n) / (slow2_n * norm_n) + d2^2 / (2 * (slow2_n + norm_n)))
  ### one row per comparison; only Mike et al. (2017) has a second slow arm
  first  <- d %>% mutate(g = d1 * j1, var_g = j1^2 * var_d1, comparison = "slow v norm")                         # 4.23-4.24
  second <- d %>% filter(!is.na(slow2_n)) %>%
    mutate(g = d2 * j2, var_g = j2^2 * var_d2, comparison = "slow2 v norm")
  bind_rows(first, second) %>% select(all_of(meta_cols), comparison, g, var_g)
}

es_within <- function(d) {
  d %>% mutate(
    Sdiff = sqrt((norm_delta_sd^2 + slow_delta_sd^2) -
                   2 * between_arm_within_subject_corr * norm_delta_sd * slow_delta_sd),                          # 4.15
    avg_corr = (norm_pre_post_cor + slow_pre_post_cor) / 2,
    Swithin = Sdiff / sqrt(2 * (1 - avg_corr)),                                                                   # 4.27
    d1 = (slow_delta_mean - norm_delta_mean) / Swithin,                                                           # 4.26
    var_d1 = (1 / slow_n) + (d1^2 / (2 * slow_n)) * 2 * (1 - avg_corr),                                           # 4.28
    J = 1 - 3 / (4 * (slow_n - 1) - 1),                                                                           # 4.22
    g = d1 * J, var_g = J^2 * var_d1, comparison = "slow v norm") %>%
    select(all_of(meta_cols), comparison, g, var_g)
}

effect_sizes <- function(sheet, which_outcome) {
  d <- sheet %>% filter(outcome == which_outcome)
  bind_rows(es_between(d %>% filter(study_design == "between")),
            es_within(d %>% filter(study_design == "within")))
}

### ------------------------------------------------------------
### The authors' model: multi-level random effects (effects nested in
### studies), sampling covariance from vcalc(), cluster-robust
### (clubSandwich) inference, 90% CI and 90% PI on the robust df.
### The authors switched to Nelder-Mead where the default optimiser
### failed to converge; the same fallback is used here.
### ------------------------------------------------------------
fit_meta <- function(dat) {
  dat <- dat %>% mutate(es_id = row_number())
  V <- vcalc(vi = var_g, cluster = descriptor, obs = es_id, rho = rho, data = dat)
  fit_with <- function(control) {
    rma.mv(yi = g, V = V, random = ~ 1 | descriptor / es_id, data = dat,
           method = "REML", test = "t", control = control)
  }
  fit <- tryCatch(fit_with(list()),
                  error = function(e) fit_with(list(optimizer = "optim",
                                                    optmethod = "Nelder-Mead",
                                                    maxit = 5000)))
  rob   <- robust(fit, cluster = dat$descriptor, clubSandwich = TRUE)
  est   <- as.numeric(rob$beta)
  se    <- as.numeric(rob$se)
  tau2  <- sum(rob$sigma2)
  tcrit <- qt(0.95, df = rob$df)
  data.frame(studies = length(unique(dat$descriptor)), k = nrow(dat),
             est = est, se = se, tau2 = tau2, df = rob$df,
             ci_lo = est - tcrit * se, ci_hi = est + tcrit * se,
             pi_lo = est - tcrit * sqrt(se^2 + tau2),
             pi_hi = est + tcrit * sqrt(se^2 + tau2))
}

### Interval for an OBSERVED effect with sampling variance v:
### the reported PI widened by the new study's own sampling error.
obs_interval <- function(m, v) {
  half <- qt(0.95, df = m$df) * sqrt(m$se^2 + m$tau2 + v)
  c(lo = m$est - half, hi = m$est + half)
}
