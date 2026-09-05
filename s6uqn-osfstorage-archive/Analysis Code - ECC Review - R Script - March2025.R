library(readxl)
library(metafor)
library(clubSandwich)
library(tidyverse)
library(dplyr)

###Read in associated datasheet available. 
df <- read_excel("C:/path/to/datasheet/ECC datasheet.xlsx")

###Assign sampling error correlation co-efficient
###rho = 0.6 is used in the main analysis
###Sensitivity analysis used rho = 0.4 and 0.8, which can be assigned here over 0.6 if the user is interested
rho <- 0.4

###Calculate effect sizes separately for within- and between designs.
###Select between-participant designs.
data_between_str <- subset(df, outcome == "strength" & study_design == "between")

###calculate effect sizes and variances
###4.xx values are equation numbers from Borenstein et al.
data_between_str <- data_between_str %>%
  mutate(
    # 4.19
    SDpooled = sqrt(((slow_n - 1) * slow_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow_n + norm_n - 2)),
    SDpooled2 = sqrt(((slow2_n - 1) * slow2_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow2_n + norm_n - 2)),
    # 4.18
    d = (slow_delta_mean - norm_delta_mean) / (SDpooled),
    d2 = (slow2_delta_mean - norm_delta_mean) / (SDpooled2),
    # 4.22
    j = (1 - (3 / (4 * (slow_n + norm_n - 2) - 1))),
    j2 = (1 - (3 / (4 * (slow2_n + norm_n - 2) - 1))),
    # 4.23
    g = (d * j),
    g2 = (d2 * j2),
    # 4.20
    var_d = (slow_n + norm_n) / (slow_n * norm_n) + (d^2 / (2 * (slow_n + norm_n))),
    var_d2 = (slow2_n + norm_n) / (slow2_n * norm_n) + (d2^2 / (2 * (slow2_n + norm_n))),
    # 4.24
    var_g = ((j^2) * (var_d)),
    var_g2 = ((j2^2) * (var_d2)),
    # 4.25
    se_g = sqrt(var_g),
    se_g2 = sqrt(var_g2)
  )

###Mutate the sheet into a longer format, separating out the effect sizes for Mike et al. which has two comparator conditions.
mike_rows <- data_between_str %>%
  filter(descriptor == "Mike et al. (2017)") %>%
  select(study, descriptor, measure, study_design, g, var_g, training_status, study_duration, match_vl_stf, multi_joint, norm_pre_mean, norm_pre_sd, norm_post_mean, norm_post_sd, norm_delta_mean, norm_delta_sd, slow_pre_mean, slow_pre_sd, slow_post_mean, slow_post_sd, slow_delta_mean, slow_delta_sd, slow2_pre_mean, slow2_pre_sd, slow2_post_mean, slow2_post_sd, slow2_delta_mean, slow2_delta_sd) %>%
  mutate(Comparison = "slow v norm")

mike_rows2 <- data_between_str %>%
  filter(descriptor == "Mike et al. (2017)") %>%
  select(study, descriptor, measure, study_design, g2, var_g2, training_status, study_duration, match_vl_stf, multi_joint, norm_pre_mean, norm_pre_sd, norm_post_mean, norm_post_sd, norm_delta_mean, norm_delta_sd, slow_pre_mean, slow_pre_sd, slow_post_mean, slow_post_sd, slow_delta_mean, slow_delta_sd, slow2_pre_mean, slow2_pre_sd, slow2_post_mean, slow2_post_sd, slow2_delta_mean, slow2_delta_sd) %>%
  rename(g = g2, var_g = var_g2) %>%
  mutate(Comparison = "slow2 v norm")

###Combine all data into long format
df_long_between_str <- data_between_str %>%
  filter(descriptor != "Mike et al. (2017)") %>%
  select(study, descriptor, measure, study_design, g, var_g, training_status, study_duration, match_vl_stf, multi_joint, norm_pre_mean, norm_pre_sd, norm_post_mean, norm_post_sd, norm_delta_mean, norm_delta_sd, slow_pre_mean, slow_pre_sd, slow_post_mean, slow_post_sd, slow_delta_mean, slow_delta_sd, slow2_pre_mean, slow2_pre_sd, slow2_post_mean, slow2_post_sd, slow2_delta_mean, slow2_delta_sd) %>%
  mutate(Comparison = "slow v norm") %>%
  bind_rows(mike_rows, mike_rows2)

###Add unique identifier
df_long_between_str <- df_long_between_str %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix for cluster-robust variance estimation.
V <- vcalc(vi = df_long_between_str$var_g,
           cluster = df_long_between_str$study,
           obs = df_long_between_str$es_id,
           rho = rho)

###Perform main multi-level meta-analysis in metafor
res_between_str <- rma.mv(yi = g, 
                            V = V,
                            random = ~ 1 | study / es_id, 
                            data = df_long_between_str,
                            method = "REML",
                            test = "t")

###Run main multi-level meta through clubSandich for cluster robust estimation
res_between_str_robust <- robust(res_between_str, cluster= study, clubSandwich=TRUE)

###Print main results
summary(res_between_str_robust)

###90%CIs
estimate <- res_between_str_robust$beta
se <- res_between_str_robust$se
deg_f <- res_between_str_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% prediction intervals
pred_int_between_str <- predict(res_between_str_robust, addx = TRUE)
print(pred_int_between_str)

###90% PIs
predicted_effect <- res_between_str_robust$b
se_predicted_effect <- res_between_str_robust$se
between_study_variance <- sum(res_between_str_robust$sigma2)
deg_f <- res_between_str_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_between_str_robust$vi)
X <- model.matrix(res_between_str_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_between_str_robust$sigma2) / (sum(res_between_str_robust$sigma2) + (res_between_str_robust$k-res_between_str_robust$p)/sum(diag(P)))
100 * res_between_str_robust$sigma2 / (sum(res_between_str_robust$sigma2) + (res_between_str_robust$k-res_between_str_robust$p)/sum(diag(P)))

###Tau2
confint <- confint(res_between_str_robust)
print(confint)

###Move on to selecting within-participant design studies
data_within_str <- subset(df, outcome == "strength" & study_design == "within")

###Assign unique ids
data_within_str <- data_within_str %>%
  mutate(es_id = 9:10)

###Calculate effect sizes with variances
data_within_str <- data_within_str %>%
  mutate(
    # 4.4
    Spooled = sqrt(((norm_pre_sd^2 + norm_post_sd^2) / 2 + (slow_pre_sd^2 + slow_post_sd^2) / 2) / 2),
    # 4.15
    Sdiff = sqrt((norm_delta_sd^2 + slow_delta_sd^2) - 2 * between_arm_within_subject_corr * norm_delta_sd * slow_delta_sd),
    # We use pre post correlations per Eq 4.27. I've elected to average the slow and norm correlations for this.
    avg_corr = (norm_pre_post_cor + slow_pre_post_cor) / 2,
    # 4.27
    Swithin = Sdiff / sqrt(2 * (1 - (avg_corr))),
    # 4.26
    d = (slow_delta_mean - norm_delta_mean) / Swithin,
    # 4.28
    var_d = (1 / slow_n) + (d^2 / (2 * slow_n)) * 2 * (1 - (avg_corr)),
    # 4.22
    J = 1 - (3 / (4 * (slow_n - 1) - 1)),
    # 4.23
    g = d * J,
    # 4.24
    var_g = (J^2 * var_d),
    # 4.25
    se_g = sqrt(var_g) 
  )

###Calculate variance covariance matrix for cluster-robust variance estimation.
V <- vcalc(vi = data_within_str$var_g,
           cluster = data_within_str$study,
           obs = data_within_str$es_id,
           rho = rho)

res_within_str <- rma.mv(yi = g, 
                             V = V, 
                             random = ~ 1 | study / es_id,  
                             data = data_within_str,
                             method = "REML",
                         test = "t")

res_within_str_robust <- robust(res_within_str, cluster= data_within_str$study, clubSandwich=TRUE)

summary(res_within_str_robust)

###Calculate 90%CIs
estimate <- res_within_str_robust$beta
se <- res_within_str_robust$se
deg_f <- res_within_str_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% prediction intervals
pred_int_within_str <- predict(res_within_str_robust, addx = TRUE)
print(pred_int_within_str)

###90% prediction intervals
predicted_effect <- res_within_str_robust$b
se_predicted_effect <- res_within_str_robust$se
between_study_variance <- sum(res_within_str_robust$sigma2)
deg_f <- res_within_str_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_within_str_robust$vi)
X <- model.matrix(res_within_str_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_within_str_robust$sigma2) / (sum(res_within_str_robust$sigma2) + (res_within_str_robust$k-res_within_str_robust$p)/sum(diag(P)))
100 * res_within_str_robust$sigma2 / (sum(res_within_str_robust$sigma2) + (res_within_str_robust$k-res_within_str_robust$p)/sum(diag(P)))

###tau2
confint <- confint(res_within_str_robust)
print(confint)

###Pooling between and within designs for strength outcome.
data_combined_str <- bind_rows(df_long_between_str, data_within_str)

###Variance covariance matrix
V <- vcalc(vi = data_combined_str$var_g,
           cluster = data_combined_str$study,
           obs = data_combined_str$es_id,
           rho = rho)

res_pooled_str <- rma.mv(yi = g, 
                         V = V, 
                         random = ~ 1 | study / es_id,  
                         data = data_combined_str,
                         method = "REML",
                         test = "t")

res_pooled_str_robust <- robust(res_pooled_str, cluster=study, clubSandwich=TRUE)

summary(res_pooled_str_robust)

###90%CIs
estimate <- res_pooled_str_robust$beta
se <- res_pooled_str_robust$se
deg_f <- res_pooled_str_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_pooled_str <- predict(res_pooled_str_robust, addx = TRUE)
print(pred_int_pooled_str)

###90% PIs
predicted_effect <- res_pooled_str_robust$b
se_predicted_effect <- res_pooled_str_robust$se
between_study_variance <- sum(res_pooled_str_robust$sigma2)
deg_f <- res_pooled_str_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_pooled_str_robust$vi)
X <- model.matrix(res_pooled_str_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_pooled_str_robust$sigma2) / (sum(res_pooled_str_robust$sigma2) + (res_pooled_str_robust$k-res_pooled_str_robust$p)/sum(diag(P)))
100 * res_pooled_str_robust$sigma2 / (sum(res_pooled_str_robust$sigma2) + (res_pooled_str_robust$k-res_pooled_str_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_str_pooled <- confint(res_pooled_str_robust)
print(confint_tau2_str_pooled)

###Sub-group analysis: training status
###Sub-group analysis for untrained (ut)
data_ut <- subset(data_combined_str, training_status == "ut")
data_ut <- data_ut %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix.
V <- vcalc(vi = data_ut$var_g,
           cluster = data_ut$study,
           obs = data_ut$es_id,
           rho = rho)

res_ut <- rma.mv(yi = g, 
                 V = V, 
                 random = ~ 1 | study / es_id, 
                 data = data_ut,
                 method = "REML",
                 test = "t")

res_ut_robust <- robust(res_ut, cluster=study, clubSandwich=TRUE)

summary(res_ut_robust)

###90%CIs
estimate <- res_ut_robust$beta
se <- res_ut_robust$se
deg_f <- res_ut_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_ut <- predict(res_ut_robust, addx = TRUE)
print(pred_int_ut)

###90% PI
predicted_effect <- res_ut_robust$b
se_predicted_effect <- res_ut_robust$se
between_study_variance <- sum(res_ut_robust$sigma2)
deg_f <- res_ut_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_ut_robust$vi)
X <- model.matrix(res_ut_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_ut_robust$sigma2) / (sum(res_ut_robust$sigma2) + (res_ut_robust$k-res_ut_robust$p)/sum(diag(P)))
100 * res_ut_robust$sigma2 / (sum(res_ut_robust$sigma2) + (res_ut_robust$k-res_ut_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_ut <- confint(res_ut_robust)
print(confint_tau2_ut)

###Subgroup analysis for trained (t) individuals
data_t <- subset(data_combined_str, training_status == "t")
data_t <- data_t %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix.
V <- vcalc(vi = data_t$var_g,
           cluster = data_t$study,
           obs = data_t$es_id,
           rho = rho)

res_t <- rma.mv(yi = g, 
                        V = V,
                       random = ~ 1 | study / es_id, 
                       data = data_t,
                       method = "REML",
                       test = "t")

res_t_robust <- robust(res_t, cluster=study, clubSandwich=TRUE)

summary(res_t_robust)

###90%CIs
estimate <- res_t_robust$beta
se <- res_t_robust$se
deg_f <- res_t_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_t <- predict(res_t_robust, addx = TRUE)
print(pred_int_t)

###90% PIs
predicted_effect <- res_t_robust$b
se_predicted_effect <- res_t_robust$se
between_study_variance <- sum(res_t_robust$sigma2)
deg_f <- res_t_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_t_robust$vi)
X <- model.matrix(res_t_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_t_robust$sigma2) / (sum(res_t_robust$sigma2) + (res_t_robust$k-res_t_robust$p)/sum(diag(P)))
100 * res_t_robust$sigma2 / (sum(res_t_robust$sigma2) + (res_t_robust$k-res_t_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_t <- confint(res_t_robust)
print(confint_tau2_t)

###Sub-group analysis: multi v single joint
###Sub-group analysis for multijoint (mj)
data_mj <- subset(data_combined_str, multi_joint == "mj")
data_mj <- data_mj %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix.
V <- vcalc(vi = data_mj$var_g,
           cluster = data_mj$study,
           obs = data_mj$es_id,
           rho = rho)

res_mj <- rma.mv(yi = g, 
                V = V,
                random = ~ 1 | study / es_id, 
                data = data_mj,
                method = "REML",
                test = "t")

res_mj_robust <- robust(res_mj, cluster=study, clubSandwich=TRUE)

summary(res_mj_robust)

###90% CIs
estimate <- res_mj_robust$beta
se <- res_mj_robust$se
deg_f <- res_mj_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_mj <- predict(res_mj_robust, addx = TRUE)
print(pred_int_mj)

###90% PIs
predicted_effect <- res_mj_robust$b
se_predicted_effect <- res_mj_robust$se
between_study_variance <- sum(res_mj_robust$sigma2)
deg_f <- res_mj_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_mj_robust$vi)
X <- model.matrix(res_mj_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_mj_robust$sigma2) / (sum(res_mj_robust$sigma2) + (res_mj_robust$k-res_mj_robust$p)/sum(diag(P)))
100 * res_mj_robust$sigma2 / (sum(res_mj_robust$sigma2) + (res_mj_robust$k-res_mj_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_mj_robust)
print(confint_tau2)

###Subgroup analysis for single-joint (sj) exercises
data_sj <- subset(data_combined_str, multi_joint == "sj")
data_sj <- data_sj %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_sj$var_g,
           cluster = data_sj$study,
           obs = data_sj$es_id,
           rho = rho)

res_sj <- rma.mv(yi = g, 
                 V = V,
                 random = ~ 1 | study / es_id, 
                 data = data_sj,
                 method = "REML",
                 test = "t")

res_sj_robust <- robust(res_sj, cluster=study, clubSandwich=TRUE)

summary(res_sj_robust)

###90% CIs
estimate <- res_sj_robust$beta
se <- res_sj_robust$se
deg_f <- res_sj_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_sj <- predict(res_sj_robust, addx = TRUE)
print(pred_int_sj)

###90% PIs
predicted_effect <- res_sj_robust$b
se_predicted_effect <- res_sj_robust$se
between_study_variance <- sum(res_sj_robust$sigma2)
deg_f <- res_sj_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_sj_robust$vi)
X <- model.matrix(res_sj_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_sj_robust$sigma2) / (sum(res_sj_robust$sigma2) + (res_sj_robust$k-res_sj_robust$p)/sum(diag(P)))
100 * res_sj_robust$sigma2 / (sum(res_sj_robust$sigma2) + (res_sj_robust$k-res_sj_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_sj_robust)
print(confint_tau2)

###Sub-group analysis: vl v stf
###Sub-group analysis for volume load (vl)
data_vl <- subset(data_combined_str, match_vl_stf == "vl")
data_vl <- data_vl %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_vl$var_g,
           cluster = data_vl$study,
           obs = data_vl$es_id,
           rho = rho)

###Issue with model convergence so Nelder-Mead optimizer is used.
res_vl <- rma.mv(yi = g, 
                 V = V,
                 random = ~ 1 | study / es_id, 
                 data = data_vl,
                 method = "REML",
                 test = "t",
                 control = list(optimizer = "optim", optmethod = "Nelder-Mead"))

res_vl_robust <- robust(res_vl, cluster=study, clubSandwich=TRUE)

summary(res_vl_robust)

###90% CIs
estimate <- res_vl_robust$beta
se <- res_vl_robust$se
deg_f <- res_vl_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_vl <- predict(res_vl_robust, addx = TRUE)
print(pred_int_vl)

###90% PI
predicted_effect <- res_vl_robust$b
se_predicted_effect <- res_vl_robust$se
between_study_variance <- sum(res_vl_robust$sigma2)
deg_f <- res_vl_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_vl_robust$vi)
X <- model.matrix(res_vl_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_vl_robust$sigma2) / (sum(res_vl_robust$sigma2) + (res_vl_robust$k-res_vl_robust$p)/sum(diag(P)))
100 * res_vl_robust$sigma2 / (sum(res_vl_robust$sigma2) + (res_vl_robust$k-res_vl_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_vl <- confint(res_vl_robust)
print(confint_tau2_vl)

###Subgroup analysis for sets to failure (stf)
data_stf <- subset(data_combined_str, match_vl_stf == "stf")
data_stf <- data_stf %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_stf$var_g,
           cluster = data_stf$study,
           obs = data_stf$es_id,
           rho = rho)

res_stf <- rma.mv(yi = g, 
                 V = V,
                 random = ~ 1 | study / es_id, 
                 data = data_stf,
                 method = "REML",
                 test = "t")

res_stf_robust <- robust(res_stf, cluster=study, clubSandwich=TRUE)

summary(res_stf_robust)

###90% CIs
estimate <- res_stf_robust$beta
se <- res_stf_robust$se
deg_f <- res_stf_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_stf <- predict(res_stf_robust, addx = TRUE)
print(pred_int_stf)

###90% PIs
predicted_effect <- res_stf_robust$b
se_predicted_effect <- res_stf_robust$se
between_study_variance <- sum(res_stf_robust$sigma2)
deg_f <- res_stf_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_stf_robust$vi)
X <- model.matrix(res_stf_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_stf_robust$sigma2) / (sum(res_stf_robust$sigma2) + (res_stf_robust$k-res_stf_robust$p)/sum(diag(P)))
100 * res_stf_robust$sigma2 / (sum(res_stf_robust$sigma2) + (res_stf_robust$k-res_stf_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_stf <- confint(res_stf_robust)
print(confint_tau2_stf)

###Moderator analysis
###Calculate variance covariance matrix for the combined dataset for strength
V_full <- vcalc(vi = data_combined_str$var_g,
                cluster = data_combined_str$study,
                obs = data_combined_str$es_id,
                rho = rho)

###Moderator analysis for training status
res_mod_str_trainingstatus <- rma.mv(yi = g, 
                                     V = V_full, 
                                     mods = ~ training_status, 
                                     random = ~ 1 | study / measure, 
                                     data = data_combined_str,
                                     method = "REML",
                                     test = "t")

res_mod_str_trainingstatus_robust <- robust(res_mod_str_trainingstatus, cluster=study, clubSandwich=TRUE)

summary(res_mod_str_trainingstatus_robust)

###90% CIs
estimate <- res_mod_str_trainingstatus_robust$beta
se <- res_mod_str_trainingstatus_robust$se
deg_f <- res_mod_str_trainingstatus_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Moderator analysis for study design
res_mod_str_design <- rma.mv(yi = g, 
                             V = V_full, 
                             mods = ~ study_design, 
                             random = ~ 1 | study / measure, 
                             data = data_combined_str,
                             method = "REML",
                             test = "t")

res_mod_str_design_robust <- robust(res_mod_str_design, cluster=study, clubSandwich=TRUE)

summary(res_mod_str_design_robust)

###90% CIs
estimate <- res_mod_str_design_robust$beta
se <- res_mod_str_design_robust$se
deg_f <- res_mod_str_design_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90
                     
###Moderator analysis for volume load/sets to failure
res_mod_str_vlstf <- rma.mv(yi = g, 
                            V = V_full, 
                            mods = ~ match_vl_stf, 
                            random = ~ 1 | study / measure, 
                            data = data_combined_str,
                            method = "REML",
                            test = "t")

res_mod_str_vlstf_robust <- robust(res_mod_str_vlstf, cluster=study, clubSandwich=TRUE)

summary(res_mod_str_vlstf_robust)

###90% CIs
estimate <- res_mod_str_vlstf_robust$beta
se <- res_mod_str_vlstf_robust$se
deg_f <- res_mod_str_vlstf_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Moderator analysis for multi-joint and single-joint
res_mod_str_mj_sj <- rma.mv(yi = g, 
                            V = V_full, 
                            mods = ~ multi_joint, 
                            random = ~ 1 | study / measure, 
                            data = data_combined_str,
                            method = "REML",
                            test = "t")

res_mod_str_mj_sj_robust <- robust(res_mod_str_mj_sj, cluster=study, clubSandwich=TRUE)

summary(res_mod_str_mj_sj_robust)

###90% CIs
estimate <- res_mod_str_mj_sj_robust$beta
se <- res_mod_str_mj_sj_robust$se
deg_f <- res_mod_str_mj_sj_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Model diagnostics for strength outcome.
###Calculating Cooks distance
k <- length(res_pooled_str_robust$beta)
n <- nrow(data_combined_str)
f_threshold <- qf(0.5, df1 = k + 1, df2 = n - k - 1)
print(f_threshold)

cooks_dist <- cooks.distance(res_pooled_str_robust, cluster = study)
plot(cooks_dist, type = "o", pch = 19, main="Cook's Distance", xlab="Study", ylab="Cook's Distance")

###Hat values
hat_values <- hatvalues(res_pooled_str_robust, cluster = study)
plot(hat_values, type = "o", pch = 19, main="Hat Values (Leverage)", xlab="Study", ylab="Hat Values")

k <- length(res_pooled_str_robust$beta)
n <- nrow(data_combined_str)
average_hat <- (k + 1) / n
threshold <- 2 * average_hat

high_leverage_points <- which(hat_values > threshold)
cat("Average hat value:", average_hat, "\n")
cat("Threshold for high leverage points:", threshold, "\n")
###No issues with Cook's distance or hat values for strength, all studies are within thresholds.

###Publication bias for strength
###Multi-level extension of funnel plot and Egger's test
residuals <- residuals(res_pooled_str_robust)
sei <- sqrt(res_pooled_str_robust$vi)
standardized_residuals <- residuals / sei

###Standard, contour enhanced funnel plot
cols <- palette.colors(length(unique(data_combined_str$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(data_combined_str$descriptor))]
funnel(res_pooled_str_robust, xlab = "Effect Size (Hedge's g)", col = cols)
funnel(res_pooled_str_robust, xlab = "Effect Size (Hedge's g)", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(-2, 2, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(res_pooled_str_robust$yi) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(res_pooled_str_robust$yi) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Standardised residuals, contour enhanced funnel plot
residuals_rma <- rma(yi = standardized_residuals, sei = sei, method = "FE")
cols <- palette.colors(length(unique(data_combined_str$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(data_combined_str$descriptor))]
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols)
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(min(standardized_residuals) - 1, max(standardized_residuals) + 1, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Multi-level extension of Egger's test for strength
egger_str <- rma.mv(yi = residuals, 
                    V = diag(sei^2),
                    mods = ~ sei, 
                    random = ~ 1 | study / g, 
                    data = data_combined_str)

summary(egger_str)

###Hypertrophy
###Select between-participant designs
data_between_hyp <- subset(df, outcome == "hypertrophy" & study_design == "between")

###Calculate effect sizes and variances
data_between_hyp <- data_between_hyp %>%
  mutate(
    # 4.19
    SDpooled = sqrt(((slow_n - 1) * slow_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow_n + norm_n - 2)),
    # 4.18
    d = (slow_delta_mean - norm_delta_mean) / (SDpooled),
    # 4.22
    j = (1 - (3 / (4 * (slow_n + norm_n - 2) - 1))),
    # 4.23
    g = (d * j),
    # 4.20
    var_d = (slow_n + norm_n) / (slow_n * norm_n) + (d^2 / (2 * (slow_n + norm_n))),
    # 4.24
    var_g = ((j^2) * (var_d)),
    # 4.25
    se_g = sqrt(var_g)
  )

###Add unique identifier
data_between_hyp <- data_between_hyp %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_between_hyp$var_g,
           cluster = data_between_hyp$study,
           obs = data_between_hyp$es_id,
           rho = rho)

res_between_hyp <- rma.mv(yi = g, 
                          V = V,
                          random = ~ 1 | study / es_id, 
                          data = data_between_hyp,
                          method = "REML",
                          test = "t")

res_between_hyp_robust <- robust(res_between_hyp, cluster=study, clubSandwich=TRUE)

summary(res_between_hyp_robust)

###90% CIs
estimate <- res_between_hyp_robust$beta
se <- res_between_hyp_robust$se
deg_f <- res_between_hyp_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_between_hyp <- predict(res_between_hyp_robust, addx = TRUE)
print(pred_int_between_hyp)

###90% PIs
predicted_effect <- res_between_hyp_robust$b
se_predicted_effect <- res_between_hyp_robust$se
between_study_variance <- sum(res_between_hyp_robust$sigma2)
deg_f <- res_between_hyp_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_between_hyp_robust$vi)
X <- model.matrix(res_between_hyp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_between_hyp_robust$sigma2) / (sum(res_between_hyp_robust$sigma2) + (res_between_hyp_robust$k-res_between_hyp_robust$p)/sum(diag(P)))
100 * res_between_hyp_robust$sigma2 / (sum(res_between_hyp_robust$sigma2) + (res_between_hyp_robust$k-res_between_hyp_robust$p)/sum(diag(P)))

###Tau2
confint_tau_between <- confint(res_between_hyp_robust)
print(confint_tau_between)

###Within studies, hypertrophy outcome
data_within_hyp <- subset(df, outcome == "hypertrophy" & study_design == "within")

###Calculate effect sizes and variances
data_within_hyp <- data_within_hyp %>%
  mutate(
    # 4.4
    Spooled = sqrt(((norm_pre_sd^2 + norm_post_sd^2) / 2 + (slow_pre_sd^2 + slow_post_sd^2) / 2) / 2),
    # 4.15
    Sdiff = sqrt((norm_delta_sd^2 + slow_delta_sd^2) - 2 * between_arm_within_subject_corr * norm_delta_sd * slow_delta_sd),
    # We use pre post correlations per Eq 4.27. I've elected to average the slow and norm correlations for this.
    avg_corr = (norm_pre_post_cor + slow_pre_post_cor) / 2,
    # 4.27
    Swithin = Sdiff / sqrt(2 * (1 - (avg_corr))),
    # 4.26
    d = (slow_delta_mean - norm_delta_mean) / Swithin,
    # 4.28
    var_d = (1 / slow_n) + (d^2 / (2 * slow_n)) * 2 * (1 - (avg_corr)),
    # 4.22
    J = 1 - (3 / (4 * (slow_n - 1) - 1)),
    # 4.23
    g = d * J,
    # 4.24
    var_g = (J^2 * var_d),
    # 4.25
    se_g = sqrt(var_g) 
  )

###Add unique identifier
data_within_hyp <- data_within_hyp %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_within_hyp$var_g,
           cluster = data_within_hyp$study,
           obs = data_within_hyp$es_id,
           rho = rho)

res_within_hyp <- rma.mv(yi = g, 
                          V = V,
                          random = ~ 1 | study / es_id, 
                          data = data_within_hyp,
                          method = "REML",
                          test = "t")

res_within_hyp_robust <- robust(res_within_hyp, cluster=study, clubSandwich=TRUE)

summary(res_within_hyp_robust)

###90% CIs
estimate <- res_within_hyp_robust$beta
se <- res_within_hyp_robust$se
deg_f <- res_within_hyp_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_within_hyp <- predict(res_within_hyp_robust, addx = TRUE)
print(pred_int_within_hyp)

###90% PIs
predicted_effect <- res_within_hyp_robust$b
se_predicted_effect <- res_within_hyp_robust$se
between_study_variance <- sum(res_within_hyp_robust$sigma2)
deg_f <- res_within_hyp_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_within_hyp_robust$vi)
X <- model.matrix(res_within_hyp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_within_hyp_robust$sigma2) / (sum(res_within_hyp_robust$sigma2) + (res_within_hyp_robust$k-res_within_hyp_robust$p)/sum(diag(P)))
100 * res_within_hyp_robust$sigma2 / (sum(res_within_hyp_robust$sigma2) + (res_within_hyp_robust$k-res_within_hyp_robust$p)/sum(diag(P)))

###Tau2
confint_tau_within <- confint(res_within_hyp_robust)
print(confint_tau_within)

###Pooling between and within designs for hypertrophy
###Reassign es_id values
data_between_hyp_selected <- select(data_between_hyp, -es_id)
data_within_hyp_selected <- select(data_within_hyp, -es_id)

###Combine the datasets
data_combined_hyp <- bind_rows(data_between_hyp_selected, data_within_hyp_selected)

###Add unique identifier again
data_combined_hyp <- data_combined_hyp %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_combined_hyp$var_g,
           cluster = data_combined_hyp$study,
           obs = data_combined_hyp$es_id,
           rho = rho)

res_pooled_hyp <- rma.mv(yi = g, 
                             V = V, 
                             random = ~ 1 | study / es_id,  
                             data = data_combined_hyp,
                             method = "REML",
                         test = "t")

res_pooled_hyp_robust <- robust(res_pooled_hyp, cluster=study, clubSandwich=TRUE)

summary(res_pooled_hyp_robust)

###90% CIs
estimate <- res_pooled_hyp_robust$beta
se <- res_pooled_hyp_robust$se
deg_f <- res_pooled_hyp_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95%PIs
pred_int_pooled_hyp <- predict(res_pooled_hyp_robust, addx = TRUE)
print(pred_int_pooled_hyp)

###90% PI
predicted_effect <- res_pooled_hyp_robust$b
se_predicted_effect <- res_pooled_hyp_robust$se
between_study_variance <- sum(res_pooled_hyp_robust$sigma2)
deg_f <- res_pooled_hyp_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_pooled_hyp_robust$vi)
X <- model.matrix(res_pooled_hyp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_pooled_hyp_robust$sigma2) / (sum(res_pooled_hyp_robust$sigma2) + (res_pooled_hyp_robust$k-res_pooled_hyp_robust$p)/sum(diag(P)))
100 * res_pooled_hyp_robust$sigma2 / (sum(res_pooled_hyp_robust$sigma2) + (res_pooled_hyp_robust$k-res_pooled_hyp_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_hyp <- confint(res_pooled_hyp_robust)
print(confint_tau2_hyp)

###Sub-group analysis: training status
### Subgroup analysis for untrained (ut) individuals
data_hyp_ut <- subset(data_combined_hyp, training_status == "ut")
data_hyp_ut <- data_hyp_ut %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_hyp_ut$var_g,
           cluster = data_hyp_ut$study,
           obs = data_hyp_ut$es_id,
           rho = rho)

res_hyp_ut <- rma.mv(yi = g, 
                 V = V, 
                 random = ~ 1 | study / es_id, 
                 data = data_hyp_ut,
                 method = "REML",
                 test = "t")

res_hyp_ut_robust <- robust(res_hyp_ut, cluster=study, clubSandwich=TRUE)

summary(res_hyp_ut_robust)

###90% CIs
estimate <- res_hyp_ut_robust$beta
se <- res_hyp_ut_robust$se
deg_f <- res_hyp_ut_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_hyp_ut <- predict(res_hyp_ut_robust, addx = TRUE)
print(pred_int_hyp_ut)

###90% PIs
predicted_effect <- res_hyp_ut_robust$b
se_predicted_effect <- res_hyp_ut_robust$se
between_study_variance <- sum(res_hyp_ut_robust$sigma2)
deg_f <- res_hyp_ut_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_hyp_ut_robust$vi)
X <- model.matrix(res_hyp_ut_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_hyp_ut_robust$sigma2) / (sum(res_hyp_ut_robust$sigma2) + (res_hyp_ut_robust$k-res_hyp_ut_robust$p)/sum(diag(P)))
100 * res_hyp_ut_robust$sigma2 / (sum(res_hyp_ut_robust$sigma2) + (res_hyp_ut_robust$k-res_hyp_ut_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_hyp_ut_robust)
print(confint_tau2)

###Subgroup analysis for trained (t) individuals
data_hyp_t <- subset(data_combined_hyp, training_status == "t")
data_hyp_t <- data_hyp_t %>%
  mutate(es_id = row_number())

###Robust variance
V <- vcalc(vi = data_hyp_t$var_g,
           cluster = data_hyp_t$study,
           obs = data_hyp_t$es_id,
           rho = rho)

res_hyp_t <- rma.mv(yi = g, 
                     V = V, 
                     random = ~ 1 | study / es_id, 
                     data = data_hyp_t,
                     method = "REML",
                     test = "t")

res_hyp_t_robust <- robust(res_hyp_t, cluster=study, clubSandwich=TRUE)

summary(res_hyp_t_robust)

###90% CIs
estimate <- res_hyp_t_robust$beta
se <- res_hyp_t_robust$se
deg_f <- res_hyp_t_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_hyp_t <- predict(res_hyp_t_robust, addx = TRUE)
print(pred_int_hyp_t)

###90% PIs
predicted_effect <- res_hyp_t_robust$b
se_predicted_effect <- res_hyp_t_robust$se
between_study_variance <- sum(res_hyp_t_robust$sigma2)
deg_f <- res_hyp_t_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_hyp_t_robust$vi)
X <- model.matrix(res_hyp_t_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_hyp_t_robust$sigma2) / (sum(res_hyp_t_robust$sigma2) + (res_hyp_t_robust$k-res_hyp_t_robust$p)/sum(diag(P)))
100 * res_hyp_t_robust$sigma2 / (sum(res_hyp_t_robust$sigma2) + (res_hyp_t_robust$k-res_hyp_t_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_hyp_t_robust)
print(confint_tau2)

###Sub-group analysis: arm v thigh
###Sub-group analysis for arm
data_hyp_arm <- subset(data_combined_hyp, muscle_group == "arm")
data_hyp_arm <- data_hyp_arm %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_hyp_arm$var_g,
           cluster = data_hyp_arm$study,
           obs = data_hyp_arm$es_id,
           rho = rho)

res_hyp_arm <- rma.mv(yi = g, 
                    V = V, 
                    random = ~ 1 | study / es_id, 
                    data = data_hyp_arm,
                    method = "REML",
                    test = "t")

res_hyp_arm_robust <- robust(res_hyp_arm, cluster=study, clubSandwich=TRUE)

summary(res_hyp_arm_robust)

###90% CIs
estimate <- res_hyp_arm_robust$beta
se <- res_hyp_arm_robust$se
deg_f <- res_hyp_arm_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_hyp_arm <- predict(res_hyp_arm_robust, addx = TRUE)
print(pred_int_hyp_arm)

###90% PIs
predicted_effect <- res_hyp_arm_robust$b
se_predicted_effect <- res_hyp_arm_robust$se
between_study_variance <- sum(res_hyp_arm_robust$sigma2)
deg_f <- res_hyp_arm_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_hyp_arm_robust$vi)
X <- model.matrix(res_hyp_arm_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_hyp_arm_robust$sigma2) / (sum(res_hyp_arm_robust$sigma2) + (res_hyp_arm_robust$k-res_hyp_arm_robust$p)/sum(diag(P)))
100 * res_hyp_arm_robust$sigma2 / (sum(res_hyp_arm_robust$sigma2) + (res_hyp_arm_robust$k-res_hyp_arm_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_hyp_arm_robust)
print(confint_tau2)

###Subgroup analysis for muscle group thigh (thigh) individuals
data_hyp_thigh <- subset(data_combined_hyp, muscle_group == "thigh")
data_hyp_thigh <- data_hyp_thigh %>%
  mutate(es_id = row_number())

###Calculate variance covariance matrix
V <- vcalc(vi = data_hyp_thigh$var_g,
           cluster = data_hyp_thigh$study,
           obs = data_hyp_thigh$es_id,
           rho = rho)

res_hyp_thigh <- rma.mv(yi = g, 
                      V = V, 
                      random = ~ 1 | study / es_id, 
                      data = data_hyp_thigh,
                      method = "REML",
                      test = "t")

res_hyp_thigh_robust <- robust(res_hyp_thigh, cluster=study, clubSandwich=TRUE)

summary(res_hyp_thigh_robust)

###90% CIs
estimate <- res_hyp_thigh_robust$beta
se <- res_hyp_thigh_robust$se
deg_f <- res_hyp_thigh_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_hyp_thigh <- predict(res_hyp_thigh_robust, addx = TRUE)
print(pred_int_hyp_thigh)

###90% PIs
predicted_effect <- res_hyp_thigh_robust$b
se_predicted_effect <- res_hyp_thigh_robust$se
between_study_variance <- sum(res_hyp_thigh_robust$sigma2)
deg_f <- res_hyp_thigh_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_hyp_thigh_robust$vi)
X <- model.matrix(res_hyp_thigh_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_hyp_thigh_robust$sigma2) / (sum(res_hyp_thigh_robust$sigma2) + (res_hyp_thigh_robust$k-res_hyp_thigh_robust$p)/sum(diag(P)))
100 * res_hyp_thigh_robust$sigma2 / (sum(res_hyp_thigh_robust$sigma2) + (res_hyp_thigh_robust$k-res_hyp_thigh_robust$p)/sum(diag(P)))

###Tau2
confint_tau2 <- confint(res_hyp_thigh)
print(confint_tau2)

###Moderator analysis for study design
V_hyp_full <- vcalc(vi = data_combined_hyp$var_g,
                    cluster = data_combined_hyp$study,
                    obs = data_combined_hyp$es_id,
                    rho = rho)

res_mod_hyp_design <- rma.mv(yi = g, 
                             V = V_hyp_full, 
                             mods = ~ study_design, 
                             random = ~ 1 | study / measure, 
                             data = data_combined_hyp,
                             method = "REML",
                             test = "t")

res_mod_hyp_design_robust <- robust(res_mod_hyp_design, cluster=study, clubSandwich=TRUE)

summary(res_mod_hyp_design_robust)

###90% CIs
estimate <- res_mod_hyp_design_robust$beta
se <- res_mod_hyp_design_robust$se
deg_f <- res_mod_hyp_design_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Moderator analysis for training status
res_mod_hyp_trainingstatus <- rma.mv(yi = g, 
                                     V = V_hyp_full, 
                                     mods = ~ training_status, 
                                     random = ~ 1 | study / measure, 
                                     data = data_combined_hyp,
                                     method = "REML",
                                     test = "t")

res_mod_hyp_trainingstatus_robust <- robust(res_mod_hyp_trainingstatus, cluster=study, clubSandwich=TRUE)

summary(res_mod_hyp_trainingstatus_robust)

###90% CIs
estimate <- res_mod_hyp_trainingstatus_robust$beta
se <- res_mod_hyp_trainingstatus_robust$se
deg_f <- res_mod_hyp_trainingstatus_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Moderator analysis for arm v thigh
res_mod_hyp_armthigh <- rma.mv(yi = g, 
                                     V = V_hyp_full, 
                                     mods = ~ muscle_group, 
                                     random = ~ 1 | study / measure, 
                                     data = data_combined_hyp,
                                     method = "REML",
                                     test = "t")

res_mod_hyp_armthigh_robust <- robust(res_mod_hyp_armthigh, cluster=study, clubSandwich=TRUE)

summary(res_mod_hyp_armthigh_robust)

###90% CIs
estimate <- res_mod_hyp_armthigh_robust$beta
se <- res_mod_hyp_armthigh_robust$se
deg_f <- res_mod_hyp_armthigh_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###Model diagnostics
###Calculate Cooks distance
k <- length(res_pooled_hyp_robust$beta)
n <- nrow(data_combined_hyp)
f_threshold <- qf(0.5, df1 = k + 1, df2 = n - k - 1)
print(f_threshold)

cooks_dist <- cooks.distance(res_pooled_hyp_robust, cluster = study)
plot(cooks_dist, type = "o", pch = 19, main="Cook's Distance", xlab="Study", ylab="Cook's Distance")
abline(h = f_threshold, col = "red", lty = 2)

###Hat values
hat_values <- hatvalues(res_pooled_hyp_robust, cluster = study)
plot(hat_values, type = "o", pch = 19, main="Hat Values (Leverage)", xlab="Study", ylab="Hat Values")

k <- length(res_pooled_hyp_robust$beta)
n <- nrow(data_combined_hyp)
average_hat <- (k + 1) / n
threshold <- 2 * average_hat

high_leverage_points <- which(hat_values > threshold)
cat("Average hat value:", average_hat, "\n")
cat("Threshold for high leverage points:", threshold, "\n")

###Publication bias for hypertrophy
residuals <- residuals(res_pooled_hyp_robust)
sei <- sqrt(res_pooled_hyp_robust$vi)
standardized_residuals <- residuals / sei

###Standard, contour enhanced funnel plot
cols <- palette.colors(length(unique(data_combined_hyp$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(data_combined_hyp$descriptor))]
funnel(res_pooled_hyp_robust, xlab = "Effect Size (Hedge's g)", col = cols)
funnel(res_pooled_hyp_robust, xlab = "Effect Size (Hedge's g)", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(-2, 2, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(res_pooled_hyp_robust$V) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(res_pooled_hyp_robust$V) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Standardised residuals, contour enhanced funnel plot
residuals_rma <- rma(yi = standardized_residuals, sei = sei, method = "FE")
cols <- palette.colors(length(unique(data_combined_hyp$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(data_combined_hyp$descriptor))]
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols)
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(min(standardized_residuals) - 1, max(standardized_residuals) + 1, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Multi-level extension of Egger's test for strength
egger_hyp <- rma.mv(yi = residuals, 
                    V = diag(sei^2),
                    mods = ~ sei, 
                    random = ~ 1 | study / g, 
                    data = data_combined_hyp)

summary(egger_hyp)

###CMJ
###Select studies for CMJ
data_pow <- subset(df, outcome == "power")

###Calculate effect sizes and variances
data_pow <- data_pow %>%
  mutate(
    # 4.19
    SDpooled = sqrt(((slow_n - 1) * slow_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow_n + norm_n - 2)),
    SDpooled2 = sqrt(((slow2_n - 1) * slow2_delta_sd^2 + (norm_n - 1) * norm_delta_sd^2) / (slow2_n + norm_n - 2)),
    # 4.18
    d = (slow_delta_mean - norm_delta_mean) / (SDpooled),
    d2 = (slow2_delta_mean - norm_delta_mean) / (SDpooled2),
    # 4.22
    j = (1 - (3 / (4 * (slow_n + norm_n - 2) - 1))),
    j2 = (1 - (3 / (4 * (slow2_n + norm_n - 2) - 1))),
    # 4.23
    g = (d * j),
    g2 = (d2 * j2),
    # 4.20
    var_d = (slow_n + norm_n) / (slow_n * norm_n) + (d^2 / (2 * (slow_n + norm_n))),
    var_d2 = (slow2_n + norm_n) / (slow2_n * norm_n) + (d2^2 / (2 * (slow2_n + norm_n))),
    # 4.24
    var_g = ((j^2) * (var_d)),
    var_g2 = ((j2^2) * (var_d2)),
    # 4.25
    se_g = sqrt(var_g),
    se_g2 = sqrt(var_g2)
  )

###Separate rows for Mike et al. study
mike_rows <- data_pow %>%
  filter(descriptor == "Mike et al. (2017)") %>%
  select(study, descriptor, measure, g, var_g, norm_delta_mean, norm_delta_sd, slow_delta_mean, slow_delta_sd, slow2_delta_mean, slow2_delta_sd) %>%
  mutate(Comparison = "slow v norm")

mike_rows2 <- data_pow %>%
  filter(descriptor == "Mike et al. (2017)") %>%
  select(study, descriptor, measure, g2, var_g2, norm_delta_mean, norm_delta_sd, slow_delta_mean, slow_delta_sd, slow2_delta_mean, slow2_delta_sd) %>%
  rename(g = g2, var_g = var_g2) %>%
  mutate(Comparison = "slow2 v norm")

###Combine data into long format
df_long_pow <- data_pow %>%
  filter(descriptor != "Mike et al. (2017)") %>%
  select(study, descriptor, measure, g, var_g, norm_delta_mean, norm_delta_sd, slow_delta_mean, slow_delta_sd, slow2_delta_mean, slow2_delta_sd) %>%
  mutate(Comparison = "slow v norm") %>%
  bind_rows(mike_rows, mike_rows2)

df_long_pow <- df_long_pow %>%
  mutate(es_id = row_number())

View(df_long_pow)

###Calculate variance covariance matrix
V <- vcalc(vi = df_long_pow$var_g,
           cluster = df_long_pow$study,
           obs = df_long_pow$es_id,
           rho = rho)

res_pow <- rma.mv(yi = g, 
                         V = V, 
                         random = ~ 1 | study / es_id,  
                         data = df_long_pow,
                         method = "REML",
                  test = "t")

###Issues with convergence with rho = 0.4 and 0.6, Nelder-Mead optimiser used.

res_pow <- rma.mv(
  yi = g, 
  V = V, 
  random = ~ 1 | study / es_id,  
  data = df_long_pow,
  method = "REML",
  test = "t",
  control = list(
    optimizer = "optim",
    optmethod = "Nelder-Mead",
    maxit = 5000,        
    trace = TRUE         
  )
)

res_pow_robust <- robust(res_pow, cluster=study, clubSandwich=TRUE)

summary(res_pow_robust)

###Calculate 90%CIs
estimate <- res_pow_robust$beta
se <- res_pow_robust$se
deg_f <- res_pow_robust$df
t_crit <- qt(0.95, df = deg_f)
ci_lower_90 <- estimate - t_crit * se
ci_upper_90 <- estimate + t_crit * se
ci_lower_90
ci_upper_90

###95% PIs
pred_int_pow <- predict(res_pow_robust, addx = TRUE)
print(pred_int_pow)

###90% PIs
predicted_effect <- res_pow_robust$b
se_predicted_effect <- res_pow_robust$se
between_study_variance <- sum(res_pow_robust$sigma2)
deg_f <- res_pow_robust$df
critical_value_90 <- qt(0.95, df=deg_f)
se_prediction <- sqrt(se_predicted_effect^2 + between_study_variance)

lower_bound_90 <- predicted_effect - critical_value_90 * se_prediction
upper_bound_90 <- predicted_effect + critical_value_90 * se_prediction
cat("90% Prediction Interval: [", lower_bound_90, ", ", upper_bound_90, "]\n")

###I2
W <- diag(1/res_pow_robust$vi)
X <- model.matrix(res_pow_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_pow_robust$sigma2) / (sum(res_pow_robust$sigma2) + (res_pow_robust$k-res_pow_robust$p)/sum(diag(P)))
100 * res_pow_robust$sigma2 / (sum(res_pow_robust$sigma2) + (res_pow_robust$k-res_pow_robust$p)/sum(diag(P)))

###Tau2
confint_tau2_pow <- confint(res_pow_robust)
print(confint_tau2_pow)

###Multi-level extension of funnel plot and Egger's test
residuals <- residuals(res_pow_robust)
sei <- sqrt(res_pow_robust$vi)
standardized_residuals <- residuals / sei

###Standard, contour enhanced funnel plot
cols <- palette.colors(length(unique(df_long_pow$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(df_long_pow$descriptor))]
funnel(res_pow_robust, xlab = "Effect Size (Hedge's g)", col = cols)
funnel(res_pow_robust, xlab = "Effect Size (Hedge's g)", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(-2, 2, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(res_pow_robust$yi) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(res_pow_robust$yi) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Standardised residuals, contour enhanced funnel plot
residuals_rma <- rma(yi = standardized_residuals, sei = sei, method = "FE")
cols <- palette.colors(length(unique(df_long_pow$descriptor)), palette="Alphabet")
cols <- cols[as.numeric(factor(df_long_pow$descriptor))]
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols)
funnel(residuals_rma, ylab = "Standard Error", xlab = "Standardized Residuals", col = cols, 
       level = c(0.95, 0.99), shade = c("lightgray", "darkgray"))
x <- seq(min(standardized_residuals) - 1, max(standardized_residuals) + 1, length = 100)
y1 <- qnorm((1 - 0.95) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
y2 <- qnorm((1 - 0.99) / 2) / sqrt(1 / (length(residuals_rma$yi) + 1) * (x^2))
lines(x, y1, col = "lightgray", lty = 2)
lines(x, -y1, col = "lightgray", lty = 2)
lines(x, y2, col = "darkgray", lty = 2)
lines(x, -y2, col = "darkgray", lty = 2)

###Multi-level extension of Egger's test for CMJ
egger_pow <- rma.mv(yi = residuals, 
                    V = diag(sei^2),
                    mods = ~ sei, 
                    random = ~ 1 | study / g, 
                    data = df_long_pow)

summary(egger_pow)