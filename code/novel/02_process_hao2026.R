### ============================================================
### Novel approach, step 2 of 4 (Processing)
###
### Turns Hao et al. (2026)'s raw participant data into rows in
### the parent paper's datasheet format, so the parent's own
### effect-size equations can be applied to it unchanged.
###
### Hao M, Qi H, Zhao L, Han W (2026). The effect of squat training
### with different eccentric contraction tempos on lower limb muscle
### strength. PLoS One 21(7): e0354725. doi:10.1371/journal.pone.0354725
### Raw data: the paper's supplementary file S1 (CC BY 4.0), copied
### to data/raw/Hao2026_pone.0354725.s001.xlsx.
###
### Coding decisions, each checked against the parent paper's rules:
###  - 4/0/X/0 group = "slow" (longer eccentric, >2 s);
###    2/0/X/0 group = "norm" (shorter eccentric, <=2 s).
###    The third group chose its own tempo, so it has no fixed
###    eccentric duration and fails the parent's eligibility
###    criterion 6 (all phases except the eccentric matched). Left out.
###  - training_status = "t": Hao required >= 2 years of strength
###    training; the parent counted 1-3 years as "trained".
###  - match_vl_stf = "stf": both groups did the same sets, with
###    repetitions "until failure to maintain tempo" (Hao, Table 2).
###  - multi_joint = "mj" (squat), study_duration = 8 (weeks).
###  - Outcomes: squat 1RM -> strength; CMJ height -> power (the
###    parent's CMJ outcome). SJ and 30-m sprint are not parent outcomes.
###  - Risk of bias was not assessed here (RoB 2 needs two reviewers);
###    it isn't used by the meta-analytic models.
###
### Run from the repo root:  Rscript code/novel/02_process_hao2026.R
### Output: data/processed/hao2026_rows.csv
### ============================================================

if (!file.exists("data/raw/Hao2026_pone.0354725.s001.xlsx"))
  stop("Run this from the repository folder: open DS340w_research-paper.Rproj in RStudio first.",
       call. = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

raw <- read_excel("data/raw/Hao2026_pone.0354725.s001.xlsx",
                  sheet = "Lower Limb Strength Indices",
                  col_names = FALSE, skip = 3, .name_repair = "minimal")
raw <- raw[, 1:9]
names(raw) <- c("id", "group", "body_mass", "height", "time",
                "squat_1rm", "cmj", "sj", "sprint_30m")

### Each participant has a T1 (pre) row followed by a T2 (post) row;
### the id and group are only written on the T1 row. The sheet also has
### two summary rows per group ("Average Value", "Standard Deviation"),
### which are dropped by keeping only rows whose time is exactly T1/T2.
participants <- raw %>%
  mutate(id = suppressWarnings(as.integer(id))) %>%
  tidyr::fill(id, group) %>%
  filter(time %in% c("T1", "T2")) %>%
  mutate(across(c(squat_1rm, cmj), as.numeric))

### one row per participant, pre and post side by side
wide <- participants %>%
  select(group, id, time, squat_1rm, cmj) %>%
  tidyr::pivot_wider(names_from = time, values_from = c(squat_1rm, cmj))

stopifnot(all(table(wide$group) == 10))

summarise_arm <- function(arm, measure) {
  d <- wide %>% filter(group == arm)
  pre  <- d[[paste0(measure, "_T1")]]
  post <- d[[paste0(measure, "_T2")]]
  delta <- post - pre
  tt <- t.test(post, pre, paired = TRUE)
  c(pre_mean = mean(pre), post_mean = mean(post), delta_mean = mean(delta),
    pre_sd = sd(pre), post_sd = sd(post), delta_sd = sd(delta),
    n = length(pre), p_value = tt$p.value, pre_post_cor = cor(pre, post),
    t_value = unname(tt$statistic))
}

make_row <- function(measure, outcome, measure_label, study_id) {
  slow <- summarise_arm("4/0/X/0", measure)
  norm <- summarise_arm("2/0/X/0", measure)
  row <- data.frame(
    descriptor = "Hao et al. (2026)", study = study_id, treatment = 1,
    outcome = outcome, measure = measure_label, study_design = "between",
    between_arm_within_subject_corr = NA_real_)
  for (nm in names(norm)) row[[paste0("norm_", nm)]] <- unname(norm[nm])
  for (nm in names(slow)) row[[paste0("slow_", nm)]] <- unname(slow[nm])
  for (nm in names(slow)) row[[paste0("slow2_", nm)]] <- NA_real_
  row$training_status <- "t"
  row$study_duration  <- 8
  row$match_vl_stf    <- "stf"
  row$multi_joint     <- "mj"
  row$muscle_group    <- NA_character_
  for (i in 1:5) row[[paste0("rob_", i)]] <- NA_character_
  row$rob_overall <- NA_character_
  row
}

### study ids continue the parent sheet's per-outcome numbering
### (strength studies are 1-9, power studies 1-3)
hao_rows <- bind_rows(
  make_row("squat_1rm", "strength", "back_squat", study_id = 10),
  make_row("cmj",       "power",    "cmj",        study_id = 4)
)

### keep the parent sheet's exact column order
parent_cols <- names(read_excel("data/raw/AK_ECC_datasheet_Kojic.xlsx", n_max = 0))
hao_rows <- hao_rows[, parent_cols]

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(hao_rows, "data/processed/hao2026_rows.csv", row.names = FALSE)

### The saved file keeps the parent datasheet's column names, because
### steps 3 and 4 read it with the parent's equations. This summary is
### the plain-English version to read: 2-s group = shorter eccentric
### phase (2/0/X/0), 4-s group = longer eccentric phase (4/0/X/0).
hao_summary <- hao_rows %>%
  transmute(`Outcome` = c(strength = "Squat 1RM (kg)", power = "CMJ height (cm)")[outcome],
            `2-s group: lifters` = norm_n,
            `2-s group: mean change` = norm_delta_mean,
            `2-s group: SD of change` = norm_delta_sd,
            `4-s group: lifters` = slow_n,
            `4-s group: mean change` = slow_delta_mean,
            `4-s group: SD of change` = slow_delta_sd)
print(hao_summary, row.names = FALSE)
cat("Wrote data/processed/hao2026_rows.csv\n")
