### ============================================================
### Orchard plot for the maximal-strength analysis, split by
### training status - an alternative to the paper's Figure 3.
###
### Same idea as an orchard plot (Nakagawa et al., 2021): every
### effect size is a dot sized by its precision, with the pooled
### mean, its confidence interval (thick bar) and its prediction
### interval (thin whisker) for each group. Unlike Figure 3, the
### prediction interval is drawn, and the ROPE is shaded.
###
### Drawn with ggplot2 (not the orchaRd package) so it can use the
### authors' own cluster-robust 90% intervals from run_analysis.R.
### Run from the repo root (or open via the .Rproj).
### ============================================================

library(ggplot2)

### Reuse every model object run_analysis.R already fits.
source("code/run_analysis.R")

### Robust 90% CI and 90% PI, computed exactly as in run_analysis.R
### (t-distribution on the robust df; PI adds the between-study variance).
robust_intervals <- function(robust_obj) {
  est   <- as.numeric(robust_obj$beta)
  se    <- as.numeric(robust_obj$se)
  tcrit <- qt(0.95, df = robust_obj$df)
  pi_se <- sqrt(se^2 + sum(robust_obj$sigma2))
  c(mean = est,
    ci_lo = est - tcrit * se,    ci_hi = est + tcrit * se,
    pi_lo = est - tcrit * pi_se, pi_hi = est + tcrit * pi_se)
}

str_dat <- data_combined_str

### Plot order, top to bottom, matching Figure 3
groups <- list(
  list(label = "Untrained", y = 3, dat = str_dat[str_dat$training_status == "ut", ], fit = res_ut_robust),
  list(label = "Trained",   y = 2, dat = str_dat[str_dat$training_status == "t", ],  fit = res_t_robust),
  list(label = "All studies", y = 1, dat = str_dat,                                  fit = res_pooled_str_robust)
)

pts <- do.call(rbind, lapply(groups, function(g) {
  data.frame(group = g$label, y = g$y, g = g$dat$g, precision = 1 / sqrt(g$dat$var_g))
}))

summ <- do.call(rbind, lapply(groups, function(g) {
  r <- robust_intervals(g$fit)
  data.frame(group = g$label, y = g$y,
             k_effects = nrow(g$dat), k_studies = length(unique(g$dat$descriptor)),
             t(r))
}))
summ$axis_label <- sprintf("%s\n%d effects, %d studies", summ$group, summ$k_effects, summ$k_studies)
summ$note <- sprintf("g = %.2f\n90%% CI [%.2f, %.2f]\n90%% PI [%.2f, %.2f]",
                     summ$mean, summ$ci_lo, summ$ci_hi, summ$pi_lo, summ$pi_hi)

rope <- 0.25   # strength ROPE, as in the paper
x_lim <- c(-3, 3)

p <- ggplot() +
  annotate("rect", xmin = -rope, xmax = rope, ymin = 0.4, ymax = 3.6, fill = "#B23A3A", alpha = 0.10) +
  geom_vline(xintercept = c(-rope, rope), colour = "#B23A3A", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(data = pts, aes(x = g, y = y, size = precision, colour = group),
             alpha = 0.55, position = position_jitter(width = 0, height = 0.14, seed = 7)) +
  geom_segment(data = summ, aes(x = pi_lo, xend = pi_hi, y = y, yend = y), colour = "#1B2A4A", linewidth = 0.55) +
  geom_segment(data = summ, aes(x = ci_lo, xend = ci_hi, y = y, yend = y), colour = "#1B2A4A", linewidth = 2.6, lineend = "butt") +
  geom_point(data = summ, aes(x = mean, y = y), shape = 21, size = 4.2, fill = "white", colour = "#1B2A4A", stroke = 1.3) +
  geom_text(data = summ, aes(x = 3.25, y = y, label = note), hjust = 0, size = 3.1, lineheight = 0.95, colour = "#1B2A4A") +
  scale_y_continuous(breaks = summ$y, labels = summ$axis_label, limits = c(0.4, 3.6)) +
  scale_colour_manual(values = c("Untrained" = "#2A9D8F", "Trained" = "#E9A23B", "All studies" = "#6B7686"), guide = "none") +
  scale_size_continuous(range = c(1.5, 5.5), name = "Precision\n(1 / SE)") +
  scale_x_continuous(breaks = -3:3) +
  coord_cartesian(xlim = x_lim, clip = "off") +
  labs(x = "Hedges' g   (← favors SEG    |    favors LEG →)",
       y = NULL,
       title = "Maximal strength by training status - orchard plot",
       subtitle = "Dots = individual effect sizes | thick bar = 90% CI | thin line = 90% PI | red band = ROPE (±0.25)") +
  theme_minimal(base_size = 11) +
  theme(plot.margin = margin(10, 120, 10, 10),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "#1B2A4A", face = "bold", lineheight = 0.95),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        legend.position = "bottom")

out_dir <- "figures and data"
dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "Orchard_strength_training_status.png"), p, width = 10.5, height = 5.4, dpi = 150, bg = "white")

cat("Orchard plot written to", file.path(out_dir, "Orchard_strength_training_status.png"), "\n")
