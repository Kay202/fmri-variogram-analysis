# =====================================================================
# Formal sign tests on Gaussian variogram residuals
# ---------------------------------------------------------------------
# formal tests on two zones defined a priori from the visual
# residual structure:
#   corner zone  = bins 2 to 4
#   plateau zone = bins 5 to 10
# Bin 1 and bins 11 to 17 are excluded with explicit physical
# interpretations.
# For each zone: exact two-sided binomial sign test, Wilcoxon
# signed-rank test, Clopper-Pearson 95% CI on the proportion of
# positive residuals.
#
# Outputs:
#   CSV
#     priority3_pathA_zone_stats.csv          per-subject zone signs
#     priority3_pathA_test_results.csv        Path A test table
#   PNG/PDF
#     figure_residual_structure.{png,pdf}     two-zone view with SE band
#     figure_residual_failure_modes.{png,pdf} four-zone presentation
#
# Requirements: dplyr, tidyr, ggplot2, readxl, readr, minpack.lm
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2)
library(readxl); library(readr)
suppressPackageStartupMessages(library(minpack.lm))

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
BINNED_XLSX  <- "all_subjects_variogram_binned.xlsx"
ALPHA        <- 0.05

# Zone definitions
CORNER_BINS   <- 2:4      # under-prediction
PLATEAU_BINS  <- 5:10     # over-prediction
EXCLUDED_BIN1 <- 1L       # short-lag over-prediction
RECOVERY_BINS <- 11:17    # long-distance recovery

# Output plot dimensions
PLOT_DPI <- 300

# Zone colours (used in all figures)
ZONE_COLS <- list(
  bin1     = "#9D6FA1",
  corner   = "#7CA37A",
  plateau  = "#C44536",
  recovery = "grey55",
  curve    = "#3D6FB5"
)
TEXT_COLS <- list(
  bin1     = "#5A3C72",
  corner   = "#3F6E3D",
  plateau  = "#7A2A22",
  recovery = "grey25"
)

# ---------------------------------------------------------------------
# 1. Gaussian fitter (self-contained, no external dependencies on the
#    main analysis script)
# ---------------------------------------------------------------------
gaussian_variogram <- function(d, nugget, sill, range) {
  nugget + sill * (1 - exp(-(d / range)^2))
}

fit_gaussian <- function(d, y) {
  if (diff(range(y, na.rm = TRUE)) < 1e-6) return(NULL)
  if (cor(d, y) < 0) return(NULL)
  starts <- list(
    list(nugget = 0,                    sill = max(y) * 0.95, range = 4),
    list(nugget = 0,                    sill = max(y) * 0.90, range = mean(d)),
    list(nugget = 0,                    sill = max(y),         range = 5),
    list(nugget = max(0, min(y) * 0.5), sill = max(y) * 0.95, range = 4)
  )
  for (s in starts) {
    fit <- tryCatch(
      suppressWarnings(nlsLM(
        y ~ gaussian_variogram(d, nugget, sill, range),
        start   = s,
        lower   = c(nugget = 0,   sill = 0,   range = 0.1),
        upper   = c(nugget = Inf, sill = Inf, range = 50),
        control = nls.lm.control(maxiter = 500))),
      error = function(e) NULL)
    if (!is.null(fit)) return(fit)
  }
  NULL
}

# ---------------------------------------------------------------------
# 2. Statistical test helpers
# ---------------------------------------------------------------------
clopper_pearson <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  a <- 1 - conf
  c(if (k == 0) 0 else qbeta(a / 2,     k,     n - k + 1),
    if (k == n) 1 else qbeta(1 - a / 2, k + 1, n - k))
}

# Run the full sign-test triple (binomial, Wilcoxon, CP CI) on a vector
# of values. Returns one row, suitable for binding into a table.
run_sign_tests <- function(values, label) {
  z <- values[is.finite(values)]
  n <- length(z); k <- sum(z > 0)
  if (n == 0) {
    return(tibble(label = label, n = 0, n_positive = NA_integer_,
                  prop_positive = NA_real_, cp_lower = NA_real_,
                  cp_upper = NA_real_, sign_test_p = NA_real_,
                  wilcoxon_p = NA_real_, mean_value = NA_real_,
                  median_value = NA_real_))
  }
  bin <- tryCatch(binom.test(k, n, p = 0.5, alternative = "two.sided"),
                  error = function(e) list(p.value = NA_real_))
  wil <- tryCatch(wilcox.test(z, mu = 0, alternative = "two.sided",
                              exact = FALSE),
                  error = function(e) list(p.value = NA_real_))
  cp <- clopper_pearson(k, n, conf = 1 - ALPHA)
  tibble(label          = label,
         n              = n,
         n_positive     = k,
         prop_positive  = k / n,
         cp_lower       = cp[1],
         cp_upper       = cp[2],
         sign_test_p    = bin$p.value,
         wilcoxon_p     = wil$p.value,
         mean_value     = mean(z),
         median_value   = median(z))
}

# ---------------------------------------------------------------------
# 3. Load empirical binned variogram and compute residuals per subject
# ---------------------------------------------------------------------
cat("Loading binned variogram and computing Gaussian residuals...\n")

emp <- read_excel(BINNED_XLSX) %>%
  arrange(Subject, Distance) %>%
  group_by(Subject) %>%
  mutate(bin = row_number()) %>%
  ungroup()

stopifnot(all(emp %>% count(Subject) %>% pull(n) == 17))

compute_residuals <- function(d_sub) {
  obj <- fit_gaussian(d_sub$Distance, d_sub$MeanGamma)
  fitted_vals <- if (is.null(obj)) rep(NA_real_, nrow(d_sub)) else fitted(obj)
  tibble(bin      = d_sub$bin,
         distance = d_sub$Distance,
         observed = d_sub$MeanGamma,
         fitted   = fitted_vals,
         residual = d_sub$MeanGamma - fitted_vals)
}

emp_resid <- emp %>%
  group_by(Subject) %>%
  group_modify(~ compute_residuals(.x)) %>%
  ungroup()

failed <- emp_resid %>%
  group_by(Subject) %>% summarise(any_na = anyNA(residual)) %>%
  filter(any_na) %>% pull(Subject)
if (length(failed))
  warning("Gaussian fit failed for subjects: ", paste(failed, collapse = ", "))

# ---------------------------------------------------------------------
# 4. Zone-level tests
# ---------------------------------------------------------------------
cat("\n=== Corner vs Plateau zone tests ===\n")

pathA_zone_means <- emp_resid %>%
  mutate(zone = case_when(
    bin %in% CORNER_BINS  ~ "corner",
    bin %in% PLATEAU_BINS ~ "plateau",
    TRUE                  ~ NA_character_)) %>%
  filter(!is.na(zone)) %>%
  group_by(Subject, zone) %>%
  summarise(zone_mean = mean(residual, na.rm = TRUE),
            n_bins    = n(),
            .groups   = "drop") %>%
  mutate(sign_pos = zone_mean > 0)

write.csv(pathA_zone_means,
          "priority3_pathA_zone_stats.csv", row.names = FALSE)

pathA_tests <- bind_rows(
  run_sign_tests(pathA_zone_means$zone_mean[pathA_zone_means$zone == "corner"],
                 "corner (bins 2-4)"),
  run_sign_tests(pathA_zone_means$zone_mean[pathA_zone_means$zone == "plateau"],
                 "plateau (bins 5-10)")
) %>% rename(zone = label)

write.csv(pathA_tests,
          "priority3_pathA_test_results.csv", row.names = FALSE)

print(pathA_tests)

cat("\nHeadline numbers:\n")
for (i in seq_len(nrow(pathA_tests))) {
  r <- pathA_tests[i, ]
  cat(sprintf(
    "  %-22s : %d/%d positive  (prop = %.3f, 95%% CP CI [%.3f, %.3f])\n",
    r$zone, r$n_positive, r$n, r$prop_positive,
    r$cp_lower, r$cp_upper))
  cat(sprintf("                          sign p = %.2e   Wilcoxon p = %.2e\n",
              r$sign_test_p, r$wilcoxon_p))
}

# ---------------------------------------------------------------------
# 5. Cohort mean residual table (used by both residual figures)
# ---------------------------------------------------------------------
cohort_resid <- emp_resid %>%
  filter(is.finite(residual)) %>%
  group_by(bin, distance) %>%
  summarise(mean_resid = mean(residual),
            se_resid   = sd(residual) / sqrt(n()),
            .groups    = "drop") %>%
  arrange(bin)

# Helper: produce x-range bounds for one set of bins, with half-bin pad
zone_bounds <- function(bins, pad = 0.45) {
  ds <- cohort_resid$distance[cohort_resid$bin %in% bins]
  c(min(ds) - pad, max(ds) + pad)
}

bin1_x    <- zone_bounds(EXCLUDED_BIN1)
corner_x  <- zone_bounds(CORNER_BINS)
plateau_x <- zone_bounds(PLATEAU_BINS)
recov_x   <- zone_bounds(RECOVERY_BINS)

# ---------------------------------------------------------------------
# 6. Figure: residual structure with zones and SE band
#    
# ---------------------------------------------------------------------
p_structure <- ggplot(cohort_resid, aes(distance, mean_resid)) +
  annotate("rect", xmin = corner_x[1],  xmax = corner_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$corner,  alpha = 0.20) +
  annotate("rect", xmin = plateau_x[1], xmax = plateau_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$plateau, alpha = 0.16) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = mean_resid - se_resid,
                  ymax = mean_resid + se_resid),
              fill = ZONE_COLS$curve, alpha = 0.25) +
  geom_line(colour  = ZONE_COLS$curve, linewidth = 0.9) +
  geom_point(colour = ZONE_COLS$curve, size = 2.2) +
  annotate("text",
           x = mean(corner_x), y = Inf, vjust = 1.6,
           label = "Corner zone\n(bins 2-4)\nunder-prediction",
           colour = TEXT_COLS$corner, size = 3.2, lineheight = 0.95) +
  annotate("text",
           x = mean(plateau_x), y = Inf, vjust = 1.6,
           label = "Plateau zone\n(bins 5-10)\nover-prediction",
           colour = TEXT_COLS$plateau, size = 3.2, lineheight = 0.95) +
  labs(title    = "Cohort mean residual from Gaussian fit, with zones",
       subtitle = sprintf("Ribbon = mean +/- 1 SE across %d subjects; bins 1 and 11-17 excluded with stated reasons",
                          length(unique(emp$Subject))),
       x = "Distance (mm)", y = "Mean residual (observed - fitted)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("figure_residual_structure.png", p_structure,
       width = 9.0, height = 5.5, dpi = PLOT_DPI)
ggsave("figure_residual_structure.pdf", p_structure,
       width = 9.0, height = 5.5, device = cairo_pdf)

# ---------------------------------------------------------------------
# 7. Figure: four-zone presentation
# ---------------------------------------------------------------------
zone_labels <- tibble(
  x_mid = c(mean(bin1_x), mean(corner_x), mean(plateau_x), mean(recov_x)),
  y     = 0.30,
  text  = c("Bin 1\nshort-lag\nover-prediction",
            "Corner zone\n(bins 2-4)\nunder-prediction",
            "Plateau zone\n(bins 5-10)\nover-prediction",
            "Long-distance recovery\n(bins 11-17)\nresidual approaches zero"),
  col   = c(TEXT_COLS$bin1, TEXT_COLS$corner,
            TEXT_COLS$plateau, TEXT_COLS$recovery),
  size  = c(3.0, 3.3, 3.3, 3.3)
)

p_modes <- ggplot(cohort_resid, aes(distance, mean_resid)) +
  annotate("rect", xmin = bin1_x[1],    xmax = bin1_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$bin1,     alpha = 0.20) +
  annotate("rect", xmin = corner_x[1],  xmax = corner_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$corner,   alpha = 0.20) +
  annotate("rect", xmin = plateau_x[1], xmax = plateau_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$plateau,  alpha = 0.16) +
  annotate("rect", xmin = recov_x[1],   xmax = recov_x[2],
           ymin = -Inf, ymax = Inf,
           fill = ZONE_COLS$recovery, alpha = 0.13) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_line(colour  = ZONE_COLS$curve, linewidth = 1.0) +
  geom_point(colour = ZONE_COLS$curve, size = 2.4) +
  geom_text(data = zone_labels,
            aes(x = x_mid, y = y, label = text),
            colour = zone_labels$col, size = zone_labels$size,
            lineheight = 0.95, inherit.aes = FALSE) +
  labs(title    = "Cohort mean residual from Gaussian fit",
       subtitle = "Four coherent regions across the binned variogram",
       x = "Distance (mm)",
       y = "Mean residual (observed - fitted)") +
  scale_x_continuous(breaks = seq(5, 20, by = 5)) +
  coord_cartesian(ylim = c(-0.22, 0.36)) +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey30", size = 10),
        panel.grid.minor = element_blank())

ggsave("figure_residual_failure_modes.png", p_modes,
       width = 10.0, height = 5.4, dpi = PLOT_DPI)
ggsave("figure_residual_failure_modes.pdf", p_modes,
       width = 10.0, height = 5.4, device = cairo_pdf)

# ---------------------------------------------------------------------
# 8. Done
# ---------------------------------------------------------------------
cat("\nOutputs saved:\n",
    "  CSV:\n",
    "    priority3_pathA_zone_stats.csv\n",
    "    priority3_pathA_test_results.csv\n",
    "  Figures (PNG + PDF):\n",
    "    figure_residual_structure\n",
    "    figure_residual_failure_modes\n", sep = "")