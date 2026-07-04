# =====================================================================
# Directional / anisotropy check
# ---------------------------------------------------------------------
# Inputs:  Per-subject directional CSVs from variogram_directional_final.m,
#          one row per (shift, run), columns:
#              Subject, run, dx, dy, dz, distance_mm, gamma_run
#
# Outputs:
#   * Per-subject directional variograms along x, y, z (axis-aligned
#     lags) plus the omnidirectional curve.
#   * Per-distance max-axis ratio of mean gamma, as a quantitative
#     anisotropy diagnostic.
#   * Plots: x/y/z/omni overlay and the ratio-vs-distance plot.
#
# Conventions:
#   * Your MATLAB shifts use dx, dy, dz in {0, 1, 2, 3, 4}; the
#     variogram is symmetric in lag direction so only non-negative
#     shifts were stored. Axis-aligned lag = exactly one of
#     (dx, dy, dz) is non-zero.
#   * Voxel size in your pipeline is read per-subject from the SPM
#     header. For the directional binning we work in millimetres
#     directly (distance_mm is already in the CSV), so this script
#     is voxel-size-agnostic.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(readr)

# ---- CONFIG ----------------------------------------------------------
# Point this at the directory written by variogram_directional_final.m
DIR_CSV   <- "directional_csvs"

# Must match SUBJECTS_FOR_DIRECTIONAL in the MATLAB driver
SUBJECTS  <- c(5, 21, 28)

# Binning for the directional curves. 1 mm bins are fine here because
# we are averaging fewer voxel pairs per bin than the omnidirectional
# fit, so a slightly wider bin keeps SE reasonable.
BIN_WIDTH <- 1.0     # mm
MAX_DIST  <- 22.0    # mm  (max possible is 4*3*sqrt(3) ~ 20.8)
# ---------------------------------------------------------------------

# ---- Load and combine -----------------------------------------------
load_subject <- function(sid) {
  f <- file.path(DIR_CSV, sprintf("sub-%02d_directional_variogram.csv", sid))
  if (!file.exists(f)) { warning("Missing CSV: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE)
}
dir_data <- bind_rows(lapply(SUBJECTS, load_subject))
stopifnot(nrow(dir_data) > 0)

cat(sprintf("Loaded %d total (shift x run) rows across %d subjects.\n",
            nrow(dir_data), length(unique(dir_data$Subject))))

# ---- Average across runs per shift (matches the original pipeline) ---
# The original variogram analysis averages gamma across runs before
# anything else (mean_gamma in variogram_directional_final.m). We do the
# same here so the directional curves are on the same footing as the
# omnidirectional one used elsewhere in the paper.
shift_data <- dir_data %>%
  group_by(Subject, dx, dy, dz, distance_mm) %>%
  summarise(gamma_mean = mean(gamma_run),
            n_runs     = n(),
            .groups = "drop")

# ---- Axis labelling --------------------------------------------------
# Axis-aligned = exactly one of (dx, dy, dz) non-zero
shift_data <- shift_data %>%
  mutate(n_nonzero = (dx != 0) + (dy != 0) + (dz != 0),
         axis = case_when(
           n_nonzero == 1 & dx != 0 ~ "x",
           n_nonzero == 1 & dy != 0 ~ "y",
           n_nonzero == 1 & dz != 0 ~ "z",
           TRUE ~ "other"))

cat("\nAxis-aligned shift counts per subject:\n")
print(shift_data %>% filter(axis %in% c("x","y","z")) %>%
        count(Subject, axis))

# ---- Bin per (Subject, axis, bin) -----------------------------------
bin_edges   <- seq(0, MAX_DIST, by = BIN_WIDTH)
bin_centres <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2

shift_data <- shift_data %>%
  mutate(bin = cut(distance_mm, breaks = bin_edges,
                   include.lowest = TRUE, labels = FALSE),
         d_centre = bin_centres[bin])

# Per-axis directional variogram
dir_axes <- shift_data %>%
  filter(axis %in% c("x","y","z")) %>%
  group_by(Subject, axis, d_centre) %>%
  summarise(gamma_mean = mean(gamma_mean),
            n_pairs    = n(),
            .groups = "drop") %>%
  rename(Direction = axis)

# Omnidirectional (all shifts pooled) -- this should match the original
# pipeline's variogram up to the binning granularity.
dir_omni <- shift_data %>%
  group_by(Subject, d_centre) %>%
  summarise(gamma_mean = mean(gamma_mean),
            n_pairs    = n(),
            .groups = "drop") %>%
  mutate(Direction = "omni")

dir_all <- bind_rows(dir_axes, dir_omni) %>%
  mutate(Direction = factor(Direction, levels = c("x","y","z","omni")))

write.csv(dir_all, "directional_variograms_binned.csv", row.names = FALSE)

# ---- Anisotropy diagnostic: max/min axis ratio ----------------------
ratio_df <- dir_axes %>%
  select(Subject, d_centre, Direction, gamma_mean) %>%
  pivot_wider(names_from = Direction, values_from = gamma_mean) %>%
  filter(!is.na(x) & !is.na(y) & !is.na(z)) %>%
  mutate(max_axis = pmax(x, y, z),
         min_axis = pmin(x, y, z),
         ratio    = max_axis / min_axis,
         range_pct = (max_axis - min_axis) / ((x + y + z)/3) * 100)

write.csv(ratio_df, "directional_anisotropy_ratios.csv", row.names = FALSE)

anisotropy_summary <- ratio_df %>%
  group_by(Subject) %>%
  summarise(median_ratio    = median(ratio),
            p90_ratio       = quantile(ratio, 0.90),
            max_ratio       = max(ratio),
            median_range_pct = median(range_pct),
            n_bins          = n(),
            .groups = "drop")

cat("\nPer-subject anisotropy summary:\n")
cat("(ratio = max(gamma_x, gamma_y, gamma_z) / min(...);  1 = perfect isotropy)\n")
print(anisotropy_summary)
write.csv(anisotropy_summary, "directional_anisotropy_summary.csv",
          row.names = FALSE)

# ---- Plot 1: x/y/z/omni overlay per subject -------------------------
p_dir <- ggplot(dir_all %>% filter(!is.na(d_centre)),
                aes(d_centre, gamma_mean,
                    colour = Direction, linetype = Direction)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ Subject, ncol = length(SUBJECTS),
             labeller = labeller(Subject = function(x) paste("Subject", x))) +
  scale_colour_manual(values = c(x = "#C44536", y = "#2A8C7D",
                                 z = "#3D6FB5", omni = "grey25")) +
  scale_linetype_manual(values = c(x = "solid", y = "solid",
                                   z = "solid", omni = "dashed")) +
  labs(title = "Directional vs omnidirectional variograms",
       subtitle = "x, y, z = axis-aligned lags only; omni = all directions pooled",
       x = "Distance (mm)", y = expression(gamma)) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("directional_variograms.png", p_dir,
       width = 4.2 * length(SUBJECTS), height = 4.8, dpi = 200)

# ---- Plot 2: anisotropy ratio vs distance ---------------------------
p_ratio <- ggplot(ratio_df, aes(d_centre, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 1.3, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ Subject, ncol = length(SUBJECTS),
             labeller = labeller(Subject = function(x) paste("Subject", x))) +
  scale_y_continuous(limits = c(1, NA)) +
  labs(title = "Anisotropy ratio across distance bins",
       subtitle = "Dashed = 1 (perfect isotropy);  dotted = 1.3 (moderate)",
       x = "Distance (mm)", y = "max/min axis ratio") +
  theme_bw()
ggsave("directional_anisotropy_ratio.png", p_ratio,
       width = 4.2 * length(SUBJECTS), height = 4.5, dpi = 200)

cat("\nPriority 2 outputs saved:\n",
    " - directional_variograms_binned.csv\n",
    " - directional_anisotropy_ratios.csv\n",
    " - directional_anisotropy_summary.csv\n",
    " - directional_variograms.png\n",
    " - directional_anisotropy_ratio.png\n", sep = "")
