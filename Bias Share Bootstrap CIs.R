# =============================================================================
#  Bootstrap CIs for bias²/variance shares
# -----------------------------------------------------------------------------
# The bootstrap analysis produces a full factorial:
#   4 noise models × 7 CV levels × 2 true-curve families × 21 subjects × 4 methods
# This script analyses a 3-CV subset of that grid (c = 0.05, 0.10, 0.15;
# see CVS below), which spans the empirically relevant low-noise regime.
#
# For each cell of that factorial, how stable is the
# mean bias² share across subjects?  Answer via subject-level bootstrap
# (resample the 21 subjects B = 1000 times, report percentile 95% CIs).
#
# Prerequisites: syn_metrics data frame already in workspace (produced by
# the metrics section of the bootstrap analysis,
# priority1_v3_synthetic_bootstrap.R), OR loaded from disk:
#   syn_metrics <- readRDS("synthetic_bootstrap_v3_metrics.rds")
#   # OR
#   syn_metrics <- read.csv("synthetic_bootstrap_v3_metrics.csv")
#
# The empirical noise level is CV ≈ 0.011–0.025 (StdGamma / (MeanGamma * sqrt(R))
# where R = 5 runs).  This sits well below the lowest CV tested (c = 0.05),
# so the c = 0.05 slice is the primary manuscript result.
# =============================================================================

library(dplyr); library(tidyr); library(ggplot2)

set.seed(20260301)
B_SUBJ      <- 1000   # subject-level bootstrap replicates
N_SUBJ      <- 21
METHODS     <- c("LocalCubic", "PTR", "NatCubicSpline", "Gaussian")
NOISE_MODELS<- c("Gaussian", "ChiSquared", "Gamma", "Lognormal")
FAMILIES    <- c("DoubleExp", "Matern")
CVS         <- c(0.05, 0.10, 0.15)

METHOD_LABELS <- c(LocalCubic     = "LP",
                   PTR            = "PTR",
                   NatCubicSpline = "NCS",
                   Gaussian       = "Gaussian")
COLS_M <- c(LocalCubic     = "#1D9E75",
            PTR            = "#BA7517",
            NatCubicSpline = "#378ADD",
            Gaussian       = "#D85A30")

# Load from disk if syn_metrics is not already in workspace
if (!exists("syn_metrics")) {
  if (file.exists("synthetic_bootstrap_v3_metrics.rds")) {
    syn_metrics <- readRDS("synthetic_bootstrap_v3_metrics.rds")
    cat("Loaded syn_metrics from RDS.\n")
  } else {
    syn_metrics <- read.csv("synthetic_bootstrap_v3_metrics.csv",
                            stringsAsFactors = FALSE)
    cat("Loaded syn_metrics from CSV.\n")
  }
}

# ------------------------------------------------------------------
# Step 1: Per-subject bias share
# ------------------------------------------------------------------
# bias_share_s = IBias2 / IMSE for each (Method, NoiseModel, CV, Family, Subject)
# Exclude rows where IMSE is NA or effectively zero (degenerate fits)

syn_bs <- syn_metrics %>%
  filter(!is.na(IMSE), IMSE > 0, !is.na(IBias2)) %>%
  mutate(bias_share = IBias2 / IMSE,
         Method_label = METHOD_LABELS[Method])

# Sanity check
cat(sprintf("Rows after filtering: %d / %d\n", nrow(syn_bs), nrow(syn_metrics)))

# ------------------------------------------------------------------
# Step 2: Subject-level bootstrap per (NoiseModel, CV, Family, Method)
# ------------------------------------------------------------------
# For each cell: resample subjects B times, compute mean bias share each time,
# return point estimate + percentile 95% CI.

bootstrap_cell <- function(shares, B = B_SUBJ) {
  n    <- length(shares)
  boot <- replicate(B, mean(sample(shares, size = n, replace = TRUE)))
  c(mean_bs  = mean(shares),        # observed mean
    lo95     = quantile(boot, 0.025, names = FALSE),
    hi95     = quantile(boot, 0.975, names = FALSE))
}

cat("\nRunning subject-level bootstrap",
    sprintf("(B = %d) across %d cells...\n", B_SUBJ,
            length(NOISE_MODELS) * length(CVS) * length(FAMILIES) * length(METHODS)))

ci_rows <- list()
for (nm in NOISE_MODELS) {
  for (cv in CVS) {
    for (fam in FAMILIES) {
      for (m in METHODS) {
        shares <- syn_bs %>%
          filter(NoiseModel == nm, CV == cv, Family == fam, Method == m) %>%
          pull(bias_share)
        if (length(shares) < 2) next
        res <- bootstrap_cell(shares)
        ci_rows[[length(ci_rows) + 1]] <- data.frame(
          NoiseModel   = nm,
          CV           = cv,
          Family       = fam,
          Method       = m,
          Method_label = METHOD_LABELS[m],
          n_subj       = length(shares),
          mean_pct     = round(res["mean_bs"]  * 100, 1),
          lo95_pct     = round(res["lo95"]     * 100, 1),
          hi95_pct     = round(res["hi95"]     * 100, 1),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

ci_df <- do.call(rbind, ci_rows)
ci_df$Method_label <- factor(ci_df$Method_label,
                             levels = c("LP", "PTR", "NCS", "Gaussian"))
cat("Bootstrap complete.\n\n")

# ------------------------------------------------------------------
# Step 3: Primary manuscript table — c = 0.05, averaged across noise models
# ------------------------------------------------------------------
# The empirical CV is ~0.011–0.025, so c = 0.05 is the closest grid point.
# Results are nearly identical across noise models (matched variance design),
# so we average across them for a compact table.

tbl_primary <- ci_df %>%
  filter(CV == 0.05) %>%
  group_by(Family, Method, Method_label) %>%
  summarise(
    mean_pct = round(mean(mean_pct), 1),
    lo95_pct = round(mean(lo95_pct), 1),
    hi95_pct = round(mean(hi95_pct), 1),
    .groups  = "drop"
  ) %>%
  mutate(Report = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                          mean_pct, lo95_pct, hi95_pct)) %>%
  arrange(Family, mean_pct)

cat("=== PRIMARY TABLE: bias² share at c = 0.05",
    "(averaged across 4 noise models, B = 1000) ===\n")
cat(sprintf("%-12s %-10s %30s\n", "Family", "Method", "Bias² share [95% CI]"))
cat(strrep("-", 54), "\n")
for (i in seq_len(nrow(tbl_primary))) {
  r <- tbl_primary[i, ]
  cat(sprintf("%-12s %-10s %30s\n", r$Family, r$Method_label, r$Report))
}
cat("\n")

# ------------------------------------------------------------------
# Step 4: How bias share CI changes across CV levels
#         (collapsed across noise models — they are nearly identical)
# ------------------------------------------------------------------

tbl_cv <- ci_df %>%
  group_by(Family, CV, Method, Method_label) %>%
  summarise(
    mean_pct = round(mean(mean_pct), 1),
    lo95_pct = round(mean(lo95_pct), 1),
    hi95_pct = round(mean(hi95_pct), 1),
    .groups  = "drop"
  )

cat("=== BIAS SHARE (%) BY CV — averaged across noise models ===\n")
for (fam in FAMILIES) {
  cat(sprintf("\nFamily: %s\n", fam))
  cat(sprintf("%-10s", "Method"))
  for (cv in CVS) cat(sprintf(" %9s", paste0("c=", cv)))
  cat("\n", strrep("-", 10 + 10*length(CVS)), "\n", sep="")
  for (m in c("LP","PTR","NCS","Gaussian")) {
    cat(sprintf("%-10s", m))
    for (cv in CVS) {
      r <- tbl_cv %>% filter(Family == fam, CV == cv, Method_label == m)
      if (nrow(r) == 1)
        cat(sprintf(" %9s", sprintf("%.0f[%.0f,%.0f]", r$mean_pct, r$lo95_pct, r$hi95_pct)))
      else cat(sprintf(" %9s", "--"))
    }
    cat("\n")
  }
}

# ------------------------------------------------------------------
# Step 5: Figures
# ------------------------------------------------------------------

# F1: Bias share + 95% CI vs CV, faceted Family × Method_label
#     One line per noise model, averaged for the CI ribbon
p_cv <- ggplot(tbl_cv,
               aes(x = CV, y = mean_pct,
                   colour = Method_label,
                   fill   = Method_label)) +
  geom_ribbon(aes(ymin = lo95_pct, ymax = hi95_pct), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 0.025, linetype = "dashed",
             colour = "gray40", linewidth = 0.7) +
  annotate("text", x = 0.027, y = 95, label = "empirical CV ≈ 0.025",
           hjust = 0, size = 3.0, colour = "gray40") +
  scale_colour_manual(values = COLS_M) +
  scale_fill_manual(values   = COLS_M) +
  scale_x_continuous(breaks  = CVS) +
  scale_y_continuous(labels  = function(x) paste0(x, "%"),
                     limits  = c(0, 100)) +
  facet_wrap(~ Family, nrow = 1) +
  labs(title    = "Priority 4: Bias² share of IMSE with 95% CIs — v3 synthetic bootstrap",
       subtitle  = "Averaged across 4 noise models | Ribbons = 95% CI (B = 1000 subject resamples) | Dashed = empirical CV",
       x        = "Coefficient of variation c",
       y        = "Bias² share of IMSE (%)",
       colour   = "Method", fill = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))
print(p_cv)
ggsave("p4_v3_bias_share_vs_cv.png", p_cv,
       width = 11, height = 5, dpi = 200)

# F2: Bar chart at c = 0.05 with CI error bars — one panel per Family,
#     bars coloured by Method, separated by Family
p_primary <- ggplot(tbl_primary,
                    aes(x    = reorder(Method_label, mean_pct),
                        y    = mean_pct,
                        fill = Method)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo95_pct, ymax = hi95_pct),
                width = 0.25, linewidth = 0.9, colour = "gray25") +
  geom_text(aes(label = sprintf("%.1f%%\n[%.1f, %.1f]",
                                mean_pct, lo95_pct, hi95_pct)),
            vjust = -0.3, size = 3.0, colour = "gray20", lineheight = 0.9) +
  scale_fill_manual(values = COLS_M, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 105)) +
  facet_wrap(~ Family, nrow = 1) +
  labs(title    = "Priority 4: Bias² share at c = 0.05 — v3 synthetic bootstrap",
       subtitle  = "Averaged across 4 noise models | Error bars = 95% CI (B = 1000 subject resamples)",
       x        = NULL,
       y        = "Bias² share of IMSE (%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_primary)
ggsave("p4_v3_bias_share_c05.png", p_primary,
       width = 10, height = 5, dpi = 200)

# F3: Per-noise-model breakdown at c = 0.05 — shows robustness
p_noise <- ci_df %>%
  filter(CV == 0.05) %>%
  ggplot(aes(x    = reorder(Method_label, mean_pct),
             y    = mean_pct,
             fill = Method)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo95_pct, ymax = hi95_pct),
                width = 0.25, linewidth = 0.8, colour = "gray30") +
  scale_fill_manual(values = COLS_M, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 105)) +
  facet_grid(Family ~ NoiseModel) +
  labs(title    = "Priority 4: Bias² share at c = 0.05 — breakdown by noise model",
       subtitle  = "Error bars = 95% CI (B = 1000) | Results are stable across all 4 noise distributions",
       x        = NULL,
       y        = "Bias² share of IMSE (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title  = element_text(face = "bold"))
print(p_noise)
ggsave("p4_v3_bias_share_c05_by_noise.png", p_noise,
       width = 12, height = 6, dpi = 200)

# ------------------------------------------------------------------
# Step 6: Save outputs
# ------------------------------------------------------------------

write.csv(ci_df,      "p4_v3_bias_share_CIs_full.csv",    row.names = FALSE)
write.csv(tbl_primary,"p4_v3_bias_share_CIs_primary.csv", row.names = FALSE)
write.csv(tbl_cv,     "p4_v3_bias_share_CIs_by_cv.csv",   row.names = FALSE)

cat("=== MANUSCRIPT NUMBERS (c = 0.05, averaged across noise models) ===\n")
for (fam in FAMILIES) {
  cat(sprintf("\nFamily: %s\n", fam))
  tbl_primary %>%
    filter(Family == fam) %>%
    arrange(mean_pct) %>%
    { for (i in seq_len(nrow(.))) {
      r <- .[i, ]
      cat(sprintf("  %s: %s\n", r$Method_label, r$Report))
    }
    }
}
cat("\nAll Priority 4 (v3) outputs saved.\n")