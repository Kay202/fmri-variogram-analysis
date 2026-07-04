# =====================================================================
# Full factorial synthetic-truth bootstrap
# ---------------------------------------------------------------------
# Design:
#   * CV grid expanded to {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35}
#   * B = 1000 replicates per cell
#   * Per-CELL checkpointing, where a cell is one
#     (noise model, CV, family, subject) combination, so a crash loses
#     at most one cell of work.
#   * Automatic resume: if synthetic_bootstrap_v3_partial.rds exists,
#     the script picks up where it left off.
#   * Wall-clock ETA printed after each cell.
#
# Total fits: 4 noise models * 7 CVs * 2 families * 21 subjects
#             * 1000 reps * 4 fitting methods = 4 704 000 fits.
#             Expect ~50 hours.
#
# Requirements (loaded before sourcing this file):
#   - locpoly3(), locpoly3_wrap(), search_fn(), kvals(), gsjs(), gcv_fn()
#   - ptr(), ptr_wrap()
#   - ncs_gcv_edf(), fit_ncs()
#   - fit_gaussian()
#   - ggplot2, dplyr, tidyr, minpack.lm
# =====================================================================

set.seed(2026)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------
# 1. Distance grid
# ---------------------------------------------------------------------
DIST_MIN_MM <- 3.4
DIST_MAX_MM <- 19.5
N_BINS      <- 17
d_mm        <- seq(DIST_MIN_MM, DIST_MAX_MM, length.out = N_BINS)
t_norm_grid <- (d_mm - DIST_MIN_MM) / (DIST_MAX_MM - DIST_MIN_MM)

# ---------------------------------------------------------------------
# 2. A priori priors
# ---------------------------------------------------------------------
prior_double_exp <- function() {
  list(sigma2 = runif(1, 0.85, 1.05),
       a      = runif(1, 0.10, 0.20),
       A      = runif(1, 0.05, 0.15),
       omega  = runif(1, 4*pi, 8*pi),
       b      = runif(1, 0.20, 0.45))
}
prior_matern <- function() {
  list(sigma2 = runif(1, 0.85, 1.05),
       rho    = runif(1, 0.10, 0.22),
       A      = runif(1, 0.05, 0.15),
       omega  = runif(1, 4*pi, 8*pi),
       b      = runif(1, 0.20, 0.45))
}
prior_nugget <- function() runif(1, 0.00, 0.10)

# ---------------------------------------------------------------------
# 3. True-curve generators
# ---------------------------------------------------------------------
true_double_exp <- function(t, p) {
  rise <- p$sigma2 * (1 - exp(-t / p$a))
  hole <- p$A * cos(p$omega * t) * exp(-t / p$b)
  rise + hole
}
matern_2p5 <- function(t, rho) {
  u <- sqrt(5) * t / rho
  (1 + u + (u^2) / 3) * exp(-u)
}
true_matern <- function(t, p) {
  rise <- p$sigma2 * (1 - matern_2p5(t, p$rho))
  hole <- p$A * cos(p$omega * t) * exp(-t / p$b)
  rise + hole
}

# ---------------------------------------------------------------------
# 4. Generate 21 true curves per family
# ---------------------------------------------------------------------
N_SUBJ       <- 21
CVS          <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)
FAMILIES     <- c("DoubleExp", "Matern")
NOISE_MODELS <- c("Gaussian", "ChiSquared", "Gamma", "Lognormal")
B            <- 1000
METHODS      <- c("LocalCubic", "PTR", "NatCubicSpline", "Gaussian")

CURVES_FILE  <- "synthetic_true_curves_v3.rds"
PARTIAL_FILE <- "synthetic_bootstrap_v3_partial.rds"
FINAL_FILE   <- "synthetic_bootstrap_v3_FINAL.rds"

if (file.exists(CURVES_FILE)) {
  cat("Loading existing true curves from", CURVES_FILE, "...\n")
  saved   <- readRDS(CURVES_FILE)
  true_curves_syn <- saved$curves
  nuggets         <- saved$nuggets
} else {
  cat("Generating new true curves...\n")
  true_curves_syn <- list(DoubleExp = list(), Matern = list())
  nuggets         <- numeric(N_SUBJ)
  for (s in seq_len(N_SUBJ)) {
    nug <- prior_nugget(); nuggets[s] <- nug
    p_de <- prior_double_exp()
    p_mt <- prior_matern()
    true_curves_syn$DoubleExp[[s]] <- list(
      params = p_de, nugget = nug,
      y_true = pmax(true_double_exp(t_norm_grid, p_de) + nug, 1e-6))
    true_curves_syn$Matern[[s]] <- list(
      params = p_mt, nugget = nug,
      y_true = pmax(true_matern(t_norm_grid, p_mt) + nug, 1e-6))
  }
  saveRDS(list(curves = true_curves_syn, nuggets = nuggets), CURVES_FILE)
}
cat("21 synthetic true curves ready for each family.\n")

# ---------------------------------------------------------------------
# 5. Noise generators (matched variance: E=gamma_k, Var=(c*gamma_k)^2)
# ---------------------------------------------------------------------
noise_gaussian  <- function(gamma, cv) {
  gamma + rnorm(length(gamma), 0, cv * gamma)
}
noise_chisq     <- function(gamma, cv) {
  nu <- 2 / cv^2
  gamma * rchisq(length(gamma), df = nu) / nu
}
noise_gamma     <- function(gamma, cv) {
  shape <- 1 / cv^2; scale <- cv^2 * gamma
  rgamma(length(gamma), shape = shape, scale = scale)
}
noise_lognormal <- function(gamma, cv) {
  s2 <- log(1 + cv^2); mu <- log(gamma) - s2 / 2
  rlnorm(length(gamma), meanlog = mu, sdlog = sqrt(s2))
}
draw_noise <- function(model, gamma, cv) {
  switch(model,
         Gaussian   = noise_gaussian(gamma, cv),
         ChiSquared = noise_chisq(gamma, cv),
         Gamma      = noise_gamma(gamma, cv),
         Lognormal  = noise_lognormal(gamma, cv),
         stop("unknown noise model: ", model))
}

# ---------------------------------------------------------------------
# 6. Pre-flight check (10k draws per model x CV)
#    Skipped if preflight summary already exists.
# ---------------------------------------------------------------------
PREFLIGHT_FILE <- "preflight_noise_summary_v3.csv"
if (!file.exists(PREFLIGHT_FILE)) {
  cat("\n--- Pre-flight noise calibration check ---\n")
  ref_gamma <- pmax(true_double_exp(t_norm_grid, list(
    sigma2 = 0.95, a = 0.15, A = 0.10, omega = 6*pi, b = 0.30)) + 0.05, 1e-6)
  N_PREFLIGHT <- 10000
  preflight <- data.frame()
  for (model in NOISE_MODELS) {
    for (cv in CVS) {
      mat <- replicate(N_PREFLIGHT, draw_noise(model, ref_gamma, cv))
      bin_sds  <- apply(mat, 1, sd)
      bin_skew <- apply(mat, 1, function(x) mean(((x-mean(x))/sd(x))^3))
      tgt_sd   <- cv * ref_gamma
      preflight <- rbind(preflight, data.frame(
        Model = model, CV = cv,
        mean_sd_ratio = mean(bin_sds / tgt_sd),
        mean_skew     = mean(bin_skew),
        frac_negative = mean(mat < 0)
      ))
    }
  }
  write.csv(preflight, PREFLIGHT_FILE, row.names = FALSE)
  cat("Pre-flight summary written.\n")
  print(preflight)
} else {
  cat("Pre-flight summary already exists; skipping.\n")
}

# ---------------------------------------------------------------------
# 7. Safe fitters
# ---------------------------------------------------------------------
mk_bad <- function(res, n_expected) {
  tag <- if (length(res) == 1 && is.na(res)) {
    paste0("error:", (attr(res, "msg") %||% "NA"))
  } else sprintf("bad_len_%d", length(res))
  list(fit = rep(NA_real_, n_expected), status = tag)
}
safe_lp <- function(t, y, n_expected) {
  res <- tryCatch({
    bw_mesh <- seq(0.05, 0.5, by = 0.05)
    lp_crit <- search_fn(t, y, bw_mesh, locpoly3_wrap)
    bw      <- bw_mesh[which.min(lp_crit[, 2])]
    locpoly3(t, y, bw, t)[[1]][, 1]
  }, error = function(e) structure(NA_real_, msg = conditionMessage(e)))
  if (length(res) == n_expected && all(is.finite(res))) list(fit = as.numeric(res), status = "ok")
  else mk_bad(res, n_expected)
}
safe_ptr <- function(t, y, n_expected) {
  res <- tryCatch({
    ptr_mesh <- 0:8
    ptr_crit <- search_fn(t, y, ptr_mesh, ptr)
    k        <- ptr_mesh[which.min(ptr_crit[, 2])]
    as.vector(ptr(t, y, k)[[1]])
  }, error = function(e) structure(NA_real_, msg = conditionMessage(e)))
  if (length(res) == n_expected && all(is.finite(res))) list(fit = as.numeric(res), status = "ok")
  else mk_bad(res, n_expected)
}
safe_ncs <- function(t, y, n_expected) {
  res <- tryCatch(fit_ncs(t, y)$fitted,
                  error = function(e) structure(NA_real_, msg = conditionMessage(e)))
  if (length(res) == n_expected && all(is.finite(res))) list(fit = as.numeric(res), status = "ok")
  else mk_bad(res, n_expected)
}
safe_gauss <- function(d_in, y, n_expected) {
  res <- tryCatch({
    obj <- fit_gaussian(d_in, y)
    if (is.null(obj) || !inherits(obj, "nls")) {
      structure(NA_real_, msg = "fit_gaussian returned non-nls")
    } else fitted(obj)
  }, error   = function(e) structure(NA_real_, msg = conditionMessage(e)),
  warning = function(w) structure(NA_real_, msg = conditionMessage(w)))
  if (length(res) == n_expected && all(is.finite(res))) list(fit = as.numeric(res), status = "ok")
  else mk_bad(res, n_expected)
}

# ---------------------------------------------------------------------
# 8. Cell key encoding and resume logic
# ---------------------------------------------------------------------
# Each cell is identified by (NoiseModel, CV, Family, Subject). To make
# checkpointing fast we store results in a flat list keyed by a string.
cell_key <- function(model, cv, fam, s) {
  sprintf("%s|cv_%.2f|%s|s_%02d", model, cv, fam, s)
}

run_one_cell <- function(family, noise_model, cv, s) {
  tc <- true_curves_syn[[family]][[s]]
  y_true <- tc$y_true
  fits_array <- array(NA_real_, dim = c(B, N_BINS, length(METHODS)),
                      dimnames = list(NULL, NULL, METHODS))
  st <- matrix("", nrow = B, ncol = length(METHODS),
               dimnames = list(NULL, METHODS))
  for (r in seq_len(B)) {
    y_sim <- draw_noise(noise_model, y_true, cv)
    out_lp    <- safe_lp(t_norm_grid, y_sim, N_BINS)
    out_ptr   <- safe_ptr(t_norm_grid, y_sim, N_BINS)
    out_ncs   <- safe_ncs(t_norm_grid, y_sim, N_BINS)
    out_gauss <- safe_gauss(d_mm,        y_sim, N_BINS)
    fits_array[r, , "LocalCubic"]     <- out_lp$fit
    fits_array[r, , "PTR"]            <- out_ptr$fit
    fits_array[r, , "NatCubicSpline"] <- out_ncs$fit
    fits_array[r, , "Gaussian"]       <- out_gauss$fit
    st[r, "LocalCubic"]     <- out_lp$status
    st[r, "PTR"]            <- out_ptr$status
    st[r, "NatCubicSpline"] <- out_ncs$status
    st[r, "Gaussian"]       <- out_gauss$status
  }
  list(fits = fits_array, status = st)
}

# ---------------------------------------------------------------------
# 9. Resume from checkpoint if available
# ---------------------------------------------------------------------
if (file.exists(PARTIAL_FILE)) {
  cat("Resume: found", PARTIAL_FILE, "- loading completed cells...\n")
  checkpoint <- readRDS(PARTIAL_FILE)
  results_flat <- checkpoint$results
  status_flat  <- checkpoint$status
  cat("Resumed with", length(results_flat), "cells already complete.\n")
} else {
  results_flat <- list()
  status_flat  <- list()
}

# Build the full task list
tasks <- expand.grid(
  Model   = NOISE_MODELS,
  CV      = CVS,
  Family  = FAMILIES,
  Subject = seq_len(N_SUBJ),
  stringsAsFactors = FALSE
)
tasks$key <- mapply(cell_key, tasks$Model, tasks$CV, tasks$Family, tasks$Subject)
tasks <- tasks[!tasks$key %in% names(results_flat), , drop = FALSE]

n_total_cells <- length(NOISE_MODELS) * length(CVS) * length(FAMILIES) * N_SUBJ
cat(sprintf("\n%d / %d cells remaining (%d already done).\n",
            nrow(tasks), n_total_cells, n_total_cells - nrow(tasks)))
cat(sprintf("Each cell = B=%d replicates x 4 methods.\n", B))

# ---------------------------------------------------------------------
# 10. Main loop with ETA reporting and per-cell checkpointing
# ---------------------------------------------------------------------
t_overall_start <- Sys.time()
fmt_eta <- function(sec) {
  if (!is.finite(sec) || sec < 0) return("--")
  h <- floor(sec / 3600); m <- floor((sec - 3600*h) / 60)
  if (h >= 24) sprintf("%dd %dh", h %/% 24, h %% 24) else sprintf("%dh %dm", h, m)
}

cells_done_this_session <- 0L
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  t0 <- Sys.time()
  out <- run_one_cell(task$Family, task$Model, task$CV, task$Subject)
  results_flat[[task$key]] <- out$fits
  status_flat[[task$key]]  <- out$status
  cells_done_this_session <- cells_done_this_session + 1L
  
  saveRDS(list(results = results_flat, status = status_flat), PARTIAL_FILE)
  
  ok_counts <- colSums(out$status == "ok")
  t_elapsed_session <- as.numeric(difftime(Sys.time(), t_overall_start, units = "secs"))
  cells_remaining   <- nrow(tasks) - i
  sec_per_cell      <- t_elapsed_session / cells_done_this_session
  eta_sec           <- sec_per_cell * cells_remaining
  cat(sprintf(
    "[%4d/%4d] %-10s c=%.2f %-9s s=%02d | LP=%d PTR=%d NCS=%d Gauss=%d | %.0fs/cell | ETA %s\n",
    i, nrow(tasks),
    task$Model, task$CV, task$Family, task$Subject,
    ok_counts["LocalCubic"], ok_counts["PTR"],
    ok_counts["NatCubicSpline"], ok_counts["Gaussian"],
    sec_per_cell, fmt_eta(eta_sec)))
}
cat(sprintf("\nDone. Total wall clock: %s\n",
            fmt_eta(as.numeric(difftime(Sys.time(), t_overall_start, units = "secs")))))

saveRDS(list(results = results_flat, status = status_flat), FINAL_FILE)

# ---------------------------------------------------------------------
# 11. Reshape flat results into the nested structure for metrics
# ---------------------------------------------------------------------
results_all <- list(); status_all <- list()
for (key in names(results_flat)) {
  parts <- strsplit(key, "\\|", fixed = FALSE)[[1]]
  model <- parts[1]
  cv_key <- parts[2]
  fam    <- parts[3]
  s_str  <- sub("^s_", "", parts[4])
  s_num  <- as.integer(s_str)
  results_all[[model]][[cv_key]][[fam]][[as.character(s_num)]] <- results_flat[[key]]
  status_all [[model]][[cv_key]][[fam]][[as.character(s_num)]] <- status_flat [[key]]
}

# ---------------------------------------------------------------------
# 12. Failure-mode aggregation
# ---------------------------------------------------------------------
cat("\n--- Failure-mode breakdown ---\n")
status_rows <- data.frame()
for (model in NOISE_MODELS) {
  for (cv in CVS) {
    cv_key <- sprintf("cv_%.2f", cv)
    for (fam in FAMILIES) {
      for (s in seq_len(N_SUBJ)) {
        st <- status_all[[model]][[cv_key]][[fam]][[as.character(s)]]
        if (is.null(st)) next
        for (m in METHODS) {
          tab <- table(st[, m])
          for (lab in names(tab)) {
            status_rows <- rbind(status_rows, data.frame(
              NoiseModel = model, CV = cv, Family = fam, Subject = s,
              Method = m, Status = lab, Count = as.integer(tab[lab])))
          }
        }
      }
    }
  }
}
fail_tab <- status_rows %>%
  group_by(NoiseModel, CV, Family, Method, Status) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  arrange(NoiseModel, CV, Family, Method, desc(n))
write.csv(fail_tab, "synthetic_bootstrap_v3_failure_modes.csv", row.names = FALSE)
non_ok <- fail_tab %>% filter(Status != "ok")
if (nrow(non_ok) == 0) cat("All cells: 100% OK across all replicates.\n")
else print(non_ok %>% head(50), n = Inf)

# ---------------------------------------------------------------------
# 13. Metrics vs known true curve
# ---------------------------------------------------------------------
compute_metrics <- function(fits_array, y_true) {
  out <- data.frame()
  for (m in METHODS) {
    mat <- fits_array[, , m]
    mat <- mat[complete.cases(mat), , drop = FALSE]
    if (nrow(mat) < 2) {
      out <- rbind(out, data.frame(Method = m, IMSE = NA, IBias2 = NA,
                                   IVar = NA, MaxAbsBias = NA,
                                   n_valid_reps = nrow(mat))); next
    }
    mean_fit <- colMeans(mat); bias_pt <- mean_fit - y_true
    var_pt   <- apply(mat, 2, var); mse_pt <- bias_pt^2 + var_pt
    out <- rbind(out, data.frame(
      Method = m, IMSE = mean(mse_pt), IBias2 = mean(bias_pt^2),
      IVar = mean(var_pt), MaxAbsBias = max(abs(bias_pt)),
      n_valid_reps = nrow(mat)))
  }
  out
}

syn_metrics <- data.frame()
for (model in NOISE_MODELS) {
  for (cv in CVS) {
    cv_key <- sprintf("cv_%.2f", cv)
    for (fam in FAMILIES) {
      for (s in seq_len(N_SUBJ)) {
        tc <- true_curves_syn[[fam]][[s]]
        fa <- results_all[[model]][[cv_key]][[fam]][[as.character(s)]]
        if (is.null(fa)) next
        m_df <- compute_metrics(fa, tc$y_true)
        m_df$NoiseModel <- model; m_df$CV <- cv
        m_df$Family <- fam; m_df$Subject <- s
        syn_metrics <- rbind(syn_metrics, m_df)
      }
    }
  }
}
saveRDS(syn_metrics, "synthetic_bootstrap_v3_metrics.rds")

syn_summary <- syn_metrics %>%
  group_by(NoiseModel, CV, Family, Method) %>%
  summarise(mean_IMSE = mean(IMSE, na.rm = TRUE),
            sd_IMSE   = sd(IMSE,   na.rm = TRUE),
            mean_IBias2 = mean(IBias2, na.rm = TRUE),
            mean_IVar   = mean(IVar,   na.rm = TRUE),
            bias2_share = mean(IBias2 / IMSE, na.rm = TRUE),
            n_subj_used = sum(!is.na(IMSE)),
            .groups = "drop") %>%
  mutate(Method     = factor(Method, levels = METHODS),
         NoiseModel = factor(NoiseModel, levels = NOISE_MODELS))
write.csv(syn_summary, "synthetic_bootstrap_v3_summary.csv", row.names = FALSE)
write.csv(syn_metrics, "synthetic_bootstrap_v3_metrics.csv", row.names = FALSE)
cat("\n--- v3 summary (head) ---\n"); print(head(syn_summary, 28), n = 28)

# ---------------------------------------------------------------------
# 14. Plots
# ---------------------------------------------------------------------
cv_lab <- function(x) paste0("c = ", x)

# Crossover line plot: IMSE vs CV for each method, faceted (NoiseModel x Family)
p_cross <- ggplot(syn_summary, aes(CV, mean_IMSE, colour = Method)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_grid(NoiseModel ~ Family, scales = "free_y") +
  scale_x_continuous(breaks = CVS) +
  labs(title = "Synthetic-truth bootstrap (v3): IMSE vs CV",
       subtitle = sprintf("Mean across 21 subjects, B=%d, matched variance", B),
       x = "Coefficient of variation c", y = "Integrated MSE") +
  theme_bw()
ggsave("syn_v3_imse_crossover.png", p_cross,
       width = 10, height = 11, dpi = 200)

# Headline grids per family (bars at each CV)
for (fam in FAMILIES) {
  d_fam <- syn_summary %>% filter(Family == fam)
  p <- ggplot(d_fam, aes(x = factor(CV), y = mean_IMSE, fill = Method)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = pmax(mean_IMSE - sd_IMSE, 0),
                      ymax = mean_IMSE + sd_IMSE),
                  position = position_dodge(width = 0.9), width = 0.25) +
    facet_wrap(~ NoiseModel, ncol = 2) +
    labs(title = sprintf("v3: IMSE by method, family = %s", fam),
         subtitle = sprintf("Mean +/- 1 SD across 21 subjects, B=%d", B),
         x = "Coefficient of variation c", y = "IMSE") +
    theme_bw()
  ggsave(sprintf("syn_v3_imse_%s.png", fam), p,
         width = 11, height = 8, dpi = 200)
}

# Bias^2 share vs CV: shows the crossover mechanism
share_df <- syn_summary %>%
  mutate(bias_share = mean_IBias2 / mean_IMSE)
p_share <- ggplot(share_df, aes(CV, bias_share, colour = Method)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_grid(NoiseModel ~ Family) +
  scale_x_continuous(breaks = CVS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "v3: Bias^2 share of IMSE vs CV",
       x = "Coefficient of variation c", y = "Bias^2 / IMSE") +
  theme_bw()
ggsave("syn_v3_bias_share.png", p_share, width = 10, height = 11, dpi = 200)

# Decomposition stacked bars per family
decomp <- syn_summary %>%
  select(NoiseModel, CV, Family, Method, mean_IBias2, mean_IVar) %>%
  pivot_longer(c(mean_IBias2, mean_IVar),
               names_to = "Component", values_to = "Value") %>%
  mutate(Component = recode(Component,
                            mean_IBias2 = "Integrated Bias^2",
                            mean_IVar   = "Integrated Variance"))
for (fam in FAMILIES) {
  d_fam <- decomp %>% filter(Family == fam)
  p <- ggplot(d_fam, aes(x = Method, y = Value, fill = Component)) +
    geom_col() +
    facet_grid(NoiseModel ~ CV, labeller = labeller(CV = cv_lab)) +
    scale_fill_manual(values = c("Integrated Bias^2"   = "#E07B39",
                                 "Integrated Variance" = "#5B8DB8")) +
    labs(title = sprintf("v3: bias^2 vs variance decomposition, family = %s", fam),
         x = NULL, y = "IMSE component") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(sprintf("syn_v3_decomp_%s.png", fam), p,
         width = 16, height = 8, dpi = 200)
}

cat("\nPriority 1 v3 outputs saved.\n")