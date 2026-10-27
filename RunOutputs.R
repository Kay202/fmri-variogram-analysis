# Regenerate every table and figure in the manuscript.
#
# Reads data/intermediates and data/out_matlab, writes outputs/tables and
# outputs/figures. R only, a few seconds, no simulation.
#
#   Rscript RunOutputs.R
#
# To regenerate the intermediates rather than trusting the shipped ones, run
# RunAnalysis.R first. The imaging stage cannot be included, because the
# residuals alone are 263 GB; see README.md.
#
# After each script this counts what appeared, and at the end it names anything
# the manuscript expects and does not have: a silent partial success is the
# failure worth guarding against.

argv <- commandArgs(trailingOnly = FALSE)
fa   <- sub("^--file=", "", argv[grep("^--file=", argv)])
HERE <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
SCR  <- file.path(HERE, "codes")
TAB  <- file.path(HERE, "outputs", "tables")
FIG  <- file.path(HERE, "outputs", "figures")
RSCRIPT <- file.path(R.home("bin"), "Rscript")

scripts <- c("Tables.R", "Figures.R")

expected_tables <- c("sim_priors.tex", "sim_imse_main.tex", "sim_stability.tex",
  "sim_amplitude.tex", "sim_sillfrac.tex", "sim_cvsweep.tex", "sim_imse_gaussian.tex",
  "sim_imse_lognormal.tex", "sim_bootci.tex", "inf_fwer.tex", "inf_power.tex",
  "inf_bound.tex", "res_profile.tex", "res_repair.tex", "res_depth.tex", "depth.tex")
expected_figures <- c(sprintf("figure%02d_%s.pdf", c(1,2,3,4,5,6),
    c("variograms","four_zones","directional","disparity","rss","residual_comparison")),
  "figure07_resampling.pdf", "figure09_ncs_edf.pdf",
  sprintf("figure%d_%s.pdf", 10:14, c("ncs_fits","ptr_fits","lp_fits","gaussian_fits","gcv_curves")),
  sprintf("figure%d_%s.pdf", 15:18, c("ncs_residuals","ptr_residuals","lp_residuals","gaussian_residuals")))

cat(sprintf("R %s, scripts in %s\n\n", getRversion(), SCR))
t0 <- Sys.time(); failed <- character(0)
for (s in scripts) {
  cat(sprintf("--- %s\n", s))
  st <- Sys.time()
  rc <- system2(RSCRIPT, shQuote(file.path(SCR, s)), stdout = NULL, stderr = NULL)
  nt <- length(list.files(TAB, pattern = "\\.tex$"))
  nf <- length(list.files(FIG, pattern = "\\.pdf$"))
  cat(sprintf("    exit %d, %.0f s, tables now %d, figures now %d\n",
              rc, as.numeric(difftime(Sys.time(), st, units = "secs")), nt, nf))
  if (rc != 0) failed <- c(failed, s)
}

cat(sprintf("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
miss_t <- setdiff(expected_tables, list.files(TAB))
miss_f <- setdiff(expected_figures, list.files(FIG))
cat(sprintf("\ntables   %d of %d present\n", length(expected_tables) - length(miss_t),
            length(expected_tables)))
if (length(miss_t)) cat("  MISSING:", paste(miss_t, collapse = ", "), "\n")
cat(sprintf("figures  %d of %d present\n", length(expected_figures) - length(miss_f),
            length(expected_figures)))
if (length(miss_f)) cat("  MISSING:", paste(miss_f, collapse = ", "), "\n")
if (length(failed)) cat("\nscripts that exited non-zero:", paste(failed, collapse = ", "), "\n")
if (!length(miss_t) && !length(miss_f) && !length(failed))
  cat("\nEverything the manuscript uses regenerated.\n")
