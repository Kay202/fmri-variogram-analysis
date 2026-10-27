# Regenerate the intermediate results, rather than trusting the shipped ones.
#
# Three things run, in this order:
#
#   1. Simulation.R   the simulation study, hours, skippable
#   2. Analysis.R     the rest of the R analysis, about 4 minutes
#   3. Inference.m    the inference experiment, MATLAB, about 33 minutes
#
# This script runs 1 and 2, then prints the command for 3 and stops. It does not
# drive MATLAB; re-invoke it afterwards and it confirms the outputs are present.
#
# Inference.m needs no imaging data and resumes by itself: each of its four
# stages is skipped when its CSV already exists.
#
# Usage
#   Rscript RunAnalysis.R
#   Rscript RunAnalysis.R --skip-simulation

argv  <- commandArgs(trailingOnly = TRUE)
afull <- commandArgs(trailingOnly = FALSE)
fa    <- sub("^--file=", "", afull[grep("^--file=", afull)])
HERE  <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
SCR   <- file.path(HERE, "codes")
RES   <- file.path(HERE, "data", "intermediates")
RSCRIPT <- file.path(R.home("bin"), "Rscript")
skip_sim <- "--skip-simulation" %in% argv

runR <- function(s, produces, note = "", args = character(0)) {
  if (all(file.exists(file.path(RES, produces)))) {
    cat(sprintf("  [have]  %-16s %s\n", s, note)); return(invisible(TRUE))
  }
  cat(sprintf("  [run ]  %-16s %s\n", s, note))
  st <- Sys.time()
  rc <- system2(RSCRIPT, c(shQuote(file.path(SCR, s)), args),
                stdout = NULL, stderr = NULL)
  cat(sprintf("          exit %d after %.0f s\n", rc,
              as.numeric(difftime(Sys.time(), st, units = "secs"))))
  # Exit 0 alone is not success: a script can run a self-test and write nothing.
  made <- file.exists(file.path(RES, produces))
  if (rc == 0 && !all(made)) {
    cat(sprintf("          exit 0 but produced nothing: %s\n",
                paste(produces[!made], collapse = ", ")))
    return(invisible(FALSE))
  }
  invisible(rc == 0)
}

needM <- function(s, produces, mins, note = "") {
  if (all(file.exists(file.path(RES, produces)))) {
    cat(sprintf("  [have]  %-16s output present\n", s)); return(invisible(TRUE))
  }
  cat(sprintf("\n  %s is a MATLAB stage. Run it, then re-invoke this script:\n\n", s))
  cat(sprintf("    matlab -batch \"cd('%s'); %s\"\n\n", SCR, sub("\\.m$", "", s)))
  cat(sprintf("  about %d minutes; writes %s\n", mins, paste(produces, collapse = ", ")))
  invisible(FALSE)
}

cat("regenerating the intermediates in data/intermediates\n\n")

cat("1. simulation study\n")
if (skip_sim) {
  cat("  [skip]  Simulation.R     --skip-simulation given\n")
} else {
  # "run" is required. Simulation.R takes a mode argument and defaults to
  # "bench", a timing benchmark that writes nothing and exits 0.
  runR("Simulation.R", c("SIM2_metrics.csv", "SIM2_curves.rds"),
       "several hours on four cores", args = "run")
}

# The files checked are the last ones Analysis.R writes, so a partial run is not
# mistaken for a complete one.
cat("\n2. shell models, interpolation signature, depth gradient, reference, sweep\n")
ok <- runR("Analysis.R",
           c("shell_model_params_v2.csv", "S44_gaussian_diagnostics.csv",
             "S17_guarded.csv", "S40_depth_vs_kernel.csv",
             "A2_depth_gradient.csv", "A6_depth_gradient_fixed.csv",
             "S25_reference_params.csv", "S27_ref_matched_params.csv",
             "S34_bands.csv", "S34_sweep_params.csv"),
           "sixteen analyses, about 4 min")

cat("\n3. inference experiment\n")
if (ok) ok <- needM("Inference.m",
                    c("S26_reference_fidelity.csv", "S32_e1e10_fwer_final.csv",
                      "S35_sweep_thresholds.csv", "S29_e1e10_power_ref.csv"), 33,
                    "reference fidelity, error rates, thresholds, power")

cat("\noptional: GeneratorCheck.m, why the reference is analytic rather than\n")
cat("shell-binned. Produces no manuscript table.\n")
cat(sprintf("    matlab -batch \"cd('%s'); GeneratorCheck\"\n", SCR))

if (ok) {
  cat("\ndone. Now run:\n  Rscript RunOutputs.R\n")
} else {
  cat("\nstopped at the MATLAB stage. Run it as printed, then re-invoke.\n")
}
