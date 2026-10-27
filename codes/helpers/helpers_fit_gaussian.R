# =============================================================================
# helpers_fit_gaussian.R
#
# Supplies fit_gaussian(), which "Boostrap analysis.R" calls but does not
# define. Source this after "Main analysis_Model Fitting.R" and before
# "Boostrap analysis.R" (step 4 of Instructions.txt).
#
# PROVENANCE. This definition is reproduced verbatim from ACT3.R line 158, the
# development script that was in the session when the published simulation was
# run. It is kept unchanged so that Table 3 and Table 4 reproduce exactly. The
# bounds matter: a looser variant of fit_gaussian() also exists in
# "Formal Sign Test.R", and substituting it changes the published numbers
#
#     source of fit_gaussian    DoubleExp IMSE   Matern IMSE
#     this file (from ACT3.R)   0.00294          0.00925
#     Formal Sign Test.R        0.00206          0.00473
#     published Table 3         0.002758         0.008809
#
# so "Formal Sign Test.R" must be sourced AFTER the simulation, not before.
# Do not "improve" the bounds here without re-running Act 2 and reporting the
# change.
#
# Requires: minpack.lm
# =============================================================================

# identical to the definition in "Main analysis_Model Fitting.R" line 122;
# repeated here so this file can be sourced on its own.
gaussian_variogram <- function(d, nugget, sill, range) {
  nugget + (sill - nugget) * (1 - exp(-(d/range)^2))
}

fit_gaussian <- function(d, y) {
  tryCatch(
    nlsLM(y ~ gaussian_variogram(d, nugget, sill, range),
          start = list(nugget = max(0, min(y) * 0.5), sill = max(y) * 0.95, range = 4),
          lower = c(nugget = 0,            sill = mean(y) * 0.5, range = 0.5),
          upper = c(nugget = min(y) * 1.5, sill = max(y) * 2,    range = 30),
          control = nls.lm.control(maxiter = 300)),
    error = function(e) NULL)
}
