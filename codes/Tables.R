# Every LaTeX table fragment in the manuscript.
#
# Reads the intermediate result tables and writes the fragments to the output
# directory. Both locations are set in the preamble below and need no editing.

# Paths, packages and the shared helpers. Each analysis below runs inside
# local() and seeds its own RNG where it needs one.
suppressPackageStartupMessages({
  library(data.table); library(splines); library(nnls)
  library(minpack.lm); library(readxl)
})
argv <- commandArgs(trailingOnly = FALSE)
fa   <- sub("^--file=", "", argv[grep("^--file=", argv)])
HERE <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
REL  <- file.path(HERE, "helpers")
DIAG <- file.path(HERE, "..", "data")
OUT  <- file.path(HERE, "..", "data", "intermediates")
TEXD <- file.path(HERE, "..", "outputs", "tables")
TAB  <- TEXD
FIG  <- file.path(HERE, "..", "outputs", "figures")
PS   <- file.path(HERE, "..", "outputs")
for (d in c(OUT, TEXD, FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
VOX <- 3

# the model-fitting functions from the published analysis
src <- function(p) { e <- parse(p)
  for (x in e) if (is.call(x) && length(x) >= 3 &&
      as.character(x[[1]])[1] %in% c("<-","=") && is.call(x[[3]]) &&
      identical(as.character(x[[3]][[1]])[1], "function")) eval(x, envir = globalenv()) }
src(file.path(REL, "Main analysis_Model Fitting.R"))

# separable Lanczos-3 kernel, and its autocorrelation averaged over
# sub-voxel offsets at integer lags 0 to 4
lanczos <- function(x, a = 3) { y <- numeric(length(x)); nz <- abs(x) < a & x != 0
  xx <- x[nz]; y[nz] <- (sin(pi*xx)/(pi*xx))*(sin(pi*xx/a)/(pi*xx/a)); y[x == 0] <- 1; y }
AK <- local({ off <- seq(0, 1, length.out = 201)[-201]; acc <- numeric(5)
  for (dd in off) { w <- lanczos((-6):6 - dd); den <- sum(w^2)
    for (i in 1:5) { h <- i - 1L; n <- seq_len(length(w) - h)
      acc[i] <- acc[i] + sum(w[n]*w[n + h])/den } }
  acc/length(off) })

# Simulation study
local({
LOG  <- file.path(OUT, "SIM2_report_log.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = LOG, append = TRUE) }
cat("", file = LOG)

S  <- readRDS(file.path(OUT, "SIM2_curves.rds")); DS <- S$design
R  <- fread(file.path(OUT, "SIM2_metrics.csv"))
NM <- length(DS$METHODS)

# reporting choices live here, not in the simulation
N_REPORT   <- c(10L, 20L, 40L, 60L, 80L, 120L)
HEAD_NOISE <- "ChiSquared"      # binned semivariance is a mean of squared differences
HEAD_CV    <- 0.05
BIAS       <- "IBias2_c"        # finite-replicate corrected, not raw
LAB <- c(LocalCubic = "LP", PTR = "PTR", NatCubicSpline = "NCS",
         ShapiroBotha = "SB", Gaussian = "Gaussian")
FAMLAB <- c(DoubleExp = "Double-exponential", Matern = "Mat\\'ern")
fmt <- function(x, d = 6) formatC(x, format = "f", digits = d)
wr  <- function(f, l) writeLines(l, file.path(TEXD, f))
mrow <- function(n, f) sprintf("\\multirow{%d}{*}{%s}", n, FAMLAB[f])

say("SIM2 report, %s\n", format(Sys.Date()))
say("%d curves per family, B = %d, %d lags, %d noise models, %d levels, %d methods\n\n",
    DS$N_CURVE, DS$B, length(DS$D_MM), length(DS$NOISE), length(DS$CVS), NM)

# 0. did anything fail to fit?
say("== 0. fitting integrity ==\n")
for (m in DS$METHODS) {
  x <- R[method == m]
  say("  %-15s lost replicates in %4d of %4d cells | worst cell kept %3d of %d | %d unusable\n",
      m, sum(x$n_ok < DS$B), nrow(x), min(x$n_ok), DS$B, sum(is.na(x$IMSE)))
}
if ("gauss_atbound" %in% names(R)) {
  g <- R[method == "Gaussian"]
  say("  Gaussian fits finishing on a constraint boundary: %.2f%%\n",
      100*sum(g$gauss_atbound, na.rm = TRUE)/sum(g$n_ok, na.rm = TRUE))
}

# 1. priors
pr <- DS$PRIOR
role <- c(s = "sill fraction reached at the longest lag",
          sigma2 = "structural sill (scale convention)",
          A = "oscillatory amplitude", omega = "oscillatory angular frequency",
          b = "damping length", tau = "nugget")
sym <- c(s = "s", sigma2 = "\\sigma^2", A = "A", omega = "\\omega", b = "b", tau = "\\tau")
L <- c("\\begin{tabular}{llc}", "\\toprule", "Parameter & Role & Prior \\\\", "\\midrule")
for (k in names(pr))
  L <- c(L, sprintf("$%s$ & %s & $U(%.2f,\\ %.3f)$ \\\\", sym[k], role[k], pr[[k]][1], pr[[k]][2]))
wr("sim_priors.tex", c(L, "\\bottomrule", "\\end{tabular}"))

# 2. main table
H <- R[noise == HEAD_NOISE & abs(cv - HEAD_CV) < 1e-9]
agg <- H[, .(IMSE = mean(IMSE, na.rm=TRUE), IB = mean(get(BIAS), na.rm=TRUE),
             IV = mean(IVar, na.rm=TRUE),
             share = 100*mean(get(BIAS)/IMSE, na.rm=TRUE)), by = .(family, method)]
L <- c("\\begin{tabular}{llcccc}", "\\toprule",
       "Family & Method & IMSE & IBias$^2$ & IVar & Bias$^2$\\,\\% \\\\", "\\midrule")
for (f in DS$FAMILIES) {
  L <- c(L, mrow(NM, f))
  for (m in DS$METHODS) {
    x <- agg[family == f & method == m]
    L <- c(L, sprintf("& %-8s & %s & %s & %s & %.1f \\\\", LAB[m],
                      fmt(x$IMSE), fmt(x$IB), fmt(x$IV), x$share))
  }
  if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
}
wr("sim_imse_main.tex", c(L, "\\bottomrule", "\\end{tabular}"))
say("\n== 2. headline cell, %s at c = %.2f ==\n", HEAD_NOISE, HEAD_CV)
for (i in seq_len(nrow(agg)))
  say("  %-10s %-15s IMSE %s  bias share %5.1f%%\n",
      agg$family[i], agg$method[i], fmt(agg$IMSE[i]), agg$share[i])
for (f in DS$FAMILIES) {
  x <- agg[family == f]
  g <- x$IMSE[x$method == "Gaussian"]; b <- min(x$IMSE[x$method != "Gaussian"])
  say("  %-10s Gaussian / best non-parametric = %.2f\n", f, g/b)
}

# 3. stability in the number of curves
L <- c("\\begin{tabular}{lr" , paste(rep("r", NM), collapse=""), "r}")
L <- c("\\begin{tabular}{lrrrrrrr}", "\\toprule",
       sprintf("Family & $N$ & %s & Ratio \\\\", paste(LAB[DS$METHODS], collapse = " & ")),
       "\\midrule")
say("\n== 3. stability in the number of curves ==\n")
for (f in DS$FAMILIES) {
  L <- c(L, mrow(length(N_REPORT), f))
  for (n in N_REPORT) {
    x <- H[family == f & curve <= n, .(v = mean(IMSE, na.rm=TRUE)), by = method]
    mv <- setNames(x$v, x$method)
    rt <- mv["Gaussian"]/min(mv[names(mv) != "Gaussian"])
    L <- c(L, sprintf(" & %3d & %s & %.2f \\\\", n,
                      paste(fmt(mv[DS$METHODS]), collapse = " & "), rt))
    say("  %-10s N=%3d  ratio %.2f\n", f, n, rt)
  }
  if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
}
wr("sim_stability.tex", c(L, "\\bottomrule", "\\end{tabular}"))

# 4. stratified by curve shape
strat <- function(var, lab, file) {
  W <- dcast(H, family + curve ~ method, value.var = "IMSE")
  W <- merge(W, unique(R[, c("curve", var), with = FALSE]), by = "curve")
  np <- setdiff(DS$METHODS, "Gaussian")
  W[, bestNP := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = np]
  W[, ratio := Gaussian/bestNP]
  # Both the mean and the median of the per-curve ratios, because the ratio is
  # heavy-tailed and the two can point in opposite directions. In the lowest
  # sill-fraction tertile three of the forty curves carry about half the summed
  # ratio, so the mean alone suggests a gradient that no rank test finds.
  # Eight columns do not fit the supplement's text block at full width, so the
  # two ratio columns share a heading rather than carrying "Mean ratio" and
  # "Median ratio", whose headers were wider than any number beneath them, and
  # the column padding is tightened inside a group so it does not leak.
  L <- c("{\\setlength{\\tabcolsep}{3pt}",
         "\\begin{tabular}{llrrrrrr}", "\\toprule",
         sprintf(paste("Family & Tertile of %s & $n$ & Mean & Best non-par. &",
                       "Gaussian & \\multicolumn{2}{c}{Ratio} \\\\"), lab),
         "\\cmidrule(lr){7-8}",
         "& & & & & & Mean & Median \\\\",
         "\\midrule")
  say("\n== stratified by %s ==\n", var)
  for (f in DS$FAMILIES) {
    x <- W[family == f]
    x[, tert := cut(get(var), quantile(get(var), c(0,1/3,2/3,1)),
                    include.lowest = TRUE, labels = c("Low","Middle","High"))]
    L <- c(L, mrow(3, f))
    for (b in levels(x$tert)) {
      y <- x[tert == b]
      L <- c(L, sprintf("& %-6s & %3d & %.3f & %s & %s & %5.2f & %5.2f \\\\", b,
                        nrow(y), mean(y[[var]]), fmt(mean(y$bestNP)),
                        fmt(mean(y$Gaussian)), mean(y$ratio), median(y$ratio)))
      say("  %-10s %-6s n=%3d  mean %s %.3f  ratio mean %.2f median %.2f\n", f, b,
          nrow(y), var, mean(y[[var]]), mean(y$ratio), median(y$ratio))
    }
    if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
    ct <- suppressWarnings(cor.test(x[[var]], x$ratio, method = "spearman"))
    say("    Spearman rho %+.3f, p = %.3g, n = %d\n", ct$estimate, ct$p.value, nrow(x))
  }
  wr(file, c(L, "\\bottomrule", "\\end{tabular}}"))
}
strat("A", "amplitude $A$",      "sim_amplitude.tex")
strat("s", "sill fraction $s$",  "sim_sillfrac.tex")

# 5. noise-level sweep
L <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
       sprintf("Family & $c$ & %s & Ratio & Gauss.\\ bias\\,\\%% \\\\",
               paste(LAB[DS$METHODS], collapse = " & ")), "\\midrule")
say("\n== 5. noise-level sweep, %s ==\n", HEAD_NOISE)
for (f in DS$FAMILIES) {
  L <- c(L, mrow(length(DS$CVS), f))
  for (cvi in DS$CVS) {
    # tolerance, not equality: the levels come from seq() so 0.15 is stored as
    # 0.15000000000000002 in the design but round-trips through the CSV as
    # exactly 0.15, and an equality test silently drops that whole level
    x <- R[noise == HEAD_NOISE & abs(cv - cvi) < 1e-9 & family == f,
           .(v = mean(IMSE, na.rm=TRUE), sh = 100*mean(get(BIAS)/IMSE, na.rm=TRUE)), by = method]
    mv <- setNames(x$v, x$method)
    rt <- mv["Gaussian"]/min(mv[names(mv) != "Gaussian"])
    L <- c(L, sprintf(" & %.2f & %s & %.2f & %.1f \\\\", cvi,
                      paste(fmt(mv[DS$METHODS]), collapse = " & "), rt,
                      x$sh[x$method == "Gaussian"]))
    say("  %-10s c=%.2f  ratio %.2f\n", f, cvi, rt)
  }
  if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
}
wr("sim_cvsweep.tex", c(L, "\\bottomrule", "\\end{tabular}"))

# 6. appendix: other noise models
for (nm in setdiff(DS$NOISE, HEAD_NOISE)) {
  X <- R[noise == nm & abs(cv - HEAD_CV) < 1e-9,
         .(IMSE = mean(IMSE, na.rm=TRUE), IB = mean(get(BIAS), na.rm=TRUE),
           IV = mean(IVar, na.rm=TRUE),
           share = 100*mean(get(BIAS)/IMSE, na.rm=TRUE)), by = .(family, method)]
  L <- c("\\begin{tabular}{llcccc}", "\\toprule",
         "Family & Method & IMSE & IBias$^2$ & IVar & Bias$^2$\\,\\% \\\\", "\\midrule")
  for (f in DS$FAMILIES) {
    L <- c(L, mrow(NM, f))
    for (m in DS$METHODS) {
      x <- X[family == f & method == m]
      L <- c(L, sprintf("& %-8s & %s & %s & %s & %.1f \\\\", LAB[m],
                        fmt(x$IMSE), fmt(x$IB), fmt(x$IV), x$share))
    }
    if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
  }
  wr(sprintf("sim_imse_%s.tex", tolower(nm)), c(L, "\\bottomrule", "\\end{tabular}"))
}

# 7. bootstrap intervals on the bias share
set.seed(11); NB <- 2000
L <- c("\\begin{tabular}{llcc}", "\\toprule",
       "Family & Method & Bias$^2$ share & 95\\% CI \\\\", "\\midrule")
say("\n== 7. bootstrap intervals, %d curve resamples ==\n", NB)
for (f in DS$FAMILIES) {
  L <- c(L, mrow(NM, f))
  for (m in DS$METHODS) {
    v <- H[family == f & method == m, get(BIAS)/IMSE]; v <- v[is.finite(v)]
    bs <- replicate(NB, mean(sample(v, length(v), TRUE)))
    ci <- quantile(bs, c(.025,.975))
    L <- c(L, sprintf("& %-8s & %.1f\\%% & $[%.1f, %.1f]\\%%$ \\\\", LAB[m],
                      100*mean(v), 100*ci[1], 100*ci[2]))
    say("  %-10s %-15s %5.1f%% [%.1f, %.1f]\n", f, m, 100*mean(v), 100*ci[1], 100*ci[2])
  }
  if (f != tail(DS$FAMILIES,1)) L <- c(L, "\\midrule")
}
wr("sim_bootci.tex", c(L, "\\bottomrule", "\\end{tabular}"))

# 8. figures
COL <- c(LocalCubic="#b2182b", PTR="#2166ac", NatCubicSpline="#1b7837",
         ShapiroBotha="#e08214", Gaussian="#7b3294")
pdf(file.path(FIG, "sim_imse_vs_cv.pdf"), width = 8.4, height = 4.0)
par(mfrow = c(1,2), mar = c(4.2,4.4,2.0,0.8), las = 1)
for (f in DS$FAMILIES) {
  X <- R[noise == HEAD_NOISE & family == f, .(v = mean(IMSE, na.rm=TRUE)), by = .(cv, method)]
  plot(range(DS$CVS), range(X$v), type="n", log="y", main = FAMLAB[f],
       xlab = "Coefficient of variation c", ylab = "Mean IMSE")
  for (m in DS$METHODS) { y <- X[method==m][order(cv)]
    lines(y$cv, y$v, col=COL[m], lwd=2); points(y$cv, y$v, col=COL[m], pch=16, cex=.7) }
  if (f == DS$FAMILIES[1])
    legend("bottomright", bty="n", cex=.7, lwd=2, col=COL[DS$METHODS], legend=LAB[DS$METHODS])
}
dev.off()

pdf(file.path(FIG, "sim_ratio_vs_shape.pdf"), width = 8.4, height = 4.0)
par(mfrow = c(1,2), mar = c(4.2,4.4,2.0,0.8), las = 1)
W <- dcast(H, family + curve ~ method, value.var = "IMSE")
W <- merge(W, unique(R[, .(curve, A, s)]), by = "curve")
np <- setdiff(DS$METHODS, "Gaussian")
W[, bestNP := do.call(pmin, c(.SD, na.rm=TRUE)), .SDcols = np][, ratio := Gaussian/bestNP]
for (v in c("A","s")) {
  plot(W[[v]], W$ratio, type="n", xlab = if (v=="A") "Amplitude A" else "Sill fraction s",
       ylab = "Gaussian / best non-parametric",
       main = if (v=="A") "by amplitude" else "by sill fraction")
  abline(h = 1, col="grey55", lty=2)
  for (f in DS$FAMILIES)
    points(W[family==f][[v]], W[family==f]$ratio, pch=16, cex=.6,
           col = adjustcolor(if (f=="DoubleExp") "#b2182b" else "#2166ac", .6))
  legend("topleft", bty="n", cex=.7, pch=16, col=c("#b2182b","#2166ac"),
         legend = FAMLAB[DS$FAMILIES])
}
dev.off()

say("\nwrote %d LaTeX fragments -> %s\n", length(list.files(TEXD, "^sim_.*tex$")), TEXD)
say("wrote 2 figures -> %s\n", FIG)
})

# Error rates, power, and the identifiability bound
local({
TXT  <- file.path(OUT, "S36_prose_numbers.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)
wr   <- function(lines, f) writeLines(lines, file.path(TAB, f), useBytes = TRUE)

F  <- fread(file.path(OUT, "S32_e1e10_fwer_final.csv"))
P  <- fread(file.path(OUT, "S29_e1e10_power_ref.csv"))
S  <- fread(file.path(OUT, "S35_sweep_thresholds.csv"))
B  <- fread(file.path(OUT, "S34_bands.csv"))
FD <- fread(file.path(OUT, "S26_reference_fidelity.csv"))
RP <- fread(file.path(OUT, "S27_ref_matched_params.csv"))
RE <- fread(file.path(OUT, "S27_ref_matched_exp.csv"))

LAB <- c(Gaussian = "Gaussian", Exponential = "Exponential", AFNImixed = "Mixed",
         SB = "Admissible", SBdamped = "Admissible$^\\dagger$", Reference = "Reference")
ORD <- names(LAB)
shells <- sort(unique(F$shell)); dmm <- sort(unique(F$depth_mm))

# Table: realised familywise error rate
M <- F[reference == "matched"]
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{Realised familywise error rate against a nominal 0.05, by assumed",
"autocorrelation model, depth shell and cluster-forming threshold. Null fields are",
"drawn from the reference covariance of Equation~\\eqref{eq:reference}, which no",
"assumed model produces; the Reference column is therefore a check on the",
"machinery rather than a result, and must read 0.05. The two admissible columns",
"differ only in how the field is drawn from the same fitted estimate, $\\dagger$",
"denoting the damped synthesis of Section~\\ref{sec:inference}. $B = 1000$",
"simulations per cell.}",
"\\label{tab:fwer}", "\\begin{tabular}{llrrrrrr}", "\\toprule",
paste0("Shell & $z >$ & ", paste(LAB[ORD], collapse = " & "), " \\\\"), "\\midrule")
for (k in shells) for (z in c(3.09, 2.30)) {
  r <- vapply(ORD, function(m)
    M[shell == k & abs(zthr - z) < 1e-9 & model_name == m]$FWER, numeric(1))
  cell <- sprintf("%.3f", r)
  cell[ORD == "Reference"] <- sprintf("\\textbf{%.3f}", r[ORD == "Reference"])
  L <- c(L, sprintf("%d\\,mm & %.2f & %s \\\\", dmm[shells == k], z,
                    paste(cell, collapse = " & ")))
}
L <- c(L, "\\midrule",
  sprintf("\\multicolumn{2}{l}{Mean} & %s \\\\",
          paste(sprintf("%.3f", vapply(ORD, function(m) mean(M[model_name == m]$FWER),
                                       numeric(1))), collapse = " & ")),
  sprintf("\\multicolumn{2}{l}{Mean, hole reference} & %s \\\\",
          paste(sprintf("%.3f", vapply(ORD, function(m)
            mean(F[reference == "hole" & model_name == m]$FWER), numeric(1))),
            collapse = " & ")),
  "\\bottomrule",
  "\\multicolumn{8}{l}{\\footnotesize The final row repeats the experiment against a reference carrying}\\\\",
  "\\multicolumn{8}{l}{\\footnotesize the deepest hole effect the data admit; no ordering changes.}",
  "\\end{tabular}", "\\end{table}")
wr(L, "inf_fwer.tex")

# Table: power
amps <- c(2, 3, 4)
PM <- P[reference == "matched"]
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{Power at three signal amplitudes under each assumed model, with the",
"realised familywise error rate alongside, because power is only interpretable",
"next to the rate at which it was bought. Cluster-forming threshold $z > 3.09$,",
"signal 15\\,mm FWHM, $B = 500$ simulations per cell. Only the Reference row",
"holds nominal error, and it is not available to a practitioner.}",
"\\label{tab:power}", "\\begin{tabular}{llrrr}", "\\toprule",
"Assumed model & FWER & \\multicolumn{3}{c}{Power at amplitude} \\\\",
"\\cmidrule(lr){3-5}",
paste0(" & & ", paste(sprintf("%d", amps), collapse = " & "), " \\\\"), "\\midrule")
for (m in ORD) {
  fw <- mean(PM[model_name == m]$FWER)
  pw <- vapply(amps, function(a) mean(PM[model_name == m & amplitude == a]$power), numeric(1))
  nm <- LAB[[m]]; if (m == "Reference") nm <- paste0("\\textbf{", nm, "}")
  L <- c(L, sprintf("%s & %.3f & %s \\\\", nm, fw,
                    paste(sprintf("%.3f", pw), collapse = " & ")))
}
L <- c(L, "\\bottomrule",
  "\\multicolumn{5}{l}{\\footnotesize Power rises with the realised error rate across the rows, which is}\\\\",
  "\\multicolumn{5}{l}{\\footnotesize what an uncalibrated threshold produces rather than detection.}",
  "\\end{tabular}", "\\end{table}")
wr(L, "inf_power.tex")

# Table: the identifiability bound
tight <- S[rmse <= 1.5*floor_rmse]
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{What the data leave undetermined. Each row spans every admissible",
"model that reproduces the observed shell variogram as well as any model of this",
"family can, within a factor of 1.5 on the residual. Those models cannot be told",
"apart from the data, yet the cluster-extent threshold they imply, and the error",
"rate that threshold realises against one fixed truth, vary as shown.}",
"\\label{tab:bound}", "\\begin{tabular}{llrrr}", "\\toprule",
"Shell & $z >$ & Models & $c^\\ast$ (voxels) & Realised FWER \\\\", "\\midrule")
for (k in shells) for (z in c(3.09, 2.30)) {
  s <- tight[shell == k & abs(zthr - z) < 1e-9]
  L <- c(L, sprintf("%d\\,mm & %.2f & %d & %.0f--%.0f & %.3f--%.3f \\\\",
                    dmm[shells == k], z, nrow(s), min(s$c_star), max(s$c_star),
                    min(s$FWER), max(s$FWER)))
}
# This range was hardcoded and disagreed with the prose, which read 4 to 423.
# Computed from the data now, so the two cannot drift again.
wide <- S[rmse <= tol]
L <- c(L, "\\bottomrule",
  "\\multicolumn{5}{l}{\\footnotesize Nominal is 0.05 throughout. Widening the tolerance to the standard}\\\\",
  sprintf("\\multicolumn{5}{l}{\\footnotesize error of the cohort curve widens $c^\\ast$ further, to %d--%d voxels.}",
          min(wide$c_star), max(wide$c_star)),
  "\\end{tabular}", "\\end{table}")
wr(L, "inf_bound.tex")

# the numbers quoted in the prose
say("Numbers used in the Section 5 prose. Every one is read from a results file\n")
say("by this script, so the text and the tables cannot drift apart.\n\n")
say("reference fit, residual sum of squares by shell : %s\n",
    paste(sprintf("%.3f", RP$sb_rss + 0), collapse = ", "))
say("  (admissible estimate against the reference, so effectively zero)\n")
say("Gaussian rss against the reference               : %s\n",
    paste(sprintf("%.3f", RP$g_rss), collapse = ", "))
say("mixed rss                                        : %s\n",
    paste(sprintf("%.3f", RP$m_rss), collapse = ", "))
say("exponential rss                                  : %s\n",
    paste(sprintf("%.3f", RE$e_rss), collapse = ", "))
say("exponential range, interior on a 0.2-80 grid     : %s\n",
    paste(sprintf("%.2f", RE$e_range), collapse = ", "))
say("mixed exponential component, grid top is 300     : %s\n",
    paste(sprintf("%.0f", RP$m_c), collapse = ", "))
say("\nsynthesis fidelity, reference against its own model, %% of curve range:\n")
say("  %s\n", paste(sprintf("%.2f", 100*FD$max_err/(FD$curve_range)), collapse = ", "))
say("\nmean realised FWER, matched reference:\n")
for (m in ORD) say("  %-22s %.3f   (hole reference %.3f)\n", LAB[[m]],
                   mean(M[model_name == m]$FWER),
                   mean(F[reference == "hole" & model_name == m]$FWER))
say("\nidentified range of total variance, ratio of largest to smallest:\n")
for (b in unique(B$basis))
  say("  %-34s %.2f to %.2f\n", b, min(B[basis == b]$ratio), max(B[basis == b]$ratio))
say("\nbound at 1.5x the best achievable fit:\n")
say("  critical cluster size spans %d to %d voxels\n", min(tight$c_star), max(tight$c_star))
say("  realised error rate spans %.3f to %.3f\n", min(tight$FWER), max(tight$FWER))
wide <- S[rmse <= tol]
say("  at the cohort standard error, %d to %d voxels\n", min(wide$c_star), max(wide$c_star))
say("\nwrote inf_fwer.tex, inf_power.tex, inf_bound.tex\n")
})

# Interpolation signature
local({
TXT  <- file.path(OUT, "S38_prose_numbers.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)
wr <- function(lines, f) writeLines(lines, file.path(TAB, f), useBytes = TRUE)

C <- fread(file.path(OUT, "S37_claims.csv"))
G <- fread(file.path(OUT, "S17_guarded.csv"))

fmtp <- function(p) ifelse(p < 1e-4, sprintf("$%.0f \\times 10^{%d}$",
  as.numeric(sub("e.*","",formatC(p, format="e", digits=0))),
  as.integer(sub(".*e","",formatC(p, format="e", digits=0)))), sprintf("%.3f", p))

sel <- C[frac %in% c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0)]
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{Three results of Section~\\ref{sec:observed} as a function of how much",
"of the fitted resampling variance is removed. The first row is the uncorrected",
"analysis and reproduces the figures reported above. The last is the largest",
"correction the data admit, at which the field's own variogram reaches zero at",
"the shortest separation. The under-prediction at the corner is unaffected",
"throughout; the fall between 6 and 9\\,mm is not.}",
"\\label{tab:resampling_profile}", "\\begin{tabular}{rrcrcrc}", "\\toprule",
"Removed & \\multicolumn{2}{c}{Fall, 6 to 9\\,mm} & \\multicolumn{2}{c}{Corner zone} & \\multicolumn{2}{c}{Plateau zone} \\\\",
"\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
" & mean & subjects & positive & $p$ & negative & $p$ \\\\", "\\midrule")
for (i in seq_len(nrow(sel)))
  L <- c(L, sprintf("%.0f\\%% & $%+.2f$\\%% & %d of 21 & %d of 21 & %s & %d of 21 & %s \\\\",
    100*sel$frac[i], sel$mean_decline[i], sel$n_pos[i], sel$corner_pos[i],
    fmtp(sel$corner_p[i]), sel$plateau_neg[i], fmtp(sel$plateau_p[i])))
L <- c(L, "\\bottomrule",
"\\multicolumn{7}{l}{\\footnotesize The resampling variance is fitted per subject-run and capped at that run's}\\\\",
"\\multicolumn{7}{l}{\\footnotesize own bound; the fitted value reaches the bound in 207 of 210 runs.}",
"\\end{tabular}", "\\end{table}")
wr(L, "res_profile.tex")

# what removing it does to smoothing-parameter selection
N <- nrow(G)
mk <- function(m, lab, unit) {
  br <- sum(G[[paste0(m, "_bound_raw")]]); ba <- sum(G[[paste0(m, "_bound_adj")]])
  sprintf("%s & %s & %d of 21 & %d of 21 & %.1f & %.1f \\\\", lab, unit, br, ba,
          median(G[[paste0(m, "_edf_raw")]]), median(G[[paste0(m, "_edf_adj")]]))
}
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{Smoothing-parameter selection before and after the resampling term is",
"removed, on the cohort curves. Each criterion is searched over the widest grid",
"it admits and is stopped only where a fit would spend more than half of the 31",
"degrees of freedom, which is the point beyond which generalised cross-validation,",
"dividing by $(1 - \\mathrm{edf}/n)^2$, ceases to be meaningful. Before the",
"correction the criteria run to that limit; afterwards they choose interior",
"values.}",
"\\label{tab:resampling_repair}", "\\begin{tabular}{llccrr}", "\\toprule",
" & & \\multicolumn{2}{c}{Stopped by the limit} & \\multicolumn{2}{c}{Median df spent} \\\\",
"\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
"Method & Parameter & raw & corrected & raw & corrected \\\\", "\\midrule",
mk("LP",  "LP",  "bandwidth"),
mk("PTR", "PTR", "$k$"),
mk("NCS", "NCS", "$\\log_{10}\\lambda$"),
"\\bottomrule",
"\\multicolumn{6}{l}{\\footnotesize Of 31 degrees of freedom available.}",
"\\end{tabular}", "\\end{table}")
wr(L, "res_repair.tex")

# the depth gradient against the correction (S40)
DP <- fread(file.path(OUT, "S40_depth_vs_kernel.csv"))
dps <- sort(unique(DP$depth_mm))
L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{The depth gradient of Table~\\ref{tab:depth} recomputed with the",
"resampling term removed, as a fraction of each shell's own admissibility",
"ceiling. At 6 and 9\\,mm the kernel products are lattice constants, identical in",
"every shell, so the correction enters a shell only through that shell's level",
"and cannot generate differences between shells. The gradient is monotone at",
"every level in both sessions, and steepens rather than flattens; what moves is",
"the depth at which the change crosses zero.}",
"\\label{tab:depth_corrected}",
paste0("\\begin{tabular}{lr", paste(rep("r", length(dps)), collapse = ""), "r}"),
"\\toprule",
sprintf("Session & Removed & %s & Slope \\\\", paste(sprintf("%d\\,mm", dps), collapse = " & ")),
"\\midrule")
for (s in 1:2) for (fr in c(0, 0.5, 1)) {
  z <- DP[session == s & abs(frac - fr) < 1e-9][order(depth_mm)]
  L <- c(L, sprintf("%s & %.0f\\%% & %s & $%+.2f$ \\\\",
    ifelse(fr == 0, sprintf("%d", s), ""), 100*fr,
    paste(sprintf("$%+.2f$", z$mean_ch), collapse = " & "), z$slope[1]))
  if (fr == 1 && s == 1) L <- c(L, "\\midrule")
}
L <- c(L, "\\bottomrule",
"\\multicolumn{8}{l}{\\footnotesize Change in semivariance from 6 to 9\\,mm, per cent, positive meaning the}\\\\",
"\\multicolumn{8}{l}{\\footnotesize curve falls. Slope is percentage points per millimetre of depth.}",
"\\end{tabular}", "\\end{table}")
wr(L, "res_depth.tex")

say("Numbers quoted in the Section 4 resampling prose.\n\n")
say("depth gradient slope per mm, session 1: %s\n",
    paste(sprintf("%.2f at %.0f%%", DP[session==1][order(frac)][, .(s=slope[1]), by=frac]$s,
                  100*unique(DP$frac)), collapse = "; "))
say("monotone at every level and session: %s\n\n", all(DP$monotone))
say("fall between 6 and 9 mm, by fraction removed:\n")
for (i in seq_len(nrow(C)))
  say("  %3.0f%%  %+6.2f%%  %2d of 21  sign p %.2e\n", 100*C$frac[i],
      C$mean_decline[i], C$n_pos[i], C$p_sign[i])
say("\ncorner zone positive: %s\n", paste(C$corner_pos, collapse = ", "))
say("plateau zone negative: %s\n", paste(C$plateau_neg, collapse = ", "))
say("\nsmoothing selection, stopped by the half-df limit:\n")
for (m in c("LP","PTR","NCS"))
  say("  %-4s raw %d of 21, corrected %d of 21; median df %.1f then %.1f\n", m,
      sum(G[[paste0(m,"_bound_raw")]]), sum(G[[paste0(m,"_bound_adj")]]),
      median(G[[paste0(m,"_edf_raw")]]), median(G[[paste0(m,"_edf_adj")]]))
say("\nwrote res_profile.tex, res_repair.tex\n")
})

# Depth gradient
local({
set.seed(20261003L)

A <- fread(file.path(OUT, "A5_depth_fixed.csv"))
W <- dcast(A[k <= 3, .(g = mean(gamma)), by = .(subject, run, session, depth_mm, k)],
           subject + run + session + depth_mm ~ k, value.var = "g")
setnames(W, c("1","2","3"), c("g3","g6","g9"))
W[, ch := 100*(g6 - g9)/g6]
S <- W[, .(ch = mean(ch)), by = .(subject, session, depth_mm)]

# the two alternative conventions, for the record
W2 <- dcast(A[k <= 3], subject + run + session + depth_mm + axis ~ k, value.var = "gamma")
setnames(W2, c("1","2","3"), c("g3","g6","g9"))
S2 <- W2[, .(ch = mean(100*(g6-g9)/g6)), by = .(subject, run, session, depth_mm)][
        , .(ch = mean(ch)), by = .(subject, session, depth_mm)]
W3 <- dcast(A[k <= 3, .(g = sum(gamma*n_base_valid)/sum(n_base_valid)),
              by = .(subject, run, session, depth_mm, k)],
            subject + run + session + depth_mm ~ k, value.var = "g")
setnames(W3, c("1","2","3"), c("g3","g6","g9"))
S3 <- W3[, .(ch = mean(100*(g6-g9)/g6)), by = .(subject, session, depth_mm)]

boot_ci <- function(x, B = 10000) {
  m <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  quantile(m, c(0.025, 0.975), names = FALSE)
}
depths <- sort(unique(S$depth_mm))
lab <- c("9--12\\,mm", "12--15\\,mm", "15--18\\,mm", "18--24\\,mm", "$\\geq 24$\\,mm")

rows <- list()
for (i in seq_along(depths)) {
  r <- list(depth = depths[i], label = lab[i])
  for (s in 1:2) {
    x <- S[session == s & depth_mm == depths[i]]$ch
    ci <- boot_ci(x)
    r[[paste0("ch", s)]] <- mean(x); r[[paste0("lo", s)]] <- ci[1]
    r[[paste0("hi", s)]] <- ci[2];   r[[paste0("nd", s)]] <- sum(x > 0)
    r[[paste0("n", s)]]  <- length(x)
  }
  r$alt2 <- mean(S2[session == 1 & depth_mm == depths[i]]$ch)
  r$alt3 <- mean(S3[session == 1 & depth_mm == depths[i]]$ch)
  rows[[i]] <- as.data.table(r)
}
D <- rbindlist(rows)
fwrite(D, file.path(OUT, "S41_depth.csv"))

cat(sprintf("%-14s %9s %9s %9s %9s %10s %10s\n", "shell", "sess1", "sess2",
            "decl 1", "decl 2", "alt conv 2", "alt conv 3"))
for (i in seq_len(nrow(D)))
  cat(sprintf("%-14s %+8.2f%% %+8.2f%% %6d/%-2d %6d/%-2d %+9.2f%% %+9.2f%%\n",
      depths[i], D$ch1[i], D$ch2[i], D$nd1[i], D$n1[i], D$nd2[i], D$n2[i],
      D$alt2[i], D$alt3[i]))
cat(sprintf("\nlargest departure of either alternative convention: %.3f pp\n",
            max(abs(c(D$alt2 - D$ch1, D$alt3 - D$ch1)))))
cat(sprintf("sessions agree to within %.2f pp at every depth\n", max(abs(D$ch1 - D$ch2))))

L <- c("\\begin{table}[htbp!]", "\\centering",
"\\caption{Change in semivariance from 6 to 9\\,mm within shells of Euclidean",
"depth from the brain boundary, positive indicating that the curve falls. Shells",
"are labelled by their inner edge. Each session is estimated independently from",
"its own five runs, so the two halves are a replication rather than a split of",
"one estimate. The change is formed within each subject-run, averaged over that",
"subject's runs, and summarised across subjects; intervals are percentile",
"bootstrap over subjects with 10{,}000 draws. The declining column counts",
"subjects in which the curve falls between 6 and 9\\,mm; 21 of 21 corresponds to",
"$p = 9.5 \\times 10^{-7}$ on a two-sided exact binomial test, which is the",
"smallest value attainable from 21 observations.}",
"\\label{tab:depth}", "\\begin{tabular}{lrrcrrc}", "\\toprule",
"& \\multicolumn{3}{c}{Session 1} & \\multicolumn{3}{c}{Session 2} \\\\",
"\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
"Depth & Change & 95\\% CI & Declining & Change & 95\\% CI & Declining \\\\",
"\\midrule")
for (i in seq_len(nrow(D)))
  L <- c(L, sprintf(paste0("%s & $%+.2f\\%%$ & $[%+.2f, %+.2f]$ & %d of %d & ",
                           "$%+.2f\\%%$ & $[%+.2f, %+.2f]$ & %d of %d \\\\"),
    D$label[i], D$ch1[i], D$lo1[i], D$hi1[i], D$nd1[i], D$n1[i],
    D$ch2[i], D$lo2[i], D$hi2[i], D$nd2[i], D$n2[i]))
L <- c(L, "\\bottomrule",
"\\multicolumn{7}{l}{\\footnotesize Forming the ratio per axis before averaging, or weighting axes by their}\\\\",
sprintf(paste0("\\multicolumn{7}{l}{\\footnotesize valid base counts, changes no cell by more than %.2f",
               " percentage points.}"), max(abs(c(D$alt2 - D$ch1, D$alt3 - D$ch1)))),
"\\end{tabular}", "\\end{table}")
writeLines(L, file.path(TAB, "depth.tex"), useBytes = TRUE)
cat("wrote depth.tex\n")
})

