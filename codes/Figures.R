# Every figure in the manuscript.
#
# Reads the per-run lag files and the intermediate result tables, and writes
# the figures to the output directory. Both locations are set in the preamble
# below and need no editing.

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

# Variograms, zones, anisotropy, fits and residuals
local({
src_funs <- function(path) {
  e <- parse(path)
  for (x in e)
    if (is.call(x) && length(x) >= 3 &&
        as.character(x[[1]])[1] %in% c("<-", "=") &&
        is.call(x[[3]]) && identical(as.character(x[[3]][[1]])[1], "function"))
      eval(x, envir = globalenv())
}
src_funs(file.path(REL, "Main analysis_Model Fitting.R"))

FACTS <- file.path(OUT, "F1_figure_facts.txt")
say <- function(...) { m <- sprintf(...); cat(m); cat(m, file = FACTS, append = TRUE) }
cat("", file = FACTS)
say("figure facts, regenerated %s\n", format(Sys.Date()))
say("every number a caption quotes is printed here so the two cannot drift\n\n")

# the corrected curves
L    <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
cur  <- L[, .(g = mean(gamma)), by = .(subject, distance_mm)][order(subject, distance_mm)]
subs <- sort(unique(cur$subject))
dgrid <- sort(unique(cur$distance_mm))
NS <- length(subs); ND <- length(dgrid)
say("basis: %d subjects, %d runs each over %d sessions, %d exact distance classes\n",
    NS, uniqueN(L$run), uniqueN(L$session), ND)
say("distance range %.1f to %.1f mm\n\n", min(dgrid), max(dgrid))

# zones, as defined in the manuscript
Z <- list(short   = c(2.9,  3.1),
          corner  = c(4.2,  6.2),
          plateau = c(7.0,  11.5),
          long    = c(11.5, 21.0))

# fit all four methods per subject
gvar <- function(d, nugget, sill, range) nugget + (sill-nugget)*(1-exp(-(d/range)^2))
fit_gauss_single <- function(d, y)
  tryCatch(nlsLM(y ~ gvar(d, nugget, sill, range),
    start = list(nugget = 0, sill = max(y)*0.95, range = 4),
    lower = c(nugget = 0, sill = max(y)*0.05, range = 0.1),
    upper = c(nugget = max(y), sill = max(y)*5, range = 100),
    control = nls.lm.control(maxiter = 1000)), error = function(e) NULL)

FIT <- list(); RSSt <- list()
cat("fitting 21 subjects ...\n")
for (s in subs) {
  x <- cur[subject == s]; d <- x$distance_mm; y <- x$g
  tn <- (d - min(d)) / (max(d) - min(d))
  ncs <- fit_ncs(tn, y)
  kk  <- (1:10)[which.min(search_fn(tn, y, 1:10, ptr_wrap)[, 2])]
  pf  <- ptr(tn, y, kk)[[1]][, 1]
  bwm <- seq(0.05, 0.60, by = 0.025)
  bw  <- bwm[which.min(search_fn(tn, y, bwm, locpoly3_wrap)[, 2])]
  lpf <- locpoly3(tn, y, bw, tn)[[1]][, 1]
  gf  <- fit_gauss_single(d, y)
  gfit <- if (is.null(gf)) rep(NA_real_, length(y)) else fitted(gf)
  FIT[[as.character(s)]] <- list(d = d, y = y, tn = tn, ncs = ncs,
    NCS = ncs$fitted, PTR = pf, LP = lpf, Gaussian = gfit,
    k = kk, bw = bw, edf = ncs$edf, lambda = ncs$lambda,
    grange = if (is.null(gf)) NA else coef(gf)["range"])
  RSSt[[length(RSSt)+1]] <- data.table(subject = s,
    NCS = sum((y-ncs$fitted)^2), PTR = sum((y-pf)^2),
    LP  = sum((y-lpf)^2),        Gaussian = sum((y-gfit)^2))
}
RS <- rbindlist(RSSt)
MTH <- c("NCS","PTR","LP","Gaussian")
COL <- c(NCS = "#1b7837", PTR = "#2166ac", LP = "#b2182b", Gaussian = "#7b3294")

# residual matrices, distance x subject
resid_mat <- function(m) sapply(subs, function(s) FIT[[as.character(s)]]$y -
                                                  FIT[[as.character(s)]][[m]])
RESID <- lapply(MTH, resid_mat); names(RESID) <- MTH

# Figure 1 - normalised empirical variograms, all 21 subjects.
# Each curve is scaled between its own minimum and maximum, and the panel is
# built around the local dip between 6 mm and about 12 mm.
Y <- sapply(subs, function(s) FIT[[as.character(s)]]$y)
Yn <- apply(Y, 2, function(v) (v - min(v)) / (max(v) - min(v)))
i6 <- which.min(abs(dgrid - 6)); i9 <- which.min(abs(dgrid - 9))
dip <- 100*(Y[i6, ] - apply(Y[dgrid > 6.5 & dgrid < 12, , drop = FALSE], 2, min))/Y[i6, ]
ch69 <- 100*(Y[i6, ] - Y[i9, ])/Y[i6, ]
say("FIGURE 1  global maximum at the final separation in %d of %d subjects\n",
    sum(apply(Y, 2, which.max) == ND), NS)
say("FIGURE 1  dip below the 6 mm value present in %d of %d, mean depth %.2f%%\n",
    sum(dip > 0), NS, mean(dip))
say("FIGURE 1  6 to 9 mm fall on the RUN-AVERAGED curve: %d of %d, mean %.2f%%\n",
    sum(ch69 > 0), NS, mean(ch69))
say("          the manuscript quotes the per-run convention instead, %s\n",
    "+4.32% and 21 of 21; see F5_basisA_facts.txt. Do not mix them in a caption.\n")

pdf(file.path(FIG, "figure01_variograms.pdf"), width = 6.5, height = 4.6)
par(mar = c(4.2, 4.2, 1.2, 1.2), las = 1)
matplot(dgrid, Yn, type = "n", xlab = "Separation (mm)",
        ylab = "Normalised semivariance", ylim = c(-0.02, 1.02))
rect(6, -0.05, 12, 1.05, col = "#eef3f8", border = NA); box()
for (j in seq_len(NS)) lines(dgrid, Yn[, j], col = adjustcolor("#2166ac", 0.5), lwd = 1)
lines(dgrid, rowMeans(Yn), col = "#b2182b", lwd = 2.6)
points(dgrid[c(i6, i9)], rowMeans(Yn)[c(i6, i9)], pch = 16, col = "#b2182b", cex = 0.9)
text(9, 0.06, "local dip", cex = 0.75, col = "grey25")
legend("bottomright", bty = "n", lwd = c(2.6, 1), cex = 0.8,
       col = c("#b2182b", adjustcolor("#2166ac", 0.6)),
       legend = c("cohort mean", sprintf("individual subjects (%d)", NS)))
dev.off()

# Figure 2 - cohort mean Gaussian residual, with the four zones
gm <- rowMeans(RESID$Gaussian); gse <- apply(RESID$Gaussian, 1, sd)/sqrt(NS)
zc <- c(short = "#f0f0f0", corner = "#d6e6f4", plateau = "#fbe3d6", long = "#f0f0f0")
pdf(file.path(FIG, "figure02_four_zones.pdf"), width = 6.0, height = 4.2)
par(mar = c(4.2, 4.6, 1.2, 1.2), las = 1)
plot(dgrid, gm, type = "n", xlab = "Separation (mm)",
     ylab = "Cohort mean Gaussian residual",
     ylim = range(gm - gse, gm + gse) * 1.15)
for (nm in names(Z))
  rect(Z[[nm]][1], par("usr")[3], Z[[nm]][2], par("usr")[4],
       col = zc[[nm]], border = NA)
box()
abline(h = 0, col = "grey40", lty = 2)
polygon(c(dgrid, rev(dgrid)), c(gm - gse, rev(gm + gse)),
        col = adjustcolor("#7b3294", 0.22), border = NA)
lines(dgrid, gm, col = "#7b3294", lwd = 2); points(dgrid, gm, pch = 16, cex = 0.6, col = "#7b3294")
text(c(3.0, 5.2, 9.2, 16), par("usr")[4]*0.9, cex = 0.72, col = "grey25",
     labels = c("short", "corner", "plateau", "long distance"))
dev.off()
say("FIGURE 2  corner zone mean residual %+.4f, plateau %+.4f (cohort means)\n",
    mean(gm[dgrid >= Z$corner[1] & dgrid <= Z$corner[2]]),
    mean(gm[dgrid >= Z$plateau[1] & dgrid <= Z$plateau[2]]))

# Figures 3 and 4 - directional variograms and the disparity ratio
# The submitted versions showed three subjects. The text now reports all 21, so
# these do too. Axis-aligned lags exist at 3, 6, 9 and 12 mm on each array axis.
ax <- L[(dx != 0) + (dy != 0) + (dz != 0) == 1]
ax[, axis := ifelse(dx != 0, "axis 1", ifelse(dy != 0, "axis 2", "axis 3"))]
AX <- ax[, .(g = mean(gamma)), by = .(subject, axis, distance_mm)]
OM <- cur[distance_mm %in% sort(unique(AX$distance_mm))]
dax <- sort(unique(AX$distance_mm))
acol <- c("axis 1" = "#b2182b", "axis 2" = "#1b7837", "axis 3" = "#2166ac")

pdf(file.path(FIG, "figure03_directional.pdf"), width = 9.0, height = 7.2)
par(mfrow = c(4, 6), mar = c(2.5, 2.9, 1.6, 0.6), mgp = c(1.7, 0.5, 0), las = 1,
    cex.axis = 0.65, cex.lab = 0.7, cex.main = 0.8)
for (s in subs) {
  a <- AX[subject == s]; o <- OM[subject == s]
  plot(dax, o$g, type = "n", ylim = range(a$g, o$g), main = sprintf("sub-%02d", s),
       xlab = "mm", ylab = expression(hat(gamma)))
  for (k in names(acol)) lines(dax, a[axis == k][order(distance_mm)]$g, col = acol[k], lwd = 1.4)
  lines(dax, o$g, col = "grey35", lwd = 1.4, lty = 2)
}
plot.new()
legend("center", bty = "n", cex = 0.8, lwd = 1.4, lty = c(1,1,1,2),
       col = c(acol, "grey35"), legend = c(names(acol), "omnidirectional"))
plot.new(); plot.new()
dev.off()

DR <- AX[, .(ratio = max(g)/min(g)), by = .(subject, distance_mm)][order(subject, distance_mm)]
DM <- DR[, .(m = mean(ratio), lo = min(ratio), hi = max(ratio)), by = distance_mm]
say("FIGURE 4  disparity ratio cohort mean by separation: %s\n",
    paste(sprintf("%.0fmm %.3f", DM$distance_mm, DM$m), collapse = "  "))
say("FIGURE 4  ratio spans %.3f to %.3f across all subject-separations\n",
    min(DR$ratio), max(DR$ratio))

pdf(file.path(FIG, "figure04_disparity.pdf"), width = 6.5, height = 4.4)
par(mar = c(4.2, 4.4, 1.2, 1.2), las = 1)
plot(DR$distance_mm, DR$ratio, type = "n", xlab = "Separation (mm)",
     ylab = "Directional disparity ratio", ylim = c(1, max(DR$ratio)*1.02))
abline(h = 1, col = "grey45", lty = 2)
for (s in subs) lines(DR[subject == s]$distance_mm, DR[subject == s]$ratio,
                      col = adjustcolor("grey40", 0.45), lwd = 0.8)
lines(DM$distance_mm, DM$m, col = "#b2182b", lwd = 2.4)
points(DM$distance_mm, DM$m, col = "#b2182b", pch = 16)
legend("topright", bty = "n", cex = 0.8, lwd = c(2.4, 0.8),
       col = c("#b2182b", adjustcolor("grey40", 0.6)),
       legend = c("cohort mean", "individual subjects"))
dev.off()

# Figure 5 - residual sum of squares per subject, by method
say("FIGURE 5  mean RSS  NCS %.4f  PTR %.4f  LP %.4f  Gaussian %.4f\n",
    mean(RS$NCS), mean(RS$PTR), mean(RS$LP), mean(RS$Gaussian))
say("FIGURE 5  Gaussian is the worst of the four in %d of %d subjects\n",
    sum(RS$Gaussian > pmax(RS$NCS, RS$PTR, RS$LP)), NS)
pdf(file.path(FIG, "figure05_rss.pdf"), width = 6.5, height = 4.2)
par(mar = c(4.2, 4.4, 1.2, 1.2), las = 1)
plot(seq_len(NS), RS$Gaussian, type = "n", xaxt = "n", log = "y",
     xlab = "Subject", ylab = "Residual sum of squares",
     ylim = range(unlist(RS[, ..MTH])))
axis(1, at = seq_len(NS), labels = subs, cex.axis = 0.6)
for (m in MTH) { lines(seq_len(NS), RS[[m]], col = COL[m], lwd = 1.6)
                 points(seq_len(NS), RS[[m]], col = COL[m], pch = 16, cex = 0.7) }
legend("topleft", bty = "n", cex = 0.8, horiz = TRUE, lwd = 1.6,
       col = COL[MTH], legend = MTH)
dev.off()

# Figure 6 - cohort mean residual by distance, all four methods
pdf(file.path(FIG, "figure06_residual_comparison.pdf"), width = 6.5, height = 4.4)
par(mar = c(4.2, 4.4, 1.2, 1.2), las = 1)
allm <- sapply(MTH, function(m) rowMeans(RESID[[m]]))
alls <- sapply(MTH, function(m) apply(RESID[[m]], 1, sd))
plot(dgrid, allm[, 1], type = "n", xlab = "Separation (mm)",
     ylab = "Cohort mean residual", ylim = range(allm - alls, allm + alls))
abline(h = 0, col = "grey45", lty = 2)
for (m in MTH) {
  polygon(c(dgrid, rev(dgrid)),
          c(allm[, m] - alls[, m], rev(allm[, m] + alls[, m])),
          col = adjustcolor(COL[m], 0.16), border = NA)
  lines(dgrid, allm[, m], col = COL[m], lwd = 2)
}
legend("bottomright", bty = "n", cex = 0.8, lwd = 2, col = COL[MTH], legend = MTH)
dev.off()
for (m in MTH)
  say("FIGURE 6  %-8s cohort mean |residual| %.4f, max %.4f\n",
      m, mean(abs(rowMeans(RESID[[m]]))), max(abs(rowMeans(RESID[[m]]))))

# Figure 9 - NCS effective degrees of freedom per subject, with the subjects
# whose smoothing parameter sits at the floor of the search grid drawn
# separately.
edf <- sapply(subs, function(s) FIT[[as.character(s)]]$edf)
atfl <- abs(log10(sapply(subs, function(s) FIT[[as.character(s)]]$lambda)) + 6) < 1e-9
say("FIGURE 9  NCS edf mean %.2f, range %.2f to %.2f; n = %d points\n",
    mean(edf), min(edf), max(edf), ND)
say("FIGURE 9  lambda at the search floor in %d of %d subjects, all at edf %.2f\n",
    sum(atfl), NS, max(edf))
pdf(file.path(FIG, "figure09_ncs_edf.pdf"), width = 6.4, height = 3.9)
par(mar = c(4.2, 4.4, 1.2, 1.2), las = 1)
plot(seq_len(NS), edf, type = "n", xaxt = "n", ylim = c(0, ND*1.05),
     xlab = "Subject", ylab = "NCS effective degrees of freedom")
abline(h = ND, col = "#b2182b", lty = 2)
points(which(!atfl), edf[!atfl], pch = 16, col = "#1b7837", cex = 1.1)
points(which(atfl),  edf[atfl],  pch = 17, col = "#b2182b", cex = 1.1)
axis(1, at = seq_len(NS), labels = subs, cex.axis = 0.6)
text(1, ND, sprintf("interpolation, df = %d", ND), pos = 3, cex = 0.7, col = "#b2182b", adj = 0)
legend("bottomleft", bty = "n", cex = 0.75, pch = c(17, 16),
       col = c("#b2182b", "#1b7837"),
       legend = c(sprintf("smoothing parameter at the grid floor (%d)", sum(atfl)),
                  sprintf("interior optimum (%d)", sum(!atfl))))
dev.off()

# Figures 10-13 - per-subject fits, one figure per method
panel_fits <- function(m, fname) {
  pdf(file.path(FIG, fname), width = 9.0, height = 7.2)
  par(mfrow = c(4, 6), mar = c(2.5, 2.9, 1.6, 0.6), mgp = c(1.7, 0.5, 0), las = 1,
      cex.axis = 0.65, cex.lab = 0.7, cex.main = 0.8)
  for (s in subs) {
    f <- FIT[[as.character(s)]]
    plot(f$d, f$y, pch = 16, cex = 0.5, col = "grey45", main = sprintf("sub-%02d", s),
         xlab = "mm", ylab = expression(hat(gamma)))
    lines(f$d, f[[m]], col = COL[m], lwd = 1.8)
  }
  plot.new(); legend("center", bty = "n", cex = 0.85,
                     lwd = c(NA, 1.8), pch = c(16, NA), pt.cex = c(0.6, NA),
                     col = c("grey45", COL[m]), legend = c("observed", m))
  plot.new(); plot.new()
  dev.off()
}
panel_fits("NCS",      "figure10_ncs_fits.pdf")
panel_fits("PTR",      "figure11_ptr_fits.pdf")
panel_fits("LP",       "figure12_lp_fits.pdf")
panel_fits("Gaussian", "figure13_gaussian_fits.pdf")

# Figure 14 - NCS generalized cross-validation curves
pdf(file.path(FIG, "figure14_gcv_curves.pdf"), width = 9.0, height = 7.2)
par(mfrow = c(4, 6), mar = c(2.5, 2.9, 1.6, 0.6), mgp = c(1.7, 0.5, 0), las = 1,
    cex.axis = 0.65, cex.lab = 0.7, cex.main = 0.8)
for (s in subs) {
  n <- FIT[[as.character(s)]]$ncs
  fl <- abs(log10(n$lambda) + 6) < 1e-9
  plot(n$coarse_log, n$gcv_c, type = "b", pch = 16, cex = 0.45, col = "#2166ac", log = "y",
       main = sprintf("sub-%02d%s", s, if (fl) " *" else ""),
       xlab = expression(log[10]~lambda), ylab = "GCV")
  points(n$fine_log, n$gcv_f, pch = 16, cex = 0.45, col = "#1b7837")
  abline(v = log10(n$lambda), col = if (fl) "#b2182b" else "grey35", lty = 2, lwd = if (fl) 1.6 else 1)
}
plot.new(); legend("center", bty = "n", cex = 0.75, pch = c(16,16,NA,NA), lty = c(NA,NA,2,2),
                   col = c("#2166ac", "#1b7837", "#b2182b", "grey35"),
                   legend = c("coarse grid", "fine grid",
                              "selected, at the grid floor (*)", "selected, interior"))
plot.new(); plot.new()
dev.off()
lam <- log10(sapply(subs, function(s) FIT[[as.character(s)]]$lambda))
say("FIGURE 14 selected log10 lambda spans %.2f to %.2f, mean %.2f\n",
    min(lam), max(lam), mean(lam))
say("FIGURE 14 selection on a grid boundary in %d of %d subjects\n",
    sum(sapply(subs, function(s) FIT[[as.character(s)]]$ncs$boundary)), NS)

# Figures 15-18 - per-subject residual profiles, one figure per method
rng_np <- range(sapply(c("NCS","PTR","LP"), function(m) range(RESID[[m]])))
panel_resid <- function(m, fname, ylim) {
  pdf(file.path(FIG, fname), width = 9.0, height = 7.2)
  par(mfrow = c(4, 6), mar = c(2.5, 2.9, 1.6, 0.6), mgp = c(1.7, 0.5, 0), las = 1,
      cex.axis = 0.65, cex.lab = 0.7, cex.main = 0.8)
  for (i in seq_along(subs)) {
    plot(dgrid, RESID[[m]][, i], type = "b", pch = 16, cex = 0.45, col = COL[m],
         ylim = ylim, main = sprintf("sub-%02d", subs[i]),
         xlab = "mm", ylab = "residual")
    abline(h = 0, col = "grey55", lty = 2)
  }
  plot.new(); legend("center", bty = "n", cex = 0.85, lwd = 1.4, col = COL[m],
                     legend = sprintf("%s residuals", m))
  plot.new(); plot.new()
  dev.off()
}
panel_resid("NCS",      "figure15_ncs_residuals.pdf",      rng_np)
panel_resid("PTR",      "figure16_ptr_residuals.pdf",      rng_np)
panel_resid("LP",       "figure17_lp_residuals.pdf",       rng_np)
panel_resid("Gaussian", "figure18_gaussian_residuals.pdf", range(RESID$Gaussian))
say("FIGURE 15-17 non-parametric residuals share the axis %.3f to %.3f\n",
    rng_np[1], rng_np[2])
say("FIGURE 18 Gaussian residuals need a wider axis, %.3f to %.3f\n",
    min(RESID$Gaussian), max(RESID$Gaussian))

fwrite(RS, file.path(OUT, "F1_rss_per_subject.csv"))
fwrite(DR, file.path(OUT, "F1_disparity.csv"))
say("\nwrote %d figures -> %s\n", length(list.files(FIG, pattern = "\\.pdf$")), FIG)
})

# Interpolation signature
local({
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, `:=`(k2 = as.integer(round((distance_mm/VOX)^2)),
         prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L],
         naxis = (dx != 0) + (dy != 0) + (dz != 0))]
V <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2, prod, naxis)][
       , .(g = mean(g)), by = .(dx, dy, dz, k2, prod, naxis)]
K <- V[, .(g = mean(g), prod = mean(prod), axis = all(naxis == 1)), by = k2][
       , d := VOX*sqrt(k2)][order(d)]

# Both panels are built from one joint fit, so the vertical offsets drawn in (a)
# are the same quantity plotted in (b).
fit <- lm(g ~ ns(d, df = 6) + prod, K)
cf  <- summary(fit)$coefficients
b   <- cf["prod", 1]
K[, sm := fitted(fit) - b*prod]     # the smooth component alone
K[, fine := g - sm]                 # partial residual for the kernel term

COL_AX <- "#B2182B"; COL_OB <- "grey35"

pdf(file.path(FIG, "figure07_resampling.pdf"), width = 9.0, height = 3.9)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.0, 1.0), las = 1, mgp = c(2.5, 0.7, 0))

# (a) the curve, with the axis-only classes marked
plot(K$d, K$g, type = "n", xlab = "Separation (mm)",
     ylab = "Semivariance", main = "(a) Cohort curve and its smooth", cex.main = 1.0)
dd <- seq(min(K$d), max(K$d), length.out = 300)
lines(K$d, K$sm, col = "grey60", lwd = 2)   # smooth component of the joint fit
segments(K$d, K$sm, K$d, K$g, col = "grey75", lwd = 1)
points(K$d[!K$axis], K$g[!K$axis], pch = 16, cex = 0.95, col = COL_OB)
points(K$d[K$axis],  K$g[K$axis],  pch = 17, cex = 1.35, col = COL_AX)
i3 <- which(abs(K$d - 3) < 1e-9); i6 <- which(abs(K$d - 6) < 1e-9)
text(K$d[i3] + 0.9, K$g[i3] - 0.02, "3 mm", col = COL_AX, cex = 0.8, adj = 0)
text(K$d[i6] + 0.5, K$g[i6] + 0.07, "6 mm", col = COL_AX, cex = 0.8, adj = 0)
legend("bottomright", bty = "n", cex = 0.8,
       legend = c("axis displacements only", "classes with oblique members", "smooth, 6 df"),
       pch = c(17, 16, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
       col = c(COL_AX, COL_OB, "grey60"))

# (b) the departure against the kernel term
plot(K$prod, K$fine, type = "n", xlab = expression(a(group("|",delta[x],"|")) *
       a(group("|",delta[y],"|")) * a(group("|",delta[z],"|"))),
     ylab = "Departure from the smooth",
     main = "(b) Departure against the kernel term", cex.main = 1.0)
abline(h = 0, col = "grey80"); abline(a = 0, b = b, col = "grey35", lwd = 2, lty = 2)
points(K$prod[!K$axis], K$fine[!K$axis], pch = 16, cex = 0.95, col = COL_OB)
points(K$prod[K$axis],  K$fine[K$axis],  pch = 17, cex = 1.35, col = COL_AX)
text(K$prod[i3], K$fine[i3], " 3 mm", col = COL_AX, cex = 0.8, adj = c(1.05, 0.4))
text(K$prod[i6], K$fine[i6], "6 mm ", col = COL_AX, cex = 0.8, adj = c(-0.05, 0.4))
legend("topright", bty = "n", cex = 0.8,
       legend = c(sprintf("slope %.2f", cf["prod",1]), sprintf("t = %.2f", cf["prod",3])))
dev.off()

cat(sprintf("kernel term: slope %.3f, se %.3f, t %.2f, p %.2e\n",
            cf["prod",1], cf["prod",2], cf["prod",3], cf["prod",4]))
cat(sprintf("axis-only classes: %s mm\n", paste(sprintf("%.1f", K$d[K$axis]), collapse = ", ")))
cat(sprintf("wrote %s\n", file.path(FIG, "figure07_resampling.pdf")))
})

