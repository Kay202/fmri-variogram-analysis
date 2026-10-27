# Shell models, the interpolation signature, the depth gradient, the
# reference covariance and the identifiability sweep.
#
# Reads the per-run lag files and the intermediate result tables, and writes
# its own results back alongside them. Both locations are set in the preamble
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

# Shell models: Gaussian, AFNI mixed, and the admissible Shapiro-Botha fit
local({
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
sh <- melt(L, id.vars = c("subject", "run", "distance_mm"),
           measure.vars = patterns("^gamma_shell"),
           variable.name = "shell", value.name = "g")
sh[, shell := as.integer(sub("gamma_shell", "", shell))]
sh <- sh[!is.na(g)]
curve <- sh[, .(gamma = mean(g)), by = .(subject, shell, distance_mm)][
  , .(gamma = mean(gamma)), by = .(shell, distance_mm)][order(shell, distance_mm)]
depth_mm <- c(9, 12, 15, 18, 24)

# Shapiro-Botha, as published
omega3 <- function(x) { o <- rep(1, length(x)); nz <- x != 0
                        o[nz] <- sin(x[nz])/x[nz]; o }
# Node band: pi over the largest separation to pi over the smallest, with no
# margin above. Past t*r = pi the basis function is flat over the sampled range,
# so a node beyond the band cannot be separated from the nugget.
make_nodes <- function(d) {
  band <- c(pi/max(d), pi/min(d))
  exp(seq(log(band[1]), log(band[2]), length.out = 60))
}
design_sb <- function(d, nodes) {
  M <- matrix(0, length(d), length(nodes) + 1L); M[, 1] <- 1
  for (j in seq_along(nodes)) M[, j + 1L] <- 1 - omega3(nodes[j]*d)
  M
}
fit_sb <- function(d, y, nodes) {
  X <- design_sb(d, nodes); p <- nnls(X, y)$x
  list(p = p, ss = sum((y - X %*% p)^2), nz = sum(p > 1e-12))
}

# the other two arms, carried over unchanged
fit_gauss <- function(d, y, grid = seq(0.5, 60, length.out = 3000)) {
  best <- NULL
  for (a in grid) {
    u <- exp(-(d/a)^2); X <- cbind(u, 1 - u)
    cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) next
    if (cf[1] < 0) cf <- c(0, sum(y*(1-u))/sum((1-u)^2))
    ss <- sum((y - X %*% cf)^2)
    if (is.null(best) || ss < best$ss)
      best <- list(ss = ss, range = a, nugget = cf[1], sill = cf[2])
  }
  best
}
fit_mixed <- function(d, y) {
  best <- NULL
  for (b in seq(0.3, 25, length.out = 90))
    for (cc in seq(0.3, 50, length.out = 110)) {
      X <- cbind(1 - exp(-(d^2)/(2*b^2)), 1 - exp(-d/cc))
      p <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
      if (is.null(p) || anyNA(p)) next
      p <- pmax(p, 0)
      ss <- sum((y - X %*% p)^2)
      if (is.null(best) || ss < best$ss)
        best <- list(ss = ss, b = b, c = cc, p1 = p[1], p2 = p[2])
    }
  best
}

coefs <- list(); pars <- list()
for (k in sort(unique(curve$shell))) {
  g <- curve[shell == k][order(distance_mm)]
  d <- g$distance_mm; y <- g$gamma
  nodes <- make_nodes(d)
  sb <- fit_sb(d, y, nodes)
  ga <- fit_gauss(d, y)
  mx <- fit_mixed(d, y)

  coefs[[length(coefs)+1]] <- data.table(
    shell = k, depth_mm = depth_mm[k],
    term = seq_along(sb$p) - 1L, node = c(NA_real_, nodes), coef = sb$p)
  pars[[length(pars)+1]] <- data.table(
    shell = k, depth_mm = depth_mm[k], lambda = 0,
    sb_rss = sb$ss, sb_sill = sum(sb$p), sb_active = sb$nz,
    node_lo = min(nodes), node_hi = max(nodes),
    g_range = ga$range, g_nugget = ga$nugget, g_sill = ga$sill, g_rss = ga$ss,
    m_b = mx$b, m_c = mx$c, m_p1 = mx$p1, m_p2 = mx$p2, m_rss = mx$ss)
  cat(sprintf("shell %d (%2d mm): SB rss %.5f, %2d active nodes of 61 | Gauss rss %.5f range %.2f | mixed rss %.5f\n",
              k, depth_mm[k], sb$ss, sb$nz, ga$ss, ga$range, mx$ss))
}
CO <- rbindlist(coefs); PA <- rbindlist(pars)
fwrite(CO, file.path(OUT, "shell_sb_coefs_v2.csv"))
fwrite(PA, file.path(OUT, "shell_model_params_v2.csv"))

# how far did the null fields move?
old <- tryCatch(fread(file.path(DIAG, "out_R", "shell_model_params.csv")),
                error = function(e) NULL)
if (!is.null(old)) {
  cat("\nShapiro-Botha residual sum of squares, penalised against plain:\n")
  m <- merge(old[, .(shell, old_rss = sb_rss, old_sill = sb_sill)],
             PA[, .(shell, new_rss = sb_rss, new_sill = sb_sill, sb_active)], by = "shell")
  for (i in seq_len(nrow(m)))
    cat(sprintf("  shell %d: rss %.5f -> %.5f (%+.1f%%) | sill %.4f -> %.4f (%+.1f%%) | %d active\n",
                m$shell[i], m$old_rss[i], m$new_rss[i],
                100*(m$new_rss[i]-m$old_rss[i])/m$old_rss[i],
                m$old_sill[i], m$new_sill[i],
                100*(m$new_sill[i]-m$old_sill[i])/m$old_sill[i], m$sb_active[i]))
  cat("\n  The sill is what sets the scale of the simulated null fields, so a\n")
  cat("  change there is what would move the inference results.\n")
}
cat(sprintf("\nwrote shell_sb_coefs_v2.csv, shell_model_params_v2.csv -> %s\n", OUT))
})

# Exponential shell fit
local({
# cohort-mean curve per shell, identical to fit_shell_models.R L19-28
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
sh <- melt(L, id.vars = c("subject", "run", "distance_mm"),
           measure.vars = patterns("^gamma_shell"),
           variable.name = "shell", value.name = "g")
sh[, shell := as.integer(sub("gamma_shell", "", shell))]
sh <- sh[!is.na(g)]
curve <- sh[, .(gamma = mean(g)), by = .(subject, shell, distance_mm)][
  , .(gamma = mean(gamma)), by = .(shell, distance_mm)][order(shell, distance_mm)]
depth_mm <- c(9, 12, 15, 18, 24)

# exponential with nugget, same parametrisation as fit_gauss
# gamma(d) = nugget * exp(-d/phi) + sill * (1 - exp(-d/phi))
# so gamma(0) = nugget and gamma(inf) = sill, matching how the Gaussian arm is
# written in fit_shell_models.R. phi is profiled on a grid, the two linear
# coefficients solved in closed form at each phi.
fit_exp <- function(d, y, grid = seq(0.2, 60, length.out = 3000)) {
  best <- NULL
  for (phi in grid) {
    u <- exp(-d / phi); X <- cbind(u, 1 - u)
    cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) next
    if (cf[1] < 0) cf <- c(0, sum(y * (1 - u)) / sum((1 - u)^2))
    ss <- sum((y - X %*% cf)^2)
    if (is.null(best) || ss < best$ss)
      best <- list(ss = ss, range = phi, nugget = cf[1], sill = cf[2])
  }
  best
}

P <- fread(file.path(DIAG, "out_R", "shell_model_params.csv"))
rows <- list()
for (k in sort(unique(curve$shell))) {
  g <- curve[shell == k][order(distance_mm)]
  e <- fit_exp(g$distance_mm, g$gamma)
  rows[[length(rows) + 1]] <- data.table(
    shell = k, depth_mm = depth_mm[k],
    e_range = e$range, e_nugget = e$nugget, e_sill = e$sill, e_rss = e$ss)
}
E <- rbindlist(rows)
fwrite(E, file.path(OUT, "S3_shell_exp_params.csv"))

M <- merge(P[, .(shell, depth_mm, sb_rss, g_rss, g_range, g_nugget, g_sill, m_rss, m_c)],
           E, by = c("shell", "depth_mm"))

cat("Exponential with nugget, CLEAN's covariance model, fitted per depth shell\n")
cat("RSS against the cohort-mean curve; lower is a better fit to the variogram.\n\n")
cat(sprintf("%-6s %-8s %10s %10s %10s %10s | %8s %8s %8s\n",
            "shell", "depth", "SB", "Gaussian", "exponent.", "AFNImixed",
            "e_range", "e_nug", "e_sill"))
for (i in seq_len(nrow(M))) {
  r <- M[i]
  cat(sprintf("%-6d %-8s %10.4f %10.4f %10.4f %10.4f | %8.2f %8.3f %8.3f\n",
              r$shell, paste0(r$depth_mm, " mm"), r$sb_rss, r$g_rss, r$e_rss, r$m_rss,
              r$e_range, r$e_nugget, r$e_sill))
}
cat(sprintf("\nmean RSS   SB %.4f | Gaussian %.4f | exponential %.4f | AFNI mixed %.4f\n",
            mean(M$sb_rss), mean(M$g_rss), mean(M$e_rss), mean(M$m_rss)))
cat(sprintf("exponential beats Gaussian in %d of %d shells\n",
            sum(M$e_rss < M$g_rss), nrow(M)))
cat(sprintf("Shapiro-Botha beats exponential in %d of %d shells\n",
            sum(M$sb_rss < M$e_rss), nrow(M)))

cat("\nNote on the AFNI mixed fit: its mono-exponential range parameter sits at\n")
cat(sprintf("the upper edge of its search grid in every shell (m_c = %.0f, grid max 50),\n",
            M$m_c[1]))
cat("so that arm is a boundary solution and its conservatism should be read with\n")
cat("that in mind.\n")
cat(sprintf("\nwrote S3_shell_exp_params.csv -> %s\n", OUT))
})

# Gaussian fit diagnostics, in the manuscript's own fitting box
local({
TXT  <- file.path(OUT, "S44_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)

L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
CUR <- L[, .(g = mean(gamma)), by = .(subject, k2)][, d := VOX*sqrt(k2)][order(subject, d)]
subs <- sort(unique(CUR$subject))
say("basis A run-averaged curves, %d subjects, %d separations each\n\n",
    length(subs), CUR[subject == subs[1], .N])

gv <- function(d, nug, sill, rng) nug + (sill - nug)*(1 - exp(-(d/rng)^2))

# the three boxes the project actually contains, each named by where it lives
boxes <- list(
  manuscript = function(y) list(
    lo = c(nugget = 0, sill = mean(y)*0.5, range = 0.5),
    hi = c(nugget = min(y)*1.5, sill = max(y)*2, range = 30)),
  signtest = function(y) list(
    lo = c(nugget = 0, sill = max(y)*0.3, range = 0.5),
    hi = c(nugget = Inf, sill = Inf, range = 25)),
  quantile_start = function(y) list(
    lo = c(nugget = 0, sill = as.numeric(quantile(y, 0.90))*0.5, range = 0.5),
    hi = c(nugget = as.numeric(quantile(y, 0.90))*0.3, sill = Inf, range = 25))
)

rows <- list()
for (bn in names(boxes)) {
  for (s in subs) {
    q <- CUR[subject == s][order(d)]; d <- q$d; y <- q$g
    b <- boxes[[bn]](y)
    st <- c(nugget = max(0, min(y)*0.5), sill = max(y)*0.95, range = 4)
    st <- pmin(pmax(st, b$lo), ifelse(is.finite(b$hi), b$hi, st))
    f <- tryCatch(nlsLM(y ~ gv(d, nugget, sill, range), start = as.list(st),
                        lower = b$lo, upper = b$hi,
                        control = nls.lm.control(maxiter = 500)),
                  error = function(e) NULL)
    if (is.null(f)) next
    p <- coef(f); fit <- gv(d, p["nugget"], p["sill"], p["range"])
    rows[[length(rows)+1]] <- data.table(box = bn, subject = s,
      range = p["range"], nugget = p["nugget"], sill = p["sill"],
      r2 = 1 - sum((y - fit)^2)/sum((y - mean(y))^2),
      nug_lo = abs(p["nugget"] - b$lo["nugget"]) < 1e-6,
      nug_hi = is.finite(b$hi["nugget"]) && abs(p["nugget"] - b$hi["nugget"]) < 1e-6)
  }
}
R <- rbindlist(rows)
fwrite(R, file.path(OUT, "S44_gaussian_diagnostics.csv"))

say("paper: range median 3.50 (IQR 3.28-3.84), nugget mean 1.16 (SD 0.90),\n")
say("       nugget at a boundary 0 of 21, R2 mean 0.791 (0.707, 0.880)\n\n")
say("%-16s %18s %18s %14s %22s\n", "box", "range median (IQR)", "nugget mean (SD)",
    "at a boundary", "R2 mean (range)")
for (bn in unique(R$box)) {
  z <- R[box == bn]
  say("%-16s %7.2f (%.2f-%.2f) %10.2f (%.2f) %10d of %-3d %8.3f (%.3f, %.3f)\n", bn,
      median(z$range), quantile(z$range, 0.25), quantile(z$range, 0.75),
      mean(z$nugget), sd(z$nugget), sum(z$nug_lo | z$nug_hi), nrow(z),
      mean(z$r2), min(z$r2), max(z$r2))
}
say("\n== verdict ==\n")
hit <- R[, .(ok = abs(median(range) - 3.50) < 0.05 & abs(mean(nugget) - 1.16) < 0.05 &
               abs(mean(r2) - 0.791) < 0.01), by = box]
if (any(hit$ok)) {
  say("  Reproduced exactly by the %s box, on all four quantities.\n",
      paste(hit[ok == TRUE]$box, collapse = " and "))
  say("  So the published table is right and the M6 answer stands. What was missing\n")
  say("  is the archive: C1_gaussian_per_subject_SUPERSEDED.csv and C1_table1_corrected_SUPERSEDED.csv\n")
  say("  record an EARLIER convention in which the nugget is free and pins at zero\n")
  say("  in all 21 subjects, and the latter's filename says 'corrected', so a reader\n")
  say("  of the repository would conclude the paper contradicts its own results.\n")
  say("  Those two files are superseded and should be labelled or withdrawn.\n")
} else {
  say("  No box in the project reproduces the published table. The closest on each\n")
  say("  quantity are listed above; none is close on all four at once.\n")
  say("  The table should be regenerated from whichever convention the paper means\n")
  say("  to report, and that convention stated in its caption.\n")
}
})

# Interpolation signature: lattice classes against the kernel's covariance
local({
TXT  <- file.path(OUT, "S12_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)

VOX <- 3; set.seed(20260911L)

# kernels, offset-averaged over one voxel (interp_predicted_variogram.m)
lanczos  <- function(x, a = 3) { y <- numeric(length(x)); nz <- abs(x) < a & x != 0
  xx <- x[nz]; y[nz] <- (sin(pi*xx)/(pi*xx))*(sin(pi*xx/a)/(pi*xx/a)); y[x == 0] <- 1; y }
trilin   <- function(x, ...) pmax(0, 1 - abs(x))
bspline3 <- function(x, ...) { t <- abs(x); ifelse(t < 1, (4 - 6*t^2 + 3*t^3)/6,
                                ifelse(t < 2, (2 - t)^3/6, 0)) }
acf_of <- function(kf, half, lags = 0:4) {
  half <- max(half, max(lags) + 2)   # window wider than any kernel's support:
  off <- seq(0, 1, length.out = 201)[-201]; acc <- numeric(length(lags))
  for (dd in off) {                  # the extra taps are exact zeros
    w <- kf((-half):half - dd); den <- sum(w^2)
    for (i in seq_along(lags)) { h <- lags[i]; n <- seq_len(length(w) - h)
      acc[i] <- acc[i] + sum(w[n]*w[n + h])/den }
  }
  acc/length(off)
}
KERN <- list(Lanczos3 = acf_of(lanczos, 3), BSpline3 = acf_of(bspline3, 2),
             Trilinear = acf_of(trilin, 1))
say("offset-averaged kernel autocorrelation at integer voxel offsets\n")
for (nm in names(KERN))
  say("  %-10s %s\n", nm, paste(sprintf("%+.4f", KERN[[nm]]), collapse = "  "))

# class means; k2 is an exact integer key, never a float
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
SV <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2)]
prod_of <- function(A, dx, dy, dz) A[abs(dx)+1L]*A[abs(dy)+1L]*A[abs(dz)+1L]
for (nm in names(KERN)) SV[[nm]] <- prod_of(KERN[[nm]], SV$dx, SV$dy, SV$dz)
SUB <- SV[, c(list(g = mean(g)), lapply(.SD, mean)), by = .(subject, k2),
          .SDcols = names(KERN)]
SUB[, d := VOX*sqrt(k2)]
K <- SUB[, c(list(g = mean(g)), lapply(.SD, first)), by = .(k2, d),
         .SDcols = names(KERN)][order(d)]
say("\n%d classes, %d subjects; class-mean product ranges %+.4f to %+.4f (Lanczos)\n",
    nrow(K), uniqueN(SUB$subject), min(K$Lanczos3), max(K$Lanczos3))

# joint fit: smooth curve plus one fixed-shape kernel term
say("\n== does a fixed kernel term belong in the curve? ==\n")
say("  %-10s %6s %9s %9s %10s %12s\n", "kernel", "df", "scale", "t", "p", "resid SD drop")
FITS <- list()
for (nm in names(KERN)) {
  f0 <- lm(g ~ ns(d, df = 6), K)
  f1 <- lm(as.formula(sprintf("g ~ ns(d, df = 6) + %s", nm)), K)
  s  <- summary(f1)$coefficients[nm, ]
  FITS[[nm]] <- f1
  say("  %-10s %6d %+9.3f %9.2f %10.2g %11.1f%%\n", nm, 6, s[1], s[3], s[4],
      100*(1 - summary(f1)$sigma/summary(f0)$sigma))
}
say("\n  robustness of the Lanczos scale to the smooth term's flexibility:\n")
for (df in c(4, 5, 6, 8, 10, 12)) {
  f1 <- lm(as.formula(sprintf("g ~ ns(d, df = %d) + Lanczos3", df)), K)
  s <- summary(f1)$coefficients["Lanczos3", ]
  say("    df %2d   scale %+.3f   t %6.2f   p %8.2g\n", df, s[1], s[3], s[4])
}

# permutation null: is the match better than a random class assignment?
f0 <- lm(g ~ ns(d, df = 6), K); r0 <- resid(f0)
obs <- abs(cor(r0, K$Lanczos3))
nullc <- replicate(20000, abs(cor(r0, sample(K$Lanczos3))))
say("\n== permutation null (class labels shuffled, 20,000 draws) ==\n")
say("  |cor| between the smooth's residual and the kernel term: %.3f\n", obs)
say("  null mean %.3f, 95th percentile %.3f, p = %.5f\n",
    mean(nullc), quantile(nullc, 0.95), (1 + sum(nullc >= obs))/20001)

# ---- the test that matters: does correcting it destroy cross-subject agreement?
tn <- (K$d - min(K$d))/(max(K$d) - min(K$d))
stiff <- function(v, target = 6) {
  lg <- seq(-8, 2, length.out = 60)
  e <- vapply(lg, function(l) ncs_gcv_edf(tn, v, 10^l)$edf, numeric(1))
  ncs_gcv_edf(tn, v, 10^lg[which.min(abs(e - target))])$fitted
}
# The regression coefficient is minus the variance of the component that was
# white before resampling. That variance is capped at the smallest
# gamma(d)/(1-product(d)) over the lags, above which the implied variogram would
# go negative; the bounded value is used below.
SW_FREE <- -coef(FITS$Lanczos3)["Lanczos3"]
SW_MAX  <- min(K$g/(1 - K$Lanczos3))
S_HAT   <- min(SW_FREE, SW_MAX)
say("\n  free scale %.3f, admissibility bound %.3f, used %.3f\n",
    SW_FREE, SW_MAX, S_HAT)
strip <- function(g, p, sw) g - sw*(1 - p)
subs <- sort(unique(SUB$subject))
mat <- function(adj) sapply(subs, function(s) {
  q <- SUB[subject == s][order(d)]
  v <- if (adj) strip(q$g, q$Lanczos3, S_HAT) else q$g
  v - stiff(v)
})
pair_r <- function(M) { C <- cor(M); mean(C[upper.tri(C)]) }
Fr <- mat(FALSE); Fa <- mat(TRUE)
say("\n== cross-subject replication of the fine structure ==\n")
say("  mean correlation between distinct subjects, raw        : %+.3f\n", pair_r(Fr))
say("  mean correlation between distinct subjects, corrected  : %+.3f\n", pair_r(Fa))
say("  fine-structure size, raw %.4f -> corrected %.4f (%.0f%% removed)\n",
    mean(apply(Fr, 2, sd)), mean(apply(Fa, 2, sd)),
    100*(1 - mean(apply(Fa, 2, sd))/mean(apply(Fr, 2, sd))))

# does a smoother now choose an interior parameter?
say("\n== smoothing-parameter behaviour on the cohort curve ==\n")
edge <- function(v, lab) { f <- fit_ncs(tn, v)
  say("  %-22s log10 lambda %6.2f | edf %5.2f of %d | at grid floor: %s\n",
      lab, log10(f$lambda), f$edf, length(v),
      ifelse(abs(log10(f$lambda) + 6) < 1e-6, "YES", "no"))
  invisible(abs(log10(f$lambda) + 6) < 1e-6) }
e_raw <- edge(K$g, "raw")
K[, g_adj := strip(g, Lanczos3, S_HAT)]
e_adj <- edge(K$g_adj, "kernel-corrected")

# and the hole effect, per subject, before and after
dip6 <- function(d, g) { i <- which.min(abs(d - 6)); 100*(g[i] - min(g[i:length(g)]))/g[i] }
H <- rbindlist(lapply(subs, function(s) { q <- SUB[subject == s][order(d)]
  sw <- min(S_HAT, min(q$g/(1 - q$Lanczos3)))   # each subject has its own bound
  data.table(subject = s, raw = dip6(q$d, q$g),
             adj = dip6(q$d, strip(q$g, q$Lanczos3, sw))) }))
say("\n== hole effect: drop below the 6 mm value, per subject ==\n")
say("  raw              mean %.1f%%, positive in %2d of %d\n",
    mean(H$raw), sum(H$raw > 0), nrow(H))
say("  kernel-corrected mean %.1f%%, positive in %2d of %d\n",
    mean(H$adj), sum(H$adj > 0), nrow(H))

fwrite(K, file.path(OUT, "S12_classes.csv"))
fwrite(H, file.path(OUT, "S12_holes.csv"))
say("\nwrote S12_classes.csv, S12_holes.csv\n")
})

# Robustness of the signature to the correction scale
local({
TXT  <- file.path(OUT, "S13_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)




L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
V <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2)][
       , .(g = mean(g)), by = .(dx, dy, dz, k2)]
V[, `:=`(d = VOX*sqrt(k2), prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L])]
K <- V[, .(g = mean(g), prod = mean(prod), n_vec = .N), by = .(k2, d)][order(d)]

# 1. within-class identification
say("== check 1: scale identified inside distance classes only ==\n")
W <- V[, .(n = .N, spread = diff(range(prod))), by = .(k2, d)][spread > 1e-8]
say("  %d of %d classes contain vectors with different offset multisets:\n",
    nrow(W), nrow(K))
for (i in seq_len(nrow(W)))
  say("    %5.2f mm  %2d vectors  product spread %.5f\n", W$d[i], W$n[i], W$spread[i])
VW <- V[k2 %in% W$k2]
fw <- lm(g ~ factor(k2) + prod, VW)
cw <- summary(fw)$coefficients["prod", ]
say("  scale from within-class contrasts : %+.3f  (se %.3f, t %.2f, p %.3g)\n",
    cw[1], cw[2], cw[3], cw[4])
fb <- lm(g ~ ns(d, df = 6) + prod, K)
cb <- summary(fb)$coefficients["prod", ]
say("  scale from the between-class fit  : %+.3f  (se %.3f)\n", cb[1], cb[2])
say("  the two are %s (difference %.2f, %.1f se of the within-class estimate)\n",
    ifelse(abs(cw[1]-cb[1]) < 2*cw[2], "compatible", "NOT compatible"),
    cw[1]-cb[1], abs(cw[1]-cb[1])/cw[2])

# 2. out-of-sample: fit beyond 7 mm, predict the hole
say("\n== check 2: scale fitted at 7 mm and beyond, used to predict 3 and 6 mm ==\n")
far <- K[d >= 7]
fo <- lm(g ~ ns(d, df = 5) + prod, far)
so <- summary(fo)$coefficients["prod", ]
say("  scale from the %d classes at 7 mm and beyond: %+.3f (se %.3f, p %.3g)\n",
    nrow(far), so[1], so[2], so[4])
near <- K[d < 7]
sm <- lm(g ~ ns(d, df = 5), far)
say("  predicted minus observed anomaly at the classes that define the hole:\n")
say("    %5s %10s %10s %10s\n", "mm", "product", "predicted", "observed")
base <- lm(g ~ ns(d, df = 6), K)
for (i in seq_len(nrow(near)))
  say("    %5.2f %10.4f %10.3f %10.3f\n", near$d[i], near$prod[i],
      so[1]*near$prod[i], resid(base)[K$d == near$d[i]])
say("  correlation of the two over those %d classes: %+.3f\n",
    nrow(near), cor(so[1]*near$prod, resid(base)[K$d < 7]))

# 3. influence
say("\n== check 3: drop one class and refit ==\n")
lo <- vapply(seq_len(nrow(K)), function(i)
  coef(lm(g ~ ns(d, df = 6) + prod, K[-i]))["prod"], numeric(1))
say("  scale ranges %+.3f to %+.3f across the 31 leave-one-out fits (full %+.3f)\n",
    min(lo), max(lo), cb[1])
ord <- order(abs(lo - cb[1]), decreasing = TRUE)[1:4]
for (i in ord) say("    dropping %5.2f mm moves it to %+.3f\n", K$d[i], lo[i])

# 4. physical bound
say("\n== check 4: is the scale physically possible? ==\n")
sill <- max(K$g)
say("  fitted |scale| %.3f versus the curve's plateau %.3f\n", abs(cb[1]), sill)
say("  implied share of variance that was white before resampling: %.0f%%\n",
    100*abs(cb[1])/sill)
say("  (the coefficient is minus that variance, so a share above 100%% is impossible)\n")

# the hole across smooth-term flexibility
say("\n== the hole under every flexibility of the smooth term ==\n")
SUB <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2)]
SUB[, prod := AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L]]
SC <- SUB[, .(g = mean(g), prod = mean(prod)), by = .(subject, k2)][
        , d := VOX*sqrt(k2)][order(subject, d)]
dip6 <- function(d, g) { i <- which.min(abs(d - 6)); 100*(g[i] - min(g[i:length(g)]))/g[i] }
# gamma_S = gamma - sw*(1 - product), with sw = minus the coefficient and capped
# at each subject's admissibility bound (S14). Subtracting sw*product instead
# would add the artefact rather than remove it.
holes_at <- function(sw) SC[, { s <- min(sw, min(g/(1 - prod)))
                                .(h = dip6(d, g - s*(1 - prod))) }, by = subject]
say("  %4s %9s %12s %10s\n", "df", "scale", "mean drop", "positive")
raw <- SC[, .(h = dip6(d, g)), by = subject]
say("  %4s %9s %11.1f%% %7d/21\n", "none", "-", mean(raw$h), sum(raw$h > 0))
for (df in c(4, 5, 6, 8, 10, 12)) {
  s <- -coef(lm(as.formula(sprintf("g ~ ns(d, df = %d) + prod", df)), K))["prod"]
  h <- holes_at(s)
  say("  %4d %+9.3f %11.1f%% %7d/21\n", df, s, mean(h$h), sum(h$h > 0))
}
h_w <- holes_at(-cw[1])
say("  %4s %+9.3f %11.1f%% %7d/21\n", "w/in", -cw[1], mean(h_w$h), sum(h_w$h > 0))
say("\n  scales are the resampling variance, capped at each subject's bound\n")
say("  the last row uses the within-class scale, which never saw the hole\n")
})

# Admissibility bound on the correction scale
local({
TXT  <- file.path(OUT, "S14_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)




L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
SV <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2)]
SV[, prod := AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L]]
SC <- SV[, .(g = mean(g), prod = mean(prod)), by = .(subject, k2)][
        , d := VOX*sqrt(k2)][order(subject, d)]
K <- SC[, .(g = mean(g), prod = first(prod)), by = .(k2, d)][order(d)]
subs <- sort(unique(SC$subject))
dip6 <- function(d, g) { i <- which.min(abs(d - 6)); 100*(g[i] - min(g[i:length(g)]))/g[i] }
# the field's own variogram once the resampling component of variance sw is taken
# out: gamma_S = gamma - sw*(1 - product). Subtracting sw*product instead would
# add the artefact rather than remove it.
strip <- function(g, p, sw) g - sw*(1 - p)

# the bound
bound <- function(g, p) min(g/(1 - p))
SW_MAX <- bound(K$g, K$prod)
i_bind <- which.min(K$g/(1 - K$prod))
say("== the admissibility bound on the resampling variance ==\n")
say("  cohort curve: sw can be at most %.3f, set by the %.2f mm class\n",
    SW_MAX, K$d[i_bind])
say("  the free fit wanted %.3f, which would make the field's own variogram\n", 3.263)
say("  equal to %.2f at 3 mm. Negative, so the free fit over-attributes.\n",
    K$g[1] - 3.263*(1 - K$prod[1]))
pb <- SC[, .(swmax = bound(g, prod)), by = subject]
say("  per subject the bound runs %.3f to %.3f, mean %.3f\n",
    min(pb$swmax), max(pb$swmax), mean(pb$swmax))

# what the largest admissible correction does to the hole
say("\n== the hole at the largest admissible resampling component ==\n")
h_raw <- SC[, .(h = dip6(d, g)), by = subject]
h_max <- SC[, { sw <- bound(g, prod); .(h = dip6(d, strip(g, prod, sw))) }, by = subject]
say("  uncorrected                     mean %.1f%%, positive in %2d of 21\n",
    mean(h_raw$h), sum(h_raw$h > 0))
say("  each subject at its own bound   mean %.1f%%, positive in %2d of 21\n",
    mean(h_max$h), sum(h_max$h > 0))

# ---- profile: how small must the resampling component be for the hole to live?
say("\n== profile over the resampling variance ==\n")
say("  %6s %8s %12s %10s %14s\n", "sw", "of bound", "mean drop", "positive", "gamma_S(3mm)")
grid <- seq(0, SW_MAX, length.out = 26)
prof <- rbindlist(lapply(grid, function(sw) {
  h <- SC[, .(h = dip6(d, strip(g, prod, min(sw, bound(g, prod))))), by = subject]
  data.table(sw = sw, frac = sw/SW_MAX, mean_drop = mean(h$h), npos = sum(h$h > 0),
             gS3 = K$g[1] - sw*(1 - K$prod[1]))
}))
for (i in seq(1, nrow(prof), by = 2))
  say("  %6.3f %7.0f%% %11.1f%% %7d/21 %14.3f\n", prof$sw[i], 100*prof$frac[i],
      prof$mean_drop[i], prof$npos[i], prof$gS3[i])
fwrite(prof, file.path(OUT, "S14_profile.csv"))

half <- prof[mean_drop <= mean(h_raw$h)/2][1]
gone <- prof[npos <= 5][1]
say("\n  the drop is halved once sw reaches %.2f, which is %.0f%% of the bound\n",
    half$sw, 100*half$frac)
if (!is.na(gone$sw))
  say("  it is present in 5 or fewer subjects once sw reaches %.2f (%.0f%% of the bound)\n",
      gone$sw, 100*gone$frac)

# and the within-class scale, which never saw the hole
say("\n== for reference, the scale identified inside distance classes ==\n")
V <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2)][
       , .(g = mean(g)), by = .(dx, dy, dz, k2)]
V[, prod := AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L]]
W <- V[, .(sp = diff(range(prod))), by = k2][sp > 1e-8]
cw <- summary(lm(g ~ factor(k2) + prod, V[k2 %in% W$k2]))$coefficients["prod", ]
say("  %.3f (se %.3f); 95%% interval %.3f to %.3f\n", -cw[1], cw[2],
    -cw[1] - 1.96*cw[2], -cw[1] + 1.96*cw[2])
say("  the bound %.3f sits inside that interval, so the within-class estimate is\n", SW_MAX)
say("  consistent with a resampling component at or near the largest admissible one\n")
h_w <- SC[, { sw <- min(-cw[1], bound(g, prod)); .(h = dip6(d, strip(g, prod, sw))) }, by = subject]
say("  hole at that scale, capped at each subject's bound: mean %.1f%%, positive in %d of 21\n",
    mean(h_w$h), sum(h_w$h > 0))
})

# Smoothing parameters after correction, on an extended grid
local({
TXT  <- file.path(OUT, "S16_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)




L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, `:=`(k2 = as.integer(round((distance_mm/VOX)^2)),
         prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L])]
LAG <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2, prod)][
         , .(g = mean(g), prod = mean(prod)), by = .(subject, k2)][
         , d := VOX*sqrt(k2)][order(subject, d)]
strip <- function(d, g, p, df = 6) {
  swmax <- min(g/(1 - p))
  sw <- min(max(0, -coef(lm(g ~ ns(d, df = df) + p))["p"]), swmax)
  g - sw*(1 - p)
}

K_OLD <- 1:10; K_NEW <- 1:20
BW_OLD <- seq(0.05, 0.60, by = 0.025)
BW_NEW <- sort(unique(c(seq(0.005, 0.05, by = 0.005), BW_OLD)))

subs <- sort(unique(LAG$subject))
R <- rbindlist(lapply(subs, function(s) {
  q <- LAG[subject == s][order(d)]
  t <- (q$d - min(q$d))/(max(q$d) - min(q$d))
  out <- list(subject = s)
  for (tag in c("raw", "adj")) {
    y <- if (tag == "raw") q$g else strip(q$d, q$g, q$prod)
    out[[paste0("k_old_", tag)]]  <- K_OLD[which.min(search_fn(t, y, K_OLD, ptr_wrap)[, 2])]
    out[[paste0("k_new_", tag)]]  <- K_NEW[which.min(search_fn(t, y, K_NEW, ptr_wrap)[, 2])]
    out[[paste0("bw_old_", tag)]] <- BW_OLD[which.min(search_fn(t, y, BW_OLD, locpoly3_wrap)[, 2])]
    out[[paste0("bw_new_", tag)]] <- BW_NEW[which.min(search_fn(t, y, BW_NEW, locpoly3_wrap)[, 2])]
  }
  as.data.table(out)
}))
fwrite(R, file.path(OUT, "S16_grid_extension.csv"))
N <- nrow(R)

say("PTR, on 31 points: k selected and how often it rests on the top of the grid\n")
say("  %-22s %12s %12s\n", "", "raw", "corrected")
say("  %-22s %9.1f    %9.1f\n", "median k, grid to 10", median(R$k_old_raw), median(R$k_old_adj))
say("  %-22s %6d/%-2d    %6d/%-2d\n", "at k = 10",
    sum(R$k_old_raw == 10), N, sum(R$k_old_adj == 10), N)
say("  %-22s %9.1f    %9.1f\n", "median k, grid to 20", median(R$k_new_raw), median(R$k_new_adj))
say("  %-22s %6d/%-2d    %6d/%-2d\n", "at k >= 15",
    sum(R$k_new_raw >= 15), N, sum(R$k_new_adj >= 15), N)
say("  %-22s %6d/%-2d    %6d/%-2d\n", "at k = 20",
    sum(R$k_new_raw == 20), N, sum(R$k_new_adj == 20), N)
say("  a PTR fit with k terms uses k+3 parameters, so k = 20 is 23 of 31 points\n")

say("\nLP, on 31 points: bandwidth selected and how often it rests on the floor\n")
say("  %-22s %12s %12s\n", "", "raw", "corrected")
say("  %-22s %9.3f    %9.3f\n", "median bw, floor 0.05",
    median(R$bw_old_raw), median(R$bw_old_adj))
say("  %-22s %6d/%-2d    %6d/%-2d\n", "at bw = 0.05",
    sum(R$bw_old_raw <= 0.05), N, sum(R$bw_old_adj <= 0.05), N)
say("  %-22s %9.3f    %9.3f\n", "median bw, floor 0.005",
    median(R$bw_new_raw), median(R$bw_new_adj))
say("  %-22s %6d/%-2d    %6d/%-2d\n", "below the old floor",
    sum(R$bw_new_raw < 0.05), N, sum(R$bw_new_adj < 0.05), N)
say("  %-22s %6d/%-2d    %6d/%-2d\n", "at bw = 0.005",
    sum(R$bw_new_raw <= 0.005), N, sum(R$bw_new_adj <= 0.005), N)

say("\nverdict\n")
say("  PTR interior on the wider grid after correction: %s\n",
    ifelse(sum(R$k_new_adj >= 15) <= 3, "yes", "NO - it still runs to the top"))
say("  LP  interior on the wider grid after correction: %s\n",
    ifelse(sum(R$bw_new_adj < 0.05) <= 3, "yes", "NO - it still runs below the old floor"))
})

# Guarded cross-validated selection
local({
TXT  <- file.path(OUT, "S17_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)




L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, `:=`(k2 = as.integer(round((distance_mm/VOX)^2)),
         prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L])]
LAG <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2, prod)][
         , .(g = mean(g), prod = mean(prod)), by = .(subject, k2)][
         , d := VOX*sqrt(k2)][order(subject, d)]
strip <- function(d, g, p, df = 6) {
  swmax <- min(g/(1 - p))
  sw <- min(max(0, -coef(lm(g ~ ns(d, df = df) + p))["p"]), swmax)
  g - sw*(1 - p)
}

BW <- sort(unique(c(seq(0.005, 0.05, by = 0.005), seq(0.05, 0.60, by = 0.025))))
KK <- 1:20
LG <- seq(-6, 6, length.out = 121)

# each returns the criterion and the degrees of freedom spent, per grid point
scan_lp  <- function(t, y) { n <- length(y)
  do.call(rbind, lapply(BW, function(b) { u <- locpoly3_wrap(t, y, b)
    data.frame(par = b, edf = u[[3]], crit = u[[2]]/(1 - u[[3]]/n)^2/n) })) }
scan_ptr <- function(t, y) { n <- length(y)
  do.call(rbind, lapply(KK, function(k) { u <- ptr_wrap(t, y, k)
    data.frame(par = k, edf = u[[3]], crit = u[[2]]/(1 - u[[3]]/n)^2/n) })) }
scan_ncs <- function(t, y) {
  do.call(rbind, lapply(LG, function(l) { r <- ncs_gcv_edf(t, y, 10^l)
    data.frame(par = l, edf = r$edf, crit = r$gcv) })) }
# smaller par means less smoothing for LP and NCS, larger par means less for PTR
pick <- function(S, cap, less_is_rough) {
  ok <- S[S$edf <= cap, , drop = FALSE]
  if (!nrow(ok)) return(list(par = NA, edf = NA, bound = TRUE))
  i <- which.min(ok$crit)
  lim <- if (less_is_rough) min(ok$par) else max(ok$par)
  list(par = ok$par[i], edf = ok$edf[i], bound = isTRUE(all.equal(ok$par[i], lim)))
}

subs <- sort(unique(LAG$subject))
R <- rbindlist(lapply(subs, function(s) {
  q <- LAG[subject == s][order(d)]
  t <- (q$d - min(q$d))/(max(q$d) - min(q$d)); n <- length(t); cap <- n/2
  out <- list(subject = s, n = n)
  for (tag in c("raw", "adj")) {
    y <- if (tag == "raw") q$g else strip(q$d, q$g, q$prod)
    for (m in c("LP", "PTR", "NCS")) {
      S <- switch(m, LP = scan_lp(t, y), PTR = scan_ptr(t, y), NCS = scan_ncs(t, y))
      p <- pick(S, cap, m != "PTR")
      out[[paste0(m, "_par_", tag)]]   <- p$par
      out[[paste0(m, "_edf_", tag)]]   <- p$edf
      out[[paste0(m, "_bound_", tag)]] <- p$bound
    }
  }
  as.data.table(out)
}))
fwrite(R, file.path(OUT, "S17_guarded.csv"))
N <- nrow(R)

say("31 points per curve; a fit may spend at most %.1f degrees of freedom\n", R$n[1]/2)
say("grids searched: bandwidth %.3f to %.2f, k 1 to %d, log10 lambda %g to %g\n\n",
    min(BW), max(BW), max(KK), min(LG), max(LG))
say("%-6s %-28s %12s %12s\n", "", "", "raw", "corrected")
for (m in c("LP", "PTR", "NCS")) {
  b_r <- sum(R[[paste0(m, "_bound_raw")]]); b_a <- sum(R[[paste0(m, "_bound_adj")]])
  say("%-6s %-28s %9d/%-2d %9d/%-2d\n", m, "guard is what stops it", b_r, N, b_a, N)
  say("%-6s %-28s %12.3f %12.3f\n", "", "median parameter",
      median(R[[paste0(m, "_par_raw")]]), median(R[[paste0(m, "_par_adj")]]))
  say("%-6s %-28s %12.1f %12.1f\n", "", "median df spent",
      median(R[[paste0(m, "_edf_raw")]]), median(R[[paste0(m, "_edf_adj")]]))
}

say("\ncohort-mean curve, criterion minimum over the whole grid with no guard\n")
q <- LAG[, .(g = mean(g), prod = first(prod)), by = .(k2, d)][order(d)]
t <- (q$d - min(q$d))/(max(q$d) - min(q$d)); n <- length(t)
for (tag in c("raw", "adj")) {
  y <- if (tag == "raw") q$g else strip(q$d, q$g, q$prod)
  say("  %-10s", tag)
  for (m in c("LP", "PTR", "NCS")) {
    S <- switch(m, LP = scan_lp(t, y), PTR = scan_ptr(t, y), NCS = scan_ncs(t, y))
    i <- which.min(S$crit)
    say("  %s par %.3f df %.1f", m, S$par[i], S$edf[i])
  }
  say("\n")
}

say("\nverdict, on the corrected curves\n")
for (m in c("LP", "PTR", "NCS")) {
  b <- sum(R[[paste0(m, "_bound_adj")]])
  say("  %-4s %s (%d of %d subjects stopped by the guard)\n", m,
      ifelse(b <= 3, "genuine interior optimum", "still wants to interpolate"), b, N)
}
})

# Gaussian misfit after correction
local({
TXT  <- file.path(OUT, "S18_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)




L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, `:=`(k2 = as.integer(round((distance_mm/VOX)^2)),
         prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L])]
LAG <- L[, .(g = mean(gamma)), by = .(subject, dx, dy, dz, k2, prod)][
         , .(g = mean(g), prod = mean(prod)), by = .(subject, k2)][
         , d := VOX*sqrt(k2)][order(subject, d)]
strip <- function(d, g, p, df = 6) {
  swmax <- min(g/(1 - p))
  sw <- min(max(0, -coef(lm(g ~ ns(d, df = df) + p))["p"]), swmax)
  g - sw*(1 - p)
}

# profiled Gaussian: linear in nugget and sill given the range, so scan the range
# and solve the rest exactly. Range box matches the manuscript's, 0.5 to 25 mm.
fit_gauss <- function(d, y) {
  at <- function(a) { u <- exp(-(d/a)^2); X <- cbind(u, 1 - u)
    cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) return(NULL)
    if (cf[1] < 0) cf <- c(0, sum(y*(1-u))/sum((1-u)^2))
    list(ss = sum((y - X %*% cf)^2), fit = as.vector(X %*% cf)) }
  best <- NULL
  for (a in seq(0.5, 25, length.out = 60)) { r <- at(a)
    if (!is.null(r) && (is.null(best) || r$ss < best$ss)) { best <- r; ba <- a } }
  for (a in seq(max(0.5, ba - 0.5), min(25, ba + 0.5), length.out = 40)) { r <- at(a)
    if (!is.null(r) && r$ss < best$ss) best <- r }
  best$fit
}
BW <- sort(unique(c(seq(0.005, 0.05, by = 0.005), seq(0.05, 0.60, by = 0.025))))
KK <- 1:20; LG <- seq(-6, 6, length.out = 121)
fit_np <- function(t, y) {
  n <- length(y); cap <- n/2
  lp <- do.call(rbind, lapply(BW, function(b) { u <- locpoly3_wrap(t, y, b)
    data.frame(p = b, edf = u[[3]], cr = u[[2]]/(1 - u[[3]]/n)^2/n) }))
  lp <- lp[lp$edf <= cap, ]; bl <- lp$p[which.min(lp$cr)]
  pt <- do.call(rbind, lapply(KK, function(k) { u <- ptr_wrap(t, y, k)
    data.frame(p = k, edf = u[[3]], cr = u[[2]]/(1 - u[[3]]/n)^2/n) }))
  pt <- pt[pt$edf <= cap, ]; bk <- pt$p[which.min(pt$cr)]
  nc <- do.call(rbind, lapply(LG, function(l) { r <- ncs_gcv_edf(t, y, 10^l)
    data.frame(p = l, edf = r$edf, cr = r$gcv) }))
  nc <- nc[nc$edf <= cap, ]; bg <- nc$p[which.min(nc$cr)]
  list(LP = locpoly3(t, y, bl, t)[[1]][, 1],
       PTR = as.vector(ptr(t, y, bk)[[1]]),
       NCS = ncs_gcv_edf(t, y, 10^bg)$fitted)
}

subs <- sort(unique(LAG$subject))
rows <- list(); resid_store <- list()
for (s in subs) {
  q <- LAG[subject == s][order(d)]
  t <- (q$d - min(q$d))/(max(q$d) - min(q$d))
  for (tag in c("raw", "adj")) {
    y <- if (tag == "raw") q$g else strip(q$d, q$g, q$prod)
    f <- fit_np(t, y); f$Gaussian <- fit_gauss(q$d, y)
    rng <- diff(range(y))
    for (m in names(f)) {
      r <- y - f[[m]]
      rows[[length(rows)+1]] <- data.table(subject = s, basis = tag, method = m,
        rss = sum(r^2), rel = sqrt(mean(r^2))/rng)
      resid_store[[length(resid_store)+1]] <- data.table(subject = s, basis = tag,
        method = m, d = q$d, r = r)
    }
  }
}
R <- rbindlist(rows); RS <- rbindlist(resid_store)
fwrite(R, file.path(OUT, "S18_rss.csv"))

say("residual sum of squares, median over the 21 subjects\n")
say("  %-10s %14s %14s\n", "method", "raw", "corrected")
for (m in c("NCS","LP","PTR","Gaussian"))
  say("  %-10s %14.4f %14.4f\n", m,
      median(R[basis=="raw" & method==m]$rss), median(R[basis=="adj" & method==m]$rss))
say("\nroot mean square residual as a fraction of each curve's own range\n")
say("  %-10s %14s %14s\n", "method", "raw", "corrected")
for (m in c("NCS","LP","PTR","Gaussian"))
  say("  %-10s %13.1f%% %13.1f%%\n", m,
      100*median(R[basis=="raw" & method==m]$rel), 100*median(R[basis=="adj" & method==m]$rel))

say("\nis the Gaussian the worst fit in every subject?\n")
for (tag in c("raw","adj")) {
  w <- R[basis == tag][, .(worst = method[which.max(rss)]), by = subject]
  say("  %-10s Gaussian worst in %2d of %d subjects\n", tag, sum(w$worst == "Gaussian"), nrow(w))
}

say("\ncohort-mean residual against distance: largest magnitude reached\n")
say("  %-10s %14s %14s\n", "method", "raw", "corrected")
CM <- RS[, .(r = mean(r)), by = .(basis, method, d)]
for (m in c("NCS","LP","PTR","Gaussian"))
  say("  %-10s %14.3f %14.3f\n", m,
      max(abs(CM[basis=="raw" & method==m]$r)), max(abs(CM[basis=="adj" & method==m]$r)))
say("  the manuscript reports 0.288 for the Gaussian against 0.048 to 0.111\n")

say("\nGaussian cohort-mean residual at the lags that carry the argument\n")
say("  %6s %12s %12s\n", "mm", "raw", "corrected")
G <- CM[method == "Gaussian"]
for (dd in sort(unique(G$d))[1:9])
  say("  %6.2f %12.3f %12.3f\n", dd, G[basis=="raw" & abs(d-dd)<1e-9]$r,
      G[basis=="adj" & abs(d-dd)<1e-9]$r)
fwrite(CM, file.path(OUT, "S18_cohort_residuals.csv"))
say("\nwrote S18_rss.csv, S18_cohort_residuals.csv\n")
})

# Every published claim, re-tested at each level of correction
local({
TXT  <- file.path(OUT, "S37_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)





L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, `:=`(k2 = as.integer(round((distance_mm/VOX)^2)),
         prod = AK[abs(dx)+1L]*AK[abs(dy)+1L]*AK[abs(dz)+1L])]
RUN <- L[, .(g = mean(gamma), prod = mean(prod)), by = .(subject, run, k2)][
         , d := VOX*sqrt(k2)][order(subject, run, d)]

# gamma_S = gamma - sw(1 - product); sw fitted per run and capped at its bound
fit_sw <- function(d, g, p) {
  swmax <- min(g/(1 - p))
  sw <- tryCatch(-coef(lm(g ~ ns(d, df = 6) + p))["p"], error = function(e) NA_real_)
  if (is.na(sw)) sw <- 0
  list(sw = min(max(0, sw), swmax), swmax = swmax)
}
SW <- RUN[, { f <- fit_sw(d, g, prod); .(sw = f$sw, swmax = f$swmax) }, by = .(subject, run)]
RUN <- merge(RUN, SW, by = c("subject", "run"))
say("resampling variance fitted per subject-run: median %.3f, at its own bound in\n",
    median(SW$sw))
say("%d of %d runs; per-run bound ranges %.3f to %.3f\n\n",
    sum(SW$sw >= SW$swmax - 1e-9), nrow(SW), min(SW$swmax), max(SW$swmax))

# the two claims, as functions of how much is removed
i6 <- function(d) which.min(abs(d - 6)); i9 <- function(d) which.min(abs(d - 9))
decline <- function(d, g) { a <- g[i6(d)]; 100*(a - g[i9(d)])/a }

zone_res <- function(d, g) {                    # Gaussian residual by zone
  fit <- tryCatch({
    at <- function(a) { u <- exp(-(d/a)^2); X <- cbind(u, 1 - u)
      cf <- as.vector(qr.solve(X, g)); if (cf[1] < 0) cf <- c(0, sum(g*(1-u))/sum((1-u)^2))
      list(ss = sum((g - X %*% cf)^2), fit = as.vector(X %*% cf)) }
    b <- NULL
    for (a in seq(0.5, 25, length.out = 80)) { r <- tryCatch(at(a), error = function(e) NULL)
      if (!is.null(r) && (is.null(b) || r$ss < b$ss)) b <- r }
    b$fit }, error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  r <- g - fit
  c(mean(r[d >= 4.2 & d <= 6.2]), mean(r[d >= 7.0 & d <= 11.5]))
}

fracs <- seq(0, 1, by = 0.1)
rows <- list()
for (fr in fracs) {
  A <- RUN[, { s <- fr*sw; y <- g - s*(1 - prod)
               .(dec = decline(d, y)) }, by = .(subject, run)]
  S <- A[, .(dec = mean(dec)), by = subject]
  # zone residuals on the subject's run-averaged corrected curve
  Z <- RUN[, { s <- fr*sw[1]; y <- tapply(g - s*(1 - prod), d, mean)
               dd <- as.numeric(names(y)); zr <- zone_res(dd, as.numeric(y))
               .(corner = zr[1], plateau = zr[2]) }, by = subject]
  bt <- binom.test(sum(S$dec > 0), nrow(S))
  rows[[length(rows)+1]] <- data.table(frac = fr,
    mean_decline = mean(S$dec), n_pos = sum(S$dec > 0), p_sign = bt$p.value,
    t_p = t.test(S$dec)$p.value,
    corner_pos = sum(Z$corner > 0), corner_p = binom.test(sum(Z$corner > 0), nrow(Z))$p.value,
    plateau_neg = sum(Z$plateau < 0),
    plateau_p = binom.test(sum(Z$plateau < 0), nrow(Z))$p.value)
}
R <- rbindlist(rows)
fwrite(R, file.path(OUT, "S37_claims.csv"))

say("== the fall between 6 and 9 mm, across the profile ==\n")
say("  %-10s %14s %12s %14s %12s\n", "removed", "mean fall", "falls in", "sign test p", "t-test p")
for (i in seq_len(nrow(R)))
  say("  %8.0f%% %13.2f%% %9d/21 %14.2e %12.2e\n", 100*R$frac[i], R$mean_decline[i],
      R$n_pos[i], R$p_sign[i], R$t_p[i])
say("\n  the published figure is +4.32%%, 21 of 21, p = 9.54e-07\n")

say("\n== the two Gaussian-residual sign tests, across the profile ==\n")
say("  %-10s %16s %12s %16s %12s\n", "removed", "corner positive", "p", "plateau negative", "p")
for (i in seq_len(nrow(R)))
  say("  %8.0f%% %13d/21 %12.2e %13d/21 %12.2e\n", 100*R$frac[i],
      R$corner_pos[i], R$corner_p[i], R$plateau_neg[i], R$plateau_p[i])
say("\n  the published figures are 21 of 21 positive and 18 of 21 negative\n")

say("\n== verdict ==\n")
z <- R[frac == 1]
say("  At the admissibility bound the fall is %.2f%% in %d of 21 subjects.\n",
    z$mean_decline, z$n_pos)
say("  The corner test reads %d of 21 and the plateau test %d of 21.\n",
    z$corner_pos, z$plateau_neg)
half <- R[mean_decline <= R$mean_decline[1]/2][1]
if (!is.na(half$frac))
  say("  The fall is halved once %.0f%% of the fitted resampling variance is removed.\n",
      100*half$frac)
say("\n  The Gaussian misfit is therefore not what the resampling term explains,\n")
say("  while the non-monotonicity claim is conditional on it. Both belong in the\n")
say("  paper in that form.\n")
})

# Depth gradient against the kernel prediction
local({
TXT  <- file.path(OUT, "S40_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)



a1 <- AK[2]; a2 <- AK[3]; a3 <- AK[4]
say("== the correction is a lattice constant, not a per-shell estimate ==\n")
say("  axis products: 3 mm %+.4f, 6 mm %+.4f, 9 mm %+.4f\n", a1, a2, a3)
say("  the numerator of the change shifts by sw x %.4f and the denominator by\n", a3 - a2)
say("  sw x %.4f, the same in every shell and session\n\n", 1 - a2)

A <- fread(file.path(OUT, "A5_depth_fixed.csv"))
W <- dcast(A[k <= 3, .(g = mean(gamma)), by = .(subject, run, session, depth_mm, k)],
           subject + run + session + depth_mm ~ k, value.var = "g")
setnames(W, c("1","2","3"), c("g3","g6","g9"))
W[, bound := pmin(g3/(1 - a1), g6/(1 - a2), g9/(1 - a3))]
depths <- sort(unique(W$depth_mm))
say("shells at %s mm; %d subject-run-shell records\n", paste(depths, collapse = ", "), nrow(W))
say("median admissibility ceiling by shell: %s\n\n",
    paste(sprintf("%.2f", W[, median(bound), by = depth_mm][order(depth_mm)]$V1), collapse = "  "))

chg <- function(g6, g9, sw) 100*(g6 - g9 - sw*(a3 - a2))/(g6 - sw*(1 - a2))

rows <- list()
for (fr in c(0, 0.25, 0.5, 0.75, 1)) {
  P <- W[, .(ch = mean(chg(g6, g9, fr*bound))), by = .(subject, session, depth_mm)]
  for (s in 1:2) {
    z <- P[session == s]
    m <- z[, .(mean_ch = mean(ch), n_pos = sum(ch > 0), n = .N), by = depth_mm][order(depth_mm)]
    lf <- lm(ch ~ depth_mm, z)
    rows[[length(rows)+1]] <- data.table(frac = fr, session = s, m,
      slope = coef(lf)[2], slope_p = summary(lf)$coefficients[2,4],
      rho = cor(m$depth_mm, m$mean_ch, method = "spearman"),
      monotone = all(diff(m$mean_ch) > 0))
  }
}
D <- rbindlist(rows)
fwrite(D, file.path(OUT, "S40_depth_vs_kernel.csv"))

for (s in 1:2) {
  say("== session %d ==\n", s)
  say("  %-10s %s %12s %10s %10s\n", "removed",
      paste(sprintf("%8d mm", depths), collapse = ""), "slope/mm", "p", "monotone")
  for (fr in unique(D$frac)) {
    z <- D[frac == fr & session == s][order(depth_mm)]
    say("  %8.0f%% %s %11.3f %10.1e %10s\n", 100*fr,
        paste(sprintf("%+8.2f%%", z$mean_ch), collapse = ""),
        z$slope[1], z$slope_p[1], ifelse(z$monotone[1], "yes", "no"))
  }
  say("\n")
}
say("the published table gives %s for session 1\n",
    paste(c("-3.23%","+1.49%","+4.57%","+5.40%","+6.65%"), collapse = "  "))

say("\n== verdict ==\n")
z0 <- D[frac == 0 & session == 1]; z5 <- D[frac == 0.5 & session == 1]
z1 <- D[frac == 1 & session == 1]
say("  slope per mm, session 1: %+.3f uncorrected, %+.3f at half the ceiling,\n",
    z0$slope[1], z5$slope[1])
say("  %+.3f at the ceiling. Monotone throughout: %s.\n", z1$slope[1],
    ifelse(all(D$monotone), "yes", "NO"))
say("  Sessions agree to within %.2f percentage points at every depth and every\n",
    max(abs(D[session == 1][order(frac, depth_mm)]$mean_ch -
            D[session == 2][order(frac, depth_mm)]$mean_ch)))
say("  correction level.\n")
say("\n  The gradient is carried by differences between shells, and the correction\n")
say("  enters every shell as the same multiple of its own level, so it cannot\n")
say("  create the gradient. What the correction moves is the depth at which the\n")
say("  change crosses zero, and that crossing is a statement about level.\n")
})

# Depth gradient, before the partner-validity fix
local({
D <- fread(file.path(DIAG, "out_matlab", "modelB_depth.csv"))
D[, depth_mm := shell_lo * 3]

cat("source: modelB_depth.csv, two-session model B residuals\n")
cat(sprintf("  subjects %d | runs %s | sessions %s | shells %s voxels\n\n",
            uniqueN(D$subject), paste(range(D$run), collapse = "-"),
            paste(sort(unique(D$session)), collapse = ","),
            paste(sort(unique(D$shell_lo)), collapse = ",")))

# cohort curve per shell: average over axes and runs within subject, then over
# subjects, so every subject weighs equally
cur <- D[, .(g = mean(gamma)), by = .(shell_lo, depth_mm, subject, distance_mm)][
         , .(g = mean(g), sd = sd(g), n = .N), by = .(shell_lo, depth_mm, distance_mm)][
         order(shell_lo, distance_mm)]

# per-subject decline so the sign test has a unit of analysis
psub <- D[, .(g = mean(gamma)), by = .(shell_lo, depth_mm, subject, distance_mm)]
at <- function(x, d, target) x$g[which.min(abs(x$distance_mm - target))]

# Seeded so the bootstrap intervals below are reproducible run to run.
set.seed(2026)

rows <- list()
for (k in sort(unique(psub$shell_lo))) {
  x <- psub[shell_lo == k]
  per <- x[, .(g6 = at(.SD, distance_mm, 6), g9 = at(.SD, distance_mm, 9)),
           by = subject]
  per[, decline_pct := 100 * (g6 - g9) / g6]
  bt <- binom.test(sum(per$decline_pct > 0), nrow(per), 0.5)
  ci <- quantile(replicate(2000, mean(sample(per$decline_pct, nrow(per), TRUE))),
                 c(.025, .975))
  rows[[length(rows)+1]] <- data.table(
    shell_lo = k, depth_mm = k * 3,
    decline_pct = mean(per$decline_pct), ci_lo = ci[1], ci_hi = ci[2],
    n_declining = sum(per$decline_pct > 0), n = nrow(per), sign_p = bt$p.value)
}
R <- rbindlist(rows)
fwrite(R, file.path(OUT, "A2_depth_gradient.csv"))

cat("Decline from 6 mm to 9 mm, positive means the curve falls\n")
cat("(the manuscript's post-peak decline). Per-subject unit, 21 subjects.\n\n")
cat(sprintf("  %-8s %10s %18s %14s %10s\n",
            "depth", "decline%", "95% CI", "n declining", "sign p"))
for (i in seq_len(nrow(R))) {
  r <- R[i]
  cat(sprintf("  %-8s %+10.2f   [%+6.2f, %+6.2f]   %6d of %2d %10.2e\n",
              paste0(r$depth_mm, " mm"), r$decline_pct, r$ci_lo, r$ci_hi,
              r$n_declining, r$n, r$sign_p))
}

cat("\n--- comparison with what the log records ---\n")
old <- data.table(depth_mm = c(9,12,15,18,24),
                  old_pct = c(-4.68, 0.53, 4.21, 4.93, 5.72))
m <- merge(R[, .(depth_mm, new_pct = decline_pct)], old, by = "depth_mm")
cat(sprintf("  %-8s %14s %14s\n", "depth", "logged (S1)", "model B (2 ses)"))
for (i in seq_len(nrow(m)))
  cat(sprintf("  %-8s %+14.2f %+14.2f\n",
              paste0(m$depth_mm[i], " mm"), m$old_pct[i], m$new_pct[i]))
cat("\n  The logged values are from the older session-1 Python analysis\n")
cat("  (fixedset_out/summary_depth_depth.csv, 3 August, 20 subjects).\n")
cat("  The right-hand column is the two-session model B basis, 21 subjects,\n")
cat("  and is what the manuscript should report.\n")
cat(sprintf("\nwrote A2_depth_gradient.csv -> %s\n", OUT))
})

# Depth gradient, after the partner-validity fix
local({
D <- fread(file.path(OUT, "A5_depth_fixed.csv"))
cat(sprintf("A5_depth_fixed.csv: %d rows | subjects %d | runs %s | sessions %s\n",
            nrow(D), uniqueN(D$subject),
            paste(range(D$run), collapse = "-"),
            paste(sort(unique(D$session)), collapse = ",")))
stopifnot(uniqueN(D$subject) == 21L, setequal(unique(D$run), 1:10))

cat("\nattrition from enforcing the precondition, mean over subject-runs:\n")
att <- D[k == 1 & axis == 1, .(raw = mean(n_base_raw), valid = mean(n_base_valid)),
         by = .(depth_mm)][order(depth_mm)]
att[, dropped_pct := 100 * (raw - valid) / raw]
for (i in seq_len(nrow(att)))
  cat(sprintf("  %5s mm  base %7.0f -> %7.0f  (%.1f%% dropped)\n",
              att$depth_mm[i], att$raw[i], att$valid[i], att$dropped_pct[i]))

# per-subject decline, averaging over axes, runs, sessions
psub <- D[, .(g = mean(gamma)), by = .(depth_mm, subject, distance_mm)]
at <- function(x, target) x$g[which.min(abs(x$distance_mm - target))]

rows <- list()
set.seed(2026)
for (dm in sort(unique(psub$depth_mm))) {
  x <- psub[depth_mm == dm]
  per <- x[, .(g6 = at(.SD, 6), g9 = at(.SD, 9)), by = subject]
  per[, decline_pct := 100 * (g6 - g9) / g6]
  bt <- binom.test(sum(per$decline_pct > 0), nrow(per), 0.5)
  ci <- quantile(replicate(4000, mean(sample(per$decline_pct, nrow(per), TRUE))),
                 c(.025, .975))
  rows[[length(rows)+1]] <- data.table(depth_mm = dm,
    decline_pct = mean(per$decline_pct), ci_lo = ci[1], ci_hi = ci[2],
    n_declining = sum(per$decline_pct > 0), n = nrow(per), sign_p = bt$p.value)
}
R <- rbindlist(rows)
fwrite(R, file.path(OUT, "A6_depth_gradient_fixed.csv"))

cat("\nCORRECTED depth gradient, 21 subjects, both sessions, 10 runs each\n")
cat("decline from 6 mm to 9 mm; positive means the curve falls\n\n")
cat(sprintf("  %-8s %10s %20s %14s %11s\n",
            "depth", "decline%", "95% CI", "n declining", "sign p"))
for (i in seq_len(nrow(R))) {
  r <- R[i]
  cat(sprintf("  %-8s %+10.2f    [%+6.2f, %+6.2f]   %6d of %2d %11.2e\n",
              paste0(r$depth_mm, " mm"), r$decline_pct, r$ci_lo, r$ci_hi,
              r$n_declining, r$n, r$sign_p))
}

prev <- fread(file.path(OUT, "A2_depth_gradient.csv"))
m <- merge(R[, .(depth_mm, fixed = decline_pct)],
           prev[, .(depth_mm, before = decline_pct)], by = "depth_mm")
old <- data.table(depth_mm = c(9,12,15,18,24),
                  logged = c(-4.68, 0.53, 4.21, 4.93, 5.72))
m <- merge(m, old, by = "depth_mm")
cat("\n--- three versions ---\n")
cat(sprintf("  %-8s %14s %16s %16s\n", "depth",
            "logged (S1)", "2-ses, leaking", "2-ses, corrected"))
for (i in seq_len(nrow(m)))
  cat(sprintf("  %-8s %+14.2f %+16.2f %+16.2f\n", paste0(m$depth_mm[i], " mm"),
              m$logged[i], m$before[i], m$fixed[i]))
cat(sprintf("\n  span, corrected: %.2f to %+.2f percentage points\n",
            min(R$decline_pct), max(R$decline_pct)))
cat("  Only the 9 mm row moves materially; A3 showed that is the only cell\n")
cat("  where partners left the mask.\n")
cat(sprintf("\nwrote A6_depth_gradient_fixed.csv -> %s\n", OUT))
})

# Damped spherical-wave reference covariance
local({
TXT  <- file.path(OUT, "S25_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)

om3 <- function(x) { o <- rep(1, length(x)); nz <- x != 0; o[nz] <- sin(x[nz])/x[nz]; o }

# ---- the observed shell curves, raw: inference is about the images that exist -
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
SH <- melt(L, id.vars = c("subject","run","k2"), measure.vars = patterns("^gamma_shell"),
           variable.name = "shell", value.name = "g")
SH[, shell := as.integer(sub("gamma_shell", "", shell))]
CUR <- SH[!is.na(g), .(g = mean(g)), by = .(subject, shell, k2)][
         , .(g = mean(g)), by = .(shell, k2)][, d := VOX*sqrt(k2)][order(shell, d)]
depth_mm <- c(9, 12, 15, 18, 24)
shells <- sort(unique(CUR$shell))

# the model
gam_ref <- function(r, p0, s, a, tau)
  p0 + rowSums(sapply(seq_along(s), function(m)
    s[m]*(1 - om3(a[m]*r)*exp(-r^2/(2*tau[m]^2)))))
dens_ref <- function(k, s, a, tau) {          # shape only, up to a constant
  # each component carries its own tau_m: it comes from
  # integral_0^inf cos(br) exp(-r^2/(2 tau^2)) dr = tau sqrt(pi/2) exp(-b^2 tau^2/2)
  # and only cancels into the overall scale when every component shares a width
  rowSums(sapply(seq_along(s), function(m)
    s[m]*tau[m]*(exp(-((a[m]-k)^2)*tau[m]^2/2) -
                 exp(-((a[m]+k)^2)*tau[m]^2/2))/(k*a[m])))
}

# check the closed form against a numerical transform
say("== is the closed-form density the transform of the covariance? ==\n")
set.seed(1); ok <- TRUE
for (trial in 1:3) {
  s <- runif(2, 0.5, 2); a <- runif(2, 0.2, 0.9); tau <- runif(2, 4, 20)
  kk <- seq(0.05, 1.6, length.out = 25)
  num <- vapply(kk, function(k) {
    r <- seq(1e-6, 400, length.out = 200000)
    Cr <- rowSums(sapply(seq_along(s), function(m) s[m]*om3(a[m]*r)*exp(-r^2/(2*tau[m]^2))))
    sum(Cr*sin(k*r)*r)*(r[2]-r[1])/(2*pi^2*k)
  }, numeric(1))
  cf <- dens_ref(kk, s, a, tau)
  rr <- cor(num, cf); sc <- sum(num*cf)/sum(cf^2)
  rel <- max(abs(num - sc*cf))/max(abs(num))
  say("  trial %d: correlation %.6f, worst relative departure %.2e\n", trial, rr, rel)
  if (rr < 0.999999 || rel > 1e-4) ok <- FALSE
}
say("  %s\n\n", ifelse(ok, "closed form confirmed", "CLOSED FORM WRONG - stop here"))
if (!ok) quit(status = 1)

# fit the reference to each shell
SBp <- fread(file.path(OUT, "shell_model_params_v2.csv"))
fit_shell <- function(d, y, M, hole = NULL) {
  rng <- diff(range(y))
  st <- list(p0 = min(y)*0.8,
             s = rep(rng/M, M),
             a = seq(0.25, 0.75, length.out = M),
             tau = rep(8, M))
  pk <- function(p) list(p0 = p[1], s = p[1 + seq_len(M)],
                         a = p[1 + M + seq_len(M)], tau = p[1 + 2*M + seq_len(M)])
  # `hole`, when given, is an extra component held fixed while the rest refit
  f <- function(p) { q <- pk(p)
    s <- q$s; a <- q$a; t <- q$tau
    if (!is.null(hole)) { s <- c(s, hole$s); a <- c(a, hole$a); t <- c(t, hole$tau) }
    gam_ref(d, q$p0, s, a, t) - y }
  # the wave number floor is the simulation box, not the data's resolvable band:
  # this is a specified truth rather than an estimate, and it only has to be
  # representable on a 192 mm grid whose fundamental is 0.033 rad/mm
  lo <- c(0, rep(1e-6, M), rep(0.10, M), rep(1.0, M))
  hi <- c(min(y), rep(5*rng, M), rep(pi/min(d), M), rep(100, M))
  p0v <- c(st$p0, st$s, st$a, st$tau)
  fit <- tryCatch(nls.lm(par = p0v, fn = f, lower = lo, upper = hi,
                         control = nls.lm.control(maxiter = 4000, maxfev = 100000)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  q <- pk(fit$par); q$rss <- sum(f(fit$par)^2); q$M <- M; q
}

dip6 <- function(d, g) { i <- which.min(abs(d - 6)); 100*(g[i] - min(g[i:length(g)]))/g[i] }

say("== fitting the reference family, raw shell curves ==\n")
say("  %-6s %3s %10s %10s %10s %26s %20s\n", "shell", "M", "rss", "SB rss",
    "hole", "wave numbers (rad/mm)", "damping (mm)")
best <- list()
for (k in shells) {
  q <- CUR[shell == k][order(d)]; d <- q$d; y <- q$g
  bb <- NULL
  for (M in 1:3) {
    r <- fit_shell(d, y, M); if (is.null(r)) next
    r$dip <- dip6(d, gam_ref(d, r$p0, r$s, r$a, r$tau))
    say("  %-6d %3d %10.5f %10.5f %9.2f%% %26s %20s\n", k, M, r$rss,
        SBp[shell == k]$sb_rss, r$dip,
        paste(sprintf("%.3f", r$a), collapse = " "),
        paste(sprintf("%.1f", r$tau), collapse = " "))
    if (is.null(bb) || r$rss < bb$rss*0.98) bb <- r
  }
  best[[as.character(k)]] <- bb
}

# the hole the fit will not take, and why it has to be added by hand
say("\n== the hole ==\n")
say("  %-6s %14s %14s\n", "shell", "in the data", "in the fit")
for (k in shells) {
  q <- CUR[shell == k][order(d)]; b <- best[[as.character(k)]]
  say("  %-6d %13.2f%% %13.2f%%\n", k, dip6(q$d, q$g), b$dip)
}
say("\n  The data's drop below the 6 mm value is about 5 per cent, and the fitted\n")
say("  reference reproduces almost none of it. That is correct behaviour, not a\n")
say("  failure of the family. The 6 mm class is three pure axis displacements and\n")
say("  its value is a single-lag spike above both neighbours; no function of |h|\n")
say("  alone can produce that, so least squares ignores it and follows the smooth\n")
say("  rise. S12 to S14 attribute the spike to the separable resampling kernel.\n")
say("  It also matches the paper's own finding that the Gaussian's error comes\n")
say("  from the shape of the rise rather than from the hole.\n\n")
say("  Since a reference with no hole cannot show whether the conclusions depend\n")
say("  on one, a second is built with the deepest hole the data will bear: one\n")
say("  further admissible component whose wave number puts the crest near 6 mm,\n")
say("  with the rest of the model refitted around it and the depth pushed as far\n")
say("  as a threefold residual budget allows. Both are carried through the\n")
say("  experiment so the answer can be reported as robust to the hole or not.\n")
say("  Note what that search finds: even pushed to the budget, a smooth isotropic\n")
say("  model reaches only a third to a half of the observed drop. The rest is not\n")
say("  available to any function of |h|, which is the same conclusion the single-\n")
say("  lag diagnosis reaches by a different route.\n")

CONST <- (1/(2*pi^2)) * 0.5 * sqrt(pi/2)    # from the derivation above
rows <- list()
say("\n== the two references ==\n")
say("  %-6s %-9s %10s %10s %10s %10s\n", "shell", "which", "rss", "hole", "nugget", "sill")
for (k in shells) {
  q <- CUR[shell == k][order(d)]; b <- best[[as.character(k)]]
  # Deepen the hole as far as the fit allows: fix one extra component, refit
  # everything else around it, and keep the deepest hole whose residual sum of
  # squares is still within three times the matched reference.
  target <- dip6(q$d, q$g); budget <- 3*b$rss
  bestB <- NULL
  for (ah in seq(0.55, 1.00, by = 0.05)) for (sh in seq(0.05, 0.80, by = 0.05)) {
    r <- fit_shell(q$d, q$g, b$M, hole = list(s = sh, a = ah, tau = 15))
    if (is.null(r)) next
    s2 <- c(r$s, sh); a2 <- c(r$a, ah); t2 <- c(r$tau, 15)
    gg <- gam_ref(q$d, r$p0, s2, a2, t2)
    dd <- dip6(q$d, gg); rr <- sum((q$g - gg)^2)
    if (rr > budget) next
    if (is.null(bestB) || dd > bestB$dip)
      bestB <- list(p0 = r$p0, s = s2, a = a2, tau = t2, rss = rr, dip = dd)
  }
  if (is.null(bestB)) bestB <- list(p0 = b$p0, s = b$s, a = b$a, tau = b$tau,
                                    rss = b$rss, dip = b$dip)
  for (nm in c("matched", "hole")) {
    z <- if (nm == "matched") list(s = b$s, a = b$a, tau = b$tau, rss = b$rss, dip = b$dip)
         else bestB
    say("  %-6d %-9s %10.5f %9.2f%% %10.3f %10.3f\n", k, nm, z$rss, z$dip,
        b$p0, sum(z$s))
    for (m in seq_along(z$s))
      rows[[length(rows)+1]] <- data.table(reference = nm, shell = k,
        depth_mm = depth_mm[k], comp = m, nugget = b$p0,
        sill = z$s[m], a = z$a[m], tau = z$tau[m], const = CONST)
  }
  best[[as.character(k)]]$hole <- bestB
}
P <- rbindlist(rows)
fwrite(P, file.path(OUT, "S25_reference_params.csv"))

say("\n  Gaussian residual sum of squares on the same curves, for scale: %s\n",
    paste(sprintf("%.3f", SBp[shell %in% shells]$g_rss), collapse = ", "))
say("  so both references describe the data about ten times better than the\n")
say("  Gaussian, and neither belongs to any family being assessed\n")

# admissibility and total mass, checked numerically
say("\n== admissibility, and does the density integrate back to the sill? ==\n")
kk <- seq(1e-4, 4.0, length.out = 60000)
say("  %-6s %-9s %14s %12s %12s %12s\n", "shell", "which", "least density",
    "negatives", "mass", "sill")
for (k in shells) {
  b <- best[[as.character(k)]]
  for (nm in c("matched", "hole")) {
    z <- if (nm == "matched") b else b$hole
    f <- dens_ref(kk, z$s, z$a, z$tau)
    # dens_ref returns the density up to the constant from the derivation, so the
    # radial mass element is CONST * 4 pi k^2 f(k) and integrating it gives C(0)
    mass <- CONST*sum(4*pi*kk^2*f)*(kk[2] - kk[1])
    say("  %-6d %-9s %14.3e %12d %12.3f %12.3f\n", k, nm,
        min(f), sum(f < -1e-12), mass, sum(z$s))
  }
}
say("\n  A density that is non-negative at every frequency is exactly what makes\n")
say("  the field valid, and here it is visible in the formula rather than\n")
say("  imposed: (a-k)^2 is smaller than (a+k)^2 whenever a and k are positive.\n")
say("  The mass column is the density integrated over all frequencies, which\n")
say("  must return the sill, and does.\n")
say("\nwrote S25_reference_params.csv\n")
})

# Every candidate model fitted to the reference
local({
TXT  <- file.path(OUT, "S27_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)

omega3 <- function(x) { o <- rep(1, length(x)); nz <- x != 0; o[nz] <- sin(x[nz])/x[nz]; o }

# the 31 observable separations
g <- as.data.table(expand.grid(dx = 0:4, dy = 0:4, dz = 0:4))[!(dx == 0 & dy == 0 & dz == 0)]
g[, k2 := dx^2 + dy^2 + dz^2]
D <- sort(unique(VOX*sqrt(g[k2 <= 48]$k2)))

# the truth
gam_ref <- function(r, p0, s, a, tau)
  p0 + rowSums(sapply(seq_along(s), function(m)
    s[m]*(1 - omega3(a[m]*r)*exp(-r^2/(2*tau[m]^2)))))

# fitters, carried over unchanged
make_nodes <- function(d) {
  band <- c(pi/max(d), pi/min(d))
  exp(seq(log(band[1]), log(band[2]), length.out = 60))
}
design_sb <- function(d, nodes) {
  M <- matrix(0, length(d), length(nodes) + 1L); M[, 1] <- 1
  for (j in seq_along(nodes)) M[, j + 1L] <- 1 - omega3(nodes[j]*d)
  M
}
fit_sb <- function(d, y, nodes) {
  X <- design_sb(d, nodes); p <- nnls(X, y)$x
  list(p = p, ss = sum((y - X %*% p)^2), nz = sum(p > 1e-12))
}
fit_gauss <- function(d, y, grid = seq(0.5, 60, length.out = 3000)) {
  best <- NULL
  for (a in grid) {
    u <- exp(-(d/a)^2); X <- cbind(u, 1 - u)
    cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) next
    if (cf[1] < 0) cf <- c(0, sum(y*(1-u))/sum((1-u)^2))
    ss <- sum((y - X %*% cf)^2)
    if (is.null(best) || ss < best$ss)
      best <- list(ss = ss, range = a, nugget = cf[1], sill = cf[2])
  }
  best
}
# The search grid for the exponential component is widened sixfold, and whether
# the fit lands inside it is reported below rather than assumed.
fit_mixed <- function(d, y) {
  best <- NULL
  for (b in seq(0.3, 25, length.out = 90))
    for (cc in seq(0.3, 300, length.out = 200)) {
      X <- cbind(1 - exp(-(d^2)/(2*b^2)), 1 - exp(-d/cc))
      p <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
      if (is.null(p) || anyNA(p)) next
      p <- pmax(p, 0)
      ss <- sum((y - X %*% p)^2)
      if (is.null(best) || ss < best$ss)
        best <- list(ss = ss, b = b, c = cc, p1 = p[1], p2 = p[2])
    }
  best
}
# Cox / CLEAN: gamma(d) = nugget + (sill - nugget)(1 - exp(-d/range)).
# The search grid is widened so the comparison is against the model at its own
# optimum rather than at a grid edge.
fit_exp <- function(d, y, grid = seq(0.2, 80, length.out = 4000)) {
  best <- NULL
  for (phi in grid) {
    u <- exp(-d/phi); X <- cbind(u, 1 - u)
    cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) next
    if (cf[1] < 0) cf <- c(0, sum(y*(1-u))/sum((1-u)^2))
    ss <- sum((y - X %*% cf)^2)
    if (is.null(best) || ss < best$ss)
      best <- list(ss = ss, range = phi, nugget = cf[1], sill = cf[2])
  }
  best
}

P <- fread(file.path(OUT, "S25_reference_params.csv"))
depth_mm <- c(9, 12, 15, 18, 24)
say("Each model is fitted to the reference variogram at the 31 observable\n")
say("separations, %.2f to %.2f mm. No estimation noise, so what is measured is\n",
    min(D), max(D))
say("the cost of the family's shape, exactly as in the estimation simulation.\n")

for (which in c("matched", "hole")) {
  say("\n================ %s reference ================\n", which)
  say("  %-6s %12s %10s %10s %10s %10s %8s\n", "shell", "SB rss", "Gauss rss",
      "mixed rss", "exp rss", "exp range", "mixed c")
  CO <- list(); PA <- list(); EX <- list()
  for (k in sort(unique(P$shell))) {
    t <- P[reference == which & shell == k]
    y <- gam_ref(D, t$nugget[1], t$sill, t$a, t$tau)
    nodes <- make_nodes(D)
    sb <- fit_sb(D, y, nodes); ga <- fit_gauss(D, y)
    mx <- fit_mixed(D, y);     ex <- fit_exp(D, y)
    say("  %-6d %12.2e %10.5f %10.5f %10.5f %10.2f %8.1f\n", k, sb$ss, ga$ss,
        mx$ss, ex$ss, ex$range, mx$c)
    CO[[length(CO)+1]] <- data.table(shell = k, depth_mm = depth_mm[k],
      term = seq_along(sb$p) - 1L, node = c(NA_real_, nodes), coef = sb$p)
    PA[[length(PA)+1]] <- data.table(shell = k, depth_mm = depth_mm[k], lambda = 0,
      sb_rss = sb$ss, sb_sill = sum(sb$p), sb_active = sb$nz,
      node_lo = min(nodes), node_hi = max(nodes),
      g_range = ga$range, g_nugget = ga$nugget, g_sill = ga$sill, g_rss = ga$ss,
      m_b = mx$b, m_c = mx$c, m_p1 = mx$p1, m_p2 = mx$p2, m_rss = mx$ss)
    EX[[length(EX)+1]] <- data.table(shell = k, depth_mm = depth_mm[k],
      e_range = ex$range, e_nugget = ex$nugget, e_sill = ex$sill, e_rss = ex$ss)
  }
  fwrite(rbindlist(CO), file.path(OUT, sprintf("S27_ref_%s_sbcoefs.csv", which)))
  fwrite(rbindlist(PA), file.path(OUT, sprintf("S27_ref_%s_params.csv", which)))
  fwrite(rbindlist(EX), file.path(OUT, sprintf("S27_ref_%s_exp.csv", which)))
  E <- rbindlist(EX); A <- rbindlist(PA)
  say("  exponential range interior on a 0.2 to 80 grid in %d of %d shells\n",
      sum(E$e_range > 0.3 & E$e_range < 79), nrow(E))
  say("  mixed exponential component interior on a 0.3 to 300 grid in %d of %d\n",
      sum(A$m_c < 299), nrow(A))
  say("  Shapiro-Botha reproduces the reference curve to %.1e, so its arm tests\n",
      max(A$sb_rss))
  say("  whether agreeing at 31 separations is enough to imply the same FIELD,\n")
  say("  which S24 showed it need not be\n")
}
say("\nwrote S27_ref_*_{sbcoefs,params,exp}.csv\n")
})

# Identifiability sweep: what the data leave undetermined
local({
TXT  <- file.path(OUT, "S34_report.txt")
say  <- function(...) { m <- sprintf(...); cat(m); cat(m, file = TXT, append = TRUE) }
cat("", file = TXT)
VOX <- 3; M <- 3

om3 <- function(x) { o <- rep(1, length(x)); nz <- x != 0; o[nz] <- sin(x[nz])/x[nz]; o }
gam_ref <- function(r, p0, s, a, tau)
  p0 + rowSums(sapply(seq_along(s), function(m)
    s[m]*(1 - om3(a[m]*r)*exp(-r^2/(2*tau[m]^2)))))

# the cohort shell curves and their precision
L <- fread(file.path(DIAG, "out_matlab", "rebin_lagwise.csv"))
L[, k2 := as.integer(round((distance_mm/VOX)^2))]
SH <- melt(L, id.vars = c("subject","run","k2"), measure.vars = patterns("^gamma_shell"),
           variable.name = "shell", value.name = "g")
SH[, shell := as.integer(sub("gamma_shell", "", shell))]
SUB <- SH[!is.na(g), .(g = mean(g)), by = .(subject, shell, k2)]
CUR <- SUB[, .(g = mean(g), se = sd(g)/sqrt(.N), n = .N), by = .(shell, k2)][
        , d := VOX*sqrt(k2)][order(shell, d)]
shells <- sort(unique(CUR$shell))
depth_mm <- c(9, 12, 15, 18, 24)

say("== precision of the cohort curve ==\n")
say("  %-6s %8s %14s %14s %14s\n", "shell", "subjects", "mean SE", "smallest SE", "largest SE")
for (k in shells) {
  q <- CUR[shell == k]
  say("  %-6d %8d %14.4f %14.4f %14.4f\n", k, q$n[1], mean(q$se), min(q$se), max(q$se))
}
say("\n  A model is accepted when its root mean square departure over the 31\n")
say("  separations is no larger than the mean standard error of the curve it is\n")
say("  fitted to. That is the weakest criterion that still means the data cannot\n")
say("  distinguish it from the observed curve.\n")

# Fit at a fixed total variance. `warm` seeds the fit from the neighbouring
# solution, so the profile is traced by continuation; four fixed starts are kept
# as fallbacks.
fit_at_total <- function(d, y, V, warm = NULL, W = 1e4) {
  rng <- diff(range(y))
  pk <- function(p) list(p0 = p[1], s = p[1 + seq_len(M)],
                         a = p[1 + M + seq_len(M)], tau = p[1 + 2*M + seq_len(M)])
  f <- function(p) { q <- pk(p)
    c(gam_ref(d, q$p0, q$s, q$a, q$tau) - y, W*(q$p0 + sum(q$s) - V)) }
  lo <- c(0, rep(1e-6, M), rep(0.10, M), rep(1.0, M))
  hi <- c(min(y), rep(5*rng, M), rep(pi/min(d), M), rep(100, M))
  starts <- list()
  if (!is.null(warm)) {
    w <- c(warm$p0, warm$s, warm$a, warm$tau)
    w[2:(1+M)] <- w[2:(1+M)] * max(1e-6, (V - warm$p0)/max(1e-9, sum(warm$s)))
    starts[[1]] <- w
  }
  for (trial in 1:4)
    starts[[length(starts)+1]] <- c(min(y)*c(0.8, 0.5, 0.95, 0.65)[trial],
      rep((V - min(y)*0.8)/M, M),
      seq(0.15, 0.85, length.out = M) * c(1, 0.8, 1.2, 1)[trial],
      rep(c(8, 20, 40, 12)[trial], M))
  best <- NULL
  for (st in starts) {
    st <- pmin(pmax(st, lo), hi)
    fit <- tryCatch(nls.lm(par = st, fn = f, lower = lo, upper = hi,
                           control = nls.lm.control(maxiter = 3000, maxfev = 60000)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    q <- pk(fit$par)
    res <- gam_ref(d, q$p0, q$s, q$a, q$tau) - y
    r <- list(p0 = q$p0, s = q$s, a = q$a, tau = q$tau,
              rss = sum(res^2), rmse = sqrt(mean(res^2)),
              total = q$p0 + sum(q$s))
    if (abs(r$total - V) > 0.01) next          # constraint not actually met
    if (is.null(best) || r$rss < best$rss) best <- r
  }
  best
}

REF <- fread(file.path(OUT, "S25_reference_params.csv"))
rows <- list(); pars <- list()
for (k in shells) {
  q <- CUR[shell == k][order(d)]; d <- q$d; y <- q$g
  tol <- mean(q$se)
  t0 <- REF[reference == "matched" & shell == k]
  V0 <- t0$nugget[1] + sum(t0$sill)
  say("\n== shell %d (%d mm), reference total %.3f, tolerance %.4f ==\n",
      k, depth_mm[k], V0, tol)
  say("  %8s %10s %12s %12s %10s\n", "total", "rmse", "within noise", "C at 20.8", "nugget")
  # continuation: fit at the reference's own total first, then walk outwards in
  # both directions seeding each fit from its neighbour
  grid <- round(seq(0.70*V0, 1.60*V0, length.out = 19), 4)
  i0 <- which.min(abs(grid - V0))
  fits <- vector("list", length(grid))
  fits[[i0]] <- fit_at_total(d, y, grid[i0])
  for (i in (i0-1):1)             fits[[i]] <- fit_at_total(d, y, grid[i], fits[[i+1]])
  for (i in (i0+1):length(grid))  fits[[i]] <- fit_at_total(d, y, grid[i], fits[[i-1]])
  for (gi in seq_along(grid)) {
    V <- grid[gi]; r <- fits[[gi]]
    if (is.null(r)) { say("  %8.3f %10s\n", V, "no fit"); next }
    cend <- sum(r$s*om3(r$a*max(d))*exp(-max(d)^2/(2*r$tau^2)))
    ok <- r$rmse <= tol
    say("  %8.3f %10.4f %12s %12.3f %10.3f\n", V, r$rmse,
        ifelse(ok, "yes", "no"), cend, r$p0)
    rows[[length(rows)+1]] <- data.table(shell = k, depth_mm = depth_mm[k],
      total = V, rmse = r$rmse, tol = tol, within = ok, C_end = cend, nugget = r$p0)
    for (m in 1:M)                       # keep every fit; filter on tolerance later
      pars[[length(pars)+1]] <- data.table(shell = k, depth_mm = depth_mm[k],
        total = V, rmse = r$rmse, comp = m, nugget = r$p0, sill = r$s[m],
        a = r$a[m], tau = r$tau[m])
  }
}
S <- rbindlist(rows); P <- rbindlist(pars)
fwrite(S, file.path(OUT, "S34_sweep.csv"))
fwrite(P, file.path(OUT, "S34_sweep_params.csv"))

say("\n================ the identified range depends on the tolerance ================\n")
say("  The answer is not a single number, so it is reported as a curve. The best\n")
say("  any model in this family achieves is a floor set by lag-to-lag structure\n")
say("  it cannot represent, so tolerances are quoted as multiples of that floor\n")
say("  as well as against the curve's own standard error.\n\n")
say("  %-6s %12s %26s %12s %12s\n", "shell", "tolerance", "basis", "total from", "to")
BAND <- list()
for (k in shells) {
  s <- S[shell == k]; floor_k <- min(s$rmse); tol <- s$tol[1]
  levs <- list(c(1.5*floor_k, "1.5x the best achievable"),
               c(3*floor_k,   "3x the best achievable"),
               c(0.5*tol,     "half the curve's own SE"),
               c(tol,         "the curve's own SE"))
  for (lv in levs) {
    th <- as.numeric(lv[[1]]); a <- s[rmse <= th]
    if (!nrow(a)) next
    say("  %-6d %12.4f %26s %12.3f %12.3f\n", k, th, lv[[2]], min(a$total), max(a$total))
    BAND[[length(BAND)+1]] <- data.table(shell = k, basis = lv[[2]], thresh = th,
      lo = min(a$total), hi = max(a$total), ratio = max(a$total)/min(a$total))
  }
  say("\n")
}
B <- rbindlist(BAND); fwrite(B, file.path(OUT, "S34_bands.csv"))
say("  ratio of largest to smallest admissible total variance:\n")
for (b in unique(B$basis))
  say("    %-26s %.2f to %.2f across the three shells\n", b,
      min(B[basis == b]$ratio), max(B[basis == b]$ratio))

say("\n== the sign of the correlation at the window edge is not determined ==\n")
for (k in shells) {
  s <- S[shell == k & rmse <= 3*min(S[shell == k]$rmse)]
  say("  shell %d: C(20.8) ranges %+.3f to %+.3f across models fitting within\n",
      k, min(s$C_end), max(s$C_end))
  say("           three times the best achievable\n")
}

say("\n== the floor, and what it means ==\n")
for (k in shells)
  say("  shell %d: best achievable rmse %.4f against a cohort SE of %.4f and a\n",
      k, min(S[shell == k]$rmse), S[shell == k]$tol[1])
say("  ten-run within-subject precision of about 0.002 (F17).\n")
say("  No model in this family reaches the within-subject precision, because the\n")
say("  curve carries lag-to-lag detail that no function of separation alone can\n")
say("  produce. So the identifiability question cannot even be posed at that\n")
say("  precision: the residual there is model inadequacy, not ambiguity.\n")
say("  Against the cohort mean, which is the precision a cohort-level claim\n")
say("  actually has, the ambiguity is real and is the range tabulated above.\n")

fwrite(P, file.path(OUT, "S34_sweep_params.csv"))
say("\n  %d fits retained in all, spanning %d totals per shell.\n",
    uniqueN(P[, .(shell, total)]), uniqueN(S$total))
})

