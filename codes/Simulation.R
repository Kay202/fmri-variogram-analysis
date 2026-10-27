# Synthetic-truth simulation: integrated mean squared error of five variogram
# estimators against curves whose shape is known.
#
# Nothing here reads the imaging data, any fitted curve, or any measured
# quantity. The lag grid comes from the voxel size and the lag radius; every
# prior range is argued from lattice geometry, the interpolation kernel, or
# identifiability.
#
# Design
#   lags        31 distinct separations reachable within a 4-voxel radius at
#               3 mm spacing, 3.000 to 20.785 mm, generated here
#   families    DoubleExp and Matern, differing only in the shape of the rise.
#               A curve index fixes every other parameter, including the sill
#               fraction, so the contrast between families is clean
#   curves      120 per family, drawn independently. Between-curve variance
#               exceeds within-curve Monte Carlo variance by about 100 to 1, so
#               the budget buys sample size rather than replicates
#   rise        uniform on the sill fraction reached at the last lag, inverted
#               per family so both span the same region. The lower bound of 0.5
#               is an identifiability statement: below it the sill cannot be
#               recovered from these lags
#   nugget      U(0, 1), parity with the structural sill
#   amplitude   U(0, 0.25); zero included so the parametric benchmark's
#               monotonicity can be tested, ceiling set by the side-lobe ratio
#               of the interpolation kernel
#   frequency   U(2, 20), a third of a cycle to a little over three
#   damping     U(0.15, 1.00), dying within one voxel to persisting across
#               the window
#   noise       Gaussian, chi-squared and lognormal at six coefficients of
#               variation to 0.30. A Gamma arm was dropped as algebraically
#               identical to the chi-squared one
#   replicates  200. Squared bias estimated from finite replicates is inflated
#               by IVar/B, and the correction is applied explicitly below
#   admissible  checked as a rejection step inside the draw loop, with the
#               rejection count reported
#
# The Shapiro-Botha fit is plain non-negative least squares, as published, with
# complexity controlled by the node set rather than a penalty. The node set spans
# exactly the band these separations resolve, with no margin above: an
# out-of-band node is flat across the sampled range and the fit uses it as a
# nugget substitute, which changes the implied random field without changing the
# error measured at the sampled lags.
#
# Usage
#   Rscript Simulation.R bank                build the curve bank and stop
#   Rscript Simulation.R bench [workers]     timing benchmark, no output files
#   Rscript Simulation.R run   [workers]     full run, checkpointed, resumable
#
# Writes SIM2_metrics.csv, SIM2_curves.rds and its checkpoint files to the
# output directory set in the preamble below.

suppressPackageStartupMessages({ library(data.table); library(parallel) })
argv <- commandArgs(trailingOnly = FALSE)
fa   <- sub("^--file=", "", argv[grep("^--file=", argv)])
HERE <- if (length(fa)) dirname(normalizePath(fa)) else getwd()
ta   <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(ta) >= 1) ta[1] else "bench"
NW   <- if (length(ta) >= 2) as.integer(ta[2]) else max(1L, detectCores(logical = FALSE) - 1L)
REL  <- file.path(HERE, "helpers")
OUT  <- file.path(HERE, "..", "data", "intermediates")

# design constants
VOX      <- 3          # mm, isotropic
KMAX     <- 4          # lag radius in voxels
B        <- 200L
N_CURVE  <- 120L
FAMILIES <- c("DoubleExp", "Matern")
NOISE    <- c("Gaussian", "ChiSquared", "Lognormal")
CVS      <- seq(0.05, 0.30, by = 0.05)
METHODS  <- c("LocalCubic", "PTR", "NatCubicSpline", "ShapiroBotha", "Gaussian")
SEED     <- 20260910L

# Shapiro-Botha non-negative spectral mixture, fitted by non-negative least
# squares over a fixed node set. There is no penalty term: complexity is set by
# the node set, which is defined below because it depends on the lag grid.

# lag grid, from the lattice alone
gg <- expand.grid(dx = 0:KMAX, dy = 0:KMAX, dz = 0:KMAX)
kk <- sort(unique(with(gg, dx^2 + dy^2 + dz^2))); kk <- kk[kk > 0]
D_MM <- VOX * sqrt(kk)
NB   <- length(D_MM)
TG   <- (D_MM - min(D_MM)) / (max(D_MM) - min(D_MM))

# The node set spans exactly the band the lag grid resolves, from pi over the
# largest separation to pi over the smallest, with no margin above. Past that
# band the basis function is flat across the sampled range and the fit cannot
# separate such a node from the nugget.
SB_BAND  <- c(pi/max(D_MM), pi/min(D_MM))
SB_NODES <- exp(seq(log(SB_BAND[1]), log(SB_BAND[2]), length.out = 60))
omega3 <- function(x) { o <- rep(1, length(x)); nz <- x != 0
                        o[nz] <- sin(x[nz])/x[nz]; o }
design_sb <- function(dd) {
  M <- matrix(0, length(dd), length(SB_NODES) + 1L); M[, 1] <- 1
  for (j in seq_along(SB_NODES)) M[, j + 1L] <- 1 - omega3(SB_NODES[j]*dd)
  M
}
safe_sb <- function(dd, yy, n_expected) {
  res <- tryCatch(as.numeric(design_sb(dd) %*% nnls::nnls(design_sb(dd), yy)$x),
                  error = function(e) NA_real_)
  if (length(res) == n_expected && all(is.finite(res)))
    list(fit = res, status = "ok") else list(fit = rep(NA_real_, n_expected), status = "bad")
}

# truth families
matern25 <- function(t, rho) { u <- sqrt(5)*t/rho; (1 + u + u^2/3) * exp(-u) }
# invert the sill fraction reached at t = 1 into each family's rise parameter
a_from_s   <- function(s) -1/log(1 - s)
rho_from_s <- function(s) uniroot(function(r) (1 - matern25(1, r)) - s,
                                  c(1e-3, 50), tol = 1e-10)$root
curve_of <- function(p, fam) {
  rise <- if (fam == "DoubleExp") p$sigma2 * (1 - exp(-TG/p$a))
          else                    p$sigma2 * (1 - matern25(TG, p$rho))
  pmax(rise + p$A * cos(p$omega * TG) * exp(-TG/p$b) + p$tau, 1e-6)
}

# admissibility in R^3, nugget included
KGRID <- seq(1e-3, 600, length.out = 40000)
cos_term <- function(k, b, w) b*(k+w)/(b^2+(k+w)^2)^2 + b*(k-w)/(b^2+(k-w)^2)^2
CM <- function(rho) gamma(4)*5^2.5/(pi^1.5*gamma(2.5)*rho^5)
admissible <- function(p, fam) {
  be <- 1/p$b
  f <- if (fam == "DoubleExp") {
    al <- 1/p$a
    (p$sigma2*2*al*KGRID/(al^2+KGRID^2)^2 - p$A*cos_term(KGRID, be, p$omega))/(2*pi^2*KGRID)
  } else {
    p$sigma2*CM(p$rho)*(5/p$rho^2+KGRID^2)^-4 - p$A*cos_term(KGRID, be, p$omega)/(2*pi^2*KGRID)
  }
  min(f) + (p$A + p$tau)/(2*pi)^3 >= 0
}

# prior, drawn independently rather than stratified, so that any initial segment
# of the sequence is an ordinary sample
PRIOR <- list(s = c(0.50, 0.999), sigma2 = c(0.85, 1.05), A = c(0.00, 0.25),
              omega = c(2, 20), b = c(0.15, 1.00), tau = c(0.00, 1.00))
draw_one <- function() lapply(PRIOR, function(rg) runif(1, rg[1], rg[2]))

build_bank <- function() {
  set.seed(SEED)
  rows <- list(); rej <- 0L
  while (length(rows) < N_CURVE) {
    p <- draw_one()
    p$a <- a_from_s(p$s); p$rho <- rho_from_s(p$s)
    # rejection leaves the accepted draws independent, from the prior
    # conditioned on admissibility, which is the prior we intend
    if (!admissible(p, "DoubleExp") || !admissible(p, "Matern")) { rej <- rej + 1L; next }
    rows[[length(rows)+1]] <- p
  }
  list(p = rows, rejected = rej)
}

# noise, three distinct laws
draw_noise <- function(model, g, cv) switch(model,
  Gaussian   = g + rnorm(length(g), 0, cv*g),
  ChiSquared = { nu <- 2/cv^2; g*rchisq(length(g), df = nu)/nu },
  Lognormal  = { s2 <- log(1+cv^2); rlnorm(length(g), log(g) - s2/2, sqrt(s2)) },
  stop("unknown noise model"))

# The parametric benchmark. Bounds are only those the model requires, set far
# outside any plausible optimum, and are the same ones used to fit the benchmark
# to the residuals.
# Fitted by profiling rather than iteration. The Gaussian variogram is linear in
# nugget and sill once the range is fixed, so the range is scanned and the other
# two solved exactly at each value. The scan is coarse-to-fine: 60 log-spaced
# values, then 40 between the neighbours of the coarse optimum.
gauss_at_range <- function(d, y, a) {
  u <- exp(-(d/a)^2); X <- cbind(u, 1 - u)
  cf <- tryCatch(as.vector(qr.solve(X, y)), error = function(e) NULL)
  if (is.null(cf) || anyNA(cf)) return(NULL)
  if (cf[1] < 0) cf <- c(0, sum(y*(1-u))/sum((1-u)^2))   # nugget cannot be negative
  list(ss = sum((y - X %*% cf)^2), range = a,
       nugget = cf[1], sill = cf[2], fit = as.vector(X %*% cf))
}
GR_LO <- 0.5; GR_HI <- 60
safe_gauss_wide <- function(d, y, n_expected) {
  scan <- function(grid) {
    best <- NULL
    for (a in grid) { r <- gauss_at_range(d, y, a)
      if (!is.null(r) && (is.null(best) || r$ss < best$ss)) best <- r }
    best
  }
  co <- exp(seq(log(GR_LO), log(GR_HI), length.out = 60))
  b1 <- scan(co)
  if (is.null(b1)) return(list(fit = rep(NA_real_, n_expected), atbound = 0L))
  i  <- which.min(abs(co - b1$range))
  b  <- scan(exp(seq(log(co[max(1, i-1)]), log(co[min(60, i+1)]), length.out = 40)))
  if (is.null(b)) b <- b1
  # the range resting on the edge of its scan, or a sill below the nugget, both
  # mean the fitted object is not a Gaussian variogram; counted, not hidden
  ab <- as.integer(b$range <= GR_LO*1.001 || b$range >= GR_HI*0.999 || b$sill < b$nugget)
  if (length(b$fit) == n_expected && all(is.finite(b$fit)))
    list(fit = b$fit, atbound = ab)
  else list(fit = rep(NA_real_, n_expected), atbound = ab)
}

# one cell
run_cell <- function(task) {
  y <- task$y
  M <- lapply(METHODS, function(m) matrix(NA_real_, B, NB)); names(M) <- METHODS
  atbound <- 0L
  for (r in seq_len(B)) {
    ys <- draw_noise(task$noise, y, task$cv)
    M$LocalCubic[r, ]     <- safe_lp(TG,   ys, NB)$fit
    M$PTR[r, ]            <- safe_ptr(TG,  ys, NB)$fit
    M$NatCubicSpline[r, ] <- safe_ncs(TG,  ys, NB)$fit
    M$ShapiroBotha[r, ]   <- safe_sb(D_MM, ys, NB)$fit
    g <- safe_gauss_wide(D_MM, ys, NB)
    M$Gaussian[r, ] <- g$fit; atbound <- atbound + g$atbound
  }
  out <- lapply(METHODS, function(m) {
    X <- M[[m]]; X <- X[complete.cases(X), , drop = FALSE]
    if (nrow(X) < 2) return(c(IMSE=NA, IBias2=NA, IBias2_c=NA, IVar=NA, n_ok=nrow(X)))
    bias <- colMeans(X) - y; v <- apply(X, 2, var)
    ib <- mean(bias^2)
    # squared bias from finite replicates is inflated by the sampling variance
    # over the number of replicates; the corrected value subtracts it
    c(IMSE = mean(bias^2 + v), IBias2 = ib,
      IBias2_c = max(ib - mean(v)/nrow(X), 0), IVar = mean(v), n_ok = nrow(X))
  })
  data.table(family = task$family, curve = task$curve, noise = task$noise,
             cv = task$cv, method = METHODS,
             IMSE = sapply(out, `[`, "IMSE"), IBias2 = sapply(out, `[`, "IBias2"),
             IBias2_c = sapply(out, `[`, "IBias2_c"), IVar = sapply(out, `[`, "IVar"),
             n_ok = sapply(out, `[`, "n_ok"),
             # constraint activity for the parametric benchmark, recorded so a
             # binding bound can never again pass unnoticed
             gauss_atbound = ifelse(METHODS == "Gaussian", atbound, NA_integer_))
}

# assemble the task list
bank <- build_bank()
CURVES <- list()
for (fam in FAMILIES) for (i in seq_len(N_CURVE))
  CURVES[[paste(fam, i)]] <- curve_of(bank$p[[i]], fam)

tasks <- list()
for (fam in FAMILIES) for (i in seq_len(N_CURVE))
  for (nm in NOISE) for (cv in CVS)
    tasks[[length(tasks)+1]] <- list(family = fam, curve = i, noise = nm, cv = cv,
                                     y = CURVES[[paste(fam, i)]])

cat(sprintf("lag grid   : %d classes, %.3f to %.3f mm, steps %.3f to %.3f\n",
            NB, min(D_MM), max(D_MM), min(diff(D_MM)), max(diff(D_MM))))
cat(sprintf("curve bank : %d per family, independent draws; %d rejected as inadmissible\n",
            N_CURVE, bank$rejected))
cat(sprintf("sill frac  : %.3f to %.3f | a %.3f to %.3f | rho %.3f to %.3f\n",
            min(sapply(bank$p, `[[`, "s")), max(sapply(bank$p, `[[`, "s")),
            min(sapply(bank$p, `[[`, "a")), max(sapply(bank$p, `[[`, "a")),
            min(sapply(bank$p, `[[`, "rho")), max(sapply(bank$p, `[[`, "rho"))))
cat(sprintf("tasks      : %d cells x %d replicates x %d methods = %s fits\n",
            length(tasks), B, length(METHODS),
            format(length(tasks)*B*length(METHODS), big.mark = ",")))

if (MODE == "bank") {
  P <- rbindlist(lapply(bank$p, as.data.table))
  cat("\nprior ranges as drawn:\n")
  for (k in names(P))
    cat(sprintf("  %-7s %8.3f to %8.3f   mean %8.3f\n", k, min(P[[k]]), max(P[[k]]), mean(P[[k]])))
  cat(sprintf("\nacceptance rate %.1f%% (%d rejected)\n",
              100*N_CURVE/(N_CURVE + bank$rejected), bank$rejected))
  quit(save = "no")
}

mk_cluster <- function(nw) {
  cl <- makeCluster(nw)
  clusterExport(cl, c("REL","TG","D_MM","NB","B","METHODS","draw_noise","run_cell",
                      "SB_NODES","omega3","design_sb",
                      "safe_sb","gauss_at_range","GR_LO","GR_HI","safe_gauss_wide"),
                envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(data.table); library(minpack.lm)
                                     library(nnls) })
    src <- function(p) { e <- parse(p)
      for (x in e) if (is.call(x) && length(x) >= 3 &&
          as.character(x[[1]])[1] %in% c("<-","=") && is.call(x[[3]]) &&
          identical(as.character(x[[3]][[1]])[1], "function")) eval(x, envir = globalenv()) }
    src(file.path(REL, "Main analysis_Model Fitting.R"))
    src(file.path(REL, "helpers_fit_gaussian.R"))
    src(file.path(REL, "Boostrap analysis.R"))
    NULL })
  cl
}

# benchmark
if (MODE == "bench") {
  suppressPackageStartupMessages({ library(minpack.lm); library(nnls) })
  src <- function(p) { e <- parse(p)
    for (x in e) if (is.call(x) && length(x) >= 3 &&
        as.character(x[[1]])[1] %in% c("<-","=") && is.call(x[[3]]) &&
        identical(as.character(x[[3]][[1]])[1], "function")) eval(x, envir = globalenv()) }
  src(file.path(REL, "Main analysis_Model Fitting.R"))
  src(file.path(REL, "helpers_fit_gaussian.R"))
  src(file.path(REL, "Boostrap analysis.R"))

  NS <- 8L
  set.seed(1); pick <- sample(seq_along(tasks), NS)
  cat(sprintf("\n--- serial benchmark, %d cells ---\n", NS))
  t0 <- Sys.time(); invisible(lapply(tasks[pick], run_cell))
  ser <- as.numeric(difftime(Sys.time(), t0, units = "secs"))/NS
  cat(sprintf("  %.2f s per cell\n", ser))
  cat(sprintf("  serial total for %d cells: %.2f hours\n",
              length(tasks), ser*length(tasks)/3600))

  for (nw in unique(c(4L, NW, 8L))) {
    NP <- nw*3L; set.seed(2); pick2 <- sample(seq_along(tasks), NP)
    cat(sprintf("\n--- parallel benchmark, %d workers, %d cells ---\n", nw, NP))
    cl <- mk_cluster(nw)
    t0 <- Sys.time(); invisible(parLapplyLB(cl, tasks[pick2], run_cell))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    stopCluster(cl)
    cat(sprintf("  %.2f s per cell wall clock | speedup %.2fx | efficiency %.0f%%\n",
                el/NP, ser/(el/NP), 100*(ser/(el/NP))/nw))
    cat(sprintf("  projected total: %.2f hours\n", (el/NP)*length(tasks)/3600))
  }
  quit(save = "no")
}

# full run
PARTIAL <- file.path(OUT, sprintf("SIM2_partial_N%d_B%d.rds", N_CURVE, B))
done <- if (file.exists(PARTIAL)) readRDS(PARTIAL) else list()
todo <- setdiff(seq_along(tasks), as.integer(names(done)))
cat(sprintf("\nresuming with %d of %d cells complete; %d workers\n",
            length(done), length(tasks), NW))

# the design travels with the curve bank so the reporting script never has to
# restate a constant and the two cannot drift apart
saveRDS(list(bank = bank, design = list(
  VOX = VOX, KMAX = KMAX, B = B, N_CURVE = N_CURVE,
  FAMILIES = FAMILIES, NOISE = NOISE, CVS = CVS, METHODS = METHODS,
  PRIOR = PRIOR, D_MM = D_MM, SEED = SEED)),
  file.path(OUT, "SIM2_curves.rds"))

cl <- mk_cluster(NW); on.exit(stopCluster(cl), add = TRUE)
CH <- split(todo, ceiling(seq_along(todo)/(NW*8)))
t0 <- Sys.time()
for (ci in seq_along(CH)) {
  res <- parLapplyLB(cl, tasks[CH[[ci]]], run_cell)
  names(res) <- as.character(CH[[ci]]); done <- c(done, res)
  saveRDS(done, PARTIAL)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[chunk %d/%d] %d cells done, %.1f min elapsed, ETA %.1f min\n",
              ci, length(CH), length(done), el/60,
              (el/length(unlist(CH[seq_len(ci)])))*(length(todo)-length(unlist(CH[seq_len(ci)])))/60))
}
R <- rbindlist(done)
P <- rbindlist(lapply(seq_len(N_CURVE), function(i) as.data.table(bank$p[[i]])))
P[, curve := seq_len(.N)]
fwrite(merge(R, P, by = "curve"), file.path(OUT, "SIM2_metrics.csv"))
cat(sprintf("\nwrote SIM2_metrics.csv (%d rows) -> %s\n", nrow(R), OUT))
