library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(gridExtra)
if (!require(minpack.lm)) install.packages("minpack.lm")
library(minpack.lm)

gsjs <- function(t, y) {
  n <- length(t); s2 <- 0; epstil <- matrix(0, n-2, 1)
  for (i in 2:(n-1)) {
    if ((t[i-1]==t[i])||(t[i]==t[i+1])) { A<-.5; B<-.5; den<-1.5 }
    else {
      den <- t[i+1]-t[i-1]
      A <- (t[i+1]-t[i])/den; B <- (t[i]-t[i-1])/den
      den <- sqrt(1+A^2+B^2)
    }
    epstil[i-1] <- (y[i]-A*y[i-1]-B*y[i+1])/den
    s2 <- s2+epstil[i-1]^2
  }
  list(epstil, s2/(n-2))
}

gcv_fn <- function(t, y, u) {
  n <- length(t); s2 <- gsjs(t,y)[[2]]; RSS <- u[[2]]
  list((RSS+2*s2*u[[3]])/n, RSS/(1-u[[3]]/n)^2/n)
}

search_fn <- function(t, y, mesh, f) {
  crit <- matrix(0, length(mesh), 2)
  for (i in seq_along(mesh)) {
    tmp <- gcv_fn(t, y, f(t, y, mesh[i]))
    crit[i,1] <- tmp[[1]]; crit[i,2] <- tmp[[2]]
  }
  crit
}

kvals <- function(u) {
  out <- matrix(0, length(u), 1)
  for (i in seq_along(u)) out[i] <- dnorm(u[i])
  out
}

locpoly3 <- function(t, y, lam, x.points) {
  n <- length(t); y <- as.matrix(y); nc <- ncol(y)
  neval <- length(x.points); fhat <- matrix(0, neval, nc)
  for (i in 1:neval) {
    u <- (x.points[i] - t) / lam; w <- as.vector(kvals(u)); sw <- sqrt(w)
    for (c in 1:nc) {
      yi <- y[, c]; X <- cbind(1, u, u^2, u^3); Xw <- X*sw; yw <- yi*sw
      beta <- tryCatch(
        qr.solve(Xw, yw),
        error = function(e) tryCatch(
          qr.solve(X[, 1:3, drop=FALSE] * sw, yw),
          error = function(e2) tryCatch(
            qr.solve(X[, 1:2, drop=FALSE] * sw, yw),
            error = function(e3) c(weighted.mean(yi, w + 1e-10), 0, 0, 0)
          )
        )
      )
      fhat[i, c] <- beta[1]
    }
  }
  list(fhat)
}

locpoly3_wrap <- function(t, y, lam) {
  n <- length(t); yhat <- locpoly3(t, y, lam, t)[[1]][, 1]
  RSS <- sum((y - yhat)^2); trH <- 0
  for (i in 1:n) {
    e_i <- rep(0, n); e_i[i] <- 1
    h_ii <- locpoly3(t, e_i, lam, t[i])[[1]][1, 1]
    trH <- trH + h_ii
  }
  list(as.matrix(yhat), RSS, trH)
}

ptr <- function(t, y, k) {
  n <- length(t); X <- matrix(0, n, k+3)
  X[,1]<-1; X[,2]<-t; X[,3]<-t^2
  if (k>0) for (j in 1:k) X[,j+3] <- sqrt(2)*cos(pi*j*t)
  qr_X <- qr(X); beta <- qr.coef(qr_X, y); beta[is.na(beta)] <- 0
  fit <- X %*% beta; RSS <- as.numeric(t(y-fit)%*%(y-fit))
  list(fit, RSS, qr_X$rank)
}

ptr_wrap <- function(t, y, k) {
  res <- tryCatch(ptr(t, y, k), error = function(e) {
    yhat <- rep(mean(y), length(y)); RSS <- sum((y-yhat)^2)
    list(as.matrix(yhat), RSS, 1)
  })
  list(as.matrix(res[[1]]), res[[2]], res[[3]])
}

ncs_gcv_edf <- function(t, y, lam) {
  n <- length(t); K <- outer(t, t, function(a,b) abs(a-b)^3)
  Tmat <- rbind(rep(1,n), t)
  A <- rbind(cbind(K+n*lam*diag(n), t(Tmat)), cbind(Tmat, matrix(0,2,2)))
  Ainv <- solve(A); S <- cbind(K, t(Tmat)) %*% Ainv[,1:n]; edf <- sum(diag(S))
  sol <- Ainv %*% c(y, rep(0,2))
  fitted <- as.vector(t(Tmat)%*%sol[(n+1):(n+2)] + K%*%sol[1:n])
  rss <- sum((y-fitted)^2)/n; gcv <- rss/(1-edf/n)^2
  list(fitted=fitted, gcv=gcv, edf=edf, c=sol[1:n], d=sol[(n+1):(n+2)])
}

fit_ncs <- function(t, y, log_min=-6, log_max=6, M_coarse=29, M_fine=10) {
  coarse_log <- seq(log_min, log_max, length.out=M_coarse)
  Delta <- coarse_log[2]-coarse_log[1]; lam_c <- 10^coarse_log
  gcv_c <- numeric(M_coarse); edf_c <- numeric(M_coarse)
  for (i in seq_along(lam_c)) {
    res <- ncs_gcv_edf(t, y, lam_c[i]); gcv_c[i] <- res$gcv; edf_c[i] <- res$edf
  }
  best_c <- which.min(gcv_c); best_log_c <- coarse_log[best_c]
  fine_log <- pmax(pmin(seq(best_log_c-Delta, best_log_c+Delta, length.out=M_fine), log_max), log_min)
  lam_f <- 10^fine_log; gcv_f <- numeric(M_fine); edf_f <- numeric(M_fine)
  for (i in seq_along(lam_f)) {
    res <- ncs_gcv_edf(t, y, lam_f[i]); gcv_f[i] <- res$gcv; edf_f[i] <- res$edf
  }
  best_f <- which.min(gcv_f); best_lam <- lam_f[best_f]; final <- ncs_gcv_edf(t, y, best_lam)
  list(fitted=final$fitted, c=final$c, d=final$d, lambda=best_lam, gcv=final$gcv,
       edf=final$edf, coarse_log=coarse_log, gcv_c=gcv_c, edf_c=edf_c,
       fine_log=fine_log, gcv_f=gcv_f, edf_f=edf_f, best_c=best_c, best_f=best_f,
       boundary=best_f==1|best_f==M_fine, Delta=Delta)
}

gaussian_variogram <- function(d, nugget, sill, range) {
  nugget + (sill-nugget) * (1-exp(-(d/range)^2))
}

fit_gaussian_constrained <- function(subj_id, data, zone_cutoff=6.5) {
  df_sub <- data %>% filter(Subject==subj_id) %>% arrange(Distance)
  d <- df_sub$Distance; y <- df_sub$MeanGamma; n <- length(d)
  fit <- tryCatch({
    nlsLM(y ~ gaussian_variogram(d, nugget, sill, range),
          start=list(nugget=max(0,min(y)*0.5), sill=max(y)*0.95, range=4),
          lower=c(nugget=0, sill=mean(y)*0.5, range=0.5),
          upper=c(nugget=min(y)*1.5, sill=max(y)*2, range=30),
          control=nls.lm.control(maxiter=500))
  }, error=function(e) { cat("nlsLM failed for subject", subj_id, ":", e$message, "\n"); NULL })
  if (is.null(fit)) return(NULL)
  params <- coef(fit); y_fitted <- fitted(fit); residuals_raw <- y - y_fitted
  RSS <- sum(residuals_raw^2); df_r <- n-3; RMSE <- sqrt(RSS/df_r)
  SS_tot <- sum((y-mean(y))^2); R2 <- 1 - RSS/SS_tot
  sigma2 <- RSS/n; loglik <- -n/2*log(2*pi*sigma2) - RSS/(2*sigma2)
  AIC_v <- -2*loglik + 2*3; BIC_v <- -2*loglik + log(n)*3
  zone <- ifelse(d<=zone_cutoff, "Rising (d <= 6.5mm)", "Post-peak (d > 6.5mm)")
  MSR_rising <- mean(residuals_raw[zone=="Rising (d <= 6.5mm)"])
  MSR_postpeak <- mean(residuals_raw[zone=="Post-peak (d > 6.5mm)"])
  d_fine <- seq(min(d), max(d), length.out=300)
  y_fine <- gaussian_variogram(d_fine, params["nugget"], params["sill"], params["range"])
  list(Subject=subj_id, params=params, RSS=RSS, RMSE=RMSE, R2=R2, AIC=AIC_v, BIC=BIC_v,
       MSR_rising=MSR_rising, MSR_postpeak=MSR_postpeak,
       Max_fail=max(abs(residuals_raw[zone=="Post-peak (d > 6.5mm)"])),
       d_obs=d, y_obs=y, y_fitted_obs=y_fitted, residuals=residuals_raw, zone=zone,
       d_fine=d_fine, y_gauss_fine=y_fine, StdGamma=df_sub$StdGamma)
}

# --- Load data ---
df_binned <- read_excel("all_subjects_variogram_binned.xlsx")
subjects <- sort(unique(df_binned$Subject)); n_subj <- length(subjects); n_obs <- 17
cat(sprintf("Binned variogram: %d subjects x %d bins\n\n", n_subj, n_obs))

# --- Step 1: Hole effect ---
COL_HOLE <- "#378ADD"; COL_MONO <- "gray65"; COL_PEAK <- "#D85A30"

hole_normalized <- df_binned %>%
  group_by(Subject) %>%
  arrange(Distance, .by_group=TRUE) %>%
  mutate(Gamma_norm=(MeanGamma-min(MeanGamma))/(max(MeanGamma)-min(MeanGamma))) %>%
  summarise(
    Peak_distance=Distance[which.max(MeanGamma)],
    Peak_gamma_raw=max(MeanGamma), Min_gamma_raw=min(MeanGamma),
    First_gamma_raw=first(MeanGamma), Last_gamma_raw=last(MeanGamma),
    Sill_raw=mean(MeanGamma[floor(0.8*n()):n()]),
    Peak_gamma_norm=Gamma_norm[which.max(MeanGamma)],
    Post_peak_min_norm={ pk <- which.max(Gamma_norm)
    if(pk<n()) min(Gamma_norm[(pk+1):n()]) else NA_real_ },
    Last_gamma_norm=last(Gamma_norm),
    Hole_depth_norm_abs=Peak_gamma_norm - Post_peak_min_norm,
    Hole_depth_norm_pct=round((Peak_gamma_norm-Post_peak_min_norm)*100,1),
    Post_peak_min_raw={ pk <- which.max(MeanGamma)
    if(pk<n()) min(MeanGamma[(pk+1):n()]) else NA_real_ },
    Hole_depth_raw_pct=round((max(MeanGamma)-Post_peak_min_raw)/max(MeanGamma)*100,1),
    N_decreases=sum(diff(MeanGamma)<0),
    Pct_decreases=round(sum(diff(MeanGamma)<0)/(n()-1)*100,1),
    .groups="drop"
  ) %>%
  mutate(Hole_5pct_raw=Hole_depth_raw_pct>5,
         Hole_5pct_norm=Hole_depth_norm_pct>5,
         Hole_10pct_norm=Hole_depth_norm_pct>10)

cat("=== Hole Effect Summary ===\n")
print(hole_normalized %>%
        select(Subject, Peak_distance, Hole_depth_raw_pct, Hole_depth_norm_pct,
               Hole_5pct_raw, Hole_5pct_norm, Hole_10pct_norm), n=Inf)

df_norm <- df_binned %>%
  group_by(Subject) %>% arrange(Distance, .by_group=TRUE) %>%
  mutate(Gamma_norm=(MeanGamma-min(MeanGamma))/(max(MeanGamma)-min(MeanGamma))) %>%
  ungroup() %>%
  left_join(hole_normalized %>%
              select(Subject, Hole_5pct_norm, Peak_distance, Hole_depth_norm_pct),
            by="Subject") %>%
  mutate(HoleGroup=ifelse(Hole_5pct_norm, "Hole effect present", "Monotone / weak"))

pC <- ggplot(df_norm, aes(x=Distance, y=Gamma_norm)) +
  geom_line(aes(color=HoleGroup), linewidth=0.9) +
  geom_point(data=hole_normalized, aes(x=Peak_distance, y=1.0),
             color=COL_PEAK, size=2.5, shape=18, inherit.aes=FALSE) +
  geom_text(data=hole_normalized, aes(x=16, y=0.15, label=paste0(Hole_depth_norm_pct, "%")),
            size=2.8, color="gray30", inherit.aes=FALSE) +
  scale_color_manual(values=c("Hole effect present"=COL_HOLE, "Monotone / weak"="gray50")) +
  facet_wrap(~Subject, nrow=4, ncol=6, labeller=label_both) +
  labs(title="Normalized Empirical Variograms — All 21 Subjects",
       subtitle="Orange diamond = peak location; gray text = post-peak drop %",
       x="Distance (mm)", y="Normalized semivariance", color=NULL) +
  theme_minimal(base_size=9) + theme(legend.position="bottom", strip.text=element_text(size=8))
print(pC)

write.csv(hole_normalized, "act1_step1_hole_normalized.csv", row.names=FALSE)
cat("Step 1 complete.\n\n")

# --- Step 2: Gaussian fit ---
ZONE_CUTOFF <- 6.5; COL_NP <- "#378ADD"; COL_FAIL <- "#FAEEDA"

gauss_results <- lapply(subjects, fit_gaussian_constrained, data=df_binned, zone_cutoff=ZONE_CUTOFF)
names(gauss_results) <- as.character(subjects)
gauss_results <- Filter(Negate(is.null), gauss_results)
cat("Successfully fitted:", length(gauss_results), "/", n_subj, "subjects\n\n")

gauss_table <- data.frame(
  Subject=sapply(gauss_results, function(r) r$Subject),
  Nugget=round(sapply(gauss_results, function(r) r$params["nugget"]),3),
  Sill=round(sapply(gauss_results, function(r) r$params["sill"]),3),
  Range_mm=round(sapply(gauss_results, function(r) r$params["range"]),2),
  RSS=round(sapply(gauss_results, function(r) r$RSS),4),
  RMSE=round(sapply(gauss_results, function(r) r$RMSE),4),
  R2=round(sapply(gauss_results, function(r) r$R2),4),
  AIC=round(sapply(gauss_results, function(r) r$AIC),2),
  MSR_rising=round(sapply(gauss_results, function(r) r$MSR_rising),4),
  MSR_postpeak=round(sapply(gauss_results, function(r) r$MSR_postpeak),4),
  Max_fail_mm=round(sapply(gauss_results, function(r) r$Max_fail),4),
  row.names=NULL, stringsAsFactors=FALSE
) %>% arrange(Subject)

cat("=== Gaussian Model — Constrained Fit ===\n"); print(gauss_table)
cat(sprintf("\nParameter summary (mean +/- SD):\n"))
cat(sprintf("  Nugget: %.3f +/- %.3f\n", mean(gauss_table$Nugget), sd(gauss_table$Nugget)))
cat(sprintf("  Sill:   %.3f +/- %.3f\n", mean(gauss_table$Sill), sd(gauss_table$Sill)))
cat(sprintf("  Range:  %.2f +/- %.2f mm\n", mean(gauss_table$Range_mm), sd(gauss_table$Range_mm)))
cat(sprintf("  RMSE:   %.4f +/- %.4f\n", mean(gauss_table$RMSE), sd(gauss_table$RMSE)))
cat(sprintf("  R2:     %.4f +/- %.4f\n", mean(gauss_table$R2), sd(gauss_table$R2)))

resid_all <- do.call(rbind, lapply(gauss_results, function(r) {
  data.frame(Subject=r$Subject, Distance=r$d_obs, Observed=r$y_obs,
             Fitted=r$y_fitted_obs, Residual=r$residuals, Zone=r$zone, StdGamma=r$StdGamma)
})) %>%
  group_by(Subject) %>%
  mutate(Residual_norm=Residual/(max(Observed)-min(Observed))) %>% ungroup()

cat("\n=== Signed residual by zone ===\n")
resid_summary <- resid_all %>%
  group_by(Zone) %>%
  summarise(Mean_resid=round(mean(Residual),4), SD_resid=round(sd(Residual),4),
            Mean_resid_norm=round(mean(Residual_norm),4),
            Pct_negative=round(mean(Residual<0)*100,1), N=n(), .groups="drop")
print(resid_summary)

postpeak_resids <- resid_all %>% filter(Zone=="Post-peak (d > 6.5mm)") %>% pull(Residual)
sign_test <- binom.test(sum(postpeak_resids<0), length(postpeak_resids), p=0.5, alternative="greater")
cat("\nSign test (post-peak residuals negative):\n")
cat("  Negative:", sum(postpeak_resids<0), "/", length(postpeak_resids),
    paste0("(", round(mean(postpeak_resids<0)*100,1), "%)\n"))
cat("  p-value:", format(sign_test$p.value, digits=4), "\n")

write.csv(gauss_table, "act1_step2_gaussian_fit_corrected.csv", row.names=FALSE)
write.csv(resid_all, "act1_step2_gaussian_residuals_corrected.csv", row.names=FALSE)
cat("Step 2 complete.\n\n")

# --- Step 3: NCS, PTR, LP vs Gaussian ---
COL_NCS <- "#378ADD"; COL_PTR <- "#BA7517"; COL_LP <- "#1D9E75"
COL_GAUSS <- "#D85A30"; COL_SEL <- "#D85A30"
k_mesh <- 1:10; bw_mesh <- seq(0.05, 0.60, by=0.025)

cat("Fitting NCS, PTR, LP to binned data...\n\n")
cat(sprintf("%-4s %6s %6s %8s %8s %8s %8s %8s\n",
            "Subj","PTR_k","LP_bw","PTR_RSS","LP_RSS","NCS_RSS","NCS_EDF","G_RSS"))
cat(strrep("-",66),"\n")

fits <- list()

for (s in subjects) {
  df_sub <- df_binned %>% filter(Subject==s) %>% arrange(Distance)
  t_mm <- df_sub$Distance; y <- df_sub$MeanGamma; n <- length(t_mm)
  t_range <- max(t_mm)-min(t_mm); t_min <- min(t_mm)
  t_norm <- (t_mm-t_min)/t_range; t_fine <- seq(0,1,length.out=300)
  t_mm_f <- t_fine*t_range+t_min
  
  ncs <- fit_ncs(t_norm, y); ncs_rss <- sum((y-ncs$fitted)^2)
  K_fine <- outer(t_fine, t_norm, function(a,b) abs(a-b)^3)
  ncs_f <- as.vector(cbind(1,t_fine)%*%ncs$d + K_fine%*%ncs$c)
  
  ptr_crit <- search_fn(t_norm, y, k_mesh, ptr_wrap); best_k <- k_mesh[which.min(ptr_crit[,2])]
  ptr_res <- ptr(t_norm, y, best_k); ptr_fit <- as.vector(ptr_res[[1]]); ptr_rss <- ptr_res[[2]]
  X_obs <- matrix(0,n,best_k+3); X_obs[,1]<-1; X_obs[,2]<-t_norm; X_obs[,3]<-t_norm^2
  if(best_k>0) for(j in 1:best_k) X_obs[,j+3]<-sqrt(2)*cos(pi*j*t_norm)
  X_fine <- matrix(0,300,best_k+3); X_fine[,1]<-1; X_fine[,2]<-t_fine; X_fine[,3]<-t_fine^2
  if(best_k>0) for(j in 1:best_k) X_fine[,j+3]<-sqrt(2)*cos(pi*j*t_fine)
  qr_obs <- qr(X_obs); beta <- qr.coef(qr_obs, y); beta[is.na(beta)] <- 0
  ptr_f <- as.vector(X_fine%*%beta)
  
  lp_crit <- search_fn(t_norm, y, bw_mesh, locpoly3_wrap); best_bw <- bw_mesh[which.min(lp_crit[,2])]
  lp_fit <- locpoly3(t_norm, y, best_bw, t_norm)[[1]][,1]; lp_rss <- sum((y-lp_fit)^2)
  lp_f <- locpoly3(t_norm, y, best_bw, t_fine)[[1]][,1]
  
  sill0 <- as.numeric(quantile(y,0.90)); rng0 <- t_mm[which(y>=0.8*sill0)[1]]
  rng0 <- ifelse(is.na(rng0)|rng0<1, mean(t_mm)/2, rng0)
  nug0 <- max(0, as.numeric(quantile(y,0.05)))
  gfit <- tryCatch(
    nlsLM(y ~ gaussian_variogram(t_mm,nugget,sill,range),
          start=list(nugget=nug0,sill=sill0,range=rng0),
          lower=c(nugget=0,sill=sill0*0.5,range=0.5),
          upper=c(nugget=sill0*0.3,sill=Inf,range=25),
          control=nls.lm.control(maxiter=500)),
    error=function(e) tryCatch(
      nlsLM(y ~ gaussian_variogram(t_mm,nugget,sill,range),
            start=list(nugget=0,sill=max(y)*0.9,range=mean(t_mm)),
            lower=c(nugget=0,sill=max(y)*0.3,range=0.5),
            upper=c(nugget=Inf,sill=Inf,range=25),
            control=nls.lm.control(maxiter=500)),
      error=function(e2) NULL))
  
  if (!is.null(gfit)) {
    gpar <- coef(gfit); gauss_fit <- predict(gfit); gauss_rss <- sum((y-gauss_fit)^2)
    gauss_f <- gaussian_variogram(t_mm_f,gpar["nugget"],gpar["sill"],gpar["range"])
    pp_idx <- t_mm>gpar["range"]
    if(sum(pp_idx)>2) {
      y_pp<-y[pp_idx]; f_pp<-gauss_fit[pp_idx]
      ss_r<-sum((y_pp-f_pp)^2); ss_t<-sum((y_pp-mean(y_pp))^2)
      g_r2pp <- if(ss_t>0) 1-ss_r/ss_t else NA
    } else g_r2pp <- NA
    g_nugget<-gpar["nugget"]; g_sill<-gpar["sill"]; g_range<-gpar["range"]
  } else {
    gauss_fit<-rep(NA,n); gauss_f<-rep(NA,300); gauss_rss<-NA; g_r2pp<-NA
    g_nugget<-NA; g_sill<-NA; g_range<-NA
  }
  
  fits[[as.character(s)]] <- list(
    t_mm=t_mm, t_norm=t_norm, y=y, n=n, t_mm_f=t_mm_f,
    ncs_fit=ncs$fitted, ncs_fine=ncs_f, ncs_lambda=ncs$lambda, ncs_gcv=ncs$gcv,
    ncs_edf=ncs$edf, ncs_rss=ncs_rss, ncs_boundary=ncs$boundary,
    ncs_coarse_log=ncs$coarse_log, ncs_gcv_c=ncs$gcv_c, ncs_edf_c=ncs$edf_c,
    ncs_fine_log=ncs$fine_log, ncs_gcv_f=ncs$gcv_f, ncs_edf_f=ncs$edf_f,
    ptr_fit=ptr_fit, ptr_fine=ptr_f, ptr_k=best_k, ptr_edf=best_k+3, ptr_rss=ptr_rss,
    ptr_crit=ptr_crit, lp_fit=lp_fit, lp_fine=lp_f,
    lp_bw=best_bw, lp_bw_mm=best_bw*t_range, lp_rss=lp_rss,
    gauss_fit=gauss_fit, gauss_fine=gauss_f, gauss_nugget=g_nugget,
    gauss_sill=g_sill, gauss_range=g_range, gauss_rss=gauss_rss, gauss_r2pp=g_r2pp)
  
  cat(sprintf("%-4d %6d %6.3f %8.5f %8.5f %8.5f %8.3f %8.5f\n",
              s, best_k, best_bw*t_range, ptr_rss, lp_rss, ncs_rss, ncs$edf,
              ifelse(is.na(gauss_rss), NA, gauss_rss)))
}
cat("\nDone.\n\n")

perf <- do.call(rbind, lapply(as.character(subjects), function(s) {
  r <- fits[[s]]
  data.frame(Subject=as.integer(s), PTR_k=r$ptr_k, PTR_EDF=r$ptr_edf,
             PTR_RSS=round(r$ptr_rss,5), LP_bw_mm=round(r$lp_bw_mm,3),
             LP_RSS=round(r$lp_rss,5), NCS_EDF=round(r$ncs_edf,3),
             NCS_lambda=signif(r$ncs_lambda,3), NCS_RSS=round(r$ncs_rss,5),
             G_RSS=round(r$gauss_rss,5), G_range=round(r$gauss_range,3),
             G_R2pp=round(r$gauss_r2pp,4), stringsAsFactors=FALSE)
}))

cat("=== RSS COMPARISON ===\n")
cat(sprintf("%-4s %5s %8s %8s %8s %8s %8s\n","Subj","PTR_k","PTR_RSS","LP_RSS","NCS_RSS","G_RSS","G_R2pp"))
cat(strrep("-",57),"\n")
for (i in seq_len(nrow(perf))) {
  r <- perf[i,]
  cat(sprintf("%-4d %5d %8.5f %8.5f %8.5f %8.5f %8.4f\n",
              r$Subject, r$PTR_k, r$PTR_RSS, r$LP_RSS, r$NCS_RSS,
              ifelse(is.na(r$G_RSS),NA,r$G_RSS), ifelse(is.na(r$G_R2pp),0,r$G_R2pp)))
}
cat(sprintf("\nNCS beats PTR: %d/21 | NCS beats LP: %d/21 | LP beats PTR: %d/21\n",
            sum(perf$NCS_RSS<perf$PTR_RSS), sum(perf$NCS_RSS<perf$LP_RSS), sum(perf$LP_RSS<perf$PTR_RSS)))
cat(sprintf("NCS beats Gaussian:  %d/21 | mean G/NCS ratio=%.3f\n",
            sum(perf$NCS_RSS<perf$G_RSS,na.rm=TRUE), mean(perf$G_RSS/perf$NCS_RSS,na.rm=TRUE)))
cat(sprintf("LP  beats Gaussian:  %d/21 | mean G/LP ratio=%.3f\n",
            sum(perf$LP_RSS<perf$G_RSS,na.rm=TRUE), mean(perf$G_RSS/perf$LP_RSS,na.rm=TRUE)))
cat(sprintf("PTR beats Gaussian:  %d/21 | mean G/PTR ratio=%.3f\n",
            sum(perf$PTR_RSS<perf$G_RSS,na.rm=TRUE), mean(perf$G_RSS/perf$PTR_RSS,na.rm=TRUE)))
cat(sprintf("Mean RSS: NCS=%.5f LP=%.5f PTR=%.5f G=%.5f\n",
            mean(perf$NCS_RSS), mean(perf$LP_RSS), mean(perf$PTR_RSS), mean(perf$G_RSS,na.rm=TRUE)))
cat(sprintf("NCS EDF: mean=%.3f SD=%.3f range=%.3f-%.3f\n",
            mean(perf$NCS_EDF), sd(perf$NCS_EDF), min(perf$NCS_EDF), max(perf$NCS_EDF)))
cat(sprintf("PTR k: mean=%.1f SD=%.1f range=%d-%d\n",
            mean(perf$PTR_k), sd(perf$PTR_k), min(perf$PTR_k), max(perf$PTR_k)))
cat(sprintf("LP bw: mean=%.2fmm SD=%.2fmm\n", mean(perf$LP_bw_mm), sd(perf$LP_bw_mm)))
cat(sprintf("Gaussian post-peak R2: mean=%.4f | negative: %d/21\n",
            mean(perf$G_R2pp,na.rm=TRUE), sum(perf$G_R2pp<0,na.rm=TRUE)))
cat(sprintf("NCS boundary selections: %d/21\n", sum(sapply(fits, function(r) r$ncs_boundary))))

fits_long <- do.call(rbind, lapply(as.character(subjects), function(s) {
  r <- fits[[s]]
  rbind(data.frame(Subject=as.integer(s), Distance=r$t_mm_f, Fitted=r$ncs_fine,   Model="NCS"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm_f, Fitted=r$ptr_fine,   Model="PTR"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm_f, Fitted=r$lp_fine,    Model="LP"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm_f, Fitted=r$gauss_fine, Model="Gaussian"))
}))

obs_long <- do.call(rbind, lapply(as.character(subjects), function(s) {
  r <- fits[[s]]; data.frame(Subject=as.integer(s), Distance=r$t_mm, Observed=r$y)
}))

resid_long <- do.call(rbind, lapply(as.character(subjects), function(s) {
  r <- fits[[s]]
  rbind(data.frame(Subject=as.integer(s), Distance=r$t_mm, Residual=r$y-r$ncs_fit,   Model="NCS"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm, Residual=r$y-r$ptr_fit,   Model="PTR"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm, Residual=r$y-r$lp_fit,    Model="LP"),
        data.frame(Subject=as.integer(s), Distance=r$t_mm, Residual=r$y-r$gauss_fit, Model="Gaussian"))
})) %>% filter(!is.na(Residual))

rss_long <- perf %>%
  dplyr::select(Subject, NCS=NCS_RSS, PTR=PTR_RSS, LP=LP_RSS, Gaussian=G_RSS) %>%
  pivot_longer(-Subject, names_to="Model", values_to="RSS") %>% filter(!is.na(RSS))

p2 <- ggplot(rss_long, aes(x=reorder(factor(Subject),RSS), y=RSS, fill=Model)) +
  geom_col(position="dodge", width=0.7, alpha=0.85) +
  scale_fill_manual(values=c(NCS=COL_NCS, PTR=COL_PTR, LP=COL_LP, Gaussian=COL_GAUSS)) +
  labs(title="RSS Comparison — All Four Methods (binned data)",
       subtitle="Lower = better fit | Gaussian is the parametric benchmark",
       x="Subject", y="RSS", fill="Method") +
  theme_minimal(base_size=11) + theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
print(p2)

edf_df <- data.frame(Subject=subjects,
                     EDF=sapply(fits, function(r) r$ncs_edf),
                     Log10Lam=log10(sapply(fits, function(r) r$ncs_lambda)))

p3 <- ggplot(edf_df, aes(x=reorder(factor(Subject),EDF), y=EDF, fill=EDF>n_obs/2)) +
  geom_col(width=0.7, alpha=0.85) +
  geom_hline(yintercept=c(2,n_obs), linetype="dashed", color=COL_SEL, linewidth=0.9) +
  geom_hline(yintercept=c(n_obs/4,n_obs/2), linetype="dotted", color="gray50", linewidth=0.8) +
  scale_fill_manual(values=c("FALSE"=COL_NCS,"TRUE"=COL_SEL),
                    labels=c("FALSE"="Within n/2","TRUE"="Above n/2")) +
  annotate("text", x=0.6, y=2.5, label="Min df=2 (linear)", hjust=0, color=COL_SEL, size=3) +
  annotate("text", x=0.6, y=n_obs-0.6, label=sprintf("Max df=%d (interpolating)",n_obs),
           hjust=0, color=COL_SEL, size=3) +
  annotate("text", x=n_subj-1, y=n_obs/4+0.4, label="n/4", hjust=1, color="gray50", size=3) +
  annotate("text", x=n_subj-1, y=n_obs/2+0.4, label="n/2", hjust=1, color="gray50", size=3) +
  labs(title="NCS Effective Degrees of Freedom",
       subtitle=sprintf("tr(S) from full hat matrix | n=17 | mean EDF=%.2f", mean(edf_df$EDF)),
       x="Subject (ordered by EDF)", y="EDF = tr(S)", fill=NULL) +
  theme_minimal(base_size=11) +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
print(p3)

coarse_df <- do.call(rbind, lapply(as.character(subjects), function(s) {
  r <- fits[[s]]
  data.frame(Subject=as.integer(s), Log10Lam=r$ncs_coarse_log, GCV_norm=r$ncs_gcv_c/min(r$ncs_gcv_c))
}))

p4 <- ggplot(coarse_df, aes(x=Log10Lam, y=GCV_norm)) +
  geom_line(color=COL_NCS, linewidth=0.9) + geom_point(color=COL_NCS, size=1.2) +
  geom_point(data=do.call(rbind, lapply(as.character(subjects), function(s) {
    r <- fits[[s]]
    data.frame(Subject=as.integer(s), Log10Lam=r$ncs_fine_log, GCV_norm=r$ncs_gcv_f/min(r$ncs_gcv_c))
  })), aes(x=Log10Lam, y=GCV_norm), color=COL_LP, size=1.8, shape=17) +
  geom_vline(data=edf_df, aes(xintercept=Log10Lam), color=COL_SEL, linetype="dashed", linewidth=0.8) +
  geom_vline(xintercept=c(-6,6), linetype="dotted", color="gray50", linewidth=0.5) +
  facet_wrap(~Subject, nrow=4, ncol=6,
             labeller=labeller(Subject=function(x) {
               vapply(x, function(xi) {
                 r <- fits[[xi]]; sprintf("S%s\nλ=%.1e df=%.1f", xi, r$ncs_lambda, r$ncs_edf)
               }, character(1))
             })) +
  coord_cartesian(ylim=c(0.9,6)) +
  labs(title="NCS GCV Curves per Subject — Coarse + Fine",
       subtitle="Blue = coarse grid | Green = fine points | Red dashed = selected lambda",
       x=expression(log[10](lambda)), y="Normalised GCV (min=1)") +
  theme_minimal(base_size=9) + theme(strip.text=element_text(size=7), panel.grid.minor=element_blank())
print(p4)

plot_resid <- function(method, col) {
  df_m <- resid_long %>% filter(Model==method)
  avg <- df_m %>% group_by(Distance) %>%
    summarise(Mean=mean(Residual,na.rm=TRUE), SD=sd(Residual,na.rm=TRUE), .groups="drop")
  ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.9) +
    geom_line(data=df_m, aes(x=Distance, y=Residual, group=factor(Subject)),
              color="gray75", linewidth=0.5, alpha=0.7) +
    geom_ribbon(data=avg, aes(x=Distance, ymin=Mean-SD, ymax=Mean+SD), fill=col, alpha=0.22) +
    geom_line(data=avg, aes(x=Distance, y=Mean), color=col, linewidth=1.5) +
    geom_point(data=avg, aes(x=Distance, y=Mean), color=col, size=3, shape=16) +
    labs(title=sprintf("%s Residuals — All 21 Subjects (binned)", method),
         subtitle="Gray lines = individual subjects | Mean +/- 1SD | Flat = good fit",
         x="Distance (mm)", y="Residual (observed - fitted)") +
    theme_minimal(base_size=12) + theme(plot.title=element_text(face="bold"))
}
print(plot_resid("NCS", COL_NCS)); print(plot_resid("PTR", COL_PTR))
print(plot_resid("LP", COL_LP)); print(plot_resid("Gaussian", COL_GAUSS))

avg_all <- resid_long %>% group_by(Model, Distance) %>%
  summarise(Mean=mean(Residual,na.rm=TRUE), SD=sd(Residual,na.rm=TRUE), .groups="drop")

p7 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.9) +
  geom_ribbon(data=avg_all, aes(x=Distance, ymin=Mean-SD, ymax=Mean+SD, fill=Model), alpha=0.12) +
  geom_line(data=avg_all, aes(x=Distance, y=Mean, color=Model), linewidth=1.4) +
  geom_point(data=avg_all, aes(x=Distance, y=Mean, color=Model), size=2.8) +
  scale_color_manual(values=c(NCS=COL_NCS, PTR=COL_PTR, LP=COL_LP, Gaussian=COL_GAUSS)) +
  scale_fill_manual(values=c(NCS=COL_NCS, PTR=COL_PTR, LP=COL_LP, Gaussian=COL_GAUSS)) +
  labs(title="Grand Mean Residuals", subtitle="Mean +/- 1SD across 21 subjects | Binned variogram",
       x="Distance (mm)", y="Mean Residual", color="Method", fill="Method") +
  theme_minimal(base_size=12) + theme(legend.position="bottom", plot.title=element_text(face="bold"))
print(p7)

make_fit_plot <- function(model_name, col, subtitle_text) {
  ggplot() +
    geom_point(data=obs_long, aes(x=Distance, y=Observed), color="gray40", size=1.5, alpha=0.75) +
    geom_line(data=fits_long %>% filter(Model==model_name, !is.na(Fitted)),
              aes(x=Distance, y=Fitted), color=col, linetype="solid", linewidth=1.0) +
    facet_wrap(~Subject, nrow=4, ncol=6, scales="free_y") +
    labs(title=paste(model_name, "Fits (Binned Variogram)"), subtitle=subtitle_text,
         x="Distance (mm)", y=expression(hat(gamma)(d))) +
    theme_minimal(base_size=8) +
    theme(strip.text=element_text(size=6.5, face="bold"),
          panel.grid.minor=element_blank(), plot.title=element_text(size=10, face="bold"))
}
print(make_fit_plot("NCS",      COL_NCS,   "Natural Cubic Spline (GCV-selected lambda)"))
print(make_fit_plot("PTR",      COL_PTR,   "Polynomial + Trigonometric Regression (GCV-selected k)"))
print(make_fit_plot("LP",       COL_LP,    "Gaussian kernel, GCV bandwidth"))
print(make_fit_plot("Gaussian", COL_GAUSS, "Parametric benchmark"))

write.csv(perf,       "act1_step3_four_methods.csv", row.names=FALSE)
write.csv(edf_df,     "act1_step3_ncs_edf.csv",      row.names=FALSE)
write.csv(resid_long, "act1_step3_residuals.csv",     row.names=FALSE)
cat("\nact1_step1_hole_normalized.csv\nact1_step2_gaussian_fit_corrected.csv\n")
cat("act1_step2_gaussian_residuals_corrected.csv\nact1_step3_four_methods.csv\n")
cat("act1_step3_ncs_edf.csv\nact1_step3_residuals.csv\n")