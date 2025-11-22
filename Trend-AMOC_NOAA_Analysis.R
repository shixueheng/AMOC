# ===============================================================
# AMOC — NOAA only (Jmax & Fmax)
# Compute & plotting are decoupled.
# Output: AMOC_NOAA_trim0.05.pdf + printed summary table
# ===============================================================

library(latex2exp)

# -------------------------
# Load data (NOAA series)
# -------------------------
load("temperature_anomalies_Sep2024.RData")

dat <- data.frame(
  year = as.numeric(Tanom_annual_df$year),
  NOAA = as.numeric(Tanom_annual_df$NOAA),
  check.names = FALSE
)

# -------------------------
# choose k to compute the statistics at k
# -------------------------
k_grid <- function(n, trim = 0.05) {
  kmin <- max(2L, ceiling(n * trim))
  kmax <- min(n - 2L, floor(n * (1 - trim)))
  if (kmin > kmax) integer(0) else kmin:kmax
}

# -------------------------
# Test statistics (robust)
# -------------------------
Jmax_stat <- function(year, y, trim = 0.05) {
  x <- as.numeric(year); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  ks <- k_grid(length(y), trim)
  if (!length(ks)) return(list(stat = NA_real_, k = max(1L, floor(length(y)/2))))
  
  Jk <- rep(NA_real_, length(ks))
  for (i in seq_along(ks)) {
    tau <- x[ks[i]]; h <- pmax(0, x - tau)
    fit <- lm(y ~ x + h)
    est <- suppressWarnings(coef(fit)[["h"]])
    se  <- suppressWarnings(summary(fit)$coef["h", "Std. Error"])
    if (is.finite(est) && is.finite(se) && se > 0) Jk[i] <- abs(est/se)
  }
  Jk2 <- Jk; Jk2[!is.finite(Jk2)] <- -Inf
  idx <- which.max(Jk2)
  if (!is.finite(Jk2[idx])) idx <- max(1L, floor(length(ks)/2))
  list(stat = Jk[idx], k = ks[idx])
}

Fmax_stat <- function(year, y, trim = 0.05) {
  x <- as.numeric(year); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  n  <- length(y)
  ks <- k_grid(n, trim)
  if (!length(ks)) return(list(stat = NA_real_, k = max(1L, floor(n/2))))
  
  Fk <- rep(NA_real_, length(ks))
  for (i in seq_along(ks)) {
    if (ks[i] < 2 || (n - ks[i]) < 2) next
    grp <- as.numeric(seq_len(n) > ks[i])
    lm0 <- lm(y ~ x)
    lm1 <- lm(y ~ x + grp + x:grp)
    Fk[i] <- anova(lm0, lm1)$`F`[2]
  }
  Fk2 <- Fk; Fk2[!is.finite(Fk2)] <- -Inf
  idx <- which.max(Fk2)
  if (!is.finite(Fk2[idx])) idx <- max(1L, floor(length(ks)/2))
  list(stat = Fk[idx], k = ks[idx])
}

# -------------------------
# 95% criticals (from your table)
# -------------------------
q95 <- list(
  J = 2.658,      # trim = 0.05
  F = 6.835       # trim = 0.05
)

# -------------------------
# Fit helpers (robust, NA-safe)
# -------------------------
# Jmax: continuous joinpoint. If hinge fit fails, fallback to left/right OLS and enforce continuity at tau.
fit_joinpoint <- function(year, y, k) {
  x <- as.numeric(year); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  n <- length(y); k <- max(1L, min(k, n))
  
  tau <- x[k]; h <- pmax(0, x - tau)
  fit <- lm(y ~ x + h)
  cf  <- suppressWarnings(coef(fit))
  
  if (all(is.finite(cf[c("(Intercept)", "x", "h")]))) {
    a <- unname(cf["(Intercept)"]); b <- unname(cf["x"]); g <- unname(cf["h"])
    b0L <- a;           b1L <- b
    b0R <- a - g * tau; b1R <- b + g
    yfit <- fitted(fit)
  } else {
    idxL <- seq_len(n) <= k; idxR <- !idxL
    if (sum(idxL) < 2 || sum(idxR) < 2)
      return(list(yfit = rep(NA_real_, n), coefs = c(b0L=NA, b1L=NA, b0R=NA, b1R=NA)))
    fitL <- lm(y[idxL] ~ x[idxL]); fitR <- lm(y[idxR] ~ x[idxR])
    aL <- unname(coef(fitL)[1]); bL <- unname(coef(fitL)[2])
    bR <- unname(coef(fitR)[2])
    aR <- (aL + bL * tau) - bR * tau  # enforce continuity
    yfit <- numeric(n)
    yfit[idxL] <- aL + bL * x[idxL]
    yfit[idxR] <- aR + bR * x[idxR]
    b0L <- aL; b1L <- bL; b0R <- aR; b1R <- bR
  }
  list(yfit = yfit, coefs = c(b0L=b0L, b1L=b1L, b0R=b0R, b1R=b1R))
}

# Fmax: discontinuous; compute left/right via separate OLS
fit_full_break <- function(year, y, k) {
  x <- as.numeric(year); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  n <- length(y); k <- max(1L, min(k, n))
  idxL <- seq_len(n) <= k; idxR <- !idxL
  if (sum(idxL) < 2 || sum(idxR) < 2)
    return(c(b0L=NA, b1L=NA, b0R=NA, b1R=NA))
  fitL <- lm(y[idxL] ~ x[idxL]); fitR <- lm(y[idxR] ~ x[idxR])
  c(b0L = unname(coef(fitL)[1]), b1L = unname(coef(fitL)[2]),
    b0R = unname(coef(fitR)[1]), b1R = unname(coef(fitR)[2]))
}

# -------------------------
# COMPUTE-ONLY: returns a structured result, no plotting
# -------------------------
compute_amoc_noaa <- function(year, y, trim = 0.05, q95 = q95) {
  x <- as.numeric(year); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  
  # stats
  Js <- Jmax_stat(x, y, trim)
  Fs <- Fmax_stat(x, y, trim)
  
  # fits
  fitJ <- fit_joinpoint(x, y, Js$k)
  fitF <- fit_full_break(x, y, Fs$k)
  
  # package
  list(
    trim = trim,
    J = list(
      k_hat    = Js$k,
      tau_year = x[Js$k],
      stat     = Js$stat,
      crit95   = q95$J,
      reject   = isTRUE(Js$stat > q95$J),
      coefs    = fitJ$coefs,
      yfit     = fitJ$yfit
    ),
    F = list(
      k_hat    = Fs$k,
      tau_year = x[Fs$k],
      stat     = Fs$stat,
      crit95   = q95$F,
      reject   = isTRUE(Fs$stat > q95$F),
      coefs    = fitF
    ),
    x = x,
    y = y
  )
}

# -------------------------
# PLOTTING (separate from compute)
# -------------------------
plot_amoc_noaa <- function(x, y, res, pdf_file = "AMOC_NOAA_trim0.05.pdf") {
  stopifnot(is.list(res), all(c("J","F") %in% names(res)))
  tauJ <- res$J$tau_year; tauF <- res$F$tau_year
  
  pdf(pdf_file, width = 10, height = 4)
  oldpar <- par(mfrow = c(1,2), mar = c(4,4,3,1), mgp = c(2,0.7,0)); on.exit({par(oldpar); dev.off()}, add = TRUE)
  
  # --- Jmax (continuous) ---
  #plot(x, y, type="l", col="gray40", lwd=0.5,
  #     xlab="Year", ylab="Temperature anomaly (°C)",
  #     main=sprintf("Jmax (τ̂ = %d)", round(tauJ)))
  plot(
    x, y, type = "l", col = "gray40", lwd = 0.5,
    xlab = "Year", ylab = "Temperature anomaly (°C)",
    main = TeX(sprintf("$J_{max}\\; (\\hat{\\tau} = %d)$", round(tauJ)))
  )
  if (any(is.finite(res$J$yfit)))
    lines(x, res$J$yfit, col="blue", lty=2, lwd=1)
  abline(v = tauJ, col = "red", lty = 3, lwd = 2)
  legend("topleft",
         legend=c("Observed","Fitted","Changepoint"),
         col=c("gray40","blue","red"), lty=c(1,2,3), lwd=2, bty="n", cex=0.8)
  
  # --- Fmax (discontinuous) ---
  plot(x, y, type="l", col="gray40", lwd=0.5,
       xlab="Year", ylab="Temperature anomaly (°C)",
       main=TeX(sprintf("$F_{max}\\; (\\hat{\\tau} = %d)$", round(tauF))))
         #sprintf("Fmax (τ̂ = %d)", round(tauF)))
  k <- res$F$k_hat
  idxL <- seq_along(x) <= k; idxR <- !idxL
  lines(x[idxL], res$F$coefs["b0L"] + res$F$coefs["b1L"]*x[idxL], col="blue", lty=2, lwd=1)
  lines(x[idxR], res$F$coefs["b0R"] + res$F$coefs["b1R"]*x[idxR], col="blue", lty=2, lwd=1)
  abline(v = tauF, col = "red", lty = 3, lwd = 2)
  legend("topleft",
         legend=c("Observed","Fitted","Changepoint"),
         col=c("gray40","blue","red"), lty=c(1,2,3), lwd=2, bty="n", cex=0.8)
}

# -------------------------
# Summary table maker (from compute result)
# -------------------------
make_summary_df <- function(res) {
  stopifnot(is.list(res), all(c("J","F") %in% names(res)))
  data.frame(
    Test         = c("Jmax","Fmax"),
    k_hat        = c(res$J$k_hat, res$F$k_hat),
    tau_year     = c(res$J$tau_year, res$F$tau_year),
    Stat         = round(c(res$J$stat, res$F$stat), 6),
    Crit95       = c(res$J$crit95, res$F$crit95),
    Reject95     = c(res$J$reject, res$F$reject),
    Intercept_L  = c(res$J$coefs["b0L"], res$F$coefs["b0L"]),
    Slope_L      = c(res$J$coefs["b1L"], res$F$coefs["b1L"]),
    Intercept_R  = c(res$J$coefs["b0R"], res$F$coefs["b0R"]),
    Slope_R      = c(res$J$coefs["b1R"], res$F$coefs["b1R"]),
    check.names  = FALSE
  )
}

# =========================
# RUN (compute → plot → summarize)
# =========================
trim_used <- 0.05
res_noaa  <- compute_amoc_noaa(dat$year, dat$NOAA, trim = trim_used, q95 = q95)
plot_amoc_noaa(res_noaa$x, res_noaa$y, res_noaa, pdf_file = sprintf("AMOC_NOAA_trim%.2f.pdf", trim_used))
summary_noaa <- make_summary_df(res_noaa)
print(summary_noaa)
