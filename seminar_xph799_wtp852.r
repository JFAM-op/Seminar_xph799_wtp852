# =============================================================================
# MS-GARCH Analysis of Danish Electricity Spot Prices DK2 (East Denmark)
# Data  : Hourly SpotPriceDKK 
# Model : Markov-Switching GARCH  (MSGARCH package, Ardia et al.)
# =============================================================================
# Please do keep in mind, that you need to change the paths for output_dir and read_excel (line 240-260)
# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
if (!require(readxl,  quietly = TRUE)) install.packages("readxl")
if (!require(MSGARCH, quietly = TRUE)) install.packages("MSGARCH")
if (!require(tseries,  quietly = TRUE)) install.packages("tseries")
if (!require(FinTS,    quietly = TRUE)) install.packages("FinTS")
if (!require(moments,  quietly = TRUE)) install.packages("moments")
if (!require(MASS,     quietly = TRUE)) install.packages("MASS")
if (!require(xtable,   quietly = TRUE)) install.packages("xtable")
if (!require(parallel, quietly = TRUE)) install.packages("parallel")

# ── Global parallel configuration ────────────────────────────────────────────
# Explicit worker count for the block bootstrap and the OOS refits.
N_CORES <- 6L
cat(sprintf("Parallel workers configured: %d\n", N_CORES))
cat(sprintf("  (detected physical cores: %d, logical: %d)\n",
            parallel::detectCores(logical = FALSE),
            parallel::detectCores(logical = TRUE)))

# ── 0B. HELPER FUNCTIONS ─────────────────────────────────────────────────────
# All reusable functions used across phases 3-5: parameter extraction,

# Extract scalar named parameter safely (returns NA if not found)
get_par <- function(p, nm) {
  v <- p[nm]
  if (length(v) == 0 || is.na(v)) NA_real_ else as.numeric(v)
}

# K x K transition matrix from MSGARCH par vector.
# the last column is 1 minus rowSums of the stored values.
build_trans_mat <- function(par, K) {
  P <- matrix(0, K, K)
  for (k in seq_len(K)) {
    for (j in seq_len(K - 1)) {
      nm <- sprintf("P_%d_%d", k, j)
      P[k, j] <- get_par(par, nm)
    }
    P[k, K] <- 1 - sum(P[k, seq_len(K - 1)])
  }
  rownames(P) <- colnames(P) <- paste0("k=", seq_len(K))
  P
}

# Ergodic distribution pi such that pi'P = pi', sum(pi) = 1.
ergodic_dist <- function(P) {
  K <- nrow(P)
  A <- t(P) - diag(K)
  A <- rbind(A, rep(1, K))
  b <- c(rep(0, K), 1)
  as.numeric(qr.solve(A, b))
}

# Compute per-regime parameter table for a fitted MSGARCH model.
regime_table <- function(fit_obj) {
  par   <- fit_obj$par
  pnms  <- names(par)
  k_idx <- as.integer(gsub(".*_", "", grep("^alpha0_", pnms, value = TRUE)))
  K     <- max(k_idx)
  gtype <- strsplit(fit_obj$spec$name[1], "_")[[1]][1]
  
  rows <- lapply(seq_len(K), function(k) {
    omega  <- get_par(par, sprintf("alpha0_%d", k))
    alpha1 <- get_par(par, sprintf("alpha1_%d", k))
    alpha2 <- get_par(par, sprintf("alpha2_%d", k))
    beta   <- get_par(par, sprintf("beta_%d",   k))
    nu     <- get_par(par, sprintf("nu_%d",     k))
    
    if (gtype == "eGARCH") {
      persist <- abs(beta)
      uncvol  <- NA_real_
    } else {
      a2_eff  <- ifelse(is.na(alpha2), 0, alpha2)
      persist <- alpha1 + 0.5 * a2_eff + beta
      uncvar  <- omega / pmax(1 - persist, 1e-6)
      uncvol  <- sqrt(uncvar)
    }
    
    data.frame(
      Regime = k, omega = omega, alpha1 = alpha1,
      alpha2 = alpha2, beta = beta, nu = nu,
      persistence = persist, uncond_vol = uncvol,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Return permutation that sorts regimes from low to high volatility.
# NOTE: eGARCH has no closed-form unconditional variance (log-variance
# parametrisation), so we order eGARCH regimes by |alpha1| as a proxy for
# ARCH responsiveness. sGARCH/GJR use the true unconditional vol.
sort_regimes <- function(rtbl, gtype) {
  if (gtype == "eGARCH") {
    key <- abs(rtbl$alpha1)
  } else if (all(is.finite(rtbl$uncond_vol))) {
    key <- rtbl$uncond_vol
  } else {
    key <- rtbl$persistence
  }
  order(key)
}

# Standard CDF for each supported distribution.
regime_cdf <- function(x, sigma, dist, df = NULL) {
  switch(dist,
         "norm" = pnorm(x, sd = sigma),
         "std"  = pt(x / (sigma * sqrt((df - 2) / df)), df = df),
         stop("Unknown distribution: ", dist)
  )
}

# Compute in-sample VaR vector for a fitted MSGARCH object at given alpha.
# Works for any K and for both "norm" and "std" error distributions.
compute_var_vec <- function(fit_obj, alpha_level) {
  par_o  <- fit_obj$par
  pnms_o <- names(par_o)
  K_o    <- max(as.integer(gsub(".*_", "", grep("^alpha0_", pnms_o, value = TRUE))))
  dt_o   <- if (!is.na(get_par(par_o, "nu_1"))) "std" else "norm"
  T_o    <- length(fit_obj$data)
  
  nu_o <- sapply(seq_len(K_o), function(k) get_par(par_o, sprintf("nu_%d", k)))
  
  st_o <- State(fit_obj)
  pp_o <- st_o$PredProb[seq_len(T_o), 1, , drop = FALSE]
  pp_o <- pp_o[, 1, ]
  if (is.null(dim(pp_o))) pp_o <- matrix(pp_o, ncol = K_o)
  
  ef_o <- ExtractStateFit(fit_obj)
  vr_o <- sapply(ef_o, Volatility)
  if (is.null(dim(vr_o))) vr_o <- matrix(vr_o, ncol = K_o)
  
  n_grid <- 4000L                        
  var_o  <- numeric(T_o)
  for (t in seq_len(T_o)) {
    w_t  <- pp_o[t, ]
    s_t  <- vr_o[t, ]
    g_t  <- 20 * max(s_t, na.rm = TRUE)  
    xg_t <- seq(-g_t, g_t, length.out = n_grid)
    mcdf <- numeric(n_grid)
    for (k in seq_len(K_o)) {
      mcdf <- mcdf + w_t[k] * regime_cdf(xg_t, s_t[k], dt_o, nu_o[k])
    }
    var_o[t] <- approx(mcdf, xg_t, xout = alpha_level, rule = 2)$y
  }
  var_o
}

# Analytic Expected Shortfall for a mixture of scaled Student-t distributions.
# Requires df_k > 1 for every regime; returns NA otherwise.
mix_es_analytic <- function(VaR_t, w, sigma_k, nu_k) {
  if (!is.finite(VaR_t)) return(NA_real_)
  integral <- 0
  for (k in seq_along(w)) {
    df_k <- nu_k[k]
    if (is.na(df_k) || !is.finite(df_k) || df_k <= 1) next
    x_star  <- VaR_t / sigma_k[k] * sqrt(df_k / (df_k - 2))
    integ_k <- sigma_k[k] * sqrt((df_k - 2) / df_k) *
      (- dt(x_star, df = df_k) * (df_k + x_star^2) / (df_k - 1))
    integral <- integral + w[k] * integ_k
  }
  integral
}

# Christoffersen (1998) Conditional Coverage test.
christoffersen_cc <- function(returns, VaR, alpha) {
  hit <- as.integer(returns < VaR)
  T   <- length(hit)
  n1  <- sum(hit)
  n0  <- T - n1
  pi_hat <- n1 / T
  
  if (n1 == 0 || n1 == T) {
    warning("All or no observations are violations; test undefined.")
    return(list(LR_UC = NA, LR_IND = NA, LR_CC = NA,
                p_UC = NA, p_IND = NA, p_CC = NA,
                hit_rate = pi_hat, expected_rate = alpha,
                n_violations = n1, n_obs = T))
  }
  
  LR_UC <- -2 * (n1 * log(alpha)    + n0 * log(1 - alpha) -
                   n1 * log(pi_hat)   - n0 * log(1 - pi_hat))
  
  n00 <- sum(hit[-T] == 0 & hit[-1] == 0)
  n01 <- sum(hit[-T] == 0 & hit[-1] == 1)
  n10 <- sum(hit[-T] == 1 & hit[-1] == 0)
  n11 <- sum(hit[-T] == 1 & hit[-1] == 1)
  
  pi01 <- n01 / max(n00 + n01, 1)
  pi11 <- n11 / max(n10 + n11, 1)
  pi2  <- (n01 + n11) / max(n00 + n01 + n10 + n11, 1)
  
  safe_log <- function(x) ifelse(x <= 0 | is.nan(x), 0, log(x))
  
  LR_IND <- -2 * (
    (n00 + n10) * safe_log(1 - pi2)  + (n01 + n11) * safe_log(pi2) -
      n00 * safe_log(1 - pi01) - n01 * safe_log(pi01) -
      n10 * safe_log(1 - pi11) - n11 * safe_log(pi11)
  )
  
  LR_CC <- LR_UC + LR_IND
  
  list(LR_UC = LR_UC, LR_IND = LR_IND, LR_CC = LR_CC,
       p_UC  = 1 - pchisq(LR_UC,  df = 1),
       p_IND = 1 - pchisq(LR_IND, df = 1),
       p_CC  = 1 - pchisq(LR_CC,  df = 2),
       hit_rate = pi_hat, expected_rate = alpha,
       n_violations = n1, n_obs = T)
}

# Engle & Manganelli (2004) Dynamic Quantile test.
dynamic_quantile_test <- function(returns, VaR, alpha, lags = 4) {
  hit <- as.integer(returns < VaR) - alpha
  T   <- length(hit)
  p   <- lags
  if (T <= p + 2) stop("Too few observations for requested lags.")
  
  idx     <- (p + 1):T
  Hit_dep <- hit[idx]
  
  lag_mat <- matrix(NA_real_, nrow = length(idx), ncol = p)
  for (k in seq_len(p)) lag_mat[, k] <- hit[idx - k]
  
  X <- cbind(1, lag_mat, VaR[idx])
  colnames(X) <- c("(Intercept)", paste0("Hit_lag", seq_len(p)), "VaR")
  
  XtX_inv  <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% Hit_dep
  
  DQ_stat <- as.numeric(
    t(Hit_dep) %*% X %*% XtX_inv %*% t(X) %*% Hit_dep / (alpha * (1 - alpha))
  )
  
  df    <- ncol(X)
  p_val <- 1 - pchisq(DQ_stat, df = df)
  
  list(DQ_stat = DQ_stat, df = df, p_value = p_val,
       beta_hat = beta_hat, lags_used = lags)
}

# =============================================================================
#  CONFIGURATION
# =============================================================================
# This script estimates three GARCH families (sGARCH, gjrGARCH, eGARCH) with
# Student-t innovations across three regime counts (K = 1, 2, 3), yielding
# nine candidate specifications. 
# VaR is computed at alpha in {0.05, 0.01}, with
# 5 lags in the Engle-Manganelli Dynamic Quantile regression. These values
# are hardcoded in phases 4a/4b rather than parameterised here.


#Note: Paths are hardcoded. Adjust as needed for your environment

# Output directory — all PDFs and .tex files are written here
output_dir <- "C:/Users/2002j/Documents/Seminar/Plots"

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
cat("\n=== 1. Loading data ===\n")
raw <- read_excel("C:/Users/2002j/Documents/Seminar/Excel/Elpriser, Øst-Danmark.xlsx")
raw <- raw[, c("HourDK", "SpotPriceDKK")]
raw <- raw[!is.na(raw$SpotPriceDKK), ]

cat(sprintf("Hourly obs (non-NA): %d\n", nrow(raw)))
cat(sprintf("Date range: %s  to  %s\n",
            format(min(raw$HourDK)), format(max(raw$HourDK))))
cat(sprintf("Negative hourly prices: %d\n", sum(raw$SpotPriceDKK < 0)))


# ── 2. AGGREGATE TO DAILY MEAN ────────────────────────────────────────────────

cat("\n=== 2. Aggregate to daily mean ===\n")

raw$Date <- as.Date(raw$HourDK)
daily    <- aggregate(SpotPriceDKK ~ Date, data = raw, FUN = mean)
daily    <- daily[order(daily$Date), ]

cat(sprintf("Hourly observations         : %d\n", nrow(raw)))
cat(sprintf("Daily observations          : %d\n", nrow(daily)))
cat(sprintf("Negative daily-mean prices  : %d (%.2f%%)\n",
            sum(daily$SpotPriceDKK < 0),
            100 * mean(daily$SpotPriceDKK < 0)))
cat(sprintf("Min daily mean price        : %.2f DKK/MWh\n",
            min(daily$SpotPriceDKK)))
cat(sprintf("Max daily mean price        : %.2f DKK/MWh\n",
            max(daily$SpotPriceDKK)))


# ── 3. HANDLE NEGATIVE DAILY-MEAN PRICES (FORWARD-FILL) ──────────────────────
# Negative or sub-1-DKK daily means are flagged and replaced with the
# previous valid price. The first observation, if flagged, falls back to
# the sample median.

cat("\n=== 3. Forward-fill non-positive daily means ===\n")

PRICE_FLOOR <- 1
bad_idx     <- which(daily$SpotPriceDKK < PRICE_FLOOR)
bad_dates   <- daily$Date[bad_idx]

cat(sprintf("Days flagged for forward-fill : %d\n", length(bad_idx)))
if (length(bad_idx) > 0) {
  cat("Affected dates:\n")
  print(bad_dates)
  
  daily$SpotPriceDKK[bad_idx] <- NA_real_
  for (i in seq_len(nrow(daily))) {
    if (is.na(daily$SpotPriceDKK[i]) && i > 1L)
      daily$SpotPriceDKK[i] <- daily$SpotPriceDKK[i - 1L]
  }
  if (is.na(daily$SpotPriceDKK[1L]))
    daily$SpotPriceDKK[1L] <- median(daily$SpotPriceDKK, na.rm = TRUE)
}

cat(sprintf("Remaining non-positive prices : %d\n",
            sum(daily$SpotPriceDKK <= 0)))
cat(sprintf("New min daily price           : %.2f DKK/MWh\n",
            min(daily$SpotPriceDKK)))


# ── 4. COMPUTE LOG-RETURNS ────────────────────────────────────────────────────

log_ret <- diff(log(daily$SpotPriceDKK))
dates   <- daily$Date[-1]

cat("\n=== 4. Log-returns ===\n")
cat(sprintf("Observations : %d\n",     length(log_ret)))
cat(sprintf("Mean         : %+.6f\n",  mean(log_ret)))
cat(sprintf("Std dev      : %.6f\n",   sd(log_ret)))
cat(sprintf("Min / Max    : %.4f  /  %.4f\n", min(log_ret), max(log_ret)))


# ── 4B. DESCRIPTIVE STATISTICS ────────────────────────────────────────────────

cat("\n=== 4B. Descriptive statistics ===\n")

# ── Summary table ──────────────────────────────────────────────────────────
make_desc <- function(x) c(
  n        = length(x),
  mean     = mean(x),
  median   = median(x),
  sd       = sd(x),
  skewness = moments::skewness(x),
  kurtosis = moments::kurtosis(x) - 3,   # excess kurtosis
  min      = min(x),
  max      = max(x)
)

desc_tbl <- rbind(
  "Daily mean price (DKK)" = make_desc(daily$SpotPriceDKK),
  "Log-returns"           = make_desc(log_ret)
)
cat("\n")
print(round(desc_tbl, 4))

# Storing for LaTeX export
desc_stats_tbl <- desc_tbl

# ── Statistical tests ───────────────────────────────────────────────────────
cat("\n── Normality, autocorrelation, ARCH, and stationarity tests ──\n")

# Jarque-Bera (H0: normal; large stat → heavy tails or skew)
jb_ret   <- tseries::jarque.bera.test(log_ret)
jb_price <- tseries::jarque.bera.test(log(daily$SpotPriceDKK))
cat(sprintf("  Jarque-Bera  log-returns  : JB = %9.1f,  p < %.4f\n",
            jb_ret$statistic,   jb_ret$p.value))
cat(sprintf("  Jarque-Bera  log-prices   : JB = %9.1f,  p < %.4f\n",
            jb_price$statistic, jb_price$p.value))

# Ljung-Box on returns (mean autocorrelation) and squared returns (ARCH)
lb_r10  <- Box.test(log_ret,    lag = 10, type = "Ljung-Box")
lb_r20  <- Box.test(log_ret,    lag = 20, type = "Ljung-Box")
lb_sq10 <- Box.test(log_ret^2,  lag = 10, type = "Ljung-Box")
lb_sq20 <- Box.test(log_ret^2,  lag = 20, type = "Ljung-Box")
cat(sprintf("  Ljung-Box  r    lag 10 : Q = %7.2f,  p = %.4f\n", lb_r10$statistic,  lb_r10$p.value))
cat(sprintf("  Ljung-Box  r    lag 20 : Q = %7.2f,  p = %.4f\n", lb_r20$statistic,  lb_r20$p.value))
cat(sprintf("  Ljung-Box  r²   lag 10 : Q = %7.2f,  p = %.4f\n", lb_sq10$statistic, lb_sq10$p.value))
cat(sprintf("  Ljung-Box  r²   lag 20 : Q = %7.2f,  p = %.4f\n", lb_sq20$statistic, lb_sq20$p.value))

# ARCH-LM 
arch_lm <- FinTS::ArchTest(log_ret, lags = 10)
cat(sprintf("  ARCH-LM  lag 10        : LM = %7.2f,  p = %.4f\n",
            arch_lm$statistic, arch_lm$p.value))

# ADF
adf_lp  <- tseries::adf.test(log(daily$SpotPriceDKK))
adf_lr  <- tseries::adf.test(log_ret)
cat(sprintf("  ADF  log-prices        : stat = %.4f,  p = %.4f\n",
            adf_lp$statistic, adf_lp$p.value))
cat(sprintf("  ADF  log-returns       : stat = %.4f,  p = %.4f\n",
            adf_lr$statistic, adf_lr$p.value))

# Storing test results for later LaTeX table
desc_tests <- list(
  jb_ret=jb_ret, jb_price=jb_price,
  lb_r10=lb_r10, lb_r20=lb_r20, lb_sq10=lb_sq10, lb_sq20=lb_sq20,
  arch_lm=arch_lm, adf_lp=adf_lp, adf_lr=adf_lr
)

# ── Fit Student-t to log-returns (used in histogram overlay) ────────────────
t_fit <- tryCatch(
  MASS::fitdistr(log_ret, "t"),
  error = function(e) NULL
)

# ── Individual PDF plots ─────────────────────────────────────────────────────

# (a) Daily mean price series
pdf(file.path(output_dir, "daily_price.pdf"), width = 7, height = 4.5)
plot(daily$Date, daily$SpotPriceDKK, type = "l", col = "steelblue", lwd = 0.7,
     main = "Daily Mean Electricity Spot Price  (DK East)",
     xlab = "Date", ylab = "DKK / MWh")
dev.off()

# (a2) Daily log-returns time series
lr_plot <- diff(log(daily$SpotPriceDKK))
d_plot  <- daily$Date[-1]

pdf(file.path(output_dir, "returns_timeseries.pdf"), width = 7, height = 4.5)
plot(d_plot, lr_plot, type = "l", col = "steelblue", lwd = 0.5,
     main = "Daily Log-Returns of Mean Electricity Spot Price (DK East)",
     xlab = "Date", ylab = "Log-Return")
abline(h = 0, col = "grey70", lty = 2)
dev.off()
cat("  Saved: returns_timeseries.pdf\n")

# (b) Histogram of log-returns + Normal + Student-t density overlay
pdf(file.path(output_dir, "returns_histogram.pdf"), width = 7, height = 4.5)
xl <- quantile(log_ret, c(0.001, 0.999))
hist(log_ret, breaks = 80, freq = FALSE, col = "grey85", border = "white",
     xlim = xl, main = "Histogram of Daily Log-Returns", xlab = "Log-Return")
xg <- seq(xl[1], xl[2], length.out = 500)
lines(xg, dnorm(xg, mean(log_ret), sd(log_ret)), col = "steelblue", lwd = 1.8)
if (!is.null(t_fit)) {
  m_t  <- t_fit$estimate["m"]
  s_t  <- t_fit$estimate["s"]
  df_t <- t_fit$estimate["df"]
  lines(xg, dt((xg - m_t) / s_t, df = df_t) / s_t, col = "firebrick",
        lwd = 1.8, lty = 2)
  legend("topright", legend = c("N(μ,σ²)", sprintf("t(df=%.1f)", df_t)),
         col = c("steelblue", "firebrick"), lty = 1:2, cex = 0.8, bg = "white")
} else {
  legend("topright", "N(μ,σ²)", col = "steelblue", lty = 1, cex = 0.8)
}
dev.off()

# (c) ACF + PACF of returns (up to lag 40)
pdf(file.path(output_dir, "acf_pacf_returns.pdf"), width = 7, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(log_ret,  lag.max = 40, main = "ACF — Log-Returns")
pacf(log_ret, lag.max = 40, main = "PACF — Log-Returns")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

# (d) ACF of squared returns
pdf(file.path(output_dir, "acf_squared_returns.pdf"), width = 7, height = 4.5)
acf(log_ret^2, lag.max = 40,
    main = "ACF of Squared Log-Returns  (evidence of ARCH effects)")
dev.off()

# (e) Normal Q-Q plot of log-returns
pdf(file.path(output_dir, "qq_returns.pdf"), width = 7, height = 4.5)
qqnorm(log_ret, pch = 1, cex = 0.3, col = "steelblue",
       main = "Normal Q-Q Plot of Log-Returns")
qqline(log_ret, col = "firebrick", lwd = 1.5)
dev.off()

cat("Descriptive plots saved to:", output_dir, "\n")


# ── 4C. PRE-FILTER LOG-RETURNS: WEEKDAY DUMMIES, HOLIDAY DUMMY, AR(1)-AR(7)

cat("\n=== 4C. Weekday / holiday / AR pre-filtering ===\n")

if (!require(timeDate, quietly = TRUE)) install.packages("timeDate")
library(timeDate)

# ── Six weekday dummies  (Monday = reference category) ───────────────────
wday_num <- as.integer(format(dates, "%u"))

# ── Danish public holiday dummy ───────────────────────────────────────────
years_in_data    <- as.integer(unique(format(dates, "%Y")))
dk_holiday_dates <- as.Date(character(0))

for (yr in years_in_data) {
  easter <- as.Date(timeDate::Easter(yr))
  
  fixed <- c(
    as.Date(sprintf("%d-01-01", yr)),   # New Year's Day
    as.Date(sprintf("%d-06-05", yr)),   # Grundlovsdag
    as.Date(sprintf("%d-12-24", yr)),   # Christmas Eve
    as.Date(sprintf("%d-12-25", yr)),   # Christmas Day
    as.Date(sprintf("%d-12-26", yr)),   # 2. Christmas Day
    as.Date(sprintf("%d-12-31", yr))    # New Year's Eve
  )
  
  movable <- c(
    easter - 3L,   # Maundy Thursday
    easter - 2L,   # Good Friday
    easter,        # Easter Sunday
    easter + 1L,   # Easter Monday
    easter + 39L,  # Kristi himmelfartsdag
    easter + 49L,  # Pentecost
    easter + 50L   # Whit Monday
  )
  
  # Store bededag
  bededag <- if (yr <= 2023) easter + 26L else as.Date(character(0))
  
  dk_holiday_dates <- c(dk_holiday_dates, fixed, movable, bededag)
}

dk_holiday_dates <- sort(unique(dk_holiday_dates))
is_holiday_full  <- dates %in% dk_holiday_dates

# ── Trim first 7 obs to give AR(7) a complete lag for every row ─ 4C onward.──
T_full   <- length(log_ret)
ok_idx   <- seq(8L, T_full)

lag_mat <- sapply(1:7, function(k) log_ret[ok_idx - k])
colnames(lag_mat) <- paste0("lag", 1:7, "_ret")

log_ret    <- log_ret[ok_idx]
dates      <- dates[ok_idx]
is_holiday <- is_holiday_full[ok_idx]
wday_w     <- wday_num[ok_idx]

is_tue <- wday_w == 2L
is_wed <- wday_w == 3L
is_thu <- wday_w == 4L
is_fri <- wday_w == 5L
is_sat <- wday_w == 6L
is_sun <- wday_w == 7L

# ── OLS pre-filter ─────────────────────────────────────────────────────────
ols_filter <- lm(log_ret ~ is_tue + is_wed + is_thu + is_fri + is_sat + is_sun +
                   is_holiday + lag_mat)

cat("\nOLS pre-filter coefficient table:\n")
print(summary(ols_filter)$coefficients)

# Save the raw trimmed series; overwrite log_ret with OLS residuals so all
# downstream code (sections 5+) automatically operates on the filtered series.
log_ret_raw <- log_ret
log_ret     <- residuals(ols_filter)

# ── Ljung-Box at lag 20: marginal effect of each filtering step ───────────
lb_raw     <- Box.test(log_ret_raw, lag = 20, type = "Ljung-Box")
lb_dummies <- Box.test(
  residuals(lm(log_ret_raw ~ is_tue + is_wed + is_thu + is_fri + is_sat + is_sun +
                 is_holiday)),
  lag = 20, type = "Ljung-Box"
)
lb_full    <- Box.test(log_ret, lag = 20, type = "Ljung-Box")

cat("\nLjung-Box (lag 20) — marginal effect of each filtering step:\n")
cat(sprintf("  (a) Raw (trimmed) series               : Q = %7.2f,  p = %.4f\n",
            lb_raw$statistic,     lb_raw$p.value))
cat(sprintf("  (b) Weekday + holiday dummies only     : Q = %7.2f,  p = %.4f\n",
            lb_dummies$statistic, lb_dummies$p.value))
cat(sprintf("  (c) Dummies + AR(1)-AR(7) residuals  : Q = %7.2f,  p = %.4f\n",
            lb_full$statistic,    lb_full$p.value))

# ── Dummy summary ──────────────────────────────────────────────────────────
n_by_day  <- tabulate(wday_w, nbins = 7L)
n_hol     <- sum(is_holiday)
n_overlap <- sum(is_holiday & wday_w %in% c(6L, 7L))

cat(sprintf("\nDummy summary (working sample: T = %d, trimmed 7 obs for AR(7)):\n",
            length(dates)))
cat(sprintf("  Mon (reference) : %d\n", n_by_day[1L]))
cat(sprintf("  Tue             : %d\n", n_by_day[2L]))
cat(sprintf("  Wed             : %d\n", n_by_day[3L]))
cat(sprintf("  Thu             : %d\n", n_by_day[4L]))
cat(sprintf("  Fri             : %d\n", n_by_day[5L]))
cat(sprintf("  Sat             : %d\n", n_by_day[6L]))
cat(sprintf("  Sun             : %d\n", n_by_day[7L]))
cat(sprintf("  Holidays (in sample)     : %d\n", n_hol))
cat(sprintf("  Holidays on a weekend    : %d  (overlap check — dummies non-collinear)\n",
            n_overlap))

# Verify all 15 coefficients are present
cat(sprintf("\nNumber of OLS coefficients fitted: %d (expected: 15)\n",
            length(coef(ols_filter))))
stopifnot(length(coef(ols_filter)) == 15L)


# ── 4D. RESIDUAL ACF/PACF DIAGNOSTICS ──────────────────────────────────────

pdf(file.path(output_dir, "acf_pacf_residuals.pdf"), width = 7, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(log_ret,  lag.max = 40, main = "ACF — OLS residuals")
pacf(log_ret, lag.max = 40, main = "PACF — OLS residuals")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

pdf(file.path(output_dir, "acf_squared_residuals.pdf"), width = 7, height = 4.5)
acf(log_ret^2, lag.max = 40,
    main = "ACF of Squared OLS Residuals (input to MS-GARCH)")
dev.off()

cat("Residual ACF/PACF plots saved.\n")


# ── 5. PHASE 3a: SINGLE-REGIME BENCHMARK ─────────────────────────────────────
# Fit three K=1 GARCH variants with Student-t errors. best becomes baseline

cat("\n=== 5. Phase 3a — Single-regime benchmark (K=1) ===\n")

bench_models <- c("sGARCH", "gjrGARCH", "eGARCH")
bench_fits   <- list()
bench_table  <- data.frame(
  GARCH  = character(),
  LogLik = numeric(),
  nPar   = integer(),
  AIC    = numeric(),
  BIC    = numeric(),
  stringsAsFactors = FALSE
)

for (gm in bench_models) {
  cat(sprintf("  Fitting %s_std_K1 ... ", gm))
  
  spec_i <- CreateSpec(
    variance.spec     = list(model        = gm),
    distribution.spec = list(distribution = "std")
  )
  
  set.seed(1234)
  fit_i <- tryCatch(
    FitML(spec = spec_i, data = log_ret),
    error = function(e) { cat("FAILED\n"); NULL }
  )
  
  if (is.null(fit_i)) next
  
  bench_fits[[gm]] <- fit_i
  bench_table <- rbind(bench_table, data.frame(
    GARCH  = gm,
    LogLik = round(fit_i$loglik, 2),
    nPar   = length(fit_i$par),
    AIC    = round(AIC(fit_i), 2),
    BIC    = round(BIC(fit_i), 2),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("LL=%.2f  BIC=%.2f\n", fit_i$loglik, BIC(fit_i)))
}

bench_table <- bench_table[order(bench_table$BIC), ]
cat("\nBenchmark ranking (sorted by BIC):\n")
print(bench_table, row.names = FALSE)

best_k1_name <- bench_table$GARCH[1]
best_k1_fit  <- bench_fits[[best_k1_name]]
cat(sprintf("\nBest K=1 model: %s\n", best_k1_name))

# Show parameters with standard errors
cat("\nParameter estimates for best K=1 model:\n")
print(summary(best_k1_fit))

# ── LB tests on K=1 standardized residuals ──
vol_k1 <- Volatility(best_k1_fit)
eta_k1 <- log_ret / pmax(vol_k1, 1e-8)
lb_eta_10  <- Box.test(eta_k1,   lag = 10, type = "Ljung-Box")
lb_eta_20  <- Box.test(eta_k1,   lag = 20, type = "Ljung-Box")
lb_eta2_10 <- Box.test(eta_k1^2, lag = 10, type = "Ljung-Box")
lb_eta2_20 <- Box.test(eta_k1^2, lag = 20, type = "Ljung-Box")

cat(sprintf("K=1 LB eta   lag 10: Q=%7.2f, p=%.4f\n", lb_eta_10$statistic,  lb_eta_10$p.value))
cat(sprintf("K=1 LB eta   lag 20: Q=%7.2f, p=%.4f\n", lb_eta_20$statistic,  lb_eta_20$p.value))
cat(sprintf("K=1 LB eta²  lag 10: Q=%7.2f, p=%.4f\n", lb_eta2_10$statistic, lb_eta2_10$p.value))
cat(sprintf("K=1 LB eta²  lag 20: Q=%7.2f, p=%.4f\n", lb_eta2_20$statistic, lb_eta2_20$p.value))

phase3a_results <- list(
  benchmark_table = bench_table,
  best_k1_name    = best_k1_name,
  best_k1_fit     = best_k1_fit,
  k1_coef         = summary(best_k1_fit)$estimate,
  k1_loglik       = best_k1_fit$loglik,
  k1_AIC          = AIC(best_k1_fit),
  k1_BIC          = BIC(best_k1_fit),
  k1_resid_tests  = list(
    lb_eta10  = lb_eta_10,
    lb_eta20  = lb_eta_20,
    lb_eta210 = lb_eta2_10,
    lb_eta220 = lb_eta2_20
  )
)

# ── 6. PHASE 3b: TWO-REGIME GRID WITH MULTIPLE STARTS ────────────────────────
# Fit three K=2 GARCH variants with multiple random seeds

cat("\n=== 6. Phase 3b — Two-regime grid (K=2) ===\n")

bench_k2_models <- c("sGARCH", "gjrGARCH", "eGARCH")
n_starts        <- 8
fits_k2         <- list()
bench_k2_table  <- data.frame(
  GARCH       = character(),
  LogLik      = numeric(),
  nPar        = integer(),
  AIC         = numeric(),
  BIC         = numeric(),
  P11         = numeric(),
  P22         = numeric(),
  Degenerate  = logical(),
  stringsAsFactors = FALSE
)

for (gm in bench_k2_models) {
  cat(sprintf("\n  Fitting %s_std_K2 with %d starts:\n", gm, n_starts))
  
  spec_i <- CreateSpec(
    variance.spec     = list(model        = rep(gm, 2)),
    distribution.spec = list(distribution = rep("std", 2))
  )
  
  best_loglik <- -Inf
  best_fit    <- NULL
  
  for (s in seq_len(n_starts)) {
    set.seed(1234 + s)
    fit_try <- tryCatch(
      FitML(spec = spec_i, data = log_ret),
      error = function(e) NULL
    )
    if (is.null(fit_try)) {
      cat(sprintf("    seed %d: FAILED\n", 1234 + s))
      next
    }
    cat(sprintf("    seed %d: LL = %.2f\n", 1234 + s, fit_try$loglik))
    if (fit_try$loglik > best_loglik) {
      best_loglik <- fit_try$loglik
      best_fit    <- fit_try
    }
  }
  
  if (is.null(best_fit)) {
    cat(sprintf("  All seeds failed for %s_std_K2\n", gm))
    next
  }
  
  # Extracts transition probabilities to check for absorbing regimes
  P11 <- as.numeric(best_fit$par["P_1_1"])
  P22 <- 1 - as.numeric(best_fit$par["P_2_1"])
  is_degenerate <- (P11 >= 0.999) || (P22 >= 0.999)
  
  fits_k2[[gm]] <- best_fit
  bench_k2_table <- rbind(bench_k2_table, data.frame(
    GARCH      = gm,
    LogLik     = round(best_fit$loglik, 2),
    nPar       = length(best_fit$par),
    AIC        = round(AIC(best_fit), 2),
    BIC        = round(BIC(best_fit), 2),
    P11        = round(P11, 4),
    P22        = round(P22, 4),
    Degenerate = is_degenerate,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  Best %s_std_K2: LL=%.2f  BIC=%.2f  P11=%.3f  P22=%.3f  %s\n",
              gm, best_fit$loglik, BIC(best_fit), P11, P22,
              ifelse(is_degenerate, "[DEGENERATE]", "[OK]")))
}

bench_k2_table <- bench_k2_table[order(bench_k2_table$BIC), ]
cat("\nK=2 ranking (sorted by BIC):\n")
print(bench_k2_table, row.names = FALSE)

# picks best non-degenerate K=2 fit
valid_k2 <- bench_k2_table[!bench_k2_table$Degenerate, ]

if (nrow(valid_k2) == 0) {
  cat("\nNO valid K=2 model converged. Phase 3b FAILS.\n")
  cat("Fall back to K=1 benchmark (gjrGARCH_std_K1).\n")
  best_k2_name <- NA
  best_k2_fit  <- NULL
  phase3b_pass <- FALSE
} else {
  best_k2_name <- valid_k2$GARCH[1]
  best_k2_fit  <- fits_k2[[best_k2_name]]
  bic_improvement <- phase3a_results$k1_BIC - valid_k2$BIC[1]
  
  cat(sprintf("\nBest non-degenerate K=2 model: %s\n", best_k2_name))
  cat(sprintf("BIC improvement over K=1: %.2f points\n", bic_improvement))
  
  if (bic_improvement >= 10) {
    cat("→ K=2 model justified. Phase 3b PASSES.\n")
    phase3b_pass <- TRUE
  } else {
    cat("→ BIC improvement < 10 points. K=2 not justified over K=1.\n")
    phase3b_pass <- FALSE
  }
}

# Show parameters with SEs for best K=2
if (!is.null(best_k2_fit)) {
  cat(sprintf("\nParameter estimates for %s_std_K2:\n", best_k2_name))
  print(summary(best_k2_fit))
}

# ── 6c. MOVING BLOCK BOOTSTRAP STANDARD ERRORS FOR K=2 EGARCH (PARALLEL) ──
cat("\n=== 6c. Moving block bootstrap SEs for K=2 eGARCH (parallel) ===\n")

library(parallel)

n_boot     <- 200
block_size <- 30
T_obs      <- length(log_ret)
n_blocks   <- ceiling(T_obs / block_size)
n_par_k2   <- length(best_k2_fit$par)
par_names_k2 <- names(best_k2_fit$par)

# ── Testing one replication ──
cat("Testing one sequential replication first...\n")
test_one <- function(b) {
  set.seed(20000 + b)
  starts <- sample.int(T_obs - block_size + 1, n_blocks, replace = TRUE)
  idx_b  <- unlist(lapply(starts, function(s) s:(s + block_size - 1)))
  idx_b  <- idx_b[seq_len(T_obs)]
  data_b <- log_ret[idx_b]
  
  spec_local <- MSGARCH::CreateSpec(
    variance.spec     = list(model        = rep("eGARCH", 2)),
    distribution.spec = list(distribution = rep("std", 2))
  )
  
  fit_b <- MSGARCH::FitML(spec = spec_local, data = data_b)
  fit_b$par
}
test_result <- tryCatch(test_one(1), error = function(e) {
  cat("Sequential test FAILED:", conditionMessage(e), "\n")
  NULL
})
if (is.null(test_result)) stop("Sequential test failed - stopping before parallel run")
cat("Sequential test OK. Starting parallel bootstrap...\n\n")

# ──parallel ──
n_cores <- N_CORES
cat(sprintf("Using %d parallel workers\n", n_cores))

cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(MSGARCH))
clusterExport(cl, c("log_ret", "block_size", "T_obs", "n_blocks", "n_par_k2"))

boot_one <- function(b) {
  set.seed(20000 + b)
  starts <- sample.int(T_obs - block_size + 1, n_blocks, replace = TRUE)
  idx_b  <- unlist(lapply(starts, function(s) s:(s + block_size - 1)))
  idx_b  <- idx_b[seq_len(T_obs)]
  data_b <- log_ret[idx_b]
  
  # Build the spec fresh on the worker
  spec_local <- MSGARCH::CreateSpec(
    variance.spec     = list(model        = rep("eGARCH", 2)),
    distribution.spec = list(distribution = rep("std", 2))
  )
  
  fit_b <- tryCatch(
    MSGARCH::FitML(spec = spec_local, data = data_b),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(rep(NA_real_, n_par_k2))
  fit_b$par
}

t0_boot <- proc.time()["elapsed"]
boot_list <- parLapplyLB(cl, seq_len(n_boot), boot_one)
stopCluster(cl)
t1_boot <- proc.time()["elapsed"]

cat(sprintf("\nBootstrap complete: %.1f minutes total\n", (t1_boot - t0_boot) / 60))

boot_par <- do.call(rbind, boot_list)
colnames(boot_par) <- par_names_k2

boot_par_valid <- boot_par[complete.cases(boot_par), , drop = FALSE]
n_valid        <- nrow(boot_par_valid)
cat(sprintf("Successful replications: %d / %d\n", n_valid, n_boot))

boot_se   <- apply(boot_par_valid, 2, sd)
boot_mean <- apply(boot_par_valid, 2, mean)

robust_inference <- data.frame(
  Estimate    = round(best_k2_fit$par, 6),
  Boot_SE     = round(boot_se, 6),
  Boot_Mean   = round(boot_mean, 6),
  t_stat      = round(best_k2_fit$par / boot_se, 3),
  CI_low_95   = round(apply(boot_par_valid, 2, quantile, 0.025), 4),
  CI_high_95  = round(apply(boot_par_valid, 2, quantile, 0.975), 4)
)
print(robust_inference)


# ── 6c-post. Label-switching correction + robust inference ──
cat("\n=== 6c-post. Label-switching correction ===\n")

# ── Helper: swap regime labels in a K=2 MSGARCH parameter vector ──
swap_regimes_k2 <- function(par_vec) {
  out <- par_vec
  for (base in c("alpha0", "alpha1", "alpha2", "beta", "nu")) {
    i1 <- paste0(base, "_1"); i2 <- paste0(base, "_2")
    out[i1] <- par_vec[i2]
    out[i2] <- par_vec[i1]
  }
  # Transition-matrix relabeling under regime-swap:
  #   new P_1_1 = old P_2_2 = 1 - old P_2_1
  #   new P_2_1 = old P_1_2 = 1 - old P_1_1
  out["P_1_1"] <- 1 - par_vec["P_2_1"]
  out["P_2_1"] <- 1 - par_vec["P_1_1"]
  out
}

id_rule <- "alpha2"   

needs_swap <- function(par_vec, rule = id_rule) {
  switch(rule,
         alpha2 = par_vec["alpha2_1"] < par_vec["alpha2_2"],  
         beta   = par_vec["beta_1"]   > par_vec["beta_2"],    
         nu     = par_vec["nu_1"]     < par_vec["nu_2"]      
  )
}

# Apply swap per replication.
boot_par_aligned <- boot_par_valid
n_swapped <- 0L
for (i in seq_len(nrow(boot_par_aligned))) {
  if (isTRUE(unname(needs_swap(boot_par_aligned[i, ])))) {
    boot_par_aligned[i, ] <- swap_regimes_k2(boot_par_aligned[i, ])
    n_swapped <- n_swapped + 1L
  }
}
cat(sprintf("Label-switch corrections: %d / %d (%.1f%%) using rule '%s'\n",
            n_swapped, nrow(boot_par_aligned),
            100 * n_swapped / nrow(boot_par_aligned), id_rule))

# ── Consistency check ──
agree <- data.frame(
  alpha2 = vapply(seq_len(nrow(boot_par_valid)),
                  function(i) isTRUE(unname(needs_swap(boot_par_valid[i, ], "alpha2"))), logical(1)),
  beta   = vapply(seq_len(nrow(boot_par_valid)),
                  function(i) isTRUE(unname(needs_swap(boot_par_valid[i, ], "beta"))),   logical(1)),
  nu     = vapply(seq_len(nrow(boot_par_valid)),
                  function(i) isTRUE(unname(needs_swap(boot_par_valid[i, ], "nu"))),     logical(1))
)
cat("\nSwap-rule agreement (TRUE = rule says swap):\n")
print(colSums(agree))
cat(sprintf("All three rules agree in %d / %d replications\n",
            sum(rowSums(agree) %in% c(0, 3)), nrow(agree)))

# ── Diagnostics: histograms before/after ──
op <- par(mfrow = c(3, 2), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
for (p in c("beta_1", "beta_2", "alpha2_1")) {
  hist(boot_par_valid[, p],   breaks = 30, main = paste(p, "(raw)"),     xlab = "", col = "grey80")
  abline(v = best_k2_fit$par[p], col = "red", lwd = 2)
  hist(boot_par_aligned[, p], breaks = 30, main = paste(p, "(aligned)"), xlab = "", col = "steelblue")
  abline(v = best_k2_fit$par[p], col = "red", lwd = 2)
}
par(op)

# ── Robust inference: median + MAD as primary, mean + SD as secondary ──
boot_med  <- apply(boot_par_aligned, 2, median)
boot_mad  <- apply(boot_par_aligned, 2, mad)         
boot_sd   <- apply(boot_par_aligned, 2, sd)
boot_q025 <- apply(boot_par_aligned, 2, quantile, 0.025)
boot_q975 <- apply(boot_par_aligned, 2, quantile, 0.975)

robust_inference_aligned <- data.frame(
  Estimate    = round(best_k2_fit$par, 4),
  Boot_Median = round(boot_med,  4),
  Boot_MAD    = round(boot_mad,  4),
  Boot_SD     = round(boot_sd,   4),
  t_MAD       = round(best_k2_fit$par / boot_mad, 2),
  CI_low_95   = round(boot_q025, 4),
  CI_high_95  = round(boot_q975, 4)
)
cat("\n── Aligned bootstrap inference (median/MAD-based) ──\n")
print(robust_inference_aligned)

# ── Flag boundary-hits on nu ──
nu_at_bound_1 <- mean(boot_par_aligned[, "nu_1"] > 99)
nu_at_bound_2 <- mean(boot_par_aligned[, "nu_2"] > 99)
cat(sprintf("\nnu_1 hitting upper bound (>99): %.1f%% of replications\n", 100 * nu_at_bound_1))
cat(sprintf("nu_2 hitting upper bound (>99): %.1f%% of replications\n", 100 * nu_at_bound_2))

# ── Build phase3b_results from current workspace ──
phase3b_results <- list(
  k2_grid_table  = bench_k2_table,
  best_k2_name   = best_k2_name,
  best_k2_fit    = best_k2_fit,
  phase3b_pass   = TRUE,
  bic_improvement_over_k1 = phase3a_results$k1_BIC - bench_k2_table$BIC[1],
  bootstrap_n_reps     = n_boot,
  bootstrap_block_size = block_size,
  bootstrap_n_valid    = n_valid,
  robust_inference         = robust_inference,
  robust_inference_aligned = robust_inference_aligned,
  label_switch_n_swapped   = n_swapped,
  label_switch_rule        = id_rule
)
cat("phase3b_results constructed.\n")

# ── 7. PHASE 3c: THREE-REGIME GRID ───────────────────────────
# Fit K=3 variants to test whether a third regime adds explanatory value.

cat("\n=== 7. Phase 3c — Three-regime grid (K=3) ===\n")

bench_k3_models <- c("sGARCH", "gjrGARCH", "eGARCH")
n_starts_k3     <- 5
fits_k3         <- list()
bench_k3_table  <- data.frame(
  GARCH       = character(),
  LogLik      = numeric(),
  nPar        = integer(),
  AIC         = numeric(),
  BIC         = numeric(),
  P11         = numeric(),
  P22         = numeric(),
  P33         = numeric(),
  Min_ergodic = numeric(),
  Degenerate  = logical(),
  stringsAsFactors = FALSE
)

# Helper: ergodic distribution for K-regime tansition matrix
ergodic_from_par <- function(par, K) {
  P <- matrix(0, K, K)
  for (k in seq_len(K)) {
    for (j in seq_len(K - 1)) {
      nm <- sprintf("P_%d_%d", k, j)
      P[k, j] <- as.numeric(par[nm])
    }
    P[k, K] <- 1 - sum(P[k, seq_len(K - 1)])
  }
  A  <- t(P) - diag(K)
  A  <- rbind(A, rep(1, K))
  b  <- c(rep(0, K), 1)
  as.numeric(qr.solve(A, b))
}

for (gm in bench_k3_models) {
  cat(sprintf("\n  Fitting %s_std_K3 with %d starts:\n", gm, n_starts_k3))
  
  spec_i <- CreateSpec(
    variance.spec     = list(model        = rep(gm, 3)),
    distribution.spec = list(distribution = rep("std", 3))
  )
  
  best_loglik <- -Inf
  best_fit    <- NULL
  
  for (s in seq_len(n_starts_k3)) {
    set.seed(2000 + s)
    fit_try <- tryCatch(
      FitML(spec = spec_i, data = log_ret),
      error = function(e) NULL
    )
    if (is.null(fit_try)) {
      cat(sprintf("    seed %d: FAILED\n", 2000 + s))
      next
    }
    cat(sprintf("    seed %d: LL = %.2f\n", 2000 + s, fit_try$loglik))
    if (fit_try$loglik > best_loglik) {
      best_loglik <- fit_try$loglik
      best_fit    <- fit_try
    }
  }
  
  if (is.null(best_fit)) {
    cat(sprintf("  All seeds failed for %s_std_K3\n", gm))
    next
  }
  
  # Compute diagonal of transition matrix
  P11 <- as.numeric(best_fit$par["P_1_1"])
  P22 <- as.numeric(best_fit$par["P_2_2"])
  P33 <- 1 - as.numeric(best_fit$par["P_3_1"]) - as.numeric(best_fit$par["P_3_2"])
  
  # Ergodic distribution
  erg <- tryCatch(ergodic_from_par(best_fit$par, 3),
                  error = function(e) rep(NA, 3))
  min_erg <- min(erg, na.rm = TRUE)
  
  # Degenerate if any diagonal element >= 0.999 OR any regime has ergodic prob < 0.02
  is_degenerate <- any(c(P11, P22, P33) >= 0.999, na.rm = TRUE) ||
    (is.finite(min_erg) && min_erg < 0.02)
  
  fits_k3[[gm]] <- best_fit
  bench_k3_table <- rbind(bench_k3_table, data.frame(
    GARCH       = gm,
    LogLik      = round(best_fit$loglik, 2),
    nPar        = length(best_fit$par),
    AIC         = round(AIC(best_fit), 2),
    BIC         = round(BIC(best_fit), 2),
    P11         = round(P11, 4),
    P22         = round(P22, 4),
    P33         = round(P33, 4),
    Min_ergodic = round(min_erg, 4),
    Degenerate  = is_degenerate,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  Best %s_std_K3: LL=%.2f  BIC=%.2f  diag=(%.3f, %.3f, %.3f)  min_erg=%.3f  %s\n",
              gm, best_fit$loglik, BIC(best_fit), P11, P22, P33, min_erg,
              ifelse(is_degenerate, "[DEGENERATE]", "[OK]")))
}

bench_k3_table <- bench_k3_table[order(bench_k3_table$BIC), ]
cat("\nK=3 ranking (sorted by BIC):\n")
print(bench_k3_table, row.names = FALSE)

# Decision logic
valid_k3 <- bench_k3_table[!bench_k3_table$Degenerate, ]

if (nrow(valid_k3) == 0) {
  cat("\nAll K=3 models degenerate. K=2 remains the preferred specification.\n")
  best_k3_name <- NA
  best_k3_fit  <- NULL
  phase3c_pass <- FALSE
} else {
  k2_bic       <- bench_k2_table$BIC[1]
  bic_gain     <- k2_bic - valid_k3$BIC[1]
  cat(sprintf("\nBest non-degenerate K=3 model: %s\n", valid_k3$GARCH[1]))
  cat(sprintf("BIC improvement over K=2: %.2f points\n", bic_gain))
  
  if (bic_gain >= 20) {
    cat("→ K=3 justifies the extra complexity. Phase 3c PASSES.\n")
    phase3c_pass <- TRUE
    best_k3_name <- valid_k3$GARCH[1]
    best_k3_fit  <- fits_k3[[best_k3_name]]
  } else {
    cat("→ BIC improvement < 20 points. K=2 remains preferred specification.\n")
    phase3c_pass <- FALSE
    best_k3_name <- valid_k3$GARCH[1]   # Keep for documentation
    best_k3_fit  <- fits_k3[[best_k3_name]]
  }
}

# Always show transition matrix of best K=3 fit for documentation
if (!is.null(best_k3_fit)) {
  cat(sprintf("\nTransition matrix for %s_std_K3 (for documentation):\n",
              best_k3_name))
  K_t <- 3
  P_mat <- matrix(0, K_t, K_t)
  for (k in seq_len(K_t)) {
    for (j in seq_len(K_t - 1)) {
      P_mat[k, j] <- as.numeric(best_k3_fit$par[sprintf("P_%d_%d", k, j)])
    }
    P_mat[k, K_t] <- 1 - sum(P_mat[k, seq_len(K_t - 1)])
  }
  rownames(P_mat) <- colnames(P_mat) <- paste0("k=", seq_len(K_t))
  print(round(P_mat, 4))
}

# Save phase 3c results
phase3c_results <- list(
  k3_grid_table  = bench_k3_table,
  best_k3_name   = best_k3_name,
  best_k3_fit    = best_k3_fit,
  phase3c_pass   = phase3c_pass
)

cat("\nParameter estimates for best K=3 model:\n")
print(summary(best_k3_fit))


# ── 7b. BLOCK BOOTSTRAP STANDARD ERRORS FOR K=3 EGARCH (PARALLEL) ──
cat("\n=== 7b. Block bootstrap SEs for K=3 eGARCH (parallel) ===\n")

library(parallel)

n_boot_k3     <- 200
block_size_k3 <- 30
T_obs         <- length(log_ret)
n_blocks_k3   <- ceiling(T_obs / block_size_k3)
n_par_k3      <- length(best_k3_fit$par)
par_names_k3  <- names(best_k3_fit$par)

# ── One test replication ──
cat("Testing one sequential replication first...\n")
test_one_k3 <- function(b) {
  set.seed(30000 + b)
  starts <- sample.int(T_obs - block_size_k3 + 1, n_blocks_k3, replace = TRUE)
  idx_b  <- unlist(lapply(starts, function(s) s:(s + block_size_k3 - 1)))
  idx_b  <- idx_b[seq_len(T_obs)]
  data_b <- log_ret[idx_b]
  
  spec_local <- MSGARCH::CreateSpec(
    variance.spec     = list(model        = rep("eGARCH", 3)),
    distribution.spec = list(distribution = rep("std", 3))
  )
  
  fit_b <- MSGARCH::FitML(spec = spec_local, data = data_b)
  fit_b$par
}
test_result_k3 <- tryCatch(test_one_k3(1), error = function(e) {
  cat("Sequential test FAILED:", conditionMessage(e), "\n")
  NULL
})
if (is.null(test_result_k3)) stop("Sequential test failed - stopping before parallel run")
cat("Sequential test OK. Starting parallel bootstrap...\n\n")

# ── Parallel ──
n_cores <- N_CORES
cat(sprintf("Using %d parallel workers\n", n_cores))

cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(MSGARCH))
clusterExport(cl, c("log_ret", "block_size_k3", "T_obs", "n_blocks_k3", "n_par_k3"))

boot_one_k3 <- function(b) {
  set.seed(30000 + b)
  starts <- sample.int(T_obs - block_size_k3 + 1, n_blocks_k3, replace = TRUE)
  idx_b  <- unlist(lapply(starts, function(s) s:(s + block_size_k3 - 1)))
  idx_b  <- idx_b[seq_len(T_obs)]
  data_b <- log_ret[idx_b]
  
  spec_local <- MSGARCH::CreateSpec(
    variance.spec     = list(model        = rep("eGARCH", 3)),
    distribution.spec = list(distribution = rep("std", 3))
  )
  
  fit_b <- tryCatch(
    MSGARCH::FitML(spec = spec_local, data = data_b),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(rep(NA_real_, n_par_k3))
  fit_b$par
}

t0_boot_k3 <- proc.time()["elapsed"]
boot_list_k3 <- parLapplyLB(cl, seq_len(n_boot_k3), boot_one_k3)
stopCluster(cl)
t1_boot_k3 <- proc.time()["elapsed"]

cat(sprintf("\nBootstrap complete: %.1f minutes total\n", (t1_boot_k3 - t0_boot_k3) / 60))

boot_par_k3 <- do.call(rbind, boot_list_k3)
colnames(boot_par_k3) <- par_names_k3

boot_par_k3_valid <- boot_par_k3[complete.cases(boot_par_k3), , drop = FALSE]
n_valid_k3        <- nrow(boot_par_k3_valid)
cat(sprintf("Successful replications: %d / %d\n", n_valid_k3, n_boot_k3))

boot_se_k3   <- apply(boot_par_k3_valid, 2, sd)
boot_mean_k3 <- apply(boot_par_k3_valid, 2, mean)

robust_inference_k3 <- data.frame(
  Estimate    = round(best_k3_fit$par, 6),
  Boot_SE     = round(boot_se_k3, 6),
  Boot_Mean   = round(boot_mean_k3, 6),
  t_stat      = round(best_k3_fit$par / boot_se_k3, 3),
  CI_low_95   = round(apply(boot_par_k3_valid, 2, quantile, 0.025), 4),
  CI_high_95  = round(apply(boot_par_k3_valid, 2, quantile, 0.975), 4)
)
print(robust_inference_k3)

# ── 7b-post. Label-switching correction for K=3 ──
cat("\n=== 7b-post. Label-switching correction for K=3 ===\n")

# ── Swap regime labels in K=3 MSGARCH parameter vector ──

permute_regimes_k3 <- function(par_vec, perm) {
  out <- par_vec
  
  # Within-regime parameters: just relabel
  for (base in c("alpha0", "alpha1", "alpha2", "beta", "nu")) {
    for (k_new in 1:3) {
      k_old <- perm[k_new]
      out[paste0(base, "_", k_new)] <- par_vec[paste0(base, "_", k_old)]
    }
  }
  
  # Transition matrix: rebuild full 3x3 first, permute, then write back
  P_old <- matrix(0, 3, 3)
  for (i in 1:3) {
    for (j in 1:2) P_old[i, j] <- par_vec[sprintf("P_%d_%d", i, j)]
    P_old[i, 3] <- 1 - sum(P_old[i, 1:2])
  }
  # P_new[i_new, j_new] = P_old[perm[i_new], perm[j_new]]
  P_new <- P_old[perm, perm]
  for (i in 1:3) {
    for (j in 1:2) out[sprintf("P_%d_%d", i, j)] <- P_new[i, j]
  }
  out
}

# ── Distance to the point estimate's key signature ──
beta_mle <- c(best_k3_fit$par["beta_1"],
              best_k3_fit$par["beta_2"],
              best_k3_fit$par["beta_3"])

# All 6 permutations of (1,2,3)
all_perms <- list(
  c(1,2,3), c(1,3,2), c(2,1,3),
  c(2,3,1), c(3,1,2), c(3,2,1)
)

# ── Apply label-switching correction ──
boot_par_k3_aligned <- boot_par_k3_valid
n_swapped_k3 <- 0L
perm_used    <- integer(nrow(boot_par_k3_aligned))

for (i in seq_len(nrow(boot_par_k3_aligned))) {
  par_i <- boot_par_k3_aligned[i, ]
  best_dist <- Inf
  best_perm <- c(1, 2, 3)
  
  for (p_idx in seq_along(all_perms)) {
    perm <- all_perms[[p_idx]]
    beta_perm <- c(par_i[paste0("beta_", perm[1])],
                   par_i[paste0("beta_", perm[2])],
                   par_i[paste0("beta_", perm[3])])
    dist <- sum((beta_perm - beta_mle)^2)
    if (dist < best_dist) {
      best_dist <- dist
      best_perm <- perm
      perm_used[i] <- p_idx
    }
  }
  
  if (!identical(best_perm, c(1L, 2L, 3L)) && !identical(best_perm, c(1, 2, 3))) {
    boot_par_k3_aligned[i, ] <- permute_regimes_k3(par_i, best_perm)
    n_swapped_k3 <- n_swapped_k3 + 1L
  }
}

cat(sprintf("Label-switch corrections: %d / %d (%.1f%%) replications permuted\n",
            n_swapped_k3, nrow(boot_par_k3_aligned),
            100 * n_swapped_k3 / nrow(boot_par_k3_aligned)))

cat("\nDistribution of chosen permutations:\n")
perm_labels <- c("(1,2,3)=identity", "(1,3,2)", "(2,1,3)",
                 "(2,3,1)", "(3,1,2)", "(3,2,1)")
print(table(factor(perm_used, levels = 1:6, labels = perm_labels)))

# ── Histograms before/after for the most exposed parameters ──
op <- par(mfrow = c(3, 2), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
for (p in c("beta_1", "beta_2", "beta_3")) {
  hist(boot_par_k3_valid[, p],   breaks = 30,
       main = paste(p, "(raw)"), xlab = "", col = "grey80")
  abline(v = best_k3_fit$par[p], col = "red", lwd = 2)
  hist(boot_par_k3_aligned[, p], breaks = 30,
       main = paste(p, "(aligned)"), xlab = "", col = "steelblue")
  abline(v = best_k3_fit$par[p], col = "red", lwd = 2)
}
par(op)

# ── Robust inference ──
boot_med_k3  <- apply(boot_par_k3_aligned, 2, median)
boot_mad_k3  <- apply(boot_par_k3_aligned, 2, mad)
boot_sd_k3   <- apply(boot_par_k3_aligned, 2, sd)
boot_q025_k3 <- apply(boot_par_k3_aligned, 2, quantile, 0.025)
boot_q975_k3 <- apply(boot_par_k3_aligned, 2, quantile, 0.975)

robust_inference_k3_aligned <- data.frame(
  Estimate    = round(best_k3_fit$par, 4),
  Boot_Median = round(boot_med_k3,  4),
  Boot_MAD    = round(boot_mad_k3,  4),
  Boot_SD     = round(boot_sd_k3,   4),
  t_MAD       = round(best_k3_fit$par / boot_mad_k3, 2),
  CI_low_95   = round(boot_q025_k3, 4),
  CI_high_95  = round(boot_q975_k3, 4)
)
cat("\n── Aligned K=3 bootstrap inference (median/MAD-based) ──\n")
print(robust_inference_k3_aligned)

# ── Boundary-hits on nu ──
nu_at_bound_k3 <- c(
  mean(boot_par_k3_aligned[, "nu_1"] > 99),
  mean(boot_par_k3_aligned[, "nu_2"] > 99),
  mean(boot_par_k3_aligned[, "nu_3"] > 99)
)
cat(sprintf("\nnu_1 hitting upper bound (>99): %.1f%% of replications\n", 100 * nu_at_bound_k3[1]))
cat(sprintf("nu_2 hitting upper bound (>99): %.1f%% of replications\n", 100 * nu_at_bound_k3[2]))
cat(sprintf("nu_3 hitting upper bound (>99): %.1f%% of replications\n", 100 * nu_at_bound_k3[3]))

# ── LB tests on K=3 standardized residuals ──
vol_k3 <- Volatility(best_k3_fit)
eta_k3 <- log_ret / pmax(vol_k3, 1e-8)
lb_eta_k3_10  <- Box.test(eta_k3,   lag = 10, type = "Ljung-Box")
lb_eta_k3_20  <- Box.test(eta_k3,   lag = 20, type = "Ljung-Box")
lb_eta2_k3_10 <- Box.test(eta_k3^2, lag = 10, type = "Ljung-Box")
lb_eta2_k3_20 <- Box.test(eta_k3^2, lag = 20, type = "Ljung-Box")

cat(sprintf("K=3 LB eta   lag 10: Q=%7.2f, p=%.4f\n", lb_eta_k3_10$statistic,  lb_eta_k3_10$p.value))
cat(sprintf("K=3 LB eta   lag 20: Q=%7.2f, p=%.4f\n", lb_eta_k3_20$statistic,  lb_eta_k3_20$p.value))
cat(sprintf("K=3 LB eta²  lag 10: Q=%7.2f, p=%.4f\n", lb_eta2_k3_10$statistic, lb_eta2_k3_10$p.value))
cat(sprintf("K=3 LB eta²  lag 20: Q=%7.2f, p=%.4f\n", lb_eta2_k3_20$statistic, lb_eta2_k3_20$p.value))

phase3c_results <- list(
  k3_grid_table       = bench_k3_table,
  best_k3_name        = best_k3_name,
  best_k3_fit         = best_k3_fit,
  phase3c_pass        = TRUE,
  bic_improvement_over_k2 = bench_k2_table$BIC[1] - bench_k3_table$BIC[1],
  bootstrap_n_reps    = 200,
  bootstrap_block_size = 30,
  bootstrap_n_valid   = 200,
  robust_inference         = robust_inference_k3,
  robust_inference_aligned = robust_inference_k3_aligned,
  k3_resid_tests      = list(
    lb_eta10  = lb_eta_k3_10,
    lb_eta20  = lb_eta_k3_20,
    lb_eta210 = lb_eta2_k3_10,
    lb_eta220 = lb_eta2_k3_20
  )
)

# ── SAVE PHASE OBJECTS TO DISK ───────────────────────────────────────────────
saveRDS(phase3a_results, file.path(output_dir, "phase3a_results.rds"))
saveRDS(phase3b_results, file.path(output_dir, "phase3b_results.rds"))
saveRDS(phase3c_results, file.path(output_dir, "phase3c_results.rds"))
cat("\nAll three phase results saved to disk.\n")

# ── 8a. PHASE 4a: COMPUTE IN-SAMPLE VaR VECTORS ──────────────────────────────
# Compute mixture-quantile VaR at alpha = 0.05 and 0.01.

cat("\n=== 8a. Phase 4a — In-sample VaR computation ===\n")

t0_var <- proc.time()["elapsed"]

cat("\n  Computing VaR for K=1 gjrGARCH ...\n")
var_k1_05 <- compute_var_vec(best_k1_fit, 0.05)
cat("    alpha=0.05 done\n")
var_k1_01 <- compute_var_vec(best_k1_fit, 0.01)
cat("    alpha=0.01 done\n")

t1_var <- proc.time()["elapsed"]
cat(sprintf("\n  VaR computation complete: %.1f seconds\n", t1_var - t0_var))

# ── Summary: hit-rate vs expected rate ───────────────────────────────────────
cat("\n--- VaR summary ---\n")
cat(sprintf("  K=1, alpha=0.05: range [%.4f, %.4f]   hit-rate: %.4f (expected 0.0500)\n",
            min(var_k1_05), max(var_k1_05), mean(log_ret < var_k1_05)))
cat(sprintf("  K=1, alpha=0.01: range [%.4f, %.4f]   hit-rate: %.4f (expected 0.0100)\n",
            min(var_k1_01), max(var_k1_01), mean(log_ret < var_k1_01)))

# ── 8b. PHASE 4a: BACKTEST K=1 gjrGARCH ──────────────────────────────────────
# Christoffersen Conditional Coverage and Engle-Manganelli Dynamic Quantile
# tests for K=1 gjrGARCH at alpha = 0.05 and 0.01.

cat("\n=== 8b. Phase 4a — Backtest K=1 gjrGARCH ===\n")

# --- alpha = 0.05 ---
cat("\n--- alpha = 0.05 ---\n")
cc_k1_05 <- christoffersen_cc(log_ret, var_k1_05, alpha = 0.05)
dq_k1_05 <- dynamic_quantile_test(log_ret, var_k1_05, alpha = 0.05, lags = 5)

cat(sprintf("  Observed hit rate    : %.4f  (expected 0.0500)\n", cc_k1_05$hit_rate))
cat(sprintf("  Total violations     : %d / %d\n", cc_k1_05$n_violations, cc_k1_05$n_obs))
cat(sprintf("  LR_UC  (df=1) = %7.3f   p = %.4f   %s\n",
            cc_k1_05$LR_UC, cc_k1_05$p_UC,
            ifelse(cc_k1_05$p_UC  < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  LR_IND (df=1) = %7.3f   p = %.4f   %s\n",
            cc_k1_05$LR_IND, cc_k1_05$p_IND,
            ifelse(cc_k1_05$p_IND < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  LR_CC  (df=2) = %7.3f   p = %.4f   %s\n",
            cc_k1_05$LR_CC, cc_k1_05$p_CC,
            ifelse(cc_k1_05$p_CC  < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  DQ stat (df=%d) = %7.3f   p = %.4f   %s\n",
            dq_k1_05$df, dq_k1_05$DQ_stat, dq_k1_05$p_value,
            ifelse(dq_k1_05$p_value < 0.05, "[REJECT H0]", "[fail to reject]")))

# --- alpha = 0.01 ---
cat("\n--- alpha = 0.01 ---\n")
cc_k1_01 <- christoffersen_cc(log_ret, var_k1_01, alpha = 0.01)
dq_k1_01 <- dynamic_quantile_test(log_ret, var_k1_01, alpha = 0.01, lags = 5)

cat(sprintf("  Observed hit rate    : %.4f  (expected 0.0100)\n", cc_k1_01$hit_rate))
cat(sprintf("  Total violations     : %d / %d\n", cc_k1_01$n_violations, cc_k1_01$n_obs))
cat(sprintf("  LR_UC  (df=1) = %7.3f   p = %.4f   %s\n",
            cc_k1_01$LR_UC, cc_k1_01$p_UC,
            ifelse(cc_k1_01$p_UC  < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  LR_IND (df=1) = %7.3f   p = %.4f   %s\n",
            cc_k1_01$LR_IND, cc_k1_01$p_IND,
            ifelse(cc_k1_01$p_IND < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  LR_CC  (df=2) = %7.3f   p = %.4f   %s\n",
            cc_k1_01$LR_CC, cc_k1_01$p_CC,
            ifelse(cc_k1_01$p_CC  < 0.05, "[REJECT H0]", "[fail to reject]")))
cat(sprintf("  DQ stat (df=%d) = %7.3f   p = %.4f   %s\n",
            dq_k1_01$df, dq_k1_01$DQ_stat, dq_k1_01$p_value,
            ifelse(dq_k1_01$p_value < 0.05, "[REJECT H0]", "[fail to reject]")))

# ── 8c. PHASE 4a: PLOT K=1 VaR WITH VIOLATIONS ──
pdf(file.path(output_dir, "var_k1_backtest.pdf"), width = 11, height = 5)

# Violations
vio_05 <- which(log_ret < var_k1_05)
vio_01 <- which(log_ret < var_k1_01)

plot(dates, log_ret, type = "l", col = "grey50", lwd = 0.5,
     main = sprintf("K=1 gjrGARCH In-Sample VaR  |  α=0.05 hit=%.3f  α=0.01 hit=%.3f",
                    mean(log_ret < var_k1_05), mean(log_ret < var_k1_01)),
     xlab = "Date", ylab = "Log-Return (OLS-residual)")
lines(dates, var_k1_05, col = "firebrick",  lwd = 1.2, lty = 2)
lines(dates, var_k1_01, col = "darkorange", lwd = 1.2, lty = 3)
points(dates[vio_05], log_ret[vio_05], col = "firebrick",  pch = 19, cex = 0.3)
points(dates[vio_01], log_ret[vio_01], col = "darkorange", pch = 17, cex = 0.5)
abline(h = 0, col = "grey80", lty = 1)
legend("topright",
       legend = c("Return", "VaR α=0.05", "VaR α=0.01", "Viol α=0.05", "Viol α=0.01"),
       col    = c("grey50", "firebrick", "darkorange", "firebrick", "darkorange"),
       lty    = c(1, 2, 3, NA, NA),
       pch    = c(NA, NA, NA, 19, 17),
       cex    = 0.75, bg = "white")
dev.off()

cat("\nPlot saved: var_k1_backtest.pdf\n")

# Gem fase 4a-resultater
phase4a_results <- list(
  var_k1_05    = var_k1_05,
  var_k1_01    = var_k1_01,
  cc_k1_05     = cc_k1_05,
  cc_k1_01     = cc_k1_01,
  dq_k1_05     = dq_k1_05,
  dq_k1_01     = dq_k1_01,
  primary_model = "K=1 gjrGARCH (Student-t)"
)

saveRDS(phase4a_results, file.path(output_dir, "phase4a_results.rds"))
cat("phase4a_results saved.\n")

# ── 8d. PHASE 4a: NUMERICAL DIAGNOSIS OF K=1 VIOLATIONS ──────────────────────
# Distribution of violations across years and worst-return diagnostic.

cat("\n=== 8d. Phase 4a — Numerical diagnosis of K=1 violations ===\n")

vio_05_dates <- dates[which(log_ret < var_k1_05)]
vio_01_dates <- dates[which(log_ret < var_k1_01)]

# Violations per year
vio_05_per_year <- table(format(vio_05_dates, "%Y"))
vio_01_per_year <- table(format(vio_01_dates, "%Y"))

cat("\n--- Violations per year (alpha = 0.05) ---\n")
print(vio_05_per_year)

cat("\n--- Violations per year (alpha = 0.01) ---\n")
print(vio_01_per_year)

# 15 most negative returns and corresponding VaR levels
worst_idx <- order(log_ret)[1:15]
worst_returns_diagnosis <- data.frame(
  Date   = as.character(dates[worst_idx]),
  Return = round(log_ret[worst_idx], 4),
  VaR_05 = round(var_k1_05[worst_idx], 4),
  VaR_01 = round(var_k1_01[worst_idx], 4),
  Hit_05 = log_ret[worst_idx] < var_k1_05[worst_idx],
  Hit_01 = log_ret[worst_idx] < var_k1_01[worst_idx],
  stringsAsFactors = FALSE
)

cat("\n--- 15 most negative log-returns ---\n")
print(worst_returns_diagnosis, row.names = FALSE)

# Update phase4a_results with diagnosis
phase4a_results$vio_05_per_year         <- vio_05_per_year
phase4a_results$vio_01_per_year         <- vio_01_per_year
phase4a_results$worst_returns_diagnosis <- worst_returns_diagnosis

# Re-save to disk
saveRDS(phase4a_results, file.path(output_dir, "phase4a_results.rds"))
cat("\nUpdated phase4a_results saved to disk.\n")

# ── 9. PHASE 4b: PERIODIC THREE-WAY OOS BACKTEST (Risk nahead=1) ──────────────

cat("\n=== 9. Phase 4b — Periodic three-way OOS (Risk nahead=1) ===\n")
library(parallel)

W_REFIT <- 90L                                   # <<< re-estimation every 90 days
T_obs   <- length(log_ret)
burn_in <- T_obs - 1095
n_oos   <- T_obs - burn_in
ret_oos   <- log_ret[(burn_in + 1):T_obs]
dates_oos <- dates[(burn_in + 1):T_obs]

run_oos_periodic <- function(model_vec, n_starts, W, seed_base, label) {
  K        <- length(model_vec)
  refit_at <- seq(burn_in, T_obs - 1L, by = W)
  cat(sprintf("  %-14s : %d refits (W=%d) + %d billige VaR-kald\n",
              label, length(refit_at), W, n_oos)); flush.console()
  t0 <- proc.time()["elapsed"]
  cl <- makeCluster(N_CORES); on.exit(stopCluster(cl))
  clusterEvalQ(cl, library(MSGARCH))
  clusterExport(cl, c("log_ret", "model_vec", "n_starts", "seed_base", "K"),
                envir = environment())
  
  ## Phase 1 — estimate parameters at each refit origin
  par_list <- parLapplyLB(cl, refit_at, function(tau) {
    sp <- MSGARCH::CreateSpec(variance.spec     = list(model = model_vec),
                              distribution.spec = list(distribution = rep("std", K)))
    best_ll <- -Inf; best <- NULL
    for (s in seq_len(n_starts)) {
      set.seed(seed_base + tau * 7L + s)
      f <- tryCatch(MSGARCH::FitML(sp, data = log_ret[1:tau]), error = function(e) NULL)
      if (!is.null(f) && is.finite(f$loglik) && f$loglik > best_ll) {
        best_ll <- f$loglik; best <- f$par
      }
    }
    best
  })
  names(par_list) <- as.character(refit_at)
  n_fitfail <- sum(vapply(par_list, is.null, logical(1)))
  
  ## Phase 2 — one-step VaR per out-of-sample day
  clusterExport(cl, c("par_list", "refit_at"), envir = environment())
  v_list <- parLapplyLB(cl, (burn_in + 1L):T_obs, function(t) {
    tau <- max(refit_at[refit_at <= t - 1L])
    p   <- par_list[[as.character(tau)]]
    if (is.null(p)) return(c(NA_real_, NA_real_))
    sp <- MSGARCH::CreateSpec(variance.spec     = list(model = model_vec),
                              distribution.spec = list(distribution = rep("std", K)))
    tryCatch(as.numeric(MSGARCH::Risk(sp, par = p, data = log_ret[1:(t - 1L)],
                                      alpha = c(0.05, 0.01), nahead = 1)$VaR),
             error = function(e) c(NA_real_, NA_real_))     # eGARCH could fail here 
  })
  mat <- do.call(rbind, v_list); colnames(mat) <- c("VaR05", "VaR01")
  cat(sprintf("    %.1f min  |  fit-fejl: %d/%d refits  |  VaR-skips: %d/%d\n",
              (proc.time()["elapsed"] - t0)/60, n_fitfail, length(refit_at),
              sum(is.na(mat[, 1])), n_oos))
  mat
}

eval_track <- function(mat, label) {
  cat(sprintf("\n--- %s ---\n", label)); out <- list()
  for (j in 1:2) {
    a <- c(0.05, 0.01)[j]; v <- mat[, j]; ok <- !is.na(v)
    cc <- christoffersen_cc(ret_oos[ok], v[ok], a)
    dq <- tryCatch(dynamic_quantile_test(ret_oos[ok], v[ok], a, 5), error = function(e) NULL)
    cat(sprintf("  a=%.2f  hit=%.4f  (n=%d, viol=%d)  p_UC=%.4f  p_CC=%.4f%s\n",
                a, cc$hit_rate, sum(ok), cc$n_violations, cc$p_UC, cc$p_CC,
                if (!is.null(dq)) sprintf("  p_DQ=%.4f", dq$p_value) else ""))
    out[[as.character(a)]] <- list(cc = cc, dq = dq)
  }
  invisible(out)
}

get_or_run <- function(path, fun) { 
  if (file.exists(path)) { cat("  cache:", basename(path), "\n"); readRDS(path) }
  else { x <- fun(); saveRDS(x, path); x }
}

# ── Validating the cheap spec engine against the fit engine at one origin
tau_chk <- T_obs - 1095
sp_chk  <- CreateSpec(variance.spec = list(model = "gjrGARCH"),
                      distribution.spec = list(distribution = "std"))
fit_chk <- FitML(sp_chk, data = log_ret[1:tau_chk])
v_fit   <- as.numeric(Risk(fit_chk, alpha = c(0.05, 0.01), nahead = 1)$VaR)
v_spec  <- as.numeric(Risk(sp_chk, par = fit_chk$par, data = log_ret[1:tau_chk],
                           alpha = c(0.05, 0.01), nahead = 1)$VaR)
cat("\nSpec-vs-fit validering (skal være identiske):\n"); print(rbind(fit = v_fit, spec = v_spec))
stopifnot(max(abs(v_fit - v_spec)) < 1e-6)

# ── K = 1 first
oos_k1 <- get_or_run(file.path(output_dir, "oos_k1_gjr_risk.rds"),
                     function() run_oos_periodic("gjrGARCH", 1L, W_REFIT, 40000L, "K=1 gjr"))
ev_k1 <- eval_track(oos_k1, "K=1 gjrGARCH (Risk nahead=1)")

# Integrate the K=1
phase4b_results <- list(
  T_obs = T_obs, burn_in = burn_in, n_oos = n_oos, n_skipped = sum(is.na(oos_k1[, 1])),
  oos_dates = dates_oos, oos_returns = ret_oos,
  VaR_oos_k1_05 = oos_k1[, 1], VaR_oos_k1_01 = oos_k1[, 2],
  cc_oos_05 = ev_k1[["0.05"]]$cc, cc_oos_01 = ev_k1[["0.01"]]$cc,
  dq_oos_05 = ev_k1[["0.05"]]$dq, dq_oos_01 = ev_k1[["0.01"]]$dq,
  primary_model = "K=1 gjrGARCH (Student-t)",
  oos_design = sprintf("Expanding window, T-1095 burn-in, Risk(nahead=1), refit hver %d. dag", W_REFIT))
saveRDS(phase4b_results, file.path(output_dir, "phase4b_results.rds"))

# ── The two expensive K = 3 tracks
DO_K3 <- TRUE
if (DO_K3) {
  oos_k3_gjr <- get_or_run(file.path(output_dir, "oos_k3_gjr_risk.rds"),
                           function() run_oos_periodic(rep("gjrGARCH", 3), 3L, W_REFIT, 50000L, "K=3 gjr"))
  ev_k3_gjr <- eval_track(oos_k3_gjr, "K=3 gjrGARCH"); print(summary(oos_k3_gjr[, 1]))
  
  oos_k3_eg <- get_or_run(file.path(output_dir, "oos_k3_egarch_risk.rds"),
                          function() run_oos_periodic(rep("eGARCH", 3), 3L, W_REFIT, 60000L, "K=3 eGARCH"))
  ev_k3_eg <- eval_track(oos_k3_eg, "K=3 eGARCH"); print(summary(oos_k3_eg[, 1]))
}

# ── 10. GENERATE LATEX TABLES FOR PAPER ──────────────────────────────────────
# Generate all .tex tables in one batch.

cat("\n=== 10. Generating LaTeX tables for paper ===\n")

library(xtable)

# ── Formatting helpers ───────────────────────────────────────────────────────

# Format a numeric value with fixed number of decimals, preserving trailing zeros.
fmt <- function(x, d = 4) {
  formatC(x, format = "f", digits = d)
}

# Format p-value
fmt_p <- function(p, threshold = 0.001) {
  out <- character(length(p))
  for (i in seq_along(p)) {
    if (is.na(p[i])) {
      out[i] <- "---"
    } else if (p[i] < 0.0001) {
      out[i] <- "$<\\!0.0001$"
    } else if (p[i] < threshold) {
      out[i] <- "$<\\!0.001$"
    } else {
      out[i] <- formatC(p[i], format = "f", digits = 4)
    }
  }
  out
}

# Map raw MSGARCH parameter names to Greek LaTeX labels.
greek_label <- function(par_name) {
  if (grepl("^alpha0_(\\d+)$", par_name))
    return(sprintf("$\\omega_{%s}$", sub("alpha0_", "", par_name)))
  if (grepl("^alpha1_(\\d+)$", par_name))
    return(sprintf("$\\alpha_{%s}$", sub("alpha1_", "", par_name)))
  if (grepl("^alpha2_(\\d+)$", par_name))
    return(sprintf("$\\gamma_{%s}$", sub("alpha2_", "", par_name)))
  if (grepl("^beta_(\\d+)$", par_name))
    return(sprintf("$\\beta_{%s}$", sub("beta_", "", par_name)))
  if (grepl("^nu_(\\d+)$", par_name))
    return(sprintf("$\\nu_{%s}$", sub("nu_", "", par_name)))
  if (grepl("^P_(\\d+)_(\\d+)$", par_name)) {
    parts <- strsplit(par_name, "_")[[1]]
    return(sprintf("$P_{%s%s}$", parts[2], parts[3]))
  }
  par_name
}

# Write xtable to .tex file with consistent settings.
xt_write <- function(xt, fname, ...) {
  path <- file.path(output_dir, fname)
  print(xt, file = path, booktabs = TRUE, floating = FALSE,
        sanitize.text.function = function(x) x,
        include.rownames = FALSE, ...)
  cat(sprintf("  Saved: %s\n", fname))
}

# ── (1) Descriptive statistics ────
desc_df <- data.frame(
  Series   = rownames(desc_stats_tbl),
  n        = format(desc_stats_tbl[, "n"], big.mark = ","),
  mean     = fmt(desc_stats_tbl[, "mean"], 4),
  median   = fmt(desc_stats_tbl[, "median"], 4),
  sd       = fmt(desc_stats_tbl[, "sd"], 4),
  skewness = fmt(desc_stats_tbl[, "skewness"], 3),
  kurtosis = fmt(desc_stats_tbl[, "kurtosis"], 3),
  min      = fmt(desc_stats_tbl[, "min"], 2),
  max      = fmt(desc_stats_tbl[, "max"], 2),
  stringsAsFactors = FALSE
)
names(desc_df) <- c("Series", "$n$", "Mean", "Median", "SD",
                    "Skew.", "Ex. kurt.", "Min", "Max")

xt_write(
  xtable(desc_df,
         caption = paste0("Descriptive statistics for daily mean electricity ",
                          "spot prices (DKK/MWh) and log-returns. Excess ",
                          "kurtosis reported (Gaussian benchmark $= 0$). ",
                          "Log-returns lose one observation to differencing."),
         label   = "tab:desc_stats",
         align   = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r")),
  "tab_desc_stats.tex"
)

# ── (2) Pre-estimation tests ─────────────────────────────────────────────────
pretests_stats <- c(
  desc_tests$jb_ret$statistic,   desc_tests$jb_price$statistic,
  desc_tests$lb_r10$statistic,   desc_tests$lb_r20$statistic,
  desc_tests$lb_sq10$statistic,  desc_tests$lb_sq20$statistic,
  desc_tests$arch_lm$statistic,
  desc_tests$adf_lp$statistic,   desc_tests$adf_lr$statistic)

pretests_pvals <- c(
  desc_tests$jb_ret$p.value,   desc_tests$jb_price$p.value,
  desc_tests$lb_r10$p.value,   desc_tests$lb_r20$p.value,
  desc_tests$lb_sq10$p.value,  desc_tests$lb_sq20$p.value,
  desc_tests$arch_lm$p.value,
  desc_tests$adf_lp$p.value,   desc_tests$adf_lr$p.value)

# ── Annotate values capped at 0.01 by tseries package ──#
pretests_pvals_fmt <- fmt_p(pretests_pvals)
pretests_pvals_fmt[8:9] <- "$\\leq\\!0.01$"  # ADF rows

pretests_df <- data.frame(
  Test = c("Jarque-Bera (log-returns)",
           "Jarque-Bera (log-prices)",
           "Ljung-Box $r_t$, lag 10",
           "Ljung-Box $r_t$, lag 20",
           "Ljung-Box $r_t^2$, lag 10",
           "Ljung-Box $r_t^2$, lag 20",
           "ARCH-LM, lag 10",
           "ADF (log-prices)",
           "ADF (log-returns)"),
  Statistic = fmt(pretests_stats, 2),
  p_value   = pretests_pvals_fmt,
  stringsAsFactors = FALSE
)
names(pretests_df) <- c("Test", "Statistic", "$p$-value")

xt_write(
  xtable(pretests_df,
         caption = paste0("Pre-estimation diagnostic tests on raw log-returns ",
                          "and log-prices. ADF $p$-values are bounded above ",
                          "at $0.01$ by the interpolation table in ",
                          "\\texttt{tseries::adf.test}."),
         label   = "tab:pretests",
         align   = c("l", "l", "r", "r")),
  "tab_pretests.tex"
)

# ── (3) Mean-model OLS coefficients ──────────────────────────────────────────
ols_coef <- summary(ols_filter)$coefficients
ols_df <- data.frame(
  Regressor = c("Intercept", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun",
                "Holiday", paste0("AR(", 1:7, ")")),
  Estimate  = fmt(ols_coef[, "Estimate"], 4),
  SE        = fmt(ols_coef[, "Std. Error"], 4),
  t_value   = fmt(ols_coef[, "t value"], 2),
  p_value   = fmt_p(ols_coef[, "Pr(>|t|)"]),
  stringsAsFactors = FALSE
)
names(ols_df) <- c("Regressor", "Estimate", "Std. Error", "$t$", "$p$-value")

xt_write(
  xtable(ols_df,
         caption = paste0("OLS pre-filter coefficients. Monday is the ",
                          "reference weekday; weekday dummies and the holiday ",
                          "dummy capture the deterministic weekly demand ",
                          "cycle, while AR(1)--AR(7) absorbs short-lag ",
                          "mean-reversion."),
         label   = "tab:mean_model",
         align   = c("l", "l", "r", "r", "r", "r")),
  "tab_mean_model.tex"
)

# ── (4) Ljung-Box filter steps ───────────────────────────────────────────────
lb_steps_df <- data.frame(
  Stage = c("(a) Raw (trimmed) series",
            "(b) Weekday + holiday dummies only",
            "(c) Dummies + AR(1)--AR(7) residuals"),
  Q_lag20 = fmt(c(lb_raw$statistic, lb_dummies$statistic, lb_full$statistic), 2),
  p_value = fmt_p(c(lb_raw$p.value, lb_dummies$p.value, lb_full$p.value)),
  stringsAsFactors = FALSE
)
names(lb_steps_df) <- c("Filtering stage", "$Q(20)$", "$p$-value")

xt_write(
  xtable(lb_steps_df,
         caption = paste0("Marginal effect of each pre-filtering step on the ",
                          "Ljung-Box statistic at lag $20$. The full model ",
                          "removes $83\\%$ of the mean autocorrelation."),
         label   = "tab:lb_filter_steps",
         align   = c("l", "l", "r", "r")),
  "tab_lb_filter_steps.tex"
)
# ── (4b) Consolidated model comparison across all nine specifications ────────.

compare_all <- rbind(
  data.frame(
    Model = paste0(bench_table$GARCH, "_std_K1"),
    K = 1, GARCH = bench_table$GARCH, Dist = "std",
    LogLik = bench_table$LogLik, nPar = bench_table$nPar,
    AIC = bench_table$AIC, BIC = bench_table$BIC,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Model = paste0(bench_k2_table$GARCH, "_std_K2"),
    K = 2, GARCH = bench_k2_table$GARCH, Dist = "std",
    LogLik = bench_k2_table$LogLik, nPar = bench_k2_table$nPar,
    AIC = bench_k2_table$AIC, BIC = bench_k2_table$BIC,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Model = paste0(bench_k3_table$GARCH, "_std_K3"),
    K = 3, GARCH = bench_k3_table$GARCH, Dist = "std",
    LogLik = bench_k3_table$LogLik, nPar = bench_k3_table$nPar,
    AIC = bench_k3_table$AIC, BIC = bench_k3_table$BIC,
    stringsAsFactors = FALSE
  )
)
compare_all <- compare_all[order(compare_all$BIC), ]

# Replace underscores in model names
compare_df <- data.frame(
  Model  = gsub("_", "\\\\_", compare_all$Model),
  K      = compare_all$K,
  GARCH  = compare_all$GARCH,
  Dist   = compare_all$Dist,
  LogLik = fmt(compare_all$LogLik, 2),
  nPar   = compare_all$nPar,
  AIC    = fmt(compare_all$AIC, 2),
  BIC    = fmt(compare_all$BIC, 2),
  stringsAsFactors = FALSE
)
names(compare_df) <- c("Model", "$K$", "GARCH", "Dist.",
                       "Log-likelihood", "nPar", "AIC", "BIC")

xt_write(
  xtable(compare_df,
         caption = paste0("Model comparison across all nine candidate ",
                          "specifications: three GARCH families ",
                          "(sGARCH, gjrGARCH, eGARCH) crossed with three ",
                          "regime counts ($K \\in \\{1, 2, 3\\}$), all with ",
                          "Student-$t$ innovations. Models are ranked by BIC ",
                          "(lower is better). The $K = 3$ eGARCH ",
                          "specification dominates on both AIC and BIC, ",
                          "and within each $K$, eGARCH outperforms gjrGARCH ",
                          "and sGARCH on log-likelihood."),
         label   = "tab:model_comparison",
         align   = c("l", "l", "r", "l", "l", "r", "r", "r", "r")),
  "tab_model_comparison.tex"
)

# ── (5) K=1 gjrGARCH parameter estimates ─────────────────────────────────────
k1_coef <- phase3a_results$k1_coef
k1_df <- data.frame(
  Parameter = sapply(rownames(k1_coef), greek_label),
  Estimate  = fmt(k1_coef[, "Estimate"], 4),
  SE        = fmt(k1_coef[, "Std. Error"], 4),
  t_value   = fmt(k1_coef[, "t value"], 2),
  p_value   = fmt_p(k1_coef[, "Pr(>|t|)"]),
  stringsAsFactors = FALSE
)
names(k1_df) <- c("Parameter", "Estimate", "Std. Error", "$t$", "$p$-value")

xt_write(
  xtable(k1_df,
         caption = paste0("Parameter estimates for the single-regime ",
                          "gjrGARCH$(1,1)$ baseline model with Student-$t$ ",
                          "errors. Hessian-based standard errors. ",
                          "Log-likelihood: ",
                          fmt(phase3a_results$k1_loglik, 2),
                          "; BIC: ",
                          fmt(phase3a_results$k1_BIC, 2), "."),
         label   = "tab:k1_params",
         align   = c("l", "l", "r", "r", "r", "r")),
  "tab_k1_params.tex"
)

# ── (6) K=2 eGARCH parameter estimates (ALIGNED bootstrap inference) ─────────
k2_inf <- phase3b_results$robust_inference_aligned
k2_df <- data.frame(
  Parameter = sapply(rownames(k2_inf), greek_label),
  Estimate  = fmt(k2_inf$Estimate, 4),
  Boot_SE   = fmt(k2_inf$Boot_SD, 4),
  t_value   = fmt(k2_inf$t_MAD, 2),
  CI_low    = fmt(k2_inf$CI_low_95, 4),
  CI_high   = fmt(k2_inf$CI_high_95, 4),
  stringsAsFactors = FALSE
)
names(k2_df) <- c("Parameter", "Estimate", "Boot. SE", "$t$",
                  "CI low $(2.5\\%)$", "CI high $(97.5\\%)$")

k2_bic <- phase3b_results$k2_grid_table$BIC[1]
xt_write(
  xtable(k2_df,
         caption = paste0("Parameter estimates for the two-regime ",
                          "eGARCH$(1,1)$ specification with Student-$t$ ",
                          "errors. Standard errors (MAD-based) and percentile ",
                          "confidence intervals from stationary block bootstrap ",
                          "($B=200$, block length $L=30$), after ",
                          "label-switching alignment. ",
                          "Log-likelihood: ",
                          fmt(phase3b_results$best_k2_fit$loglik, 2),
                          "; BIC: ",
                          fmt(k2_bic, 2),
                          " ($\\Delta\\text{BIC}$ vs $K=1$: ",
                          fmt(phase3b_results$bic_improvement_over_k1, 2),
                          ")."),
         label   = "tab:k2_params",
         align   = c("l", "l", "r", "r", "r", "r", "r")),
  "tab_k2_params.tex"
)

# ── (7) K=3 eGARCH parameter estimates (ALIGNED bootstrap inference) ─────────
k3_inf <- phase3c_results$robust_inference_aligned
k3_df <- data.frame(
  Parameter = sapply(rownames(k3_inf), greek_label),
  Estimate  = fmt(k3_inf$Estimate, 4),
  Boot_SE   = fmt(k3_inf$Boot_SD, 4),
  t_value   = fmt(k3_inf$t_MAD, 2),
  CI_low    = fmt(k3_inf$CI_low_95, 4),
  CI_high   = fmt(k3_inf$CI_high_95, 4),
  stringsAsFactors = FALSE
)
names(k3_df) <- c("Parameter", "Estimate", "Boot. SE", "$t$",
                  "CI low $(2.5\\%)$", "CI high $(97.5\\%)$")

k3_bic <- phase3c_results$k3_grid_table$BIC[1]
xt_write(
  xtable(k3_df,
         caption = paste0("Parameter estimates for the three-regime ",
                          "eGARCH$(1,1)$ specification with Student-$t$ ",
                          "errors. Standard errors (MAD-based) and percentile ",
                          "confidence intervals from stationary block bootstrap ",
                          "($B=200$, block length $L=30$), after ",
                          "label-switching alignment. ",
                          "Log-likelihood: ",
                          fmt(phase3c_results$best_k3_fit$loglik, 2),
                          "; BIC: ",
                          fmt(k3_bic, 2),
                          " ($\\Delta\\text{BIC}$ vs $K=2$: ",
                          fmt(phase3c_results$bic_improvement_over_k2, 2),
                          ")."),
         label   = "tab:k3_params",
         align   = c("l", "l", "r", "r", "r", "r", "r")),
  "tab_k3_params.tex"
)

# ── (8) K=3 transition matrix + ergodic distribution ─────────────────────────
P_k3 <- build_trans_mat(phase3c_results$best_k3_fit$par, 3)
erg_k3 <- ergodic_dist(P_k3)
dur_k3 <- 1 / (1 - diag(P_k3))

P_df <- data.frame(
  From = paste0("$k=", 1:3, "$"),
  to_1 = fmt(P_k3[, 1], 4),
  to_2 = fmt(P_k3[, 2], 4),
  to_3 = fmt(P_k3[, 3], 4),
  Duration = fmt(dur_k3, 1),
  Ergodic  = fmt(erg_k3, 4),
  stringsAsFactors = FALSE
)
names(P_df) <- c("From regime", "To $k=1$", "To $k=2$", "To $k=3$",
                 "Exp. duration (days)", "Ergodic prob.")

xt_write(
  xtable(P_df,
         caption = paste0("Estimated transition probability matrix $\\hat{P}$, ",
                          "expected regime durations ($1/(1-\\hat{P}_{kk})$), ",
                          "and ergodic probabilities for the $K=3$ eGARCH ",
                          "model. Regime $3$ is the most persistent ",
                          "(crisis-like) state but not absorbing ",
                          "($\\hat{P}_{33} = ",
                          fmt(P_k3[3, 3], 4), "$)."),
         label   = "tab:k3_transition",
         align   = c("l", "l", "r", "r", "r", "r", "r")),
  "tab_k3_transition.tex"
)

# ── Backtest helper: format a single backtest row with significance marks ────
fmt_backtest_row <- function(cc, dq, alpha) {
  # Star marker for hits where p < 0.05
  star <- function(p) ifelse(is.na(p) | p >= 0.05, "", "$^{*}$")
  
  data.frame(
    Alpha      = sprintf("$%.2f$", alpha),
    Hit_rate   = fmt(cc$hit_rate, 4),
    Violations = format(cc$n_violations, big.mark = ","),
    LR_UC      = paste0(fmt(cc$LR_UC, 2), star(cc$p_UC)),
    p_UC       = fmt_p(cc$p_UC),
    LR_IND     = paste0(fmt(cc$LR_IND, 2), star(cc$p_IND)),
    p_IND      = fmt_p(cc$p_IND),
    LR_CC      = paste0(fmt(cc$LR_CC, 2), star(cc$p_CC)),
    p_CC       = fmt_p(cc$p_CC),
    DQ         = paste0(fmt(dq$DQ_stat, 2), star(dq$p_value)),
    p_DQ       = fmt_p(dq$p_value),
    stringsAsFactors = FALSE
  )
}

# ── (9) In-sample VaR backtest (K=1 gjrGARCH) ────────────────────────────────
bt_is_df <- rbind(
  fmt_backtest_row(phase4a_results$cc_k1_05, phase4a_results$dq_k1_05, 0.05),
  fmt_backtest_row(phase4a_results$cc_k1_01, phase4a_results$dq_k1_01, 0.01)
)
names(bt_is_df) <- c("$\\alpha$", "Hit rate", "Viol.",
                     "$LR_{UC}$", "$p_{UC}$",
                     "$LR_{IND}$", "$p_{IND}$",
                     "$LR_{CC}$", "$p_{CC}$",
                     "$DQ$", "$p_{DQ}$")

xt_write(
  xtable(bt_is_df,
         caption = paste0("In-sample VaR backtest for $K=1$ gjrGARCH ",
                          "($T = ", format(length(log_ret), big.mark = ","),
                          "$ observations). Christoffersen (1998) Conditional ",
                          "Coverage tests ($LR_{UC}$, $LR_{IND}$, $LR_{CC}$) ",
                          "and Engle--Manganelli (2004) Dynamic Quantile ",
                          "test with $5$ lags. ",
                          "Asterisks ($^{*}$) mark statistics with $p<0.05$."),
         label   = "tab:backtest_insample",
         align   = c("l", "l", rep("r", ncol(bt_is_df) - 1))),
  "tab_backtest_insample.tex"
)

# ── (10) Out-of-sample rolling backtest (K=1 gjrGARCH) ───────────────────────
bt_oos_df <- rbind(
  fmt_backtest_row(phase4b_results$cc_oos_05, phase4b_results$dq_oos_05, 0.05),
  fmt_backtest_row(phase4b_results$cc_oos_01, phase4b_results$dq_oos_01, 0.01)
)
names(bt_oos_df) <- c("$\\alpha$", "Hit rate", "Viol.",
                      "$LR_{UC}$", "$p_{UC}$",
                      "$LR_{IND}$", "$p_{IND}$",
                      "$LR_{CC}$", "$p_{CC}$",
                      "$DQ$", "$p_{DQ}$")

xt_write(
  xtable(bt_oos_df,
         caption = paste0("Out-of-sample VaR backtest for $K=1$ gjrGARCH ",
                          "(expanding window, burn-in $T - 1095$, ",
                          format(phase4b_results$n_oos, big.mark = ","),
                          " OOS observations covering ",
                          format(min(phase4b_results$oos_dates), "%b %Y"),
                          "--",
                          format(max(phase4b_results$oos_dates), "%b %Y"),
                          "). Asterisks ($^{*}$) mark statistics with $p<0.05$."),
         label   = "tab:backtest_oos",
         align   = c("l", "l", rep("r", ncol(bt_oos_df) - 1))),
  "tab_backtest_oos.tex"
)

cat("\nAll LaTeX tables written to:", output_dir, "\n")

# ── (11) Three-way OOS comparison: K=1 gjr / gjr-K3 / eGARCH-K3 ──
oos_row <- function(ev, model, a) {
  cc <- ev[[as.character(a)]]$cc; dq <- ev[[as.character(a)]]$dq
  cbind(Model = model, fmt_backtest_row(cc, dq, a))
}
bt_3way_df <- rbind(
  oos_row(ev_k1,     "K=1 gjrGARCH", 0.05), oos_row(ev_k1,     "K=1 gjrGARCH", 0.01),
  oos_row(ev_k3_gjr, "K=3 gjrGARCH", 0.05), oos_row(ev_k3_gjr, "K=3 gjrGARCH", 0.01),
  oos_row(ev_k3_eg,  "K=3 eGARCH",   0.05), oos_row(ev_k3_eg,  "K=3 eGARCH",   0.01)
)
names(bt_3way_df) <- c("Model", "$\\alpha$", "Hit rate", "Viol.",
                       "$LR_{UC}$", "$p_{UC}$", "$LR_{IND}$", "$p_{IND}$",
                       "$LR_{CC}$", "$p_{CC}$", "$DQ$", "$p_{DQ}$")
xt_write(
  xtable(bt_3way_df,
         caption = paste0("Out-of-sample VaR backtest: single-regime baseline ",
                          "(K=1 gjrGARCH) vs.\\ two three-regime specifications, ",
                          "all Student-$t$, via one-step predictive-mixture VaR ",
                          "(\\texttt{Risk}, $h{=}1$) on an expanding window with ",
                          "quarterly re-estimation. Violation rates fall ",
                          "monotonically from K=1 to K=3 at both levels, with ",
                          "eGARCH-K3 closest to nominal, though all three still ",
                          "reject unconditional coverage. eGARCH-K3 evaluated over ",
                          "1,091 days (4 excluded for forecast failure). ",
                          "Asterisks mark $p<0.05$."),
         label   = "tab:backtest_oos_3way",
         align   = c("l", "l", "l", rep("r", 10))),
  "tab_backtest_oos_3way.tex"
)

# ── 11. PAPER PLOTS: REGIME STRUCTURE AND DIAGNOSTICS ────────────────────────
# Five plots:
#   (a) Smoothed regime probabilities for K=3 eGARCH (3 panels)
#   (b) Viterbi regime sequence with shaded background (K=3)
#   (c) Conditional volatility from K=1 gjrGARCH
#   (d) Standardised residual diagnostics (histogram + Q-Q) for K=1
#   (e) OOS VaR backtest plot (K=1 gjrGARCH)

cat("\n=== 11. Paper plots ===\n")

# ── (a) Smoothed regime probabilities (K=3) ──────────────────────────────────

states_k3   <- State(best_k3_fit)
T_k3        <- length(log_ret)
smooth_k3   <- states_k3$SmoothProb[seq_len(T_k3), 1, , drop = FALSE]
smooth_k3   <- smooth_k3[, 1, ]
if (is.null(dim(smooth_k3))) smooth_k3 <- matrix(smooth_k3, ncol = 3)

regime_labels <- c("Regime 1 (calm)", "Regime 2 (reactive)", "Regime 3 (crisis)")
regime_cols   <- c("steelblue", "darkorange", "firebrick")

pdf(file.path(output_dir, "regime_smoothed_probs_k3.pdf"),
    width = 10, height = 7)
par(mfrow = c(3, 1), mar = c(3, 4, 2.5, 1))
for (k in 1:3) {
  plot(dates, smooth_k3[, k], type = "l", col = regime_cols[k],
       ylim = c(0, 1), lwd = 0.7,
       main = sprintf("Smoothed Probability -- %s   (mean = %.3f)",
                      regime_labels[k], mean(smooth_k3[, k])),
       xlab = "", ylab = sprintf("P(s_t = %d | data)", k))
  abline(h = 0.5, col = "grey60", lty = 2)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
cat("  Saved: regime_smoothed_probs_k3.pdf\n")


# ── (b) Viterbi regime sequence with shaded background (K=3) ─────────────────

viterbi_k3 <- states_k3$Viterbi[seq_len(T_k3)]

pdf(file.path(output_dir, "regime_viterbi_k3.pdf"),
    width = 11, height = 4.5)
shade_cols <- c("lightblue", "navajowhite", "lightpink")

plot(dates, log_ret, type = "n",
     main = "Daily Log-Returns with Viterbi Regime Sequence (K=3 eGARCH)",
     xlab = "Date", ylab = "Log-Return (OLS-residual)")

# Shade background by contiguous regime runs
rle_v   <- rle(as.integer(viterbi_k3))
cum_end <- cumsum(rle_v$lengths)
cum_beg <- c(1, cum_end[-length(cum_end)] + 1)
for (seg in seq_along(rle_v$values)) {
  k_s <- rle_v$values[seg]
  x0  <- dates[cum_beg[seg]]
  x1  <- dates[cum_end[seg]]
  rect(x0, par("usr")[3], x1, par("usr")[4],
       col = shade_cols[k_s], border = NA)
}
lines(dates, log_ret, col = "grey25", lwd = 0.4)
abline(h = 0, col = "grey80", lty = 1)
legend("topright", legend = regime_labels,
       fill = shade_cols, cex = 0.8, bg = "white", inset = 0.02)
dev.off()
cat("  Saved: regime_viterbi_k3.pdf\n")


# ── (c) Conditional volatility from K=1 gjrGARCH ─────────────────────────────
vol_k1_plot <- Volatility(best_k1_fit)

pdf(file.path(output_dir, "conditional_volatility_k1.pdf"),
    width = 10, height = 4)
plot(dates, vol_k1_plot, type = "l", col = "darkgreen", lwd = 0.6,
     main = "Conditional Volatility -- K=1 gjrGARCH(1,1) with Student-t",
     xlab = "Date", ylab = expression(hat(sigma)[t]))
abline(h = median(vol_k1_plot), col = "grey60", lty = 2)
dev.off()
cat("  Saved: conditional_volatility_k1.pdf\n")


# ── (d) Standardised residual diagnostics from K=1 gjrGARCH ──────────────────

eta_k1_plot <- log_ret / pmax(vol_k1_plot, 1e-8)
nu_k1       <- as.numeric(best_k1_fit$par["nu_1"])

pdf(file.path(output_dir, "residuals_diag_k1.pdf"),
    width = 11, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Histogram
xl <- quantile(eta_k1_plot, c(0.001, 0.999))
hist(eta_k1_plot, breaks = 80, freq = FALSE, col = "grey85", border = "white",
     xlim = xl,
     main = "Standardised Residuals (K=1 gjrGARCH)",
     xlab = expression(hat(eta)[t]))
xg <- seq(xl[1], xl[2], length.out = 400)
lines(xg, dnorm(xg), col = "steelblue", lwd = 1.8)
# Student-t with model's estimated df
lines(xg, dt(xg, df = nu_k1), col = "firebrick", lwd = 1.8, lty = 2)
legend("topright",
       legend = c("N(0,1)", sprintf("t(df=%.2f)", nu_k1)),
       col = c("steelblue", "firebrick"), lty = 1:2, lwd = 1.8,
       cex = 0.8, bg = "white")

# Q-Q against Student-t (not Normal, since model assumes t)
qqplot(qt(ppoints(length(eta_k1_plot)), df = nu_k1),
       eta_k1_plot,
       pch = 1, cex = 0.3, col = "steelblue",
       main = sprintf("Q-Q vs Student-t(df=%.2f) -- K=1", nu_k1),
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "firebrick", lwd = 1.5)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
cat("  Saved: residuals_diag_k1.pdf\n")


# ── (e) OOS VaR backtest plot (K=1 gjrGARCH) ─────────────────────────────────

ret_oos_plot <- phase4b_results$oos_returns
dates_oos_plot <- phase4b_results$oos_dates
var05_oos    <- phase4b_results$VaR_oos_k1_05
var01_oos    <- phase4b_results$VaR_oos_k1_01

vio_05_oos <- which(ret_oos_plot < var05_oos)
vio_01_oos <- which(ret_oos_plot < var01_oos)

pdf(file.path(output_dir, "oos_var_plot_k1.pdf"),
    width = 11, height = 5)
plot(dates_oos_plot, ret_oos_plot, type = "l", col = "grey50", lwd = 0.5,
     main = sprintf(paste0("OOS Returns vs VaR -- K=1 gjrGARCH  |  ",
                           "alpha=0.05 hit=%.4f  alpha=0.01 hit=%.4f"),
                    mean(ret_oos_plot < var05_oos, na.rm = TRUE),
                    mean(ret_oos_plot < var01_oos, na.rm = TRUE)),
     xlab = "Date", ylab = "Log-Return (OLS-residual)")
lines(dates_oos_plot, var05_oos, col = "firebrick",  lwd = 1.2, lty = 2)
lines(dates_oos_plot, var01_oos, col = "darkorange", lwd = 1.2, lty = 3)
points(dates_oos_plot[vio_05_oos], ret_oos_plot[vio_05_oos],
       col = "firebrick", pch = 19, cex = 0.4)
points(dates_oos_plot[vio_01_oos], ret_oos_plot[vio_01_oos],
       col = "darkorange", pch = 17, cex = 0.6)
abline(h = 0, col = "grey80", lty = 1)
legend("bottomleft",
       legend = c("OOS Return", "VaR alpha=0.05", "VaR alpha=0.01",
                  "Violation alpha=0.05", "Violation alpha=0.01"),
       col = c("grey50", "firebrick", "darkorange", "firebrick", "darkorange"),
       lty = c(1, 2, 3, NA, NA),
       pch = c(NA, NA, NA, 19, 17),
       cex = 0.75, bg = "white")
dev.off()
cat("  Saved: oos_var_plot_k1.pdf\n")

cat("\nAll paper plots saved to:", output_dir, "\n")
