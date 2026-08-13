devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

# DGP: 3-regime bubble (unit root -> explosive -> stationary collapse), with
# a volatility burst concentrated in the FIRST 20% of regime 1 (sd = high),
# dropping to a much smaller sd for the rest of regime 1. This mirrors the
# KS(2023) WLS paper's own finding that the correction helps most when a
# volatility break sits near the start/end of the sample: the noisy early
# segment inflates OLS's unweighted objective and can bias the origination
# split; WLS should downweight it via the estimated spot variance.
run_once <- function(seed) {
  set.seed(seed)
  n1_len <- 150; n2_len <- 80; n3_len <- 100
  burst_len <- round(0.2 * n1_len)

  e1 <- c(rnorm(burst_len, sd = 4), rnorm(n1_len - burst_len, sd = 0.3))
  regime1 <- cumsum(e1)
  regime2 <- regime1[n1_len] * 1.07^(1:n2_len) + cumsum(rnorm(n2_len, sd = 0.15))
  peak <- regime2[n2_len]
  rho3 <- 0.5
  regime3 <- numeric(n3_len)
  regime3[1] <- rho3 * peak + rnorm(1, sd = 0.5)
  for (t in 2:n3_len) regime3[t] <- rho3 * regime3[t - 1] + rnorm(1, sd = 0.5)

  y <- c(regime1, regime2, regime3)
  true_origination <- n1_len
  true_collapse <- n1_len + n2_len

  out_ols <- dating_pdc(y, regimes = 3L, trim = 0.05, type = "ols")
  out_wls <- dating_pdc(y, regimes = 3L, trim = 0.05, type = "wls")

  c(
    ols_orig_err = out_ols$origination - true_origination,
    ols_coll_err = out_ols$collapse - true_collapse,
    wls_orig_err = out_wls$origination - true_origination,
    wls_coll_err = out_wls$collapse - true_collapse
  )
}

res <- t(sapply(1:40, run_once))
cat("Per-seed errors (first 10 rows):\n")
print(head(res, 10))

mae <- colMeans(abs(res))
cat("\nMean absolute error across 40 seeds:\n")
print(mae)

cat(sprintf(
  "\nOrigination MAE: OLS=%.2f  WLS=%.2f  (WLS better if smaller)\n",
  mae["ols_orig_err"], mae["wls_orig_err"]
))
cat(sprintf(
  "Collapse MAE:     OLS=%.2f  WLS=%.2f  (WLS better if smaller)\n",
  mae["ols_coll_err"], mae["wls_coll_err"]
))
