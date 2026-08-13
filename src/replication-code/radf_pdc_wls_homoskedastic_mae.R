devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

run_once <- function(seed) {
  set.seed(seed)
  n1_len <- 150; n2_len <- 80; n3_len <- 100
  regime1 <- cumsum(rnorm(n1_len, sd = 0.5))
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
mae <- colMeans(abs(res))
cat("Homoskedastic case -- mean absolute error across 40 seeds:\n")
print(mae)
cat(sprintf(
  "\nOrigination MAE: OLS=%.2f  WLS=%.2f  (should be close -- no volatility signal to exploit)\n",
  mae["ols_orig_err"], mae["wls_orig_err"]
))
cat(sprintf(
  "Collapse MAE:     OLS=%.2f  WLS=%.2f\n",
  mae["ols_coll_err"], mae["wls_coll_err"]
))
