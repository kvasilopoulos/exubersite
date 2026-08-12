# Validation script for radf_knp() (Kejriwal, Nguyen & Perron 2025, "An
# Improved Procedure for Retrospectively Dating the Emergence and
# Collapse of Bubbles"). See docs/enhancements/dating-and-root-
# inference.md, "SSR/BIC dating vs. PSY recursive dating", for the full
# write-up. Run from the exuber/ package root (or adjust the
# devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Formula-exact: knp_find_break() vs brute-force nested lm() search ===\n")
set.seed(3)
n <- 26
y <- cumsum(rnorm(n))
n1 <- n - 1L
x <- y[1:n1]; z <- y[2:(n1 + 1)] - y[1:n1]
k_min <- max(2L, ceiling(0.1 * n1))

for (omit in c(FALSE, TRUE)) {
  fit <- exuber:::knp_find_break(y, trim = 0.1, omit = omit)
  best <- list(ssr = Inf)
  for (tau1 in k_min:(n1 - 2 * k_min)) {
    for (tau2 in (tau1 + k_min):(n1 - k_min)) {
      idx_mid <- (tau1 + 1):tau2
      ssr <- sum(z[1:tau1]^2) + sum(resid(lm(z[idx_mid] ~ x[idx_mid]))^2) + sum(z[(tau2 + 1):n1]^2)
      if (omit) ssr <- ssr - z[tau2 + 1]^2
      if (ssr < best$ssr) best <- list(tau1 = tau1, tau2 = tau2, ssr = ssr)
    }
  }
  cat(sprintf("omit=%s: vectorized=(%d,%d,%.4f)  brute=(%d,%d,%.4f)\n",
              omit, fit$tau1, fit$tau2, fit$ssr, best$tau1, best$tau2, best$ssr))
}

cat("\n=== 2. Reproducing Theorem 1 (naive inconsistency) vs Theorem 2\n")
cat("     (omission-corrected consistency) ===\n")
cat("(KNP's own DGP: unit root -> no-intercept explosive AR(1) -> an\n")
cat("instantaneous collapse back near the pre-bubble level -> fresh unit\n")
cat("root. Theorem 1 proves plain OLS's origination-date estimate\n")
cat("converges to the TRUE COLLAPSE date, not the true origination date;\n")
cat("Theorem 2 proves the single-residual omission fixes this.)\n\n")

sim_knp <- function(seed, T1 = 50, T2 = 90, T = 200, delta = 1.05) {
  set.seed(seed)
  y <- numeric(T)
  y[1] <- 0
  for (t in 2:T1) y[t] <- y[t - 1] + rnorm(1)
  for (t in (T1 + 1):T2) y[t] <- delta * y[t - 1] + rnorm(1)
  y[T2 + 1] <- y[T1] + rnorm(1)
  if (T2 + 2 <= T) for (t in (T2 + 2):T) y[t] <- y[t - 1] + rnorm(1)
  list(y = y, T1 = T1, T2 = T2)
}
run <- function(seed, omit) {
  sim <- sim_knp(seed)
  fit <- exuber:::knp_find_break(sim$y, trim = 0.05, omit = omit)
  c(tau1 = fit$tau1, tau2 = fit$tau2, T1 = sim$T1, T2 = sim$T2)
}
res_naive <- t(sapply(1:30, run, omit = FALSE))
res_om <- t(sapply(1:30, run, omit = TRUE))

cat(sprintf("Naive (omit=FALSE):     mean|tau1-T1|=%.1f   mean|tau1-T2|=%.1f  (tau1 should track T2, not T1)\n",
            mean(abs(res_naive[, "tau1"] - res_naive[, "T1"])),
            mean(abs(res_naive[, "tau1"] - res_naive[, "T2"]))))
cat(sprintf("Omission-corrected:     mean|tau1-T1|=%.2f   mean|tau2-T2|=%.2f\n\n",
            mean(abs(res_om[, "tau1"] - res_om[, "T1"])),
            mean(abs(res_om[, "tau2"] - res_om[, "T2"]))))

cat("=== 3. delta_hat accuracy under omission correction (true delta=1.05) ===\n")
run_delta <- function(seed) {
  sim <- sim_knp(seed)
  out <- radf_knp(sim$y, trim = 0.05, omit = TRUE)
  unname(out$delta[["series1"]])
}
deltas <- sapply(1:30, run_delta)
cat(sprintf("mean delta_hat = %.3f (true = 1.05), mean|bias| = %.3f\n", mean(deltas), mean(abs(deltas - 1.05))))
