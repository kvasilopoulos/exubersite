# Validation script for lbi_test() (Breitung & Diegel 2025's static
# locally best invariant test). See docs/enhancements/monitoring.md for
# the full write-up. Run from the exuber/ package root (or adjust the
# devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Breitung & Diegel's eq. 4 telescoping identity ===\n")
cat("(2*sum(Delta y_t * y_{t-1}) = y_T^2 - T*sigma_tilde^2, y_1=0 case)\n\n")
set.seed(1)
T <- 100
y <- c(0, cumsum(rnorm(T)))
dy <- diff(y)
ylag <- y[1:T]
lhs <- 2 * sum(dy * ylag)
sigma2_tilde <- mean(dy^2)
rhs <- y[T + 1]^2 - T * sigma2_tilde
cat("LHS:", lhs, " RHS:", rhs, " match:", isTRUE(all.equal(lhs, rhs)), "\n\n")

cat("=== 2. Basic run ===\n")
set.seed(1)
y2 <- cumsum(rnorm(100))
out <- lbi_test(y2)
print(out)

cat("\n=== 3. H0 calibration: does the statistic actually follow N(0,1)? ===\n")
cat("(the paper claims a standard normal null distribution directly --\n")
cat("this checks that claim, not just an approximately-correct size)\n\n")
run_stat <- function(seed) {
  set.seed(seed)
  y <- cumsum(rnorm(100))
  lbi_test(y)$stat[["series1"]]
}
stats <- sapply(1:500, run_stat)
cat(sprintf("mean=%.3f (expect ~0), sd=%.3f (expect ~1)\n", mean(stats), sd(stats)))
cat(sprintf("empirical P(stat > qnorm(0.95)=%.3f): %.3f (expect ~0.05)\n", qnorm(0.95), mean(stats > qnorm(0.95))))
cat(sprintf("KS test vs N(0,1): p-value = %.3f (expect not tiny)\n", ks.test(stats, "pnorm")$p.value))

cat("\n=== 4. Power under a genuine explosive alternative, vs standard SADF ===\n")
run_power_lbi <- function(seed) {
  set.seed(seed)
  n1 <- 60
  y <- 100 * 1.03^(1:n1) + cumsum(rnorm(n1, sd = 1))
  lbi_test(y)$detected[["series1"]]
}
run_power_sadf <- function(seed) {
  set.seed(seed)
  n1 <- 60
  y <- 100 * 1.03^(1:n1) + cumsum(rnorm(n1, sd = 1))
  full <- radf(y, minw = 20)
  cv <- radf_mc_cv(n = n1, minw = 20, nrep = 200, seed = 1)
  full$sadf > cv$sadf_cv["95%"]
}
cat(sprintf("Detection rate, LBI:  %.3f\n", mean(sapply(1:60, run_power_lbi))))
cat(sprintf("Detection rate, SADF: %.3f\n", mean(sapply(1:60, run_power_sadf))))
