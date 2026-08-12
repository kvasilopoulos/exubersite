# Validation script for radf_quantile() (Wu, Shi & Wu 2025, "Quantile
# analysis for financial bubble detection and surveillance", the
# "global test" of their Section 3.1). See docs/enhancements/
# alternative-paradigms.md, "Quantile-based detection", for the full
# write-up. Run from the exuber/ package root (or adjust the
# devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Structural check: the critical value's Q functional is\n")
cat("exactly radf()'s own single-shot adf t-statistic, not a new\n")
cat("computation (bit-for-bit) ===\n\n")
set.seed(1)
y <- cumsum(rnorm(100))
full <- radf(y, minw = 90)
q_manual <- exuber:::quantile_adf_tstat(y)
cat("radf()$adf:", full$adf, " manual Q:", q_manual, " match:",
    isTRUE(all.equal(unname(full$adf), q_manual, tolerance = 1e-8)), "\n\n")

cat("=== 2. Empirical size under H0 (pure random walk), tau='optimal' ===\n")
run_null <- function(seed) {
  set.seed(seed)
  y <- cumsum(rnorm(100))
  radf_quantile(y, nrep = 200, seed = 1)$detected[["series1"]]
}
rate_h0 <- mean(sapply(1:100, run_null))
cat(sprintf("False-detection rate (100 reps, nominal 5%%): %.3f\n\n", rate_h0))

cat("=== 3. Empirical size under H0, fixed tau=0.5 ===\n")
run_null_fixed <- function(seed) {
  set.seed(seed)
  y <- cumsum(rnorm(100))
  radf_quantile(y, tau = 0.5, nrep = 200, seed = 1)$detected[["series1"]]
}
rate_h0_fixed <- mean(sapply(1:100, run_null_fixed))
cat(sprintf("False-detection rate at tau=0.5 (nominal 5%%): %.3f\n\n", rate_h0_fixed))

cat("=== 4. Power under a genuine explosive alternative, vs standard SADF ===\n")
run_power <- function(seed) {
  set.seed(seed)
  n1 <- 60
  y <- 100 * 1.03^(1:n1) + cumsum(rnorm(n1, sd = 1))
  radf_quantile(y, nrep = 200, seed = 1)$detected[["series1"]]
}
run_power_sadf <- function(seed) {
  set.seed(seed)
  n1 <- 60
  y <- 100 * 1.03^(1:n1) + cumsum(rnorm(n1, sd = 1))
  full <- radf(y, minw = 20)
  cv <- radf_mc_cv(n = n1, minw = 20, nrep = 200, seed = 1)
  full$sadf > cv$sadf_cv["95%"]
}
cat(sprintf("Detection rate, radf_quantile: %.3f\n", mean(sapply(1:60, run_power))))
cat(sprintf("Detection rate, standard SADF (same DGP, rough cross-check): %.3f\n\n",
            mean(sapply(1:60, run_power_sadf))))

cat("=== 5. Optimal-tau selection lands sensibly within the search grid ===\n")
taus <- sapply(1:20, function(s) {
  set.seed(s)
  y <- cumsum(rnorm(120))
  radf_quantile(y, nrep = 100, seed = 1)$tau[["series1"]]
})
cat("selected taus:", paste(taus, collapse = ", "), "\n")
cat("all within [0.2, 0.8]:", all(taus >= 0.2 & taus <= 0.8), "\n")
