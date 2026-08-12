# Replication script for radf_kp() (kernel-purge test, Harvey, Leybourne,
# Taylor & Zu 2024). Archived retroactively -- see
# docs/enhancements/volatility-robustness.md, "Kernel-purge test",
# "Independent validation (2026-08-09)" for the narrative this reproduces.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== Critical-value replication vs Table I (T = 400, PSY_sigma) ===\n")
cat("Published: 1.712 / 1.935 / 2.296 (10%/5%/1%)\n\n")

set.seed(31415)
gsadf_kp <- replicate(800, radf_kp(cumsum(rnorm(400)))$gsadf)
own <- quantile(gsadf_kp, c(0.9, 0.95, 0.99))
cat("Independent run (n=400, nrep=800, seed=31415):", round(own, 3), "\n")
cat("Abs. difference from published:", round(abs(own - c(1.712, 1.935, 2.296)), 3), "\n\n")

cat("=== Consistency with radf_mc_cv() (Remark 3.2: T=Inf row should match\n")
cat("    the standard homoskedastic GSADF null) ===\n")
set.seed(2)
gsadf_kp_300 <- replicate(500, radf_kp(cumsum(rnorm(300)))$gsadf)
cv_mc_300 <- radf_mc_cv(300, nrep = 500, seed = 2)
cat("radf_kp() null quantiles (n=300):", round(quantile(gsadf_kp_300, c(0.9, 0.95, 0.99)), 3), "\n")
cat("radf_mc_cv(300) quantiles:       ", round(cv_mc_300$gsadf_cv, 3), "\n\n")

cat("=== Full test-kp.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-kp.R",
  reporter = "summary"
)
