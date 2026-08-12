# Replication script for radf_common()/radf_common_cv() (common-bubble
# detection via PCA + PSY, Chen, Phillips & Shi 2023). Archived
# retroactively -- see docs/enhancements/multivariate.md, "Common-bubble
# detection via PCA + PSY", "Independent validation (2026-08-09)" for the
# narrative this reproduces: Theorem 4.3 claims the PSY-on-PC1 statistic's
# null is asymptotically identical to the plain univariate GSADF null
# (independent of panel width N), but this doesn't hold at practical N --
# the true null quantile *grows* with N (opposite of what a naive
# "converges as N grows" reading of the paper's own finite-sample section
# would suggest), which is why radf_common_cv() simulates its own,
# N-dependent null rather than reusing radf_mc_cv().
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== radf_common()'s null quantiles grow with panel width N ===\n")
cat("(sharpest possible null: N independent random walks, no true common factor)\n")
cat("T = 250 throughout, seed = 271828\n\n")

cv_univariate <- radf_mc_cv(250, nrep = 2000, seed = 1)
cat("radf_mc_cv(250)'s 95% gsadf_cv (the univariate benchmark):",
  round(cv_univariate$gsadf_cv["95%"], 3), "\n\n"
)

set.seed(271828)
results <- list()
for (N in c(6, 20, 50, 100)) {
  gsadf_null <- replicate(400, {
    panel <- replicate(N, cumsum(rnorm(250)))
    radf_common(panel)$gsadf
  })
  q <- quantile(gsadf_null, c(0.9, 0.95, 0.99))
  results[[as.character(N)]] <- q
  cat(sprintf(
    "N=%3d: 90%%=%.3f 95%%=%.3f 99%%=%.3f  (diff from radf_mc_cv 95%%: %+.2f)\n",
    N, q[1], q[2], q[3], q[2] - cv_univariate$gsadf_cv["95%"]
  ))
}

cat("\nPublished (doc): N=6 -> (2.420,2.662,3.057); N=20 -> (3.225,3.451,4.144);\n")
cat("                 N=50 -> (4.057,4.365,4.921); N=100 -> (4.880,5.203,5.836)\n\n")

cat("=== radf_common_cv() reproduces this N-dependence directly ===\n")
cv_small <- radf_common_cv(n = 250, N = 4, minw = NULL, nrep = 500, seed = 1)
cv_large <- radf_common_cv(n = 250, N = 30, minw = NULL, nrep = 500, seed = 1)
cat("radf_common_cv(N=4)  gsadf_cv:", round(cv_small$gsadf_cv, 3), "\n")
cat("radf_common_cv(N=30) gsadf_cv:", round(cv_large$gsadf_cv, 3), "\n")
cat("N=30's null quantile is higher than N=4's:",
  all(cv_large$gsadf_cv > cv_small$gsadf_cv), "\n\n"
)

cat("=== Full test-common.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-common.R",
  reporter = "summary"
)
