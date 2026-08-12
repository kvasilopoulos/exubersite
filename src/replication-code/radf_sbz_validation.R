# Replication script for radf_sbz_cv() (SBZ: WLS + kernel volatility,
# Harvey, Leybourne & Zu 2019). Archived retroactively -- see
# docs/enhancements/volatility-robustness.md, "SBZ (WLS + kernel
# volatility)", "Independent validation (2026-08-09) -- found and fixed a
# real bug" for the narrative this reproduces: an off-by-one bootstrap
# -indexing bug (`pointer <- length(ystar) - 1L - minw`, missing the -1L)
# that silently extracted `adf` instead of `sadf` from rls_gsadf()'s flat
# result vector, badly oversizing supDF/U. Fixed in radf_sbz.R; this script
# reproduces the empirical-size-under-H0 check that caught it, on the
# already-fixed code (so it should reproduce the "after fix" column, not
# the original 0.640/0.573 broken numbers).
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== Empirical size under H0 (pure random walk, no bubble) ===\n")
cat("Target: ~0.05 nominal for each of supDF/supBZ/U\n")
cat("After-fix numbers reported in the doc: supDF=0.033, supBZ=0.080, U=0.060\n\n")

set.seed(13579)
n <- 150
nrep <- 150
nboot <- 199

p_supDF <- p_supBZ <- p_U <- numeric(nrep)
for (i in seq_len(nrep)) {
  y <- cumsum(rnorm(n))
  res <- radf_sbz_cv(y, minw = 20, nboot = nboot, seed = NULL)
  p_supDF[i] <- res$p_supDF
  p_supBZ[i] <- res$p_supBZ
  p_U[i] <- res$p_U
}

cat(sprintf("supDF empirical rejection rate: %.3f\n", mean(p_supDF < 0.05)))
cat(sprintf("supBZ empirical rejection rate: %.3f\n", mean(p_supBZ < 0.05)))
cat(sprintf("U     empirical rejection rate: %.3f\n", mean(p_U < 0.05)))

cat("\n=== Full test-sbz.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-sbz.R",
  reporter = "summary"
)
