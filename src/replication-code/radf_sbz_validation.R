# Replication script for radf_sbz_union() (SBZ: WLS + kernel volatility,
# Harvey, Leybourne & Zu 2019). See docs/volatility-robustness.md,
# "SBZ (WLS + kernel volatility)". The script checks the empirical size under
# a pure random walk, where supDF, supBZ and U should reject at about the
# nominal 5% rate.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== Empirical size under H0 (pure random walk, no bubble) ===\n")
cat("Target: ~0.05 nominal for each of supDF/supBZ/U\n")
cat("After-fix numbers reported in the doc: supDF=0.033, supBZ=0.060, U=0.053\n\n")

set.seed(13579)
n <- 150
nrep <- 150
nboot <- 199

p_supDF <- p_supBZ <- p_U <- numeric(nrep)
for (i in seq_len(nrep)) {
  y <- cumsum(rnorm(n))
  res <- radf_sbz_union(y, minw = 20, nboot = nboot, seed = NULL)
  p_supDF[i] <- res$p_supDF
  p_supBZ[i] <- res$p_supBZ
  p_U[i] <- res$p_U
}

cat(sprintf("supDF empirical rejection rate: %.3f\n", mean(p_supDF < 0.05)))
cat(sprintf("supBZ empirical rejection rate: %.3f\n", mean(p_supBZ < 0.05)))
cat(sprintf("U     empirical rejection rate: %.3f\n", mean(p_U < 0.05)))

cat("\n=== Full test-sbz.R suite ===\n")
testthat::test_file(
  "exuber/tests/testthat/test-sbz.R",
  reporter = "summary"
)
