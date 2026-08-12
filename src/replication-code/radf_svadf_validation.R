# Validation of radf_svadf() -- Sarkar & Wells (2026, preprint)'s SV-ADF
# asymmetric-threshold bubble dating. See docs/enhancements/
# volatility-robustness.md, "SV-ADF", for the full writeup.
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Structural: badf reused bit-for-bit from radf() ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
out <- radf_svadf(y)
r <- radf(y, minw = attr(out, "minw"), lag = 0)
cat("max|badf diff|:", max(abs(out$badf[, 1] - r$badf[, 1])), "\n")

cat("\n=== 2. Threshold formulas ===\n")
cat("svadf_threshold(100,'origination'):", exuber:::svadf_threshold(100, "origination"),
  " expect log(100)/10=", log(100) / 10, "\n")
cat("svadf_threshold(100,'collapse'):", exuber:::svadf_threshold(100, "collapse"),
  " expect log(100)/2=", log(100) / 2, "\n")

cat("\n=== 3. Collapse never dated before origination (20 reps) ===\n")
set.seed(3)
ok <- TRUE
for (i in 1:20) {
  yy <- cumsum(rnorm(150))
  oo <- radf_svadf(yy)
  if (!is.na(oo$origination) && !is.na(oo$collapse) && oo$collapse <= oo$origination) ok <- FALSE
}
cat("collapse always after origination when both detected:", ok, "\n")

cat("\n=== 4. Dating accuracy on a synthetic bubble+collapse episode (20 reps) ===\n")
orig_err <- coll_err <- c()
detected <- 0
nrep <- 20
for (i in 1:nrep) {
  set.seed(100 + i)
  n1 <- 60
  yy1 <- 100 + cumsum(rnorm(n1))
  n2 <- 40
  bubble <- yy1[n1] * 1.04^(1:n2) + cumsum(rnorm(n2, sd = 1))
  n3 <- 40
  coll <- bubble[n2] - cumsum(abs(rnorm(n3, mean = 3, sd = 1)))
  yy <- c(yy1, bubble, coll)
  out <- radf_svadf(yy)
  if (!is.na(out$origination)) {
    detected <- detected + 1
    orig_err <- c(orig_err, abs(out$origination - n1))
    if (!is.na(out$collapse)) coll_err <- c(coll_err, abs(out$collapse - (n1 + n2)))
  }
}
cat("detection rate:", detected / nrep, "\n")
cat("mean |origination error|:", mean(orig_err), "\n")
cat("mean |collapse error|:", mean(coll_err), "\n")

cat("\n=== 5. False-alarm rate under H0 (pure random walk, 60 reps) ===\n")
set.seed(5)
fa <- 0
nrep2 <- 60
for (i in 1:nrep2) {
  set.seed(1000 + i)
  yy <- cumsum(rnorm(150))
  out <- radf_svadf(yy)
  if (!is.na(out$origination)) fa <- fa + 1
}
cat("false origination-alarm rate:", fa / nrep2, "\n")

cat("\ndone\n")
