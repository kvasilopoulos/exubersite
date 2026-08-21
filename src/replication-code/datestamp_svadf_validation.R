# Validation of datestamp(option = "svadf") -- Sarkar & Wells (2026,
# preprint)'s SV-ADF asymmetric-threshold bubble dating. See docs/
# enhancements/volatility-robustness.md, "SV-ADF", for the full writeup.
# Originally shipped as its own radf_svadf() entry point (2026-08-11);
# folded into datestamp() as an option (2026-08-18) -- see exuber/CLAUDE.md,
# "Naming: not everything is radf_* anymore".
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Structural: badf reused bit-for-bit from radf() ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
r <- radf(y, lag = 0)
r2 <- radf(y, minw = attr(r, "minw"), lag = 0)
cat("max|badf diff|:", max(abs(r$badf[, 1] - r2$badf[, 1])), "\n")

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
  rr <- radf(yy, lag = 0)
  out <- datestamp(rr, option = "svadf", min_duration = psy_ds(150))
  if (length(out) > 0 && out[[1]]$End[1] <= out[[1]]$Start[1]) ok <- FALSE
}
cat("collapse always after origination when detected:", ok, "\n")

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
  rr <- radf(yy, lag = 0)
  out <- datestamp(rr, option = "svadf", min_duration = psy_ds(length(yy)))
  if (length(out) > 0) {
    detected <- detected + 1
    orig_err <- c(orig_err, abs(out[[1]]$Start[1] - n1))
    if (!isTRUE(out[[1]]$Ongoing[1])) coll_err <- c(coll_err, abs(out[[1]]$End[1] - (n1 + n2)))
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
  rr <- radf(yy, lag = 0)
  out <- datestamp(rr, option = "svadf", min_duration = psy_ds(150))
  if (length(out) > 0) fa <- fa + 1
}
cat("false origination-alarm rate:", fa / nrep2, "\n")

cat("\ndone\n")
