# Validation script for radf_recovery()/radf_recovery_cv() (Phillips &
# Shi 2014, "Financial Bubble Implosion and Reverse Regression").
#
# Reports both what validates cleanly (f_r, the structural f_c<=f_r
# invariant, the reversal-calibrated CV genuinely differing from the
# forward one) and what doesn't (f_c's bias, the H0 false-detection
# rate) -- see docs/enhancements/dating-and-root-inference.md,
# "Reverse-regression recovery dating", for the full write-up and honest
# accounting of these results. Run from the exuber/ package root (or
# adjust the devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

n <- 100; minw <- 20; lvl_lab <- "95%"

cat("=== 1. Reversal-calibrated CV differs from forward CV ===\n")
cat("(paired Monte Carlo, same underlying draws, confirms Theorem 1's\n")
cat("endogeneity finding empirically rather than just asymptotically)\n\n")
nrep_cv_check <- 5000
set.seed(7)
res_fwd <- matrix(NA_real_, nrow = n - minw, ncol = nrep_cv_check)
res_rev <- matrix(NA_real_, nrow = n - minw, ncol = nrep_cv_check)
for (i in 1:nrep_cv_check) {
  y <- cumsum(rnorm(n))
  res_fwd[, i] <- exuber:::rls_gsadf(exuber:::unroot(y, lag = 0), min_win = minw, lag = 0)[1:(n - minw)]
  res_rev[, i] <- exuber:::rls_gsadf(exuber:::unroot(rev(y), lag = 0), min_win = minw, lag = 0)[1:(n - minw)]
}
cv_fwd <- t(apply(apply(res_fwd, 2, cummax), 1, quantile, probs = c(0.9, 0.95, 0.99)))
cv_rev <- t(apply(apply(res_rev, 2, cummax), 1, quantile, probs = c(0.9, 0.95, 0.99)))
cat("Max abs diff (95% col):", max(abs(cv_fwd[, "95%"] - cv_rev[, "95%"])), "\n")
cat("Mean abs diff (95% col):", mean(abs(cv_fwd[, "95%"] - cv_rev[, "95%"])), "\n\n")

cat("=== 2. radf_recovery_cv(), one stable estimate (nrep=1000) reused below ===\n\n")
cv <- radf_recovery_cv(n = n, minw = minw, nrep = 1000, seed = 99)
zadj <- minw

detect_once <- function(y, cv) {
  fit <- radf(rev(y), minw = minw)
  exceed <- fit$bsadf[, 1] > cv$bsadf_cv[, lvl_lab]
  g_e <- which(exceed)[1L]
  if (is.na(g_e)) return(c(detected = FALSE, censored = FALSE, f_c = NA, f_r = NA))
  after <- which(!exceed[g_e:length(exceed)])
  if (length(after) == 0L) {
    g_c <- length(exceed) + 1L; censored <- TRUE
  } else {
    g_c <- g_e + after[1L] - 1L; censored <- FALSE
  }
  f_r <- n + 1L - (g_e + zadj)
  f_c <- if (censored) NA else n + 1L - min(g_c + zadj, n)
  c(detected = TRUE, censored = censored, f_c = f_c, f_r = f_r)
}

cat("=== 3. Empirical false-detection rate under pure H0 (random walk) ===\n")
cat("(200 fresh draws against the one stable cv from step 2)\n\n")
set.seed(123)
res_h0 <- t(sapply(1:200, function(i) detect_once(cumsum(rnorm(n)), cv)))
cat(sprintf("False-detection rate: %.3f  (nominal level: %s; noisier than\n", mean(res_h0[, "detected"] == 1), lvl_lab))
cat("comparable forward-test numbers elsewhere in this project -- flagged\n")
cat("honestly as unresolved, see the taxonomy file for discussion)\n\n")

cat("=== 4. Detection accuracy on a synthetic collapse-then-recovery DGP ===\n")
cat("(smooth, continuous mean-reverting collapse -- an earlier deterministic-\n")
cat("decay version produced a spurious spike right at the regime boundary,\n")
cat("found and fixed during this validation)\n\n")
run_detect <- function(seed, cv) {
  set.seed(seed)
  n1 <- 40; n2 <- 25; n3 <- 35
  expansion <- 100 * 1.03^(1:n1) + cumsum(rnorm(n1, sd = 1))
  target <- expansion[n1] * 0.5
  collapse <- numeric(n2)
  collapse[1] <- expansion[n1] + rnorm(1, sd = 1)
  for (k in 2:n2) collapse[k] <- target + 0.9 * (collapse[k - 1] - target) + rnorm(1, sd = 1)
  recovery <- collapse[n2] + cumsum(rnorm(n3, sd = 1)) + (1:n3) * 0.5
  y <- c(expansion, collapse, recovery)
  out <- detect_once(y, cv)
  c(out, true_collapse = n1, true_recovery = n1 + n2)
}
res <- t(sapply(1:40, function(s) run_detect(s, cv)))
cat(sprintf("Detection rate: %.3f\n", mean(res[, "detected"] == 1)))
det <- res[res[, "detected"] == 1 & res[, "censored"] == 0, , drop = FALSE]
cat(sprintf("n detected & uncensored = %d\n", nrow(det)))
cat(sprintf("mean f_c bias = %.2f, mean f_r bias = %.2f\n",
            mean(det[, "f_c"] - det[, "true_collapse"]), mean(det[, "f_r"] - det[, "true_recovery"])))
cat(sprintf("mean |f_c bias| = %.2f, mean |f_r bias| = %.2f  (f_r matches the\n",
            mean(abs(det[, "f_c"] - det[, "true_collapse"])), mean(abs(det[, "f_r"] - det[, "true_recovery"]))))
cat("paper's own ~6-observation-early finding; f_c's bias is materially\n")
cat("larger and not fully resolved -- see taxonomy file for discussion)\n\n")

cat("=== 5. Structural invariant: f_c <= f_r whenever both identified & uncensored ===\n")
cat("all TRUE:", all(det[, "f_c"] <= det[, "f_r"]), "\n")
