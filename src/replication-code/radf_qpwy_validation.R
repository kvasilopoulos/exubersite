# Validation of radf_qpwy() -- Wu, Shi & Wu (2025)'s QPWY recursive
# quantile monitoring strategy (single-recursion subset; QPSY's double
# recursion is not implemented). See docs/enhancements/
# alternative-paradigms.md, "Quantile-based detection", "Implementation
# (QPWY)", for the full writeup -- including a real bug (per-r marginal
# quantile used as boundary instead of a supremum-calibrated one) found
# and fixed via the Monte Carlo false-alarm-rate check in this script.
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Point statistic: qpwy_stat_path at the full sample matches
radf_quantile()'s own per-window QR t-ratio formula exactly ===\n")
set.seed(3)
y <- cumsum(rnorm(60))
full_stat <- exuber:::qpwy_stat_path(y, 0.5, 60)

dy <- diff(y)
ylag <- y[1:59]
yresp <- y[2:60]
qr_fit <- quantreg::rq(yresp ~ ylag, tau = 0.5)
alpha_hat <- unname(coef(qr_fit)["ylag"])
f_hat <- exuber:::quantile_check_density(dy, 0.5)$f_hat
yPzy <- sum((ylag - mean(ylag))^2)
manual <- (f_hat / sqrt(0.5 * 0.5)) * sqrt(yPzy) * (alpha_hat - 1)
cat(sprintf("fast=%.8f manual=%.8f |diff|=%.2e\n", full_stat, manual, abs(full_stat - manual)))

cat("\n=== 2. qpwy_boundary_sim reuses radf()'s own badf sequence (Corollary 2) ===\n")
n <- 80
minw <- 20
Q <- exuber:::qpwy_boundary_sim(n, minw, 5, seed = 1)
set.seed(1)
ysim1 <- cumsum(rnorm(n))
r1 <- radf(ysim1, minw = minw, lag = 0)
cat("ncol(Q):", ncol(Q), " length(badf):", length(r1$badf[, 1]), " match:", ncol(Q) == length(r1$badf[, 1]), "\n")

cat("\n=== 3. THE BUG: per-r marginal quantile vs. supremum-calibrated boundary ===\n")
Q2 <- exuber:::qpwy_boundary_sim(80, minw, 200, seed = 7)
set.seed(9)
z <- rnorm(200)
delta_j <- 0.4
U <- sqrt(1 - delta_j^2) * z + delta_j * Q2
marginal_boundary_lastcol <- unname(quantile(U[, ncol(U)], probs = 0.95, names = FALSE))
sup_boundary <- unname(quantile(apply(U, 1, max), probs = 0.95, names = FALSE))
cat(sprintf(
  "marginal (per-r) boundary at last r: %.3f   supremum-calibrated boundary: %.3f\n",
  marginal_boundary_lastcol, sup_boundary
))
cat("(the fix uses the supremum-calibrated boundary; using the marginal one at every\n")
cat(" r instead was the original bug, giving a ~50%% false-alarm rate against nominal 5%%)\n")

cat("\n=== 4. Empirical false-alarm rate under H0 (post-fix) ===\n")
set.seed(2)
nrep_mc <- 60
n <- 150
fa <- 0
for (i in seq_len(nrep_mc)) {
  set.seed(2000 + i)
  yy <- cumsum(rnorm(n))
  out <- radf_qpwy(yy, tau = 0.5, nrep = 150, seed = i)
  if (!is.na(out$alarm)) fa <- fa + 1
}
cat(sprintf("false-alarm rate (nominal 5%%): %.3f\n", fa / nrep_mc))

cat("\n=== 5. Detection power on a genuine explosive DGP, vs. standard SADF ===\n")
set.seed(3)
det_qpwy <- det_sadf <- 0
for (i in seq_len(nrep_mc)) {
  set.seed(3000 + i)
  n1 <- 100
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.03^(1:50) + cumsum(rnorm(50, sd = 1))
  yy <- c(normal_part, expl_part)
  out <- radf_qpwy(yy, tau = 0.5, nrep = 150, seed = i)
  r <- radf(yy)
  cv <- radf_mc_cv(length(yy))
  if (!is.na(out$alarm)) det_qpwy <- det_qpwy + 1
  if (r$sadf > cv$sadf_cv[2]) det_sadf <- det_sadf + 1
}
cat(sprintf("QPWY power: %.3f   SADF power: %.3f\n", det_qpwy / nrep_mc, det_sadf / nrep_mc))

cat("\ndone\n")
