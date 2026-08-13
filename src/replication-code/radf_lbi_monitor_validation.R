# Validation of monitor_lbi() -- Breitung & Diegel (2025)'s sequential
# (constant-boundary mCUSUM/wCUSUM) monitoring extension of lbi_test().
# See docs/enhancements/monitoring.md, "Breitung & Diegel (2025) -- static
# LBI test AND sequential extension both done", for the full writeup.
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Weight normalization (eq. 12) ===\n")
Tm <- 500
for (cb in c(0, 1, 2, 5)) {
  w <- exuber:::bd_cusum_weights(Tm, cb)
  cat(sprintf(
    "c_bar=%.1f  sum(w^2)=%.6f (expect exactly 1 at c_bar=0, ~1 otherwise)  length=%d\n",
    cb, sum(w^2), length(w)
  ))
}

cat("\n=== 2. mCUSUM final-point statistic is formula-exact ===\n")
cat("(matches a hand-telescoped computation using training-window sigma_tilde)\n")
set.seed(1)
n <- 300
T_star <- 150
y <- cumsum(rnorm(n))
out <- monitor_lbi(y, r_star = T_star, c_bar = 0, level = 0.95)
final_stat <- out$stat[nrow(out$stat), 1]
dy <- diff(y)
sigma2_tilde <- mean(dy[seq_len(T_star - 1L)]^2)
manual <- (y[n] - y[T_star]) / sqrt(sigma2_tilde * (n - T_star))
cat(sprintf(
  "monitor_lbi: %.8f   manual telescoped: %.8f   |diff|: %.2e\n",
  final_stat, manual, abs(final_stat - manual)
))

cat("\n=== 3. Table 1 lookup (exact, both mCUSUM and wCUSUM share it) ===\n")
for (lv in c(0.90, 0.95, 0.975, 0.99, 0.995)) {
  cat(sprintf("level=%.3f -> b_alpha=%.2f\n", lv, exuber:::bd_cusum_q(lv)))
}
res <- tryCatch(exuber:::bd_cusum_q(0.80), error = function(e) "ERROR (expected)")
cat("untabulated level=0.80:", res, "\n")

cat("\n=== 4. Alarms never fire before the training window ends ===\n")
set.seed(3)
ok <- TRUE
for (i in 1:50) {
  y <- cumsum(rnorm(200))
  om <- monitor_lbi(y, r_star = 100, c_bar = 0, level = 0.95)
  if (!is.na(om$alarm) && om$alarm <= 100) ok <- FALSE
}
cat("all alarms strictly after T_star (50 reps):", ok, "\n")

cat("\n=== 5. Empirical false-alarm rate under H0 (pure random walk, 1,000 reps) ===\n")
set.seed(42)
nrep <- 1000
n <- 200
T_star <- 100
fa_mcusum <- fa_wcusum <- 0
for (i in seq_len(nrep)) {
  y <- cumsum(rnorm(n))
  om <- monitor_lbi(y, r_star = T_star, c_bar = 0, level = 0.95)
  ow <- monitor_lbi(y, r_star = T_star, c_bar = 2, level = 0.95)
  if (!is.na(om$alarm)) fa_mcusum <- fa_mcusum + 1
  if (!is.na(ow$alarm)) fa_wcusum <- fa_wcusum + 1
}
cat(sprintf("mCUSUM (c_bar=0) false-alarm rate: %.4f (nominal 0.05)\n", fa_mcusum / nrep))
cat(sprintf("wCUSUM (c_bar=2) false-alarm rate: %.4f (nominal 0.05)\n", fa_wcusum / nrep))

cat("\n=== 6. Detection power vs. monitor_cusum(type = 'standard') on the same DGP ===\n")
set.seed(4)
nrep <- 100
n <- 200
T_star <- 100
make_bubble_series <- function(n, T_star, bubble_start_frac = 0.65, rho = 1.03) {
  y <- numeric(n)
  y[seq_len(T_star)] <- cumsum(rnorm(T_star))
  bstart <- T_star + round((n - T_star) * (bubble_start_frac - T_star / n) / (1 - T_star / n))
  bstart <- max(bstart, T_star + 5)
  for (t in (T_star + 1):n) {
    y[t] <- if (t < bstart) y[t - 1] + rnorm(1) else rho * y[t - 1] + rnorm(1)
  }
  y
}

det_mcusum <- det_wcusum <- det_cusum_std <- 0
for (i in seq_len(nrep)) {
  y <- make_bubble_series(n, T_star)
  om <- monitor_lbi(y, r_star = T_star, c_bar = 0, level = 0.95)
  ow <- monitor_lbi(y, r_star = T_star, c_bar = 2, level = 0.95)
  oc <- monitor_cusum(y, r_star = T_star / n, b_alpha = 4.6)
  if (!is.na(om$alarm)) det_mcusum <- det_mcusum + 1
  if (!is.na(ow$alarm)) det_wcusum <- det_wcusum + 1
  if (!is.na(oc$alarm)) det_cusum_std <- det_cusum_std + 1
}
cat(sprintf("mCUSUM (c_bar=0) power: %.3f\n", det_mcusum / nrep))
cat(sprintf("wCUSUM (c_bar=2) power: %.3f\n", det_wcusum / nrep))
cat(sprintf("monitor_cusum(type='standard') power (same DGP): %.3f\n", det_cusum_std / nrep))

cat("\ndone\n")
