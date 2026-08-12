# Validation of radf_monitor(..., boundary = "kurozumi", s0 = 0.4/0.8) --
# Kurozumi (2020)'s GSADF_{s0} monitoring detector, re-triaged and shipped
# after initially being scoped out as needing new recursion code.
# See docs/enhancements/monitoring.md, "Kurozumi (2020, 2021) -- SADF and
# GSADF cases both implemented", for the full writeup.
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Backward compatibility: s0 = 0 (default) unchanged ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
o_default <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = "kurozumi", level = 0.95)
o_explicit <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = "kurozumi", s0 = 0, level = 0.95)
cat("stat identical:", identical(o_default$stat, o_explicit$stat), "\n")
cat("boundary identical:", identical(o_default$boundary, o_explicit$boundary), "\n")

cat("\n=== 2. kurozumi_gsadf_stat() matches radf()$badf exactly at the s0 -> 0 limit ===\n")
minw <- 20
badf <- radf(y, minw = minw, lag = 0)$badf[, 1]
stat_limit <- exuber:::kurozumi_gsadf_stat(y, T_star = minw, s0 = 1 / minw)
cat("max|badf - stat(k1_max=1)|:", max(abs(unname(stat_limit) - badf)), "\n")

cat("\n=== 3. Formula-exact vs. brute-force lm() search over the restricted band ===\n")
T_star <- 75
s0 <- 0.4
k1_max <- floor(T_star * s0)
n <- length(y)
stat <- exuber:::kurozumi_gsadf_stat(y, T_star, s0)
for (k_check in c(5, 30, 60)) {
  t_check <- T_star + k_check
  brute <- sapply(seq_len(k1_max), function(k1) {
    yy <- y[k1:t_check]
    fit <- lm(diff(yy) ~ yy[-length(yy)])
    summary(fit)$coefficients[2, 3]
  })
  cat(sprintf(
    "k=%d  brute_max=%.6f  fast=%.6f  |diff|=%.2e\n",
    k_check, max(brute), stat[k_check], abs(max(brute) - stat[k_check])
  ))
}

cat("\n=== 4. Table 1 GSADF lookups (q04_df/q08_df columns) ===\n")
for (lv in c(0.90, 0.95, 0.99)) {
  cat(sprintf("level=%.2f s_bar=1 s0=0.4 -> %.4f\n", lv, exuber:::kurozumi_gsadf_q(lv, 1, 0.4)))
}
cat("level=0.95 s_bar=1 s0=0.8 ->", exuber:::kurozumi_gsadf_q(0.95, 1, 0.8), "(expect 2.3330)\n")
cat("level=0.95 s_bar=1 s0=0.6 (tie snap to 0.4) ->", exuber:::kurozumi_gsadf_q(0.95, 1, 0.6), "\n")

cat("\n=== 5. Alarms never fire before T*+1 ===\n")
set.seed(3)
ok <- TRUE
for (i in 1:30) {
  yy <- cumsum(rnorm(150))
  oo <- radf_monitor(yy, r_star = 0.5, boundary = "kurozumi", s0 = 0.4, level = 0.95)
  if (!is.na(oo$alarm) && oo$alarm <= 75) ok <- FALSE
}
cat("all alarms strictly after T_star (30 reps):", ok, "\n")

cat("\n=== 6. Empirical false-alarm rate under H0 ===\n")
set.seed(10)
nrep <- 300
n <- 150
T_star <- 75
fa_sadf <- fa_gsadf04 <- fa_gsadf08 <- 0
for (i in seq_len(nrep)) {
  yy <- cumsum(rnorm(n))
  o_sadf <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0, level = 0.95)
  o_g04 <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0.4, level = 0.95)
  o_g08 <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0.8, level = 0.95)
  if (!is.na(o_sadf$alarm)) fa_sadf <- fa_sadf + 1
  if (!is.na(o_g04$alarm)) fa_gsadf04 <- fa_gsadf04 + 1
  if (!is.na(o_g08$alarm)) fa_gsadf08 <- fa_gsadf08 + 1
}
cat(sprintf("SADF (s0=0)   FA rate: %.3f (nominal 0.05)\n", fa_sadf / nrep))
cat(sprintf("GSADF s0=0.4  FA rate: %.3f\n", fa_gsadf04 / nrep))
cat(sprintf("GSADF s0=0.8  FA rate: %.3f\n", fa_gsadf08 / nrep))

cat("\n=== 7. Detection power on a post-training-bubble DGP ===\n")
set.seed(20)
nrep <- 60
make_bubble_series <- function(n, T_star, bubble_start_frac = 0.65, rho = 1.03) {
  yy <- numeric(n)
  yy[1:T_star] <- cumsum(rnorm(T_star))
  bstart <- T_star + round((n - T_star) * (bubble_start_frac - T_star / n) / (1 - T_star / n))
  bstart <- max(bstart, T_star + 5)
  for (t in (T_star + 1):n) {
    yy[t] <- if (t < bstart) yy[t - 1] + rnorm(1) else rho * yy[t - 1] + rnorm(1)
  }
  yy
}
det_sadf <- det_g04 <- det_g08 <- 0
for (i in seq_len(nrep)) {
  yy <- make_bubble_series(n, T_star)
  o_sadf <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0, level = 0.95)
  o_g04 <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0.4, level = 0.95)
  o_g08 <- radf_monitor(yy, r_star = T_star, boundary = "kurozumi", s0 = 0.8, level = 0.95)
  if (!is.na(o_sadf$alarm)) det_sadf <- det_sadf + 1
  if (!is.na(o_g04$alarm)) det_g04 <- det_g04 + 1
  if (!is.na(o_g08$alarm)) det_g08 <- det_g08 + 1
}
cat(sprintf("SADF (s0=0)   power: %.3f\n", det_sadf / nrep))
cat(sprintf("GSADF s0=0.4  power: %.3f\n", det_g04 / nrep))
cat(sprintf("GSADF s0=0.8  power: %.3f\n", det_g08 / nrep))

cat("\ndone\n")
