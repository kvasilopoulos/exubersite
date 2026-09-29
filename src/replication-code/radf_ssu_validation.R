# Validation of ssu_test() and cusum_test() -- Kurozumi & Nishi (2025)'s
# SSU/GSSU stochastic-unit-root bubble tests, their UR/GUR union with
# SADF/GSADF, and the CS/GCS/CSSQ/GCSSQ CUSUM-type statistics.
# See docs/volatility-robustness.md, "Stochastic explosive\n# -coefficient test", for the full writeup.
#
# Run from the exuber-project/ root.

devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Formula-exact: t^{omega,c} vs. brute-force lm() + manual cross-moment ===\n")
set.seed(2)
n <- 150
y <- cumsum(rnorm(n))
ps <- exuber:::ssu_prefix_sums(y)

brute_force_stat <- function(hi, lo = 0) {
  win <- (lo + 1):hi
  x1 <- y[win]
  d1 <- y[win + 1] - y[win]
  x2 <- x1^2
  d2 <- d1^2

  fit6 <- lm(d1 ~ x1)
  fit7 <- lm(d2 ~ x2)
  eps_hat <- residuals(fit6)
  eta_hat <- residuals(fit7)
  L <- length(win)
  sigma2_eps <- sum(eps_hat^2) / (L - 2)
  sigma2_eta <- sum(eta_hat^2) / (L - 2)
  sigma2_epseta <- sum(eps_hat * eta_hat) / (L - 2)
  sigma_eps <- sqrt(sigma2_eps)
  sigma_eta <- sqrt(sigma2_eta)
  psi_hat <- sigma2_epseta / (sigma_eps * sigma_eta)

  omega_hat <- unname(coef(fit7)[2])
  Sxx2_c <- sum((x2 - mean(x2))^2)
  t_omega <- omega_hat / sqrt(sigma2_eta / Sxx2_c)

  num_corr <- sum((x2 - mean(x2)) * d1)
  den_corr <- sqrt(Sxx2_c)
  correction <- (psi_hat / sigma_eps) * num_corr / den_corr
  (t_omega - correction) / sqrt(1 - psi_hat^2)
}

for (hi_check in c(50, 80, 120, 149)) {
  fast <- exuber:::ssu_stat_path(ps, hi_check)
  manual <- brute_force_stat(hi_check)
  cat(sprintf("hi=%d fast=%.8f manual=%.8f |diff|=%.2e\n", hi_check, fast, manual, abs(fast - manual)))
}

cat("\n=== 2. Table I lookup ===\n")
for (lv in c(90, 95, 99)) {
  cat(sprintf("sig_lvl=%g -> crit=%.2f\n", lv, exuber:::ssu_q(lv)))
}
res <- tryCatch(exuber:::ssu_q(80), error = function(e) "ERROR (expected)")
cat("untabulated sig_lvl=80:", res, "\n")

cat("\n=== 3. minw matches psy_minw() (SSU's own r0 formula) ===\n")
set.seed(1)
y2 <- cumsum(rnorm(120))
out <- ssu_test(y2)
cat("minw:", attr(out, "minw"), " psy_minw(120):", psy_minw(120), "\n")

cat("\n=== 4. Empirical false-alarm rate under H0 ===\n")
set.seed(1)
nrep <- 300
n <- 200
fa_10 <- fa_05 <- fa_01 <- 0
for (i in seq_len(nrep)) {
  yy <- cumsum(rnorm(n))
  out10 <- ssu_test(yy, sig_lvl = 90)
  out05 <- ssu_test(yy, sig_lvl = 95)
  out01 <- ssu_test(yy, sig_lvl = 99)
  if (out10$detected) fa_10 <- fa_10 + 1
  if (out05$detected) fa_05 <- fa_05 + 1
  if (out01$detected) fa_01 <- fa_01 + 1
}
cat(sprintf("sig_lvl=90%%  FA rate: %.3f (nominal 0.10)\n", fa_10 / nrep))
cat(sprintf("sig_lvl=95%%  FA rate: %.3f (nominal 0.05)\n", fa_05 / nrep))
cat(sprintf("sig_lvl=99%%  FA rate: %.3f (nominal 0.01)\n", fa_01 / nrep))

cat("\n=== 5. Detection power: stochastic-explosive-coefficient DGP (KN's own eq. 2 style) ===\n")
set.seed(2)
nrep <- 60
n <- 200
make_stochastic_bubble <- function(n, te_frac = 0.5, c1 = 3, a = 4) {
  yy <- numeric(n)
  yy[1] <- rnorm(1)
  Te <- round(te_frac * n)
  for (t in 2:n) {
    if (t <= Te) {
      yy[t] <- yy[t - 1] + rnorm(1)
    } else {
      rho_t <- 1 + c1 / n + a * rnorm(1) / sqrt(n)
      yy[t] <- rho_t * yy[t - 1] + rnorm(1)
    }
  }
  yy
}
det_ssu <- 0
for (i in seq_len(nrep)) {
  yy <- make_stochastic_bubble(n)
  out <- ssu_test(yy, sig_lvl = 95)
  if (out$detected) det_ssu <- det_ssu + 1
}
cat(sprintf("SSU power on stochastic-coefficient bubble DGP: %.3f\n", det_ssu / nrep))

cat("\n=== 6. Detection power on a deterministic explosive DGP, vs. standard SADF ===\n")
set.seed(3)
det_ssu2 <- det_sadf <- 0
for (i in seq_len(nrep)) {
  n1 <- 100
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.03^(1:100) + cumsum(rnorm(100, sd = 1))
  yy <- c(normal_part, expl_part)
  out <- ssu_test(yy, sig_lvl = 95)
  r <- radf(yy)
  cv <- radf_mc_cv(length(yy))
  if (out$detected) det_ssu2 <- det_ssu2 + 1
  if (r$sadf > cv$sadf_cv[2]) det_sadf <- det_sadf + 1
}
cat(sprintf("SSU power: %.3f  SADF power: %.3f (deterministic explosive DGP)\n", det_ssu2 / nrep, det_sadf / nrep))

cat("\n=== 7. GSSU: windows off the sample start, and the sup over starts ===\n")
d_lo <- max(vapply(c(1, 17, 60), function(l) abs(exuber:::ssu_stat_path(ps, 120, l) - brute_force_stat(120, l)), 0))
g <- exuber:::gssu_stat_path(ps, 60:70, 25)
bg <- vapply(60:70, function(h) max(vapply(0:(h - 25), function(l) brute_force_stat(h, l), 0)), 0)
cat(sprintf("max|diff| lo > 0: %.2e   gssu path max|diff|: %.2e\n", d_lo, max(abs(g - bg))))

cat("\n=== 8. CUSUM family: running max/min vs brute-force double loop ===\n")
set.seed(7)
yc <- cumsum(rnorm(60)) + c(rep(0, 40), cumsum(1.1^(1:20)))
dd <- diff(yc)
nd <- length(dd)
worst <- 0
for (type in c("cs", "gcs", "cssq", "gcssq")) {
  vals <- unlist(lapply(1:nd, function(k) {
    vapply(if (type %in% c("cs", "cssq")) 0 else 0:(k - 1), function(j) {
      w <- dd[(j + 1):k]
      if (type %in% c("cs", "gcs")) {
        sum(w) / sqrt(mean(dd^2) * nd)
      } else {
        (sum(w^2) - (k - j) / nd * sum(dd^2)) / sqrt((mean(dd^4) - mean(dd^2)^2) * nd)
      }
    }, 0)
  }))
  out <- cusum_test(yc, type = type)
  worst <- max(worst, abs(out$sup - max(vals)), if (!is.null(out$inf)) abs(out$inf - min(vals)) else 0)
}
cat(sprintf("max|diff| over cs/gcs/cssq/gcssq sup and inf: %.2e\n", worst))

cat("\n=== 9. Empirical size at 5%, n = 200, 300 reps, all statistics ===\n")
n <- 200
cv200 <- radf_mc_cv(n, nrep = 2000, seed = 1)
decide <- function(yy) {
  c(
    SSU = unname(ssu_test(yy)$detected),
    GSSU = unname(ssu_test(yy, type = "gssu")$detected),
    UR = unname(ssu_test(yy, union = TRUE, cv = cv200)$union_detected),
    GUR = unname(ssu_test(yy, type = "gssu", union = TRUE, cv = cv200)$union_detected),
    CS = unname(cusum_test(yy, type = "cs")$detected),
    GCS = unname(cusum_test(yy, type = "gcs")$detected),
    CSSQ = unname(cusum_test(yy, type = "cssq")$detected),
    GCSSQ = unname(cusum_test(yy, type = "gcssq")$detected)
  )
}
size <- rowMeans(vapply(1:300, function(i) {
  set.seed(9000 + i)
  decide(cumsum(rnorm(n)))
}, logical(8)))
print(round(size, 3))

cat("\n=== 10. Power at 5%, n = 200, 100 reps: stochastic (a = 4) vs deterministic (a = 0) coefficient ===\n")
kn_dgp <- function(n, a, c1 = 3, te_frac = 0.5) {
  yy <- numeric(n)
  for (t in 2:n) {
    rho <- if (t > te_frac * n) 1 + c1 / n + a * rnorm(1) / sqrt(n) else 1
    yy[t] <- rho * yy[t - 1] + rnorm(1)
  }
  yy
}
for (a in c(4, 0)) {
  pw <- rowMeans(vapply(1:100, function(i) {
    set.seed(9500 + i)
    decide(kn_dgp(n, a, c1 = if (a == 0) 10 else 3))
  }, logical(8)))
  cat(sprintf("a = %g:\n", a))
  print(round(pw, 3))
}

cat("\ndone\n")
