# Validation of monitor_quantile() -- Wu, Shi & Wu (2025)'s QPWY and QPSY
# recursive quantile monitoring strategies. See docs/
# alternative-paradigms.md, "Quantile-based detection", for the full
# writeup -- including two real bugs found via Monte Carlo size checks:
# (a) a per-r marginal quantile used as boundary instead of a
# supremum-calibrated one (2026-08-11), and (b) the limiting process's
# independent-BM component Z simulated as ONE z per replicate instead of a
# process over windows (2026-09-29).
#
# Run from the exuber-project/ root. Sections 4-5 take ~30 minutes.

devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Point statistics vs quantreg::rq() brute force ===\n")
set.seed(3)
y <- cumsum(rnorm(60))
rq_stat <- function(yy, tau) {
  m <- length(yy)
  ylag <- yy[1:(m - 1)]
  yresp <- yy[2:m]
  a <- unname(coef(quantreg::rq(yresp ~ ylag, tau = tau))["ylag"])
  f <- exuber:::quantile_check_density(yresp - ylag, tau)$f_hat
  (f / sqrt(tau * (1 - tau))) * sqrt(sum((ylag - mean(ylag))^2)) * (a - 1)
}
d1 <- abs(exuber:::qpwy_stat_path(y, 0.5, 60) - rq_stat(y, 0.5))
minw <- 15
r_idx <- (minw + 1):60
qpsy <- exuber:::qpsy_stat_path(y, 0.7, r_idx, minw)
brute <- vapply(r_idx, function(r) max(vapply(1:(r - minw), function(r1) rq_stat(y[r1:r], 0.7), 0)), 0)
cat(sprintf("QPWY |diff| = %.2e   QPSY max|diff| = %.2e\n", d1, max(abs(qpsy - brute))))

cat("\n=== 2. quantile_boundary_sim vs per-window brute force ===\n")
n <- 30
minw <- 8
delta <- c(0.3, 0.9)
worst <- 0
for (type in c("qpwy", "qpsy")) {
  sim <- exuber:::quantile_boundary_sim(n, minw, 3, delta, type = type, seed = 11)
  set.seed(11)
  for (i in 1:3) {
    e <- rnorm(n - 1)
    v <- rnorm(n - 1)
    x <- c(0, cumsum(e))[1:(n - 1)]
    U <- NULL
    for (hi in minw:(n - 1)) {
      for (lo in if (type == "qpwy") 0 else 0:(hi - minw)) {
        k <- (lo + 1):hi
        xb <- x[k] - mean(x[k])
        U <- rbind(U, (delta * sum(xb * e[k]) + sqrt(1 - delta^2) * sum(xb * v[k])) / sqrt(sum(xb^2)))
      }
    }
    worst <- max(worst, abs(sim[i, ] - apply(U, 2, max)))
  }
}
cat(sprintf("max|diff| over qpwy/qpsy, 3 reps, 2 deltas = %.2e\n", worst))

cat("\n=== 3. THE 2026-09-29 BUG: one z per replicate vs Z as a process ===\n")
# QPWY slice (lo = 0), n = 200, psy_minw, 4000 reps: 95% boundary under the
# old single-z construction, and that boundary's true size under the
# correct limiting process
n <- 200
minw <- psy_minw(n)
nrep <- 4000
set.seed(1)
hi <- minw:(n - 1)
Qm <- Zm <- matrix(NA_real_, nrep, length(hi))
for (i in seq_len(nrep)) {
  e <- rnorm(n - 1)
  v <- rnorm(n - 1)
  x <- c(0, cumsum(e))[1:(n - 1)]
  cs <- function(w) cumsum(w)[hi]
  sxx <- cs(x^2) - cs(x)^2 / hi
  Qm[i, ] <- (cs(x * e) - cs(x) * cs(e) / hi) / sqrt(sxx)
  Zm[i, ] <- (cs(x * v) - cs(x) * cs(v) / hi) / sqrt(sxx)
}
z1 <- rnorm(nrep)
for (d in c(0.8, 0.5, 0.2)) {
  right <- apply(d * Qm + sqrt(1 - d^2) * Zm, 1, max)
  single <- apply(d * Qm + sqrt(1 - d^2) * z1, 1, max)
  b_single <- quantile(single, 0.95)
  cat(sprintf(
    "delta = %.1f: single-z boundary %.3f, correct %.3f, true size of single-z boundary %.3f\n",
    d, b_single, quantile(right, 0.95), mean(right > b_single)
  ))
}

cat("\n=== 4. Empirical false-alarm rate under H0 (nominal 5%) ===\n")
fa_rate <- function(type, n, reps, innov, tau) {
  mean(vapply(seq_len(reps), function(i) {
    set.seed(5000 + i)
    yy <- cumsum(innov(n))
    !is.na(suppressMessages(monitor_quantile(yy, tau = tau, nrep = 300, seed = i, type = type))$alarm)
  }, logical(1)))
}
gauss <- function(n) rnorm(n)
t3 <- function(n) rt(n, df = 3)
for (cfg in list(
  list("qpwy", 150, 200, "gaussian", gauss, 0.5),
  list("qpwy", 150, 200, "t3", t3, 0.5),
  list("qpwy", 150, 200, "t3", t3, 0.9),
  list("qpwy", 150, 200, "t3", t3, 0.8),
  list("qpwy", 150, 200, "t3", t3, 0.2),
  list("qpsy", 100, 100, "gaussian", gauss, 0.5),
  list("qpsy", 100, 100, "t3", t3, 0.9),
  list("qpsy", 100, 80, "t3", t3, 0.5),
  list("qpsy", 100, 80, "t3", t3, 0.8),
  list("qpsy", 100, 80, "gaussian", gauss, 0.9)
)) {
  cat(sprintf(
    "%s n=%d reps=%d %s tau=%.1f: %.3f\n", cfg[[1]], cfg[[2]], cfg[[3]], cfg[[4]], cfg[[6]],
    fa_rate(cfg[[1]], cfg[[2]], cfg[[3]], cfg[[5]], cfg[[6]])
  ))
}

cat("\n=== 5. Detection power, explosive from t = 71 of n = 100, vs SADF ===\n")
reps <- 60
det <- c(qpwy = 0, qpsy = 0, sadf = 0)
cv <- radf_mc_cv(100, nrep = 2000, seed = 1)
for (i in seq_len(reps)) {
  set.seed(3000 + i)
  normal_part <- cumsum(rnorm(70))
  yy <- c(normal_part, normal_part[70] * 1.04^(1:30) + cumsum(rnorm(30)))
  det["qpwy"] <- det["qpwy"] + !is.na(monitor_quantile(yy, nrep = 300, seed = i)$alarm)
  det["qpsy"] <- det["qpsy"] + !is.na(monitor_quantile(yy, nrep = 300, seed = i, type = "qpsy")$alarm)
  det["sadf"] <- det["sadf"] + (radf(yy)$sadf > cv$sadf_cv[2])
}
print(round(det / reps, 3))

cat("\ndone\n")
