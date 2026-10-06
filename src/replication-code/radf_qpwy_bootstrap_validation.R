# Validation of monitor_quantile(boundary = "bootstrap"), the i.i.d. residual
# bootstrap of Algorithm 1 of Wu, Shi & Wu (2025) applied to the QPWY and
# QPSY paths. See docs/alternative-paradigms.md, "Implementation: bootstrap
# boundary".
#
# Sections 1-2 check the code against brute force. Section 3 measures the
# false-alarm rate under H0 of the bootstrap boundary and of the asymptotic
# boundary on the same series, and compares the bootstrap boundary with the
# exact finite-sample boundary (the 95% quantile of the path supremum under
# the true data generating process). Section 4 compares detection rates with
# the OLS monitors PWY and PSY.
#
# Run from the exuber-project/ root. Sections 3-4 use a worker cluster and take
# about two hours on 13 cores. PILOT=1 runs a reduced version in a few
# minutes. NCORES sets the number of workers (default 13).

devtools::load_all("exuber", quiet = TRUE)
pilot <- nzchar(Sys.getenv("PILOT"))
ncores <- as.integer(Sys.getenv("NCORES", "13"))

cat("=== 1. Bootstrap path maxima vs a brute-force rq() computation ===\n")
# Same seed, same resample, but every window statistic from quantreg::rq().
# quantile_boundary_boot() draws each replicate in a foreach loop with its own
# random-number stream, so the check calls the single-replicate draw directly.
rq_stat <- function(yy, tau) {
  m <- length(yy)
  ylag <- yy[1:(m - 1)]
  yresp <- yy[2:m]
  a <- unname(coef(quantreg::rq(yresp ~ ylag, tau = tau))["ylag"])
  f <- exuber:::quantile_check_density(yresp - ylag, tau)$f_hat
  (f / sqrt(tau * (1 - tau))) * sqrt(sum((ylag - mean(ylag))^2)) * (a - 1)
}
set.seed(9)
y <- cumsum(rt(45, df = 3))
n <- length(y)
minw <- 12
u <- diff(y)
u <- u - mean(u)
for (type in c("qpwy", "qpsy")) {
  set.seed(21)
  got <- vapply(1:4, function(i) exuber:::quantile_boot_max(u, 0.7, minw, type), numeric(1))
  set.seed(21)
  want <- vapply(1:4, function(i) {
    ys <- cumsum(c(0, sample(u, n - 1L, replace = TRUE)))
    max(vapply((minw + 1L):n, function(r) {
      if (type == "qpwy") rq_stat(ys[1:r], 0.7) else max(vapply(1:(r - minw), function(r1) rq_stat(ys[r1:r], 0.7), numeric(1)))
    }, numeric(1)))
  }, numeric(1))
  cat(sprintf("%s: max|diff| over 4 replicates = %.2e
", type, max(abs(got - want))))
}

cat("\n=== 2. The b = 100 burn-in of Algorithm 1 changes nothing ===\n")
# The statistic regresses with an intercept, so adding a constant to the
# series leaves every window statistic unchanged. With i.i.d. draws the last
# T values of a T + 100 resample have the same law as a fresh T draw.
set.seed(2)
ys <- cumsum(rnorm(60))
r_idx <- 15:60
shift <- exuber:::qpsy_stat_path(ys + 37.5, 0.8, r_idx, 14) - exuber:::qpsy_stat_path(ys, 0.8, r_idx, 14)
cat(sprintf("QPSY path, series shifted by 37.5: max|diff| = %.2e\n", max(abs(shift))))

# Worker cluster (sections 3 and 4) ------------------------------------------
cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {
  pkgload::load_all("exuber", compile = FALSE, quiet = TRUE)
  options(exuber.parallel = FALSE)
  NULL
})

gauss <- function(n) rnorm(n)
t3 <- function(n) rt(n, df = 3)
path_max <- function(y, tau, type) {
  n <- length(y)
  minw <- psy_minw(n)
  r_idx <- (minw + 1L):n
  max(if (type == "qpwy") {
    exuber:::qpwy_stat_path(y, tau, r_idx)
  } else {
    exuber:::qpsy_stat_path(y, tau, r_idx, minw)
  })
}
one_series <- function(i, type, n, innov, tau, B, seed0) {
  set.seed(seed0 + i)
  y <- cumsum(innov(n))
  minw <- psy_minw(n)
  dy <- diff(y)
  psi <- tau - as.numeric(dy < quantile(dy, tau, names = FALSE))
  delta <- max(min(cor(dy, psi), 1), -1)
  asym <- quantile(exuber:::quantile_boundary_sim(n, minw, 300, delta, type = type), 0.95, names = FALSE)
  boot <- quantile(exuber:::quantile_boundary_boot(y, tau, minw, B, type = type), 0.95, names = FALSE)
  c(stat = path_max(y, tau, type), asym = asym, boot = boot)
}
one_oracle <- function(k, type, n, innov, tau, seed0) {
  set.seed(seed0 + k)
  path_max(cumsum(innov(n)), tau, type)
}
parallel::clusterExport(cl, c("gauss", "t3", "path_max", "one_series", "one_oracle", "rq_stat"))

cat("\n=== 3. False-alarm rate under H0, bootstrap vs asymptotic boundary (nominal 5%) ===\n")
cat("exact = 95% quantile of the path supremum under the true DGP\n")
# type, n, reps, B, oracle reps, innovation label, innovation, tau
cfgs <- list(
  list("qpwy", 100, 200, 199, 1000, "gaussian", gauss, 0.5),
  list("qpwy", 100, 200, 199, 1000, "t3", t3, 0.5),
  list("qpwy", 100, 200, 199, 1000, "t3", t3, 0.2),
  list("qpwy", 100, 200, 199, 1000, "t3", t3, 0.8),
  list("qpwy", 100, 200, 199, 1000, "t3", t3, 0.9),
  list("qpwy", 100, 200, 199, 1000, "gaussian", gauss, 0.9),
  list("qpsy", 60, 120, 99, 500, "gaussian", gauss, 0.5),
  list("qpsy", 60, 120, 99, 500, "t3", t3, 0.8),
  list("qpsy", 60, 120, 99, 500, "gaussian", gauss, 0.9),
  list("qpsy", 60, 120, 99, 500, "t3", t3, 0.9)
)
if (pilot) cfgs <- lapply(cfgs[c(5, 9)], function(c) { c[[3]] <- 26; c[[4]] <- 49; c[[5]] <- 80; c })
size_rows <- list()
for (k in seq_along(cfgs)) {
  cfg <- cfgs[[k]]
  type <- cfg[[1]]; n <- cfg[[2]]; reps <- cfg[[3]]; B <- cfg[[4]]; K <- cfg[[5]]
  innov <- cfg[[7]]; tau <- cfg[[8]]
  t0 <- Sys.time()
  res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(reps), one_series, type = type, n = n, innov = innov, tau = tau, B = B, seed0 = 10000L * k))
  orc <- unlist(parallel::parLapplyLB(cl, seq_len(K), one_oracle, type = type, n = n, innov = innov, tau = tau, seed0 = 900000L + 10000L * k))
  exact <- quantile(orc, 0.95, names = FALSE)
  row <- data.frame(
    type = type, n = n, innov = cfg[[6]], tau = tau, reps = reps, B = B,
    size_asym = mean(res[, "stat"] > res[, "asym"]),
    size_boot = mean(res[, "stat"] > res[, "boot"]),
    size_exact = mean(res[, "stat"] > exact),
    cv_asym = mean(res[, "asym"]), cv_boot = mean(res[, "boot"]),
    cv_boot_sd = sd(res[, "boot"]), cv_exact = exact
  )
  size_rows[[k]] <- row
  cat(sprintf(
    "%s n=%d %s tau=%.1f (%d reps, B=%d): size asymptotic %.3f, bootstrap %.3f, exact-boundary %.3f | boundary asym %.2f, boot %.2f (sd %.2f), exact %.2f [%.0f min]\n",
    type, n, cfg[[6]], tau, reps, B, row$size_asym, row$size_boot, row$size_exact,
    row$cv_asym, row$cv_boot, row$cv_boot_sd, row$cv_exact, as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}

cat("\n=== 4. Detection of an end-of-sample bubble, bootstrap QPWY/QPSY vs OLS PWY/PSY ===\n")
# Table V design of the paper: bubble from 0.8 T with root 1 + 1 / T^0.6 to the
# sample end. Detection is the first alarm at or after the true origination,
# the PWY/PSY monitors use the flat sadf_cv/gsadf_cv of radf_mc_cv(). tau = 0.5
# and tau = 0.8 are fixed, not the optimal tau of the paper, and the boundary
# is flat, so the numbers are comparable in pattern to Table V and not in level.
sim_bubble <- function(n, innov) {
  te <- floor(0.8 * n)
  e <- innov(n)
  y <- numeric(n)
  y[1] <- e[1]
  for (t in 2:n) y[t] <- (if (t > te) 1 + 1 / n^0.6 else 1) * y[t - 1] + e[t]
  y
}
one_power <- function(i, n, innov, B_wy, B_sy, seed0, sadf_cv, gsadf_cv, do_qpsy) {
  set.seed(seed0 + i)
  y <- sim_bubble(n, innov)
  minw <- psy_minw(n)
  r_idx <- (minw + 1L):n
  first <- function(path, cv) {
    b <- which(path > cv)
    if (length(b)) minw + b[1] else NA_integer_
  }
  full <- radf(y)
  out <- c(pwy = first(full$badf, sadf_cv), psy = first(full$bsadf, gsadf_cv))
  for (tau in c(0.5, 0.8)) {
    cv <- quantile(exuber:::quantile_boundary_boot(y, tau, minw, B_wy, type = "qpwy"), 0.95, names = FALSE)
    out[paste0("qpwy", tau)] <- first(exuber:::qpwy_stat_path(y, tau, r_idx), cv)
    if (do_qpsy) {
      cv <- quantile(exuber:::quantile_boundary_boot(y, tau, minw, B_sy, type = "qpsy"), 0.95, names = FALSE)
      out[paste0("qpsy", tau)] <- first(exuber:::qpsy_stat_path(y, tau, r_idx, minw), cv)
    }
  }
  out
}
parallel::clusterExport(cl, c("sim_bubble", "one_power"))
# QPWY and PWY at n = 100; QPSY and PSY at n = 60, where the bootstrap is affordable
setups <- list(
  list(n = 100, reps = 200, B_wy = 199, B_sy = 0, do_qpsy = FALSE),
  list(n = 60, reps = 100, B_wy = 99, B_sy = 49, do_qpsy = TRUE)
)
if (pilot) setups <- lapply(setups, function(z) { z$reps <- 26; z$B_wy <- 49; z$B_sy <- 19; z })
for (st in setups) {
  cv_ols <- radf_mc_cv(st$n, nrep = 2000, seed = 1)
  for (cfg in list(list("gaussian", gauss, 1), list("t3", t3, 2))) {
    res <- do.call(rbind, parallel::parLapplyLB(
      cl, seq_len(st$reps), one_power,
      n = st$n, innov = cfg[[2]], B_wy = st$B_wy, B_sy = st$B_sy,
      seed0 = 70000L + 1000L * cfg[[3]], sadf_cv = cv_ols$sadf_cv[2], gsadf_cv = cv_ols$gsadf_cv[2],
      do_qpsy = st$do_qpsy
    ))
    te <- floor(0.8 * st$n)
    detected <- apply(res, 2, function(a) mean(!is.na(a) & a > te))
    early <- apply(res, 2, function(a) mean(!is.na(a) & a <= te))
    cat(sprintf("%s, n = %d, %d reps\n", cfg[[1]], st$n, st$reps))
    print(rbind(detected = round(detected, 3), `alarm before origination` = round(early, 3)))
  }
}

cat("\ndone\n")
