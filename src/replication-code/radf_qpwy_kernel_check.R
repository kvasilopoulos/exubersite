# Kernel check for monitor_quantile(): the Gaussian kernel of
# quantile_check_density() against the Epanechnikov kernel that Wu, Shi & Wu
# (2025, Section 4) use for the density f(b_tau). The bandwidth is the same,
# h = 0.9 T^(-1/5) min(sd, IQR / 1.34), which is stats::bw.nrd0(). See
# docs/alternative-paradigms.md, "Implementation: bootstrap boundary".
#
# Section 1 compares the two density estimates on random-walk windows,
# section 2 the QPWY path maximum on the same series, and section 3 the size
# of the bootstrap test with each kernel (its own bootstrap boundary). The
# comparison is small (60 series, B = 49), so section 3 shows whether the
# conclusion survives the kernel and does not estimate a size.
#
# Run from the exuber-project/ root. It takes about 30 minutes on 4 workers.

suppressMessages(devtools::load_all("exuber", quiet = TRUE))
ns <- asNamespace("exuber")
gauss_den <- get("quantile_check_density", ns)
# Epanechnikov kernel, h = 0.9 T^(-1/5) min(sd, IQR/1.34) (WSW Section 4); bw.nrd0 is the same h
epan_den <- function(u, tau) {
  b <- quantile(u, probs = tau, names = FALSE)
  h <- stats::bw.nrd0(u)
  z <- (b - u) / h
  list(b_tau = b, f_hat = mean(ifelse(abs(z) <= 1, 0.75 * (1 - z^2), 0)) / h)
}
set_den <- function(f) { unlockBinding("quantile_check_density", ns); assign("quantile_check_density", f, ns) }
pmax_stat <- function(y, tau) max(exuber:::qpwy_stat_path(y, tau, (psy_minw(length(y)) + 1L):length(y)))

cat("--- 2. path maximum, same series, both kernels (QPWY n = 100) ---\n")
for (dist in c("gaussian", "t3")) for (tau in c(0.5, 0.9)) {
  set.seed(2)
  m <- t(replicate(60, { y <- cumsum(if (dist == "t3") rt(100, 3) else rnorm(100)); set_den(gauss_den); g <- pmax_stat(y, tau); set_den(epan_den); e <- pmax_stat(y, tau); c(g, e) }))
  cat(sprintf("%s tau=%.1f: mean max gauss %.3f, epan %.3f, cor %.3f\n", dist, tau, mean(m[, 1]), mean(m[, 2]), cor(m[, 1], m[, 2])))
}
set_den(gauss_den)

cat("--- 3. size, QPWY n = 100, 60 series, B = 49: Gaussian vs Epanechnikov kernel ---
")
one <- function(i, dist, tau) {
  ns <- asNamespace("exuber")
  set.seed(500 + i)
  y <- cumsum(if (dist == "t3") rt(100, 3) else rnorm(100))
  minw <- psy_minw(100)
  out <- c()
  for (k in c("gauss", "epan")) {
    f <- if (k == "gauss") gauss_den else epan_den
    unlockBinding("quantile_check_density", ns); assign("quantile_check_density", f, ns)
    cv <- quantile(exuber:::quantile_boundary_boot(y, tau, minw, 49, seed = i), 0.95, names = FALSE)
    out <- c(out, setNames(c(pmax_stat(y, tau), cv), paste0(k, c("_stat", "_cv"))))
  }
  out
}
cl <- parallel::makeCluster(4)
parallel::clusterExport(cl, c("gauss_den", "epan_den", "pmax_stat", "one"))
invisible(parallel::clusterEvalQ(cl, { suppressMessages(pkgload::load_all("exuber", compile = FALSE, quiet = TRUE)); NULL }))
for (cfg in list(c("t3", 0.9), c("gaussian", 0.5))) {
  r <- do.call(rbind, parallel::parLapplyLB(cl, 1:60, one, dist = cfg[1], tau = as.numeric(cfg[2])))
  cat(sprintf("%s tau=%s: size gauss %.3f (cv %.2f), epan %.3f (cv %.2f)
", cfg[1], cfg[2],
    mean(r[, "gauss_stat"] > r[, "gauss_cv"]), mean(r[, "gauss_cv"]), mean(r[, "epan_stat"] > r[, "epan_cv"]), mean(r[, "epan_cv"])))
}
parallel::stopCluster(cl)
