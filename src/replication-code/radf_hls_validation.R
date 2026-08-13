# Validation script for dating_hls() (Harvey, Leybourne & Sollis 2017,
# "Improving the accuracy of asset price bubble start and end date
# estimators"). See docs/enhancements/dating-and-root-inference.md,
# "SSR/BIC dating vs. PSY recursive dating", for the full write-up. Run
# from the exuber/ package root (or adjust the devtools::load_all() path
# below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Formula-exact: hls_segment_ssr() vs brute-force lm() ===\n")
set.seed(1)
y <- cumsum(rnorm(40))
ps <- exuber:::hls_prefix_sums(y)
n1 <- length(y) - 1L
x_all <- y[1:n1]; z_all <- y[2:(n1 + 1)] - y[1:n1]
for (seg in list(c(1, 10), c(11, 25), c(5, 39))) {
  lo <- seg[1] - 1; hi <- seg[2]
  manual <- exuber:::hls_segment_ssr(ps, lo, hi, TRUE)
  idx <- (lo + 1):hi
  brute <- sum(resid(lm(z_all[idx] ~ x_all[idx]))^2)
  cat(sprintf("segment (%d,%d]: manual=%.6f brute=%.6f match=%s\n",
              lo, hi, manual, brute, isTRUE(all.equal(manual, brute, tolerance = 1e-8))))
}

cat("\n=== 2. Full Model 4 (3-breakpoint) grid search vs brute-force nested lm() ===\n")
set.seed(5)
n3 <- 24
y3 <- cumsum(rnorm(n3))
ps3 <- exuber:::hls_prefix_sums(y3)
m4 <- exuber:::hls_model4(y3, ps3, trim = 0.1)
n1c <- n3 - 1L
x3 <- y3[1:n1c]; z3 <- y3[2:(n1c + 1)] - y3[1:n1c]
k_min3 <- max(2L, ceiling(0.1 * n1c))
best4 <- list(ssr = Inf)
for (tau1 in k_min3:(n1c - 3 * k_min3)) {
  for (tau2 in (tau1 + k_min3):(n1c - 2 * k_min3)) {
    if (y3[tau2 + 1] <= y3[tau1 + 1]) next
    for (tau3 in (tau2 + k_min3):(n1c - k_min3)) {
      if (y3[tau2 + 1] <= y3[tau3 + 1]) next
      ssr <- sum(z3[1:tau1]^2) +
        sum(resid(lm(z3[(tau1 + 1):tau2] ~ x3[(tau1 + 1):tau2]))^2) +
        sum(resid(lm(z3[(tau2 + 1):tau3] ~ x3[(tau2 + 1):tau3]))^2) +
        sum(z3[(tau3 + 1):n1c]^2)
      if (ssr < best4$ssr) best4 <- list(tau1 = tau1, tau2 = tau2, tau3 = tau3, ssr = ssr)
    }
  }
}
cat("vectorized:", m4$tau1, m4$tau2, m4$tau3, m4$ssr, "\n")
cat("brute force:", best4$tau1, best4$tau2, best4$tau3, best4$ssr, "\n")

cat("\n=== 3. Performance at realistic sample sizes ===\n")
for (n in c(100, 200, 400)) {
  set.seed(1)
  yn <- cumsum(rnorm(n))
  t0 <- Sys.time()
  dating_hls(yn, trim = 0.05)
  cat(sprintf("n=%d: %.2f sec\n", n, as.numeric(Sys.time() - t0)))
}

cat("\n=== 4. Monte Carlo: model-selection accuracy and breakpoint bias by DGP ===\n")
cat("(bubble/collapse regimes built on a large positive base (100) so the\n")
cat("explosive signal is a genuinely large absolute-scale departure from\n")
cat("noise -- an early attempt using a small mean-zero random-walk base\n")
cat("gave a systematic Model-1 selection bias not because of a code bug\n")
cat("but because that DGP's signal-to-noise ratio was too weak; this was\n")
cat("found and fixed during this validation.)\n\n")

sim_model4 <- function(seed, n1 = 60, n2 = 25, n3 = 25, n4 = 40, base = 100, c_bubble = 1.05) {
  set.seed(seed)
  unit1 <- base + cumsum(rnorm(n1))
  bubble <- unit1[n1] * c_bubble^(1:n2) + cumsum(rnorm(n2))
  target <- bubble[n2] * 0.5
  collapse <- numeric(n3)
  collapse[1] <- bubble[n2] + rnorm(1)
  for (k in 2:n3) collapse[k] <- target + 0.85 * (collapse[k - 1] - target) + rnorm(1)
  recovery <- collapse[n3] + cumsum(rnorm(n4))
  list(y = c(unit1, bubble, collapse, recovery), true_tau1 = n1, true_tau2 = n1 + n2)
}
run4 <- function(seed) {
  sim <- sim_model4(seed)
  out <- dating_hls(sim$y, trim = 0.05)
  list(model = out$model[["series1"]],
       orig_bias = as.numeric(out$origination[["series1"]]) - sim$true_tau1,
       coll_bias = if (!is.na(out$collapse[["series1"]])) as.numeric(out$collapse[["series1"]]) - sim$true_tau2 else NA)
}
res4 <- lapply(1:30, run4)
models4 <- sapply(res4, `[[`, "model")
cat("Model 4 DGP -- selection freq:", paste(sapply(1:4, function(m) sprintf("M%d=%.2f", m, mean(models4 == m))), collapse = ", "), "\n")
cat(sprintf("  origination mean|bias|=%.2f; collapse mean|bias| (M3/4 only)=%.2f\n",
            mean(abs(sapply(res4, `[[`, "orig_bias"))),
            mean(abs(sapply(res4[models4 %in% c(3, 4)], `[[`, "coll_bias")))))

sim_model2 <- function(seed, n1 = 60, n2 = 30, n3 = 60, base = 100, c_bubble = 1.05) {
  set.seed(seed)
  unit1 <- base + cumsum(rnorm(n1))
  bubble <- unit1[n1] * c_bubble^(1:n2) + cumsum(rnorm(n2))
  unit2 <- bubble[n2] + cumsum(rnorm(n3))
  list(y = c(unit1, bubble, unit2), true_tau1 = n1)
}
run2 <- function(seed) {
  sim <- sim_model2(seed)
  out <- dating_hls(sim$y, trim = 0.05)
  list(model = out$model[["series1"]], orig_bias = as.numeric(out$origination[["series1"]]) - sim$true_tau1)
}
res2 <- lapply(1:30, run2)
models2 <- sapply(res2, `[[`, "model")
cat("Model 2 DGP -- selection freq:", paste(sapply(1:4, function(m) sprintf("M%d=%.2f", m, mean(models2 == m))), collapse = ", "),
    sprintf("; origination mean|bias|=%.2f\n", mean(abs(sapply(res2, `[[`, "orig_bias")))))

sim_model1 <- function(seed, n1 = 80, n2 = 60, base = 100, c_bubble = 1.05) {
  set.seed(seed)
  unit1 <- base + cumsum(rnorm(n1))
  bubble <- unit1[n1] * c_bubble^(1:n2) + cumsum(rnorm(n2))
  list(y = c(unit1, bubble), true_tau1 = n1)
}
run1 <- function(seed) {
  sim <- sim_model1(seed)
  out <- dating_hls(sim$y, trim = 0.05)
  list(model = out$model[["series1"]], orig_bias = as.numeric(out$origination[["series1"]]) - sim$true_tau1)
}
res1 <- lapply(1:30, run1)
models1 <- sapply(res1, `[[`, "model")
cat("Model 1 DGP -- selection freq:", paste(sapply(1:4, function(m) sprintf("M%d=%.2f", m, mean(models1 == m))), collapse = ", "),
    sprintf("; origination mean|bias|=%.2f\n", mean(abs(sapply(res1, `[[`, "orig_bias")))))

run_h0 <- function(seed) {
  set.seed(seed)
  y0 <- 100 + cumsum(rnorm(150))
  dating_hls(y0, trim = 0.05)$model[["series1"]]
}
models0 <- sapply(1:30, run_h0)
cat("Pure H0 (no bubble) -- selection freq:", paste(sapply(1:4, function(m) sprintf("M%d=%.2f", m, mean(models0 == m))), collapse = ", "), "\n")
