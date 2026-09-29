# Replication script for the 2026-08 innovation-generator DGP extensions:
# sim_vol_break(), sim_vol_garch() (+ TGARCH via gamma > 0), sim_vol_cir(),
# sim_vol_sv(), sim_fi(), sim_innov() -- R/sim.R. No prior validation
# script existed for these; written from scratch alongside the pyexuber
# port (2026-09-17).
#
# All of these are otherwise-deterministic transformations of Gaussian
# noise, so each is checked by fixing the underlying rnorm() draws (one
# call: set.seed(1); rnorm(10)) and hand-tracing the recursion -- exactly
# the technique pyexuber/tests/test_sim.py already uses for sim_psy1's
# regime-switching branches (a _FakeRNG feeding the same fixed pool to
# np.random.default_rng()'s call sites, in the same order R draws them).
devtools::load_all("exuber", quiet = TRUE)

set.seed(1)
eps10 <- rnorm(10)
cat("eps10:", paste(sprintf("%.8f", eps10), collapse = ","), "\n\n")

cat("=== sim_vol_break(n=10, tau=0.5, ratio=3, sigma=2) ===\n")
n <- 10; tau <- 0.5; ratio <- 3; sigma <- 2
sd_t <- sigma * ifelse(seq_len(n) > floor(tau * n), ratio, 1)
z_break <- sd_t * eps10
cat("sd_t:", sd_t, "\n")
cat("z:", paste(sprintf("%.8f", z_break), collapse = ","), "\n\n")

cat("=== sim_vol_garch(n=5, omega=0.1, alpha=0.1, beta=0.8, gamma=0) ===\n")
omega <- 0.1; alpha <- 0.1; beta <- 0.8
eps5 <- eps10[1:5]
h <- numeric(5); z <- numeric(5); h_prev <- 0; z_prev <- 0
for (t in 1:5) {
  h[t] <- omega + alpha * z_prev ^ 2 + beta * h_prev
  z[t] <- sqrt(h[t]) * eps5[t]
  h_prev <- h[t]; z_prev <- z[t]
}
cat("h:", paste(sprintf("%.8f", h), collapse = ","), "\n")
cat("z:", paste(sprintf("%.8f", z), collapse = ","), "\n\n")

cat("=== sim_vol_garch(..., gamma=0.2) [TGARCH] ===\n")
gamma2 <- 0.2
h2 <- numeric(5); z2 <- numeric(5); h_prev <- 0; z_prev <- 0
for (t in 1:5) {
  h2[t] <- omega + alpha * z_prev ^ 2 + beta * h_prev + gamma2 * z_prev ^ 2 * (z_prev < 0)
  z2[t] <- sqrt(h2[t]) * eps5[t]
  h_prev <- h2[t]; z_prev <- z2[t]
}
cat("h:", paste(sprintf("%.8f", h2), collapse = ","), "\n")
cat("z:", paste(sprintf("%.8f", z2), collapse = ","), "\n\n")

cat("=== sim_vol_cir(n=5, kappa=0.03, theta=0.25, xi=0.1) ===\n")
kappa <- 0.03; theta <- 0.25; xi <- 0.1; sigma0_sq <- theta
dt <- 1 / 5
db <- eps10[1:4] * sqrt(dt)
mult <- eps10[6:10]
sig2 <- numeric(5); sig2[1] <- sigma0_sq
for (i in 2:5) {
  prev <- max(sig2[i - 1], 0)
  sig2[i] <- max(prev + kappa * (theta - prev) * dt + xi * sqrt(prev) * db[i - 1], 0)
}
cat("sig2:", paste(sprintf("%.8f", sig2), collapse = ","), "\n")
cat("z:", paste(sprintf("%.8f", sqrt(sig2) * mult), collapse = ","), "\n\n")

cat("=== sim_vol_sv(n=5, phi=0.98, tau=0.1) ===\n")
phi <- 0.98; tau_sv <- 0.1
eta <- eps10[1:4] * tau_sv
mult_sv <- eps10[6:10]
log_sig2 <- numeric(5); log_sig2[1] <- 0
for (i in 2:5) log_sig2[i] <- phi * log_sig2[i - 1] + eta[i - 1]
cat("log_sig2:", paste(sprintf("%.8f", log_sig2), collapse = ","), "\n")
cat("z:", paste(sprintf("%.8f", exp(log_sig2 / 2) * mult_sv), collapse = ","), "\n\n")

cat("=== sim_fi(): psi recursion (d=0.2, m=5) + hand convolution (n=3) ===\n")
d <- 0.2; m <- 5
psi <- numeric(m + 1); psi[1] <- 1
for (j in 2:(m + 1)) psi[j] <- psi[j - 1] * (j - 2 + d) / (j - 1)
cat("psi:", paste(sprintf("%.8f", psi), collapse = ","), "\n")
epsfi <- eps10[1:8]
u <- vapply(1:3, function(k) sum(psi * rev(epsfi[k:(k + m)])), numeric(1))
cat("u:", paste(sprintf("%.8f", u), collapse = ","), "\n\n")

cat("=== sim_innov(dist='normal', sigma=2) ===\n")
cat("z:", paste(sprintf("%.8f", eps10 * 2), collapse = ","), "\n\n")

cat("=== sim_innov: t/skew_t closed-form constants (df=5) ===\n")
df <- 5
cat("t rescale 1/sqrt(df/(df-2)):", sprintf("%.10f", 1 / sqrt(df / (df - 2))), "\n")
xi <- -0.75
delta <- xi / sqrt(1 + xi ^ 2)
e_abs_t0 <- (2 * sqrt(df) / ((df - 1) * beta(df / 2, 0.5))) / sqrt(df / (df - 2))
mean_raw <- delta * e_abs_t0
sd_raw <- sqrt(max(1 - delta ^ 2 * e_abs_t0 ^ 2, .Machine$double.eps))
cat("delta:", sprintf("%.10f", delta), " e_abs_t0:", sprintf("%.10f", e_abs_t0),
    " mean_raw:", sprintf("%.10f", mean_raw), " sd_raw:", sprintf("%.10f", sd_raw), "\n")

# All numbers above are reproduced bit-for-bit in
# docs/replication/simulation-dgps/sim_vol_innovations_validation.py and
# in pyexuber/tests/test_sim.py (same reference values, duplicated as
# literals there per docs/replication/README.md's convention).
