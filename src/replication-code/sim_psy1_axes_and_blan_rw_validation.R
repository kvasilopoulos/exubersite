# Replication script for sim_psy1()'s optional e/shifts/coef_noise/coef_a
# axes and sim_blan(type = "rotermann_wilfling") -- R/sim.R. No prior
# validation script existed for these; written from scratch alongside the
# pyexuber port (2026-09-17).
devtools::load_all("exuber", quiet = TRUE)

cat("=== sim_psy1(e = ...): custom innovations replace rnorm(n-1, sd=sigma) ===\n")
set.seed(42)
y_plain <- sim_psy1(50, seed = 42)
y_zero_e <- sim_psy1(50, seed = 42, e = rep(0, 49))
cat("y_zero_e finite:", all(is.finite(y_zero_e)), "\n")
cat("y_zero_e != y_plain (e overrides rnorm):", !isTRUE(all.equal(y_zero_e, y_plain)), "\n\n")

cat("=== sim_psy1(shifts = ...): one-period deterministic bump ===\n")
y_shift <- sim_psy1(50, seed = 42, shifts = list(date = 10, size = 100))
cat("y_shift[10] - y_plain[10] (expect 100):", y_shift[10] - y_plain[10], "\n")
cat("y_shift[1:9] == y_plain[1:9] (no effect before shift date):",
    isTRUE(all.equal(y_shift[1:9], y_plain[1:9])), "\n\n")

cat("=== sim_blan(type = 'rotermann_wilfling') ===\n")
set.seed(7)
n_b <- 5; pi_b <- 0.7; delta_b <- 0.984; rw_sigma <- 0.05; b0 <- 0.1
theta_b <- rbinom(n_b - 1, 1, pi_b)
u_b <- rlnorm(n_b - 1, meanlog = -rw_sigma ^ 2 / 2, sdlog = rw_sigma)
cat("theta:", theta_b, "\n")
cat("u:", paste(sprintf("%.8f", u_b), collapse = ","), "\n")
b <- b0
for (i in 1:(n_b - 1)) {
  b[i + 1] <- if (theta_b[i] == 1) {
    b[i] * u_b[i] / delta_b
  } else {
    ((1 - pi_b * delta_b) / (1 - pi_b)) * b[i] * u_b[i]
  }
}
cat("b:", paste(sprintf("%.8f", b), collapse = ","), "\n\n")

cat("=== sim_blan(type = 'rotermann_wilfling') via the real function ===\n")
set.seed(7)
b_real <- sim_blan(n_b, pi = pi_b, type = "rotermann_wilfling", delta = delta_b,
                    rw_sigma = rw_sigma, b0 = b0)
cat("b (real fn call):", paste(sprintf("%.8f", b_real), collapse = ","), "\n")
cat("matches hand-traced b:", isTRUE(all.equal(as.numeric(b_real), b)), "\n")

# All numbers above are reproduced (the deterministic/formula parts
# bit-for-bit, the RNG-driven parts structurally) in
# docs/replication/simulation-dgps/sim_psy1_axes_and_blan_rw_validation.py
# and in pyexuber/tests/test_sim.py.
