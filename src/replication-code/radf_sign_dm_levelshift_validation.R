# Replicates Harvey, Leybourne, Tatlow & Zu (2025)'s Table 1, Case 1 row
# k = 2, mu = 5, T = 400, p = 0.8 (their worst PSY over-size case at this T):
# published PSY/sPSY/sbarPSY empirical sizes = 0.337 / 0.121 / 0.050.
#
# DGP (Section 5.1): alpha_n = 0.5, alpha_mu = 0, n_T = floor(k*sqrt(T))
# level shifts, n_T^+ = round(p*n_T) positive shifts of magnitude mu_T = mu,
# n_T^- = n_T - n_T^+ negative shifts of magnitude -mu_T. Shift locations
# t_i iid from floor(T*U[0,1])+1, no repeats. y_1 = eps_1, mu_0 = 0. pi =
# 0.1 trimming (minw = floor(0.1*T)), tests at nominal 0.05, critical
# values simulated under the null model with no shifts, eps_t ~ iid N(0,1).
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

set.seed(20260811)
Tn <- 400
k <- 2
mu <- 5
p <- 0.8
nrep <- 500
minw <- floor(0.1 * Tn)

n_T <- floor(k * sqrt(Tn))
n_pos <- round(p * n_T)
n_neg <- n_T - n_pos

simulate_path <- function() {
  eps <- rnorm(Tn)
  shift_locs <- sample.int(Tn, n_T, replace = FALSE)
  shift_sign <- c(rep(1, n_pos), rep(-1, n_neg))
  mu_t <- numeric(Tn)
  mu_t[shift_locs] <- shift_sign * mu
  level <- cumsum(mu_t)
  y <- cumsum(eps) + level
  y
}

cat("Simulating critical values under the null (no shifts) ...\n")
cv_mc <- radf_mc_cv(Tn, minw = minw, nrep = nrep, seed = 1)
cv_sign <- radf_sign_cv(Tn, minw = minw, nrep = nrep, seed = 2)
cv_sign_dm <- radf_sign_dm_cv(Tn, minw = minw, nrep = nrep, seed = 3)

cat("Simulating rejection rates under Case 1 (k=2, mu=5, p=0.8, T=400) ...\n")
reject_psy <- reject_spsy <- reject_sbarpsy <- logical(nrep)
for (r in seq_len(nrep)) {
  y <- simulate_path()
  reject_psy[r] <- radf(y, minw = minw)$gsadf > cv_mc$gsadf_cv["95%"]
  reject_spsy[r] <- radf_sign(y, minw = minw)$gsadf > cv_sign$gsadf_cv["95%"]
  reject_sbarpsy[r] <- radf_sign_dm(y, minw = minw)$gsadf > cv_sign_dm$gsadf_cv["95%"]
}

cat("\n--- Empirical size (nominal 0.05), Case 1: k=2, mu=5, p=0.8, T=400 ---\n")
cat(sprintf("PSY:      %.3f  (published: 0.337)\n", mean(reject_psy)))
cat(sprintf("sPSY:     %.3f  (published: 0.121)\n", mean(reject_spsy)))
cat(sprintf("sbarPSY:  %.3f  (published: 0.050)\n", mean(reject_sbarpsy)))
