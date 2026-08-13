# Validation of contagion_reg() -- Greenaway-McGrevy & Phillips (2016)'s
# bubble contagion regression, minimum-viable subset (fixed-window AR(1)
# sequence, single-delay Nadaraya-Watson regression, LOOCV bandwidth).
# See docs/enhancements/multivariate.md, "Contagion regression",
# "Implementation -- done (2026-08-10), minimum-viable subset", for the
# full writeup, including two real bugs found and fixed here (a
# window-width off-by-one and a matrix-orientation bug in the LOOCV SSE
# helper).
#
# No published numeric table exists to validate against -- the source
# paper's own results are Figures 7-8, not tabulated numbers. Validated
# instead via brute-force cross-checks of each closed-form piece, plus a
# directional sensible-behavior check.
#
# Run from the exuber/ package root.

devtools::load_all(".")

cat("=== 1. Fixed-window AR(1) coefficient sequence (eq. 1) vs. brute-force lm() ===\n")
set.seed(1)
n <- 150
S <- 50
core <- cumsum(rnorm(n))
beta_core <- exuber:::contagion_fixed_window_beta(core, S)
for (t_check in c(60, 80, 100, 130, 150)) {
  win <- (t_check - S + 1):t_check
  fit <- lm(core[win[-1]] ~ core[win[-length(win)]])
  cf <- beta_core[as.character(t_check)]
  cat(sprintf(
    "t=%d closed-form=%.8f lm()=%.8f |diff|=%.2e\n",
    t_check, cf, coef(fit)[2], abs(cf - coef(fit)[2])
  ))
}

cat("\n=== 2. Nadaraya-Watson ratio (eq. 6) vs. manual weighted-least-squares ===\n")
y <- 0.5 * core + cumsum(rnorm(n, sd = 0.5))
bc <- exuber:::contagion_fixed_window_beta(core, S)
bj <- exuber:::contagion_fixed_window_beta(y, S)
d <- 2
r_test <- 0.5
h_test <- 0.2
fast <- exuber:::contagion_nw_delta2(bc, bj, n, r_test, h_test, d)
s <- as.integer(names(bj))
core_shift <- (bc - mean(bc))[as.character(s - d)]
valid <- !is.na(core_shift)
s2 <- s[valid]
bjc <- (bj - mean(bj))[valid]
csh <- core_shift[valid]
w <- dnorm((s2 / n - r_test) / h_test) / h_test
manual <- sum(w * bjc * csh) / sum(w * csh^2)
cat(sprintf("fast=%.10f manual=%.10f |diff|=%.2e\n", fast, manual, abs(fast - manual)))

cat("\n=== 3. LOOCV SSE (eq. 7) vs. manual leave-one-out double loop ===\n")
fast_sse <- exuber:::contagion_loocv_sse(h_test, bc, bj, n, d)
m <- length(s2)
manual_sse <- 0
for (i in seq_len(m)) {
  r_i <- s2[i] / m
  num <- 0
  den <- 0
  for (p in seq_len(m)) {
    if (p == i) next
    wp <- dnorm((s2[p] / n - r_i) / h_test) / h_test
    num <- num + wp * bjc[p] * csh[p]
    den <- den + wp * csh[p]^2
  }
  pred <- (num / den) * csh[i]
  manual_sse <- manual_sse + (bjc[i] - pred)^2
}
cat(sprintf("fast=%.10f manual=%.10f |diff|=%.2e\n", fast_sse, manual_sse, abs(fast_sse - manual_sse)))

cat("\n=== 4. Bandwidth CV: interior optimum, SSE no worse than either H_T endpoint ===\n")
h_opt <- exuber:::contagion_bandwidth_cv(bc, bj, n, d)
H_T <- c(m^(-1 / 2), m^(-1 / 10))
cat("H_T:", H_T, " h_opt:", h_opt, "\n")
cat(
  "SSE at h_opt:", exuber:::contagion_loocv_sse(h_opt, bc, bj, n, d),
  " SSE at H_T[1]:", exuber:::contagion_loocv_sse(H_T[1], bc, bj, n, d),
  " SSE at H_T[2]:", exuber:::contagion_loocv_sse(H_T[2], bc, bj, n, d), "\n"
)

cat("\n=== 5. Directional sensible-behavior check: planted contagion vs. independent series ===\n")
set.seed(42)
nrep <- 15
planted_range <- indep_range <- numeric(nrep)
for (i in seq_len(nrep)) {
  set.seed(1000 + i)
  core_i <- cumsum(rnorm(n))
  y_planted <- numeric(n)
  y_planted[1:10] <- rnorm(10)
  for (t in 11:n) {
    local_rho <- 0.5 + 0.4 * tanh((core_i[max(t - 3, 1)] - core_i[max(t - 13, 1)]) / 5)
    y_planted[t] <- local_rho * y_planted[t - 1] + rnorm(1)
  }
  y_indep <- cumsum(rnorm(n))

  out_planted <- contagion_reg(y_planted, core_i, S = S, d = 3, h = 0.3)
  out_indep <- contagion_reg(y_indep, core_i, S = S, d = 3, h = 0.3)
  planted_range[i] <- diff(range(out_planted$delta2))
  indep_range[i] <- diff(range(out_indep$delta2))
}
cat("mean range(delta2), planted contagion: ", mean(planted_range), "\n")
cat("mean range(delta2), independent series:", mean(indep_range), "\n")

cat("\ndone\n")
