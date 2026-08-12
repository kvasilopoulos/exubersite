devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Formula check: coexplosive_stat() vs independent brute-force loop ===\n")
set.seed(1)
Tn <- 100
x <- rnorm(Tn)
y <- 2 + 0.5 * x + rnorm(Tn)
lag <- 2L

res <- exuber:::coexplosive_stat(y, x, lag)

# Independent brute-force recomputation: manual loop, manual OLS via lm(),
# manual cumulative sum loop (not vectorized the same way as the package code)
lo <- max(lag, 0L) + 1L
hi <- Tn + min(lag, 0L)
yy <- y[lo:hi]
xx <- x[(lo:hi) - lag]
fit_lm <- lm(yy ~ xx)
e_brute <- residuals(fit_lm)
n_brute <- length(e_brute)
sigma2_brute <- sum(e_brute^2) / n_brute
S_manual_sum <- 0
running <- 0
for (t in seq_along(e_brute)) {
  running <- running + e_brute[t]
  S_manual_sum <- S_manual_sum + running^2
}
S_brute <- S_manual_sum / (sigma2_brute * n_brute^2)

cat(sprintf("package S=%.8f  brute S=%.8f  diff=%.2e\n", res$S, S_brute, abs(res$S - S_brute)))
stopifnot(abs(res$S - S_brute) < 1e-8)
cat("OK: exact match\n\n")

cat("=== 2. Size under H0, homoskedastic errors (should be near nominal 5%) ===\n")
run_h0_homo <- function(seed) {
  set.seed(seed)
  Tn <- 150
  # x contains a genuine explosive episode (Model 1: unit root -> explosive)
  Te <- 90
  ex <- c(cumsum(rnorm(Te)), NA)
  ex <- ex[1:Te]
  expl <- ex[Te] * 1.05^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  x <- c(ex, expl)
  # y co-explosive with x at lag 0: y = a + b*x + iid noise
  y <- 1 + 0.8 * x + rnorm(Tn, sd = 1)
  out <- exuber:::radf_cobubble(y, x, lag = 0L, nboot = 199L, seed = 1)
  out$reject
}
rejections <- sapply(1:100, run_h0_homo)
cat(sprintf("Empirical size (homoskedastic, 100 reps): %.3f (nominal 0.05)\n\n", mean(rejections)))

cat("=== 3. Size under H0, heteroskedastic errors (tests the wild bootstrap's whole point) ===\n")
run_h0_hetero <- function(seed) {
  set.seed(seed)
  Tn <- 150
  Te <- 90
  ex <- cumsum(rnorm(Te))
  expl <- ex[Te] * 1.05^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  x <- c(ex, expl)
  # heteroskedastic errors: sd jumps from 1 to 4 partway through
  sd_pattern <- c(rep(1, Tn %/% 2), rep(4, Tn - Tn %/% 2))
  y <- 1 + 0.8 * x + rnorm(Tn, sd = sd_pattern)
  out <- exuber:::radf_cobubble(y, x, lag = 0L, nboot = 199L, seed = 1)
  out$reject
}
rejections_hetero <- sapply(1:100, run_h0_hetero)
cat(sprintf("Empirical size (heteroskedastic, 100 reps): %.3f (nominal 0.05)\n\n", mean(rejections_hetero)))

cat("=== 4. Power under H1: y and x are NOT co-explosive (independent explosive episodes) ===\n")
run_h1 <- function(seed) {
  set.seed(seed)
  Tn <- 150
  Te <- 90
  ex <- cumsum(rnorm(Te))
  expl_x <- ex[Te] * 1.05^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  x <- c(ex, expl_x)
  # y has its OWN independent explosive episode (not driven by x at all)
  ey <- cumsum(rnorm(Te))
  expl_y <- ey[Te] * 1.05^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  y <- ey_full <- c(ey, expl_y)
  out <- exuber:::radf_cobubble(y, x, lag = 0L, nboot = 199L, seed = 1)
  out$reject
}
rejections_h1 <- sapply(1:60, run_h1)
cat(sprintf("Empirical power (H1, independent explosive episodes, 60 reps): %.3f (should be high)\n\n", mean(rejections_h1)))

cat("=== 5. Lag recovery: true lag = 3, does coexplosive_select_lag() find it? ===\n")
run_lag_recovery <- function(seed) {
  set.seed(seed)
  Tn <- 200
  Te <- 120
  true_lag <- 3L
  ex <- cumsum(rnorm(Te))
  expl <- ex[Te] * 1.06^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  x <- c(ex, expl)
  # y co-explosive with x lagged by true_lag: y_t = a + b*x_{t-true_lag} + noise
  # build y such that y[t] depends on x[t - true_lag]
  y <- rep(NA_real_, Tn)
  for (t in (true_lag + 1):Tn) y[t] <- 1 + 0.8 * x[t - true_lag] + rnorm(1, sd = 0.5)
  y[1:true_lag] <- x[1:true_lag] + rnorm(true_lag, sd = 0.5)
  est_lag <- exuber:::coexplosive_select_lag(y, x, lags = -6:6)
  est_lag
}
lags_est <- sapply(1:20, run_lag_recovery)
cat("Estimated lags across 20 seeds (true = 3):\n")
print(table(lags_est))
