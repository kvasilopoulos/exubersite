Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Formula check: one_sided_kernel_spot_vol() vs brute-force loop ===\n")
set.seed(1)
n <- 60
dy <- rnorm(n)
N <- 10
res <- exuber:::one_sided_kernel_spot_vol(dy, N = N, kernel = "gaussian")

sigma2_brute <- numeric(n)
kern <- function(u) dnorm(u)
w_full <- kern((0:N) / N); w_full <- w_full / sum(w_full)
for (j in seq_len(n)) {
  if (j <= N) {
    sigma2_brute[j] <- 1  # paper's own convention
  } else {
    idx <- (j - N):j  # s = 0..N -> dy[j], dy[j-1], ..., dy[j-N]
    ww <- rev(w_full)  # align so w_full[1] (s=0) multiplies dy[j] (last in idx)
    sigma2_brute[j] <- sum(ww * dy[idx]^2)
  }
}
cat("max abs diff (j > N only):", max(abs(res[(N+1):n] - sigma2_brute[(N+1):n])), "\n")
cat("all j<=N equal to 1:", all(res[1:N] == 1), "\n\n")

cat("=== 2. Formula check: cusum_stat_path_kernel() vs independent brute-force ===\n")
set.seed(2)
n2 <- 100
T_star <- 50
y <- cumsum(rnorm(n2))
b_alpha <- 4.6
res2 <- exuber:::cusum_stat_path_kernel(y, T_star, b_alpha, N = 20, kernel = "gaussian")

dy2 <- diff(y)
sigma2_dy <- exuber:::one_sided_kernel_spot_vol(dy2, N = 20, kernel = "gaussian")
SV_brute <- numeric(n2 - T_star)
for (k in seq_len(n2 - T_star)) {
  t <- T_star + k
  js <- (T_star + 1):t  # Delta y_j for j = T_star+1..t -> dy index (j-1)
  SV_brute[k] <- sum(dy2[js - 1] / sqrt(sigma2_dy[js - 1]))
}
cat("max abs diff SV:", max(abs(res2$S - SV_brute)), "\n\n")

cat("=== 3. Under HOMOSKEDASTIC H0: standard vs kernel CUSUM should have",
    "comparable (not necessarily identical) empirical false-alarm rates ===\n")
run_null_homo <- function(seed, type) {
  set.seed(seed)
  y <- cumsum(rnorm(150))
  out <- monitor_cusum(y, r_star = 0.5, b_alpha = 4.6, type = type)
  !is.na(out$alarm)
}
rate_std_homo <- mean(sapply(1:60, function(s) run_null_homo(s, "standard")))
rate_ker_homo <- mean(sapply(1:60, function(s) run_null_homo(s, "kernel")))
cat(sprintf("Standard CUSUM false-alarm rate (homoskedastic): %.3f\n", rate_std_homo))
cat(sprintf("Kernel CUSUMV false-alarm rate (homoskedastic):  %.3f\n\n", rate_ker_homo))

cat("=== 4. Under HETEROSKEDASTIC H0: standard CUSUM should become",
    "OVERSIZED, kernel CUSUMV should stay controlled -- this is the",
    "paper's whole selling point ===\n")
run_null_hetero <- function(seed, type) {
  set.seed(seed)
  n <- 150
  # sharp one-time volatility jump partway through the monitoring region
  vol <- c(rep(1, 90), rep(8, 60))
  y <- cumsum(rnorm(n) * vol)
  out <- monitor_cusum(y, r_star = 0.5, b_alpha = 4.6, type = type)
  !is.na(out$alarm)
}
rate_std_hetero <- mean(sapply(1:60, function(s) run_null_hetero(s, "standard")))
rate_ker_hetero <- mean(sapply(1:60, function(s) run_null_hetero(s, "kernel")))
cat(sprintf("Standard CUSUM false-alarm rate (heteroskedastic, vol jump 1->8): %.3f\n", rate_std_hetero))
cat(sprintf("Kernel CUSUMV false-alarm rate (heteroskedastic, vol jump 1->8):  %.3f\n\n", rate_ker_hetero))

cat("=== 5. Detection power under a genuine post-training bubble (homoskedastic) ===\n")
run_detect <- function(seed, type) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- monitor_cusum(y, r_star = n1 / length(y), type = type)
  !is.na(out$alarm)
}
power_std <- mean(sapply(1:30, function(s) run_detect(s, "standard")))
power_ker <- mean(sapply(1:30, function(s) run_detect(s, "kernel")))
cat(sprintf("Detection rate, standard: %.3f\n", power_std))
cat(sprintf("Detection rate, kernel:   %.3f\n", power_ker))
