Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Formula check: cusum_stat_path() vs independent brute-force loop ===\n")
set.seed(1)
n <- 80
T_star <- 40
y <- cumsum(rnorm(n))
b_alpha <- 4.6
res <- exuber:::cusum_stat_path(y, T_star, b_alpha)

S_brute <- numeric(n - T_star)
bnd_brute <- numeric(n - T_star)
for (k in seq_len(n - T_star)) {
  t <- T_star + k
  dy <- diff(y[1:t])
  sigma2_t <- sum(dy^2) / (t - 1)
  S_brute[k] <- (y[t] - y[T_star]) / sqrt(sigma2_t)
  c_t <- sqrt(b_alpha + log(t / T_star))
  bnd_brute[k] <- c_t * sqrt(t)
}
cat("max abs diff S:", max(abs(res$S - S_brute)), "\n")
cat("max abs diff boundary:", max(abs(res$boundary - bnd_brute)), "\n\n")

cat("=== 2. Empirical size under H0 (pure random walk, no bubble anywhere) ===\n")
run_null <- function(seed) {
  set.seed(seed)
  n <- 150
  y <- cumsum(rnorm(n))
  out <- monitor_cusum(y, r_star = 0.5, b_alpha = 4.6)
  !is.na(out$alarm)
}
false_alarm_rate <- mean(sapply(1:100, run_null))
cat(sprintf("Cumulative false-alarm rate (100 reps, H0 throughout, 75 monitoring obs): %.3f\n",
            false_alarm_rate))
cat("(This is HB's own asymptotic UPPER BOUND on the cumulative FPR at b_alpha=4.6,\n",
    " i.e. their eq 28's exp(-b/2) bound for the whole monitoring horizon, not a\n",
    " per-point 5% level -- so a rate notably below 5% would be consistent with a\n",
    " genuinely conservative bound; well above 5% would be a red flag.)\n\n")

cat("=== 3. Detection power under a genuine post-training bubble ===\n")
run_detect <- function(seed) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- monitor_cusum(y, r_star = n1 / length(y), b_alpha = 4.6)
  c(alarm = unname(out$alarm), true_origination = n1)
}
res2 <- t(sapply(1:30, run_detect))
detected <- !is.na(res2[, "alarm"])
cat("Detection rate:", mean(detected), "\n")
if (any(detected)) {
  cat("Alarm delay (alarm - true origination) among detections:\n")
  print(summary(res2[detected, "alarm"] - res2[detected, "true_origination"]))
}
