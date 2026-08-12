Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Basic run: no bubble anywhere -- should mostly not alarm ===\n")
set.seed(1)
n <- 150
y_null <- cumsum(rnorm(n))
out_null <- radf_monitor(y_null, r_star = 0.5, minw = 20, nboot = 199, seed = 1)
print(out_null)

cat("\n=== 2. Empirical false-alarm rate under H0 (no bubble anywhere),",
    "should be roughly <= nominal level over the monitoring horizon ===\n")
run_null <- function(seed) {
  set.seed(seed)
  n <- 150
  y <- cumsum(rnorm(n))
  out <- radf_monitor(y, r_star = 0.5, minw = 20, nboot = 199, seed = 1)
  !is.na(out$alarm)
}
false_alarm_rate <- mean(sapply(1:40, run_null))
cat(sprintf("False alarm rate (40 reps, H0 throughout, T-T*=75 monitoring obs): %.3f\n",
            false_alarm_rate))
cat("(NOTE: this is a CUMULATIVE false-alarm probability over 75 monitoring points\n",
    " at a per-point 95% threshold, so a much higher rate than 5% is expected and\n",
    " not itself a red flag -- Family A/PSY-style monitoring's FPR grows with the\n",
    " monitoring horizon, exactly as flagged in monitoring.md's own AHLST/Whitehouse\n",
    " discussion; this is descriptive, not a strict pass/fail bound.)\n\n")

cat("=== 3. A genuine bubble starting AFTER T* should be detected, and",
    "the alarm date should be close to the true origination date ===\n")
run_detect <- function(seed) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_monitor(y, r_star = n1 / length(y), minw = 20, nboot = 199, seed = 1)
  c(alarm = unname(out$alarm), true_origination = n1)
}
res <- t(sapply(1:15, run_detect))
cat("Detection rate:", mean(!is.na(res[, "alarm"])), "\n")
cat("Alarm delay (alarm - true origination) among detections:\n")
print(summary(res[!is.na(res[, "alarm"]), "alarm"] - res[!is.na(res[, "alarm"]), "true_origination"]))

cat("\n=== 4. No false alarm strictly WITHIN the training window",
    "(monitoring only starts after T*) ===\n")
# check monitoring rows never include anything before T*
set.seed(5)
n <- 150; T_star <- 75
y <- cumsum(rnorm(n))
out <- radf_monitor(y, r_star = T_star, minw = 20, nboot = 99, seed = 1)
cat("T_star:", out$T_star, " n bsadf rows:", nrow(out$bsadf), "\n")
