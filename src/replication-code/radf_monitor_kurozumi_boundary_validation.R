Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Table lookup sanity checks ===\n")
cat("level=0.95, s_bar=1 (expect 1.0381):", exuber:::kurozumi_sadf_q(0.95, 1), "\n")
cat("level=0.95, s_bar=1.2 (snap to 1, expect 1.0381):", exuber:::kurozumi_sadf_q(0.95, 1.2), "\n")
cat("level=0.95, s_bar=3 (expect 1.3330):", exuber:::kurozumi_sadf_q(0.95, 3), "\n")
cat("level=0.95, s_bar=5 (expect 1.4255):", exuber:::kurozumi_sadf_q(0.95, 5), "\n")
cat("level=0.90, s_bar=1 (expect 0.6946):", exuber:::kurozumi_sadf_q(0.90, 1), "\n")
cat("level=0.99, s_bar=1 (expect 1.6474):", exuber:::kurozumi_sadf_q(0.99, 1), "\n")
tryCatch(exuber:::kurozumi_sadf_q(0.93, 1), error = function(e) cat("level=0.93 correctly errors:", conditionMessage(e), "\n"))
cat("\n")

cat("=== 2. Basic run with boundary='kurozumi' ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
out <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = "kurozumi", level = 0.95)
print(out)
cat("\n")

cat("=== 3. Empirical false-alarm rate under H0, boundary='kurozumi' vs 'bootstrap' ===\n")
run_null <- function(seed, boundary) {
  set.seed(seed)
  y <- cumsum(rnorm(150))
  out <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = boundary, nboot = 99, seed = 1)
  !is.na(out$alarm)
}
rate_kuro <- mean(sapply(1:100, function(s) run_null(s, "kurozumi")))
rate_boot <- mean(sapply(1:100, function(s) run_null(s, "bootstrap")))
cat(sprintf("False-alarm rate, kurozumi boundary (100 reps, target ~10%% since s_bar=(150-75)/75=1 -> nearest table row): %.3f\n", rate_kuro))
cat(sprintf("False-alarm rate, bootstrap boundary (100 reps): %.3f\n\n", rate_boot))

cat("=== 4. Detection power under a genuine post-training bubble ===\n")
run_detect <- function(seed, boundary) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_monitor(y, r_star = n1 / length(y), minw = 20, boundary = boundary, nboot = 99, seed = 1)
  c(alarm = unname(out$alarm), true_origination = n1)
}
res_kuro <- t(sapply(1:30, function(s) run_detect(s, "kurozumi")))
res_boot <- t(sapply(1:30, function(s) run_detect(s, "bootstrap")))
cat(sprintf("Detection rate, kurozumi:  %.3f\n", mean(!is.na(res_kuro[, "alarm"]))))
cat(sprintf("Detection rate, bootstrap: %.3f\n", mean(!is.na(res_boot[, "alarm"]))))

cat("\n=== 5. Structural check: alarm never before T_star, for kurozumi boundary ===\n")
run_check <- function(seed) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_monitor(y, r_star = n1 / length(y), minw = 20, boundary = "kurozumi")
  alarm <- unname(out$alarm)
  if (is.na(alarm)) NA else alarm > out$T_star
}
cat("all TRUE or NA:", all(na.omit(sapply(1:10, run_check))), "\n")
