Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== Power check with dist_skew=TRUE bootstrap, ORDINARY normal innovations ===\n")
run_power <- function(seed) {
  set.seed(seed)
  Tn <- 100
  Te <- 60
  normal_part <- cumsum(rnorm(Te))
  expl_part <- normal_part[Te] * 1.06^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.3))
  y <- c(normal_part, expl_part)
  obs <- radf(y, minw = 20)$sadf
  cv <- radf_wb_cv(y, minw = 20, nboot = 199, dist_skew = TRUE, seed = 1)
  obs > cv$sadf_cv[1, "95%"]
}
power <- mean(sapply(1:30, run_power))
cat(sprintf("Empirical power (30 reps): %.3f\n\n", power))

cat("=== Size under H0 with the paper's own right-skewed innovation distribution
     (negative log-chi-square(1)), no heteroskedasticity added on top ===\n")
rskew_innov <- function(n) {
  z <- rnorm(n)
  e <- -log(z^2)
  (e - mean(e)) / sd(e)
}
run_size <- function(seed, dist_skew) {
  set.seed(seed)
  Tn <- 150
  y <- cumsum(rskew_innov(Tn))  # pure random walk under H0
  obs <- radf(y, minw = 20)$sadf
  cv <- radf_wb_cv(y, minw = 20, nboot = 199, dist_skew = dist_skew, seed = 1)
  obs > cv$sadf_cv[1, "95%"]
}
rej_normal <- mean(sapply(1:60, function(s) run_size(s, FALSE)))
rej_skew   <- mean(sapply(1:60, function(s) run_size(s, TRUE)))
cat(sprintf("Empirical size, normal multiplier: %.3f\n", rej_normal))
cat(sprintf("Empirical size, skewed multiplier: %.3f\n", rej_skew))
cat("(nominal 0.05)\n")
