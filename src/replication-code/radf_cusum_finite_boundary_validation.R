# Validation script for radf_cusum(..., boundary = "finite") (Homm &
# Breitung 2012's finite-sample CUSUM boundary, their Table 8). See
# docs/enhancements/monitoring.md, "Implementation (CUSUM)", for the
# full write-up. Run from the exuber/ package root (or adjust the
# devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Table lookup sanity checks (Homm & Breitung Table 8(i)) ===\n")
cat("level=0.95, n=100, k=2 (expect 1.51):", exuber:::hb_cusum_finite_q(0.95, 100, 2), "\n")
cat("level=0.95, n=100, k=10 (expect 3.36):", exuber:::hb_cusum_finite_q(0.95, 100, 10), "\n")
cat("level=0.90, n=20, k=2 (expect 0.81):", exuber:::hb_cusum_finite_q(0.90, 20, 2), "\n")
tryCatch(exuber:::hb_cusum_finite_q(0.93, 100, 2), error = function(e) cat("level=0.93 correctly errors:", conditionMessage(e), "\n"))

cat("\n=== 2. Basic run with boundary='finite' ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
out <- radf_cusum(y, r_star = 0.5, boundary = "finite", level = 0.95)
print(out)
cat("b_alpha used:", attr(out, "b_alpha"), "(vs asymptotic default 4.6)\n")

cat("\n=== 3. False-alarm rate and detection power: asymptotic vs finite ===\n")
run_null <- function(seed, boundary) {
  set.seed(seed)
  y <- cumsum(rnorm(150))
  out <- radf_cusum(y, r_star = 0.5, boundary = boundary)
  !is.na(out$alarm)
}
run_detect <- function(seed, boundary) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_cusum(y, r_star = n1 / length(y), boundary = boundary)
  !is.na(out$alarm)
}
cat(sprintf("False-alarm rate, asymptotic: %.3f\n", mean(sapply(1:100, function(s) run_null(s, "asymptotic")))))
cat(sprintf("False-alarm rate, finite:     %.3f\n\n", mean(sapply(1:100, function(s) run_null(s, "finite")))))
cat(sprintf("Detection rate, asymptotic: %.3f\n", mean(sapply(1:30, function(s) run_detect(s, "asymptotic")))))
cat(sprintf("Detection rate, finite:     %.3f\n", mean(sapply(1:30, function(s) run_detect(s, "finite")))))

cat("\n=== 4. boundary='finite' also works with type='kernel' (CUSUMV), ===\n")
cat("    extending Corollary 1's shared-boundary result to the table ===\n")
out_k <- radf_cusum(y, r_star = 0.5, boundary = "finite", type = "kernel")
print(out_k)
