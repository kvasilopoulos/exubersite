# Validation script for radf_monitor(..., boundary = "fluc") (Homm &
# Breitung 2012's FLUC monitoring detector). See docs/enhancements/
# monitoring.md for the full write-up. Run from the exuber/ package root
# (or adjust the devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. Table lookup sanity checks (Homm & Breitung Table 7(i)) ===\n")
cat("level=0.95, n=100, k=2 (expect 4.50):", exuber:::hb_fluc_q(0.95, 100, 2), "\n")
cat("level=0.95, n=100, k=10 (expect 6.26):", exuber:::hb_fluc_q(0.95, 100, 10), "\n")
cat("level=0.95, n=50, k=4 (expect 5.11):", exuber:::hb_fluc_q(0.95, 50, 4), "\n")
cat("level=0.90, n=20, k=2 (expect 2.49):", exuber:::hb_fluc_q(0.90, 20, 2), "\n")
cat("level=0.99, n=100, k=8 (expect 9.79):", exuber:::hb_fluc_q(0.99, 100, 8), "\n")
cat("level=0.95, n=73, k=7 (snaps n->50, k->6, expect 5.50):", exuber:::hb_fluc_q(0.95, 73, 7), "\n")
tryCatch(exuber:::hb_fluc_q(0.93, 100, 2), error = function(e) cat("level=0.93 correctly errors:", conditionMessage(e), "\n"))

cat("\n=== 2. Basic run with boundary='fluc' ===\n")
set.seed(1)
y <- cumsum(rnorm(150))
out <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = "fluc", level = 0.95)
print(out)

cat("\n=== 3. False-alarm rate under H0: fluc vs kurozumi vs bootstrap ===\n")
run_null <- function(seed, boundary) {
  set.seed(seed)
  y <- cumsum(rnorm(150))
  out <- radf_monitor(y, r_star = 0.5, minw = 20, boundary = boundary, nboot = 99, seed = 1)
  !is.na(out$alarm)
}
cat(sprintf("fluc:      %.3f\n", mean(sapply(1:100, function(s) run_null(s, "fluc")))))
cat(sprintf("kurozumi:  %.3f\n", mean(sapply(1:100, function(s) run_null(s, "kurozumi")))))
cat(sprintf("bootstrap: %.3f\n", mean(sapply(1:100, function(s) run_null(s, "bootstrap")))))

cat("\n=== 4. Detection power under a genuine post-training bubble ===\n")
cat("(HB's own paper: FLUC/CUSUM monitoring generally has LESS power than\n")
cat("a supDF-style test, though FLUC beats CUSUM -- a lower detection rate\n")
cat("here than kurozumi/bootstrap is consistent with that, not a defect)\n\n")
run_detect <- function(seed, boundary) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_monitor(y, r_star = n1 / length(y), minw = 20, boundary = boundary, nboot = 99, seed = 1)
  !is.na(out$alarm)
}
cat(sprintf("fluc:      %.3f\n", mean(sapply(1:30, function(s) run_detect(s, "fluc")))))
cat(sprintf("kurozumi:  %.3f\n", mean(sapply(1:30, function(s) run_detect(s, "kurozumi")))))
cat(sprintf("bootstrap: %.3f\n", mean(sapply(1:30, function(s) run_detect(s, "bootstrap")))))

cat("\n=== 5. Structural check: alarm never before T_star, for fluc boundary ===\n")
run_check <- function(seed) {
  set.seed(seed)
  n1 <- 75; n2 <- 40
  normal_part <- cumsum(rnorm(n1))
  expl_part <- normal_part[n1] * 1.05^(1:n2) + cumsum(rnorm(n2, sd = 0.3))
  y <- c(normal_part, expl_part)
  out <- radf_monitor(y, r_star = n1 / length(y), minw = 20, boundary = "fluc")
  alarm <- unname(out$alarm)
  if (is.na(alarm)) NA else alarm > out$T_star
}
cat("all TRUE or NA:", all(na.omit(sapply(1:10, run_check))), "\n")
