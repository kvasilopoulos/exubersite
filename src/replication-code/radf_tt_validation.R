# Replication script for radf_tt()/radf_tt_cv() (STADF/GSTADF, Kurozumi,
# Skrobotov & Tsarev 2024). Archived retroactively -- this item predates the
# replication/ convention; see docs/enhancements/volatility-robustness.md,
# "Time-transformed test (STADF / GSTADF)", "Independent validation
# (2026-08-09)" for the narrative this reproduces.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Formula check: gls_dfstat_grid() vs brute-force lm() ===\n")
set.seed(4242)
y <- cumsum(rnorm(65))
minw <- 15
res <- exuber:::gls_dfstat_grid(y, minw)

n1 <- length(y) - 1L
dy <- diff(y - y[1])
ylag <- (y - y[1])[1:n1]
b_idx <- minw:n1
badf_lm <- vapply(b_idx, function(b) {
  fit <- lm(dy[1:b] ~ ylag[1:b] - 1)
  summary(fit)$coefficients[1, "t value"]
}, numeric(1))

cat("max|badf_formula - badf_lm| =", max(abs(res$badf - badf_lm)), "\n\n")

cat("=== 2. Critical-value replication vs Whitehouse (2019)'s published STADF triple ===\n")
cat("Published (r0 = 0.1): 2.319, 2.626, 3.223\n")
set.seed(99001)
n <- 300
minw2 <- 30
sadf <- vapply(replicate(4000, exuber:::gls_dfstat_grid(cumsum(rnorm(n)), minw2),
  simplify = FALSE
), `[[`, numeric(1), "sadf")
own_mc <- quantile(sadf, c(0.9, 0.95, 0.99))
cat("Own independent MC (n=300, nrep=4000, seed=99001):", round(own_mc, 3), "\n")

cv <- radf_tt_cv(n = 300, minw = 30, nrep = 4000, seed = 555)
cat("radf_tt_cv(n=300, minw=30, nrep=4000, seed=555) sadf_cv:", round(cv$sadf_cv, 3), "\n\n")

cat("=== 3. Behavioral sanity check: unit-root -> explosive -> collapse series ===\n")
set.seed(1)
normal_part <- cumsum(rnorm(90))
expl_part <- normal_part[90] * 1.04^(1:60) + cumsum(rnorm(60, sd = 0.5))
y_bubble <- c(normal_part, expl_part)
res_bubble <- radf_tt(y_bubble, minw = 20)
cat("gsadf =", round(res_bubble$gsadf, 3), "vs critical-value range ~2.0-3.3 above\n\n")

cat("=== 4. Full test-tt.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-tt.R",
  reporter = "summary"
)
