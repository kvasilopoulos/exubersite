# Replication script for explosive_root()/root_ci()/root_ci_datestamp()
# (root inference: Guo, Sun & Wang 2019 normal-t CI + Phillips-Magdalinos
# 2007 Cauchy CI). Archived retroactively -- see
# docs/enhancements/dating-and-root-inference.md, "Root inference",
# "Independent validation (2026-08-09)" for the narrative this reproduces.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Cauchy percentiles (Skrobotov 2023 review's footnote 17) ===\n")
cat("Published: C_0.10=6.315, C_0.05=12.7, C_0.01=63.65674\n")
cat("qt(0.95, df=1) =", qt(0.95, df = 1), "\n")
cat("qt(0.975, df=1) =", qt(0.975, df = 1), "\n")
cat("qt(0.995, df=1) =", qt(0.995, df = 1), "\n\n")

cat("=== 2. Point estimate: super-consistency at rho=1.05, n=200 ===\n")
set.seed(8675309)
n <- 200
y <- numeric(n)
e <- rnorm(n)
for (t in 2:n) y[t] <- 1.05 * y[t - 1] + e[t]
er <- explosive_root(y, 1, n)
ci <- root_ci(er)
cat("rho_hat =", round(ci$rho, 4), "vs rho_true = 1.0500\n")
cat("95% CI:", round(ci$rho_ci, 6), "\n\n")

cat("=== 3. Coverage (rho=1.05, n=200, 800 reps, seed=24601) ===\n")
cat("Published: 94.6% coverage of nominal 95%\n")
set.seed(24601)
covered <- replicate(800, {
  y <- numeric(200)
  e <- rnorm(200)
  for (t in 2:200) y[t] <- 1.05 * y[t - 1] + e[t]
  ci <- root_ci(explosive_root(y, 1, 200))
  ci$rho_ci[1] <= 1.05 && 1.05 <= ci$rho_ci[2]
})
cat("Own coverage:", mean(covered), "\n\n")

cat("=== 4. Cauchy-type CI (eq. 27, Phillips-Magdalinos 2007) formula-exact check ===\n")
set.seed(11)
y2 <- numeric(150)
e2 <- rnorm(150)
for (t in 2:150) y2[t] <- 1.04 * y2[t - 1] + e2[t]
er2 <- explosive_root(y2, 1, 150)
ci_cauchy <- root_ci(er2, type = "cauchy", level = 0.95)
q <- qt(0.975, df = 1)
rho_hat <- er2$rho
n2 <- er2$n
half_width_formula <- q * (rho_hat^2 - 1) / rho_hat^n2
cat("root_ci(type='cauchy') half-width:", ci_cauchy$rho_ci[2] - ci_cauchy$rho, "\n")
cat("eq. 27 formula half-width:        ", half_width_formula, "\n\n")

cat("=== 5. root_ci_datestamp() end-to-end vs calling explosive_root()/root_ci() directly ===\n")
set.seed(42)
normal_part <- cumsum(rnorm(90))
expl_part <- normal_part[90] * 1.04^(1:60) + cumsum(rnorm(60, sd = 0.5))
y3 <- c(normal_part, expl_part)
r <- radf(y3, minw = 20)
cv <- radf_mc_cv(length(y3), minw = 20, nrep = 500, seed = 1)
ds <- datestamp(r, cv = cv)
rc <- root_ci_datestamp(r, ds)
cat("root_ci_datestamp() episodes found:", nrow(rc[[1]]), "\n")
print(rc)

cat("\n=== Full test-root-ci.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-root-ci.R",
  reporter = "summary"
)
