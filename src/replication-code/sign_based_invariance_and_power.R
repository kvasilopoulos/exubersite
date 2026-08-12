Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Cross-check radf_sign_cv() vs paper's Table 1 asymptotic (T=Inf) values ===\n")
cat("Paper: sPWY (10%,5%,1%) = (2.410, 2.734, 3.248); sPSY = (2.933, 3.180, 3.655)\n")
n <- 300
minw <- round(0.1 * n)
cv <- radf_sign_cv(n, minw = minw, nrep = 1500, seed = 1)
cat("Simulated sadf_cv (-> sPWY):", cv$sadf_cv, "\n")
cat("Simulated gsadf_cv (-> sPSY):", cv$gsadf_cv, "\n\n")

cat("=== 2. Exact invariance to heteroskedasticity: same sign-path, different scale ===\n")
set.seed(7)
n2 <- 150
Te <- 90
base_incr <- rnorm(Te)  # normal-part increments
expl_incr <- c(rnorm(1, mean = 3), rnorm(n2 - Te - 1))  # some explosive-ish increments
raw_dy <- c(base_incr, expl_incr)
y_homo <- cumsum(raw_dy)  # constant volatility

# now scale by a wild heteroskedastic pattern -- SAME SIGNS, different magnitude
vol_pattern <- c(rep(0.1, 40), rep(10, 60), rep(1, n2 - 100))
y_hetero <- cumsum(raw_dy * vol_pattern)

r_homo <- radf_sign(y_homo, minw = 20)
r_hetero <- radf_sign(y_hetero, minw = 20)
cat(sprintf("homoskedastic:   sadf=%.6f gsadf=%.6f\n", r_homo$sadf, r_homo$gsadf))
cat(sprintf("heteroskedastic: sadf=%.6f gsadf=%.6f\n", r_hetero$sadf, r_hetero$gsadf))
cat("identical:", isTRUE(all.equal(r_homo$sadf, r_hetero$sadf)) &&
                    isTRUE(all.equal(r_homo$gsadf, r_hetero$gsadf)), "\n\n")

cat("=== 3. Power check: does radf_sign correctly flag a clear bubble via its own cv? ===\n")
run_power <- function(seed) {
  set.seed(seed)
  Tn <- 150
  Te <- 90
  normal_part <- cumsum(rnorm(Te))
  expl_part <- normal_part[Te] * 1.05^(1:(Tn - Te)) + cumsum(rnorm(Tn - Te, sd = 0.5))
  y <- c(normal_part, expl_part)
  obs <- radf_sign(y, minw = 20)
  obs$gsadf > cv2$gsadf_cv["95%"]
}
cv2 <- radf_sign_cv(150, minw = 20, nrep = 500, seed = 2)
power <- mean(sapply(1:30, run_power))
cat(sprintf("Empirical power (30 reps, gsadf vs simulated 95%% cv at n=150): %.3f\n", power))
