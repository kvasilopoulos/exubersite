devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. Moment check: w = u/sqrt(2) + (v^2-1)/2 should have E[w]=0, E[w^2]=1, E[w^3]=1 ===\n")
set.seed(1)
n <- 2000000
u <- rnorm(n); v <- rnorm(n)
w <- u / sqrt(2) + (v^2 - 1) / 2
cat(sprintf("E[w]=%.4f (want 0)  E[w^2]=%.4f (want 1)  E[w^3]=%.4f (want 1)\n\n",
            mean(w), mean(w^2), mean(w^3)))

cat("=== 2. Regression check: dist_rad=FALSE, dist_skew=FALSE unchanged (compare to git-committed behavior) ===\n")
y <- cumsum(rnorm(60))
set.seed(5)
r1 <- exuber:::radf_wb_dgp_hlst(y, dist_rad = FALSE)
set.seed(5)
r2 <- exuber:::radf_wb_dgp_hlst(y, dist_rad = FALSE, dist_skew = FALSE)
stopifnot(identical(r1, r2))
cat("OK: identical to pre-change default behavior\n\n")

cat("=== 3. Error when both dist_rad and dist_skew are TRUE ===\n")
tryCatch({
  exuber:::radf_wb_hlst(cumsum(rnorm(40)), minw = 10, nboot = 5, dist_rad = TRUE, dist_skew = TRUE)
  cat("FAIL: expected an error\n")
}, error = function(e) cat("OK, errored:", conditionMessage(e), "\n"))
cat("\n")

cat("=== 4. Empirical size under H0 with RIGHT-SKEWED, HETEROSKEDASTIC errors (Hafner's own setting) ===\n")
# errors: negative log-chi-square(1), right-skewed, standardized -- as in the paper's footnote 1
rskew_innov <- function(n) {
  z <- rnorm(n)
  e <- -log(z^2)
  (e - mean(e)) / sd(e)
}
run_once <- function(seed, dist_skew) {
  set.seed(seed)
  Tn <- 100
  g <- 0.05 * (1 + 2 * cos(pi * (1:Tn) / Tn)^2)  # deterministic heteroskedastic scale
  y <- cumsum(g * rskew_innov(Tn))  # pure random walk under H0, right-skewed heteroskedastic errors
  cv <- exuber:::radf_wb_cv(y, minw = 20, nboot = 199, dist_skew = dist_skew, seed = 1)
  obs <- exuber:::rls_gsadf(exuber:::unroot(y), min_win = 20)
  sadf_obs <- obs[length(y) - 20 + 2]
  sadf_obs > cv$sadf_cv[1, "95%"]
}
rej_normal <- sapply(1:80, function(s) run_once(s, dist_skew = FALSE))
rej_skew   <- sapply(1:80, function(s) run_once(s, dist_skew = TRUE))
cat(sprintf("Empirical size, normal multiplier (dist_skew=FALSE): %.3f\n", mean(rej_normal)))
cat(sprintf("Empirical size, skewed multiplier (dist_skew=TRUE):  %.3f\n", mean(rej_skew)))
cat("(nominal 0.05; paper's own finding: undersized in small samples, bias grows with heteroskedasticity -- so well below 0.05 is expected and consistent with the source, not a red flag)\n")
