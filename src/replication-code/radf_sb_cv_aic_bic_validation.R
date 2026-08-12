# Replication script for radf_sb_cv(type = "aic"/"bic") (Pedersen &
# Schuette 2020's automatic lag-order selection for the sieve bootstrap,
# Bundle 1, 2026-08-09). Archived retroactively -- see
# docs/enhancements/volatility-robustness.md, "Pedersen & Schuette sieve
# bootstrap" for the narrative this reproduces. Different seeds/sample
# sizes from test-sb.R throughout.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== 1. type = 'fixed' (pre-existing default) unaffected by the new arguments ===\n")
set.seed(101)
y <- cumsum(rnorm(90))
a <- radf_sb_cv(y, lag = 2, nboot = 60, seed = 22)
b <- radf_sb_cv(y, lag = 2, type = "fixed", nboot = 60, seed = 22)
cat("identical gsadf_panel_cv:", isTRUE(all.equal(a$gsadf_panel_cv, b$gsadf_panel_cv)), "\n\n")

cat("=== 2. type = 'aic'/'bic' picks up a nonzero lag on AR(2)-autocorrelated data ===\n")
set.seed(202)
n <- 200
e <- arima.sim(list(ar = c(0.35, 0.25)), n = n)
y2 <- cumsum(as.numeric(e))
sb_bic <- radf_sb_cv(y2, type = "bic", max_lag = 5, nboot = 60, seed = 22)
sb_aic <- radf_sb_cv(y2, type = "aic", max_lag = 5, nboot = 60, seed = 22)
cat("class:", class(sb_bic), "\n")
cat("selected lag (bic):", attr(sb_bic, "lag"), "\n")
cat("selected lag (aic):", attr(sb_aic, "lag"), "\n\n")

cat("=== 3. type = 'bic' modal lag on pure random-walk data (no true\n")
cat("    autocorrelation) is 0, across 10 independent draws ===\n")
lags <- vapply(1:10, function(s) {
  set.seed(s + 900)
  y <- cumsum(rnorm(120))
  attr(radf_sb_cv(y, type = "bic", max_lag = 6, nboot = 25, seed = 22), "lag")
}, integer(1))
cat("selected lags across 10 draws:", lags, "\n")
cat("modal lag:", as.integer(names(sort(table(lags), decreasing = TRUE))[1]), "\n\n")

cat("=== Full test-sb.R suite ===\n")
testthat::test_file(
  "c:/Users/User/Documents/05-R/exuber-project/exuber/tests/testthat/test-sb.R",
  reporter = "summary"
)
