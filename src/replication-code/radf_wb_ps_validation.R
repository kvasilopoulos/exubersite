# Replication script for radf_wb_ps_cv() and radf_wb_ps_distr() (the
# Phillips & Shi (2020) wild bootstrap variant in R/radf_wb.R) and the shared
# OLS lag-selection and AR-fit routines it needs (adf_res() and lag_select(),
# also in R/radf_wb.R).
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. lag_select()/adf_res(): deterministic (no RNG), bit-for-bit reference ===\n")
set.seed(11)
y <- cumsum(rnorm(40))

cat("lag_select(y, 'aic', max_lag=5):", exuber:::lag_select(y, criterion = "aic", max_lag = 5), "\n")
cat("lag_select(y, 'bic', max_lag=5):", exuber:::lag_select(y, criterion = "bic", max_lag = 5), "\n")

ar <- exuber:::adf_res(y, adflag = 2, type = "fixed")
cat("adf_res(y, adflag=2, type='fixed')$beta:", sprintf("%.10f", ar$beta), "\n")
cat("adf_res(...)$res[1:5]:", sprintf("%.10f", ar$res[1:5]), "\n")
cat("adf_res(...)$res length:", length(ar$res), "\n\n")

cat("=== 2. radf_wb_ps_cv(): structural + shape checks ===\n")
cat("(bootstrap draws use R's own RNG, not comparable bit-for-bit to a\n")
cat(" numpy Generator port -- this checks shape/monotonicity/sanity only)\n")
minw <- psy_minw(length(y))
wb <- radf_wb_ps_cv(y, minw = minw, nboot = 80, adflag = 0, seed = 5)
cat("class:", paste(class(wb), collapse = ","), "\n")
cat("adf_cv dim:", dim(wb$adf_cv), "\n")
cat("gsadf_cv dim:", dim(wb$gsadf_cv), "\n")
cat("bsadf_cv dim:", dim(wb$bsadf_cv), "\n")
cat("gsadf_cv monotonic (90 <= 95 <= 99):", all(diff(as.vector(wb$gsadf_cv)) >= 0), "\n\n")

cat("=== 3. tb (training-window) mode: badf/bsadf collapse to repeated sadf/gsadf ===\n")
tb <- minw + 10
wbt <- radf_wb_ps_cv(y, minw = minw, nboot = 60, adflag = 0, tb = tb, seed = 6)
pointer_full <- length(y) - minw
cat("badf_cv dim (expect", pointer_full, "x 3 x 1):", dim(wbt$badf_cv), "\n")
cat(
  "badf_cv[1,,1] == sadf_cv[1,]:",
  isTRUE(all.equal(as.vector(wbt$badf_cv[1, , 1]), as.vector(wbt$sadf_cv[1, ]))), "\n"
)
cat(
  "bsadf_cv[1,,1] == gsadf_cv[1,]:",
  isTRUE(all.equal(as.vector(wbt$bsadf_cv[1, , 1]), as.vector(wbt$gsadf_cv[1, ]))), "\n"
)

# Reference output (R 4.6.1, seed as above):
#
# lag_select(y, 'aic', max_lag=5): 5
# lag_select(y, 'bic', max_lag=5): 4
# adf_res(y, adflag=2, type='fixed')$beta: -0.2766891117 -0.0510771847 0.0954335649
# adf_res(...)$res[1:5]: -1.1659634957 1.5303078394 -0.4672254328 1.4401135170 1.0583623423
# adf_res(...)$res length: 37
#
# adf_cv dim: 1 3 ; gsadf_cv dim: 1 3 ; bsadf_cv dim: 29 3 1
# gsadf_cv monotonic: TRUE
#
# badf_cv dim: 29 3 1 ; badf_cv[1,,1] == sadf_cv[1,]: TRUE ;
# bsadf_cv[1,,1] == gsadf_cv[1,]: TRUE
