Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("c:/Users/User/Documents/05-R/exuber-project/exuber", quiet = TRUE)

cat("=== Cross-check radf_sign_cv() at EXACT finite T=100, vs paper's own T=100 table row ===\n")
cat("(avoids T->Inf convergence ambiguity -- paper documents sPSY converges slowly)\n")
cat("Paper T=100: sPWY (10%,5%,1%) = (2.470, 2.859, 3.656); sPSY = (4.381, 5.578, 13.056)\n\n")

n <- 100
minw <- round(0.1 * n)
cv <- radf_sign_cv(n, minw = minw, nrep = 2000, seed = 1)
cat("Simulated sadf_cv (-> sPWY):", cv$sadf_cv, "\n")
cat("Simulated gsadf_cv (-> sPSY):", cv$gsadf_cv, "\n\n")

cat("=== Also T=200 for a second data point ===\n")
cat("Paper T=200: sPWY = (2.405, 2.735, 3.434); sPSY = (3.469, 3.901, 4.957)\n")
n2 <- 200
minw2 <- round(0.1 * n2)
cv2 <- radf_sign_cv(n2, minw = minw2, nrep = 1500, seed = 1)
cat("Simulated sadf_cv (-> sPWY):", cv2$sadf_cv, "\n")
cat("Simulated gsadf_cv (-> sPSY):", cv2$gsadf_cv, "\n")
