# Replication script for the radf_obj methods of tidy() and augment()
# (R/radf-tidiers.R). Only the radf_obj methods are ported to pyexuber, not
# the radf_cv and radf_distr tidiers, tidy_join/augment_join or
# summary()/diagnostics(). See the tidy.py module docstring in pyexuber.
Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

set.seed(5)
dta <- data.frame(a = cumsum(rnorm(30)), b = cumsum(rnorm(30)))
rfd <- radf(dta, minw = 10, lag = 1)

cat("adf:", sprintf("%.8f", rfd$adf), "\n")
cat("sadf:", sprintf("%.8f", rfd$sadf), "\n")
cat("gsadf:", sprintf("%.8f", rfd$gsadf), "\n")
cat("gsadf_panel:", sprintf("%.8f", rfd$gsadf_panel), "\n")
cat("badf[,1] (series a):", sprintf("%.8f", rfd$badf[, 1]), "\n")
cat("badf[,2] (series b):", sprintf("%.8f", rfd$badf[, 2]), "\n")
cat("bsadf[,1] (series a):", sprintf("%.8f", rfd$bsadf[, 1]), "\n")
cat("bsadf[,2] (series b):", sprintf("%.8f", rfd$bsadf[, 2]), "\n")
cat("bsadf_panel:", sprintf("%.8f", rfd$bsadf_panel), "\n")
cat("minw:", attr(rfd, "minw"), " lag:", attr(rfd, "lag"), " n:", attr(rfd, "n"), "\n\n")

cat("=== tidy(rfd) [wide] ===\n")
print(tidy(rfd))
cat("=== tidy(rfd, format='long') ===\n")
print(tidy(rfd, format = "long"))
cat("=== tidy(rfd, panel=TRUE) ===\n")
print(tidy(rfd, panel = TRUE))
cat("=== tidy(rfd, panel=TRUE, format='long') ===\n")
print(tidy(rfd, panel = TRUE, format = "long"))

cat("=== augment(rfd) [wide, trunc=TRUE] -- first 4 rows ===\n")
print(head(augment(rfd), 4))
cat("nrow:", nrow(augment(rfd)), "(expect (n - minw - lag) * ncol =", (30 - 10 - 1) * 2, ")\n")
cat("=== augment(rfd, format='long') -- first 4 rows ===\n")
print(head(augment(rfd, format = "long"), 4))
cat("=== augment(rfd, panel=TRUE) -- first 3 rows ===\n")
print(head(augment(rfd, panel = TRUE), 3))
cat("=== augment(rfd, trunc=FALSE) -- first 4 rows (pre-minw+lag: NA) ===\n")
print(head(augment(rfd, trunc = FALSE), 4))

# pyexuber's tidy()/augment() are exercised against these exact printed
# numbers by docs/replication/core-workflow/tidy_validation.py (via a
# hand-built RadfResult carrying the same values -- pyexuber's radf()
# itself can't be run outside CI, see pyexuber/CLAUDE.md), and by
# pyexuber/tests/test_tidy.py.
