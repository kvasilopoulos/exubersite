# Validation script for radf_hlw() (Harvey, Leybourne & Whitehouse 2020,
# "Date-stamping multiple bubble regimes" -- the two-step wrapper around
# radf_hls() that extends single-bubble SSR/BIC dating to series with
# more than one explosive episode). See docs/enhancements/
# dating-and-root-inference.md, "SSR/BIC dating vs. PSY recursive
# dating", for the full write-up. Run from the exuber/ package root (or
# adjust the devtools::load_all() path below).

Sys.setenv(NOT_CRAN = "true")
options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
devtools::load_all("exuber", quiet = TRUE)

cat("=== 1. hlw_local_to_global() arithmetic (window-local breakpoint\n")
cat("     -> global i-index / date position) ===\n\n")
g <- exuber:::hlw_local_to_global(local_tau = 5L, s = 21L)
cat("local_tau=5, s=21 -> i_index=", g$i_index, " position=", g$position,
    " (expect 25, 26)\n\n")

cat("=== 2. Two genuine bubbles: window count + breakpoint accuracy (20 reps) ===\n")
cat("(datestamp()'s own step-1 PSY detection can fragment a true episode\n")
cat("into extra spurious windows under noise -- this is expected and is\n")
cat("HLW's own documented issue, not a radf_hlw() bug; we report the full\n")
cat("window-count distribution rather than hiding the fragmentation)\n\n")
sim_two_bubbles <- function(seed, n1a = 50, n2a = 20, n3a = 30, n1b = 50, n2b = 20, n3b = 30) {
  set.seed(seed)
  e1 <- 100 + cumsum(rnorm(n1a))
  b1 <- e1[n1a] * 1.05^(1:n2a) + cumsum(rnorm(n2a))
  u1 <- b1[n2a] + cumsum(rnorm(n3a))
  e2 <- u1[n3a] + cumsum(rnorm(n1b))
  b2 <- e2[n1b] * 1.05^(1:n2b) + cumsum(rnorm(n2b))
  u2 <- b2[n2b] + cumsum(rnorm(n3b))
  y <- c(e1, b1, u1, e2, b2, u2)
  list(y = y, true1 = c(n1a, n1a + n2a), true2 = c(n1a + n2a + n3a + n1b, n1a + n2a + n3a + n1b + n2b))
}
run_two <- function(seed) {
  sim <- sim_two_bubbles(seed)
  out <- radf_hlw(sim$y, trim = 0.1, min_duration = psy_ds(length(sim$y)), nboot = 199, seed = 1)
  df <- out[["series1"]]
  list(n_windows = nrow(df), df = df, true1 = sim$true1, true2 = sim$true2)
}
res_two <- lapply(1:20, run_two)
n_windows <- sapply(res_two, `[[`, "n_windows")
cat("Window count distribution:", paste(names(table(n_windows)), table(n_windows), sep = "=", collapse = ", "), "\n")
two_win <- res_two[n_windows == 2]
cat(sprintf("Reps with exactly 2 windows: %d/20\n", length(two_win)))
if (length(two_win) > 0) {
  orig1_bias <- sapply(two_win, function(r) as.numeric(r$df$origination[1]) - r$true1[1])
  coll1_bias <- sapply(two_win, function(r) as.numeric(r$df$collapse[1]) - r$true1[2])
  orig2_bias <- sapply(two_win, function(r) as.numeric(r$df$origination[2]) - r$true2[1])
  coll2_bias <- sapply(two_win, function(r) as.numeric(r$df$collapse[2]) - r$true2[2])
  cat(sprintf("Bubble 1: origination mean|bias|=%.2f, collapse mean|bias|=%.2f\n",
              mean(abs(orig1_bias)), mean(abs(coll1_bias))))
  cat(sprintf("Bubble 2: origination mean|bias|=%.2f, collapse mean|bias|=%.2f\n\n",
              mean(abs(orig2_bias)), mean(abs(coll2_bias))))
}

cat("=== 3. Windows are always ordered (never overlapping) ===\n")
ok <- sapply(res_two, function(r) {
  if (nrow(r$df) < 2) return(TRUE)
  all(diff(as.numeric(r$df$origination)) > 0)
})
cat("all ordered:", all(ok), "\n\n")

cat("=== 4. Pure H0 (no bubble at all): no error, and how often 0 windows ===\n")
run_h0 <- function(seed) {
  set.seed(seed)
  y0 <- 100 + cumsum(rnorm(150))
  out <- radf_hlw(y0, trim = 0.1, nboot = 199, seed = 1)
  nrow(out[["series1"]])
}
n0 <- sapply(1:20, run_h0)
cat("Window counts under H0:", paste(names(table(n0)), table(n0), sep = "=", collapse = ", "), "\n\n")

cat("=== 5. Single clean bubble: final window vs standalone radf_hls() ===\n")
run_single <- function(seed, n1 = 60, n2 = 25, n3 = 25, n4 = 40) {
  set.seed(seed)
  unit1 <- 100 + cumsum(rnorm(n1))
  bubble <- unit1[n1] * 1.05^(1:n2) + cumsum(rnorm(n2))
  target <- bubble[n2] * 0.5
  collapse <- numeric(n3)
  collapse[1] <- bubble[n2] + rnorm(1)
  for (k in 2:n3) collapse[k] <- target + 0.85 * (collapse[k - 1] - target) + rnorm(1)
  recovery <- collapse[n3] + cumsum(rnorm(n4))
  y <- c(unit1, bubble, collapse, recovery)
  hls_out <- radf_hls(y, trim = 0.1)
  hlw_out <- radf_hlw(y, trim = 0.1, min_duration = psy_ds(length(y)), nboot = 199, seed = 1)
  df <- hlw_out[["series1"]]
  last <- df[nrow(df), ]
  list(
    n_windows = nrow(df),
    match = identical(unname(hls_out$model[["series1"]]), last$model) &&
      identical(unname(hls_out$origination[["series1"]]), last$origination) &&
      identical(unname(hls_out$collapse[["series1"]]), last$collapse)
  )
}
res_single <- lapply(1:15, run_single)
cat(sprintf("Final-window match rate vs standalone radf_hls(): %.2f\n",
            mean(sapply(res_single, `[[`, "match"))))
cat("Window count distribution:",
    paste(names(table(sapply(res_single, `[[`, "n_windows"))),
          table(sapply(res_single, `[[`, "n_windows")), sep = "=", collapse = ", "), "\n")
