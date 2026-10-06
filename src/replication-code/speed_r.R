# Timing of exuber::radf() in R: median elapsed time per call on a random walk
# of length n (median of five block medians), with minw = 30 and lag = 1, for n = 100, 200, ..., 1000.
# speed_py.py runs the same grid with the same settings through pyexuber, and
# speed_plot.R draws both series.
#
# Time the installed package, not devtools::load_all(). A development load
# compiles without the optimisation flags of a release build and would
# understate the speed of R. From the project root:
#   R CMD INSTALL exuber
#   Rscript docs/replication/speed/speed_r.R
suppressMessages({
  library(exuber)
  library(microbenchmark)
})
options(exuber.show_progress = FALSE, exuber.parallel = FALSE)

minw <- 30L
lag <- 1L
sample_size <- seq(100L, 1000L, by = 100L)
times <- 60L # calls per block
blocks <- 5L # the reported time is the median of the block medians

res <- do.call(rbind, lapply(sample_size, function(n) {
  set.seed(123)
  rw <- cumsum(rnorm(n))
  radf(rw, minw = minw, lag = lag) # warm-up call
  t_ms <- median(replicate(blocks, {
    median(microbenchmark(
      radf(rw, minw = minw, lag = lag),
      times = times
    )$time) / 1e6
  }))
  cat(sprintf("n = %4d   R exuber %s   %9.3f ms\n", n, packageVersion("exuber"), t_ms))
  data.frame(software = paste0("exuber ", packageVersion("exuber"), " (R)"), n = n, time_ms = t_ms)
}))

out <- file.path("docs", "replication", "speed", "speed-r.csv")
write.csv(res, out, row.names = FALSE)
cat("written:", out, "\n")
