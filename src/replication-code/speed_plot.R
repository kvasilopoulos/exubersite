# Figure comparing exuber::radf() in R and Python (pyexuber), from the two
# tables written by speed_r.R and speed_py.py. Both use minw = 30, lag = 1 and
# a random walk of length n. The y axis is on a log scale. Run from the
# project root:
#   Rscript docs/replication/speed/speed_plot.R
suppressMessages({
  library(ggplot2)
})

dir <- file.path("docs", "replication", "speed")
tbl <- rbind(
  read.csv(file.path(dir, "speed-r.csv")),
  read.csv(file.path(dir, "speed-py.csv"))
)
# Slowest (top of legend) to fastest (bottom), ranked by the n = 1000 time.
at_max <- tbl[tbl$n == max(tbl$n), ]
tbl$software <- factor(tbl$software, levels = at_max$software[order(-at_max$time_ms)])

p <- ggplot(tbl, aes(n, time_ms, colour = software, linetype = software)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(100, 1000, 100)) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x = "Sample size", y = "Milliseconds (log scale)", colour = NULL, linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.09, 0.86),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = "grey80"),
    legend.key.width = unit(1.4, "lines")
  )

out <- file.path(dir, "speed-r-vs-python.png")
ggsave(out, p, width = 8, height = 5.5, dpi = 120)
cat("written:", out, "\n")
print(reshape(tbl, idvar = "n", timevar = "software", direction = "wide"), digits = 4)
