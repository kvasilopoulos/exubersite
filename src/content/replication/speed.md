---
title: "Speed"
blurb: "Computation time of radf() across software, and of the R and Python implementations of exuber on the same grid."
order: 7
---
This page reports how long `radf()` takes. The first comparison puts exuber next to the other software that computes the same recursive tests. The second compares the R and the Python implementation of exuber on identical settings.

## exuber against other software

The figure reproduces the comparison in Section 4 of Vasilopoulos, Pavlidis and Martinez-Garcia (2022, *Journal of Statistical Software*). It covers the R packages MultipleBubbles and psymonitor, EViews (`rtadf`), MATLAB (`PSY.m`) and Stata, together with exuber as it was benchmarked in the paper (0.4.1) and the current version. Every run uses `minw = 30`, `lag = 1` and a random walk of length `n`, and each point is the median elapsed time over repeated runs. The vertical axis is on a log scale.

![Elapsed time of radf() and of the alternative software, by sample size](/speed/benchmark-plot-1.png)

All series except exuber 2.0.0 are the archived runs from the paper, so they do not depend on the current implementation. MultipleBubbles and psymonitor are `O(T^2)` loops written in pure R and already take minutes per run at `n = 1000`. exuber 2.0.0 is faster than 0.4.1 at every sample size shown. The speed comes from the recursive least-squares algorithm in `radf()`, which uses the matrix inversion lemma and never inverts a matrix for each window.

## R and Python

The R package and pyexuber call the same C++ core (exubercore v0.3.1), so the statistic costs about the same in either language. We timed both on one machine with identical settings: `minw = 30`, `lag = 1`, a random walk of length `n = 100, 200, ..., 1000`, and the median of five block medians of 60 calls each. R is exuber 2.0.0 installed with `R CMD INSTALL`, since a development load compiles without release optimisation flags. Python is pyexuber 0.1.0 from PyPI.

![Elapsed time of radf() in R and in Python, by sample size](/speed/speed-r-vs-python.png)

| n | exuber 2.0.0 (R), ms | pyexuber 0.1.0 (Python), ms |
|---|---|---|
| 100 | 2.0 | 0.9 |
| 200 | 3.5 | 2.8 |
| 300 | 5.9 | 5.9 |
| 400 | 12.3 | 11.5 |
| 500 | 17.1 | 18.4 |
| 600 | 23.1 | 30.0 |
| 700 | 38.1 | 44.8 |
| 800 | 45.4 | 54.1 |
| 900 | 59.7 | 57.0 |
| 1000 | 66.7 | 95.6 |

Neither implementation is consistently faster. Python is quicker at the smallest sizes, where the call overhead is a larger share of the time, and the two cross several times after that. The ranking at `n = 1000` changed between repeated runs on this machine (Python took between 49 and 96 ms across three runs, R between 55 and 67 ms), so differences of this size should not be read as a property of either language. The difference that matters for users is outside `radf()`. The Monte Carlo and bootstrap critical values repeat the call many times, and there the setup of the loop and the number of parallel workers decide the time, not the language.

The figure and the table come from three scripts that you can rerun: [speed_r.R](#script-speed_r) times the R package, `speed_py.py` times pyexuber, and [speed_plot.R](#script-speed_plot) draws the figure from the two tables they write.
