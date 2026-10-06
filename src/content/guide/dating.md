---
title: "Dating and root inference"
blurb: "Date a bubble you already believe is there, and estimate how fast it grows."
group: "Beyond PSY"
order: 2
family: "dating-and-root-inference"
---
## Dating Methods: Alternatives to datestamp()

``` r
library(exuber)
```

### Why these exist alongside `datestamp()`

`datestamp()` applies the rule of Phillips, Shi and Yu to a `radf()` result: a bubble runs from the first point at which the recursive statistic crosses its critical value to the first point at which it falls back below it. The rule is simple and well understood, but it is not the only way to date a bubble once you already believe that one exists. The `dating_*()` family fits an explicit regime model to the raw series: a unit root, then an explosive regime, then a unit root again. It chooses the break dates that minimize the residual sum of squares (SSR). These functions need no critical value. Given a window that contains at most one bubble, they tell you where the bubble starts and ends, and they do not tell you whether there is one.

| Function | Paper | Idea |
|---|---|---|
| `dating_hls()` | Harvey, Leybourne & Sollis (2017) | Fits four candidate regime-dummy models (with or without a distinct collapse regime, with or without recovery) using closed-form segment SSR, and the BIC picks among them. |
| `dating_knp()` | Kejriwal, Nguyen & Perron (2025) | Uses the same model as Model 2 of HLS. The authors prove that the plain SSR minimizer is inconsistent, because it converges to the collapse date and not to the origination date, and they correct this by omitting one squared residual from the objective. |
| `dating_pdc()` | Pang, Du & Chong (2021); Kurozumi & Skrobotov (2023) | Assumes a fixed structure of three or four regimes and finds each breakpoint sequentially in closed form, starting with the collapse because it is stochastically dominant. There is no BIC step. |
| `dating_hlw()` | Harvey, Leybourne & Whitehouse (2020) | A wrapper that first runs `radf()` and `datestamp()` to find how many episodes there are and roughly where they lie, and then applies HLS-style fitting separately within each detected window. |

None of the four takes a `radf_cv` object, so they do not work with `summary()`, `tidy()` or `autoplot()`. Each prints its own dating table instead. [The results, tidying and plotting page](/guide/pipeline) describes the full pipeline.

### A single bubble, four estimates

We use one series simulated with `sim_ps1()`, the single-bubble data generating process of Phillips & Shi (2018). It has a unit-root run-up, an explosive regime that starts at 40, a mildly integrated collapse regime that starts at 61, and a return to a unit root at 70. This is the regime structure that these estimators are built for.

``` r
y <- sim_ps1(n = 100, seed = 1)
```

The true origination date is 40 and the true collapse date is 61. We run all four estimators:

``` r
dating_hls(y, trim = 0.05)
#> 
#> ── dating_hls (n = 100, trim = 0.05) ───────────────────────────────────────────
#> 
#>    series  model  origination  collapse  recovery
#>   series1      4           39        60        70
```

``` r
dating_knp(y, trim = 0.05)
#> 
#> ── dating_knp (n = 100, trim = 0.05, omit = TRUE, breaks = 2 ───────────────────
#> 
#>    series  bubble  origination  collapse   delta
#>   series1       1           60        70  0.9178
```

``` r
dating_pdc(y, regimes = 3, trim = 0.05)
#> 
#> ── dating_pdc (n = 100, regimes = 3, type = ols) ───────────────────────────────
#> 
#>    series  origination  collapse
#>   series1           38        59
```

``` r
dating_hlw(y, trim = 0.1, nboot = 199, seed = 1)
#> 
#> ── dating_hlw (n = 100, trim = 0.1) ────────────────────────────────────────────
#> 
#> series1:
#>  model origination collapse recovery
#>      4          39       60       70
```

On this draw, `dating_hls()` and `dating_hlw()` both select the four-regime model and land within one period of every true date (39, 60 and 70 against 40, 61 and 70). `dating_pdc()` is one or two periods early on both dates (38 and 59). `dating_knp()` performs worst here. Its model assumes an instantaneous collapse, with the unit root resuming from a shifted level, so it has no room for the ten-period mildly integrated crash (61 to 70) that `sim_ps1()` generates. It dates that crash as the episode instead and reports 60 and 70. This follows from the mismatch between the model and the data and not from the estimator's bias correction.

We run all four side by side to show that they are different estimators with different failure modes, and not to single out one as correct. When they disagree on real data, the disagreement is informative, and it should not be settled by picking a favorite.

### Which to reach for

- If you believe the window contains exactly one bubble and you want the best-fitting regime model, with or without a distinct collapse and recovery regime, use `dating_hls()`.
- If you have the same setting but the accuracy of the origination date matters more than that of the collapse date, use `dating_knp()`. Its authors find that the naive estimator is biased toward the collapse date.
- If you want closed-form estimates with no BIC model search and can fix the number of regimes in advance, use `dating_pdc()`. Add `weights` for the volatility-corrected variant.
- If you do not know how many episodes there are, or you want the dating to start from an actual `radf()` and `datestamp()` detection, use `dating_hlw()`.

## Root Inference: How Fast Is the Bubble Growing

``` r
library(exuber)
```

### A different question from "is there a bubble"

`radf()`, `datestamp()` and the `dating_*()`, `_test()`, `monitor()` and `monitor_*()` families all answer some version of the question of whether there is a bubble and when it happened. None of them says anything about its magnitude. Once an explosive episode is dated, how fast is the underlying autoregressive root growing? `rootstamp()` (Phillips & Magdalinos 2007; Guo, Sun & Wang 2019) answers that question. It is a follow-up step after detection and dating, and it replaces neither.

The function fits a no-intercept AR(1), `y_t = rho * y_{t-1} + e_t`, over a given sub-sample. It reports the estimate of `rho` with a confidence interval and the implied doubling time, `log(2) / log(rho)`, which is the number of periods the bubble needs to double in size at the estimated growth rate. There are two methods for two starting points. The default method takes a numeric sub-sample and fits it once, with one confidence interval. The `radf_obj` method takes a `radf_obj` together with its `datestamp()` result and fits every episode at once, with no manual loop. Neither method returns the `radf_obj` class, so `rootstamp()` does not work with `summary()`, `tidy()` or `autoplot()` (see [the results, tidying and plotting page](/guide/pipeline)).

### Detect, date, then estimate the root

The series has a unit-root run-up followed by an explosive regime with `rho = 1.04`:

``` r
y <- sim_psy1(n = 100, te = 60, tf = 100, c = 0.04, alpha = 0, sigma = 1, seed = 2026)
```

We first detect and date the episode in the usual way:

``` r
r <- radf(y, minw = 20)
cv <- radf_mc_cv(length(y), minw = 20, nrep = 300, seed = 4)
ds <- datestamp(r, cv = cv, min_duration = 3)
ds
#> 
#> ── Datestamp (min_duration = 3) ───────────────────────────────── Monte Carlo ──
#> 
#> series1 :
#>   Start Peak End Duration   Signal Ongoing
#> 1    63  100 100       38 positive    TRUE
```

Then we estimate the root over the detected episode. The default method takes the sub-sample directly, which we slice with the `Start` and `End` of the episode:

``` r
ep <- ds[["series1"]]
rootstamp(y[ep$Start[1]:ep$End[1]]) # normal-t interval (Guo, Sun & Wang 2019), true rho = 1.04
#> 
#> ── rootstamp (n = 37, sig_lvl = 95%, type = normal) ────────────────────────────
#> 
#>    rho         se  t_stat  rho_lower  rho_upper  doubling_time  dt_lower
#>   1.04  0.0007295    54.2      1.038      1.041          17.87     17.26
#>   dt_upper
#>      18.53
rootstamp(y[ep$Start[1]:ep$End[1]], type = "cauchy") # fixed-root Cauchy interval (Phillips & Magdalinos 2007)
#> 
#> ── rootstamp (n = 37, sig_lvl = 95%, type = cauchy) ────────────────────────────
#> 
#>    rho         se  t_stat  rho_lower  rho_upper  doubling_time  dt_lower
#>   1.04  0.0007295    54.2     0.7955      1.284          17.87     2.776
#>   dt_upper
#>      -3.03
```

The estimate of `rho` is close to the true value of 1.04. The output also reports `rho_ci` and the implied `doubling_time` and `doubling_time_ci`. The two interval types answer slightly different questions. `type = "normal"`, the default, is the safer choice under drift or weak dependence, and it gives a noticeably tighter interval here. `type = "cauchy"` assumes a fixed root that does not drift. A Cauchy distribution has much fatter tails than a normal one, so this interval is visibly wider even at the same nominal level.

### Every episode at once

When a `datestamp()` result contains more than one episode, the `radf_obj` method runs the default method on each of them without a manual loop. Pass the original `radf()` result and the `datestamp()` result together:

``` r
rootstamp(r, ds)
#> 
#> ── rootstamp (sig_lvl = 95%, type = normal) ────────────────────────────────────
#> 
#> series1 :
#>   Start End  rho rho_lower rho_upper doubling_time doubling_time_lower
#> 1    63 100 1.04     1.038     1.041         17.87               17.26
#>   doubling_time_upper
#> 1               18.53
```

Root inference on a very short episode is close to meaningless, because there are too few points to estimate an AR(1) coefficient precisely. Filter the episodes with `datestamp(..., min_duration = ...)` before passing them in. The method does not decide for you what counts as too short.

## Experimental Methods: radf_recovery() and datestamp(option = 'svadf')

``` r
library(exuber)
```

### What "experimental" means here

Most methods in exuber implement the procedure of a peer-reviewed paper and pass the package's standard validation. That validation consists of a formula-exact check against a brute-force reimplementation, a lookup against published tables, a Monte Carlo check of size, and a check of power against a true alternative. `radf_recovery()` and `datestamp(option = "svadf")` went through the same validation and both give useful results, but each has one disclosed gap that keeps it below the standard. For that reason they print an "Experimental" badge and emit a caveat message when called. Treat their output as a guide to where episodes lie, and do not assume it is as well calibrated as the rest of the package.

### `radf_recovery()`: dating a collapse and a recovery

This function uses the reverse-regression idea of Phillips & Shi (2014). We reverse the series in time, run the BSADF recursion that `radf()` already computes, and map the crossing dates back to the original time axis. In the reversed series a collapse followed by a recovery turns the collapse into an explosive regime and the recovery into the end of that regime, so the forward machinery run backwards dates both.

``` r
# sim_ps1(): unit root -> explosive (40-60) -> collapse (61-70) -> recovery (71+)
y <- sim_ps1(n = 100, seed = 2)
res <- radf_recovery(y, minw = 15, nrep = 200, seed = 1)
res
#> 
#> ── radf_recovery (n = 100, minw = 15, level = 95%) ─────────────────────────────
#> 
#> ℹ Experimental. f_c and the overall false-detection rate are exploratory pending further validation; see ?radf_recovery, Caveats section.
#> 
#>    series  f_c  f_r  detected  censored
#>   series1   57   67      TRUE     FALSE
```

The estimate `f_c` (crisis onset, 57) falls just before the true collapse start (61). The estimate `f_r` (recovery, 67) falls after it, inside the collapse regime and before the true recovery date (70). The two dates come out in the right order by construction, because the down-crossing search only starts at the up-crossing. The disclosed gap is that `f_c` and the overall false-detection rate are exploratory until they are validated further (see the Caveats section of `?radf_recovery`). The ordering of the dates is reliable, but the calibration of false alarms is not yet.

### `datestamp(option = "svadf")`: a preprint

This option implements Sarkar & Wells (2026), an arXiv preprint that has not been peer reviewed. Every other paper implemented in the package has been, so the evidence behind this method is weaker. Its statistic is the `badf` sequence that `radf()` already computes, compared against two closed-form thresholds that depend only on the sample size and come from the applied methodology of the paper. It is an option of `datestamp()` and not a separate function.

``` r
res <- radf(sim_data, lag = 0)
datestamp(res, option = "svadf", min_duration = psy_ds(nrow(sim_data)))
#> 
#> ── Datestamp (min_duration = 5) ──────────────── SV-ADF (Sarkar & Wells 2026) ──
#> 
#> ℹ Experimental. Sarkar & Wells (2026) is a non-peer-reviewed preprint; see ?datestamp, Caveats section.
#> 
#> psy1 :
#>   Start Peak End Duration   Signal Ongoing
#> 1    48   48  49        1 positive   FALSE
#> 
#> psy2 :
#>   Start Peak End Duration   Signal Ongoing
#> 1    23   23  24        1 positive   FALSE
```

`psy1` and `psy2` receive clear origination and collapse dates, while `evans`, `div` and `blan` never cross the threshold in this panel.

### Using them responsibly

Both methods are worth using. The date ordering from `radf_recovery()` and the point statistic from `datestamp(option = "svadf")` are reliable. Neither should be the only basis for a claim about false-alarm rates or exact calibration, though. When that matters, prefer `radf()` and `datestamp()`, or one of the peer-reviewed alternatives in [the alternative paradigms page](/guide/alternatives) and [the dating and root inference page](/guide/dating). Treat these two methods as a second opinion until their caveats are resolved. The caveats are listed in the replication notes for dating and volatility-robustness.
