---
title: "Volatility-robust tests"
blurb: "Tests robust to time-varying innovation variance: time-transformed, kernel-purged, WLS, sign-based, and stochastic-coefficient routes."
order: 1
---
The GSADF test of Phillips, Wu & Yu (2011) and Phillips, Shi & Yu (2015) assumes that the innovation variance is constant. When the variance moves over time, whether deterministically or stochastically, the test over-rejects. Every method in this file is a different remedy for that problem. The status column uses these labels. `done` means the method is implemented and checked against a published number. `evaluated` means we have read the source and not implemented the method.

| Method | Paper | Status |
|---|---|---|
| [Time-transformed test (STADF/GSTADF)](#time-transformed-test-stadf--gstadf) | Kurozumi, Skrobotov & Tsarev (2024) | done |
| [Kernel-purge test](#kernel-purge-test) | Harvey, Leybourne, Taylor & Zu (2024/2025) | done (with-intercept variant) |
| [SBZ (WLS + kernel volatility)](#sbz-wls--kernel-volatility) | Harvey, Leybourne & Zu (2019) | done |
| [Sieve bootstrap (autocorrelated innovations)](#pedersen--schütte-sieve-bootstrap) | Pedersen & Montes Schütte (2020) | done |
| [Skewness-corrected wild bootstrap](#hafner-skewness-corrected-wild-bootstrap) | Hafner (2020) | done |
| [Sign-based sGSADF](#sign-based-sgsadf) | Harvey, Leybourne & Zu (2020); level shifts: Harvey, Leybourne, Tatlow & Zu (2025) | done |
| [Stochastic explosive-coefficient test](#stochastic-explosive-coefficient-test) | Kurozumi & Nishi (2025) | done (`ssu_test()`, `cusum_test()`) |
| [SV-ADF](#sv-adf) | Sarkar & Wells (2026, preprint) | done (`datestamp(option = "svadf")`; preprint, not peer-reviewed) |

All papers are listed in [references.md](/replication/references#volatility-robustness).

---

## Time-transformed test (STADF / GSTADF)

Status: done. The test is in `exuber/R/radf_tt.R` (`radf_tt()`, `radf_tt_cv()`) and is checked in `exuber/tests/testthat/test-tt.R`.

### Source

Kurozumi, E., Skrobotov, A. & Tsarev, A. (2024). Time-Transformed Test for Bubbles under Non-stationary Volatility. *Journal of Financial Econometrics*, `doi:10.1093/jjfinec/nbae026`. We worked from the arXiv version (2012.13937v2, 15 Nov 2021), and equation and theorem numbers below refer to it.

### Idea

If the variance path were known, we could stretch calendar time so that the rescaled series has constant variance, and then run an ordinary GLS-demeaned SADF or GSADF test on it. Theorem 1 shows that the null limit of the statistic computed on the time-transformed series coincides with the homoskedastic GLS-demeaned SADF/GSADF distribution of Whitehouse (2019), whatever the volatility path. Theorem 2 shows that this still holds when the variance profile is estimated from the data. The null distribution is therefore pivotal, and `radf_tt_cv()` simulates critical values once from a plain random walk, with no bootstrap and no dependence on the dataset.

The paper reports one anchor value. Footnote 4 gives, for $r_0 = 0.1$, critical values $2.319$, $2.626$ and $3.223$ at the 10%, 5% and 1% levels. These belong to STADF, the single-supremum statistic $\sup_{r_2} \mathrm{ADF}^{r_2}_0$, and not to the double-supremum GSTADF. Whitehouse (2019) studies the GLS-demeaned version of the PWY (2011) test, which is the single-supremum one. GSTADF critical values are obtained with `radf_tt_cv()`, because the paper gives no published value for them.

### Formulas

Let $\check y_t = y_t - y_0$ be the GLS-demeaned series. We subtract the first observation and fit no intercept. This differs from `radf()`, which always fits an intercept.

The recursive statistic of eq. (9) is the $t$-statistic on $\beta$ in the no-intercept regression $\Delta \check y_t = \beta\, \check y_{t-1} + e_t$ over the window $[r_1, r_2]$:

$$
\mathrm{ADF}^{r_2}_{r_1} = \frac{\sum \check y_{t-1}\,\Delta \check y_t}{\sqrt{\hat\sigma^2(r_1, r_2) \sum \check y_{t-1}^2}} .
$$

`gls_dfstat_grid()` computes it over the whole $(r_1, r_2)$ grid from cumulative sums, so no per-window regression is run.

The variance profile is estimated in three steps.

- Eq. (18) is a local kernel (Nadaraya–Watson type) estimate $\hat\delta_t$ of the time-varying AR(1) coefficient. The default kernel is uniform, as in the paper's Monte Carlo.
- Eq. (19) defines the variance profile $\hat\eta(s)$ as a normalized cumulative sum of squares of the truncated local-regression residuals. It is piecewise linear on the observation grid, so we invert it exactly by linear interpolation to get $\hat g(s)$, which is used to resample and time-transform the series.
- Footnote 6 sets the truncation threshold to $\psi_T = \bar c\, T^{1/7}$, where $\bar c$ is the largest residual standard deviation over rolling windows of 10% of the sample.

### Choices that differ from the paper

- **Bandwidth.** The paper selects $h \in [T^{-0.5}, T^{-0.3}]$ by leave-one-out cross-validation. The default here is the plug-in $h = T^{-2/5}$, the midpoint of that range on a log scale, and `h` is a user argument. Cross-validation would repeat the $O(T \cdot Th)$ kernel fit about ten times. It is the first thing to add if empirical size control turns out to matter.
- **Lags.** The paper's Monte Carlo uses no augmentation lags, and `radf_tt()` likewise supports only $p = 0$.

### Why the statistic is not computed in exubercore

The recursive statistic is a closed-form ratio of cumulative sums, so `gls_dfstat_grid()` is fully vectorized in R. The C++ `radf()` always fits an intercept when `lag == 0` (see `exubercore/src/radf.cpp`), so it could not compute this no-intercept statistic without modification. A port becomes worthwhile only if the $O(T^2)$ grid is a bottleneck for $T$ in the thousands. Runs up to $T \approx 2000$ showed no such problem.

### Validation

1. **Formula.** `gls_dfstat_grid()` matches a brute-force per-window `lm(dy ~ ylag - 1)` fit on a series of length 65 with `minw = 15`, with a maximum absolute difference of $6.7 \times 10^{-16}$.
2. **Critical values.** A Monte Carlo with $n = 300$, `minw = 30` and 4000 replications gives the STADF critical values below. The deviation from the published triple has the same direction and size as in the package test, which uses a tolerance of 0.15. It is consistent with Monte Carlo noise around a $T \to \infty$ target. `radf_tt_cv(n = 300, minw = 30, nrep = 4000)` gives $(2.399, 2.737, 3.346)$ on another seed.

   | | 10% | 5% | 1% |
   |---|---|---|---|
   | Published (Whitehouse 2019) | 2.319 | 2.626 | 3.223 |
   | Simulated | 2.355 | 2.715 | 3.435 |
   | Absolute difference | 0.036 | 0.089 | 0.212 |

   For comparison, the GSTADF (double-supremum) statistic under the same setup has critical values $(3.157, 3.436, 4.302)$.
3. **Power.** On a series that follows a unit root, turns explosive ($\rho = 1.04$) and then collapses, `radf_tt()` gives a GSADF statistic of 5.781, well above every critical value in the table.

Replication script: [replication/volatility-robustness/radf_tt_validation.R](#script-radf_tt_validation).

### `datestamp()` and `autoplot()` support (done)

`radf_tt_cv()` returns the scalar critical values `adf_cv`, `sadf_cv` and `gsadf_cv`, and also the time-varying boundaries `badf_cv` and `bsadf_cv` that `datestamp()` and `autoplot()` need. `gls_dfstat_grid()` already returns the supremum over all window starts at each end point, so the boundaries are the per-time-point quantiles of `badf` and `bsadf` across replications. This is the same construction as in `radf_mc_cv()`.

Checks:

1. The last row of `badf_cv` is identical to `adf_cv`, as it must be, since `adf` is the last point of `badf` in each replication.
2. Under a pure random walk ($n = 100$, `minw = 20`, 2000 critical-value replications, 300 test replications) the false-alarm rate is 3.3% at a nominal 5%, for both `option = "gsadf"` and `option = "sadf"`.
3. On a short, moderate bubble (60 pre-bubble observations, 30 explosive observations at $\rho = 1.03$, 10 post-collapse observations, 50 replications), `datestamp()` on `radf_tt()` detects the bubble in 18% of runs, against 16% for the `radf()` and `radf_mc_cv()` pipeline on the same series.
4. On the bundled `psy2` series, `datestamp()` finds two episodes (Start 21 / Peak 27 / End 35, and Start 55 / Peak 55 / End 73).

`radf_sign_cv()` and `radf_sign_dm_cv()` share the construction and are described in [Sign-based sGSADF](#sign-based-sgsadf).

### Open items

- Provide critical values over a grid of $r_0$. Only $r_0 = 0.1$ is anchored to a published value. The authors point to R code at https://sites.google.com/site/antonskrobotov/.

---

## Kernel-purge test

Status: done (with-intercept variant). The test is in `exuber/R/radf_kp.R` (`radf_kp()`) and is checked in `exuber/tests/testthat/test-kp.R`.

### Source

Harvey, D. I., Leybourne, S. J., Taylor, A. M. R. & Zu, Y. (2024). A new heteroskedasticity-robust test for explosive bubbles. *Journal of Time Series Analysis*, `doi:10.1111/jtsa.12784`. Open access.

### Idea

Instead of resampling (`radf_wb_cv()`) or deforming time (`radf_tt()`), the test removes the unconditional heteroskedasticity directly. It estimates the spot volatility $\hat\sigma_t$ with a Gaussian kernel (eq. 4), divides each first difference by it, and cumulates the result (eq. 5) into a volatility-standardized series $x_t$. The unmodified PSY/GSADF test is then run on $x_t$.

Theorem 1 and Remark 3.2 show that the null limit of the purged statistic equals the standard homoskedastic GSADF limit. The Monte Carlo critical values of `radf_mc_cv()` therefore apply directly, and `tidy()`, `autoplot()` and `datestamp()` work on the result. `radf_kp()` calls `kernel_spot_vol()`, which is shared with [SBZ](#sbz-wls--kernel-volatility), and then the unmodified `radf()` on the purged series.

The paper also proposes a without-intercept variant $\mathrm{PSY}^*_\sigma$ and a union-of-rejections test $\mathrm{UPSY}_\sigma$ that combines the two. Neither is implemented. The without-intercept variant needs a no-intercept regression on non-demeaned data, which differs slightly from the GLS-demeaned `gls_dfstat_grid()`. The union scaling has the same form as the union in `radf_sbz_cv()`.

### Published critical values

Table I gives critical values for a minimum window of 0.1, a Gaussian kernel, $h = 0.1\,T^{-0.25}$ and 2000 replications.

| $T$ | $\mathrm{PSY}_\sigma$ (10% / 5% / 1%) | $\mathrm{PSY}^*_\sigma$ (10% / 5% / 1%) | $\mathrm{UPSY}_\sigma$ (10% / 5% / 1%) |
|---|---|---|---|
| 100 | 1.629 / 1.828 / 2.392 | 3.637 / 4.158 / 5.553 | 3.950 / 4.527 / 6.129 |
| 200 | 1.608 / 1.789 / 2.140 | 3.226 / 3.595 / 4.330 | 3.468 / 3.804 / 4.589 |
| 400 | 1.712 / 1.935 / 2.296 | 3.167 / 3.446 / 4.007 | 3.361 / 3.598 / 4.145 |
| $\infty$ | 1.875 / 2.094 / 2.486 | 2.978 / 3.296 / 3.859 | 3.186 / 3.486 / 3.951 |

The implemented statistic is the $\mathrm{PSY}_\sigma$ column. By Remark 3.2 its $T = \infty$ row should equal the asymptotic GSADF critical values of PSY (2015) at the same minimum window.

### Validation

`test-kp.R` simulates the null GSADF distribution of `radf_kp()` at $n = 300$ and compares it with `radf_mc_cv(300)` and with the published $T = 400$ row. One run gave $(1.685, 1.897, 2.337)$ for `radf_kp`, against $(1.873, 2.121, 2.429)$ for `radf_mc_cv(300)` and $(1.712, 1.935, 2.296)$ in the table. Dividing by an estimated volatility path, and not the true one, plausibly shrinks the effective variance a little in finite samples.

A run at $n = 400$ (the table row, so no interpolation), 800 replications, gives:

| | 10% | 5% | 1% |
|---|---|---|---|
| Published (Table I, $T = 400$) | 1.712 | 1.935 | 2.296 |
| Simulated ($n = 400$, 800 replications) | 1.752 | 1.944 | 2.324 |
| Absolute difference | 0.040 | 0.009 | 0.028 |

The core regression is the one `radf()` already computes, applied to transformed data, so the check concentrates on the transform `kernel_purge()`.

Replication script: [replication/volatility-robustness/radf_kp_validation.R](#script-radf_kp_validation).

---

## SBZ (WLS + kernel volatility)

Status: done. The test is in `exuber/R/radf_sbz.R` (`radf_sbz_cv()`, `radf_sbz_union()`) and is checked in `exuber/tests/testthat/test-sbz.R`.

### Source

Harvey, D. I., Leybourne, S. J. & Zu, Y. (2019). Testing explosive bubbles with time-varying volatility. *Econometric Reviews*, 38(10), 1131–1151. We used the open working paper, Granger Centre Discussion Paper 18/05, University of Nottingham.

### Idea

SBZ is a weighted-least-squares variant of the PWY/PSY sup-ADF test. The weights come from a nonparametric kernel estimate of the time-varying volatility (eq. 6: Gaussian kernel, leave-one-out cross-validated bandwidth with $h \in [1/(2T), 1/6]$, footnote 2). The null distribution still depends on the volatility path, so size is controlled with a wild bootstrap, run jointly for $\sup \mathrm{DF}$ (the classic PWY/PSY test) and $\sup \mathrm{BZ}$ (the WLS test).

The paper's second contribution is a union of rejections that combines the two tests. We reject if

$$
U = \max\!\left(\sup \mathrm{DF},\ \frac{q_{\mathrm{DF}}}{q_{\mathrm{BZ}}} \sup \mathrm{BZ}\right) > q_{\mathrm{DF}},
$$

where $q_{\mathrm{DF}}$ and $q_{\mathrm{BZ}}$ are bootstrap critical values at the chosen level. The scaling constant makes the union asymptotically correctly sized (Theorem 3).

### Published numbers

Table 1 reports bootstrap $p$-values ($M = 499$) for FTSE (December 1985 to December 1999) and S&P 500 (January 1980 to March 2000).

| Series | $\sup\mathrm{DF}$ | $\sup\mathrm{BZ}$ | $U$ |
|---|---|---|---|
| FTSE Daily | 0.288 | 0.016 | 0.046 |
| FTSE Weekly | 0.275 | 0.146 | 0.201 |
| FTSE Monthly | 0.477 | 0.279 | 0.315 |
| SP500 Daily | 0.267 | 0.000 | 0.003 |
| SP500 Weekly | 0.170 | 0.002 | 0.011 |
| SP500 Monthly | 0.153 | 0.044 | 0.071 |

These agree with the paper's text: $\sup\mathrm{DF}$ rejects for no series, $\sup\mathrm{BZ}$ rejects at 5% for daily FTSE and all S&P 500 frequencies, and $U$ keeps those rejections at a slightly weaker level for monthly S&P 500.

The paper publishes no table of fixed critical values, because every $p$-value comes from a wild bootstrap on the authors' own series. Table 1 therefore cannot be reproduced bit for bit without their data and random draws. The check we can run is the empirical size under the null.

### Implementation

Three pieces are needed.

1. A Gaussian-kernel volatility estimator with leave-one-out bandwidth selection.
2. A WLS recursive Dickey–Fuller statistic with an intercept and heteroskedasticity weights, computed per window. It differs from `radf()` (OLS) and from `radf_tt()` (no intercept).
3. The wild bootstrap of `radf_wb.R`, extended to run $\sup\mathrm{DF}$ and $\sup\mathrm{BZ}$ on the same bootstrap draws, as the union requires, and to compute the scaling constant.

### Validation

Size under the null: random walks with $n = 150$, 150 replications, 199 bootstrap draws, nominal 5%.

| Statistic | Rejection rate |
|---|---|
| $\sup\mathrm{DF}$ | 0.033 |
| $\sup\mathrm{BZ}$ | 0.060 |
| $U$ | 0.053 |

All three rates are within Monte Carlo noise of the nominal level.

Replication script: [replication/volatility-robustness/radf_sbz_validation.R](#script-radf_sbz_validation).

---

## Pedersen & Schütte sieve bootstrap

Status: done. The sieve bootstrap of exuber (`radf_sb_cv()`) covers the method, and the lag selection of Pedersen & Schütte is available through its `type` argument.

### Source

Pedersen, T. Q. & Montes Schütte, E. C. (2020). Testing for explosive bubbles in the presence of autocorrelated innovations. *Journal of Empirical Finance*, 58, 207–225. We used the open working paper, CREATES Research Paper 2017-9.

### What is implemented

exuber's sieve bootstrap (`R/radf_sb.R`) follows the construction that Pedersen & Schütte propose for autocorrelated innovations. It fits an AR($p$) sieve to the first differences, resamples the residuals and rebuilds bootstrap paths with `stats::filter(..., "rec")`.

Their contribution is to show that a fixed lag order distorts size under autocorrelated innovations, and that BIC lag selection with a maximum lag `kmax` repairs it (Section 4, for example Table 4.8). `radf_sb_cv(type = "aic" / "bic", max_lag = ...)` selects the lag for each series by AIC or BIC with the `lag_select()` routine shared with `radf_wb_ps_cv()`, and uses the maximum across the panel. The lag is common to the whole panel, in line with the single `lag` argument of `radf()`. `type = "fixed"` keeps the original behaviour.

### Validation

`type = "bic"` picks a nonzero lag on AR(2)-autocorrelated data. On pure random-walk data the modal selection is 0 across 8 independent draws. A single draw is not asserted on, because BIC can pick a nonzero lag by chance in any one finite sample.

Replication script: [replication/volatility-robustness/radf_sb_cv_aic_bic_validation.R](#script-radf_sb_cv_aic_bic_validation).

---

## Hafner skewness-corrected wild bootstrap

Status: done, as `radf_wb_cv(..., dist_skew = TRUE)` and `radf_wb_distr(..., dist_skew = TRUE)`.

### Source

Hafner, C. M. (2020). Testing for bubbles in cryptocurrencies with time-varying volatility. *Journal of Financial Econometrics*, 18(2), 233–249. We used the IRTG 1792 Discussion Paper 2018-005 (Humboldt University Berlin).

### Idea

Hafner changes the multiplier distribution of the wild bootstrap of Harvey et al. (2016), so that the bootstrap approximates the PWY null distribution better when returns have time-varying volatility and are right-skewed, as cryptocurrency returns are. The multiplier is

$$
w_t = \frac{u_t}{\sqrt 2} + \frac{v_t^2 - 1}{2}, \qquad u_t, v_t \overset{\text{iid}}{\sim} N(0, 1) \text{ independent},
$$

which has $E[w_t] = 0$, $E[w_t^2] = 1$ and $E[w_t^3] = 1$. It is a fixed right-skewed multiplier, not matched to the skewness of each series, and it replaces the $N(0,1)$ or Rademacher multiplier in the otherwise unchanged wild bootstrap $y^*_t = w_t \hat e_t$.

### Implementation

A `dist_skew` argument passes through the wild-bootstrap code in `exuber/R/radf_wb.R`. `radf_wb_dgp_hlst(y, dist_rad, dist_skew = FALSE)` has a third multiplier branch, next to the default $N(0,1)$ branch and the Rademacher branch. `radf_wb_hlst()`, `radf_wb_cv()` and `radf_wb_distr()` pass the argument along (default `FALSE`) and reject `dist_rad = TRUE` together with `dist_skew = TRUE`.

### Validation

1. **Moments.** A simulation of $w_t$ with more than 500,000 draws gives $E[w] \approx 0$, $E[w^2] \approx 1$ and $E[w^3] \approx 1$.
2. **Default path.** With `dist_skew = FALSE` the bootstrap DGP is identical, draw for draw, to the one without the option.
3. **Power.** Under a mildly explosive alternative with ordinary innovations, `dist_skew = TRUE` rejects 86.7% of the time (30 replications).
4. **Size.** With the paper's own innovation distribution, negative log-$\chi^2(1)$ (footnote 1), the rejection rate at a nominal 5% is 3.3% without added heteroskedasticity (60 replications) and 0% with a deterministic volatility pattern added (80 replications). The test is conservative in small samples, which agrees with the paper's finding that it is undersized and that the effect grows with global heteroskedasticity.
5. **Heavy tails.** Combining log-$\chi^2$ innovations with a mild explosive alternative ($\rho = 1.06$) gives no power in short samples. The distribution of $-\log Z^2$ produces occasional extreme single-step outliers, and in a short series one of them can dominate the sum of squared residuals in both the observed and the bootstrap statistics. The same alternative with normal innovations has strong power (item 3), so this is a property of the noise process and not of `dist_skew`.

Replication scripts: [replication/volatility-robustness/hafner_dist_skew_moments_and_regression.R](#script-hafner_dist_skew_moments_and_regression), [hafner_dist_skew_power_and_size.R](#script-hafner_dist_skew_power_and_size).

---

## Sign-based sGSADF

Status: done, as `radf_sign()` / `radf_sign_cv()` (Harvey, Leybourne & Zu 2020) and `radf_sign_dm()` / `radf_sign_dm_cv()` (the recursively demeaned variant of Harvey, Leybourne, Tatlow & Zu 2025).

### Source

Harvey, D. I., Leybourne, S. J. & Zu, Y. (2020). Sign-based unit root tests for explosive financial bubbles in the presence of deterministically time-varying volatility. *Econometric Theory*, 36(1), 122–169, `doi:10.1017/S0266466619000057`.

Harvey, D. I., Leybourne, S. J., Tatlow, D. & Zu, Y. (2025). Unit root tests for explosive financial bubbles in the presence of deterministic level shifts. *Oxford Bulletin of Economics and Statistics*, 87(5), 879–901, `doi:10.1111/obes.12668`.

### Idea

The sign of $\Delta y_t$ depends only on which side of zero the innovation falls, so it does not depend on the innovation volatility. A test built on cumulated signs is therefore exactly invariant to any volatility pattern and needs no bootstrap. Let

$$
C_t = \sum_{i \le t} \operatorname{sign}(\Delta y_i).
$$

The statistic $s\mathrm{PSY}$, and its single-supremum case $s\mathrm{PWY}$ (eq. 4), is the double-supremum recursive Dickey–Fuller construction of PSY applied to $C_t$ instead of $y_t$, fitted without an intercept: $C_t = \rho(r_1, r_2)\, C_{t-1} + e_t$. Theorem 2 shows that its null limit does not depend on the volatility process $\sigma(s)$. Remark 1 adds that this follows from excluding the intercept and using only signs. A no-intercept version of PSY on the raw series would be invariant only asymptotically, whereas $s\mathrm{PSY}$ is invariant in finite samples.

### Implementation

The code is in `exuber/R/radf_sign.R`, modelled on `radf_tt.R`.

- `sign_transform(y)` returns `c(0, cumsum(sign(diff(y))))`.
- `radf_sign(data, minw)` transforms each series and calls `gls_dfstat_grid()`, the no-intercept recursive Dickey–Fuller routine of [STADF](#time-transformed-test-stadf--gstadf). It returns a `radf_sign_obj` that inherits from `radf_obj`, so the existing S3 methods apply.
- `radf_sign_cv(n, minw, nrep, seed)` simulates critical values from a random walk. The null distribution is pivotal (Theorem 2), so this is done once and not per dataset.

The `unroot()` and `rls_gsadf()` route of `radf()` cannot be reused. Despite the column names of `unroot(y, lag = 0)`, it fits a regression with an intercept, which is the wrong form here.

Not implemented: the paper's union of rejections with the standard PSY/PWY test through a wild bootstrap (Section 4), which has the same structure as the union in [SBZ](#sbz-wls--kernel-volatility), and the sign-based dating extension of Section 6.

### Validation

1. **Invariance.** A series and the same series multiplied by a strongly time-varying volatility pattern ($0.1$, $10$ and $1$ over three segments) give identical `sadf` and `gsadf` statistics, bit for bit.
2. **Published critical values.** At $T = 200$ and $\mathrm{minw}/T = 0.1$ the simulated critical values for $s\mathrm{PSY}$ (`gsadf_cv`) are $(3.482, 3.900, 4.925)$, against $(3.469, 3.901, 4.957)$ in Table 1 for 10%, 5% and 1%. For $s\mathrm{PWY}$ (`sadf_cv`) they are $(2.337, 2.665, 3.403)$, against $(2.405, 2.735, 3.434)$. At $T = 100$ the match is looser, and the published 1% value for $s\mathrm{PSY}$ (13.056) reflects the extreme-value behaviour that the paper notes at small samples and high quantiles.
3. **Finite-$T$ comparison.** The comparison has to use the matching finite $T$. The finite-sample critical values of $s\mathrm{PSY}$ converge to their asymptotic limit slowly (the paper's text says convergence is "fairly slow (particularly for $s\mathrm{PSY}$)"), and a $T = 300$ simulation lies between the paper's $T = 200$ and $T = 400$ rows.
4. **Power.** The rejection rate is 96.7% (30 replications) under a mildly explosive alternative, with the simulated 95% critical value of `radf_sign_cv()`.

Replication scripts: [replication/volatility-robustness/sign_based_invariance_and_power.R](#script-sign_based_invariance_and_power), [sign_based_finite_T_crosscheck.R](#script-sign_based_finite_T_crosscheck).

### Level-shift robustness (HLTZ 2025)

Status: done. `radf_sign()` and `radf_sign_cv()` carry a documentation section on level shifts, and `radf_sign_dm()` and `radf_sign_dm_cv()` implement the demeaned variant.

The 2025 paper evaluates the same sign-based statistics under a different departure from the null, deterministic level shifts in the series, with $n_T = O(T^{\alpha_n})$ shifts. Theorem 1 gives the null limit of the standard PSY statistic. Theorem 2 covers the cumulated-sign $s\mathrm{PWY}$/$s\mathrm{PSY}$, which is `radf_sign()`. Theorem 3 covers the recursively demeaned analogue $\bar s\mathrm{PWY}$/$\bar s\mathrm{PSY}$.

PSY needs a joint restriction on the number and the size of the shifts (Assumption 3). In the paper's Table 1, PSY is essentially never correctly sized once the number of shifts grows at rate $\sqrt T$ (Case 1), with empirical size up to 0.425 at a nominal 0.05. Both sign-based statistics need only a restriction on the number of shifts (Assumption 4, $\alpha_n < 1/2$) and none on their magnitude. At the boundary rate $\alpha_n = 1/2$ they pick up a term that depends on the shifts, but the resulting over-sizing is bounded independently of the shift magnitude (Remark 4).

The demeaned series of Theorem 3 uses an expanding-window mean:

$$
\tilde C_t = \sum_{i=2}^{t} \left\{ \operatorname{sign}(\Delta y_i) - \frac{1}{i-1} \sum_{j=2}^{i} \operatorname{sign}(\Delta y_j) \right\}.
$$

The inner sum is the running sum $C_i$ that `radf_sign()` already computes, so the demeaning is a one-time $O(T)$ transform, `sign_demean_transform()`, which feeds the same `gls_dfstat_grid()`.

Not implemented: the serial-correlation correction of Remark 6, which augments regression (4) with lagged $\Delta C_t$. The no-shift implementation does not handle serial correlation either.

`radf_sign_cv()`, `radf_sign_dm_cv()` and `radf_tt_cv()` carry the class tag `"mc_cv"`, so their `print()` and `summary()` methods fall through to the existing `tidy_radf_cv.mc_cv`.

#### Validation

1. **Formula.** `sign_demean_transform()` matches a brute-force loop that recomputes the recursive mean for each $i$, to about $1.8 \times 10^{-15}$.
2. **Table 1 of the paper.** For Case 1 ($\alpha_n = 0.5$) with $k = 2$, $\mu = 5$, $p = 0.8$ and $T = 400$, the published empirical size at a nominal 5% is 0.337 for PSY, 0.121 for $s\mathrm{PSY}$ and 0.050 for $\bar s\mathrm{PSY}$. Our replication with 500 replications, critical values simulated under the no-shift null, the paper's trimming $\mathrm{minw} = \lfloor 0.1T \rfloor$ and the same construction of shift count, size and location gives 0.302, 0.090 and 0.036. The ordering and magnitudes agree. PSY is badly oversized, $s\mathrm{PSY}$ mildly so, and $\bar s\mathrm{PSY}$ is closest to nominal and slightly conservative. The differences are within Monte Carlo noise at 500 replications against the paper's 2000.
3. **Invariance.** `radf_sign_dm()` keeps the exact invariance to volatility rescaling.
4. **Power.** `radf_sign_dm()` rejects a mildly explosive alternative with its own simulated critical value.

Replication script: [replication/volatility-robustness/radf_sign_dm_levelshift_validation.R](#script-radf_sign_dm_levelshift_validation).

### `datestamp()` and `autoplot()` support (done)

`radf_sign_cv()` and `radf_sign_dm_cv()` return the time-varying boundaries `badf_cv` and `bsadf_cv`, built in the same way as for [STADF/GSTADF](#time-transformed-test-stadf--gstadf): `gls_dfstat_grid()` on the transformed series already returns the supremum over window starts at each end point.

Checks:

1. The last row of `badf_cv` equals `adf_cv` for both functions.
2. Under a pure random walk ($n = 100$, `minw = 20`, 2000 critical-value replications, 200 test replications) the false-alarm rate is 5.5% for `radf_sign` and 3.5% for `radf_sign_dm`, at a nominal 5%.
3. On the synthetic bubble used for STADF/GSTADF (the `radf()` baseline detects 16%), `radf_sign` detects 20% and `radf_sign_dm` detects 8%. The lower rate for the demeaned variant reflects the trade-off between invariance and power in this family.
4. `datestamp()` (`option = "gsadf"` and `"sadf"`) and `autoplot()` run on `radf_sign(sim_data)` and `radf_sign_dm(sim_data)`.

---

## Stochastic explosive-coefficient test

Status: done. SSU is `ssu_test()`. GSSU, the UR/GUR union and the four CUSUM-type statistics are `ssu_test(type = "gssu", union = TRUE)` and `cusum_test()`.

### Source

Kurozumi, E. & Nishi, M. (2025). Testing for a bubble with a stochastically varying explosive coefficient. *Journal of Time Series Analysis*, 46(5), 945–965. Open access.

### Idea

The paper treats the explosive AR(1) coefficient itself as random:

$$
\rho_t = 1 + \frac{c_1}{T} + \frac{a\, u_t}{\sqrt T}, \qquad (2)
$$

where $u_t$ is i.i.d. with mean zero and unit variance, independent of the innovations. Every other method in this project assumes the deterministic form $1 + c/T^\alpha$. The motivation (Figure 1, rolling AR(1) estimates on Japanese daily stock prices) is that the speed of explosion looks unstable within a bubble episode. The innovation variance is constant (Assumption 1a), so this is not a volatility fix. It sits in this file because it appeared in the same *JTSA* special issue as the other heteroskedasticity papers.

The paper proposes three kinds of statistic.

1. **SSU/GSSU** (eq. 7–8) is a stochastic-unit-root test in the style of Lee (1998) and Nagakura (2009). It regresses squared differences on squared lagged levels, $(\Delta y_t)^2 = \mu^2 + \omega\, y_{t-1}^2 + \eta_t$. The raw $t$-statistic is not pivotal, and a bias-correction term $\hat\rho(r_1, r_2)$, a recursively estimated cross-moment between the residuals of the two regressions, gives the corrected statistic $t^{c}_{r_1, r_2}$ (following Nishi & Kurozumi 2024).
2. **CUSUM and CUSUM-SQ**, the parameter-constancy tests of Brown et al. (1975), applied to the same detection problem.
3. **A union of rejections** of SADF/GSADF with SSU/GSSU, which the paper recommends in practice because neither family dominates. SSU/GSSU is better when the coefficient is stochastic ($a \ne 0$) and SADF/GSADF when it is deterministic ($a = 0$).

### SSU

`ssu_test(data, minw = NULL, level = 0.95)` is in `exuber/R/ssu_test.R`. The regression of eq. 7 is a two-variable OLS over a window, so the closed-form window sums of `hls_prefix_sums()` and `hls_segment_coef()` apply, with $x = y_{t-1}^2$ and $z = (\Delta y_t)^2$. Expanding the residual cross-moment $\sum \hat\varepsilon_t \hat\eta_t$ gives a bilinear combination of window sums of twelve per-observation products. `ssu_prefix_sums()` builds the twelve cumulative sums and `ssu_stat_path()` evaluates the bias-corrected statistic $t^{\omega,c}_{0, r_2}$ at every end point in $O(1)$ per window. The statistic is the running maximum. Page 6 of the paper divides all three moments ($\sigma^2_\varepsilon$, $\sigma^2_\eta$, $\sigma_{\varepsilon\eta}$) by the window count minus 2, and the code does the same.

Table I of the paper gives the asymptotic critical values $2.90$, $3.30$ and $4.20$ at the 10%, 5% and 1% levels. They are scalars, because SSU fixes $r_1 = 0$. The default minimum window is the $r_0 = 0.01 + 1.8/\sqrt T$ that the paper recommends, which is the existing `psy_minw()`.

### GSSU, union and CUSUM-type statistics

**GSSU** (`ssu_test(type = "gssu")`). `ssu_stat_path()` takes a window start as well as an end, so the same twelve prefix sums give any window $(\mathrm{lo}, \mathrm{hi}]$ in $O(1)$. The statistic is the maximum over end points of the supremum over starts. The minimum window is the paper's $r_0 = -0.004 + 2.24/\sqrt T$, because the Table I note says that the `psy_minw()` formula oversizes GSSU. The critical values are $4.83$, $5.37$ and $6.81$.

**Union** (`ssu_test(union = TRUE)`). The UR statistic is

$$
\mathrm{UR} = \max\!\left( \frac{\mathrm{SADF}}{cv_{\mathrm{SADF}}},\ \frac{\mathrm{SSU}}{cv_{\mathrm{SSU}}} \right)
$$

and is compared with $1.16$, $1.13$ and $1.09$. GUR uses GSADF and GSSU and is compared with $1.11$, $1.10$ and $1.08$. A scaling constant is valid only at the level it was built for. The SADF/GSADF side is `radf(x, lag = 0)`, compared with the precomputed critical values at the sample size of the data.

**CUSUM-type statistics** (`cusum_test(type = "cs" | "gcs" | "cssq" | "gcssq")`, page 7). With $\sigma^2 = \operatorname{mean}((\Delta y)^2)$, not demeaned,

$$
S_k = \frac{\sum_{t \le k} \Delta y_t}{\sigma \sqrt T}.
$$

CS is $\max_k S_k$ and GCS is $\max_{j<k}(S_k - S_j)$, a running-minimum drawup. For CSSQ, with $\sigma_\eta^2 = \operatorname{mean}((\Delta y)^4) - \sigma^4$,

$$
D_k = \frac{\sum_{t \le k} (\Delta y_t)^2 - \tfrac{k}{T} \sum_{t \le T} (\Delta y_t)^2}{\sigma_\eta \sqrt T}.
$$

CSSQ is $\max_k D_k$ and $\min_k D_k$, and GCSSQ is the drawup and drawdown of $D$. The CUSUM-SQ tests are two-sided with each tail at $\alpha/2$. The CSSQ columns of Table I sit at the Brownian-bridge supremum quantiles for $\alpha/2$ (for example $1.32$ at 5%, against $\sqrt{-\log(0.025)/2} = 1.36$ before discretization), so a "level $\alpha$" row is already the two-sided test. All of these statistics are sups or infs of a partial-sum process and take $O(T)$ through a running max or min. Table I also gives the union constants, so no joint simulation is needed.

### Validation

- **Exact checks.** `ssu_stat_path()` matches a brute-force computation (two separate `lm()` fits and a manual cross-moment) to $2.4 \times 10^{-15}$ for windows with $\mathrm{lo} > 0$, and to $6.1 \times 10^{-15}$ for the GSSU path. All four CUSUM-type statistics, sup and inf, match a brute-force double loop to $8.9 \times 10^{-16}$. `test-ssu.R` checks every Table I column value by value. The default `minw` of `ssu_test()` equals `psy_minw()` exactly.
- **Size** at 5%, $n = 200$, 300 replications: SSU 0.067, GSSU 0.060, UR 0.050, GUR 0.057, CS 0.043, GCS 0.047, CSSQ 0.023, GCSSQ 0.037. A separate SSU run gives $12.0\%$, $9.0\%$ and $2.7\%$ at the 10%, 5% and 1% levels, so SSU is mildly oversized.
- **Power** at 5%, $n = 200$, 100 replications, with the bubble over the second half. With a stochastic coefficient ($c_1 = 3$, $a = 4$): SSU 0.90, GSSU 0.97, UR 0.87, GUR 0.96, CS 0.01, GCS 0.01, CSSQ 0.85, GCSSQ 0.78. This reproduces the paper's main finding (Theorem 2, Figure 2): CUSUM-type tests lose essentially all power once $a \ne 0$, SSU and CUSUM-SQ keep it, and the union stays close to the better of its two parts. With a deterministic coefficient ($c_1 = 10$, $a = 0$), every SSU, union and CUSUM-SQ statistic reaches 1.00, and CS and GCS reach 0.53.
- **Alternatives.** On a stochastic-coefficient DGP (the alternative of eq. 2, 60 replications) SSU has 85.0% power. On a deterministic explosive DGP it has 80.0%, against 90.0% for SADF, which is the trade-off of Theorem 2.

Tests are in `test-ssu.R` and `test-cusum-test.R`. The functions are also in pyexuber (`ssu_test(type=, union=)`, `cusum_test()`) and agree with the R values on a shared input.

Replication script: [replication/volatility-robustness/radf_ssu_validation.R](#script-radf_ssu_validation).

---

## SV-ADF

Status: done, as `datestamp(option = "svadf")`. The source is a preprint and has not been peer-reviewed, which is a lower bar than for the other sources in this project.

### Source

Sarkar, A. & Wells, M. T. (2026). Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance. arXiv:2604.12062. The theory is in the same authors' Double Local-to-Unity: Inference under Nearly Nonstationary Volatility, arXiv:2512.06823.

### Idea

SV-ADF extends the recursive right-tailed ADF test to highly persistent stochastic volatility, for example

$$
\log \sigma_t^2 = \phi_n \log \sigma_{t-1}^2 + \eta_t, \qquad \phi_n \to 1 \text{ at an iterated-logarithmic rate.}
$$

The other methods in this file need $\sigma(s)$ to be a fixed, non-stochastic function of calendar time. SV-ADF allows the volatility itself to be near-unit-root persistent, as in a GARCH with $\alpha + \beta$ close to 1.

The feasible statistics $\mathrm{SV\text{-}ADF}_r$ and $\mathrm{SV\text{-}ADF}_{rt}$ (eq. 3–4) are built from the same recursive OLS estimator and within-window residual variance $\hat\tau^2 = \tau^{-1} \sum \hat e^2$ as the recursive ADF statistic of `radf()` (eq. A.13–A.14). The contribution of Theorem 3.1 is an asymptotic justification under much weaker volatility conditions, with the limit coinciding with the homoskedastic one (Phillips & Yu 2009) once normalized.

What differs in practice is the threshold. Section 5.1 gives the calibration. For origination, the authors simulate the statistic under $H_0$ at $n \in \{100, 200, \dots, 1000\}$ (1000 replications each) and find that the 90th percentiles "are well approximated by $\log(n)/10$", which they adopt as the origination threshold. For collapse, they average the 10th-percentile threshold over random nuisance-parameter configurations and find it "most closely approximated by $\log(n)/2$". Both are closed-form in the sample size and need no estimated nuisance parameter. Origination and collapse use different thresholds (Remark 1).

### Implementation

`datestamp(data, option = "svadf", min_duration = NULL)` (helper `datestamp_svadf()` in `exuber/R/svadf.R`) reuses the `badf` sequence of `radf()`. `min_duration` defaults to `psy_ds(n)`, the existing $\log(n)$ rule of exuber, in place of the paper's data-frequency requirement of two consecutive months or one month. Origination is dated at the first run of at least `min_duration` consecutive points with `badf` above $\log(t)/10$. Collapse is dated at the first run of at least `min_duration` consecutive points with `badf` below $\log(t)/2$, searched only after the origination date.

### Validation

- The `badf` field of `datestamp(option = "svadf")` equals a direct `radf()` call, bit for bit, and the thresholds equal $\log(t)/10$ and $\log(t)/2$ exactly.
- Collapse is never dated before origination, which holds by construction and was confirmed over 20 replications.
- On a synthetic bubble and collapse (large-base bubble DGP) the detection rate is 100% over 20 replications. The mean absolute origination-date error is 4.95 periods and the mean absolute collapse-date error is 20.25 periods. Collapse is dated less precisely because the expanding `badf` window dilutes a post-collapse downward signal more than the upward signal at origination, a general property of single-recursion statistics.
- Under $H_0$ (60 replications, pure random walk) the false-alarm rate is 13.3%. This counts any origination crossing in a 150-period path, which is a much larger compound opportunity than a single-point 10% test.

`test-datestamp-svadf.R` has 7 tests. Replication script: [replication/volatility-robustness/datestamp_svadf_validation.R](#script-datestamp_svadf_validation).
