---
title: "Multivariate bubble tests"
blurb: "Panel and cross-series tests -- common bubbles, co-bubbles, and bubble contagion."
order: 4
---
This file covers tests for many series at once. They ask whether a bubble is shared across series, whether series explode together or one after another, and whether an explosive episode in one market passes into another. Three methods have their own code: common-bubble detection (Chen, Phillips & Shi), the co-bubble test (Evripidou, Harvey, Leybourne & Sollis) and the contagion regression (Greenaway-McGrevy & Phillips). The fourth, bubble migration, is an analysis that uses `radf()` and `datestamp()` series by series and needs no new code.

- `radf_common()` is in `exuber/R/radf_common.R` and is tested in `exuber/tests/testthat/test-common.R`. Its critical value depends on the panel width $N$ and is provided by `radf_common_cv()`.
- `cobubble_test()` is in `exuber/R/cobubble_test.R` and is tested in `exuber/tests/testthat/test-cobubble.R`.
- `contagion_reg()` is in `exuber/R/contagion_reg.R` and is tested in `exuber/tests/testthat/test-contagion.R`.

`radf()` returns `bsadf_panel` and `gsadf_panel` (`exuber/R/radf_.R`), but these are a panel average of independently estimated univariate BSADF sequences (`apply(bsadf, 1, mean)` followed by `max()`). They detect bubbles somewhere in the panel on average. They do not detect a shared latent bubble, lead and lag co-movement, migration or contagion, which is what the methods below add.

| Method | Paper | Status |
|---|---|---|
| [Common-bubble detection (PCA + PSY)](#common-bubble-detection-via-pca--psy) | Chen, Phillips & Shi (2020/2023) | done (including `radf_common_cv()`) |
| [Bubble migration](#bubble-migration) | Phillips & Yu (2011) | evaluated, no new code needed |
| [Co-bubble test](#co-bubble-test) | Evripidou, Harvey, Leybourne & Sollis (2022) | done |
| [Contagion regression](#contagion-regression) | Greenaway-McGrevy & Phillips (2015/2016) | done (the automatic delay search of eq. 8 is not included) |

All papers are listed in [references.md](/replication/references#multivariate).

## Common-bubble detection via PCA + PSY

### Source

Chen, Y., Phillips, P. C. B. & Shi, S. Common Bubble Detection in Large Dimensional Financial Systems. Cowles Foundation Discussion Paper 2251, Yale University, August 2020. Published in *Journal of Financial Econometrics*, 2023, 21(4), 989–1063. Open.

### Idea

The method is a two-step procedure for a bubble that is common to a panel of $N$ series. The motivating case in the paper is real-estate prices in 89 Chinese cities.

1. **PCA.** Estimate the dominant common factor of the panel by principal components, solving $\min_{\Lambda, F} \frac{1}{NT} \sum (X_{it} - \lambda_i f_t)^2$ subject to $\frac1N \Lambda'\Lambda = I_r$ (eq. 3.1–3.2). The estimated loadings are $\sqrt N$ times the eigenvectors of $X'X$ for the $r$ largest eigenvalues, and the factor estimate is $\tilde F = X \tilde\Lambda / N$. Only the first component $\tilde y_t$ is used for bubble detection ("sufficient... for the purpose of bubble identification", footnote 3), so $r$ need not be chosen by an information criterion.
2. **PSY on the factor.** Apply the PSY (2015a,b) recursive GSADF procedure to $\tilde y_t$, with the same ADF regression as `radf()`, $\tilde y_t = \mu + \rho\, \tilde y_{t-1} + v_t$ (eq. 3.3, OLS-demeaned, with an intercept).

This targets a different null and alternative from the panel average of exuber. The factor model (eq. 2.3, 2.7) mixes an I(1) factor (normal times), a mildly explosive factor (the bubble) and a stationary factor (after the collapse), with idiosyncratic errors for each series. The alternative is that a subsample of the panel shares the explosive factor.

### Theorems

- Theorem 4.2 (p. 12) says that the DF statistic computed on the estimated factor $\tilde y_t$ (eq. 4.2) has a "limit distribution... unaffected by factor estimation and is identical to that of the DF statistic computed from the original data, as in Phillips et al. (2015a)".
- Theorem 4.3 (p. 13, eq. 4.3) states that the limit of the resulting PSY-factor statistic is "then identical to that of the original PSY statistic (i.e., $F_{r_2}(W, r_0)$ in Phillips et al. (2015a))". $F_{r_2}(W, r_0)$ is the object that `radf_mc_cv()` simulates for `gsadf`.
- The identity is asymptotic ($N, T \to \infty$). Section 5 reports that in finite samples "the finite sample distribution lies to the left of the asymptotic, which implies slight undersizing if asymptotic critical values are employed" for small $N$ (Figure 1, $T = 60, 100, 140$, $N = 20$–$100$), under a DGP with a genuine shared explosive factor. The authors' own application therefore uses finite-sample simulated critical values.
- Theorems 4.5–4.8 give the alternative-hypothesis asymptotics, divergence rates and consistency of the origination and collapse dates. They are proved for the factor-model DGP, and the date-stamping rule of exuber (the $\log T$ duration filter) is not covered by them.
- The empirical application covers 89 Chinese cities, monthly from January 2003 to March 2013. The authors find three common-bubble episodes in Tier 1 and 2 cities and none in Tier 3. The underlying data are not available here.

### Implementation

`radf_common(data, r = 1, ...)` takes the first principal component of the panel matrix built by `parse_data()` (with `stats::prcomp()`) and calls `radf()` on the $T$-vector of scores. The result is a standard `radf_obj` on a single series, so `datestamp()` works on it unchanged.

### Critical values depend on the panel width

A panel with no common factor, $N$ independent random walks, is the sharpest null, because the first component has nothing to detect legitimately. Simulating `radf_common()` on such panels with $T = 250$ and 400 replications gives these null quantiles, to be compared with `radf_mc_cv()` at 95% (2.133):

| $N$ | 90% | 95% | 99% | difference from `radf_mc_cv()` at 95% |
|---|---|---|---|---|
| 6 | 2.420 | 2.662 | 3.057 | +0.53 |
| 20 | 3.239 | 3.484 | 3.960 | +1.35 |
| 50 | 4.014 | 4.333 | 4.922 | +2.20 |
| 100 | 4.890 | 5.201 | 5.714 | +3.07 |

The gap grows with $N$. At $N = 100$ the 95% critical value (5.2) is more than twice that of `radf_mc_cv()`. With more independent unit-root series, the leading component increasingly picks up whatever transient co-movement arises by chance, and its apparent persistence grows with $N$. The undersizing in Figure 1 of the paper is measured under a DGP that contains a genuine explosive factor, so it says nothing about this case. Using the univariate `radf_mc_cv(n, minw)`, which has no argument for $N$, would make the test badly oversized at realistic panel sizes.

`radf_common_cv(n, N, minw, nrep, seed)` simulates this null (an $N$-column panel of independent random walks, PCA, GSADF). It returns the same shape as `radf_mc_cv()` (`adf_cv`, `sadf_cv`, `gsadf_cv`, `badf_cv`, `bsadf_cv`), so it works as the `cv` argument of `datestamp()`, `tidy()` and `autoplot()`. The exuber join machinery matches the `method` attribute against `"Monte Carlo"` exactly, so the object carries that label and the panel width sits in a separate `N` attribute.

`test-common.R` checks that the function returns a `radf_cv`/`mc_cv`-shaped object that `datestamp()` accepts, and that the null quantile for $N = 30$ is significantly higher than for $N = 4$.

Replication script: [replication/multivariate/radf_common_validation.R](#script-radf_common_validation).

## Co-bubble test

### Source

Evripidou, A. C., Harvey, D. I., Leybourne, S. J. & Sollis, R. (2022). Testing for Co-explosive Behaviour in Financial Time Series. *Oxford Bulletin of Economics and Statistics*, 84(3), 624–650, `doi:10.1111/obes.12487`.

Status: done, as `cobubble_test()`.

### Idea

The test asks whether two series, each with an explosive episode, are related through co-explosive behaviour: a linear combination of the two is integrated of order zero although each series is locally explosive. It is the explosive analogue of cointegration.

The DGP (eq. 2) is $y_t = \mu_y + \beta_x\, x_{t-i} + \beta_z\, z_t + e_{y,t}$, where $x_t$ is observed and $z_t$ is an unobserved explosive process. Under $H_0\colon \beta_x > 0,\ \beta_z = 0$, the series $y_t$ and $x_{t-i}$ are co-explosive. The statistic (eq. 3) is a KPSS-type LM statistic on the OLS residuals of $y_t$ regressed on a constant and $x_{t-i}$:

$$
\hat e_{y,t} = y_t - \hat\alpha - \hat\beta\, x_{t-i}, \qquad
S = \hat\sigma_y^{-2}\,(T - |i|)^{-2} \sum_t \Bigl(\sum_{s \le t} \hat e_{y,s}\Bigr)^2 .
$$

The sum runs over the overlapping valid range $t = \max(i, 0) + 1, \dots, T + \min(i, 0)$. The testing direction is opposite to PSY/ADF tests, because the null is stationarity (co-explosivity) and not a unit root. Theorem 1 shows that the null limit does not depend on the properties of the regressor $x$, because its mild explosivity is asymptotically negligible for this statistic. It does depend on the pattern of heteroskedasticity in $e_{y,t}$, which rules out a fixed table of critical values. A wild bootstrap ($y^*_t = w_t \hat e_{y,t}$ with $w_t \sim \mathrm{IIDN}(0,1)$, refitted on the same regressor $x_{t-i}$ per Remark 2) reproduces the heteroskedasticity pattern and gives asymptotically size-controlled critical values (Theorem 2).

When the lead or lag $i$ is unknown (Section VI), it is estimated as $\hat i = \arg\min_j \hat\sigma_y^2(j)$ over a set of candidate lags. A misspecified $j \ne i$ leaves a neglected explosive term in the residuals that inflates their variance, so the variance-minimising $j$ consistently recovers $i$.

### Implementation

`cobubble_test(y, x, lag = NULL, lags = -6:6, nboot = 499L, level = 0.05, seed = NULL)` is in `exuber/R/cobubble_test.R`.

- `coexplosive_stat_aligned(y, xreg)` computes the statistic of eq. 3 on two aligned vectors of equal length.
- `coexplosive_stat(y, x, lag)` builds the pair $(y_t, x_{t-\mathrm{lag}})$ over the overlapping range and calls the aligned core.
- `coexplosive_select_lag(y, x, lags)` runs the $\hat i$ search of Section VI.
- `cobubble_test()` selects `lag` if none is given, computes the statistic and runs the wild bootstrap. Each bootstrap sample is regressed on the same $x_{t-\mathrm{lag}}$ regressor (Remark 2), because the paper found that omitting this makes the bootstrap a worse finite-sample match. It returns the statistic, the critical value, the $p$-value and the decision.

The recursive wild-bootstrap functions of `radf_wb.R` are not reused, because they are built around the recursive window structure of the ADF family. Here there is one static OLS fit per bootstrap draw.

### Validation

1. `coexplosive_stat()` matches a brute-force computation (a separate `lm()` call and an explicit loop for the cumulative sum) to floating-point precision (difference about $10^{-17}$).
2. Empirical size under $H_0$ with homoskedastic errors is 6.0% at a nominal 5% (100 replications).
3. Empirical size under $H_0$ with a volatility jump partway through the sample is 5.0% at a nominal 5%. This is the central claim of the paper (Theorem 2) and the reason for the wild bootstrap.
4. Power, where $y$ and $x$ have independent, unrelated explosive episodes, is a 100% rejection rate over 60 replications.
5. With a true lag of 3, `coexplosive_select_lag()` recovers it exactly in 20 of 20 seeds.

`test-cobubble.R` checks size and power against loose bounds. Replication script: [replication/multivariate/radf_cobubble_validation.R](#script-radf_cobubble_validation).

## Bubble migration

### Source

Phillips, P. C. B. & Yu, J. (2011). Dating the Timeline of Financial Bubbles During the Subprime Crisis. *Quantitative Economics*, 2(3), 455–491, `doi:10.3982/QE82`. Working paper: Cowles Foundation DP 1770. Open.

### Migration is an analysis, not a joint test

The paper has no separate migration statistic with its own null hypothesis, test equation or critical values. It does the following.

1. It applies the single-series recursive right-tailed unit-root test of Phillips, Wu & Yu (2011), the single-supremum predecessor of GSADF and the statistic `radf()` computes as `sadf` and `badf`, separately to seven unrelated financial series (the Nasdaq index, a home price index, asset-backed commercial paper, crude oil, platinum, the Baa bond rate and Pound/USD).
2. Each series gets its own origination and collapse date from the standard PWY rule ($\max \mathrm{DF}_r$ and $\max \mathrm{DF}_{r,t}$, with a $\log n$ minimum duration). Table 4 follows this format, for example "Heating oil: max DFr 6.9092, max DFrt 2.2416, origination March/08, collapse August/08".
3. "Migration" is a qualitative comparison of these independently estimated dates across series. The collapse of the equity-market bubble roughly precedes the origination of the housing-market bubble, which roughly precedes the mortgage-market bubble. The paper calls this a migration mechanism and matches it informally against the prediction of Caballero, Farhi & Gourinchas (2008).

From the Conclusion (pp. 34–35): "The dates are matched against the onset date for the subprime crisis as well as a specific sequential hypothesis concerning bubble migrations that are predicted in the theoretical model proposed by CFG (2008a)." The migration hypothesis is therefore compared with an external theoretical timeline and not tested statistically.

### Published numbers

Table 4 (search for additional series):

| Series | max DFr | max DFrt | origination | collapse |
|---|---|---|---|---|
| Heating oil | 6.9092 | 2.2416 | March/08 | August/08 |
| Coffee | -1.6035 | -0.7002 | NA | NA |
| Cotton | -0.2466 | -0.0866 | NA | NA |
| Cocoa | 2.4876 | 0.9872 | NA | NA |
| Sugar | -0.7408 | -0.2220 | NA | NA |
| Feeder cattle | 1.0336 | 0.4327 | NA | NA |
| Euro/USD | 0.4091 | 0.3311 | NA | NA |
| Yen/USD | 3.8949 | 1.4247 | NA | NA |
| Cnd/USD | 4.0494 | 2.6956 | Sep/21/07 | Nov/23/07 |

The abstract describes a bubble that "first emerged in the equity market during mid-1995 lasting to the end of 2000, followed by a bubble in the real estate market between January 2001 and July 2007 and in the mortgage market between November 2005 and August 2007", and that after the crisis erupted migrated "selectively into the commodity market and the foreign exchange market".

### In exuber

Reproducing the analysis means running `radf()` and `datestamp()` on each series and comparing the date ranges, which is what the paper does. A helper that lays several `datestamp()` outputs on a shared timeline would be a plotting convenience and not a statistical addition. See practitioner-guidance.md for the same idea on a more recent dataset.

## Contagion regression

### Source

Greenaway-McGrevy, R. & Phillips, P. C. B. Hot Property in New Zealand: Empirical Evidence of Housing Bubbles in the Metropolitan Centres. Cowles Foundation DP 2004, Yale University (2015). Published in *New Zealand Economic Papers*, 50(1), 88–113 (2016), `doi:10.1080/00779954.2015.1065903`. Open.

### Idea

Section 2.6 ("Bubble Contagion") defines a functional (time-varying) coefficient regression that tests and quantifies the contagion of bubble behaviour from a hypothesised core region to other regions. The paper applies it to New Zealand regional house prices, with Auckland City as the core. Equation numbers below follow the paper (pages 17 and 26).

1. **Rolling AR(1) coefficients.** For every region $i$ and every subsample-ending date $s$, estimate a fixed-width rolling-window OLS AR(1) regression (eq. 1, $y_t = \delta + \beta\, y_{t-1} + e_t$, a plain Dickey–Fuller regression with an intercept and one lag) to get slope estimates $\hat\beta_{i,s}$. The window width is $S = \lfloor 0.33\,T \rfloor$ in the paper's application. This is a moving-window version of the expanding-window `badf` construction, computed from the prefix-sum pattern of `hls_segment_ssr()` and `hls_prefix_sums()`.
2. **Functional regression** (eq. 4, p. 17):
   $$\hat\beta_{j,s} = \delta_{1j} + \delta_{2j}\,\frac{s}{T - S + 1}\,\hat\beta_{\mathrm{core}, s-d} + \text{error}_s,$$
   where $d \in \{0, \dots, 12\}$ is an integer delay. The text calls it months, but the empirical section discusses delays in quarters against quarterly data, so $d$ is read as native sampling periods of the input series.
3. **Time-varying coefficient** (eq. 6, p. 26). The local-constant Nadaraya–Watson estimate with a Gaussian kernel is
   $$\hat\delta_{2j}(r; h, d) = \frac{\sum_s K_{hs}(r)\, \tilde\beta_{j,s}\, \tilde\beta_{\mathrm{core}, s-d}}{\sum_s K_{hs}(r)\, \tilde\beta^2_{\mathrm{core}, s-d}},$$
   with $\tilde\beta_{j,s} = \hat\beta_{j,s} - \operatorname{mean}(\hat\beta_{j,\cdot})$ (centred) and $K_{hs}(r) = h^{-1} K\bigl((s/T - r)/h\bigr)$. It is the closed-form solution of a no-intercept, single-regressor weighted least squares fit at each $r$, a ratio of two weighted sums. It vectorises as one $T \times T$ Gaussian weight matrix times two length-$T$ vectors.
4. **Bandwidth** (eq. 7, p. 26), by leave-one-out cross-validation:
   $$\check h_{jT}(d) = \arg\min_{h \in H_T} \sum_s \bigl\{\tilde\beta_{j,s} - \check\delta_{2j}\bigl(s/(T-S+1); h, d\bigr)\, \tilde\beta_{\mathrm{core}, s-d}\bigr\}^2, \qquad H_T = \bigl[(T-S+1)^{-1/2},\ (T-S+1)^{-1/10}\bigr],$$
   where $\check\delta_{2j}$ is the leave-one-out version of eq. 6 (it excludes $p = s$). This is a bounded one-dimensional optimisation, so `stats::optimize()` applies.
5. **Delay** (eq. 8, p. 26). Select $d \in \{0, \dots, 12\}$ by minimising the cross-validated SSE. The prose says to choose $d$ by NLS with the largest $R^2$, and eq. 8 says to minimise the SSE. These are one criterion, because $R^2 = 1 - \mathrm{SSE}/\mathrm{TSS}$ and the TSS of the centred $\tilde\beta_{j,s}$ sequence is constant across $d$.

The paper reports no numeric table of estimated $d$, $h$ or $\hat\delta_{2j}$. Its results are Figures 7–8 (time-varying contagion coefficients from Auckland City), which are graphical. It performs no formal inference on $\hat\delta_{2j}(r)$, and consists of point estimation, cross-validated tuning and plots.

### Implementation

`contagion_reg(y, core, S, d, h, r_grid)` is in `exuber/R/contagion_reg.R`. It computes the fixed-window AR(1) sequence (eq. 1), the Nadaraya–Watson regression at a caller-supplied delay `d` (eq. 6) and the cross-validated bandwidth (eq. 7) when `h` is not supplied. The automatic delay search of eq. 8 is not included. It can be done by calling `contagion_reg` once per candidate $d$ and comparing the cross-validated SSE.

A window "$\{t = s-S+1, \dots, s\}$" has $S$ levels and so $S - 1$ regression pairs. The Nadaraya–Watson ratio in `contagion_nw_delta2()` uses `crossprod(K, v)` to sum over $s$, and the cross-validation helper `contagion_loocv_sse()` uses the same orientation, since the kernel weight matrix is not symmetric.

### Validation

No published numeric table exists to match.

- `contagion_fixed_window_beta()` (eq. 1) matches a brute-force `lm()` fit (`< 1e-8`) at five window-end dates.
- `contagion_nw_delta2()` (eq. 6) matches a manual Gaussian-kernel weighted least squares ratio (`< 1e-10`).
- `contagion_loocv_sse()` (eq. 7) matches a manual leave-one-out double loop (`< 1e-8`).
- `contagion_bandwidth_cv()` picks a bandwidth strictly inside the interval $H_T$, with cross-validated SSE no worse than at either endpoint.
- A synthetic series whose local persistence tracks that of the core series (with a known delay) shows a wider range of estimated $\delta_2(r)$ than an independent series (mean range 0.50 against 0.42 over 15 replications), so the estimator responds to a genuine time-varying relationship.

`test-contagion.R` has 19 tests. Replication script: [replication/multivariate/radf_contagion_validation.R](#script-radf_contagion_validation).
