---
title: "Multivariate bubble tests"
blurb: "Panel and cross-series tests -- common bubbles, co-bubbles, and bubble contagion."
order: 4
---
This file covers tests for many series at once. They ask whether a bubble is
shared across series, whether series explode together or one after another, and
whether an explosive episode in one market passes into another. Three of the
four methods needed new code: common-bubble detection (Chen, Phillips & Shi),
the co-bubble test (Evripidou, Harvey, Leybourne & Sollis) and the contagion
regression (Greenaway-McGrevy & Phillips). The fourth, bubble migration, needs
no new code.

- `radf_common()` is in `exuber/R/radf_common.R` and is tested in
  `exuber/tests/testthat/test-common.R`. The detection statistic was nearly free
  to implement. Independent validation (2026-08-09, below) showed that the
  critical value we had first recommended, `radf_mc_cv(n, minw)`, is not valid
  at realistic panel sizes. The correct, `N`-dependent critical value ships as
  `radf_common_cv()`. The "Independent validation" and "Update" subsections
  under Common-bubble detection give the numbers.
- `cobubble_test()` (2026-08-09) is in `exuber/R/cobubble_test.R` and is tested
  in `exuber/tests/testthat/test-cobubble.R`. See "Implementation" under
  Co-bubble test.
- `contagion_reg()` (2026-08-10, minimum-viable subset) is in
  `exuber/R/contagion_reg.R` and is tested in
  `exuber/tests/testthat/test-contagion.R`. See "Implementation" under
  Contagion regression. Re-triaging this item showed that our earlier
  assessment, "most expensive, no reuse", was wrong on two of its three cost
  drivers.

Before these additions exuber had no proper multivariate bubble test.
`radf()` returns `bsadf_panel` and `gsadf_panel` (`exuber/R/radf_.R`, lines
98-99), but these are a panel average of independently estimated univariate
BSADF sequences (`apply(bsadf, 1, mean)` followed by `max()`). They are not a
common-factor or cross-series test. They detect "bubbles somewhere in the
panel on average", and they do not detect a shared latent bubble, lead/lag
co-movement, migration or contagion. All four items below were new
capabilities relative to what existed.

| Method | Paper | Status |
|---|---|---|
| [Common-bubble detection (PCA + PSY)](#common-bubble-detection-via-pca--psy) | Chen, Phillips & Shi (2020/2023) | done (incl. `radf_common_cv()`) |
| [Bubble migration](#bubble-migration) | Phillips & Yu (2011) | evaluated, no new code needed |
| [Co-bubble test](#co-bubble-test) | Evripidou, Harvey, Leybourne & Sollis (2022) | done |
| [Contagion regression](#contagion-regression) | Greenaway-McGrevy & Phillips (2015/2016) | done (2026-08-10, minimum-viable subset, eq. 8's automatic delay search not included) |

All papers: [references.md](/replication/references#multivariate).

## Common-bubble detection via PCA + PSY

### Source

Chen, Y., Phillips, P.C.B. & Shi, S. "Common Bubble Detection in Large
Dimensional Financial Systems." Cowles Foundation Discussion Paper No. 2251,
Yale University, August 2020. Published version: *Journal of Financial
Econometrics*, 2023, 21(4), 989-1063.

We accessed it through the open Cowles Foundation PDF, which has no paywall.
We converted the PDF with `pdftotext -layout` and read the theorem statements
quoted below directly from that extraction. We did not render PNG pages,
because the equations in question are short and the `pdftotext -layout`
output was unambiguous on inspection. This is a lighter verification standard
than most other files in this project use, and we note it for honesty.

### What it is

The method is a two-step procedure for detecting a bubble that is common to a
panel of N series. The motivating case in the paper is real-estate prices in
89 Chinese cities.

1. Step 1 is PCA. We estimate the dominant common factor of the panel by
   principal components, solving `min_{Λ,F} (1/NT) Σ (X_it − λ_i f_t)²`
   subject to `(1/N)Λ'Λ = I_r` (their eq. 3.1-3.2). The estimated loadings
   are `√N` times the eigenvectors of `X'X` for the `r` largest eigenvalues,
   and the factor estimate is `F̃ = XΛ̃/N`. Only the first component is
   used for bubble detection (`ỹ_t`, "sufficient... for the purpose of
   bubble identification," footnote 3), so we do not need to determine `r`
   through an information criterion for this application.
2. Step 2 is PSY on the factor. We apply the standard PSY (2015a,b)
   recursive-evolving GSADF procedure to `ỹ_t`. It uses the same ADF
   regression form as the existing `radf()` of exuber:
   `ỹ_t = μ + ρ ỹ_{t-1} + v_t` (their eq. 3.3, OLS-demeaned, intercept
   included). It does not use the GLS/no-intercept form of `radf_tt()`.

This targets a different null and alternative than the panel average of
exuber. The factor model (their eq. 2.3, 2.7) explicitly mixes an I(1)
factor (normal times), a mildly explosive factor (the bubble) and a
stationary factor (after the collapse), with idiosyncratic errors for each
series. The alternative is that a subsample of the panel shares this
explosive factor.

### Exact numbers / theorem reproduced

We verified the following directly from the theorem statements and not from
a secondary citation.

- Theorem 4.2 (p.12 of the PDF) says that the DF statistic computed on the
  estimated factor `ỹ_t` (eq. 4.2, a ratio of stochastic integrals of
  standard Brownian motion `W₁`) has a "limit distribution... unaffected by
  factor estimation and is identical to that of the DF statistic computed
  from the original data, as in Phillips et al. (2015a)."
- Theorem 4.3 (p.13, eq. 4.3) states verbatim that the limit of the resulting
  "PSY-factor" statistic is *"then identical to that of the original PSY
  statistic (i.e., `F_{r2}(W, r0)` in Phillips et al. (2015a))"*. We confirmed
  the claim in the theorem of the paper itself and did not infer it.
  `F_{r2}(W, r0)` is exactly the object that `radf_mc_cv()` already simulates
  for `gsadf`.
- The identity is asymptotic (N,T → ∞), and a casual read does not make that
  obvious. Section 5 (Simulations) reports that in finite samples "the finite
  sample distribution lies to the left of the asymptotic, which implies slight
  undersizing if asymptotic critical values are employed" for small N (their
  Figure 1, T = 60/100/140, N = 20-100). The distribution converges rapidly as
  N and T grow, under their own DGP, which includes a genuine shared explosive
  factor. Their own empirical application therefore uses finite-sample
  simulated critical values and not the published PSY asymptotic table.
  Correction (2026-08-09, from the independent validation below): plain
  `radf_mc_cv(n, minw)` is not a safe substitute. It has no dependence on `N`
  at all, and the validation section below shows that the true null quantiles
  diverge further from those of `radf_mc_cv()` as `N` grows, and do not
  converge. The first draft of this section treated `radf_mc_cv(n, minw)` as
  "the correct match", which was wrong. The validation section below shows
  what needs simulating instead.
- We located Theorems 4.5-4.8 (alternative-hypothesis asymptotics, divergence
  rates, origination and collapse date consistency) but did not read them in
  depth, because they were outside this feasibility pass. A faithful
  implementation that needs the date-stamping consistency proof, beyond the
  detection statistic, would have to return to them.
- The empirical application covers 89 Chinese cities, January 2003 to March
  2013 monthly, with factor loadings drawn or estimated in [0.3, 1.7] and an
  idiosyncratic error SD of about 0.1. The authors detect three common-bubble
  episodes in Tier 1/2 cities and none in Tier 3. We report this for
  completeness. We cannot re-derive it without the underlying data, so we did
  not cross-check it bit for bit.

### Cost/feasibility for exuber, near-free, confirmed

This is the cheapest item in this file, and it holds up on inspection.

- The PCA step uses base R `stats::prcomp()` (or a two-line manual
  eigen-decomposition of `X'X`) on the panel matrix that `radf()` already
  builds via `parse_data()`. It needs no new dependency and no new C++.
- The GSADF step is `radf(factor_series)`, the existing function unmodified,
  called on a length-T vector (the first principal-component score). The code
  needs no changes for this, because `radf()` and `radf_mc_cv()` do not care
  where their input vector came from. The caveat above, and the independent
  validation below, still apply. The critical values of `radf_mc_cv(n, minw)`
  are valid for the output of `radf_common()` only at the specific `N` (panel
  width) they were calibrated at, and in practice that is never, because
  `radf_mc_cv()` has no `N` argument at all. Using them off the shelf, as this
  file first suggested, is not sound at realistic panel sizes.
- The new code is small. (a) A thin wrapper function, for example
  `radf_common(data, r = 1, ...)`, calls `prcomp()` on the panel, takes the
  scores of PC1 as a `T`-vector, calls `radf()`, and prints or returns the
  result tagged as a common-bubble result. (b) Optionally, date-stamping
  helpers can reuse `datestamp()` unchanged, because the output is a standard
  `radf_obj` on a single series once the factor is extracted.
- A first pass should leave one thing out. The consistency guarantees of
  Theorems 4.5-4.8 for origination and collapse dates are proved for the
  factor-model DGP specifically. `datestamp()` would run mechanically on the
  factor series regardless, but to claim that the consistency theorem of the
  paper applies to the exact date-stamping rule of exuber (the
  `log(T)`-duration filter and so on), we would need a closer read of Section
  4.3, which we have not done.
- It ranks first because it is a one-to-two file addition that reuses all of
  the existing panel-parsing (`parse_data()`), estimation (`radf()`) and
  critical-value (`radf_mc_cv()`/`radf_wb_cv()`) machinery. It needs no new
  statistic, no new asymptotic theory and no new bootstrap. The only new code
  is the PCA call itself.

### Independent validation (2026-08-09), Theorem 4.3's asymptotic identity does not hold at practical N, and the gap grows with N

The direct test of the claim of Theorem 4.3 simulates a panel with no true
common factor at all, namely `N` independent random walks. This is the
sharpest possible null, because PC1 has nothing to detect legitimately. We run
`radf_common()` many times and compare its null quantiles with those of plain
`radf_mc_cv()`, the univariate GSADF benchmark that the theorem says it should
coincide with asymptotically. We use `T=250` throughout, with seed=271828.

```r
for (N in c(6, 20, 50, 100)) {
  gsadf_null <- replicate(400, {
    panel <- replicate(N, cumsum(rnorm(250)))   # independent walks, no factor
    radf_common(panel)$gsadf
  })
  quantile(gsadf_null, c(0.9, 0.95, 0.99))
}
```

| N | 90% | 95% | 99% | diff from `radf_mc_cv()`'s 95% (2.133) |
|---|---|---|---|---|
| 6 | 2.420 | 2.662 | 3.057 | +0.53 |
| 20 | 3.239 | 3.484 | 3.960 | +1.35 |
| 50 | 4.014 | 4.333 | 4.922 | +2.20 |
| 100 | 4.890 | 5.201 | 5.714 | +3.07 |

The result is the opposite of what a "converges as N grows" reading of
Section 5 of the paper would predict, and the gap is large. At `N=100` the
true 95% critical value (5.2) is more than double that of `radf_mc_cv()`
(2.1). The mechanism differs from the one that Figure 1 of the paper
documents. Their finite-sample undersizing result is measured under their own
DGP, which contains a genuine shared explosive factor, so it says nothing
about behavior when the "common factor" premise is false. When we extract PC1
from a panel of purely independent (non-cointegrated) I(1) series, it does not
behave like a single random walk once `N` grows. With more independent
unit-root series to draw from, the leading principal component increasingly
picks up whatever transient co-movement exists by chance, and the persistence
and apparent explosiveness of that component grow with `N`. This is a
recognized general hazard of applying PCA or factor extraction to
non-cointegrated nonstationary panels. We have not verified it against a
specific citation, so treat it as a claim about the mechanism and not as a
sourced fact. We keep to the usual standard of this file and assert no
unverified citations.

The practical consequence is that using `radf_mc_cv(n, minw)`, the plain
univariate critical value without an `N` argument, as the critical value of
`radf_common()` would make the test badly oversized. It would detect "common
bubbles" that are not there for any realistically sized panel, and the problem
gets worse as the panel grows. This corrects the earlier claim in this section
that "the correct match on the exuber side is `radf_mc_cv(n, minw)`", made
when the section was first drafted, before this validation pass. That claim is
wrong at any `N` beyond a handful of series.

**Update (2026-08-09, Bundle 1): implemented.** `radf_common_cv(n, N, minw,
nrep, seed)` now ships. It simulates the null exactly as above (an
`N`-column panel of independent random walks, PCA, GSADF) and does not defer
to `radf_mc_cv()`. It returns the same shape as `radf_mc_cv()`
(`adf_cv`/`sadf_cv`/`gsadf_cv`/`badf_cv`/`bsadf_cv`), so it works as a drop-in
`cv` argument for `datestamp()`, `tidy()` and `autoplot()` on a
`radf_common()` result. One implementation wrinkle remains. The generic
`is_mc()` and join machinery of exuber (used by `datestamp()` and others) does
an exact string match of the `method` attribute against `"Monte Carlo"`.
`radf_common_cv()` therefore has to set `method = "Monte Carlo"`, and the
panel-specific information lives in a separate `N` attribute instead of in the
method label.

`test-common.R` tests two things. The function runs and returns a
`radf_cv`/`mc_cv`-shaped object that `datestamp()` can use directly. And, as a
direct re-confirmation of the independent-validation finding above, the null
quantile for `N=30` is significantly higher than for `N=4`, where a naive
reading of Theorem 4.3 might expect them to match. The test suite now
exercises the dependence on panel width that validation found, in addition to
the one-off simulation script.

Replication script:
[replication/multivariate/radf_common_validation.R](#script-radf_common_validation).

## Co-bubble test

### Source

Evripidou, A.C., Harvey, D.I., Leybourne, S.J. & Sollis, R. (2022). "Testing
for Co-explosive Behaviour in Financial Time Series." *Oxford Bulletin of
Economics and Statistics*, 84(3), 624-650. `doi:10.1111/obes.12487`.

Status: done (2026-08-09), shipped as `cobubble_test()`. We read the full PDF,
from the abstract through Section VI, including the model, Theorems 1-2, the
wild bootstrap algorithm and the lag-selection procedure. We did not stop at
the abstract.

### What it is

The test asks whether two series, each containing an explosive (bubble)
episode, are related through co-explosive behaviour. This means that a linear
combination of the two is integrated of order zero even while each series is
locally explosive. It is analogous to cointegration, but for explosive and
not unit-root regimes, hence "co-bubble".

The DGP (their eq. 2) is `y_t = mu_y + beta_x*x_{t-i} + beta_z*z_t +
e_y,t`, where `x_t` is observed and `z_t` is an unobserved explosive process.
Under `H0: beta_x > 0, beta_z = 0`, `y_t` and `x_{t-i}` are co-explosive. The
test statistic (eq. 3) is a KPSS-type LM statistic on the OLS residuals of
`y_t` regressed on a constant and `x_{t-i}`:

```
e_y,t = y_t - alpha_hat - beta_hat * x_{t-i}          (OLS residual)
S = sigma_y_hat^-2 * (T-|i|)^-2 * sum_t (cumsum of e_y up to t)^2
```

The sum runs over the overlapping valid range `t = max(i,0)+1, ..., T+min(i,0)`.
The testing direction is the opposite of PSY/ADF-style tests, because the null
is stationarity (co-explosivity) and not the presence of a unit root. The
paper proves (Theorem 1) that the limit null distribution does not depend on
the properties of the `x` regressor at all, since its mild explosivity is
asymptotically negligible for this statistic. It does depend on the pattern of
heteroskedasticity in `e_y,t`. That rules out a fixed table of critical values,
so a wild bootstrap (`y*_t = w_t * e_y,t`, `w_t ~ IIDN(0,1)`, refit on the same
`x_{t-i}` regressor per Remark 2) reproduces the heteroskedasticity pattern and
gives asymptotically size-controlled critical values (Theorem 2).

When the lead or lag `i` is unknown (Section VI), we estimate it as
`i_hat = argmin_j sigma_y_hat(j)^2` over a candidate set of lags. A
misspecified `j != i` leaves a neglected explosive term in the residuals that
inflates their variance, so the variance-minimizing `j` consistently recovers
`i`.

### Implementation

The function ships as `cobubble_test(y, x, lag = NULL, lags = -6:6, nboot =
499L, level = 0.05, seed = NULL)` in `exuber/R/cobubble_test.R`.

- `coexplosive_stat_aligned(y, xreg)` computes the core KPSS-type statistic
  (eq. 3) on two already aligned vectors of equal length.
- `coexplosive_stat(y, x, lag)` builds the `(y_t, x_{t-lag})` pair over the
  overlapping valid range of the paper and calls the aligned core.
- `coexplosive_select_lag(y, x, lags)` runs the `i_hat` search of Section VI.
- `cobubble_test()` ties these together. It selects `lag` if none is given,
  computes the observed statistic, and runs the wild bootstrap. Each bootstrap
  sample is regressed on the same `x_{t-lag}` regressor, per Remark 2, because
  the paper found that omitting it makes the bootstrap distribution a worse
  finite-sample match. The function returns the statistic, the critical value,
  the p-value and the reject/no-reject decision.

The function does not reuse the existing wild-bootstrap DGP functions of
`radf_wb.R` (`radf_wb_dgp_ps()` and `radf_wb_dgp_hlst()`). Those are built
around the recursive ADF-family ptr/window structure and do not apply to a
single-shot KPSS statistic on a fixed regression. The idea of resampling with
wild multipliers is the same, but the mechanics (one static OLS fit and not a
recursive scan) are simple enough that no shared scaffolding is needed.

Independent validation:

1. Formula-exact check: `coexplosive_stat()` matches an independently written
   brute-force computation (a separate `lm()` call, and a manual loop-based
   cumulative sum instead of the vectorized `cumsum()`) to floating-point
   precision (`diff ~ 1e-17`).
2. Empirical size under H0 with homoskedastic errors is 6.0% at the nominal 5%
   level (100 Monte Carlo reps), within sampling noise.
3. Empirical size under H0 with heteroskedastic errors (a volatility jump
   partway through the sample) is 5.0% at the nominal 5% level. This is the
   central claim of the paper (Theorem 2) and the reason for the wild
   bootstrap in place of a fixed KPSS table, and it holds up.
4. Power under H1, where `y` and `x` have independent, unrelated explosive
   episodes and are not co-explosive, is a 100% rejection rate over 60 reps.
5. Lag recovery: with a true lag of 3 in the DGP,
   `coexplosive_select_lag()` recovered the exact true lag in 20 of 20 seeds.

The new tests in `test-cobubble.R` check the size and power Monte Carlo
results against loose bounds. This matches the usual convention of this
project, which tests stochastic simulation results against a tolerance and not
against an exact number.

Replication script:
[replication/multivariate/radf_cobubble_validation.R](#script-radf_cobubble_validation).

### Cost/feasibility note (as originally scoped, before implementation)

Of the four items in this file, this was the one structurally most different
from anything already in exuber. It is a KPSS-type (LM/stationarity) statistic
and not an ADF/DF-family recursive one, and it is inherently bivariate and not
panel-shaped. It turned out to be tractable because the statistic is
closed-form for each fixed `(lag, i)`, needing one OLS fit and one cumulative
sum, with no nonparametric bandwidth selection and no recursive window scan.
The wild bootstrap loop is therefore just "refit and recompute the same
closed-form statistic `nboot` times". It has the same computational shape as
the loop in `radf_wb_cv()`, even though the underlying statistic and DGP
function are unrelated to it.

## Bubble migration

### Source

Phillips, P.C.B. & Yu, J. (2011). "Dating the Timeline of Financial Bubbles
During the Subprime Crisis." *Quantitative Economics*, 2(3), 455-491.
`doi:10.3982/QE82`. Working paper: Cowles Foundation Discussion Paper No.
1770.
fully open, with no paywall.

### What it is, important finding: "migration" is not a distinct joint test

Reading the primary source directly overturns the naive framing that the
title alone suggests. The paper contains no separate "migration test"
statistic with its own null hypothesis, test equation or critical values. What
the paper does is this.

1. It applies the same single-series recursive right-tailed unit-root
   methodology as Phillips, Wu & Yu (2010) and Phillips & Yu (2009), the
   PWY-style sup-ADF test. This is the predecessor of the GSADF of PSY, and
   exactly the statistic that `radf()` already computes (`sadf`/`badf`, the
   single-sup PWY branch and not even the double-sup GSADF). The paper applies
   it separately and independently to seven unrelated financial series (the
   Nasdaq index, a home price index, asset-backed commercial paper, crude oil,
   platinum, the Baa bond rate and Pound/USD).
2. Each series gets its own origination and collapse date, estimated by the
   standard PWY dating rule (`max DF_r`/`max DF_{r,t}`, a `log(n)` minimum
   duration). Table 4 of the paper reproduces exactly this format, for example
   "Heating oil: max DFr 6.9092, max DFrt 2.2416, origination March/08,
   collapse August/08".
3. "Migration" is a narrative, qualitative comparison of these independently
   estimated dates across series after the fact. The paper observes that the
   collapse date of the equity-market bubble roughly precedes the origination
   date of the housing-market bubble, which roughly precedes the
   mortgage-market bubble, and so on. It calls this sequence a "migration
   mechanism" and matches it informally against the theoretical prediction in
   Caballero, Farhi & Gourinchas (2008), without a formal joint statistic.

Direct quote (Conclusion, p.34-35 of the PDF): *"The dates are matched
against the onset date for the subprime crisis as well as a specific
sequential hypothesis concerning bubble migrations that are predicted in
the theoretical model proposed by CFG (2008a)."* The "test" of the migration
hypothesis is therefore a match against an external theoretical timeline and
not a statistical hypothesis test internal to the econometric methodology.

### Exact numbers reproduced

Table 4 (search for additional series, secondary dataset), reproduced
exactly from the PDF text:

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

The primary-series migration narrative (Abstract, verbatim) is that a bubble
"first emerged in the equity market during mid-1995 lasting to the end of
2000, followed by a bubble in the real estate market between January 2001 and
July 2007 and in the mortgage market between November 2005 and August 2007,"
and then, after the crisis erupted, migrated "selectively into the commodity
market and the foreign exchange market."

### Cost/feasibility note, effectively zero marginal statistical cost, but nothing new to build

The per-series test is exactly what `radf()` already computes (or even just
its `sadf`/`badf` output, ignoring the `gsadf`/`bsadf` fields that PSY added
later), and "migration" is a narrative overlay and not a formal statistic.

- There is no new test to implement. Reproducing the "migration" analysis of
  Phillips & Yu with exuber today means running `radf()` and `datestamp()` on
  each of several series independently (already fully supported) and then
  comparing the resulting date ranges by eye, which is exactly what the paper
  does.
- The only concrete deliverable this item could motivate is a vignette or a
  small helper that lays out several `datestamp()` outputs on a shared
  timeline, for visual comparison of the "migration". That is a plotting and
  reporting convenience and no statistical addition. See
  practitioner-guidance.md for the same idea
  applied to a more current dataset.
- It ranks below common-bubble detection for a different reason than the
  others. It is not expensive. There is simply nothing statistically new to
  build, so implementing it would be documentation or vignette work that
  demonstrates an existing capability. That is a decision about product scope
  and not a question of feasibility.

## Contagion regression

### Source

Greenaway-McGrevy, R. & Phillips, P.C.B. "Hot Property in New Zealand:
Empirical Evidence of Housing Bubbles in the Metropolitan Centres." Cowles
Foundation Discussion Paper No. 2004, Yale University (2015). Published
version: *New Zealand Economic Papers*, 50(1), 88-113 (2016).
`doi:10.1080/00779954.2015.1065903`.
fully open, with no paywall.

### What it is

Section 2.6 ("Bubble Contagion") defines a functional (time-varying)
coefficient regression to test and quantify the contagion of bubble behaviour
from a hypothesized "core" region to other regions. The paper applies it to
New Zealand regional house prices, with Auckland City as the core. The steps
below use the equation numbers and notation of the paper. We re-verified them
on 2026-08-10 against rendered PDF pages 17 and 26. Our earlier transcription
used `γ̂`, `β_{1j}` and `β_{2j}` throughout, whereas the paper's own symbols
are `β̂` for the AR coefficient and `δ_{1j}`/`δ_{2j}` for the
functional-regression coefficients. That was a notation slip in our earlier
write-up and not a substantive error, and it is corrected below.

1. For every region `i` and every subsample-ending date `s`, estimate a
   fixed-width rolling-window OLS AR(1) regression (their eq. 1,
   `y_t = δ + β·y_{t-1} + e_t`, a plain, non-augmented Dickey-Fuller
   regression with an intercept and one lag) to get a sequence of
   slope-coefficient estimates `β̂_{i,s}` (window width `S = ⌊0.33×T⌋` in
   their application). This is structurally close to the existing `badf`
   construction of `radf()`, a recursive AR-coefficient-type sequence per
   window, but with a fixed width in place of an expanding one. `radf()` does
   not currently support a fixed-width variant. The `O(1)`-per-window
   closed-form OLS machinery that `dating_hls()` already uses in
   `hls_segment_ssr()` and `hls_prefix_sums()` extends directly, with the same
   prefix-sum trick applied to a moving window and not an expanding one.
2. Fit the functional regression (eq. 4, p.17):
   `β̂_{j,s} = δ_{1j} + δ_{2j}·(s/(T−S+1))·β̂_{core,s−d} + error_s`, where
   `d` is an integer delay parameter in `{0,...,12}`. The text of the paper
   calls this "months", but the empirical section discusses delays in
   quarters against quarterly data. This is an inconsistency of units inside
   the paper. We confirmed it on the rendered page, so it is no transcription
   artifact. Read `d` as "native sampling periods of the input series".
3. Estimate the time-varying coefficient `δ̂_{2j}(r)` by Nadaraya-Watson
   local-constant kernel regression (eq. 6, p.26, Gaussian kernel):
   `δ̂_{2j}(r;h,d) = [Σ_s K_{hs}(r)·β̃_{j,s}·β̃_{core,s−d}] / [Σ_s
   K_{hs}(r)·β̃²_{core,s−d}]`, with `β̃_{j,s} := β̂_{j,s} − mean(β̂_{j,·})`
   (centered) and `K_{hs}(r) = (1/h)·K((s/T − r)/h)`. This is not a general
   nonparametric estimator. It is the closed-form solution of a no-intercept,
   single-regressor weighted least squares fit at each `r`, that is, a ratio of
   two weighted sums. It vectorizes fully as one `T×T` Gaussian weight matrix
   times two length-`T` vectors (`outer()` and `dnorm()`), the same
   closed-form-ratio pattern used throughout `exuber/R/`, with kernel weights
   in place of a hard window indicator.
4. Select the bandwidth `h` by leave-one-out cross-validation (eq. 7, p.26):
   `ȟ_{jT}(d) := argmin_{h∈H_T} Σ_s {β̃_{j,s} −
   δ̌_{2j}(s/(T−S+1);h,d)·β̃_{core,s−d}}²`, with `H_T =
   [(T−S+1)^(−1/2}, (T−S+1)^(−1/10)]`. This is a bounded interval and not a
   grid, and `δ̌_{2j}` is the leave-one-out version of eq. 6 (it excludes
   `p=s` from the sums). It is a textbook one-dimensional bounded
   optimization, so `stats::optimize()` works directly and no custom search is
   needed.
5. Select the delay `d ∈ {0,...,12}` by minimizing the CV-bandwidth-conditional
   SSE (eq. 8, p.26). There are 13 candidate delays, and each needs one
   `optimize()` call for `h` and then a single SSE evaluation, so this is no
   nested `bandwidth-grid × 13 × T` search. An earlier pass flagged an apparent
   inconsistency, because the prose says "choose `d` by NLS, largest R²" while
   eq. 8 says to minimize SSE. It resolves on inspection. R² = 1 − SSE/TSS, and
   TSS (the variance of the centered `β̃_{j,s}` sequence) is constant across
   candidate `d`. Maximizing R² and minimizing the SSE of eq. 8 therefore pick
   the same `d`, and the two descriptions are one criterion stated two ways.

### Exact numbers reproduced

We re-verified equations 1, 4, 6, 7 and 8 above on 2026-08-10 against rendered
PDF pages 17 and 26 (PyMuPDF), and no longer rely on the raw
`pdftotext -layout` extraction that an earlier pass had flagged for this
re-check. The structure (local-constant kernel regression, the LOOCV bandwidth
search range, and joint bandwidth and delay selection) was already right. Only
the notation (`β̂`/`δ_{1j}`/`δ_{2j}`, and not `γ̂`/`β_{1j}`/`β_{2j}`) needed
correcting, as noted above. The paper contains no numeric table of estimated
`d`, `h` or `δ̂_{2j}` values. Its results are Figures 7-8, "Time-varying
Contagion Coefficients from the Auckland City Real...", which are graphical and
not tabular. We confirmed that by grepping the full paper for every plausible
inference keyword ("confidence," "standard error," "significan*," "critical
value," "bootstrap," "asymptotic distribution"). Every hit is for the
unrelated GSADF/BSADF bubble-detection statistic (Section 2.6/5.1.1), and none
concerns the contagion coefficient. The paper performs no formal inference on
`δ̂_{2j}(r)`, and consists of point estimation, CV-based tuning and plotting
only.

### Cost/feasibility note, re-triaged 2026-08-10, materially cheaper than originally scoped

The original verdict, "most expensive item, no reuse", was wrong on two of its
three cost drivers. Only the CV-bandwidth step was assessed correctly.

1. The rolling-window AR(1) coefficient sequence is a small, mechanical
   extension of the closed-form prefix-sum window pattern that `dating_hls()`
   already uses in `hls_segment_ssr()` and `hls_prefix_sums()`. Going from an
   expanding to a moving window changes which cumulative-sum differences we
   take, in one line, and needs no new machinery.
2. The Nadaraya-Watson step is not "from-scratch nonparametric machinery".
   Re-reading eq. 6 directly shows that it is a closed-form ratio of two
   Gaussian-kernel-weighted sums. It is a WLS solution, structurally the same
   closed-form-ratio idea used throughout this project, with kernel weights in
   place of a 0/1 window indicator. It vectorizes in roughly 10-15 lines.
3. The bandwidth and delay selection does need new but standard code, which is
   the one place where the original assessment was right. It is a
   one-dimensional bounded `optimize()` call repeated 13 times (once per
   candidate delay), and not a nested two-dimensional grid search. It is cheap
   at the scale of the paper itself (`T` in the tens of quarterly
   observations).
4. No inference or critical-value machinery is needed at all, because the
   source paper does none. The earlier worry about "new simulation machinery"
   for this item was moot from the start, since there was never anything to
   port.

The minimum-viable subset mirrors the "SSU alone, no GSSU" scoping used
elsewhere in this project. It has three parts. (a) The fixed-window AR(1)
coefficient sequence. (b) A single-delay Nadaraya-Watson regression (eq. 6
only, with `d` supplied by the caller and the bandwidth chosen either by
`optimize()` over `H_T` of eq. 7 or supplied by the user as `bw`), shipped as
the point-estimate and plotting pipeline. (c) The joint delay search of eq. 8
as a thin follow-on wrapper (repeat (b) for `d = 0..12` and keep the best),
which is separable and not a prerequisite for (a)-(b). The item is tractable
now and is no longer the most expensive one in this file.

### Implementation, done (2026-08-10), minimum-viable subset

The function ships as `contagion_reg(y, core, S, d, h, r_grid)` in
`exuber/R/contagion_reg.R`. It computes the fixed-window AR(1) coefficient
sequence (eq. 1), the Nadaraya-Watson regression at a caller-supplied delay
`d` (eq. 6), and leave-one-out cross-validated bandwidth selection (eq. 7)
when `h` is not supplied. The automatic delay search of eq. 8 is not
implemented. If we need it later, we can call `contagion_reg` once per
candidate `d` and compare the results. This is the same "thin, separable
follow-on" scoping that we used for the multi-bubble DP algorithm of KNP and
the fragmentation-joining heuristic of HLW elsewhere in this project.

Brute-force validation caught two real bugs during implementation. Trusting
the transcribed formulas would not have caught either of them.

1. A window-width off-by-one. The text of GMP ("window width `S`... data
   window `{t=s-S+1,...,s}`") means `S` levels per window, which is `S-1`
   regression pairs. An initial implementation used `S` pairs (`S+1` levels)
   instead. A brute-force `lm()` cross-check at five separate window-end dates
   caught it immediately, and the results did not match to machine precision
   until we fixed it.
2. A matrix-orientation bug in the SSE of the LOOCV bandwidth. The
   Nadaraya-Watson ratio (`contagion_nw_delta2()`) correctly used
   `crossprod(K, v)` (a sum over kernel-weight rows and positions for each
   evaluation-point column) to implement the sum over `s` in eq. 6. The LOOCV
   SSE helper (`contagion_loocv_sse()`) initially used plain `K %*% v`
   instead. That is a different and wrong orientation, because the kernel
   weight matrix is not symmetric (weighing position `p` at evaluation point
   `i` differs from weighing position `i` at point `p`). A brute-force
   double-loop cross-check disagreed with the closed-form version by about 1%,
   which is more than machine-precision noise, until we traced the cause and
   fixed it.

Validation. No published numeric table exists, because the results of GMP are
Figures 7-8 and not tabulated numbers. The subagent-assisted re-triage above
grepped the full paper for every inference-related keyword and found none for
the contagion coefficient specifically. The situation is the same "figures
only" one as for
[SBZ's](/replication/volatility-robustness#why-this-cant-be-a-bit-exact-numeric-cross-check-test)
own size and power results.

- `contagion_fixed_window_beta()` (eq. 1) matches a brute-force `lm()` fit
  exactly (`< 1e-8`) at five separate window-end dates.
- `contagion_nw_delta2()` (eq. 6) matches a manual Gaussian-kernel
  weighted-least-squares ratio exactly (`< 1e-10`).
- `contagion_loocv_sse()` (eq. 7) matches a manual leave-one-out double loop
  exactly (`< 1e-8`, after the matrix-orientation fix above).
- `contagion_bandwidth_cv()` picks a bandwidth strictly inside the `H_T`
  interval of eq. 7, with LOOCV SSE no worse than at either endpoint.
- A sensible-behavior check has no ground truth to match, but Figures 7-8 of
  the paper motivate a directional check. A synthetic series whose local
  persistence tracks that of the core series (with a known delay) shows a
  visibly wider range of estimated `delta_2(r)` than an independent series
  with no relationship to the core (mean range 0.50 against 0.42 over 15
  replications). The estimator responds to a genuine time-varying
  relationship and does not return noise whatever the input.

`test-contagion.R` is new (19 tests), and the full package suite passes after
the addition. Replication script:
[replication/multivariate/radf_contagion_validation.R](#script-radf_contagion_validation).

## Summary ranking (implementation cost, cheapest first)

| Rank | Item | New statistical code needed | Source access |
|---|---|---|---|
| 1 | Chen, Phillips & Shi (PCA + PSY) | PCA call + thin wrapper; `radf()`/`radf_mc_cv()` reused as-is | Open (Cowles), theorem verified directly |
| 2 | Phillips & Yu (migration) | None, already fully covered by existing `radf()`/`datestamp()`; "migration" is narrative, not a statistic | Open (Cowles), read directly |
| 3 | Evripidou et al. (co-bubble) | New KPSS-type bivariate statistic + lead/lag search; wild-bootstrap *pattern* reusable, code is not | Now open (institutional access), not yet re-read |
| 4 | Greenaway-McGrevy & Phillips (contagion) | Re-triaged 2026-08-10: fixed-window AR sequence (small extension of `hls_segment_ssr()`'s pattern) + closed-form NW kernel ratio (~10-15 lines) + 13× 1-D `optimize()` calls for bandwidth/delay; no inference machinery needed (the paper does none) | Open (Cowles), read directly |

The ranking by engineering cost differs from the order of discovery and
citation. Once we read Phillips & Yu, it required no new statistical code,
since it is a documentation and vignette opportunity and no implementation
gap. The cost of Evripidou et al. was uncertain pending source access, and it
is now unblocked. If we order strictly by "cheapest to ship first", the order
is Chen-Phillips-Shi, then Phillips-Yu (as a vignette and not new R code),
then Evripidou (now that access is open) or Greenaway-McGrevy, whichever we
read first.
