---
title: "Multivariate bubble tests"
blurb: "Panel and cross-series tests -- common bubbles, co-bubbles, and bubble contagion."
order: 4
---
**Status: common-bubble (Chen, Phillips & Shi), co-bubble (Evripidou,
Harvey, Leybourne & Sollis), and contagion regression (Greenaway-McGrevy
& Phillips) all done — every item in this file is now implemented in
some form.** `radf_common()` is implemented in `exuber/R/radf_common.R`,
cross-checked in `exuber/tests/testthat/test-common.R`, and this file's
research (below) confirmed the *detection statistic* as genuinely
near-free to implement. Independent validation (2026-08-09, below) then
found that the critical value this file originally recommended pairing
it with (`radf_mc_cv(n, minw)`) is not actually valid at realistic panel
sizes — the correct, `N`-dependent critical value now ships as
`radf_common_cv()` (Bundle 1) — see the "Independent validation"/"Update"
subsections under Common-bubble detection for the numbers.
`radf_cobubble()` (Bundle 4, 2026-08-09) is implemented in
`exuber/R/radf_cobubble.R`, cross-checked in
`exuber/tests/testthat/test-cobubble.R` — see "Implementation" under
Co-bubble test below. `radf_contagion()` (2026-08-10, minimum-viable
subset) is implemented in `exuber/R/radf_contagion.R`, cross-checked in
`exuber/tests/testthat/test-contagion.R` — see "Implementation" under
Contagion regression below; re-triaging this item first found the
earlier "most expensive, no reuse" assessment was wrong on two of its
three cost drivers. exuber currently has no multivariate bubble test
proper *before* these; `radf()`'s `bsadf_panel`/`gsadf_panel`
(`exuber/R/radf_.R`, lines 98-99) is a *panel-average* of independently
estimated univariate BSADF sequences (`apply(bsadf, 1, mean)` then
`max()`), not a common-factor or cross-series test — it detects "bubbles
somewhere in the panel on average," not a shared latent bubble, lead/lag
co-movement, migration, or contagion. All four items below were
genuinely new capabilities relative to what existed.

| Method | Paper | Status |
|---|---|---|
| [Common-bubble detection (PCA + PSY)](#common-bubble-detection-via-pca--psy) | Chen, Phillips & Shi (2020/2023) | **done** (incl. `radf_common_cv()`) |
| [Bubble migration](#bubble-migration) | Phillips & Yu (2011) | evaluated — no new code needed |
| [Co-bubble test](#co-bubble-test) | Evripidou, Harvey, Leybourne & Sollis (2022) | **done** |
| [Contagion regression](#contagion-regression) | Greenaway-McGrevy & Phillips (2015/2016) | **done** (2026-08-10, minimum-viable subset — eq. 8's automatic delay search not included) |

All papers: [references.md](/replication/references#multivariate).

## Common-bubble detection via PCA + PSY

### Source

Chen, Y., Phillips, P.C.B. & Shi, S. "Common Bubble Detection in Large
Dimensional Financial Systems." Cowles Foundation Discussion Paper No. 2251,
Yale University, August 2020. Published version: *Journal of Financial
Econometrics*, 2023, 21(4), 989-1063.

Accessed via the open Cowles Foundation PDF — fully open, no paywall.
Downloaded and converted with `pdftotext -layout`; the theorem statements
quoted below were read directly from that extraction (no PNG rendering
needed — the equations in question are short enough that the
`pdftotext -layout` output was unambiguous on inspection; this is a
lighter verification standard than most other files in this project, noted
here for honesty).

### What it is

A two-step procedure for detecting a bubble that is common to a panel of N
series (their motivating case: real-estate prices in 89 Chinese cities):

1. **Step 1 — PCA.** Estimate the dominant common factor of the panel by
   principal components: solve `min_{Λ,F} (1/NT) Σ (X_it − λ_i f_t)²`
   subject to `(1/N)Λ'Λ = I_r` (their eq. 3.1-3.2); the estimated loadings
   are `√N` times the eigenvectors of `X'X` for the `r` largest eigenvalues,
   and the factor estimate is `F̃ = XΛ̃/N`. Only the **first** component is
   used for bubble detection (`ỹ_t`, "sufficient... for the purpose of
   bubble identification," footnote 3) — no need to determine `r` via an
   information criterion for this application.
2. **Step 2 — PSY on the factor.** Apply the standard PSY (2015a,b)
   recursive-evolving GSADF procedure to `ỹ_t`, using the *same* ADF
   regression form as exuber's existing `radf()`: `ỹ_t = μ + ρ ỹ_{t-1} + v_t`
   (their eq. 3.3 — OLS-demeaned, intercept included, not the GLS/no-intercept
   form used by `radf_tt()`).

This targets a genuinely different null/alternative than exuber's panel
average: the factor model (their eq. 2.3, 2.7) explicitly mixes an I(1)
factor (normal times), a mildly-explosive factor (bubble), and a stationary
factor (post-collapse), with idiosyncratic errors per series; the alternative
is that a *subsample* of the panel shares this explosive factor.

### Exact numbers / theorem reproduced

Verified directly from the theorem statements (not a secondary citation):

- **Theorem 4.2** (p.12 of the PDF): the DF statistic computed on the
  estimated factor `ỹ_t` (eq. 4.2, a ratio of stochastic integrals of
  standard Brownian motion `W₁`) has "limit distribution... unaffected by
  factor estimation and is identical to that of the DF statistic computed
  from the original data, as in Phillips et al. (2015a)."
- **Theorem 4.3** (p.13, eq. 4.3): the resulting "PSY-factor" statistic's
  limit is stated verbatim as *"then identical to that of the original PSY
  statistic (i.e., `F_{r2}(W, r0)` in Phillips et al. (2015a))"* — i.e. the
  claim is confirmed directly from the paper's own theorem, not inferred.
  `F_{r2}(W, r0)` is exactly the object exuber's `radf_mc_cv()` already
  simulates for `gsadf`.
- **Important caveat, not obvious from a casual read**: this identity is
  *asymptotic* (N,T → ∞). Section 5 (Simulations) reports that in finite
  samples "the finite sample distribution lies to the left of the
  asymptotic, which implies slight undersizing if asymptotic critical values
  are employed" for small N (their Figure 1, T = 60/100/140, N = 20-100),
  converging rapidly as N and T grow — **under their own DGP, which
  includes a genuine shared explosive factor**. Their own empirical
  application therefore uses **finite-sample simulated critical values**,
  not the published PSY asymptotic table. **Correction (2026-08-09,
  independent validation below): plain `radf_mc_cv(n, minw)` is *not* a
  safe substitute** — it has no dependence on `N` at all, and the
  independent-validation section below shows the true null quantiles
  diverge further from `radf_mc_cv()`'s as `N` *grows*, not shrinks. The
  originally-drafted framing here (treating `radf_mc_cv(n, minw)` as "the
  correct match") was wrong; see below for what actually needs simulating.
- Theorems 4.5-4.8 (alternative-hypothesis asymptotics, divergence rates,
  origination/collapse date consistency) were located but not read in
  depth — out of scope for this feasibility pass; flag for later if a
  faithful implementation needs the date-stamping consistency proof, not
  just the detection statistic.
- Empirical application: 89 Chinese cities, Jan 2003-Mar 2013 monthly, factor
  loadings drawn/estimated in [0.3, 1.7], idiosyncratic error SD ≈ 0.1;
  three common-bubble episodes detected in Tier 1/2 cities, none in Tier 3.
  (Reported for completeness; not independently re-derivable without the
  underlying data, so not cross-checked bit-for-bit.)

### Cost/feasibility for exuber — near-free, confirmed

This is the cheapest item in this file, and it holds up on inspection:

- **PCA step**: base R `stats::prcomp()` (or a 2-line manual eigen-decomposition
  of `X'X`) on the panel matrix `radf()` already builds via `parse_data()` —
  no new dependency, no new C++.
- **GSADF step**: literally `radf(factor_series)` — the *existing* function,
  unmodified, called on a length-T vector (the first principal component
  score). The *code* needs no changes for this — `radf()`/`radf_mc_cv()`
  don't care where their input vector came from. **But** (see caveat
  above and the independent-validation section below) `radf_mc_cv(n,
  minw)`'s critical values are only valid for `radf_common()`'s output at
  the specific `N` (panel width) they were implicitly calibrated at (in
  practice: never, since `radf_mc_cv()` has no `N` argument at all) — using
  them off-the-shelf as this file originally suggested is not sound at
  realistic panel sizes.
- **What's new, concretely**: (a) a thin wrapper function, e.g.
  `radf_common(data, r = 1, ...)`, that calls `prcomp()` on the panel, takes
  PC1's scores as a `T`-vector, and calls `radf()` + prints/returns it
  tagged as a common-bubble result; (b) optionally, date-stamping helpers
  reusing `datestamp()` unchanged since the output is a standard `radf_obj`
  on a single series once the factor is extracted.
- **What's explicitly out of scope for a first pass**: Theorem 4.5-4.8's
  origination/collapse-date *consistency* guarantees are proved for the
  factor-model DGP specifically; `datestamp()` would run mechanically on the
  factor series regardless, but claiming the paper's consistency theorem
  applies to exuber's exact date-stamping rule (`log(T)`-duration filter
  etc.) would need a closer read of Section 4.3, not done here.
- **Why it ranks first**: this is a 1-2 file addition reusing 100% of
  existing panel-parsing (`parse_data()`), estimation (`radf()`), and
  critical-value (`radf_mc_cv()`/`radf_wb_cv()`) machinery — no new
  statistic, no new asymptotic theory to implement, no new bootstrap. The
  only genuinely new code is the PCA call itself.

### Independent validation (2026-08-09) — Theorem 4.3's asymptotic identity does not hold at practical N, and the gap grows with N

The direct test of Theorem 4.3's claim: simulate a panel with **no true
common factor at all** (`N` independent random walks — the sharpest
possible null, since there's nothing for PC1 to legitimately detect), run
`radf_common()` many times, and compare its null quantiles to plain
`radf_mc_cv()` (the univariate GSADF benchmark the paper's theorem says it
should coincide with asymptotically). `T=250` throughout, seed=271828:

```r
for (N in c(6, 20, 50, 100)) {
  gsadf_null <- replicate(400, {
    panel <- replicate(N, cumsum(rnorm(250)))   # independent walks, no factor
    radf_common(panel)$gsadf
  })
  quantile(gsadf_null, c(0.9, 0.95, 0.99))
}
```

| N | 90% | 95% | 99% | diff from `radf_mc_cv()`'s 95% (2.307) |
|---|---|---|---|---|
| 6 | 2.420 | 2.662 | 3.057 | +0.36 |
| 20 | 3.225 | 3.451 | 4.144 | +1.14 |
| 50 | 4.057 | 4.365 | 4.921 | +2.06 |
| 100 | 4.880 | 5.203 | 5.836 | +2.90 |

**This is the opposite of what the "converges as N grows" reading of the
paper's Section 5 would predict**, and it isn't small: at `N=100`, the
true 95% critical value (5.2) is more than *double* `radf_mc_cv()`'s
(2.3). The mechanism is a known pitfall distinct from what the paper's own
Figure 1 documents: their finite-sample undersizing result is measured
*under their own DGP*, which bakes in a genuine shared explosive factor —
it says nothing about behavior when the "common factor" premise is false.
Extracting PC1 from a panel of purely independent (non-cointegrated) I(1)
series doesn't behave like a single random walk once `N` grows: with more
independent unit-root series to draw from, the leading principal
component increasingly picks up whatever transient co-movement exists by
chance, and that component's own persistence/apparent-explosiveness grows
with `N`. (This is a recognized general hazard of applying PCA/factor
extraction to non-cointegrated nonstationary panels — not verified against
a specific citation here, flagged as a claim about the *mechanism*, not a
sourced fact, to keep this file's usual standard of not asserting
unverified citations.)

**Practical consequence**: using `radf_mc_cv(n, minw)` — the plain
univariate critical value, which has no `N` argument — as `radf_common()`'s
critical value would make the test **badly oversized** (over-detect
"common bubbles" that aren't there) for any realistically-sized panel,
and the problem gets *worse*, not better, as the panel grows. This
directly corrects the "the correct match on the exuber side is
`radf_mc_cv(n, minw)`" claim made earlier in this section when it was
first drafted (before this validation pass) — that claim is now shown to
be wrong at any `N` beyond a handful of series.

**Update (2026-08-09, Bundle 1): implemented.** `radf_common_cv(n, N,
minw, nrep, seed)` now ships, simulating the null exactly as above (an
`N`-column panel of independent random walks, PCA, GSADF) rather than
deferring to `radf_mc_cv()` — same return shape as `radf_mc_cv()`
(`adf_cv`/`sadf_cv`/`gsadf_cv`/`badf_cv`/`bsadf_cv`), so it's a drop-in
`cv` argument for `datestamp()`/`tidy()`/`autoplot()` on a `radf_common()`
result. One implementation wrinkle: exuber's generic `is_mc()`/join
machinery (used by `datestamp()` etc.) does an *exact string match* on the
`method` attribute against `"Monte Carlo"` — `radf_common_cv()` has to set
`method = "Monte Carlo"` for that to work, so the panel-specific
distinguishing info lives in a separate `N` attribute instead of in the
method label.

Tested in `test-common.R`: runs and returns a `radf_cv`/`mc_cv`-shaped
object usable directly by `datestamp()`; and, directly re-confirming the
independent-validation finding above, `N=30`'s null quantile is
significantly higher than `N=4`'s (not the same, as a naive reading of
Theorem 4.3 might expect) — the dependence on panel width validation found
is exercised by the test suite now, not just a one-off simulation script.

Replication script:
[replication/multivariate/radf_common_validation.R](#script-radf_common_validation).

## Co-bubble test

### Source

Evripidou, A.C., Harvey, D.I., Leybourne, S.J. & Sollis, R. (2022). "Testing
for Co-explosive Behaviour in Financial Time Series." *Oxford Bulletin of
Economics and Statistics*, 84(3), 624-650. `doi:10.1111/obes.12487`.

**Status: done (2026-08-09).** Shipped as `radf_cobubble()`. Full PDF read
(abstract through Section VI — model, Theorems 1-2, the wild bootstrap
algorithm, and the lag-selection procedure), not just the abstract.

### What it is

A test for whether two series, each containing an explosive (bubble)
episode, are related via **co-explosive behaviour**: a linear combination
of the two is integrated of order zero even while each series
individually is locally explosive (analogous to cointegration, but for
explosive rather than unit-root regimes — hence "co-bubble").

The DGP (their eq. 2) is `y_t = mu_y + beta_x*x_{t-i} + beta_z*z_t +
e_y,t`, where `x_t` is observed and `z_t` is an unobserved explosive
process; under `H0: beta_x > 0, beta_z = 0`, `y_t` and `x_{t-i}` are
co-explosive. The test statistic (eq. 3) is a **KPSS-type LM statistic**
on the OLS residuals of `y_t` regressed on a constant and `x_{t-i}`:

```
e_y,t = y_t - alpha_hat - beta_hat * x_{t-i}          (OLS residual)
S = sigma_y_hat^-2 * (T-|i|)^-2 * sum_t (cumsum of e_y up to t)^2
```

over the overlapping valid range `t = max(i,0)+1, ..., T+min(i,0)`. This
is the *opposite* testing direction from PSY/ADF-style tests — the null
is stationarity (co-explosivity), not the presence of a unit root — and
the paper proves (Theorem 1) the limit null distribution does **not**
depend on the properties of the `x` regressor at all (its mild
explosivity is asymptotically negligible to this statistic), but **does**
depend on the pattern of heteroskedasticity in `e_y,t`. Since that rules
out a fixed table of critical values, a **wild bootstrap**
(`y*_t = w_t * e_y,t`, `w_t ~ IIDN(0,1)`, refit on the same `x_{t-i}`
regressor per Remark 2) reproduces the heteroskedasticity pattern and
gives asymptotically size-controlled critical values (Theorem 2).

When the lead/lag `i` is unknown (Section VI), it is estimated as
`i_hat = argmin_j sigma_y_hat(j)^2` over a candidate set of lags — a
misspecified `j != i` leaves a neglected explosive term in the residuals
that inflates their variance, so the variance-minimizing `j` consistently
recovers `i`.

### Implementation

Shipped as `radf_cobubble(y, x, lag = NULL, lags = -6:6, nboot = 499L,
level = 0.05, seed = NULL)` in `exuber/R/radf_cobubble.R`:

- `coexplosive_stat_aligned(y, xreg)` — the core KPSS-type statistic (eq.
  3) on two already-aligned, equal-length vectors.
- `coexplosive_stat(y, x, lag)` — builds the `(y_t, x_{t-lag})` pair over
  the paper's overlapping valid range and calls the aligned core.
- `coexplosive_select_lag(y, x, lags)` — Section VI's `i_hat` search.
- `radf_cobubble()` — ties it together: selects `lag` if not given, fits
  the observed statistic, runs the wild bootstrap (regressing each
  bootstrap sample on the *same* `x_{t-lag}` regressor, per Remark 2 —
  the paper explicitly found omitting it makes the bootstrap distribution
  a worse finite-sample match), and returns the statistic, critical
  value, p-value, and reject/no-reject decision.

No reuse of `radf_wb.R`'s existing wild-bootstrap DGP functions
(`radf_wb_dgp_ps()`/`radf_wb_dgp_hlst()`) — those are built around the
recursive ADF-family ptr/window structure and don't apply to a
single-shot KPSS statistic on a fixed regression; the wild-multiplier
resampling *idea* is the same, but the mechanics (one static OLS fit,
not a recursive scan) are simple enough to not need shared scaffolding.

**Independent validation**:

1. **Formula-exact check**: `coexplosive_stat()` matches an independently
   written brute-force computation (separate `lm()` call, manual
   loop-based cumulative sum instead of vectorized `cumsum()`) to
   floating-point precision (`diff ~ 1e-17`).
2. **Empirical size under H0, homoskedastic errors**: 6.0% at the nominal
   5% level (100 Monte Carlo reps) — within sampling noise.
3. **Empirical size under H0, heteroskedastic errors** (a volatility jump
   partway through the sample): **5.0%** at the nominal 5% level — this
   is the paper's central claim (Theorem 2) and the whole reason the
   wild bootstrap exists rather than a fixed KPSS table, and it holds up.
4. **Power under H1** (`y` and `x` have independent, unrelated explosive
   episodes — not co-explosive): **100%** rejection rate over 60 reps.
5. **Lag recovery**: with a true lag of 3 baked into the DGP,
   `coexplosive_select_lag()` recovered the exact true lag in 20/20 seeds.

New tests in `test-cobubble.R` (regression-tested against loose bounds on
the size/power Monte Carlo checks, matching this project's usual
tolerance-not-exact-number pattern for stochastic simulation results).

Replication script:
[replication/multivariate/radf_cobubble_validation.R](#script-radf_cobubble_validation).

### Cost/feasibility note (as originally scoped, before implementation)

This was structurally the most different of the four items in this file
from anything already in exuber — a KPSS-type (LM/stationarity)
statistic rather than an ADF/DF-family recursive one, and inherently
bivariate rather than panel-shaped. It turned out to be tractable exactly
because the statistic is closed-form per fixed `(lag, i)` — one OLS fit
and one cumulative sum, no nonparametric bandwidth selection, no
recursive window scan — so the wild bootstrap loop is just "refit and
recompute the same closed-form statistic `nboot` times," the same
computational shape as `radf_wb_cv()`'s loop even though the underlying
statistic and DGP function are unrelated to it.

## Bubble migration

### Source

Phillips, P.C.B. & Yu, J. (2011). "Dating the Timeline of Financial Bubbles
During the Subprime Crisis." *Quantitative Economics*, 2(3), 455-491.
`doi:10.3982/QE82`. Working paper: Cowles Foundation Discussion Paper No.
1770. —
fully open, no paywall.

### What it is — important finding: "migration" is not a distinct joint test

Reading the primary source directly overturns the naive framing you'd get
from the title alone. There is **no separate "migration test" statistic**
with its own null hypothesis, test equation, or critical values in this
paper. What the paper actually does:

1. Applies the same single-series recursive right-tailed unit-root
   methodology as **Phillips, Wu & Yu (2010) / Phillips & Yu (2009)** — the
   PWY-style sup-ADF test, i.e. the *predecessor* of PSY's GSADF and
   precisely the statistic exuber's `radf()` already computes (`sadf`/`badf`,
   the single-sup PWY branch, not even the double-sup GSADF) — separately
   and independently to seven unrelated financial series (Nasdaq index,
   home price index, asset-backed commercial paper, crude oil, platinum,
   Baa bond rate, Pound/USD).
2. Each series gets its own origination/collapse date estimated by the
   standard PWY dating rule (`max DF_r`/`max DF_{r,t}`, `log(n)` minimum
   duration — Table 4 in the paper reproduces exactly this format, e.g.
   "Heating oil: max DFr 6.9092, max DFrt 2.2416, origination March/08,
   collapse August/08").
3. "Migration" is a **narrative/qualitative comparison** of these
   independently-estimated dates across series after the fact: the paper
   observes that the equity-market bubble's collapse date roughly precedes
   the housing-market bubble's origination date, which roughly precedes the
   mortgage-market bubble, etc., and calls this sequence a "migration
   mechanism" — matched informally against the theoretical prediction in
   Caballero, Farhi & Gourinchas (2008), not tested with a formal joint
   statistic.

Direct quote (Conclusion, p.34-35 of the PDF): *"The dates are matched
against the onset date for the subprime crisis as well as a specific
sequential hypothesis concerning bubble migrations that are predicted in
the theoretical model proposed by CFG (2008a)."* — i.e. the "test" of the
migration hypothesis is a match against an external theoretical timeline,
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

Primary-series migration narrative (Abstract, verbatim): bubble "first
emerged in the equity market during mid-1995 lasting to the end of 2000,
followed by a bubble in the real estate market between January 2001 and
July 2007 and in the mortgage market between November 2005 and August
2007," then after the crisis erupted, migrated "selectively into the
commodity market and the foreign exchange market."

### Cost/feasibility note — effectively zero marginal statistical cost, but nothing new to build

Because the underlying per-series test is exactly what `radf()` (or even
just its `sadf`/`badf` output, ignoring the `gsadf`/`bsadf` fields PSY added
later) already computes, and "migration" is a narrative overlay rather than
a formal statistic:

- **There is no new test to implement.** Reproducing Phillips & Yu's
  "migration" analysis with exuber today would mean: run `radf()` +
  `datestamp()` on each of several series independently (already fully
  supported), then eyeball/compare the resulting date ranges — exactly what
  the paper itself does.
- **The only concrete deliverable this item could motivate** is a vignette
  or small helper for laying out multiple `datestamp()` outputs on a shared
  timeline for visual/narrative "migration" comparison (a plotting/reporting
  convenience, not a statistical addition) — see
  practitioner-guidance.md for the same idea
  applied to a more current dataset.
- Ranked below common-bubble detection for a different reason than the
  others: not because it's expensive, but because there is nothing
  statistically new to build — implementing it would just be
  documentation/vignette work demonstrating an existing capability, which
  is a product-scope decision, not a feasibility question.

## Contagion regression

### Source

Greenaway-McGrevy, R. & Phillips, P.C.B. "Hot Property in New Zealand:
Empirical Evidence of Housing Bubbles in the Metropolitan Centres." Cowles
Foundation Discussion Paper No. 2004, Yale University (2015). Published
version: *New Zealand Economic Papers*, 50(1), 88-113 (2016).
`doi:10.1080/00779954.2015.1065903`. —
fully open, no paywall.

### What it is

Section 2.6 ("Bubble Contagion") defines a **functional (time-varying)
coefficient regression** to test/quantify contagion of bubble behaviour from
a hypothesized "core" region to other regions, applied empirically to New
Zealand regional house prices with Auckland City as the core. Steps, with
equation numbers and notation as in the paper (re-verified 2026-08-10
against rendered PDF pages 17 and 26 — the earlier pass's transcription
used `γ̂`/`β_{1j}`/`β_{2j}` throughout; the paper's own symbols are `β̂`
for the AR coefficient and `δ_{1j}`/`δ_{2j}` for the functional-regression
coefficients, a notation slip in the earlier write-up, not a substantive
error — corrected below):

1. For every region `i` and every subsample-ending date `s`, estimate a
   **fixed-width rolling-window** OLS AR(1) regression (their eq. 1,
   `y_t = δ + β·y_{t-1} + e_t` — a plain, non-augmented Dickey-Fuller
   regression, intercept + one lag) to get a sequence of slope-coefficient
   estimates `β̂_{i,s}` (window width `S = ⌊0.33×T⌋` in their application).
   This is structurally close to `radf()`'s existing `badf` construction
   (a recursive AR-coefficient-type sequence per window), but fixed-width
   rather than expanding — `radf()` doesn't currently support a
   fixed-width variant, but the underlying `O(1)`-per-window closed-form
   OLS machinery already used for `radf_hls()`'s `hls_segment_ssr()`/
   `hls_prefix_sums()` extends directly (same prefix-sum trick, just a
   moving rather than expanding window).
2. Fit the functional regression (**eq. 4**, p.17):
   `β̂_{j,s} = δ_{1j} + δ_{2j}·(s/(T−S+1))·β̂_{core,s−d} + error_s`, where
   `d` is an integer delay parameter in `{0,...,12}` (the paper's own text
   calls this "months," but the empirical section discusses delays in
   *quarters* against quarterly data — a paper-internal unit
   inconsistency, confirmed on the rendered page, not a transcription
   artifact; read `d` as "native sampling periods of the input series").
3. Estimate the **time-varying coefficient** `δ̂_{2j}(r)` by Nadaraya-Watson
   local-constant kernel regression (**eq. 6**, p.26, Gaussian kernel):
   `δ̂_{2j}(r;h,d) = [Σ_s K_{hs}(r)·β̃_{j,s}·β̃_{core,s−d}] / [Σ_s
   K_{hs}(r)·β̃²_{core,s−d}]`, `β̃_{j,s} := β̂_{j,s} − mean(β̂_{j,·})`
   (centered), `K_{hs}(r) = (1/h)·K((s/T − r)/h)`. This is **not** a
   general nonparametric estimator — it is the closed-form solution of a
   no-intercept, single-regressor *weighted* least squares fit at each
   `r`, i.e. a ratio of two weighted sums, fully vectorizable as one
   `T×T` Gaussian weight matrix times two length-`T` vectors (`outer()` +
   `dnorm()`) — the same closed-form-ratio pattern already used
   throughout `exuber/R/`, just with kernel weights instead of a hard
   window indicator.
4. Select the bandwidth `h` by **leave-one-out cross-validation** (**eq.
   7**, p.26): `ȟ_{jT}(d) := argmin_{h∈H_T} Σ_s {β̃_{j,s} −
   δ̌_{2j}(s/(T−S+1);h,d)·β̃_{core,s−d}}²`, `H_T =
   [(T−S+1)^(−1/2}, (T−S+1)^(−1/10)]` — a **bounded interval**, not a
   grid; `δ̌_{2j}` is eq. 6's leave-one-out version (excludes `p=s` from
   the sums). A textbook 1-D bounded optimization, `stats::optimize()`
   directly, no custom search needed.
5. Select the delay `d ∈ {0,...,12}` by minimizing the CV-bandwidth
   -conditional SSE (**eq. 8**, p.26) — 13 candidate delays, each needing
   *one* `optimize()` call for `h` then a single SSE evaluation, not a
   nested `bandwidth-grid × 13 × T` search. The apparent internal
   inconsistency an earlier pass flagged (prose says "choose `d` by NLS,
   largest R²"; eq. 8 says minimize SSE) **resolves on inspection, not
   just a flag for later**: R² = 1 − SSE/TSS, and TSS (variance of the
   centered `β̃_{j,s}` sequence) is constant across candidate `d`, so
   maximizing R² and minimizing eq. 8's SSE pick the identical `d` — the
   two descriptions are the same criterion stated two ways, not competing
   ones.

### Exact numbers reproduced

Equations 1, 4, 6, 7, 8 above re-verified 2026-08-10 against rendered PDF
pages 17 and 26 (PyMuPDF, not the raw `pdftotext -layout` extraction an
earlier pass had flagged as needing this re-check) — the structure
(local-constant kernel regression, LOOCV bandwidth search range, joint
bandwidth/delay selection) was already right, but the notation (`β̂`/
`δ_{1j}`/`δ_{2j}`, not `γ̂`/`β_{1j}`/`β_{2j}`) is corrected above. No
numeric table of estimated `d`/`h`/`δ̂_{2j}` values exists in the paper
(results are Figures 7-8, "Time-varying Contagion Coefficients from the
Auckland City Real..." — graphical, not tabular) — confirmed by grepping
the full paper for every plausible inference keyword ("confidence,"
"standard error," "significan*," "critical value," "bootstrap,"
"asymptotic distribution"): every hit is for the unrelated GSADF/BSADF
bubble-detection statistic (Section 2.6/5.1.1), **none for the contagion
coefficient at all**. The paper performs no formal inference on
`δ̂_{2j}(r)` — point estimation, CV-based tuning, and plotting only.

### Cost/feasibility note — re-triaged 2026-08-10, materially cheaper than
originally scoped

**The original "most expensive item, no reuse" verdict was wrong on two
of its three cost drivers; only the CV-bandwidth step was assessed
correctly.**

1. **Rolling-window AR(1) coefficient sequence** — a small, mechanical
   extension of the closed-form prefix-sum window pattern already used
   for `radf_hls()`'s `hls_segment_ssr()`/`hls_prefix_sums()` (expanding
   → moving window is a one-line change to which cumulative-sum
   differences get taken), not new machinery.
2. **The Nadaraya-Watson step is not "from-scratch nonparametric
   machinery"** — re-reading eq. 6 directly shows it is a closed-form
   ratio of two Gaussian-kernel-weighted sums (a WLS solution, structurally
   the same closed-form-ratio idea used throughout this project, just
   with kernel weights instead of a 0/1 window indicator), vectorizable
   in roughly 10-15 lines.
3. **The bandwidth/delay selection genuinely needs new (but standard) code
   — the one place the original assessment was right** — but it is a
   1-D bounded `optimize()` call repeated 13 times (once per candidate
   delay), not a nested 2-D grid search; cheap at the paper's own scale
   (`T` in the tens of quarterly observations).
4. **No inference/critical-value machinery is needed at all**, because
   the source paper itself does none — the earlier assessment's implicit
   worry about "new simulation machinery" for this item was moot from
   the start; there is nothing to port because there was never anything
   there.

**Minimum-viable subset** (mirrors this project's "SSU alone, no GSSU"
scoping elsewhere): (a) the fixed-window AR(1) coefficient sequence, (b)
a single-delay Nadaraya-Watson regression (eq. 6 only, `d` supplied by
the caller, bandwidth either via `optimize()` over eq. 7's `H_T` or a
user-supplied `bw`), shipped as the point-estimate/plotting pipeline; (c)
the eq. 8 joint delay search as a thin follow-on wrapper (repeat (b) for
`d = 0..12`, keep the best) — separable, not a prerequisite for (a)-(b).
Ranked as a genuinely tractable item now, not the most expensive one in
this file.

### Implementation — done (2026-08-10), minimum-viable subset

Shipped as `radf_contagion(y, core, S, d, h, r_grid)` in
`exuber/R/radf_contagion.R`: the fixed-window AR(1) coefficient sequence
(eq. 1), the Nadaraya-Watson regression at a caller-supplied delay `d`
(eq. 6), and leave-one-out cross-validated bandwidth selection (eq. 7)
when `h` isn't supplied. Eq. 8's automatic delay search is not
implemented — call `radf_contagion` once per candidate `d` and compare
if an automatic search is needed later, the same "thin, separable
follow-on" scoping used for KNP's multi-bubble DP algorithm and HLW's
fragmentation-joining heuristic elsewhere in this project.

**Two real bugs found and fixed during implementation, both caught by
brute-force validation, neither by just trusting the transcribed
formulas**:

1. **Window-width off-by-one.** GMP's own text ("window width `S`...
   data window `{t=s-S+1,...,s}`") means `S` *levels* per window, which
   is `S-1` regression pairs — an initial implementation used `S` pairs
   (`S+1` levels) instead, caught immediately by a brute-force `lm()`
   cross-check at five separate window-end dates that didn't match to
   machine precision until fixed.
2. **Matrix-orientation bug in the LOOCV bandwidth's SSE.** The
   Nadaraya-Watson ratio (`contagion_nw_delta2()`) correctly used
   `crossprod(K, v)` (sum over kernel-weight rows/positions for each
   evaluation-point column) to implement eq. 6's sum over `s`, but the
   LOOCV SSE helper (`contagion_loocv_sse()`) initially used plain
   `K %*% v` instead — a different, wrong orientation, since the kernel
   weight matrix isn't symmetric (weighing position `p` at evaluation
   point `i` differs from weighing position `i` at point `p`). A
   brute-force double-loop cross-check disagreed with the closed-form
   version by about 1% (not machine-precision noise) until traced to
   this and fixed.

**Validated** (no published numeric table exists — GMP's own results are
Figures 7-8, not tabulated numbers, confirmed by the subagent-assisted
re-triage above grepping the full paper for every inference-related
keyword and finding none for the contagion coefficient specifically —
the same "figures only" situation as
[SBZ's](/replication/volatility-robustness#why-this-cant-be-a-bit-exact-numeric-cross-check-test)
own size/power results):

- `contagion_fixed_window_beta()` (eq. 1) matches a brute-force `lm()`
  fit exactly (`< 1e-8`) at five separate window-end dates.
- `contagion_nw_delta2()` (eq. 6) matches a manual Gaussian-kernel
  weighted-least-squares ratio exactly (`< 1e-10`).
- `contagion_loocv_sse()` (eq. 7) matches a manual leave-one-out
  double-loop exactly (`< 1e-8`, after the matrix-orientation fix above).
- `contagion_bandwidth_cv()` picks a bandwidth strictly inside eq. 7's
  own `H_T` interval, with LOOCV SSE no worse than either endpoint.
- Sensible-behavior check (no ground truth to match, but a directional
  check the paper's own Figures 7-8 motivate): a synthetic series whose
  local persistence genuinely tracks the core series' own (with a known
  delay) shows a visibly wider range of estimated `delta_2(r)` than an
  independent series with no relationship to the core (mean range 0.50
  vs. 0.42 over 15 replications) — the estimator responds to a genuine
  time-varying relationship rather than returning noise regardless of
  input.

New `test-contagion.R` (19 tests). Full package suite green after
addition. Replication script:
[replication/multivariate/radf_contagion_validation.R](#script-radf_contagion_validation).

## Summary ranking (implementation cost, cheapest first)

| Rank | Item | New statistical code needed | Source access |
|---|---|---|---|
| 1 | Chen, Phillips & Shi (PCA + PSY) | PCA call + thin wrapper; `radf()`/`radf_mc_cv()` reused as-is | Open (Cowles), theorem verified directly |
| 2 | Phillips & Yu (migration) | None — already fully covered by existing `radf()`/`datestamp()`; "migration" is narrative, not a statistic | Open (Cowles), read directly |
| 3 | Evripidou et al. (co-bubble) | New KPSS-type bivariate statistic + lead/lag search; wild-bootstrap *pattern* reusable, code is not | Now open (institutional access), not yet re-read |
| 4 | Greenaway-McGrevy & Phillips (contagion) | Re-triaged 2026-08-10: fixed-window AR sequence (small extension of `hls_segment_ssr()`'s pattern) + closed-form NW kernel ratio (~10-15 lines) + 13× 1-D `optimize()` calls for bandwidth/delay; no inference machinery needed (the paper does none) | Open (Cowles), read directly |

Note the ranking by *engineering cost* differs from discovery/citation
order: Phillips & Yu turns out to require zero new statistical code once
actually read (it's a documentation/vignette opportunity, not an
implementation gap), while Evripidou et al. was genuinely uncertain in
cost pending source access and is now unblocked. If ordering strictly by
"cheapest to ship first," the order is Chen-Phillips-Shi, then Phillips-Yu
(as a vignette, not new R code), then Evripidou (now that access is
open) or Greenaway-McGrevy, whichever gets read first.
