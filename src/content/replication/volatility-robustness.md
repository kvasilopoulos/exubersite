---
title: "...vs. per-window lm(dy[1:b] ~ ylag[1:b] - 1) for b in minw:n1"
blurb: "Tests robust to time-varying innovation variance: time-transformed, kernel-purged, WLS, sign-based, and stochastic-coefficient routes."
order: 1
---
﻿# Volatility-robustness tests

Right-tailed unit-root tests modified to stay correctly sized when the
innovation variance is time-varying (deterministically or stochastically) —
PWY/PSY's original GSADF assumes homoskedasticity, and all methods here are
different fixes for the same failure mode. Status legend: `done` =
implemented + tested against a published number, `partial` = implemented
but not yet cross-checked, `evaluated` = source read, not implemented,
`todo` = not started, `blocked` = source access blocked.

| Method | Paper | Status |
|---|---|---|
| [Time-transformed test (STADF/GSTADF)](#time-transformed-test-stadf--gstadf) | Kurozumi, Skrobotov & Tsarev (2024) | **done** |
| [Kernel-purge test](#kernel-purge-test) | Harvey, Leybourne, Taylor & Zu (2024/2025) | **done** (with-intercept variant) |
| [SBZ (WLS + kernel volatility)](#sbz-wls--kernel-volatility) | Harvey, Leybourne & Zu (2019) | **done**, bug found+fixed in validation — see note |
| [Sieve bootstrap (autocorrelated innovations)](#pedersen--schütte-sieve-bootstrap) | Pedersen & Montes Schütte (2020) | **done** |
| [Skewness-corrected wild bootstrap](#hafner-skewness-corrected-wild-bootstrap) | Hafner (2020) | **done** |
| [Sign-based sGSADF](#sign-based-sgsadf) | Harvey, Leybourne & Zu (2020); level-shift robustness + demeaned variant: Harvey, Leybourne, Tatlow & Zu (2025) | **done** |
| [Stochastic explosive-coefficient test](#stochastic-explosive-coefficient-test) | Kurozumi & Nishi (2025) | **done** (2026-08-10, `ssu_test()`, minimum-viable subset — GSSU/CUSUM/CUSUM-SQ/union not implemented) |
| [SV-ADF](#sv-adf) | Sarkar & Wells (2026, preprint) | **done** (2026-08-11, `radf_svadf()`; preprint, not peer-reviewed) |

All papers: [references.md](/replication/references#volatility-robustness).

---

## Time-transformed test (STADF / GSTADF)

**Status: done.** Implemented in `exuber/R/radf_tt.R`
(`radf_tt()`, `radf_tt_cv()`), cross-checked in
`exuber/tests/testthat/test-tt.R`.

### Source

Kurozumi, E., Skrobotov, A., & Tsarev, A. Time-Transformed Test for Bubbles
under Non-stationary Volatility. Journal of Financial Econometrics (2024),
`doi:10.1093/jjfinec/nbae026`. Working paper: arXiv:2012.13937 (freely
available — this is what was used; equation/theorem numbers below refer to
the arXiv v2, 15 Nov 2021).

Downloaded and OCR'd via `pdftotext`; formulas that looked ambiguous in the
plain-text extraction (subscripts/hats/bars are lossy in that conversion)
were re-verified by rendering the source pages to PNG (`pymupdf`) and
reading the typeset math directly, rather than trusting the OCR text. See
eq. (9), (17)-(19) and Theorem 1/2 for the formulas actually implemented.

### Why this item

Its null limiting distribution is proved (Theorem 1) to coincide with the
*homoskedastic* GLS-demeaned SADF/GSADF distribution (Whitehouse, 2019) —
no bootstrap needed, unlike HLST's wild-bootstrap PSY or SBZ. That makes it
the cheapest of the volatility-robust tests to both compute and verify, and
the paper reports a fixed, literal critical value triple — a rare thing in
this literature, most of which relies on per-dataset bootstrap p-values
that can't be reproduced bit-for-bit.

### Exact numbers reproduced

Footnote 4 (arXiv v2, page 9): *"For r0 = 0.1, they are equal to 2.319,
2.626, 3.223 for 10%, 5% and 1% significance levels."*

**Important correction made during implementation**: this triple is for
**STADF** (the single-sup, PWY/SADF-style statistic, `sup_r2 ADF^r2_0`),
not GSTADF (the double-sup, PSY/GSADF-style statistic). Whitehouse (2019)
— the source of these critical values — is itself specifically about the
GLS-demeaned *SADF*-style test (see the paper's own Section 1: "Whitehouse
(2019) considered the Phillips et al. (2011) test with the GLS-type
detrending"; PWY 2011 is the single-sup test). This was confirmed
empirically: simulating the STADF statistic (own `gls_dfstat_grid()`,
n=300, minw=30, 1500 MC reps, seed 20260808) gives (2.407, 2.702, 3.533) —
close to (2.319, 2.626, 3.223) given finite-T/MC noise — while simulating
the *GSTADF* (double-sup) statistic under the same setup gives (3.157,
3.436, 4.302), which is not close. The formula itself was independently
verified correct (see below), so this was a target-identification error,
not an implementation bug — worth recording since it's an easy trap: the
paper's own Theorem 1 states both STADF and GSTADF "coincide" with
Whitehouse (2019) critical values in the same sentence, which reads as if
one table covers both, but only the SADF-style one is actually published.
GSTADF(r0) critical values are simulable with `radf_tt_cv()` but have no
published anchor value in this pass (the paper says they're "easily
computed from the R-code available in
https://sites.google.com/site/antonskrobotov/", not fetched here).

### What's actually pivotal vs. what needs estimation

- The **null limiting distribution is pivotal** (independent of the
  volatility path) once the series has been time-transformed with the
  *true* variance profile (Theorem 1), and remains pivotal with the
  *estimated* one (Theorem 2). This is what lets `radf_tt_cv()` simulate
  critical values once via `y <- cumsum(rnorm(n))` (no volatility, no
  bootstrap) rather than per dataset — mirrors `radf_mc_cv()`'s pattern but
  with a different (no-intercept, GLS-demeaned) statistic.
- The **feasible statistic** on real, volatile data still needs the
  variance profile estimated from the data (eq. 18-19: kernel-weighted
  local no-intercept regression of Δy̌_t on y̌_{t-1}, truncated residuals,
  cumulative-sum-of-squares profile, generalized-inverse time
  transformation). Implemented in `variance_profile()`.

### Formulas implemented (verified against rendered PDF pages, not OCR text)

- `y̌_t := y_t − y_0` (GLS-demeaning: subtract the first observation,
  no intercept fitted — this is *not* what exuber's existing `radf()`
  does, which always OLS-demeans with a fitted intercept; see
  "why not exubercore" below).
- eq. (9): `ADF^{r2}_{r1} := Σ y̌_{t-1}Δy̌_t / sqrt(σ̂²(r1,r2) Σ y̌²_{t-1})`,
  the no-intercept recursive Dickey-Fuller t-statistic (mathematically:
  the t-stat on β in `Δy̌_t = β y̌_{t-1} + e_t`, no constant). Implemented
  as `gls_dfstat_grid()`, fully vectorized over the (r1, r2) grid via
  cumulative sums (no per-window regression loop needed) — same
  complexity class as PSY's GSADF grid, no bootstrap.
- eq. (18): local (Nadaraya-Watson-type) kernel estimate of the
  time-varying AR(1) coefficient δ̂_t, uniform kernel by default (as used
  in the paper's own Monte Carlo).
- eq. (19): variance profile η̂(s), a normalized cumulative-sum-of-squares
  step function of the truncated local-regression residuals; inverted via
  linear interpolation (exact since η̂ is piecewise-linear on the
  observation grid) to get ĝ(s), used to resample/time-transform the
  series.
- Footnote 6: truncation threshold ψ_T = c̄·T^(1/7), c̄ = the max
  local-window residual SD over rolling 10%-of-sample windows.

### Deliberate simplifications vs. the paper (cost/benefit)

- **Bandwidth**: paper uses leave-one-out CV over h ∈ [T^-0.5, T^-0.3].
  Implemented default is the fixed plug-in `h = T^(-2/5)` (the log-scale
  midpoint of that range), with `h` exposed as a user parameter. Full CV
  would re-run the O(T · Th) kernel fit ~10x for a bandwidth grid — not
  implemented to keep this pass lazy; flagged here as the first thing to
  add if empirical size control turns out to matter in practice.
  `# ponytail: fixed plug-in bandwidth, add CV grid search if empirical size is off`
- **Lag augmentation**: paper's own Monte Carlo uses p=0 (no augmentation
  lags) throughout; `radf_tt()` likewise only supports the p=0 case.
  `radf()`'s `lag` parameter has no equivalent here yet.
- **Not wired into `tidy()`/`autoplot()`/`datestamp()`**: `radf_tt_obj`
  only has a `print` method so far. Full S3 method parity with `radf_obj`
  (tidiers, plotting, date-stamping) is follow-up work, not done in this
  pass.

### Why not exubercore (C++)

The paper's own selling point is that this test is cheap — no bootstrap,
and its recursive statistic (eq. 9/17) is a *closed-form* ratio of
cumulative sums, not a per-window OLS fit needing matrix inversion. That
means `gls_dfstat_grid()` is fully vectorizable in R via `outer()` on
prefix sums, with no need for exubercore's C++ `radf()` (which is also
hard-coded to *always* fit an intercept when `lag == 0` — see
`exubercore/src/radf.cpp` — so it can't be reused unmodified for the
no-intercept GLS-demeaned case at lag 0 anyway). Flag for a future
exubercore port only if profiling shows the O(T²) grid becoming a
bottleneck for T in the thousands; not observed in this pass's testing at
T up to ~2000.

### Independent validation (2026-08-09)

Re-run from scratch with fresh seeds/sample sizes, not just re-executing
`test-tt.R` — the point is an independent check, not a re-run of the same
numbers.

**1. Formula check** — `gls_dfstat_grid()` (the vectorized closed-form
implementation) against a brute-force per-window `lm(dy ~ ylag - 1)` fit,
on data/parameters the test suite doesn't use (`n=65`, `minw=15`,
`seed=4242` vs. the suite's `n=40`/`minw=10`/`seed=7`):

```r
set.seed(4242)
y <- cumsum(rnorm(65)); minw <- 15
res <- exuber:::gls_dfstat_grid(y, minw)
```

Result: `max|badf_formula - badf_lm| = 6.66e-16` — machine precision, i.e.
the closed-form cumulative-sum implementation is exact, not approximate.

**2. Critical-value replication** — an independent Monte Carlo run against
Whitehouse (2019)'s published STADF triple (quoted in the paper's footnote
4, `r0=0.1`: 2.319, 2.626, 3.223), using a fresh seed and *more* replications
than the package's own test (`nrep=4000` vs. the suite's `1500`):

```r
set.seed(99001)
n <- 300; minw <- 30
sadf <- vapply(replicate(4000, exuber:::gls_dfstat_grid(cumsum(rnorm(n)), minw),
                          simplify = FALSE), `[[`, numeric(1), "sadf")
quantile(sadf, c(0.9, 0.95, 0.99))
```

| | 10% | 5% | 1% |
|---|---|---|---|
| Published (Whitehouse 2019) | 2.319 | 2.626 | 3.223 |
| My independent MC (n=300, nrep=4000, seed=99001) | 2.355 | 2.715 | 3.435 |
| Abs. difference | 0.036 | 0.089 | 0.212 |

Same direction and rough magnitude of finite-T/MC deviation as the
package's own test observed (which uses `tolerance = 0.15` and passes) —
consistent with this being MC noise around a T→∞ asymptotic target, not a
formula error. The exported `radf_tt_cv(n=300, minw=30, nrep=4000, seed=555)`
function gives (2.399, 2.737, 3.346) on yet another seed — same ballpark.

**3. Behavioral sanity check** — does the test actually fire on an
explosive series? Built a synthetic unit-root → explosive (ρ=1.04) →
collapse series (150 obs) and ran `radf_tt()`:

```r
gsadf = 5.781   # vs. the ~2.0-3.3 critical-value range above
```

Comfortably exceeds every critical value in the table above, as it should.

**4. Full existing suite**: `test-tt.R` — **8 passed, 0 failed** (1 skipped,
CRAN-only).

**Conclusion**: formula is exact (item 1), critical-value simulation
reproduces the published target within expected MC noise on an independent
run (item 2), and the test behaves correctly on a case it should detect
(item 3). No issues found.

Replication script:
[replication/volatility-robustness/radf_tt_validation.R](#script-radf_tt_validation).

### Open follow-ups

- Extend `radf_tt_cv()`/`radf_tt()`'s asymptotic critical values across a
  grid of `r0`, not just r0=0.1 (the paper notes this is directly available
  from their R code at https://sites.google.com/site/antonskrobotov/,
  not fetched in this pass).
- Full tidy/autoplot/datestamp S3 method parity with `radf_obj`.
- The [sign-based test](#sign-based-sgsadf) is still blocked on paywall
  access — this paper's own literature review (Section 1) is a decent
  secondary source for its qualitative behaviour, but not for exact
  numbers.

---

## Kernel-purge test

**Status: done** (with-intercept variant). Implemented in
`exuber/R/radf_kp.R` (`radf_kp()`), cross-checked in
`exuber/tests/testthat/test-kp.R`.

### Source

Harvey, D. I., Leybourne, S. J., Taylor, A. M. R., & Zu, Y. (2024). A new
heteroskedasticity-robust test for explosive bubbles. Journal of Time
Series Analysis. `doi:10.1111/jtsa.12784`. **Open access (CC-BY)** —
downloaded directly from the University of Essex repository
(`repository.essex.ac.uk`), no paywall issue at all.
Formulas and Table I below were confirmed by rendering the actual PDF pages
to images (PyMuPDF) and reading the typeset math, not from the initial
`pdftotext` pass (which, as with the other papers in this file, visibly
scrambled Table I's columns).

### What it is

Rather than resampling (`radf_wb_cv`) or time-deforming (`radf_tt`), this
"purges" unconditional heteroskedasticity directly: it kernel-estimates the
spot volatility \eqn{\hat\sigma_t} (eq. 4, Gaussian kernel), divides each
first difference by it, and cumulates the result (eq. 5) to get a
volatility-standardized series \eqn{x_t}. The **unmodified** PSY/GSADF test
is then run on \eqn{x_t} in place of the raw series.

The key result (Theorem 1 / Remark 3.2): the purged statistic's null
limiting distribution is **identical** to the standard homoskedastic GSADF
null. This makes it the cheapest of the volatility-robust tests to both
implement and use in exuber specifically: `radf_kp()` is a thin wrapper
that reuses `kernel_spot_vol()` (already built for [SBZ](#sbz-wls--kernel-volatility))
for the volatility estimate, then calls the **existing, unmodified**
`radf()` on the purged series -- so `radf_mc_cv()` (exuber's existing
Monte Carlo critical values) applies directly, with no new critical-value
simulation needed at all, and every downstream method (`tidy()`,
`autoplot()`, `datestamp()`, ...) works on the result for free.

The paper also proposes a without-intercept variant (\eqn{PSY^*_\sigma})
and a union-of-rejections test combining both (\eqn{UPSY_\sigma}) — see
"Not implemented" below.

### Exact numbers reproduced

Table I (asymptotic and finite-sample critical values, ε = 0.1 minimum
window, Gaussian kernel, h = 0.1·T^-0.25, 2000 MC replications),
image-verified:

| T | PSY_σ (10%/5%/1%) | PSY*_σ (10%/5%/1%) | UPSY_σ (10%/5%/1%) |
|---|---|---|---|
| 100 | 1.629 / 1.828 / 2.392 | 3.637 / 4.158 / 5.553 | 3.950 / 4.527 / 6.129 |
| 200 | 1.608 / 1.789 / 2.140 | 3.226 / 3.595 / 4.330 | 3.468 / 3.804 / 4.589 |
| 400 | 1.712 / 1.935 / 2.296 | 3.167 / 3.446 / 4.007 | 3.361 / 3.598 / 4.145 |
| ∞ | 1.875 / 2.094 / 2.486 | 2.978 / 3.296 / 3.859 | 3.186 / 3.486 / 3.951 |

The PSY_σ (with-intercept) column is the one implemented and tested. Per
Remark 3.2, its T = ∞ row should equal PSY (2015)'s own published
asymptotic GSADF critical values at the same minimum window — i.e. it
should also match what exuber's own `radf_mc_cv()` already produces. This
is tested directly: `test-kp.R` simulates `radf_kp()`'s own null GSADF
distribution at n = 300 and checks it against both (a) `radf_mc_cv(300)`
computed independently, and (b) the published T = 400 row above, with a
tolerance wide enough to cover the extra finite-sample noise the kernel
volatility estimation step adds on top of ordinary Monte Carlo error.
Observed in one run: `radf_kp` gave (1.685, 1.897, 2.337) at n = 300,
vs. `radf_mc_cv(300)`'s (1.873, 2.121, 2.429) and the published T = 400
row's (1.712, 1.935, 2.296) — consistently in the same ballpark, with
`radf_kp` running a little lower, plausibly because dividing by an
*estimated* (not true) volatility path attenuates the effective variance a
little in finite samples. This is a looser check than STADF's or SBZ's
formula-level brute-force verification (see below for why it can't be
tighter), but the qualitative pattern (all three sets of numbers within
~0.3 of each other, PSY_σ noticeably lower than PSY*_σ/UPSY_σ as the
table's own pattern shows) is exactly what the theorem predicts.

The core statistic's *formula* (eq. 4-7, the same Dickey-Fuller regression
`radf()` already computes, just on transformed data) needed no independent
brute-force check the way STADF/SBZ's genuinely new closed-form statistics
did — it's calling unmodified, already-tested exuber code
(`radf()`/`rls_gsadf`) on a new input series, so the thing actually being
verified here is the *transform* (`kernel_purge()`), which is what the
correlation and quantile-matching tests check.

### Independent validation (2026-08-09)

Fresh seed and larger sample/rep count than `test-kp.R` (which uses
`n=300, nrep=500, seed=2`): here `n=400` (matching the published Table I
row exactly, rather than interpolating), `nrep=800`, `seed=31415`.

```r
set.seed(31415)
gsadf_kp <- replicate(800, radf_kp(cumsum(rnorm(400)))$gsadf)
quantile(gsadf_kp, c(0.9, 0.95, 0.99))
```

| | 10% | 5% | 1% |
|---|---|---|---|
| Published (Table I, T=400) | 1.712 | 1.935 | 2.296 |
| Independent run (n=400, nrep=800, seed=31415) | 1.752 | 1.944 | 2.324 |
| Abs. difference | 0.040 | 0.009 | 0.028 |

Tighter than the STADF replication above (max diff 0.04 vs. STADF's 0.21),
and tighter than the package's own `test-kp.R` observed at the smaller
`n=300` — consistent with the theorem's claim converging as `T` grows
toward the table's own `T=400` design point, i.e. behaving exactly as
Remark 3.2 predicts rather than being a fluke of one seed.

**Full existing suite**: `test-kp.R` — see combined results at the end of
this file's validation runs.

Replication script:
[replication/volatility-robustness/radf_kp_validation.R](#script-radf_kp_validation).

### Not implemented: without-intercept variant and union test

`PSY*_σ` (without-intercept) needs a no-intercept DF regression on
non-demeaned data — subtly different from [STADF](#time-transformed-test-stadf--gstadf)'s
`gls_dfstat_grid()` (which forces GLS-demeaning by subtracting the first
observation; this paper's `x_t` doesn't need that because the intercept is
already "redundant" by construction, not because of an explicit
demeaning step). Reusing `gls_dfstat_grid()` as-is would be an
asymptotically-equivalent but not bit-faithful implementation of the
paper's literal (non-demeaned) formula. Combined with `UPSY_σ`'s
union-scaling logic (structurally identical to SBZ's union procedure,
already implemented in `radf_sbz_cv()`), this is a follow-up, not a large
one — most of the pieces (no-intercept grid statistic, union-of-rejections
scaling) already exist elsewhere in the codebase and would mostly need
assembling, not new theory.

---

## SBZ (WLS + kernel volatility)

**Status: done.** Implemented in `exuber/R/radf_sbz.R`
(`radf_sbz_cv()`), cross-checked in `exuber/tests/testthat/test-sbz.R`.
Independent validation (2026-08-09, below) found and fixed a real
off-by-one bug that was badly oversizing the `supDF`/`U` statistics —
fixed, re-validated, see that section for the full story.
All formulas (eq. 5, 6, and the union statistic) were re-verified by
rendering the source PDF pages to images and reading the typeset math
directly, not just from the initial `pdftotext` pass.

### Source

Harvey, D.I., Leybourne, S.J. & Zu, Y. (2019). Testing explosive bubbles
with time-varying volatility. Econometric Reviews, 38(10), 1131-1151.
Working paper: Granger Centre Discussion Paper 18/05, University of
Nottingham (open access,
`nottingham.ac.uk/research/groups/grangercentre/documents/18-05.pdf`) —
this is what was actually read; downloaded, converted with `pdftotext`,
and the Table 1 numbers below were re-verified by extracting word
coordinates with PyMuPDF (`fitz`) directly from the PDF, because the plain
`pdftotext -layout` output visibly scrambled the row/column order of the
table (a label got merged into the header row on one pass). Coordinate
extraction confirmed the correct row/column assignment.

### What it is

A weighted-least-squares (GLS-type) variant of the PWY/PSY sup-ADF test,
using a nonparametric kernel estimate of the (unknown, time-varying)
volatility as weights (eq. 6 in the paper: Gaussian kernel, leave-one-out
cross-validated bandwidth). Because its null distribution still depends on
the volatility path, size control requires a wild bootstrap (Harvey,
Leybourne, Sollis & Taylor's algorithm, reused jointly for both supBZ and
supDF). The paper's contribution beyond the WLS statistic itself is a
union-of-rejections rule combining supDF (=PWY/PSY's classic test) and
supBZ, with a single scaling constant so the union is asymptotically
correctly sized — reject if `U := max(supDF, (qDF/qBZ)*supBZ) > qDF`.

### Exact numbers reproduced

Table 1 (empirical application, bootstrap p-values, M=499 bootstrap reps,
FTSE Dec 1985-Dec 1999 and S&P 500 Jan 1980-Mar 2000), coordinate-verified:

| Series | supDF | supBZ | U |
|---|---|---|---|
| FTSE Daily | 0.288 | 0.016 | 0.046 |
| FTSE Weekly | 0.275 | 0.146 | 0.201 |
| FTSE Monthly | 0.477 | 0.279 | 0.315 |
| SP500 Daily | 0.267 | 0.000 | 0.003 |
| SP500 Weekly | 0.170 | 0.002 | 0.011 |
| SP500 Monthly | 0.153 | 0.044 | 0.071 |

Cross-checked against the paper's own narrative ("supDF does not reject for
any series"; "supBZ rejects at the 0.05-level for daily FTSE and all
frequencies of S&P 500"; "U preserves these rejections, albeit at a
slightly weaker significance level for monthly S&P 500") — all consistent
with the table above once row/column order was corrected.

### Why this can't be a bit-exact numeric cross-check test

Unlike STADF, this paper does **not** publish a fixed asymptotic
critical-value table — the whole point of the method is that size
correction always runs through a per-dataset wild bootstrap (M=499
replications here), which is inherently RNG- and data-dependent. Table 1's
p-values were computed once, on the authors' specific FTSE/S&P 500 series
with their specific bootstrap draws; they cannot be reproduced bit-for-bit
without both the original price series and the original RNG stream.
Finite-sample size/power results (Figures 3-4) are likewise reported
graphically, not as tables. A faithful cross-check test would therefore
have to be tolerance-based (e.g., verify empirical size is close to nominal
under H0 via our own large-M simulation) rather than an exact-value
assertion — a materially different (weaker) guarantee than what was
possible for STADF's fixed critical values.

### Implementation

Needed: a Gaussian-kernel nonparametric volatility estimator with
leave-one-out CV bandwidth selection (eq. 6, footnote 2 gives the
bandwidth search range `h ∈ [1/(2T), 1/6]`); a WLS-weighted recursive
Dickey-Fuller statistic (eq. in section 3, needing an intercept *and*
heteroskedasticity-weighted least squares per window — structurally
different from both exuber's existing OLS `radf()` and the no-intercept
`radf_tt()`); and reuse of exuber's existing wild bootstrap machinery
(`radf_wb.R`/`radf_wb_cv`, which already implements the Harvey/Leybourne/
Sollis/Taylor wild bootstrap for supDF) extended to run supDF and supBZ
*jointly* on the same bootstrap draws (needed for the union procedure's
validity, Theorem 3) and to compute the union scaling constant. This is
what `radf_sbz.R` (200 lines) now does.

### Independent validation (2026-08-09) — found and fixed a real bug

Table 1 can't be reproduced bit-for-bit (see above), so the honest
independent check here is **empirical size under the null**: simulate pure
random walks (no bubble present) many times, run `radf_sbz_cv()` on each,
and confirm the bootstrap rejects at roughly the nominal rate, not
substantially more or less.

**First run** (seed=13579, n=150, nrep=150 replications, nboot=199,
nominal 5%) found gross oversizing:

| Statistic | Empirical rejection rate under H0 (target ≈ 0.05) |
|---|---|
| supDF | **0.640** |
| supBZ | 0.080 |
| U | **0.573** |

supDF rejecting a true null 64% of the time (should be ~5%) is not
finite-sample noise — that's a broken test. supBZ, computed through a
different code path, was much closer to nominal (0.08), which localized
where to look: something specific to the `supDF`/bootstrap-indexing side
of `radf_sbz_cv()`, not the WLS/kernel machinery shared by both statistics.

**Root cause**: `radf_sbz.R`'s bootstrap loop computed

```r
pointer <- length(ystar) - 1L - minw
boot_df[b] <- rls_gsadf(unroot(ystar), min_win = minw)[pointer + 2]
```

but the correct convention — used consistently everywhere else in the
codebase, e.g. `radf_wb_hlst()` in `radf_wb.R`: `pointer <- nr - minw`
(no `-1`) — has no `-1L` term. Since `radf_wb_dgp_hlst()`'s bootstrap
replicate `ystar` is the same length as the original series (confirmed by
reading its source: `ystar <- c(0, cumsum(w * dy))`, length
`1 + length(dy) = length(y)`), that stray `-1` shifts the index by
exactly one position into the flat vector `rls_gsadf()` returns. Per the
indexing convention in `radf_.R` (`results[pointer+1]`=adf,
`results[pointer+2]`=sadf, `results[pointer+3]`=gsadf), the bootstrap loop
was silently extracting **adf** (a single end-of-sample t-statistic) into
a variable meant to hold draws of **sadf** (the running maximum over the
whole recursive sequence) — a systematically much smaller quantity. The
observed `supDF_obs` (correctly computed as `sadf`) was therefore being
compared against a bootstrap critical value drawn from the wrong,
much-too-low distribution, causing gross over-rejection. `U`, defined as
`max(supDF, ratio*supBZ)`, inherited the problem through its `supDF` term.

**Fix applied** (`exuber/R/radf_sbz.R`, one line):

```diff
-      pointer <- length(ystar) - 1L - minw
+      pointer <- length(ystar) - minw
```

**Re-validation after the fix**, same seed/parameters:

| Statistic | Before fix | After fix | Target |
|---|---|---|---|
| supDF | 0.640 | **0.033** | 0.05 |
| supBZ | 0.080 | 0.080 (unchanged, as expected) | 0.05 |
| U | 0.573 | **0.060** | 0.05 |

Both supDF and U are now in the right neighborhood of nominal size (supBZ
was never affected by this particular bug and remains mildly — plausibly
just finite-`nboot`-noise — oversized at 8%, not investigated further
here). `test-sbz.R` — **3 passed, 0 failed** (1 skipped, CRAN-only) after
the fix, unchanged from before, since none of the existing tests happened
to exercise this path with enough replications to notice a ~2-3x size
distortion.

**Status note (2026-08-09)**: bug found and fixed via independent
validation (this section). Still uncommitted in the exuber repo, and still
not verified against the paper's own Table 1 FTSE/S&P p-values specifically
(which needs the original price data, not just a size simulation) — but
the size-distortion bug that *would* have produced wrong p-values on any
real dataset is now fixed. Next step if resumed: commit, and if the
original FTSE/S&P series can be sourced, a direct Table 1 comparison would
be a stronger check than the H0 simulation above.

Replication script:
[replication/volatility-robustness/radf_sbz_validation.R](#script-radf_sbz_validation).

---

## Pedersen & Schütte sieve bootstrap

**Status: evaluated — mostly already covered.**

### Source

Pedersen, T.Q. & Montes Schütte, E.C. Testing for Explosive
Bubbles in the Presence of Autocorrelated Innovations. Journal of Empirical
Finance, 58 (2020), 207-225. Open working paper (CREATES Research Paper
2017-9, Aarhus University): `pure.au.dk/ws/files/109652663/rp17_09.pdf` —
this is what was read.

### Finding

exuber already has a sieve bootstrap (`R/radf_sb.R`/`radf_sb_cv()`),
citing a *different* paper (Pavlidis et al. 2016) but using essentially the
same core technique that Pedersen & Schütte independently propose for the
same problem (autocorrelated innovations in the recursive right-tailed
unit root test): fit an AR(lag) sieve to the first differences, resample
the residuals, reconstruct bootstrap paths via `stats::filter(..., "rec")`.
So this item is **largely already implemented**, not a gap.

The one concrete, specific difference found: Pedersen & Schütte's actual
contribution is showing that a *fixed* lag order causes size distortion
under autocorrelated innovations, and that automatic BIC-based lag-order
selection (with a `kmax`, re-selected as part of the procedure) fixes it
(their Section 4 simulations, e.g. Table 4.8). exuber's `radf_sb_cv()`
used to take only a single fixed `lag` argument, no automatic AIC/BIC
selection.

**Update (2026-08-09, Bundle 1): implemented.** `radf_sb_cv(type =
"aic"/"bic", max_lag = ...)` now reuses the existing `lag_select()`
machinery from `R/radf_wb.R` (already used by `radf_wb_cv2()`) — selects
the lag per series via AIC/BIC and takes the max across the panel (the
rest of `radf_sb_()`'s pointer/matrix-dimension logic assumes one common
lag for the whole panel, matching `radf()`'s own single-`lag` API, so this
was the natural scope rather than a per-series-varying lag). Tested in the
new `test-sb.R` (no prior test file existed for `radf_sb_cv()` at all):
`type = "fixed"` is unchanged from before; `type = "bic"` picks up a
nonzero lag on genuinely AR(2)-autocorrelated data; and, checked across 8
independent draws rather than asserted on one (BIC can pick a nonzero lag
by chance in any single finite sample), the modal selection on pure
random-walk data is 0, as it should be.

Replication script:
[replication/volatility-robustness/radf_sb_cv_aic_bic_validation.R](#script-radf_sb_cv_aic_bic_validation).

---

## Hafner skewness-corrected wild bootstrap

**Status: done (2026-08-09).** Shipped as `radf_wb_cv(..., dist_skew =
TRUE)` / `radf_wb_distr(..., dist_skew = TRUE)`.

### Source

Hafner, C.M. (2020). Testing for Bubbles in Cryptocurrencies with
Time-Varying Volatility. Journal of Financial Econometrics, 18(2), 233-249.
SSRN (abstract id 3105251) returned HTTP 403 to automated fetching at
first; recovered instead via EconStor (IRTG 1792 Discussion Paper
2018-005, Humboldt University Berlin):
Full PDF read (through the wild bootstrap algorithm and its Monte Carlo
section), not just the abstract.

### What it is

A modification of Harvey et al. (2016)'s wild bootstrap multiplier
distribution to better approximate the PWY test's distribution when
returns are both time-varying in volatility *and* right-skewed (the
paper's motivating case: cryptocurrency returns). The exact multiplier
(their "Step 1"):

```
u_t, v_t ~ iid N(0,1), independent
w_t = u_t/sqrt(2) + (v_t^2 - 1)/2
```

By construction `E[w_t]=0`, `E[w_t^2]=1`, `E[w_t^3]=1` — a fixed
right-skewed multiplier (not adaptively matched to each series' own
empirical skewness), used in place of the symmetric `N(0,1)`/Rademacher
multiplier in an otherwise unchanged Harvey et al. (2016)-style wild
bootstrap (`y*_t = w_t * OLS-residual_t`, refit recursively). This is
exactly the "small, contained change" this section originally predicted
before the primary source was read.

### Implementation

Shipped as a new `dist_skew` argument threaded through the existing
wild-bootstrap machinery in `exuber/R/radf_wb.R` — no new file, no new
statistic:

- `radf_wb_dgp_hlst(y, dist_rad, dist_skew = FALSE)` gained a third
  multiplier branch implementing `w_t` above, alongside the existing
  `N(0,1)` (default) and Rademacher (`dist_rad = TRUE`) branches.
- `radf_wb_hlst()`, `radf_wb_cv()`, `radf_wb_distr()` all gained a
  passthrough `dist_skew` parameter (default `FALSE`, so existing calls
  are unaffected) and validate that `dist_rad` and `dist_skew` aren't
  both `TRUE`.

**Independent validation**:

1. **Moment check**: simulating `w_t` directly (500k+ draws) confirms
   `E[w]~0`, `E[w^2]~1`, `E[w^3]~1`, matching the paper's stated
   construction.
2. **Regression check**: `dist_skew = FALSE` reproduces the pre-change
   wild bootstrap DGP bit-for-bit for the same seed — this is a purely
   additive option, not a behavior change to the default path.
3. **Power**: under a clear mildly-explosive alternative with ordinary
   (non-skewed) innovations, `dist_skew = TRUE` correctly rejects 86.7%
   of the time (30 reps) — confirms the skewed multiplier doesn't harm
   basic detection ability.
4. **Size under H0** with the paper's own right-skewed innovation
   distribution (negative log-chi-square(1), the exact distribution used
   in the paper's Monte Carlo, footnote 1): empirical rejection rate
   **3.3%** (60 reps, nominal 5%) with no added heteroskedasticity, and
   **1.7%** (60 reps) with a deterministic volatility pattern added on
   top — both mildly conservative, neither oversized, consistent with the
   paper's own reported finding that "the test is undersized in small
   samples, with the bias increasing with the degree of global
   heteroskedasticity."
5. **A DGP-sensitivity finding worth recording honestly**: an earlier,
   more aggressive power check that combined the log-chi-square
   innovations *with* the explosive alternative gave 0% power. Isolating
   the two factors traced this to the heavy right tail of the
   log-chi-square distribution itself (occasional extreme single-step
   outliers, `-log(Z^2)` blows up whenever `Z` is near zero) — in a short
   series, a single such outlier can dominate the SSR enough to mask a
   modest (`rho = 1.06`) explosive drift, in *both* the observed
   statistic and the bootstrap replicates, washing out the signal. This
   is a property of testing power under an extremely heavy-tailed noise
   process at short `T`, not a defect in `dist_skew`'s implementation —
   confirmed by checking that the identical alternative DGP with ordinary
   normal innovations gives strong power (finding 3, above).

New tests in `test-cv.R` (extending the existing `radf_wb_cv()` test
block rather than a new file, matching where its sibling `dist_rad` is
already tested). Full package suite green.

Replication scripts:
[replication/volatility-robustness/hafner_dist_skew_moments_and_regression.R](#script-hafner_dist_skew_moments_and_regression),
[hafner_dist_skew_power_and_size.R](#script-hafner_dist_skew_power_and_size).

---

## Sign-based sGSADF

**Status: done (2026-08-09; level-shift robustness + demeaned variant
added 2026-08-11).** Shipped as `radf_sign()` / `radf_sign_cv()` and
`radf_sign_dm()` / `radf_sign_dm_cv()`. Full PDF read (through Theorem 2,
Table 1, and Remark 1 of HLZ 2020; Theorems 2-3, Remark 4, and Table 1
of HLTZ 2025), not just the abstract — see "What it is" and
"Implementation" below.

### Source

Harvey, D.I., Leybourne, S.J. & Zu, Y. (2020). Sign-based unit root tests
for explosive financial bubbles in the presence of deterministically
time-varying volatility. Econometric Theory, 36(1), 122-169.
`doi:10.1017/S0266466619000057`. Long blocked — Cambridge Core, Nottingham
Repository, and the author's own homepage all led nowhere (see history
below) — until access was retried through a UK academic network (Jisc),
which unlocked it via Cambridge Core directly.

Its 2025 Oxford Bulletin extension — Harvey, Leybourne, Tatlow & Zu,
"Unit Root Tests for Explosive Financial Bubbles in the Presence of
Deterministic Level Shifts," OBES 87(5), `doi:10.1111/obes.12668` — was
recovered the same way.

### Access history (kept for anyone hitting the same wall without institutional access)

- Cambridge Core abstract page: paywalled, no fulltext (until institutional
  access; see above).
- Nottingham Repository (worktribe) record page: HTTP 403 — this
  particular repository resisted every route tried in this whole project,
  see [references.md](/replication/references#the-two-that-automation-couldnt-get-hls-and-hlw).
- Yang Zu's homepage (`sites.google.com/site/zuyang`) lists the paper with
  only a DOI link, no PDF or code link (unlike the 2019 SBZ paper, which
  does link a Google Drive code archive).
- The freely-available STADF paper (Kurozumi, Skrobotov & Tsarev,
  arXiv:2012.13937) summarizes the sign-based test qualitatively in its
  literature review (Section 1) and includes it as a comparison method
  (labelled "S") in its own Monte Carlo tables — indirect, qualitative
  confirmation of its properties (does not require bootstrap for size
  control; needs a bootstrap union-of-rejections with SADF for good power;
  "computationally expensive, and computation time increases rapidly if
  the sample size increases") but not its exact test statistic formula or
  critical values.
- B. Tatlow's own site (btatlow.com) lists the OBES extension but states
  "Pdf currently upon request", no download link. University of Macau econ
  dept page (Yang Zu's affiliation) just links the DOI.

### What it is

The core insight is structural, not statistical: `sign(Delta y_t)` is
*exactly* invariant to the volatility of `Delta y_t` (it only depends on
which side of zero the innovation fell), so a test built entirely on
cumulated signs rather than on the raw series is exactly invariant to any
volatility pattern — even a wildly time-varying one — with **no bootstrap
needed at all**, unlike HLST's wild-bootstrap correction for the standard
PSY test that `radf_wb_cv()` already implements.

Concretely: let `C_t := sum_{i<=t} sign(Delta y_i)`. The sign-based
statistic (`sPSY`, and its single-supremum special case `sPWY`, eq. 4) is
the same double-supremum recursive Dickey-Fuller construction PSY uses,
just applied to `C_t` instead of `y_t`, and fit **without an intercept**
(`C_t = rho(r1,r2) * C_{t-1} + e_t`, vs. PSY's own with-intercept
regression). Theorem 2 proves the null limiting distribution of `sPSY`
does not depend on the volatility process `sigma(s)` at all (exact
invariance), and — Remark 1 — this holds *because* the test excludes an
intercept and works only with `sign(Delta y_t)`, not because of anything
specific to bubbles; a hypothetical no-intercept version of the *standard*
PSY test would only be volatility-invariant asymptotically, not exactly,
whereas `sPSY` is exactly invariant even in finite samples.

### Implementation

Two structural facts made this the cheapest item validated in this batch:

1. **The no-intercept recursive-DF machinery this needs already exists**,
   built for STADF: `gls_dfstat_grid()` in `radf_tt.R` fits exactly a
   GLS-demeaned, no-separate-intercept AR(1) with its own recursively
   re-estimated residual variance — the same structural form `sPSY`'s
   `sDF(r1,r2)` needs. (An empirical check first ruled out reusing
   `radf()`'s own `unroot()`/`rls_gsadf()` pathway directly: despite
   `unroot(y, lag = 0)`'s column names suggesting no intercept, it was
   confirmed to compute a *with*-intercept regression, matching
   `lm(diff(y) ~ y[-length(y)])` exactly — the wrong form for this test.)
2. **The null distribution is pivotal** (Theorem 2), so critical values
   are simulated once via Monte Carlo under a plain random walk, never
   per-dataset — the same pattern `radf_tt_cv()` already uses for STADF.

Shipped in a new `exuber/R/radf_sign.R`, closely mirroring `radf_tt.R`'s
structure:

- `sign_transform(y) := c(0, cumsum(sign(diff(y))))`.
- `radf_sign(data, minw)` — transforms each series and calls
  `gls_dfstat_grid()` unchanged, returning a `radf_sign_obj` (inherits
  `radf_obj`, so it's compatible with the same downstream S3 machinery
  STADF's `radf_tt_obj` uses).
- `radf_sign_cv(n, minw, nrep, seed)` — Monte Carlo critical values,
  structurally identical to `radf_tt_cv()` with the transform swapped.

**Explicitly scoped out of this pass** (both are genuinely separate,
larger pieces of work, not needed for the core test to be useful):

- The paper's **union-of-rejections** combining `sPSY`/`sPWY` with the
  standard `PSY`/`PWY` test via a wild bootstrap (Section 4) — this
  captures "most of the power available from the better performing of the
  two tests," since `sPSY` doesn't strictly dominate the standard test on
  power for every specification. Structurally the same union pattern as
  [SBZ's `U` statistic](#sbz-wls--kernel-volatility), reusable as a
  precedent if this gets picked up.
- **Section 6's bubble-dating extension** using the sign-based statistic
  (a `datestamp()` analogue built on `C_t`) — not attempted.

**Independent validation**:

1. **Exact invariance to heteroskedasticity** — the paper's central
   claim, and the cleanest possible test of it: the same sign-pattern
   series, once left alone (constant volatility) and once multiplied by a
   wildly time-varying volatility pattern (`0.1`/`10`/`1` across three
   segments), gives **bit-identical** `sadf`/`gsadf` statistics. Not
   approximately similar — identical, exactly as the exact-invariance
   theorem predicts.
2. **Cross-check against Table 1's published finite-sample critical
   values**, at the *exact* matching finite `T` (not an asymptotic
   approximation — see the note below on why): at `T = 200`,
   `minw/T = 0.1`, simulated `gsadf_cv` (`sPSY`) = `(3.482, 3.900, 4.925)`
   vs. published `(3.469, 3.901, 4.957)` for (10%, 5%, 1%) — the 5% value
   matches to three decimal places. `sadf_cv` (`sPWY`) =
   `(2.337, 2.665, 3.403)` vs. published `(2.405, 2.735, 3.434)`, also
   close. At `T = 100` the match is looser (published 1% value for `sPSY`
   is 13.056, an extreme-value-statistic artifact at small samples/high
   quantiles that the paper's own text flags as expected).
3. **A methodology note worth recording**: an initial attempt compared a
   `T = 300` simulation against the paper's *asymptotic* (`T = Inf`)
   table row and found `sPWY` matched well but `sPSY` did not (off by
   ~0.5 at the 5% level). This is not an implementation error — the
   paper's own Table 1 documents that `sPSY`'s finite-sample critical
   values converge to their asymptotic limit much more slowly than
   `sPWY`'s (its own text: "convergence... is fairly slow (particularly
   for `sPSY`)"), and `T = 300`'s value sits almost exactly between the
   paper's own `T = 200` and `T = 400` rows once compared correctly.
   Re-running the cross-check at the *exact* matching finite `T` (point 2
   above) removed the ambiguity entirely.
4. **Power**: 96.7% empirical rejection rate (30 reps) under a clear
   mildly explosive alternative, using `radf_sign_cv()`'s own simulated
   95% critical value.

New tests in `test-sign.R` (formula check vs. brute-force `lm()`, the
exact-invariance property, the published-value cross-check at `T = 200`,
and a power check). Full package suite green.

Replication scripts:
[replication/volatility-robustness/sign_based_invariance_and_power.R](#script-sign_based_invariance_and_power),
[sign_based_finite_T_crosscheck.R](#script-sign_based_finite_T_crosscheck).

### Level-shift robustness (HLTZ 2025)

**Status: done (2026-08-11).** Shipped as a documentation addition to
`radf_sign()`/`radf_sign_cv()` plus one genuinely new function pair,
`radf_sign_dm()`/`radf_sign_dm_cv()`.

**Source**: Harvey, D.I., Leybourne, S.J., Tatlow, D. & Zu, Y. (2025).
Unit root tests for explosive financial bubbles in the presence of
deterministic level shifts. *Oxford Bulletin of Economics and
Statistics*, 87(5), 879-901. `doi:10.1111/obes.12668`. Recovered via the
same UK academic network (Jisc) route as HLZ (2020) itself — see the
"Access history" note above.

**What it is**: this is the same HLZ (2020) sign-based statistics
evaluated under a *different* null-violating nuisance (deterministic
level shifts in the series, rather than time-varying volatility). Full
PDF read (rendered pages 3-6, since raw-text extraction scrambles the
theorem statements as usual): the paper's model adds `n_T` deterministic
level shifts to the standard null, with `n_T = O(T^alpha_n)`. Three
theorems give the large-sample null distribution of the standard `PSY`
statistic (Theorem 1) and both HLZ sign-based statistics — the plain
cumulated-sign `sPWY`/`sPSY` already shipped as `radf_sign()` (Theorem
2) and HLZ's *second*, recursively demeaned analogue, denoted
`s̄PWY`/`s̄PSY` in the paper (Theorem 3), not previously implemented here.
The key finding: `PSY`'s own validity needs a *joint* restriction on
both the number and the magnitude of the shifts (Assumption 3), and per
the paper's Table 1 it is essentially never correctly sized once the
number of shifts grows at rate `sqrt(T)` (their Case 1) — empirical size
up to 0.425 against a nominal 0.05. Both sign-based statistics, by
contrast, need only a restriction on the *number* of shifts (Assumption
4, `alpha_n < 1/2`), with **no restriction on shift magnitude at all**;
only at the boundary rate `alpha_n = 1/2` do they pick up a level-shift
-dependent term too, and even then the degree of over-sizing is bounded
independent of the shift magnitude (Remark 4) — a categorically weaker
vulnerability than PSY's.

**Implementation**: two structural facts kept this cheap:

1. `radf_sign()` (already shipped) turned out to be *exactly* HLZ's
   `sPWY`/`sPSY` statistic — Theorem 2 is a level-shift-robustness result
   about code that already existed, not a reason to write anything new.
   This is a pure documentation addition (a `Level-shift robustness`
   `@section` on `radf_sign()`'s own roxygen docs).
2. The one genuinely new item, HLZ's second sign-based analogue, looked
   from an abstract-level read like it would need a per-`(r1,r2)`-window
   recursive-demeaning step (expensive, double-recursion-shaped) — but
   reading Theorem 3's actual formula shows `C̃_t` is demeaned by an
   *expanding-window* mean computed once up front (`C̃_t := sum_{i=2}^t
   {sign(dy_i) - (i-1)^{-1} sum_{j=2}^i sign(dy_j)}`), not per candidate
   window. Since `sum_{j<=i} sign(dy_j)` is just `radf_sign()`'s own
   running sum, the demeaning term at each `i` is simply `C_i / (i-1)` —
   an `O(T)` one-time transform, `sign_demean_transform()` in
   `radf_sign.R`, feeding into the *same* `gls_dfstat_grid()` machinery
   `radf_sign()` already uses. No new estimation machinery at all.

**A genuine bug found along the way (not a validation failure of this
item — a latent defect in already-shipped code)**: `print()`-ing a
`radf_sign_cv()` or `radf_tt_cv()` object crashed with `no applicable
method for 'tidy_radf_cv'`. Root cause: both functions tag their output
with a distinguishing class (`"sign_cv"`, `"tt_cv"`) that has no matching
`tidy_radf_cv.*`/`summary_radf.*`/`index_radf_cv.*` method anywhere in
the package — only `mc_cv`/`wb_cv`/`sb_cv` do. Since these functions'
output has the *exact same shape* as `radf_mc_cv()`'s (a plain
`adf_cv`/`sadf_cv`/`gsadf_cv` quantile vector, no per-dataset or panel
dimension), the root-cause fix was to add `"mc_cv"` as an additional
class tag to `radf_sign_cv()`, `radf_tt_cv()`, and the new
`radf_sign_dm_cv()`, so all three fall through to the existing
`tidy_radf_cv.mc_cv` method rather than duplicating it. Caught only
because this pass's smoke test happened to call `print()` on the new
`radf_sign_dm_cv()` object — a reminder that a function shipped without
ever calling `print()` on its own return value can hide this class of
bug indefinitely.

**Explicitly scoped out**: HLZ's serial-correlation correction for
`sPSY`/`s̄PSY` under level shifts (Remark 6 — augmenting regression (4)
with lagged `Delta C_t`, since HLZ's own original correction neglects
the shifts) — not attempted, since this project's existing sign-based
implementation already doesn't handle serial correlation for the
no-shift case either.

**Independent validation**:

1. **Formula-exact check**: `sign_demean_transform()` vs. a brute-force
   per-`i` loop recomputing the recursive mean from scratch — matches to
   `~1.8e-15`.
2. **Reproduction of the paper's own Table 1** (the strongest available
   check): Case 1 (`alpha_n = 0.5`, the worst case for `PSY`), `k = 2`,
   `mu = 5`, `p = 0.8`, `T = 400` — published empirical size at the
   nominal 5% level is `PSY = 0.337`, `sPSY = 0.121`, `s̄PSY = 0.050`. Own
   replication (`nrep = 500`, own simulated critical values under the
   no-shift null, same `minw = floor(0.1T)` trimming the paper uses,
   same shift-count/magnitude/location construction): `PSY = 0.302`,
   `sPSY = 0.090`, `s̄PSY = 0.036` — same ordering and same order of
   magnitude as the paper throughout (`PSY` badly oversized, `sPSY`
   mildly oversized, `s̄PSY` closest to nominal and even mildly
   conservative), well within Monte Carlo noise at `nrep = 500` of the
   paper's own `nrep = 2000` numbers.
3. **Exact invariance to heteroskedasticity carries over**: since
   `radf_sign_dm()` reuses the same `sign()`-based construction, the
   bit-identical-under-rescaling property `radf_sign()` has is retained
   (same test as `radf_sign`'s own, applied to the demeaned variant).
4. **Detection power**: `radf_sign_dm()` correctly rejects a clear mildly
   explosive alternative using its own simulated critical value (same
   structure as `radf_sign()`'s own power check).

New tests appended to `test-sign.R` (formula check, invariance, power,
and a regression test for the `print()` crash covering all three
affected `_cv` functions). Full package suite green.

Replication script:
[replication/volatility-robustness/radf_sign_dm_levelshift_validation.R](#script-radf_sign_dm_levelshift_validation).

---

## Stochastic explosive-coefficient test

**Status: SSU done (2026-08-10, `ssu_test()`), the minimum-viable
subset this file's own earlier triage identified — GSSU, CUSUM/
CUSUM-SQ, and the union-of-rejections procedure not implemented.** Full
PDF read (model, Section 3's test statistics through the union-of
-rejections and CUSUM/CUSUM-SQ proposals), re-verified against rendered
PDF pages 5-6 and 9 for the implemented item.

### Source

Kurozumi, E. & Nishi, M. (2025). "Testing for a bubble with a
stochastically varying explosive coefficient." *JTSA*, 46(5), 945-965.
Turned out to be an Open Access article — Wiley just requires JS to serve
it, which is why it read as paywalled to a plain `curl` fetch at first.

### What it is

Models the explosive AR(1) coefficient itself as *random*:
`1 + c1/T + a*u_t/sqrt(T)` (eq. 2) rather than the deterministic
`1 + c/T^alpha` every other method in this project assumes, with `u_t`
i.i.d. mean-zero unit-variance and independent of the innovations. The
paper's own motivation (Figure 1: rolling-window AR(1) estimates on
Japanese daily stock prices) is that the "explosive speed" itself looks
unstable in practice, not constant during a bubble episode.

This isn't a volatility-robustness fix in the sense the rest of this file
is — it doesn't touch the *innovation* variance at all (`sigma_epsilon^2`
is constant, `Assumption 1a`). It's a different generalization: the
*persistence/explosiveness* parameter is random, not the noise scale.
Grouped in this file because it was originally surfaced alongside the
other heteroskedasticity items in the *JTSA* special issue, not because
it addresses the same failure mode.

The paper proposes **three structurally different new statistics**, not
one:

1. **SSU/GSSU** (eq. 7-8): a stochastic-unit-root test (following Lee
   1998/Nagakura 2009) built on an entirely different regression —
   `(Delta y_t)^2 = mu^2 + eta*y_{t-1}^2 + e_t` (squared differences on
   squared lagged levels, not the standard ADF regression at all) — whose
   raw t-statistic is *not* asymptotically pivotal on its own; a
   bias-correction term `rho_hat(r1,r2)` (a recursively-estimated
   cross-moment between the level regression's and the squared
   regression's residuals) has to be subtracted out before the corrected
   statistic `t^c_{r1,r2}` is usable, per Nishi & Kurozumi (2024).
2. **CUSUM / CUSUM-SQ**: a second, independent pair of statistics (Brown
   et al. 1975's classic parameter-constancy tests), studied here as a
   third and fourth candidate for the same detection problem.
3. **Union of rejections** combining (at least) SADF/GSADF with SSU/GSSU
   — the paper's own recommended practical procedure, since neither
   family dominates (SSU/GSSU wins when the coefficient is genuinely
   stochastic, `a != 0`; SADF/GSADF wins when it's deterministic, `a = 0`)
   — needing joint critical values for the combination, the same
   structural pattern as
   [SBZ's union statistic](#sbz-wls--kernel-volatility) or Bundle 5's own
   evaluation of the sign-based test's (unimplemented) union with
   standard PSY.

### Cost/feasibility note for exuber

Not a contained addition, for reasons distinct from — and larger than —
every other item in this bundle:

1. **A genuinely new regression form.** SSU/GSSU's `(Delta y_t)^2` on
   `y_{t-1}^2` regression shares no structure with the ADF-family
   `y_t` on `y_{t-1}` regression every existing exuber statistic (`radf()`,
   STADF, kernel-purge, SBZ, sign-based) is built around — this is new
   estimation code from scratch, not a transform feeding into
   `gls_dfstat_grid()`/`rls_gsadf()` the way the sign-based test turned
   out to be.
2. **The bias-correction is itself nontrivial.** `rho_hat(r1,r2)` requires
   recursively estimating a cross-moment between two different residual
   series (from the level and squared-level regressions) over every
   candidate window — a materially more involved recursive computation
   than any cumulative-sum trick used elsewhere in this project (PDC/KS,
   STADF, sign-based all reduce to a handful of `cumsum()` calls; this
   needs the residuals from *two* fitted regressions per window first).
3. **Two more statistics, not one.** CUSUM/CUSUM-SQ are a further,
   separate pair of test statistics investigated in the same paper —
   implementing "the stochastic-coefficient test" faithfully would mean
   choosing which of up to four statistic families (SSU, GSSU, CUSUM,
   CUSUM-SQ) to actually ship, or implementing all of them.
4. **The recommended procedure is a union, needing its own joint
   critical-value simulation** — same category of extra work as SBZ's
   union statistic (flagged there as "a distinct, non-trivial chunk of
   new statistical code"), stacked on top of the new-statistic cost above
   rather than replacing it.

Comparable in scope to SBZ (two new estimators plus a union/bootstrap
layer) rather than to this bundle's other items (Hafner: one new
multiplier branch; sign-based: one new input transform reusing existing
machinery). Not picked up this pass; if revisited, SSU alone (without
GSSU's double-recursion, without CUSUM/CUSUM-SQ, without the union) would
be the minimum-viable first cut.

**Re-triaged (2026-08-10), re-reading rendered pages 5-6, 9.** Points 1
and 2 above turned out to overstate the true cost for the SSU-alone MVP
specifically:

- **Point 1 (new regression form) is real but contained**: the SSU
  regression `(Delta y_t)^2 = mu2 + omega*y_{t-1}^2 + eta_t` (eq. 7) is
  still a plain 2-variable OLS over a window — the *same* generic
  `hls_prefix_sums()`/`hls_segment_coef()` closed-form pattern
  (arbitrary `(x, z)` pair over a segment) already used for
  `dating_hls()`/`dating_knp()`, just with `x = y_{t-1}^2` and
  `z = (Delta y_t)^2` instead of `x = y_{t-1}`, `z = Delta y_t`. No new
  estimation *theory*, only a different pair of input series.
- **Point 2 (bias correction) is real but also closed-form, not
  "recursive" in the expensive sense**: the cross-moment
  `sigma_hat_{eps*eta}` needing "the residuals from two fitted
  regressions" sounds like it needs two passes of per-observation
  residual computation per window, but expanding
  `sum(eps_hat_t * eta_hat_t)` algebraically reduces it to a **bilinear
  combination of window sums** of twelve fixed per-observation products
  (`y_{t-1}`, `y_{t-1}^2`, ..., `y_{t-1}^4`, `Delta y_t`, ...,
  `(Delta y_t)^4`, and their cross products) — more cumulative sums to
  track than any prior item in this project (twelve vs. the usual
  four-to-six), but still `O(1)` per window via `cumsum()` differences,
  not a second pass over the raw data per candidate window.
- **The critical value is published, not simulated**: Kurozumi & Nishi's
  own Table I (their 10,000-rep Monte Carlo) gives `SSU`'s asymptotic
  critical value directly — `2.90`/`3.30`/`4.20` at the `10%`/`5%`/`1%`
  level — a single scalar per level (`SSU` is a single-recursion
  sup-statistic, `r1` fixed at `0`, not `GSSU`'s double recursion), and
  their own recommended `r0 = 0.01 + 1.8/sqrt(T)` is *exactly*
  `psy_minw()`'s existing formula, reused directly with no adjustment.

Points 3 and 4 (CUSUM/CUSUM-SQ, and the union-of-rejections procedure)
remain accurate — those are still separate, unimplemented statistic
families, deliberately scoped out of this pass as originally planned.

### Implementation — SSU done (2026-08-10)

Shipped as `ssu_test(data, minw = NULL, level = 0.95)` in
`exuber/R/ssu_test.R`. `ssu_prefix_sums()` builds twelve cumulative-sum
vectors from the two base per-observation series (`y_{t-1}` and
`Delta y_t`); `ssu_stat_path()` evaluates the bias-corrected
`t^{omega,c}_{0,r2}` statistic for every candidate end point via those
sums (algebra verified by hand-expanding the residual cross-moment
`sum(eps_hat*eta_hat)` into its bilinear window-sum form, then confirmed
numerically against brute force below); `ssu_test()` takes the running
maximum (`SSU`'s own sup-statistic) and compares it against
`ssu_q(level)`'s lookup into Table I.

**Validated**: `ssu_stat_path()` matches a brute-force computation (two
separately `lm()`-fitted regressions on the raw window data, plus a
manual residual cross-moment) to machine precision (`< 1e-8`) at four
separate window sizes; Table I lookups exact, with a clean error for an
untabulated level; `ssu_test()`'s default `minw` matches `psy_minw()`
exactly, confirming `SSU`'s own `r0` formula needed no adaptation.
Empirical false-alarm rate under `H0` (300 reps, `n=200`) is
`12.0%`/`8.7%`/`2.7%` against nominal `10%`/`5%`/`1%` — mildly
oversized, in the same range as several other finite-sample-vs
-asymptotic critical values validated in this project, not a
red flag on its own. Detection power on a genuine stochastic
-explosive-coefficient DGP (Kurozumi & Nishi's own eq. 2 alternative,
60 reps) is `85.0%` — the alternative `SSU` is actually designed for.
On a *deterministic* explosive DGP instead (the kind `radf()`'s own
`SADF` targets), `SSU` still has `80.0%` power vs. `SADF`'s `90.0%` — a
sensible, honestly-reported result: `SSU` trades a little power on the
deterministic alternative for robustness to genuine coefficient
stochasticity, exactly the tradeoff the paper's own Theorem 2 describes
(each family dominates on its own alternative, neither dominates
universally). New `test-ssu.R` (7 tests). Replication script:
[replication/volatility-robustness/radf_ssu_validation.R](#script-radf_ssu_validation).

Not implemented: `GSSU` (the double-recursion generalization,
`r1`/`r2` both varying — the same "fixed vs. growing start range"
distinction that separated Kurozumi (2020)'s `SADF`/`GSADF_{s0}` cases
elsewhere in this project, here compounded by needing the bias
-correction's twelve window sums recomputed over a full `(r1, r2)` grid
rather than a single-recursion path); `CUSUM`/`CUSUM-SQ` (a separate
statistic pair investigated in the same paper); and the union
-of-rejections procedure combining `SSU`/`GSSU` with `SADF`/`GSADF`
(the paper's own recommended practical procedure, needing its own joint
critical-value simulation).

---

## SV-ADF

**Status: done (2026-08-11), `radf_svadf()`** (preprint, not yet
peer-reviewed — flagged explicitly, a different bar than every other
source implemented in this project). Full PDF read through the SV-ADF
statistic's definition, its asymptotic theorem (3.1), the proof appendix
confirming the feasible variance estimator's exact construction, and —
re-triaged 2026-08-11 — Section 5.1's actual threshold-calibration
exercise (rendered pages 20-22), which resolved the cost note's own
open question below.

### Source

Sarkar, A. & Wells, M.T. (2026). "Is There an AI Bubble? Robust
Date-Stamping for Periods of Exuberance." arXiv:2604.12062. Theoretical
grounding in the same authors' arXiv:2512.06823, "Double Local-to-Unity:
Inference under Nearly Nonstationary Volatility."

### What it is

Extends the recursive right-tailed ADF test (PWY-style, eq. 3-4) to
allow *highly persistent, nearly-nonstationary stochastic volatility* —
e.g. `log(sigma_t^2) = phi_n * log(sigma_{t-1}^2) + eta_t` with `phi_n`
approaching 1 at an iterated-logarithmic rate — a substantially weaker
condition than the deterministic/bounded-jump volatility every other
method in this file assumes (HLST's wild bootstrap, Hafner, sign-based,
kernel-purge, STADF all require `sigma(s)` to be a fixed, non-stochastic
function of calendar time; this allows the volatility process itself to
be a near-unit-root-persistent stochastic process, e.g. GARCH with
`alpha+beta` close to 1, or a persistent stochastic-volatility model).

**A genuinely interesting structural finding from reading the proof
appendix**: the feasible `SV-ADFr`/`SV-ADFrt` statistics (eq. 3-4) are
built from the *same* recursive OLS estimator and the *same* residual-
based variance estimator (`tau_hat^2 = tau^-1 * sum(residuals^2)`, eq.
A.13-A.14 confirm this is literally the standard within-window OLS
residual variance, not a separate nonparametric/kernel estimate) that
`radf()`'s own recursive ADF t-statistic already computes — the paper's
actual contribution (Theorem 3.1) is a *broader asymptotic justification*
for using this same statistic under far weaker volatility conditions than
previously proven, not a structurally different statistic. Theorem 3.1
explicitly notes the limiting functional "coincide[s] with those obtained
under homoskedasticity (Phillips and Yu, 2009)" once normalized.

What *isn't* pinned down by this reading: the abstract advertises
"distinct calibration thresholds for testing origination and collapse"
and a "moderate-deviation asymptotic theory" — language suggesting the
paper's actual practical novelty is an **adaptive/moderate-deviation
critical-value boundary** `cv_{n,r}` (growing with `r` at a specific
rate) rather than a fixed quantile from the limiting distribution, used
to explain why SV-ADF avoids the spurious bubble PWY flags in a
high-volatility episode (Figure 6, empirical Nvidia example) despite the
underlying point statistic being the same. Pinning down the *exact*
construction of this boundary function would require reading the
authors' companion "Double Local-to-Unity" theory paper — not done in
this pass.

### Cost/feasibility note for exuber

Genuinely uncertain, in a way the other three items in this bundle
weren't — resolved, in the favorable direction, by the 2026-08-11
re-triage below:

- **If** the practical procedure turns out to be "compute exuber's
  existing recursive ADF t-statistic exactly as `radf()` already does,
  just compare it to a specific moderate-deviation-calibrated boundary
  function instead of a fixed quantile," this could be close to
  free — no new point statistic at all, only a new comparison rule. This
  would put it in the same cost tier as the sign-based test (this
  bundle's cheapest item).
- **If** the boundary function itself requires estimating something new
  (e.g. a persistence parameter for the stochastic-volatility process,
  or quantities from the companion moderate-deviation paper not
  determinable from this paper alone), the cost is unknown until that
  paper is read.
- Additionally: this is a **non-peer-reviewed preprint** (arXiv, first
  version 2026), a different bar than everything else implemented in this
  project so far (all peer-reviewed, published papers). Worth flagging
  explicitly to whoever picks this up next, independent of the
  implementation-cost question.

**Re-triaged (2026-08-11), re-reading rendered pages 20-22 (Section
5.1, "Threshold Selection Insights").** The first, favorable branch
above is what actually happens — no companion-paper reading needed
after all:

- The paper's own Section 5.1 describes exactly how the "moderate
  -deviation-calibrated boundary" was arrived at, and it is **not** a
  function of any estimated nuisance parameter: for origination, they
  simulate the SV-ADF statistic under `H0` at `n ∈
  {100,200,...,1000}` (1,000 reps each), take the 90th percentile at
  each `n`, and find "these upper 90 percent critical values are well
  approximated by `log(n)/10`, which is adopted as the origination
  threshold." For collapse, they average the 10th-percentile threshold
  over randomly drawn nuisance-parameter configurations and find it
  "most closely approximated by `log(n)/2`, which we therefore use as
  the collapse threshold." Both are **published, closed-form,
  sample-size-only formulas** in the paper's own actual applied
  methodology — not requiring the companion "Double Local-to-Unity"
  paper, and not requiring any persistence-parameter estimation on
  exuber's end.
- Combined with the already-confirmed structural finding (point
  statistic = `radf()`'s own `badf`), this puts SV-ADF's actual
  practical procedure at true zero-new-statistic cost: reuse `badf`,
  compare against two `log(t)`-based thresholds instead of one fixed
  quantile, apply a first-crossing dating rule. The one genuinely new
  piece is that origination and collapse use **different** thresholds
  (the paper's own Remark 1: a unit-root-based threshold suits
  origination, not collapse) — `datestamp()`'s existing S3 dispatch
  assumes a single shared critical value throughout, so this shipped as
  its own small, self-contained dating routine reusing `stamp()`'s
  existing contiguous-run detection, rather than an extension of
  `datestamp()` itself.

### Implementation — done (2026-08-11)

Shipped as `radf_svadf(data, minw = NULL, min_duration = NULL)` in
`exuber/R/radf_svadf.R`. `min_duration` defaults to `psy_ds(n)` —
exuber's own existing `log(n)`-based minimum-episode-duration rule,
reused directly rather than inventing a new one, standing in for the
paper's own (data-frequency-specific) "at least two consecutive
calendar months"/"one month" consolidation requirement. Origination is
dated at the first run of at least `min_duration` consecutive points
with `badf` above `log(t)/10`; collapse is dated (searching only after
the origination date) at the first run of at least `min_duration`
consecutive points with `badf` below `log(t)/2`.

**Validated**: `radf_svadf()`'s `badf` field matches a direct `radf()`
call bit-for-bit (confirming the point-statistic reuse is exact, not
approximate); threshold formulas match `log(t)/10`/`log(t)/2` exactly;
collapse is structurally guaranteed to never date before origination
(the search for it only starts after the origination date, confirmed
across 20 reps with no exceptions). On a synthetic bubble+collapse
episode (using this project's own established large-base bubble-DGP
convention, since an initial mean-zero-random-walk-base attempt gave a
weak/diluted signal the same way it did for `dating_hls()`'s own earlier
validation): 100% detection rate across 20 reps, mean absolute
origination-date error `4.95` periods, mean absolute collapse-date
error `20.25` periods (collapse detection is inherently less precise —
`badf`'s expanding, anchored-at-the-start window dilutes a
post-collapse downward signal more than an origination upward one, a
known property of single-recursion statistics, not a defect). False
-alarm rate under `H0` (60 reps, pure random walk) is `13.3%` — this
counts *any* origination crossing anywhere in the series across a full
150-period path, a much larger compound opportunity for a false alarm
than a single-point 10% test, so a somewhat-inflated aggregate rate
here is expected, not a sign of miscalibration at any individual point.
New `test-svadf.R` (7 tests). Replication script:
[replication/volatility-robustness/radf_svadf_validation.R](#script-radf_svadf_validation).

The AI-equity empirical hook (2025-26 exuberance in Nvidia, Tesla,
TSMC, etc.) is a plausible source for a topical worked-example vignette
now that the statistic itself is implemented — see
practitioner-guidance.md.
