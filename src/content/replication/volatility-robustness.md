---
title: "Volatility-robustness tests"
blurb: "Tests robust to time-varying innovation variance: time-transformed, kernel-purged, WLS, sign-based, and stochastic-coefficient routes."
order: 1
---
This file covers right-tailed unit-root tests that are modified to stay
correctly sized when the innovation variance changes over time, whether
deterministically or stochastically. The original GSADF of PWY and PSY
assumes homoskedasticity, and every method here is a different fix for
that same failure. The status legend is as follows. `done` means the method
is implemented and tested against a published number. `partial` means it is
implemented but not yet cross-checked. `evaluated` means we have read the
source and have not implemented the method. `todo` means we have not
started, and `blocked` means we could not get access to the source.

| Method | Paper | Status |
|---|---|---|
| [Time-transformed test (STADF/GSTADF)](#time-transformed-test-stadf--gstadf) | Kurozumi, Skrobotov & Tsarev (2024) | **done** |
| [Kernel-purge test](#kernel-purge-test) | Harvey, Leybourne, Taylor & Zu (2024/2025) | **done** (with-intercept variant) |
| [SBZ (WLS + kernel volatility)](#sbz-wls--kernel-volatility) | Harvey, Leybourne & Zu (2019) | **done**, bug found and fixed in validation, see note |
| [Sieve bootstrap (autocorrelated innovations)](#pedersen--schütte-sieve-bootstrap) | Pedersen & Montes Schütte (2020) | **done** |
| [Skewness-corrected wild bootstrap](#hafner-skewness-corrected-wild-bootstrap) | Hafner (2020) | **done** |
| [Sign-based sGSADF](#sign-based-sgsadf) | Harvey, Leybourne & Zu (2020); level-shift robustness and demeaned variant: Harvey, Leybourne, Tatlow & Zu (2025) | **done** |
| [Stochastic explosive-coefficient test](#stochastic-explosive-coefficient-test) | Kurozumi & Nishi (2025) | **done** (SSU 2026-08-10; GSSU, UR/GUR union and CS/GCS/CSSQ/GCSSQ 2026-09-29, `ssu_test(type =, union =)`, `cusum_test()`) |
| [SV-ADF](#sv-adf) | Sarkar & Wells (2026, preprint) | **done** (2026-08-11, `datestamp(option = "svadf")`; preprint, not peer-reviewed) |

All papers: [references.md](/replication/references#volatility-robustness).

---

## Time-transformed test (STADF / GSTADF)

Status: done. The test is implemented in `exuber/R/radf_tt.R`
(`radf_tt()`, `radf_tt_cv()`) and cross-checked in
`exuber/tests/testthat/test-tt.R`.

### Source

Kurozumi, E., Skrobotov, A., & Tsarev, A. Time-Transformed Test for Bubbles
under Non-stationary Volatility. Journal of Financial Econometrics (2024),
`doi:10.1093/jjfinec/nbae026`. The working paper is arXiv:2012.13937, which
is freely available and is the version we used. Equation and theorem
numbers below refer to arXiv v2, 15 Nov 2021.

We downloaded the paper and converted it with `pdftotext`. Subscripts,
hats and bars are lossy in that conversion, so we rendered the source pages
of every ambiguous formula to PNG (`pymupdf`) and read the typeset math
directly and did not trust the OCR text. The formulas we implemented are
eq. (9), (17)-(19) and Theorems 1 and 2.

### Why this item

Theorem 1 proves that the null limiting distribution coincides with the
homoskedastic GLS-demeaned SADF/GSADF distribution (Whitehouse, 2019). No
bootstrap is needed, unlike HLST's wild-bootstrap PSY or SBZ. That makes
this the cheapest of the volatility-robust tests to compute and to verify.
The paper also reports a fixed, literal triple of critical values, which is
rare in this literature, because most of it relies on per-dataset bootstrap
p-values that cannot be reproduced bit for bit.

### Exact numbers reproduced

Footnote 4 (arXiv v2, page 9): *"For r0 = 0.1, they are equal to 2.319,
2.626, 3.223 for 10%, 5% and 1% significance levels."*

We corrected one point during implementation. This triple belongs to
STADF, the single-sup, PWY/SADF-style statistic `sup_r2 ADF^r2_0`, and not
to GSTADF, the double-sup, PSY/GSADF-style statistic. Whitehouse (2019),
the source of these critical values, is itself about the GLS-demeaned
SADF-style test. The paper says so in its Section 1: "Whitehouse (2019)
considered the Phillips et al. (2011) test with the GLS-type detrending,"
and PWY 2011 is the single-sup test. We confirmed this by simulation. The
STADF statistic (our `gls_dfstat_grid()`, n=300, minw=30, 1500 Monte Carlo
replications, seed 20260808) gives (2.407, 2.702, 3.533), which is close to
(2.319, 2.626, 3.223) once finite-T and Monte Carlo noise are allowed for.
The GSTADF (double-sup) statistic under the same setup gives (3.157, 3.436,
4.302), which is not close. We verified the formula independently (see
below), so the problem was in identifying the target and not an
implementation bug. It is worth recording because it is an easy trap. The
paper's Theorem 1 says in one sentence that both STADF and GSTADF
"coincide" with the Whitehouse (2019) critical values, which reads as if
one table covered both, but only the SADF-style one is published.
GSTADF(r0) critical values can be simulated with `radf_tt_cv()`. We have no
published anchor value for them in this pass. The paper says they are
"easily computed from the R-code available in
https://sites.google.com/site/antonskrobotov/", and we did not fetch that
code.

### What's actually pivotal vs. what needs estimation

- The null limiting distribution is pivotal, meaning it does not depend on
  the volatility path, once the series has been time-transformed with the
  true variance profile (Theorem 1). It remains pivotal when we use the
  estimated profile (Theorem 2). This is why `radf_tt_cv()` can simulate
  critical values once, with `y <- cumsum(rnorm(n))` (no volatility, no
  bootstrap), instead of once per dataset. It follows the pattern of
  `radf_mc_cv()` but uses a different statistic (no intercept,
  GLS-demeaned).
- The feasible statistic on real, volatile data still needs the variance
  profile to be estimated from the data (eq. 18-19). That involves a
  kernel-weighted local no-intercept regression of Δy̌_t on y̌_{t-1},
  truncated residuals, a cumulative-sum-of-squares profile and a
  generalized-inverse time transformation. It is implemented in
  `variance_profile()`.

### Formulas implemented (verified against rendered PDF pages, not OCR text)

- `y̌_t := y_t − y_0` is the GLS-demeaning. We subtract the first
  observation and fit no intercept. This differs from what exuber's
  `radf()` does, because `radf()` always OLS-demeans with a fitted
  intercept (see "why not exubercore" below).
- Eq. (9): `ADF^{r2}_{r1} := Σ y̌_{t-1}Δy̌_t / sqrt(σ̂²(r1,r2) Σ y̌²_{t-1})`
  is the no-intercept recursive Dickey-Fuller t-statistic. It is the
  t-statistic on β in `Δy̌_t = β y̌_{t-1} + e_t`, with no constant. We
  implemented it as `gls_dfstat_grid()`, which is fully vectorized over the
  (r1, r2) grid through cumulative sums, so no per-window regression loop
  is needed. It is in the same complexity class as PSY's GSADF grid and
  needs no bootstrap.
- Eq. (18): a local (Nadaraya-Watson-type) kernel estimate of the
  time-varying AR(1) coefficient δ̂_t. The default is a uniform kernel, as
  in the paper's Monte Carlo.
- Eq. (19): the variance profile η̂(s) is a normalized cumulative-sum-of-
  squares step function of the truncated local-regression residuals. We
  invert it by linear interpolation (exact, since η̂ is piecewise linear on
  the observation grid) to obtain ĝ(s), which we use to resample and
  time-transform the series.
- Footnote 6: the truncation threshold is ψ_T = c̄·T^(1/7), where c̄ is the
  maximum local-window residual SD over rolling windows of 10% of the
  sample.

### Deliberate simplifications vs. the paper (cost/benefit)

- Bandwidth. The paper uses leave-one-out cross-validation over h ∈
  [T^-0.5, T^-0.3]. Our default is the fixed plug-in `h = T^(-2/5)`, the
  midpoint of that range on the log scale, and `h` is exposed as a user
  parameter. Full cross-validation would rerun the O(T · Th) kernel fit
  about ten times for a bandwidth grid. We left it out to keep this pass
  light. If empirical size control matters in practice, it is the first
  thing to add.
  `# ponytail: fixed plug-in bandwidth, add CV grid search if empirical size is off`
- Lag augmentation. The paper's Monte Carlo uses p=0 (no augmentation
  lags) throughout, and `radf_tt()` likewise supports only p=0. `radf()`'s
  `lag` parameter has no equivalent here yet.
- Not wired into `tidy()`, `autoplot()` or `datestamp()`. At first,
  `radf_tt_obj` had only a `print` method. Full S3 method parity with
  `radf_obj` (tidiers, plotting and date-stamping) was follow-up work,
  which we completed later (see the 2026-08-18 note below).

### Why not exubercore (C++)

The paper's own selling point is that this test is cheap. It needs no
bootstrap, and its recursive statistic (eq. 9/17) is a closed-form ratio of
cumulative sums and not a per-window OLS fit that needs a matrix
inversion. `gls_dfstat_grid()` can therefore be fully vectorized in R with
`outer()` on prefix sums, and we do not need exubercore's C++ `radf()`. The
C++ `radf()` is also hard-coded to fit an intercept whenever `lag == 0`
(see `exubercore/src/radf.cpp`), so it could not be reused unmodified for
the no-intercept GLS-demeaned case at lag 0 in any event. We should port
it to exubercore only if profiling shows the O(T²) grid becoming a
bottleneck for T in the thousands. We did not see that in this pass's
testing, which went up to T of about 2000.

### Independent validation (2026-08-09)

We reran everything from scratch with fresh seeds and sample sizes, so
this is an independent check and not a rerun of the same numbers in
`test-tt.R`.

1. Formula check. We compared `gls_dfstat_grid()`, the vectorized
   closed-form implementation, with a brute-force per-window `lm(dy ~ ylag
   - 1)` fit on data and parameters that the test suite does not use
   (`n=65`, `minw=15`, `seed=4242`, whereas the suite uses `n=40`,
   `minw=10`, `seed=7`):

   ```r
   set.seed(4242)
   y <- cumsum(rnorm(65)); minw <- 15
   res <- exuber:::gls_dfstat_grid(y, minw)
   # ...vs. per-window lm(dy[1:b] ~ ylag[1:b] - 1) for b in minw:n1
   ```

   The result is `max|badf_formula - badf_lm| = 6.66e-16`, which is
   machine precision. The closed-form cumulative-sum implementation is
   exact and not an approximation.

2. Critical-value replication. We ran an independent Monte Carlo against
   Whitehouse (2019)'s published STADF triple (quoted in footnote 4 of the
   paper, `r0=0.1`: 2.319, 2.626, 3.223), with a fresh seed and more
   replications than the package's own test (`nrep=4000`, against the
   suite's `1500`):

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

   The deviation has the same direction and roughly the same size as in
   the package's own test, which uses `tolerance = 0.15` and passes. This
   is consistent with Monte Carlo noise around a T→∞ asymptotic target,
   and it does not point to a formula error. The exported
   `radf_tt_cv(n=300, minw=30, nrep=4000, seed=555)` gives (2.399, 2.737,
   3.346) on yet another seed, in the same range.

3. Behavioral check. We wanted to know whether the test fires on an
   explosive series. We built a synthetic series that follows a unit root,
   then turns explosive (ρ=1.04) and then collapses (150 observations),
   and ran `radf_tt()` on it:

   ```r
   gsadf = 5.781   # vs. the ~2.0-3.3 critical-value range above
   ```

   The statistic exceeds every critical value in the table above with
   room to spare, as it should.

4. Full existing suite. `test-tt.R`: 8 passed, 0 failed (1 skipped,
   CRAN-only).

Conclusion: the formula is exact (item 1), the critical-value simulation
reproduces the published target within the expected Monte Carlo noise on
an independent run (item 2), and the test behaves correctly on a case it
should detect (item 3). We found no issues.

Replication script:
[replication/volatility-robustness/radf_tt_validation.R](#script-radf_tt_validation).

### `datestamp()`/`autoplot()` parity (done, 2026-08-18)

The "Open follow-ups" below used to list full tidy, autoplot and datestamp
S3 method parity with `radf_obj`, and that item is now closed.
`radf_tt_cv()` computed only the three scalar critical values
(`adf_cv`/`sadf_cv`/`gsadf_cv`) that `summary()` and `tidy()` need. It
never computed the time-varying boundaries `badf_cv` and `bsadf_cv` that
`datestamp()` and `autoplot()` need. The gap surfaced when a user reported
the equivalent gap in `radf_sign()`, which led us to check all three
functions of the GLS-demeaned family (`radf_sign()`, `radf_sign_dm()` and
`radf_tt()`). The fix turned out to be small. `gls_dfstat_grid()`, the
shared no-intercept recursive-DF machinery that `radf_tt()` already calls,
returns the sup-over-all-window-starts `bsadf` at each point directly. It
does not need the cummax-of-single-start-series shortcut that
`radf_mc_cv()`'s `bsadf_cv` uses to adapt the output shape of the base C++
engine. We needed no new recursion and no new shortcut derivation. We only
had to collect `badf` and `bsadf` across replicates, which were already
computed and then discarded, and take the per-time-point quantile. This is
the same construction that `radf_mc_cv()` uses.

Validation:

1. Formula-exact check. The last row of `badf_cv` is bit-identical to
   `adf_cv`. This is an identity and not an approximation, because `adf` is
   by definition the last point of `badf` in each replicate, so their
   quantiles across replicates must match exactly.
2. Monte Carlo size. The false-alarm rate under `H0` (pure random walk,
   n=100, minw=20, nrep=2000, 300 replications) is 3.3% against a nominal
   5%, for both `option = "gsadf"` (`bsadf_cv`) and `option = "sadf"`
   (`badf_cv`). The test is conservative and shows no size distortion.
3. Power, cross-checked against the established baseline. On an identical
   synthetic bubble (60 pre-bubble observations, 30 explosive observations
   at ρ=1.03 and 10 post-collapse observations, 50 replications), the new
   `datestamp()` for `radf_tt()` detects the bubble 18% of the time, while
   the long-established `radf()`/`radf_mc_cv()` pipeline detects it 16% of
   the time. The low absolute number reflects a DGP with a short, moderate
   bubble at this sample size. We know that because the already trusted
   baseline gives the same number on the identical series, and so the new
   critical-value code is not at fault.
4. On `sim_data`'s own bundled bubble series (`psy2`), `datestamp()` finds
   two sensible episodes (Start=21/Peak=27/End=35 and
   Start=55/Peak=55/End=73).

`radf_sign_cv()` and `radf_sign_dm_cv()` had the identical gap. We checked
it immediately afterwards and found that the same shortcut applies
unchanged. Both call `gls_dfstat_grid()` as well, on a sign-transformed
series, so the same fix went in during the same pass. See the note on
`datestamp()`/`autoplot()` parity in the [Sign-based sGSADF](#sign-based-sgsadf)
section for that validation record.

### Open follow-ups

- Extend `radf_tt_cv()` and `radf_tt()`'s asymptotic critical values across
  a grid of `r0`, not only r0=0.1. The paper says these are directly
  available from their R code at https://sites.google.com/site/antonskrobotov/,
  which we did not fetch in this pass.
- ~~`radf_sign_cv()` and `radf_sign_dm_cv()` still lack `badf_cv` and
  `bsadf_cv`~~. This was done on 2026-08-18, in the same pass. See the
  `datestamp()`/`autoplot()` parity note in the
  [Sign-based sGSADF](#sign-based-sgsadf) section below.
- The [sign-based test](#sign-based-sgsadf) was blocked on paywall access
  when this item was written. This paper's literature review (Section 1)
  is a decent secondary source for its qualitative behaviour, but not for
  exact numbers.

---

## Kernel-purge test

Status: done (with-intercept variant). The test is implemented in
`exuber/R/radf_kp.R` (`radf_kp()`) and cross-checked in
`exuber/tests/testthat/test-kp.R`.

### Source

Harvey, D. I., Leybourne, S. J., Taylor, A. M. R., & Zu, Y. (2024). A new
heteroskedasticity-robust test for explosive bubbles. Journal of Time
Series Analysis. `doi:10.1111/jtsa.12784`. The paper is open access
(CC-BY), and we downloaded it directly from the University of Essex
repository (`repository.essex.ac.uk`) without any paywall problem.
We checked the formulas and Table I below by rendering the PDF pages to
images (PyMuPDF) and reading the typeset math. The first `pdftotext` pass
scrambled the columns of Table I, as it did for other papers in this file,
so we did not rely on it.

### What it is

Instead of resampling (`radf_wb_cv`) or deforming time (`radf_tt`), this
test "purges" unconditional heteroskedasticity directly. It estimates the
spot volatility \eqn{\hat\sigma_t} with a kernel (eq. 4, Gaussian kernel),
divides each first difference by it, and cumulates the result (eq. 5) to
obtain a volatility-standardized series \eqn{x_t}. The unmodified PSY/GSADF
test is then run on \eqn{x_t} in place of the raw series.

The key result (Theorem 1 and Remark 3.2) is that the null limiting
distribution of the purged statistic is identical to the standard
homoskedastic GSADF null. That makes the test the cheapest of the
volatility-robust tests to implement and to use in exuber. `radf_kp()` is
a thin wrapper. It reuses `kernel_spot_vol()`, which we had already built
for [SBZ](#sbz-wls--kernel-volatility), for the volatility estimate, and
then calls the existing, unmodified `radf()` on the purged series. The
Monte Carlo critical values that exuber already has (`radf_mc_cv()`) apply
directly, so no new critical-value simulation is needed. Every downstream
method (`tidy()`, `autoplot()`, `datestamp()` and so on) works on the
result without further effort.

The paper also proposes a without-intercept variant (\eqn{PSY^*_\sigma})
and a union-of-rejections test that combines both (\eqn{UPSY_\sigma}). See
"Not implemented" below.

### Exact numbers reproduced

Table I gives asymptotic and finite-sample critical values (minimum window
ε = 0.1, Gaussian kernel, h = 0.1·T^-0.25, 2000 Monte Carlo replications).
We verified it against the rendered images.

| T | PSY_σ (10%/5%/1%) | PSY*_σ (10%/5%/1%) | UPSY_σ (10%/5%/1%) |
|---|---|---|---|
| 100 | 1.629 / 1.828 / 2.392 | 3.637 / 4.158 / 5.553 | 3.950 / 4.527 / 6.129 |
| 200 | 1.608 / 1.789 / 2.140 | 3.226 / 3.595 / 4.330 | 3.468 / 3.804 / 4.589 |
| 400 | 1.712 / 1.935 / 2.296 | 3.167 / 3.446 / 4.007 | 3.361 / 3.598 / 4.145 |
| ∞ | 1.875 / 2.094 / 2.486 | 2.978 / 3.296 / 3.859 | 3.186 / 3.486 / 3.951 |

The PSY_σ (with-intercept) column is the one we implemented and tested.
By Remark 3.2, its T = ∞ row should equal PSY (2015)'s published
asymptotic GSADF critical values at the same minimum window, and so it
should also match what exuber's `radf_mc_cv()` already produces. We test
this directly. `test-kp.R` simulates the null GSADF distribution of
`radf_kp()` at n = 300 and checks it against (a) an independently computed
`radf_mc_cv(300)` and (b) the published T = 400 row above. The tolerance is
wide enough to cover the extra finite-sample noise that the kernel
volatility estimation step adds on top of ordinary Monte Carlo error. In
one run, `radf_kp` gave (1.685, 1.897, 2.337) at n = 300, compared with
`radf_mc_cv(300)`'s (1.873, 2.121, 2.429) and the published T = 400 row's
(1.712, 1.935, 2.296). The numbers are in the same range, with `radf_kp`
a little lower. A plausible reason is that dividing by an estimated
volatility path, and not the true one, attenuates the effective variance a
little in finite samples. This check is looser than the formula-level
brute-force verification of STADF or SBZ (see below for why it cannot be
tighter). The qualitative pattern matches the theorem, though: all three
sets of numbers lie within about 0.3 of each other, and PSY_σ is clearly
lower than PSY*_σ and UPSY_σ, as the table's own pattern shows.

The formula of the core statistic (eq. 4-7) is the same Dickey-Fuller
regression that `radf()` already computes, applied to transformed data. It
did not need an independent brute-force check of the kind that the new
closed-form statistics of STADF and SBZ required. It calls unmodified,
already tested exuber code (`radf()` and `rls_gsadf`) on a new input
series, so what we actually verify is the transform (`kernel_purge()`).
The correlation and quantile-matching tests do that.

### Independent validation (2026-08-09)

This run uses a fresh seed and a larger sample and replication count than
`test-kp.R`, which uses `n=300, nrep=500, seed=2`. We used `n=400`, which
matches the published Table I row exactly and so avoids interpolation,
with `nrep=800` and `seed=31415`.

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

This is tighter than the STADF replication above (maximum difference 0.04
against STADF's 0.21), and tighter than what the package's own
`test-kp.R` observed at the smaller `n=300`. It agrees with the theorem's
claim that the statistic converges as `T` grows toward the table's design
point of `T=400`. The result therefore behaves as Remark 3.2 predicts and
does not look like an accident of one seed.

Full existing suite: `test-kp.R`, see the combined results at the end of
this file's validation runs.

Replication script:
[replication/volatility-robustness/radf_kp_validation.R](#script-radf_kp_validation).

### Not implemented: without-intercept variant and union test

`PSY*_σ` (the without-intercept variant) needs a no-intercept DF regression
on non-demeaned data. That is subtly different from the
`gls_dfstat_grid()` of [STADF](#time-transformed-test-stadf--gstadf),
which forces GLS-demeaning by subtracting the first observation. The `x_t`
of this paper does not need that, because the intercept is already
redundant by construction, and not because of an explicit demeaning step.
Reusing `gls_dfstat_grid()` as it stands would be asymptotically
equivalent to the paper's literal (non-demeaned) formula but not
bit-faithful to it. The union-scaling logic of `UPSY_σ` has the same
structure as SBZ's union procedure, which is already implemented in
`radf_sbz_cv()`. Both variants are follow-up work, and not a large piece
of it. Most of the pieces (the no-intercept grid statistic and the
union-of-rejections scaling) already exist elsewhere in the codebase and
mainly need assembling, without any new theory.

---

## SBZ (WLS + kernel volatility)

Status: done. The test is implemented in `exuber/R/radf_sbz.R`
(`radf_sbz_cv()`) and cross-checked in `exuber/tests/testthat/test-sbz.R`.
Independent validation on 2026-08-09 (below) found a real off-by-one bug
that badly oversized the `supDF` and `U` statistics. We fixed it and
revalidated, and the section on that validation tells the full story. We
checked all formulas (eq. 5, 6 and the union statistic) by rendering the
source PDF pages to images and reading the typeset math, and not only from
the first `pdftotext` pass.

### Source

Harvey, D.I., Leybourne, S.J. & Zu, Y. (2019). Testing explosive bubbles
with time-varying volatility. Econometric Reviews, 38(10), 1131-1151. We
read the working paper, Granger Centre Discussion Paper 18/05, University
of Nottingham (open access,
`nottingham.ac.uk/research/groups/grangercentre/documents/18-05.pdf`). We
downloaded it and converted it with `pdftotext`. We then re-verified the
Table 1 numbers below by extracting word coordinates directly from the PDF
with PyMuPDF (`fitz`), because the plain `pdftotext -layout` output
scrambled the row and column order of the table (on one pass a label was
merged into the header row). The coordinate extraction confirmed the
correct row and column assignment.

### What it is

SBZ is a weighted-least-squares (GLS-type) variant of the PWY/PSY sup-ADF
test. It uses a nonparametric kernel estimate of the unknown, time-varying
volatility as weights (eq. 6 in the paper: Gaussian kernel,
leave-one-out cross-validated bandwidth). The null distribution still
depends on the volatility path, so size control requires a wild bootstrap
(the algorithm of Harvey, Leybourne, Sollis & Taylor, applied jointly to
supBZ and supDF). Beyond the WLS statistic itself, the paper contributes a
union-of-rejections rule that combines supDF (the classic PWY/PSY test) and
supBZ, with a single scaling constant that makes the union asymptotically
correctly sized. We reject if `U := max(supDF, (qDF/qBZ)*supBZ) > qDF`.

### Exact numbers reproduced

Table 1 shows the empirical application (bootstrap p-values, M=499
bootstrap replications, FTSE Dec 1985-Dec 1999 and S&P 500 Jan 1980-Mar
2000). We verified it by coordinates.

| Series | supDF | supBZ | U |
|---|---|---|---|
| FTSE Daily | 0.288 | 0.016 | 0.046 |
| FTSE Weekly | 0.275 | 0.146 | 0.201 |
| FTSE Monthly | 0.477 | 0.279 | 0.315 |
| SP500 Daily | 0.267 | 0.000 | 0.003 |
| SP500 Weekly | 0.170 | 0.002 | 0.011 |
| SP500 Monthly | 0.153 | 0.044 | 0.071 |

The numbers agree with the paper's own narrative ("supDF does not reject
for any series"; "supBZ rejects at the 0.05-level for daily FTSE and all
frequencies of S&P 500"; "U preserves these rejections, albeit at a
slightly weaker significance level for monthly S&P 500") once we corrected
the order of rows and columns.

### Why this can't be a bit-exact numeric cross-check test

Unlike STADF, this paper does not publish a table of fixed asymptotic
critical values. The purpose of the method is that size correction always
goes through a wild bootstrap for each dataset (M=499 replications here),
and the bootstrap depends on the random number generator and on the data.
The authors computed the p-values in Table 1 once, on their FTSE and S&P
500 series with their bootstrap draws. We cannot reproduce them bit for bit
without the original price series and the original RNG stream. The paper
reports its finite-sample size and power results (Figures 3-4) graphically
and not in tables. A faithful cross-check test therefore has to be
tolerance-based. For example, it can verify that the empirical size under
H0 is close to nominal in our own large-M simulation. That is a weaker
guarantee than an exact-value assertion, which was possible for the fixed
critical values of STADF.

### Implementation

We needed three pieces. The first is a Gaussian-kernel nonparametric
volatility estimator with leave-one-out CV bandwidth selection (eq. 6;
footnote 2 gives the bandwidth search range `h ∈ [1/(2T), 1/6]`). The
second is a WLS-weighted recursive Dickey-Fuller statistic (the equation in
section 3), which needs an intercept and heteroskedasticity-weighted least
squares for each window, and so differs structurally from exuber's existing
OLS `radf()` and from the no-intercept `radf_tt()`. The third is exuber's
existing wild-bootstrap machinery (`radf_wb.R` and `radf_wb_cv`, which
already implements the Harvey/Leybourne/Sollis/Taylor wild bootstrap for
supDF). We extended it to run supDF and supBZ jointly on the same bootstrap
draws, which the validity of the union procedure requires (Theorem 3), and
to compute the union scaling constant. This is what `radf_sbz.R` (200
lines) now does.

### Independent validation (2026-08-09), found and fixed a real bug

We cannot reproduce Table 1 bit for bit (see above), so the independent
check we can run is the empirical size under the null. We simulate many
pure random walks with no bubble, run `radf_sbz_cv()` on each, and check
that the bootstrap rejects at roughly the nominal rate, neither much more
nor much less.

The first run (seed=13579, n=150, nrep=150 replications, nboot=199,
nominal 5%) showed gross oversizing:

| Statistic | Empirical rejection rate under H0 (target ≈ 0.05) |
|---|---|
| supDF | **0.640** |
| supBZ | 0.080 |
| U | **0.573** |

A supDF that rejects a true null 64% of the time, where it should reject
about 5% of the time, is not finite-sample noise, and the test was broken.
The supBZ rate of 0.08 came from a different code path and was much closer
to nominal. That told us where to look: something specific to the supDF and
bootstrap-indexing side of `radf_sbz_cv()`, and not the WLS and kernel
machinery that both statistics share.

The root cause was in the bootstrap loop of `radf_sbz.R`, which computed

```r
pointer <- length(ystar) - 1L - minw
boot_df[b] <- rls_gsadf(unroot(ystar), min_win = minw)[pointer + 2]
```

The correct convention, used consistently elsewhere in the codebase (for
example `radf_wb_hlst()` in `radf_wb.R`), is `pointer <- nr - minw`, with
no `-1L` term. The bootstrap replicate `ystar` from `radf_wb_dgp_hlst()` has
the same length as the original series. We confirmed this by reading its
source: `ystar <- c(0, cumsum(w * dy))`, which has length `1 + length(dy) =
length(y)`. The stray `-1` therefore shifted the index by exactly one
position into the flat vector that `rls_gsadf()` returns. The indexing
convention in `radf_.R` is `results[pointer+1]`=adf,
`results[pointer+2]`=sadf and `results[pointer+3]`=gsadf. The bootstrap
loop was silently extracting adf, a single end-of-sample t-statistic, into
a variable meant to hold draws of sadf, the running maximum over the whole
recursive sequence. That quantity is systematically much smaller. We had
correctly computed the observed `supDF_obs` as `sadf`, but compared it with
a bootstrap critical value drawn from the wrong, much too low,
distribution, and this caused gross over-rejection. `U`, defined as
`max(supDF, ratio*supBZ)`, inherited the problem through its `supDF` term.

The fix in `exuber/R/radf_sbz.R` is one line:

```diff
-      pointer <- length(ystar) - 1L - minw
+      pointer <- length(ystar) - minw
```

Revalidation after the fix, with the same seed and parameters:

| Statistic | Before fix | After fix | Target |
|---|---|---|---|
| supDF | 0.640 | **0.033** | 0.05 |
| supBZ | 0.080 | 0.080 (unchanged, as expected) | 0.05 |
| U | 0.573 | **0.060** | 0.05 |

supDF and U are now in the right neighborhood of the nominal size. The bug
never affected supBZ, and its 8%, a mild oversizing, is plausibly
finite-`nboot` noise. We did not investigate it further. `test-sbz.R` gave
3 passed and 0 failed (1 skipped, CRAN-only) after the fix, unchanged from
before, because none of the existing tests exercised this path with enough
replications to notice a distortion of size by a factor of two or three.

Rerun (2026-09-14): `radf_sbz_cv()` was split on 2026-08-22, and the
p-values now come from `radf_sbz_union()`, whose bootstrap RNG draws are
ordered differently. The archived `radf_sbz_validation.R` (which now calls
`radf_sbz_union()`) gives supDF 0.033, supBZ 0.060 and U 0.053 under the
same seed. The supDF value is unchanged. The 8% for supBZ above is
therefore `nboot`/`nrep` noise, and all three rates are within Monte Carlo
noise of the nominal level.

Status note (2026-08-09): we found and fixed the bug through the
independent validation described here. The fix was still uncommitted in the
exuber repo at the time, and we had not yet compared our results with the
FTSE and S&P p-values in the paper's Table 1, which needs the original
price data and cannot be done with a size simulation alone. The
size-distortion bug that would have produced wrong p-values on any real
dataset is now fixed. If this work resumes, the next steps are to commit
the fix and, if we can source the original FTSE and S&P series, to compare
with Table 1 directly, which would be a stronger check than the H0
simulation above.

Replication script:
[replication/volatility-robustness/radf_sbz_validation.R](#script-radf_sbz_validation).

---

## Pedersen & Schütte sieve bootstrap

Status: evaluated, and mostly covered already.

### Source

Pedersen, T.Q. & Montes Schütte, E.C. Testing for Explosive Bubbles in the
Presence of Autocorrelated Innovations. Journal of Empirical Finance, 58
(2020), 207-225. We read the open working paper (CREATES Research Paper
2017-9, Aarhus University): `pure.au.dk/ws/files/109652663/rp17_09.pdf`.

### Finding

exuber already has a sieve bootstrap (`R/radf_sb.R`, `radf_sb_cv()`). It
cites a different paper (Pavlidis et al. 2016) but uses essentially the
technique that Pedersen & Schütte propose independently for the same
problem, autocorrelated innovations in the recursive right-tailed unit root
test. We fit an AR(lag) sieve to the first differences, resample the
residuals and reconstruct bootstrap paths with `stats::filter(..., "rec")`.
So this item is largely implemented already and is not a gap.

We found one specific difference. Pedersen & Schütte's contribution is to
show that a fixed lag order causes size distortion under autocorrelated
innovations, and that automatic BIC-based lag-order selection (with a
`kmax`, re-selected as part of the procedure) fixes it (their Section 4
simulations, for example Table 4.8). `radf_sb_cv()` used to take only a
single fixed `lag` argument, with no automatic AIC or BIC selection.

Update (2026-08-09, Bundle 1): implemented. `radf_sb_cv(type =
"aic"/"bic", max_lag = ...)` now reuses the existing `lag_select()`
machinery in `R/radf_wb.R`, which `radf_wb_ps_cv()` already uses. It selects
the lag for each series by AIC or BIC and takes the maximum across the
panel. The rest of the pointer and matrix-dimension logic of `radf_sb_()`
assumes one common lag for the whole panel, in line with the single-`lag`
API of `radf()`, so this was the natural scope and we did not allow the lag
to vary by series. We tested it in a new `test-sb.R`, since no test file for
`radf_sb_cv()` existed before. `type = "fixed"` is unchanged. `type = "bic"`
picks up a nonzero lag on AR(2)-autocorrelated data. On pure random-walk
data the modal selection is 0, as it should be. We checked this over 8
independent draws and did not assert it on one, because BIC can pick a
nonzero lag by chance in any single finite sample.

Replication script:
[replication/volatility-robustness/radf_sb_cv_aic_bic_validation.R](#script-radf_sb_cv_aic_bic_validation).

---

## Hafner skewness-corrected wild bootstrap

Status: done (2026-08-09). Shipped as `radf_wb_cv(..., dist_skew = TRUE)`
and `radf_wb_distr(..., dist_skew = TRUE)`.

### Source

Hafner, C.M. (2020). Testing for Bubbles in Cryptocurrencies with
Time-Varying Volatility. Journal of Financial Econometrics, 18(2), 233-249.
SSRN (abstract id 3105251) returned HTTP 403 to automated fetching at
first. We recovered the paper from EconStor (IRTG 1792 Discussion Paper
2018-005, Humboldt University Berlin).
We read the full PDF, through the wild bootstrap algorithm and its Monte
Carlo section, and not only the abstract.

### What it is

Hafner modifies the wild-bootstrap multiplier distribution of Harvey et al.
(2016) so that it approximates the distribution of the PWY test better when
returns have time-varying volatility and are also right-skewed. The paper's
motivating case is cryptocurrency returns. The multiplier is (their "Step
1")

```
u_t, v_t ~ iid N(0,1), independent
w_t = u_t/sqrt(2) + (v_t^2 - 1)/2
```

By construction `E[w_t]=0`, `E[w_t^2]=1` and `E[w_t^3]=1`. It is a fixed
right-skewed multiplier and is not adaptively matched to the empirical
skewness of each series. It replaces the symmetric `N(0,1)` or Rademacher
multiplier in an otherwise unchanged wild bootstrap in the style of Harvey
et al. (2016) (`y*_t = w_t * OLS-residual_t`, refit recursively). This is
the small, contained change that this section predicted before we read the
primary source.

### Implementation

Hafner's bootstrap shipped as a new `dist_skew` argument threaded through
the existing wild-bootstrap machinery in `exuber/R/radf_wb.R`. It needed no
new file and no new statistic.

- `radf_wb_dgp_hlst(y, dist_rad, dist_skew = FALSE)` gained a third
  multiplier branch that implements `w_t` above, alongside the existing
  `N(0,1)` branch (the default) and the Rademacher branch (`dist_rad =
  TRUE`).
- `radf_wb_hlst()`, `radf_wb_cv()` and `radf_wb_distr()` all gained a
  passthrough `dist_skew` parameter (default `FALSE`, so existing calls are
  unaffected). They check that `dist_rad` and `dist_skew` are not both
  `TRUE`.

Independent validation:

1. Moment check. We simulated `w_t` directly (500k+ draws) and found
   `E[w]~0`, `E[w^2]~1` and `E[w^3]~1`, which matches the construction
   stated in the paper.
2. Regression check. With `dist_skew = FALSE`, the code reproduces the
   pre-change wild bootstrap DGP bit for bit for the same seed. The option
   is purely additive and does not change the default path.
3. Power. Under a clear, mildly explosive alternative with ordinary
   (non-skewed) innovations, `dist_skew = TRUE` rejects 86.7% of the time
   (30 replications). The skewed multiplier does not harm basic detection
   ability.
4. Size under H0 with the paper's own right-skewed innovation distribution
   (negative log-chi-square(1), the exact distribution used in the Monte
   Carlo of the paper, footnote 1). The empirical rejection rate is 3.3%
   (60 replications, nominal 5%) with no added heteroskedasticity, and 0%
   (80 replications) with a deterministic volatility pattern added on top.
   Both are conservative and neither is oversized. This agrees with the
   paper's finding that "the test is undersized in small samples, with the
   bias increasing with the degree of global heteroskedasticity."
5. A DGP-sensitivity finding that we record as it came out. An earlier,
   more aggressive power check combined log-chi-square innovations with the
   explosive alternative and gave 0% power. Isolating the two factors
   showed that the heavy right tail of the log-chi-square distribution
   itself was responsible. Occasional extreme single-step outliers occur
   (`-log(Z^2)` blows up whenever `Z` is near zero), and in a short series
   one such outlier can dominate the SSR enough to mask a modest (`rho =
   1.06`) explosive drift, in both the observed statistic and the bootstrap
   replicates, and so wash out the signal. This is a property of power
   under an extremely heavy-tailed noise process at short `T`, and not a
   defect in the implementation of `dist_skew`. We confirmed that by
   checking that the identical alternative DGP with normal innovations
   gives strong power (finding 3, above).

The new tests are in `test-cv.R`. We extended the existing `radf_wb_cv()`
test block and did not add a new file, matching where its sibling
`dist_rad` is already tested. The full package suite passes.

Replication scripts:
[replication/volatility-robustness/hafner_dist_skew_moments_and_regression.R](#script-hafner_dist_skew_moments_and_regression),
[hafner_dist_skew_power_and_size.R](#script-hafner_dist_skew_power_and_size).

---

## Sign-based sGSADF

Status: done (2026-08-09; level-shift robustness and the demeaned variant
added 2026-08-11). Shipped as `radf_sign()` / `radf_sign_cv()` and
`radf_sign_dm()` / `radf_sign_dm_cv()`. We read the full PDFs, through
Theorem 2, Table 1 and Remark 1 of HLZ 2020 and Theorems 2-3, Remark 4 and
Table 1 of HLTZ 2025, and not only the abstracts. See "What it is" and
"Implementation" below.

### Source

Harvey, D.I., Leybourne, S.J. & Zu, Y. (2020). Sign-based unit root tests
for explosive financial bubbles in the presence of deterministically
time-varying volatility. Econometric Theory, 36(1), 122-169.
`doi:10.1017/S0266466619000057`. Access was blocked for a long time:
Cambridge Core, the Nottingham Repository and the author's own homepage all
led nowhere (see the history below). It opened when we retried through a UK
academic network (Jisc), which gave access to Cambridge Core directly.

We recovered its 2025 Oxford Bulletin extension in the same way: Harvey,
Leybourne, Tatlow & Zu, "Unit Root Tests for Explosive Financial Bubbles in
the Presence of Deterministic Level Shifts," OBES 87(5),
`doi:10.1111/obes.12668`.

### Access history (kept for anyone hitting the same wall without institutional access)

- The Cambridge Core abstract page was paywalled and gave no full text
  (until we had institutional access, as above).
- The Nottingham Repository (worktribe) record page returned HTTP 403.
  This repository resisted every route we tried in this project. See
  [references.md](/replication/references#the-two-that-automation-couldnt-get-hls-and-hlw).
- Yang Zu's homepage (`sites.google.com/site/zuyang`) lists the paper with
  only a DOI link and no PDF or code link. The 2019 SBZ paper, in contrast,
  links a Google Drive code archive.
- The freely available STADF paper (Kurozumi, Skrobotov & Tsarev,
  arXiv:2012.13937) summarizes the sign-based test qualitatively in its
  literature review (Section 1) and includes it as a comparison method
  (labelled "S") in its own Monte Carlo tables. This gives indirect,
  qualitative confirmation of its properties: size control does not require
  a bootstrap, good power needs a bootstrap union-of-rejections with SADF,
  and the test is "computationally expensive, and computation time
  increases rapidly if the sample size increases." It does not give the
  exact test statistic formula or the critical values.
- B. Tatlow's own site (btatlow.com) lists the OBES extension but states
  "Pdf currently upon request" and gives no download link. The University
  of Macau economics department page (Yang Zu's affiliation) links only
  the DOI.

### What it is

The core insight is structural and not statistical. `sign(Delta y_t)` is
exactly invariant to the volatility of `Delta y_t`, because it depends only
on which side of zero the innovation fell. A test built entirely on
cumulated signs and not on the raw series is therefore exactly invariant to
any volatility pattern, even a wildly time-varying one, and needs no
bootstrap at all. HLST's wild-bootstrap correction for the standard PSY
test, which `radf_wb_cv()` already implements, does need one.

Concretely, let `C_t := sum_{i<=t} sign(Delta y_i)`. The sign-based
statistic (`sPSY`, and its single-supremum special case `sPWY`, eq. 4) is
the same double-supremum recursive Dickey-Fuller construction that PSY
uses, applied to `C_t` and not to `y_t`, and fit without an intercept
(`C_t = rho(r1,r2) * C_{t-1} + e_t`, whereas PSY's own regression includes
an intercept). Theorem 2 proves that the null limiting distribution of
`sPSY` does not depend on the volatility process `sigma(s)` at all, which
is exact invariance. Remark 1 adds that this holds because the test
excludes an intercept and works only with `sign(Delta y_t)`, and not
because of anything specific to bubbles. A hypothetical no-intercept
version of the standard PSY test would be volatility-invariant only
asymptotically, whereas `sPSY` is exactly invariant even in finite
samples.

### Implementation

Two structural facts made this the cheapest item that we validated in this
batch.

1. The no-intercept recursive-DF machinery that it needs already exists.
   We built it for STADF: `gls_dfstat_grid()` in `radf_tt.R` fits a
   GLS-demeaned AR(1) with no separate intercept and a recursively
   re-estimated residual variance, which is the structural form that the
   `sDF(r1,r2)` of `sPSY` needs. An empirical check first ruled out
   reusing the `unroot()` and `rls_gsadf()` pathway of `radf()` directly.
   The column names of `unroot(y, lag = 0)` suggest that there is no
   intercept, but we confirmed that it computes a with-intercept
   regression, which matches `lm(diff(y) ~ y[-length(y)])` exactly. That is
   the wrong form for this test.
2. The null distribution is pivotal (Theorem 2), so we simulate critical
   values once by Monte Carlo under a plain random walk, and never per
   dataset. `radf_tt_cv()` already uses the same pattern for STADF.

The code is in a new file, `exuber/R/radf_sign.R`, which closely mirrors the
structure of `radf_tt.R`.

- `sign_transform(y) := c(0, cumsum(sign(diff(y))))`.
- `radf_sign(data, minw)` transforms each series and calls
  `gls_dfstat_grid()` unchanged. It returns a `radf_sign_obj`, which
  inherits `radf_obj` and is therefore compatible with the same downstream
  S3 machinery that the `radf_tt_obj` of STADF uses.
- `radf_sign_cv(n, minw, nrep, seed)` computes Monte Carlo critical values.
  Its structure is identical to that of `radf_tt_cv()`, with the transform
  swapped.

We scoped two pieces out of this pass. Both are separate and larger pieces
of work, and the core test is useful without them.

- The paper's union-of-rejections, which combines `sPSY`/`sPWY` with the
  standard `PSY`/`PWY` test through a wild bootstrap (Section 4). It
  captures "most of the power available from the better performing of the
  two tests," since `sPSY` does not dominate the standard test on power in
  every specification. It has the same union structure as [SBZ's `U`
  statistic](#sbz-wls--kernel-volatility), which can serve as a precedent
  if someone picks this up.
- The bubble-dating extension of Section 6 that uses the sign-based
  statistic (a `datestamp()` analogue built on `C_t`). We did not attempt
  it.

Independent validation:

1. Exact invariance to heteroskedasticity, the paper's central claim and
   the cleanest possible test of it. We took a series with a given sign
   pattern, left it alone once (constant volatility) and multiplied it once
   by a wildly time-varying volatility pattern (`0.1`/`10`/`1` across three
   segments). The two versions give bit-identical `sadf` and `gsadf`
   statistics. They are not approximately similar. They are identical, as
   the exact-invariance theorem predicts.
2. Cross-check against the published finite-sample critical values of
   Table 1, at the exactly matching finite `T` (and not an asymptotic
   approximation, see the note below on why). At `T = 200` and
   `minw/T = 0.1`, the simulated `gsadf_cv` (`sPSY`) is `(3.482, 3.900,
   4.925)`, against the published `(3.469, 3.901, 4.957)` for (10%, 5%,
   1%). The 5% value matches to three decimal places. The simulated
   `sadf_cv` (`sPWY`) is `(2.337, 2.665, 3.403)`, against the published
   `(2.405, 2.735, 3.434)`, which is also close. At `T = 100` the match is
   looser (the published 1% value for `sPSY` is 13.056, an artifact of the
   extreme-value statistic at small samples and high quantiles that the
   paper's text flags as expected).
3. A methodological note worth recording. Our first attempt compared a
   `T = 300` simulation with the asymptotic (`T = Inf`) row of the paper's
   table. `sPWY` matched well and `sPSY` did not (off by about 0.5 at the
   5% level). This was not an implementation error. The paper's Table 1
   documents that the finite-sample critical values of `sPSY` converge to
   their asymptotic limit much more slowly than those of `sPWY` (its text:
   "convergence... is fairly slow (particularly for `sPSY`)"), and the
   `T = 300` value sits almost exactly between the paper's `T = 200` and
   `T = 400` rows once we compare it correctly. Repeating the cross-check at
   the exactly matching finite `T` (point 2 above) removed the ambiguity.
4. Power. The empirical rejection rate is 96.7% (30 replications) under a
   clear, mildly explosive alternative, using the simulated 95% critical
   value of `radf_sign_cv()`.

The new tests are in `test-sign.R`. They cover a formula check against a
brute-force `lm()`, the exact-invariance property, the cross-check against
the published values at `T = 200`, and a power check. The full package
suite passes.

Replication scripts:
[replication/volatility-robustness/sign_based_invariance_and_power.R](#script-sign_based_invariance_and_power),
[sign_based_finite_T_crosscheck.R](#script-sign_based_finite_T_crosscheck).

### Level-shift robustness (HLTZ 2025)

Status: done (2026-08-11). This shipped as a documentation addition to
`radf_sign()` and `radf_sign_cv()`, and as one new function pair,
`radf_sign_dm()` and `radf_sign_dm_cv()`.

Source: Harvey, D.I., Leybourne, S.J., Tatlow, D. & Zu, Y. (2025). Unit
root tests for explosive financial bubbles in the presence of
deterministic level shifts. *Oxford Bulletin of Economics and Statistics*,
87(5), 879-901. `doi:10.1111/obes.12668`. We recovered it through the same
UK academic network (Jisc) route as HLZ (2020) itself (see the "Access
history" note above).

What it is: the paper evaluates the same HLZ (2020) sign-based statistics
under a different violation of the null, deterministic level shifts in the
series, and not time-varying volatility. We read the full PDF (rendered
pages 3-6, since raw-text extraction scrambles the theorem statements, as
usual). The paper's model adds `n_T` deterministic level shifts to the
standard null, with `n_T = O(T^alpha_n)`. Three theorems give the
large-sample null distribution of the standard `PSY` statistic (Theorem 1)
and of both HLZ sign-based statistics. The first is the plain
cumulated-sign `sPWY`/`sPSY` that we already shipped as `radf_sign()`
(Theorem 2). The second is HLZ's recursively demeaned analogue, which the
paper denotes `s̄PWY`/`s̄PSY` (Theorem 3) and which we had not implemented.
The key finding is that `PSY`'s validity needs a joint restriction on both
the number and the magnitude of the shifts (Assumption 3). By the paper's
Table 1, `PSY` is essentially never correctly sized once the number of
shifts grows at rate `sqrt(T)` (their Case 1), with empirical size up to
0.425 against a nominal 0.05. Both sign-based statistics, in contrast,
need only a restriction on the number of shifts (Assumption 4, `alpha_n <
1/2`) and no restriction on the magnitude of the shifts at all. Only at the
boundary rate `alpha_n = 1/2` do they pick up a term that depends on the
level shifts, and even then the degree of over-sizing is bounded
independently of the shift magnitude (Remark 4). That is a categorically
weaker vulnerability than that of PSY.

Implementation: two structural facts kept the work small.

1. `radf_sign()`, already shipped, turned out to be exactly HLZ's
   `sPWY`/`sPSY` statistic. Theorem 2 is a level-shift robustness result
   about code that already existed, so it gave us no reason to write
   anything new. The change is a pure documentation addition, a
   `Level-shift robustness` `@section` in the roxygen docs of
   `radf_sign()`.
2. The one new item, HLZ's second sign-based analogue, looked from an
   abstract-level read as though it needed a recursive-demeaning step for
   each `(r1,r2)` window, which would be expensive and shaped like a double
   recursion. The actual formula in Theorem 3 shows that `C̃_t` is demeaned
   by an expanding-window mean computed once up front (`C̃_t :=
   sum_{i=2}^t {sign(dy_i) - (i-1)^{-1} sum_{j=2}^i sign(dy_j)}`), and not
   per candidate window. The sum `sum_{j<=i} sign(dy_j)` is the running sum
   that `radf_sign()` already computes, so the demeaning term at each `i`
   is simply `C_i / (i-1)`. This is a one-time `O(T)` transform,
   `sign_demean_transform()` in `radf_sign.R`, and it feeds into the same
   `gls_dfstat_grid()` machinery that `radf_sign()` already uses. We needed
   no new estimation machinery.

We found a real bug along the way. It was not a validation failure of this
item and was a latent defect in code that had already shipped. Calling
`print()` on a `radf_sign_cv()` or `radf_tt_cv()` object crashed with `no
applicable method for 'tidy_radf_cv'`. Both functions tag their output with
a distinguishing class (`"sign_cv"`, `"tt_cv"`), and no matching
`tidy_radf_cv.*`, `summary_radf.*` or `index_radf_cv.*` method exists
anywhere in the package. Only `mc_cv`, `wb_cv` and `sb_cv` have them. The
output of these functions has exactly the same shape as that of
`radf_mc_cv()`, a plain `adf_cv`/`sadf_cv`/`gsadf_cv` quantile vector with
no per-dataset or panel dimension. The root-cause fix was therefore to add
`"mc_cv"` as an additional class tag to `radf_sign_cv()`, `radf_tt_cv()` and
the new `radf_sign_dm_cv()`, so that all three fall through to the existing
`tidy_radf_cv.mc_cv` method and we do not duplicate it. We caught this only
because the smoke test of this pass happened to call `print()` on the new
`radf_sign_dm_cv()` object. A function that is shipped without ever calling
`print()` on its own return value can hide this class of bug indefinitely.

Scoped out: HLZ's serial-correlation correction for `sPSY` and `s̄PSY`
under level shifts (Remark 6, which augments regression (4) with lagged
`Delta C_t`, because the original HLZ correction neglects the shifts). We
did not attempt it, because the existing sign-based implementation in this
project does not handle serial correlation for the no-shift case either.

Independent validation:

1. Formula check. `sign_demean_transform()` matches a brute-force per-`i`
   loop that recomputes the recursive mean from scratch, to about `1.8e-15`.
2. Reproduction of the paper's own Table 1, which is the strongest check
   available. In Case 1 (`alpha_n = 0.5`, the worst case for `PSY`), with
   `k = 2`, `mu = 5`, `p = 0.8` and `T = 400`, the published empirical size
   at the nominal 5% level is `PSY = 0.337`, `sPSY = 0.121` and `s̄PSY =
   0.050`. Our own replication (`nrep = 500`, our own simulated critical
   values under the no-shift null, the same `minw = floor(0.1T)` trimming
   that the paper uses, and the same construction of shift count,
   magnitude and location) gives `PSY = 0.302`, `sPSY = 0.090` and `s̄PSY =
   0.036`. The ordering and order of magnitude agree with the paper
   throughout. `PSY` is badly oversized, `sPSY` is mildly oversized, and
   `s̄PSY` is closest to nominal and even mildly conservative. The
   differences are well within the Monte Carlo noise at `nrep = 500`
   relative to the paper's own `nrep = 2000` numbers.
3. Exact invariance to heteroskedasticity carries over. `radf_sign_dm()`
   reuses the same `sign()`-based construction, so it keeps the property
   of `radf_sign()` that rescaling leaves the statistic bit-identical. The
   test is the same as the one for `radf_sign` and is applied to the
   demeaned variant.
4. Power. `radf_sign_dm()` rejects a clear, mildly explosive alternative
   with its own simulated critical value (same structure as the power check
   of `radf_sign()`).

We appended new tests to `test-sign.R`: a formula check, the invariance
check, a power check, and a regression test for the `print()` crash that
covers all three affected `_cv` functions. The full package suite passes.

Replication script:
[replication/volatility-robustness/radf_sign_dm_levelshift_validation.R](#script-radf_sign_dm_levelshift_validation).

### `datestamp()`/`autoplot()` parity (done, 2026-08-18)

This is the same gap and the same fix as for
[STADF/GSTADF](#time-transformed-test-stadf--gstadf). `radf_sign_cv()` and
`radf_sign_dm_cv()` computed only the three scalar critical values that
`summary()` and `tidy()` need. They discarded the `badf`/`bsadf` path that
their own calls, `gls_dfstat_grid(sign_transform(y), minw)` and
`gls_dfstat_grid(sign_demean_transform(y), minw)`, already compute in each
replicate. We checked this right after confirming the fix for
`radf_tt_cv()`. All three functions call `gls_dfstat_grid()` on a
transformed series, and the identical construction applies unchanged. We
needed no cummax shortcut, because the `bsadf` of each replicate is already
the full sup-over-all-window-starts statistic, whichever series feeds the
grid.

We ran the same battery as for STADF/GSTADF on each function.

1. Formula-exact check. The last row of `badf_cv` is bit-identical to
   `adf_cv` for both `radf_sign_cv()` and `radf_sign_dm_cv()`. This is a
   hard identity (`adf <- badf[length(badf)]` inside `gls_dfstat_grid()`)
   and holds regardless of the sign transform.
2. Monte Carlo size (pure random walk, n=100, minw=20, nrep=2000, 200
   replications). The false-alarm rate is 5.5% for `radf_sign` and 3.5% for
   `radf_sign_dm`, against a nominal 5%, so we see no size distortion.
3. Power, cross-checked against the same baseline. On the identical
   synthetic bubble used for the STADF/GSTADF check (the `radf()` and
   `radf_mc_cv()` baseline detects 16%), `radf_sign` detects 20% and
   `radf_sign_dm` detects 8%. The lower number for `radf_sign_dm` agrees
   with the trade-off between heteroskedasticity invariance and power that
   is documented above for this whole family of tests, and it is not a
   validation concern.
4. Both `datestamp()` (`option = "gsadf"` and `"sadf"`, which exercise
   `bsadf_cv` and `badf_cv` respectively) and `autoplot()` run without
   error on `radf_sign(sim_data)` and `radf_sign_dm(sim_data)`.

We added new tests to `test-sign.R` (a shape and identity check, and a
full-pipeline smoke test for both functions). The full package suite passes
(724 assertions, up from 700 before this `badf_cv`/`bsadf_cv` pass).

---

## Stochastic explosive-coefficient test

Status: done. SSU was implemented on 2026-08-10 (`ssu_test()`). GSSU, the
UR/GUR union-of-rejections procedure and the four CUSUM-type statistics
followed on 2026-09-29 (`ssu_test(type = "gssu", union = TRUE)`,
`cusum_test()`). We read the full PDF (the model, and Section 3's test
statistics through the union-of-rejections and CUSUM/CUSUM-SQ proposals)
and re-verified the implemented items against rendered PDF pages 5-6 and 9.

### Source

Kurozumi, E. & Nishi, M. (2025). "Testing for a bubble with a
stochastically varying explosive coefficient." *JTSA*, 46(5), 945-965. The
article is open access. Wiley requires JavaScript to serve it, which is why
it first looked paywalled to a plain `curl` fetch.

### What it is

The paper models the explosive AR(1) coefficient itself as random,
`1 + c1/T + a*u_t/sqrt(T)` (eq. 2). Every other method in this project
assumes the deterministic form `1 + c/T^alpha`. Here `u_t` is i.i.d. with
mean zero and unit variance and is independent of the innovations. The
motivation (Figure 1: rolling-window AR(1) estimates on Japanese daily
stock prices) is that the explosive speed itself looks unstable in
practice and is not constant during a bubble episode.

This is not a volatility-robustness fix in the sense of the rest of this
file. It does not touch the innovation variance, since `sigma_epsilon^2` is
constant (`Assumption 1a`). It generalizes the model in a different
direction: the persistence, or explosiveness, parameter is random, and not
the noise scale. We group it in this file because we first found it next to
the other heteroskedasticity items in the *JTSA* special issue, and not
because it addresses the same failure.

The paper proposes three structurally different new statistics and not one.

1. SSU/GSSU (eq. 7-8) is a stochastic-unit-root test in the style of Lee
   (1998) and Nagakura (2009). It is built on an entirely different
   regression, `(Delta y_t)^2 = mu^2 + eta*y_{t-1}^2 + e_t`, with squared
   differences on squared lagged levels, and not the standard ADF
   regression. Its raw t-statistic is not asymptotically pivotal on its
   own. We have to subtract a bias-correction term `rho_hat(r1,r2)`, a
   recursively estimated cross-moment between the residuals of the level
   regression and those of the squared regression, before the corrected
   statistic `t^c_{r1,r2}` is usable, following Nishi & Kurozumi (2024).
2. CUSUM and CUSUM-SQ are a second, independent pair of statistics (the
   classic parameter-constancy tests of Brown et al. 1975), studied here as
   a third and fourth candidate for the same detection problem.
3. A union of rejections combines (at least) SADF/GSADF with SSU/GSSU. This
   is the paper's own recommended practical procedure, because neither
   family dominates. SSU/GSSU wins when the coefficient is stochastic (`a
   != 0`), and SADF/GSADF wins when it is deterministic (`a = 0`). The
   union needs joint critical values. It has the same structure as [SBZ's
   union statistic](#sbz-wls--kernel-volatility) or the evaluation in
   Bundle 5 of a union of the sign-based test (not implemented) with
   standard PSY.

### Cost/feasibility note for exuber

This is not a contained addition, for reasons that differ from those of
every other item in this bundle and that are larger.

1. It needs a new regression form. The `(Delta y_t)^2` on `y_{t-1}^2`
   regression of SSU/GSSU shares no structure with the `y_t` on `y_{t-1}`
   ADF-family regression that every existing exuber statistic (`radf()`,
   STADF, kernel-purge, SBZ, sign-based) builds on. We would have to write
   new estimation code from scratch, and could not feed a transform into
   `gls_dfstat_grid()` or `rls_gsadf()` the way we could for the
   sign-based test.
2. The bias correction is itself nontrivial. `rho_hat(r1,r2)` requires a
   recursive estimate of a cross-moment between two different residual
   series (from the level and squared-level regressions) for every
   candidate window. That is a more involved recursive computation than
   any cumulative-sum trick used elsewhere in this project. PDC/KS, STADF
   and the sign-based test all reduce to a handful of `cumsum()` calls,
   whereas this one needs the residuals of two fitted regressions for each
   window first.
3. There are two more statistics, and not one. CUSUM and CUSUM-SQ are a
   further, separate pair of test statistics in the same paper.
   Implementing "the stochastic-coefficient test" faithfully would mean
   either choosing which of up to four statistic families (SSU, GSSU,
   CUSUM, CUSUM-SQ) to ship, or implementing all of them.
4. The recommended procedure is a union, which needs its own joint
   critical-value simulation. This is the same category of extra work as
   the union statistic of SBZ (flagged there as "a distinct, non-trivial
   chunk of new statistical code"), and it comes on top of the cost of the
   new statistic and does not replace it.

The scope is comparable to SBZ (two new estimators plus a union and
bootstrap layer) and not to the other items in this bundle (Hafner: one new
multiplier branch; sign-based: one new input transform that reuses
existing machinery). We did not take it up in this pass. If it is
revisited, SSU alone, without the double recursion of GSSU, without
CUSUM/CUSUM-SQ and without the union, would be the minimum viable first
cut.

We re-triaged this on 2026-08-10 by re-reading rendered pages 5-6 and 9.
Points 1 and 2 above overstate the cost for the SSU-alone first cut.

- Point 1 (a new regression form) is real but contained. The SSU
  regression `(Delta y_t)^2 = mu2 + omega*y_{t-1}^2 + eta_t` (eq. 7) is a
  plain two-variable OLS over a window. It is covered by the generic
  closed-form pattern of `hls_prefix_sums()` and `hls_segment_coef()`
  (an arbitrary `(x, z)` pair over a segment) that `dating_hls()` and
  `dating_knp()` already use. The only change is that `x = y_{t-1}^2` and
  `z = (Delta y_t)^2`, in place of `x = y_{t-1}` and `z = Delta y_t`. No new
  estimation theory is needed, only a different pair of input series.
- Point 2 (the bias correction) is real but also closed-form, and it is
  not "recursive" in the expensive sense. The cross-moment
  `sigma_hat_{eps*eta}` needs "the residuals from two fitted regressions,"
  which sounds as if it needed two passes of per-observation residual
  computation for each window. Expanding `sum(eps_hat_t * eta_hat_t)`
  algebraically reduces it to a bilinear combination of window sums of
  twelve fixed per-observation products (`y_{t-1}`, `y_{t-1}^2`, ...,
  `y_{t-1}^4`, `Delta y_t`, ..., `(Delta y_t)^4`, and their cross
  products). That is more cumulative sums to track than in any previous
  item in this project (twelve, against the usual four to six), but it
  still takes `O(1)` per window through `cumsum()` differences and needs no
  second pass over the raw data for each candidate window.
- The critical value is published and does not need simulation. Table I of
  Kurozumi & Nishi (their 10,000-replication Monte Carlo) gives the
  asymptotic critical value of `SSU` directly, `2.90`/`3.30`/`4.20` at the
  `10%`/`5%`/`1%` level. It is a single scalar for each level, because
  `SSU` is a single-recursion sup-statistic with `r1` fixed at `0`, unlike
  the double recursion of `GSSU`. The paper recommends `r0 = 0.01 +
  1.8/sqrt(T)`, which is exactly the existing formula of `psy_minw()` and
  can be reused without adjustment.

Points 3 and 4 (CUSUM/CUSUM-SQ and the union-of-rejections procedure) were
left out of the SSU pass. We re-triaged them on 2026-09-29, using rendered
pages 6-9, and they too turned out to be overstated. Table I publishes
critical values for every statistic in the paper (GSSU, CS, GCS, and both
tails of CSSQ and GCSSQ) and also the union scaling constants `ur` and
`gur`, so the union needs no joint simulation. Every CUSUM-type statistic is
a sup or inf of a partial-sum process, which takes `O(T)` through a running
max or min. See [the 2026-09-29 implementation](#implementation-gssu-union-cusum-done-2026-09-29)
below.

### Implementation, SSU done (2026-08-10)

SSU shipped as `ssu_test(data, minw = NULL, level = 0.95)` in
`exuber/R/ssu_test.R`. `ssu_prefix_sums()` builds twelve cumulative-sum
vectors from the two base per-observation series (`y_{t-1}` and `Delta
y_t`). `ssu_stat_path()` evaluates the bias-corrected `t^{omega,c}_{0,r2}`
statistic for every candidate end point from those sums. We verified the
algebra by hand-expanding the residual cross-moment `sum(eps_hat*eta_hat)`
into its bilinear window-sum form and then confirmed it numerically against
brute force, as described below. `ssu_test()` takes the running maximum
(the sup-statistic of `SSU`) and compares it with the lookup of
`ssu_q(level)` into Table I.

Validation: `ssu_stat_path()` matches a brute-force computation (two
separately `lm()`-fitted regressions on the raw window data plus a manual
residual cross-moment) to machine precision (`< 1e-8`) at four window
sizes. The Table I lookups are exact, with a clean error for an untabulated
level. The default `minw` of `ssu_test()` matches `psy_minw()` exactly,
which shows that the `r0` formula of `SSU` needed no adaptation. The
false-alarm rate under `H0` (300 replications, `n=200`) is
`12.0%`/`9.0%`/`2.7%` (after the cross-moment denominator fix of
2026-09-29, below), against a nominal `10%`/`5%`/`1%`. The test is mildly
oversized, in the same range as several other finite-sample-against-
asymptotic critical values validated in this project, and this alone is not
a cause for concern. Detection power on a stochastic explosive-coefficient
DGP (the alternative of eq. 2 in Kurozumi & Nishi, 60 replications) is
`85.0%`, which is the alternative `SSU` is designed for. On a deterministic
explosive DGP instead, the kind that the `SADF` of `radf()` targets, `SSU`
still has `80.0%` power, against `90.0%` for `SADF`. This is a sensible
result: `SSU` gives up a little power on the deterministic alternative for
robustness to stochasticity in the coefficient, which is the trade-off
described in the paper's Theorem 2 (each family dominates on its own
alternative, and neither dominates universally). The new `test-ssu.R` has 7
tests. Replication script:
[replication/volatility-robustness/radf_ssu_validation.R](#script-radf_ssu_validation).

### Implementation, GSSU, union, CUSUM (done 2026-09-29)

GSSU (`ssu_test(type = "gssu")`): `ssu_stat_path()` now takes a window start
as well, so the same twelve prefix sums give any window `(lo, hi]` in
`O(1)`. The GSSU path is the sup over starts at each end point (the shape
of `bsadf`), and its maximum is the statistic. The minimum window is the
paper's `r0 = -0.004 + 2.24/sqrt(T)`. The Table I note says that
`psy_minw()`'s formula oversizes GSSU. The Table I critical values are
`4.83`/`5.37`/`6.81`.

Union of rejections (`ssu_test(union = TRUE)`): `UR = max(SADF / cv_sadf,
SSU / cv_ssu)` is compared with `ur` = `1.16`/`1.13`/`1.09`. GUR, with
GSADF and GSSU, is compared with `gur` = `1.11`/`1.10`/`1.08`. The constant
is valid only at the level for which the statistic was built. The
SADF/GSADF side is `radf(x, lag = 0)` against `cv`. The default is the
precomputed store at the `n` of the data, that is, the finite-sample
critical values at `psy_minw(n)`, which approach the asymptotic ones that
Table I's `ur` and `gur` were calibrated with.

CUSUM-type statistics (`cusum_test(type = "cs" | "gcs" | "cssq" |
"gcssq")`), page 7: `S_k = sum_{t<=k} Delta y_t / (sigma sqrt(T))`, with
`sigma^2 = mean((Delta y)^2)` (not demeaned). CS is `max_k S_k`. GCS is
`max_{j<k} (S_k - S_j)`, a running-minimum drawup. For CSSQ, `D_k =
(sum_{t<=k} (Delta y)^2 - k/T sum (Delta y)^2) / (sigma_eta sqrt(T))`, with
`sigma_eta^2 = mean((Delta y)^4) - sigma^4`. CSSQ is `max_k D_k` and
`min_k D_k`, and GCSSQ is the drawup and drawdown of `D`. The CUSUM-SQ
tests are two-sided, with each tail at `alpha/2`. The CSSQ columns of Table
I sit at the Brownian-bridge sup quantiles for `alpha/2` (for example
`1.32` at 5%, against `sqrt(-log(0.025)/2) = 1.36` before discretization).
This confirms that a "level alpha" row already is the two-sided test at
alpha.

We corrected one point in the shipped SSU. Page 6 divides all three
moments (`sigma_eps^2`, `sigma_eta^2` and the cross-moment `sigma_{eps
eta}`) by the same `floor(T r2) - floor(T r1) - 1`, which is the window
count minus 2. `ssu_stat_path()` used the count minus 1 for the
cross-moment. The brute-force check had copied the same choice, so it could
not catch the error. We changed both to the count minus 2. The effect is
`O(1/T)` (for example, the 5% size of SSU below moved from 8.7% to 9.0%).

Validated (sections 7-10 of the replication script):

- Exact checks. Windows with `lo > 0` match the brute-force `lm()` with a
  manual cross-moment (`2.4e-15`), and so does the GSSU sup path
  (`6.1e-15`). All four CUSUM-type statistics, sup and inf, match a
  brute-force double loop over every window (`8.9e-16`). `test-ssu.R`
  checks every Table I column value by value.
- Size at 5%, `n = 200`, 300 replications: SSU `0.067`, GSSU `0.060`, UR
  `0.050`, GUR `0.057`, CS `0.043`, GCS `0.047`, CSSQ `0.023`, GCSSQ
  `0.037`. Section 4's separate SSU run gives `12.0%`/`9.0%`/`2.7%` at
  10/5/1%, so SSU is mildly oversized, as we recorded above.
- Power at 5%, `n = 200`, 100 replications, with the bubble over the second
  half. With a stochastic coefficient (`c1 = 3`, `a = 4`): SSU `0.90`, GSSU
  `0.97`, UR `0.87`, GUR `0.96`, CS `0.01`, GCS `0.01`, CSSQ `0.85`,
  GCSSQ `0.78`. This reproduces the paper's main qualitative finding
  (Theorem 2 and Figure 2): CUSUM-type tests lose essentially all their
  power once `a != 0`, while the SSU and CUSUM-SQ types keep it, and the
  union stays close to the better of its two parts. With a deterministic
  coefficient (`c1 = 10`, `a = 0`), every SSU, union and CUSUM-SQ statistic
  reaches `1.00`, and CS and GCS reach `0.53`.

The tests are in `test-ssu.R` and `test-cusum-test.R`. Replication script:
[replication/volatility-robustness/radf_ssu_validation.R](#script-radf_ssu_validation).
We ported the functions to pyexuber (`ssu_test(type=, union=)` and
`cusum_test()`) and cross-checked them against the R values on a shared
input.

---

## SV-ADF

Status: done (2026-08-11), `datestamp(option = "svadf")`. The source is a
preprint and has not been peer-reviewed, which is a lower bar than for every
other source implemented in this project, and we flag it explicitly. We read
the full PDF through the definition of the SV-ADF statistic, its asymptotic
theorem (3.1) and the proof appendix, which confirms the exact construction
of the feasible variance estimator. On 2026-08-11 we also re-triaged
Section 5.1's threshold-calibration exercise (rendered pages 20-22), which
resolved the open question in the cost note below.

### Source

Sarkar, A. & Wells, M.T. (2026). "Is There an AI Bubble? Robust
Date-Stamping for Periods of Exuberance." arXiv:2604.12062. The theoretical
grounding is in the same authors' arXiv:2512.06823, "Double Local-to-Unity:
Inference under Nearly Nonstationary Volatility."

### What it is

SV-ADF extends the recursive right-tailed ADF test (PWY-style, eq. 3-4) to
highly persistent, nearly nonstationary stochastic volatility, for example
`log(sigma_t^2) = phi_n * log(sigma_{t-1}^2) + eta_t` with `phi_n`
approaching 1 at an iterated-logarithmic rate. That condition is much
weaker than the deterministic or bounded-jump volatility that every other
method in this file assumes. HLST's wild bootstrap, Hafner, the sign-based
test, kernel-purge and STADF all require `sigma(s)` to be a fixed,
non-stochastic function of calendar time. SV-ADF allows the volatility
process itself to be a near-unit-root-persistent stochastic process, such
as a GARCH with `alpha+beta` close to 1 or a persistent stochastic-volatility
model.

Reading the proof appendix gave us an interesting structural finding. The
feasible `SV-ADFr` and `SV-ADFrt` statistics (eq. 3-4) are built from the
same recursive OLS estimator and the same residual-based variance estimator
that the recursive ADF t-statistic of `radf()` already computes. The
variance estimator is `tau_hat^2 = tau^-1 * sum(residuals^2)`. Eq.
A.13-A.14 confirm that it is the standard within-window OLS residual
variance and not a separate nonparametric or kernel estimate. The paper's
real contribution (Theorem 3.1) is a broader asymptotic justification for
using this same statistic under much weaker volatility conditions than were
previously proven. It is not a structurally different statistic. Theorem 3.1
notes that the limiting functional "coincide[s] with those obtained under
homoskedasticity (Phillips and Yu, 2009)" once it is normalized.

This reading did not settle one point. The abstract advertises "distinct
calibration thresholds for testing origination and collapse" and a
"moderate-deviation asymptotic theory." That language suggests that the
practical novelty of the paper is an adaptive, moderate-deviation
critical-value boundary `cv_{n,r}`, which grows with `r` at a specific
rate, and not a fixed quantile of the limiting distribution. It would
explain why SV-ADF avoids the spurious bubble that PWY flags in a
high-volatility episode (Figure 6, the empirical Nvidia example) although
the underlying point statistic is the same. To pin down the exact
construction of this boundary function, we would have to read the authors'
companion "Double Local-to-Unity" theory paper, which we did not do in this
pass.

### Cost/feasibility note for exuber

The cost was uncertain, more so than for the other three items in this
bundle. The 2026-08-11 re-triage below resolved it in the favorable
direction.

- If the practical procedure is to compute exuber's existing recursive ADF
  t-statistic exactly as `radf()` already does and then compare it with a
  specific moderate-deviation-calibrated boundary function and not a fixed
  quantile, the cost could be close to zero. We would need no new point
  statistic, only a new comparison rule. That would put it in the same cost
  tier as the sign-based test, the cheapest item in this bundle.
- If the boundary function itself requires estimating something new (for
  example a persistence parameter for the stochastic-volatility process, or
  quantities from the companion moderate-deviation paper that this paper
  does not determine), the cost is unknown until we read that paper.
- This is a preprint on arXiv (first version 2026) and has not been
  peer-reviewed. Everything else that we implemented in this project so far
  came from peer-reviewed, published papers. Whoever picks this up next
  should know that, independently of the question of implementation cost.

On 2026-08-11 we re-triaged the item by re-reading rendered pages 20-22
(Section 5.1, "Threshold Selection Insights"). The first, favorable branch
above is what actually happens, and we did not need to read the companion
paper after all.

- Section 5.1 of the paper describes how the "moderate-deviation-calibrated
  boundary" was obtained, and the boundary is not a function of any
  estimated nuisance parameter. For origination, the authors simulate the
  SV-ADF statistic under `H0` at `n ∈ {100,200,...,1000}` (1,000
  replications each), take the 90th percentile at each `n`, and find "these
  upper 90 percent critical values are well approximated by `log(n)/10`,
  which is adopted as the origination threshold." For collapse, they
  average the 10th-percentile threshold over randomly drawn
  nuisance-parameter configurations and find it "most closely approximated
  by `log(n)/2`, which we therefore use as the collapse threshold." Both
  are published, closed-form formulas in the sample size only, and they
  are part of the authors' own applied methodology. They do not require the
  companion "Double Local-to-Unity" paper, and exuber does not have to
  estimate any persistence parameter.
- Together with the structural finding already confirmed (the point
  statistic is the `badf` of `radf()`), this puts the practical procedure
  of SV-ADF at zero new-statistic cost. We reuse `badf`, compare it with
  two `log(t)`-based thresholds in place of one fixed quantile, and apply
  a first-crossing dating rule. One piece is new. Origination and collapse
  use different thresholds (the paper's Remark 1: a unit-root-based
  threshold suits origination and not collapse). The existing S3 dispatch of
  `datestamp()` assumes a single shared critical value throughout, so this
  shipped as its own small, self-contained dating routine that reuses the
  existing contiguous-run detection of `stamp()`, and not as an extension
  of `datestamp()` itself.

### Implementation, done (2026-08-11)

SV-ADF shipped as `datestamp(data, option = "svadf", min_duration = NULL)`
(helper `datestamp_svadf()` in `exuber/R/svadf.R`). It started as its own
`datestamp(option = "svadf")` entry point and was folded into
`datestamp()` on 2026-08-18. `min_duration` defaults to `psy_ds(n)`, the
existing `log(n)`-based minimum-episode-duration rule of exuber. We reused
it directly and did not invent a new one. It stands in for the paper's
data-frequency-specific requirement of "at least two consecutive calendar
months" or "one month" for consolidation. Origination is dated at the first
run of at least `min_duration` consecutive points with `badf` above
`log(t)/10`. Collapse is dated at the first run of at least `min_duration`
consecutive points with `badf` below `log(t)/2`, and we search for it only
after the origination date.

Validated: the `badf` field of `datestamp(option = "svadf")` matches a
direct `radf()` call bit for bit, which confirms that the reuse of the
point statistic is exact and not approximate. The threshold formulas match
`log(t)/10` and `log(t)/2` exactly. Collapse is guaranteed by construction
never to be dated before origination, because the search for it starts
after the origination date, and we confirmed this across 20 replications
with no exceptions. We used a synthetic bubble-and-collapse episode, built
with this project's established large-base bubble-DGP convention, because
a first attempt with a mean-zero random-walk base gave a weak, diluted
signal, as it did in the earlier validation of `dating_hls()`. The
detection rate is 100% across 20 replications. The mean absolute
origination-date error is `4.95` periods and the mean absolute
collapse-date error is `20.25` periods. Collapse detection is inherently
less precise, because the expanding `badf` window, anchored at the start,
dilutes a post-collapse downward signal more than an upward signal at
origination. This is a known property of single-recursion statistics and
not a defect. The false-alarm rate under `H0` (60 replications, pure random
walk) is `13.3%`. This counts any origination crossing anywhere in a full
150-period path, which is a much larger compound opportunity for a false
alarm than a single-point 10% test, so a somewhat inflated aggregate rate is
expected, and it does not mean the test is miscalibrated at any individual
point. The new `test-datestamp-svadf.R` has 7 tests. Replication script:
[replication/volatility-robustness/datestamp_svadf_validation.R](#script-datestamp_svadf_validation).

The AI-equity empirical application (the 2025-26 exuberance in Nvidia,
Tesla, TSMC and others) is a plausible source for a topical worked-example
vignette now that the statistic itself is implemented. See
practitioner-guidance.md.
