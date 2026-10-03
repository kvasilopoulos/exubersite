---
title: "Dating and root inference"
blurb: "Origination, collapse and recovery dates, plus confidence intervals on the explosive root itself."
order: 2
---
Two related problems come up once `radf()` or `datestamp()` has flagged an
episode as explosive. **Dating** asks when the episode started, ended, or
recovered. **Root inference** asks how explosive it was: what the
autoregressive root `rho` is, what confidence interval it carries, and how
fast the series was doubling. We keep them in one file because they are the
same step in exuber's workflow (`radf()` → `datestamp()` → sub-sample →
refine). The methods themselves share little statistical machinery. The
status legend is the one in [volatility-robustness.md](/replication/volatility-robustness).

| Method | Paper | Status |
|---|---|---|
| [SSR/BIC dating, PDC/KS route](#ssrbic-dating-vs-psy-recursive-dating) | Pang/Du/Chong (2021), Kurozumi/Skrobotov (2023) | **done** |
| [SSR/BIC dating, HLS/HLW route](#ssrbic-dating-vs-psy-recursive-dating) | Harvey/Leybourne/Sollis (2017), Harvey/Leybourne/Whitehouse (2020) | **done** (both HLS single-bubble and HLW multi-bubble routes) |
| [Root inference (Cauchy CI + normal-t CI)](#root-inference) | Phillips & Magdalinos (2007), Guo, Sun & Wang (2019) | **done** |
| [Confidence sets for bubble dates](#confidence-sets-for-bubble-dates) | Kurozumi & Skrobotov (2025) | re-triaged 2026-08-10. The critical values are cheap (closed form or a published response surface) but the statistic is still multi-step. Not implemented |
| [Improved retrospective dating](#improved-retrospective-dating) | Kejriwal, Nguyen & Perron (2025) | **done** (single-bubble omission fix 2026-08-10; multi-bubble dynamic programme `dating_knp(breaks = )` 2026-09-29) |
| [WLS dating under time-varying volatility](#wls-dating-under-time-varying-volatility) | Kurozumi & Skrobotov (2023) | **done** |
| [Reverse-regression recovery dating](#reverse-regression-recovery-dating) | Phillips & Shi (2014/2019) | **done, shipped with caveats** |

All papers: [references.md](/replication/references#dating-and-root-inference).

---

## SSR/BIC dating vs. PSY recursive dating

**Status: the PDC/KS route is done (2026-08-09). HLS (2017)'s single-bubble
SSR+BIC route is done (2026-08-10), and so is HLW (2020)'s multi-bubble
two-step wrapper around it (2026-08-10).** The PDC/KS sequential
sample-splitting estimator is implemented as `dating_pdc()`; see
[Implementation (PDC/KS route)](#implementation-pdcks-route) below. HLS's
four-model grid-search BIC approach is implemented as `dating_hls()`; see
[Implementation (HLS route)](#implementation-hls-route) below, which also
records a sign-constraint bug that we found and fixed shortly after the
first release. HLW's two-step extension is implemented as `dating_hlw()`;
see [Implementation (HLW route)](#implementation-hlw-route) below.

### Source

1. Harvey, D.I., Leybourne, S.J. & Sollis, R. (2017). "Improving the accuracy
   of asset price bubble start and end date estimators." *Journal of
   Empirical Finance*, 40, 121-138. `doi:10.1016/j.jempfin.2016.11.001`.
   ("HLS")
2. Harvey, D.I., Leybourne, S.J. & Whitehouse, E.J. (2020). "Date-stamping
   multiple bubble regimes." *Journal of Empirical Finance*, 58, 226-246.
   `doi:10.1016/j.jempfin.2020.06.004`. ("HLW")
3. Pang, T., Du, L. & Chong, T.T.L. "Estimating multiple breaks in
   nonstationary autoregressive models." *Journal of Econometrics*, 221(1),
   277-311 (2021 publication; working paper circulated 2018-2019). ("PDC")
4. Kurozumi, E. & Skrobotov, A. "On the asymptotic behavior of bubble date
   estimators." Published version: *Journal of Time Series Analysis*, 44(4),
   359-373 (2023). ("KS")

We read all four papers as full PDFs and not through secondary summaries.
HLS and HLW were for a long time the two hardest papers in this project to
obtain. Both are paywalled at Elsevier/ScienceDirect, and their open-access
mirror at the University of Nottingham's repository
(`nottingham-repository.worktribe.com`) sits behind a Cloudflare bot check
that blocked every automated route we tried, including a real Playwright
browser session over an institutional network. [references.md](/replication/references#the-two-that-automation-couldnt-get-hls-and-hlw)
has the full history. We once took two Wayback Machine snapshots
(`.../preview/829483/bubble_dates.pdf` and
`.../preview/4418036/bubble_dates.pdf`) for recovered copies. The working
copies now in the paper library follow the same filename convention, but they were
supplied directly and not retrieved by any tool in this project. PDC's
working paper is open (MPRA 92074), and so is KS (arXiv:2110.04500).

We checked every formula quoted below against the typeset math. These are
HLS Models 1-4 and BICj, HLW's two-step Models and BICj, and PDC's Step 1
and Step 2 estimators. We rendered the PDF pages to PNG (`pymupdf`,
`matrix=fitz.Matrix(2.5,2.5)`) and read them directly, because
`pdftotext -layout` scrambles stacked multi-line formulas and subscripts in
all four PDFs. HLS's DGP definition and its SSR order-of-magnitude table on
p.4-5, for example, come out as out-of-order fragments. The Model and BIC
formulas on p.7-8 are not scrambled and match the OCR text once laid out,
but we checked that by rendering them. We also
checked HLW's Table 1 (BIC model-selection frequencies) against the image.
In the plain-text dump of that table the row labels `A`, `B`, `C`... look
shifted relative to their `j=1`/`j=2` data rows, and a first read could
easily pair them wrongly. The rendered image gives the correct pairing,
which differs from the naive text-order reading.

### What it is

All three lines of work replace PSY's **threshold-crossing** dating rule
(find where the recursive BSADF/ADF t-statistic sequence first crosses, or
re-crosses, a critical value) with a **model-based SSR-minimisation + BIC**
rule. They differ in how they handle *multiple* bubbles in one series.

#### 1. HLS (2017), single-bubble SSR+BIC dating

Defines four candidate regime-structure models for a series `y_t`. Each is
a piecewise OLS regression of `Δy_t` on regime-indicator dummies and
dummy-interacted `y_{t-1}` (own notation, checked against the page image,
p.7):

```
Model 1: Δyt = μ1·Dt(τ1,1)      + δ1·Dt(τ1,1)yt-1                                    + v1t   (unit root → bubble to sample end)
Model 2: Δyt = μ1·Dt(τ1,τ2)     + δ1·Dt(τ1,τ2)yt-1                                    + v2t   (unit root → bubble → unit root)
Model 3: Δyt = μ1·Dt(τ1,τ2) + μ2·Dt(τ2,1)  + δ1·Dt(τ1,τ2)yt-1 + δ2·Dt(τ2,1)yt-1      + v3t   (unit root → bubble → collapse to end)
Model 4: Δyt = μ1·Dt(τ1,τ2) + μ2·Dt(τ2,τ3) + δ1·Dt(τ1,τ2)yt-1 + δ2·Dt(τ2,τ3)yt-1     + v4t   (unit root → bubble → collapse → unit root)
```

with `Dt(a,b) = 1(⌊aT⌋ < t ≤ ⌊bT⌋)`. For each model, the break-fraction(s)
`(τ̂1, τ̂2, τ̂3)` are the values that jointly minimise the model's residual
sum of squares (`SSRj`) over all candidate dates satisfying ordering/sign
constraints (e.g. `y⌊τ2T⌋ > y⌊τ1T⌋` to force the bubble phase to be upward).
Theorem 1 proves `⌊τ̂iT⌋ − ⌊τi,0T⌋ →p 0` for each *correctly paired*
DGP/Model. Under a fixed-magnitude bubble the estimator is therefore
consistent for the exact date, and not only for the break fraction.

Model selection across the four candidates uses a BIC whose penalty is the
number of fitted coefficients *plus* the number of estimated break dates,
times `ln(T)` (checked against the image, p.8):

```
BIC1 = T·ln{T⁻¹SSR1(τ̂1,1)}        + (2+1)ln(T)
BIC2 = T·ln{T⁻¹SSR2(τ̂1,τ̂2)}       + (2+2)ln(T)
BIC3 = T·ln{T⁻¹SSR3(τ̂1,τ̂2,1)}     + (4+2)ln(T)
BIC4 = T·ln{T⁻¹SSR4(τ̂1,τ̂2,τ̂3)}    + (4+3)ln(T)
jopt = argmin_j BICj
```

In practice (§5) HLS impose minimum regime durations (`τ1 ≥ s`,
`τ2−τ1 ≥ s`, `τ3−τ2 ≥ s/2`, with `s = 0.1` in the T=200 simulations and
`s = 0.05` in the T=389 empirical application). Otherwise the method is a
**brute-force grid search** over 1, 2, or 3 breakpoints depending on the
model. The paper describes no dynamic-programming (Bai-Perron-style)
speedup.

#### 2. HLW (2020), two-step extension to multiple bubbles

HLW address PSY's known late end-date bias directly. HLS's four-model set
generalises badly to N bubbles, because the model set grows
combinatorially. HLW therefore propose a **two-step** procedure instead
(checked against the image, p.9; the formulas match the `pdftotext` text
closely):

- **Step 1**: run PSY's GSADF/BSADF detection and dating exactly as it
  stands, to get preliminary start and end fractions `τ̂P SY_j1, τ̂P SY_j2`
  for each of the `N̂` detected bubbles. Use them to carve the sample into
  `N̂` disjoint sub-sample "date windows" `[sj, ej]`. Windows are split at
  the midpoint between consecutive PSY-detected regimes, with a rule that
  nudges the split so a window always starts inside a fitted
  post-explosive (unit-root) regime and never mid-bubble.
- **Step 2**: apply HLS's SSR+BIC Model 1-4 procedure *independently within
  each date window*. For every window but the last, restrict to Models 2
  and 4, since a window boundary is by construction a unit-root point and
  not a sample end.

HLW present this as bolting the HLS refinement onto PSY's own output, and
not as a replacement detection method: "we propose a dating methodology
based on minimum sum of squared residual estimators and BIC model
selection, but using prior information gleaned from the PSY dating
procedure as a means of reducing the dimensionality" (HLW, §1).

#### 3. PDC (2021, journal) / KS (2023, journal), sequential sample-splitting

A structurally different and computationally much cheaper approach to the
*same* SSR-minimisation idea, for a single bubble episode (3- or 4-regime
model). HLW themselves suggest extending it to the multiple-bubble case by
running it inside each PSY date window instead of using HLS's joint fit.

PDC's model is unit-root → explosive → stationary-collapse (3 regimes, 2
breakpoints `τ1_0 < τ2_0`). HLS minimise a 2-breakpoint SSR surface
jointly. PDC instead show, through a stochastic-order argument (their
"Example 3" and Lemmas A.2-A.4), that under the bubble DGP the **collapse
date is always identified first**: the SSR drop at the collapse breakpoint
dominates, in stochastic order, the drop at the origination breakpoint.
The two breaks can therefore be estimated **sequentially** and not jointly:

```
Step 1: τ̂2 = argmin_{τ∈(0,1)} RSS2,T(τ),  RSS2,T(τ) = Σ_{t≤⌊τT⌋}(yt − β̂x(τ)yt-1)² + Σ_{t>⌊τT⌋}(yt − β̂3(τ)yt-1)²
        where β̂x(τ) = Σ_{t≤⌊τT⌋} yt·yt-1 / Σ_{t≤⌊τT⌋} yt-1²   (no-intercept AR(1) slope, full-sample split at τ)
Step 2: on the left subsample [1, τ̂2T] only, repeat the same one-break RSS minimisation to get τ̂1.
```

(checked against the image, PDC p.9). This is a **no-intercept model with a
single AR(1) coefficient per regime**, unlike HLS's intercept+AR(1)
dummies. Each `β̂(τ)` is a closed-form ratio of two prefix sums
(`Σy_t y_{t-1}`, `Σy_{t-1}²`), so the whole `RSS(τ)` curve over all
candidate `τ` is computable in `O(T)` via cumulative sums. There is no
joint grid search and no model-selection BIC step: PDC's algorithm always
assumes the 3-regime structure holds and estimates its two breaks one at a
time.

KS (2023) extend PDC's 3-regime model to 4 regimes by adding a final
unit-root "recovery" regime after the stationary collapse, and reuse PDC's
sequential logic for the extra breakpoint. They also contrast their own
cost with HLS's: "we perform the three SSR minimization with one break each
with O(T) computations, while Harvey et al. (2017) requires minimizing the
three break model over all possible combinations of these breaks" (KS, §1,
p.3). KS thus read HLS's Model 4 grid search as an un-sped-up
multi-dimensional combinatorial search, which agrees with our reading of
HLS §5 above. KS also say their method is meant to slot into HLW's
per-window second step: "Harvey et al. (2020) proposed to initially
identify the bubble regimes based on [the] PSY approach ... one can use our
approach in the second step" (KS, §1).

### Exact numbers reproduced

**HLS Table 1** (Nasdaq composite real price index, PWY's own 1973:2-2005:6
series). We checked the table against the image and against `pdftotext`;
the two matched exactly, with none of the scrambling that affects HLW's
Table 1:

| Sample | PSY test | PSY start | PSY end | BICopt model | BICopt start | BICopt end |
|---|---|---|---|---|---|---|
| 1973:2-2005:6 (full) | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2-2000:9 (pseudo-real-time) | 3.07*** | 1998:11 | 2000:9 | 1 | 2000:1 | 2000:9 |
| 1973:2-2000:10 | 3.07*** | 1998:11 | 2000:10 | 1 | 2000:1 | 2000:10 |
| 1973:2-2000:11 | 3.07*** | 1998:11 | 2000:11 | 1 | 1999:12 | 2000:11 |
| 1973:2-2000:12 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2-2001:1 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |

This is the single concrete quantitative comparison HLS publish between the
two dating rules. Both agree on the 1998:11 start, but PSY's end date
(2000:12) is 3 months later than BICopt's (2000:9). In this worked example,
then, PSY's original dating strategy is late on the *end* date only.

**What could not be reduced to an exact number.** HLS's Monte Carlo
dating-accuracy comparison (§6, Figures 1-3) is reported **only as plots**
of frequency against bubble magnitude. That is the head-to-head "BICopt vs
PSY, % correct within k observations" evidence, and it is not tabulated.
The text states the qualitative conclusion ("BICopt out-performs PSY... in
finite samples, particularly with respect to the bubble's end date") but
gives no percentage or RMSE in prose either. HLW §4 is in the same
position. Its Figures 1-6 are histograms of PSY and BIC start/end date
estimates, and the prose describes patterns ("PSY estimates ... typically
fall somewhat later than the true date", "BIC estimated end dates equal the
true end date in almost every replication") with no numeric table. This
mirrors the [SBZ finding](/replication/volatility-robustness#why-this-cant-be-a-bit-exact-numeric-cross-check-test):
a Monte Carlo comparison reported as a figure and not as a table cannot be
cross-checked bit for bit, so we report it here as a qualitative claim
only.

**HLW Table 1** (BIC model-selection frequencies across 6 DGPs A-F). We
checked it against the image. In the `pdftotext -layout` dump the row
labels (`A`,`B`,`C`,`D`,`E`,`F`) look offset from their `j=1/j=2/j=3` data
rows, and the values below come from the rendered image and not from the
text dump:

| DGP | Regime | Model 1 | Model 2 | Model 3 | Model 4 | True model |
|---|---|---|---|---|---|---|
| A | j=1 | 0.007 | 0.133 | 0.128 | **0.732** | 4 |
| A | j=2 | 0.033 | 0.062 | 0.451 | **0.454** | 4 |
| B | j=1 | 0.013 | 0.127 | 0.134 | **0.726** | 4 |
| B | j=2 | 0.026 | 0.060 | **0.633** | 0.281 | 3 |
| C | j=1 | 0.049 | **0.765** | 0.080 | 0.106 | 2 |
| C | j=2 | 0.089 | 0.176 | 0.299 | **0.436** | 4 |
| D | j=1 | 0.285 | **0.684** | 0.009 | 0.022 | 2 |
| D | j=2 | **0.627** | 0.343 | 0.016 | 0.014 | 1 |
| E | j=1 | 0.006 | 0.027 | 0.001 | **0.966** | 4 |
| E | j=2 | 0.156 | 0.039 | 0.022 | **0.782** | 4 |
| E | j=3 | **0.642** | 0.015 | 0.002 | 0.342 | 1 |
| F | j=1 | 0.002 | 0.088 | 0.027 | **0.883** | 4 |
| F | j=2 | 0.018 | 0.095 | 0.252 | **0.635** | 4 |
| F | j=3 | 0.014 | 0.061 | **0.524** | 0.401 | 3 |

BIC picks the true model most often in every row. The weakest cases are
DGP A/j=2 and DGP C/j=2, where the correct-model rate is about 44-45%. In
both, a post-collapse reversion to a unit root close to the end of the
fitted window is easily missed, as the paper's own caveat text says. This
number bears on dating, since it is the frequency of correct model
identification, but it is **not** dating accuracy. The paper is explicit
that choosing the wrong model (for example Model 3 over Model 4) still
typically gives accurate `τ̂1`/`τ̂2` estimates.

**KS empirical application** (§6, stated in prose and not tabulated).
NASDAQ Composite, Jan 1985-Aug 2013 monthly: BIC4 (3387.727) beats BIC3
(3404.268) and BIC2 (3409.296), so the 4-regime model is selected. The
collapse is dated Feb 2000, origination (from the left-subsample re-split)
Aug 1998, and recovery (right subsample) Sep 2001. US house price index
(FHFA, real), Jan 1991-Dec 2012: BIC4 (-2918.223) is again selected, with
collapse Nov 2006, origination Sep 1997 and recovery May 2011.

**KS Monte Carlo** (§5) is reported as histograms (Figures 1-4) and not as
tables. The only numbers in the prose are approximate: "the frequency of
selecting the true break date is not very high (it is around 30% for T=400
and 65% for T=800)... increases to almost 100%..."; "close to 100%";
"approximately 75% and 100% for T=400 and 800, respectively". We quote
these in the paper's own words and have not computed them independently, so
read them as orders of magnitude and not as exact figures.

### Cost/feasibility note for exuber

exuber's `datestamp()` (`exuber/R/radf-methods.R`, `datestamp.radf_obj`,
plus helpers `stamp()`/`stamp_to_index()`/`add_peak()` in the same file) is
**pure post-processing on statistics that `radf()` has already computed**.
It takes the recursive `bsadf`/`badf` sequences (from
`augment_join(object, cv)`), compares them with the simulated or
wild-bootstrap critical values at a chosen significance level
(`tstat > crit`), and finds contiguous runs of exceedance (`stamp()` is
literally `which(diff(x) != 1)` on the indices where the inequality holds).
No new regression is fit anywhere in this path. It is threshold-crossing on
numbers that already exist.

Describing a `datestamp(method = "bic")` as a "contained, high-value
addition" undersells the work for the HLS/HLW route, though it is roughly
right for the PDC/KS route. The two need different amounts of work:

1. **HLS-style BIC dating is not a post-processing step.** We built it
   anyway (2026-08-10). Per bubble episode it required (a) fitting up to 4
   new OLS regression *specifications* (regime-dummy-interacted AR(1)
   models with intercept), none of which match exuber's existing `radf()`
   machinery; (b) a grid search over 1-3 breakpoints jointly per model; and
   (c) a BIC comparison across the 4 fitted models. We worried that the
   `O(T^m)` scaling would force a Bai-Perron-style dynamic-programming
   rewrite, and it did not. The four models' regime dummies never overlap,
   so any candidate partition's total SSR decomposes exactly into
   independent per-segment closed-form OLS fits. Each is an `O(1)` lookup
   into precomputed cumulative sums, the same trick that `dating_pdc()`'s
   (differently specified) breakpoint search uses. The grid search
   therefore stays fast (about 2 seconds at `T=400` for the full 4-model
   search), even though in absolute terms it is still an `O(T)`, `O(T^2)`
   or `O(T^3)` search for the four models respectively; see
   [Implementation (HLS route)](#implementation-hls-route) below. HLW's
   step 1 (splitting into date windows) can reuse exuber's existing
   `radf()`/`datestamp()` output almost as it stands, since it is exactly
   PSY's detected start and end dates, which `datestamp()` already returns.
   Step 2 (applying `dating_hls()` inside each window) is a thin wrapper
   around code that had already shipped, and not the "100% new estimation
   code" we first expected. The wrapper itself (window construction and the
   Models-2-and-4-only restriction for non-final windows) was not built at
   the time of this note. It came later as `dating_hlw()`; see
   [Implementation (HLW route)](#implementation-hlw-route).

2. **PDC/KS-style sequential sample-splitting is a meaningfully smaller,
   "contained" addition.** Each breakpoint estimate is a single `O(T)` scan
   of a *closed-form* ratio of cumulative sums, with no intercept, no joint
   multi-dimensional grid and no BIC loop across 4 specifications. It sits
   in the same complexity class as
   [STADF's `gls_dfstat_grid()`](/replication/volatility-robustness#formulas-implemented-verified-against-rendered-pdf-pages-not-ocr-text),
   so it is plausibly implementable in R via `cumsum()`/`outer()` without
   new C++. The price is a narrower model. PDC/KS take the regime count (3
   or 4) as known, where HLS select among 4 alternatives by BIC. The PDC/KS
   route therefore cannot tell a bubble that collapses from one that does
   not, or from one still running at the sample end, and that distinction
   is what HLS/HLW's model-selection step buys in practice. We recommended
   this route as the next step at the time, and it has since shipped as
   `dating_pdc()`.

3. **Either way this is a new statistical layer and not a flag on the
   existing rule.** Both routes need new regime-dummy or no-intercept
   AR(1) fitting code, which editing `stamp()`/`datestamp()` cannot supply.
   Both need an explicit minimum-regime-duration trimming parameter (`s` in
   HLS/HLW, the `ψ`-style hyperparameters in PDC/KS), analogous to exuber's
   `psy_ds()`/`minw` but serving a different purpose. The HLS/HLW route
   also needs a BIC computation and a 1-3 dimensional grid search with no
   exuber analogue to build on. A `datestamp(..., method = "bic")` would be
   a self-contained R-only feature (PDC/KS variant) or a two-tier R feature
   that reuses existing `datestamp()` output for step 1 (HLW variant).
   Neither is a small diff on top of `stamp()`/`add_peak()`. Each is new
   estimation code, comparable in size to the
   [SBZ effort](/replication/volatility-robustness#implementation) for HLS/HLW and
   somewhat smaller for PDC/KS, and nowhere near the near-zero cost of
   STADF's formula reuse.

### Implementation (HLS route)

Shipped as `dating_hls(data, trim = 0.05)` in `exuber/R/dating_hls.R`,
tested in `exuber/tests/testthat/test-hls.R`. It fits all four of HLS's
regime-dummy models, each by exact SSR minimisation over its candidate
breakpoint(s), subject to the paper's minimum-regime-duration trim and
upward-bubble sign constraint (`y_{tau2} > y_{tau1}`). It selects among
them with the paper's BIC formula (`n*log(SSR/n) + df*log(n)`, with `df`
in `{3,4,6,7}` for Models 1-4 respectively). The function returns the
selected model number and its origination, collapse and recovery dates
(`NA` for any the selected model does not have), plus every candidate
model's BIC value so you can see how close the selection was.

**The key implementation move is a change of algorithm.** HLS's paper
describes a brute-force grid search with no speedup, and the
cost/feasibility note above worried that it would scale as `O(T^m)` for an
`m`-breakpoint model, with a large constant. Because the four models'
regime dummies never overlap, any candidate partition's SSR is exactly the
*sum* of independent per-segment OLS fits. A segment with no active dummy
has zero fitted parameters (`SSR = sum(Delta y_t^2)`), and a segment with an
active dummy is a plain intercept+slope fit. Both reduce to a closed-form
ratio of cumulative sums (`Sx, Sxx, Sz, Szz, Sxz`), so every candidate
breakpoint (or pair or triple of them) evaluates in `O(1)` given
precomputed prefix sums, with no repeated `lm()` calls. `dating_pdc()`'s
breakpoint search uses the same trick with a different (no-intercept)
specification. The search stays fast: `dating_hls()` returns in well under
2 seconds at `T=400`, HLW's own sample size, in plain R with no C++.

**Validation.** The formula checks are exact. `hls_segment_ssr()` matches a
brute-force `lm()` fit's SSR (tolerance `1e-8`) for arbitrary segments. The
*full joint 3-breakpoint grid search* for Model 4 also matches an
exhaustive brute-force nested-`lm()` search bit for bit on a small
synthetic series, so the check covers more than the per-segment formula in
isolation. The Monte Carlo runs use synthetic regime DGPs built on a large
positive base level (`100 + cumsum(rnorm(...))`), so that the explosive
signal is a large departure from noise on the absolute scale:

- **A true Model 4 DGP** (unit root, bubble, mean-reverting collapse,
  unit-root recovery): BIC selects Model 3 or 4 in 100% of 30 replications,
  so it always identifies that a distinct collapse regime exists. The split
  between the two is 80%/20%. This is consistent with HLW's finding that
  telling a final recovery regime from continued collapse is the hardest
  case: their Table 1 reports a correct-model rate as low as about 44-45%
  in their weakest DGP. The origination date's mean absolute bias is about
  0 observations, and the collapse date's mean absolute bias (among Model
  3/4 selections) is about 1 observation.
- **A true Model 2 DGP** (bubble fully reverting to a unit root, no
  distinct collapse regime): BIC selects Model 2 in 100% of 30
  replications, with origination bias exactly 0 in every replication.
- **A true Model 1 DGP** (bubble ongoing at the sample end): BIC selects
  Model 1 in 100% of 30 replications, with origination bias exactly 0 in
  every replication.
- **Pure `H0`** (no bubble at all): BIC never selects the most complex
  model (0% Model 4 across 30 replications). It splits between the three
  simpler ones (30% Model 1, 57% Model 2, 13% Model 3), so complex models
  do not run away under pure noise. HLS's method has no "no bubble" option,
  since it is designed to run downstream of a PSY-detected episode, so this
  is a curiosity check and not a core requirement.

**A DGP construction pitfall, fixed during this validation (not a code
bug).** An early version of the synthetic bubble DGP built the explosive
regime as `unit1[n1] * c^i`, where `unit1[n1]` is itself a small, mean-zero
random-walk value that can even be negative. That gave a weak and sometimes
backwards "explosive" signal relative to the cumulative noise, and under it
BIC rationally preferred the simplest model (Model 1) about 80% of the time
whatever the true DGP. We diagnosed this by comparing the selected models'
actual SSR and BIC values and not only the selection frequencies, which
traced the problem to DGP signal strength and not to the search or the BIC
code. Rebuilding the DGP on a large positive base level fixed it. The same
construction was already working in the Monte Carlo checks for
`quantile_test()`.

**A sign-constraint gap, found and fixed shortly after the first release
while we were reading HLW (2020) for the multi-bubble extension.** HLW
restate HLS's Models 1/3/4 with sign constraints that our first version had
abbreviated. Model 1 requires `y_T > y_{tau1}`, so the series must end above
where the bubble started. Models 3 and 4 require the fitted peak `y_{tau2}`
to exceed not only the bubble's starting level but also the level the
series reaches after the fitted collapse regime (`y_T` for Model 3,
`y_{tau3}` for Model 4); otherwise `y_{tau2}` is not a peak. The initial
`dating_hls()` checked only `y_{tau2} > y_{tau1}` and so missed the second
half of each constraint. We fixed all three affected search functions and
updated the brute-force validation to match. Re-running the Monte Carlo
checks improved accuracy further: Model 1 and Model 2 selection went from
97% to 100% correct, with bias dropping to exactly 0 in every replication.
The numbers above already reflect the fix. Replication script:
[replication/dating-and-root-inference/radf_hls_validation.R](#script-radf_hls_validation).

### Implementation (HLW route)

Shipped as `dating_hlw(data, cv = NULL, minw = NULL, trim = 0.1, min_duration = NULL, nboot = 199L, seed = NULL)`
in `exuber/R/dating_hlw.R`, tested in `exuber/tests/testthat/test-hlw.R`. It
is a two-step wrapper around `dating_hls()`, following HLW's paper:

1. **Step 1**: run PSY's detection and dating (`radf()`/`datestamp()`, with
   `min_duration` defaulting to `psy_ds(n)`, HLW's own `ln(T)`
   minimum-episode-length rule) to get a preliminary start and end position
   for each detected explosive episode.
2. **Step 2**: carve the sample into disjoint date windows using HLW's
   formula. The end of window `j` is the midpoint between that episode's
   PSY-detected end and the next episode's PSY-detected start
   (`e_j = tau2_psy[j] + floor((tau1_psy[j+1] - tau2_psy[j])/2)`), and the
   last window runs to the sample end. Then fit each window with
   `hls_fit_series()`, a shared internal helper factored out of
   `dating_hls()` for this reuse. Every window but the last is restricted to
   Models 2 and 4, since a window boundary is by construction a unit-root or
   collapse point and not a sample end. After fitting window `j`, HLW's
   sequential-adjustment rule sets the *next* window's start to the first
   observation of the just-fitted post-explosive regime
   (`s_{j+1} = s_j + tau2_local` for a Model 2 fit, `s_j + tau3_local` for
   Model 4). This prevents a later window from starting mid-bubble, a
   failure mode that HLW's paper discusses explicitly.

**Two bugs found and fixed while building this** (neither was a design
choice):

1. `datestamp()` raises a hard error ("Cannot reject H0 at the 5%
   significance level"), not a warning, when no series has any detected
   episode. An initial version of `dating_hlw()` caught only `warning`
   conditions from this call, so a series with no bubble crashed instead of
   returning the intended empty result. It now catches `error` as well.
2. HLW's restatement of HLS (2017)'s Models 1/3/4 sign constraints turned
   out to be more complete than what the just-shipped `dating_hls()`
   checked; see the sign-constraint fix above. We fixed it in
   `dating_hls()`, which it affects directly, and `dating_hlw()` inherits
   the fix by reuse.

**Validation** (structural checks, Monte Carlo, and a strong equivalence
check). We checked `hlw_local_to_global()`'s index arithmetic directly. On a
synthetic two-bubble DGP (20 reps), PSY step-1 detection found exactly 2
windows in 13 of 20 reps. In the rest it fragmented into 3 or more spurious
sub-windows. HLW's paper discusses this failure mode and proposes a
run-joining heuristic for it, which we did not implement; see "Not
implemented" below. Among the 13 clean reps, origination and collapse date
bias for *both* bubbles was exactly 0 in every replication. Windows were
correctly ordered and non-overlapping in all 20 reps. Under a pure `H0` null
(no bubble), `dating_hlw()` never errors and returned 0 detected windows in
all 20 reps tried, so there were no false positives. On a single clean
bubble episode (15 reps), the wrapper's *final* window matched standalone
`dating_hls()` on the whole series exactly (model, origination and collapse
all identical) in every rep that had a final window. This is a strong
structural check, since HLW's paper states that the two-step procedure
reduces to plain HLS when there is only one episode. Replication script:
[replication/dating-and-root-inference/radf_hlw_validation.R](#script-radf_hlw_validation).

**Not implemented.** HLW's run-joining heuristic for step-1 fragmentation
("if up to 3 non-rejections are surrounded on either side by an explosive
regime of length `ln(T)`, treat them as a single episode"). `dating_hlw()`
uses `datestamp()`'s regimes as detected, un-joined, so a single true bubble
can occasionally surface as several windows when PSY's step-1 detection
fragments it. We quantified this above: 13/20 clean vs. 7/20 fragmented on
the two-bubble DGP, and 12/15 vs. 3/15 on the single-bubble DGP. The
fragmentation comes from noise in PSY's step-1 detection, and not from the
window construction or the per-window fitting, both of which validated
exactly. We scoped the heuristic out as its own small, well-defined
follow-on.

### Implementation (PDC/KS route)

Shipped as `dating_pdc(data, regimes = 3L, trim = 0.05)` in
`exuber/R/dating_pdc.R`, tested in `exuber/tests/testthat/test-pdc.R`. The
internal helper `pdc_find_break(y, trim)` implements the single-breakpoint
no-intercept AR(1) RSS minimiser exactly as specified above (§"3. PDC/KS").
It computes the closed-form `β̂(τ)` from prefix sums of `Σy_t y_{t-1}` and
`Σy_{t-1}²`, so the whole `RSS(τ)` curve costs `O(T)`. `dating_pdc()` calls
it sequentially following PDC's Step 1 and Step 2 (collapse first, then
origination on the left subsample), and once more on the right subsample
for KS's 4-regime recovery date when `regimes = 4`.

**Validation methodology.** This function was already present in the
repository, uncommitted and with no test coverage, so we validated it from
scratch before trusting it:

1. **Formula-exact check.** `pdc_find_break()`'s closed-form break index and
   RSS match a brute-force scan that fits `lm(y[t] ~ y[t-1] - 1)`
   separately on every candidate left/right split, bit for bit
   (`tolerance = 1e-8`). The cumulative-sum algebra is therefore both fast
   and correct.
2. **Consistency in the low-noise, long-series, strong-effect limit.** On a
   synthetic 3-regime series (unit root → explosive → collapse) and a
   4-regime series (... → recovery), `dating_pdc()` recovers the true break
   dates to within 1-2 observations. We use the same technique elsewhere in
   this project (SBZ, common-bubble) to separate "the estimator is
   asymptotically correct" from "the estimator happens to pass one
   finite-sample tolerance check". A check that holds only at moderate
   sample sizes cannot tell those apart, whereas convergence as noise
   shrinks and sample size grows can.
3. **A DGP-construction pitfall worth recording.** Our first 4-regime
   consistency check used a *deterministic* exponential decay
   (`target + (peak - target) * exp(-k*t)`) for the collapse regime. That
   decay saturates numerically well before the regime ends, and its noisy
   "flat" tail near `target` is then statistically indistinguishable from
   the random walk of the following recovery regime, since both look like
   small increments near a level. `pdc_find_break()` correctly finds a
   break, but the wrong one, *inside* the collapse regime and not at its
   boundary with the recovery regime. The bug was in the test and not in the
   estimator: `pdc_find_break()` fits a no-intercept AR(1), and a
   deterministic trend plus noise is not one. We fixed it by generating the
   collapse regime as a stationary AR(1) recursion
   (`y_t = ρ·y_{t-1} + ε_t`, `ρ = 0.5`). The dynamics are then homogeneous
   throughout the regime, and cleanly distinguishable from the `ρ = 1`
   regimes on either side.
4. **A plain account of finite-sample accuracy.** A separate moderate-`T`
   check (`T = 350`, 30 seeds, §4 of `radf_pdc_validation.R`) reports the
   exact-date recovery rate (origination 3.3%, collapse 0%) and mean
   absolute error (origination 5.5, collapse 1.0), without asserting tight
   recovery. This matches KS's own Monte Carlo (§5, quoted above), which
   reports only about 30% exact-date recovery at `T = 400`, rising to about
   65% at `T = 800`. A synthetic test that demanded tight accuracy at small
   `T` would either be cherry-picked or silently contradict what the source
   paper reports about its method. We record this in the same way as the
   [root-CI coverage finding](#root-inference) below.

**Result**: 13/13 assertions pass (`devtools::test()`, full suite: 306
passed, 0 failed). We found no implementation bug. The only fix was to the
test's synthetic data-generating process, and `dating_pdc.R` itself did not
change.

Replication script:
[replication/dating-and-root-inference/radf_pdc_validation.R](#script-radf_pdc_validation).

---

## Root inference

**Status: done.** Guo, Sun & Wang's normal-t CI with doubling time, the
Phillips-Magdalinos Cauchy CI, and per-episode datestamp integration are
all implemented in `exuber/R/rootstamp.R` as a single S3 generic,
`rootstamp()`. The `default` method handles a single sub-sample and the
`radf_obj` method handles every `datestamp()` episode at once. Tests are in
`exuber/tests/testthat/test-rootstamp.R`. **Update (2026-08-18):** the code
first shipped as three separate functions (`explosive_root()`, `root_ci()`,
`root_ci_datestamp()`). We consolidated them into `rootstamp()` before
release, because the three-call handoff and the backwards-reading name
`root_ci_datestamp()` were friction worth fixing early, and the old API was
not stable enough to be worth preserving. Mentions of the three old names
below are left as they were, since they describe what was true at each
dated entry. Read them as `rootstamp()` in current terms.

### Source

- Phillips, P. C. B., & Magdalinos, T. (2007). Limit theory for moderate
  deviations from a unit root. Journal of Econometrics, 136(1), 115-130.
  Working paper: Cowles Foundation DP 1471 (July 2004), **open** at
  `cowles.yale.edu/sites/default/files/2022-08/d1471.pdf`. It is not
  paywalled; the original search simply missed it. We verified Theorem 4.3
  by rendering the page to PNG on 2026-08-09.
- Guo, G., Sun, Y., & Wang, S. (2019). Testing for moderate explosiveness.
  The Econometrics Journal, 22(3), 279-303. The paper is paywalled and we
  first could not access it; we later recovered it through institutional
  access.
  We implemented `root_ci()` from a secondary restatement before finding it
  (see below), so it is worth cross-checking against the primary source now
  that it is available.

The first implementation drew on a secondary review that restates both
papers' results: Skrobotov, A. (2023), "Testing for explosive bubbles: a
review", open on arXiv (`arxiv.org/pdf/2207.08249`). It is the same review
this whole project started from.

### What's implemented, and why only part of it

The review states two distinct results:

1. **Phillips-Magdalinos (2007a)**: for a mildly explosive AR(1)
   \eqn{\rho_T = 1 + c/k_T}, a specific normalization of
   \eqn{(\hat\rho_T - \rho_T)} converges to a **standard Cauchy**
   distribution, and this holds under non-Gaussian errors too. The review
   gives the resulting CI as
   \eqn{\hat\rho_T \pm \frac{\sqrt{\hat\rho_T^{2n}-1}}{\hat\rho_T^n} C_\xi}.
2. **Guo, Sun & Wang (2019)**: generalizing to allow a drift term, they show
   that the **ordinary regression t-statistic** for \eqn{\rho_T} is
   asymptotically **standard normal** under i.i.d. errors (Student's-t/HAR
   under dependence). This is practically simpler, because it does not need
   to know \eqn{c}, \eqn{k_T}, or the exact rate.

At first we implemented only (2). Result (1)'s exact formula came through
`pdftotext` with its subscripts and exponents visibly mangled. The review
is a secondary source restating a third paper's theorem, so we had no
primary-source PDF to render as an image and check. Shipping the formula
with an unverified exponent risked a silently wrong confidence interval,
which is worse than no interval. Result (2) survived as unambiguous prose
("the t-statistic... is asymptotically standard normal"), which was enough
to implement it correctly at low risk: `root_ci()` is just
`rho_hat +/- qnorm(...) * se`, with the standard error from an ordinary
no-intercept OLS fit (`explosive_root()`). The one unusual point is the
justification for treating this as normal despite `rho > 1`, which is
Guo/Sun/Wang's CLT and not the classical stationary-AR one.

### Exact number verified

The review's footnote 17 gives the two-sided Cauchy percentiles used in
result (1): C_0.10 = 6.315 (≈6.314), C_0.05 = 12.7, C_0.01 = 63.65674. The
standard Cauchy distribution is identical to Student's *t* with 1 degree of
freedom, a mathematical fact independent of whichever paper states it. The
percentiles can therefore be checked exactly, with no simulation and no
primary-source risk: `qt(0.95, df = 1) = 6.313752`,
`qt(0.975, df = 1) = 12.7062`, `qt(0.995, df = 1) = 63.65674`. These match
the footnote to the precision given and are tested in `test-rootstamp.R`.
The check validates the quoted percentile table only. It says nothing about
the *exponent* in result (1)'s CI formula, which is why that formula did
not ship immediately.

### What `root_ci()` gives you

Given a sub-sample (for example an explosive episode already identified by
`datestamp()` and converted to row positions), `explosive_root()` fits the
no-intercept AR(1) OLS regression that Phillips-Magdalinos's model
specifies. `root_ci()` then reports:
- `rho`, `rho_ci`: point estimate and Wald interval,
- `doubling_time`, `doubling_time_ci`: `log(2)/log(rho)` (periods for the
  bubble to double at the estimated rate) with its own interval, obtained
  by transforming the `rho` interval's endpoints (doubling time decreases
  in `rho`, so the bounds flip).

### Empirical coverage

A 500-replication simulation at a modest sample size (explosive episode of
149 observations, ρ = 1.03) gave empirical coverage of about 90% for the
nominal 95% CI. This is finite-sample undercoverage and not a bug. We
verified that the formula reduces to the textbook Wald interval, and the
CLT it relies on is a T → ∞ result for an estimator known in the
literature to converge slowly. The coverage test in `test-rootstamp.R`
checks against a deliberately loose bound (>80%) for this reason, and does
not assert the nominal rate.

### Independent validation (2026-08-09)

We used different parameters and a different seed from `test-rootstamp.R`
throughout (that suite uses `rho=1.03, n=150`; here `rho=1.05, n=200`), to
check that the finding is not an artifact of one parameterization.

**Point estimate.** On a fresh simulated explosive AR(1) (`seed=8675309`),
`explosive_root()` recovered `rho_hat = 1.0500` against `rho_true = 1.0500`,
matching to 4 decimal places, with a 95% CI so tight that it is
indistinguishable from the point estimate at that precision. This is not a
bug. It is the expected "super-consistency" of explosive-root estimation:
the regressor `y_{t-1}` itself grows geometrically, so
`sxx = Σy_{t-1}²` (the denominator of the OLS standard error) explodes at
rate `ρ^{2n}`. By `n=200` with `ρ=1.05` that denominator is astronomically
large, and the standard error collapses far faster than in the unit-root or
stationary case. This agrees with the theory the CI is built on (Guo, Sun &
Wang's CLT for an explosive root).

**Coverage** (the more informative check): 800 independent replications at
`rho=1.05, n=200, seed=24601`:

```r
covered <- replicate(800, {
  y <- numeric(200); e <- rnorm(200)
  for (t in 2:200) y[t] <- 1.05 * y[t-1] + e[t]
  ci <- rootstamp(y)
  ci$rho_ci[1] <= 1.05 && 1.05 <= ci$rho_ci[2]
})
mean(covered)
```

Result: **94.6%** coverage of the nominal 95% interval, a marked
improvement on the roughly 90% that the package's own test observes at
`rho=1.03, n=150`. The direction is the expected one. A larger `n` and a
root further from 1 (`1.05` vs `1.03`) both bring the estimator into its
asymptotic (`T→∞`) regime faster, so coverage should approach nominal as
either increases, and this run shows that it does. The result supports the
known finite-sample undercoverage disclosure. The undercoverage shrinks in
the direction the theory predicts, which also indicates that the CLT
justification is the right one and the match is not a coincidence.

**Full existing suite**: `test-rootstamp.R`, **8 passed, 0 failed** (1
skipped, CRAN-only).

**Conclusion**: we found no issues. The point estimate and CI behave as the
asymptotic theory predicts, including the direction of convergence as
`(n, ρ)` move further from the boundary.

### Update: primary source found, exponent now verified

The Phillips-Magdalinos (2007) primary source turned out to be open after
all; see Source above. `pdftotext` mangles the formula in the same way it
mangled SBZ's Table 1, with exponents and fraction bars merging into runs
like "nnn". We therefore rendered page 14 to a PNG with PyMuPDF (the method
used for STADF/SBZ) and read it directly. **Theorem 4.3** (their eq. 26),
for the moderate-deviations model \eqn{\rho_n = 1 + c/n^\alpha},
\eqn{c>0}, \eqn{\alpha \in (0,1)}, reads:

\deqn{\frac{n^\alpha \rho_n^n}{2c}(\hat\rho_n - \rho_n) \Rightarrow C}

with `C` standard Cauchy. Remark (i)/eq. 27 also gives the simpler
**fixed-root exact-explosive case** (White 1958, restated there). It needs
no \eqn{\alpha} or \eqn{c}, which makes it the more usable form for a
plug-in CI:

\deqn{\frac{\rho^n}{\rho^2-1}(\hat\rho_n - \rho) \Rightarrow C}

Plugging \eqn{\hat\rho} in for the unknown \eqn{\rho} in the normalization
(standard practice for a self-normalized pivot of this kind), a two-sided
CI is \eqn{\hat\rho \pm q_{\alpha/2} \cdot (\hat\rho^2-1)/\hat\rho^n},
where \eqn{q_{\alpha/2}} is a standard-Cauchy quantile (`qcauchy()`, or
equivalently `qt(., df = 1)`, already verified exactly above). This
formula is **simpler than, and different from, the review's garbled
restatement** quoted earlier (`sqrt(rho^(2n)-1)/rho^n`). Now that we have
read the primary source directly, treat that restatement as wrong and
superseded, and not as an equivalent alternative form.

**Update (2026-08-09, Bundle 1): implemented.** `root_ci(x, type = "cauchy")`
now ships this fixed-root form (eq. 27) as a second CI type alongside the
default `"normal"` one. The roxygen docs note that the Cauchy interval
assumes a *fixed* explosive root, while the default normal-t interval
(Guo/Sun/Wang) allows drift and dependence and stays the safer default.
Tests in `test-rootstamp.R` check that it brackets the point estimate and
that it matches the closed-form eq. 27 formula exactly. The true moderate-deviations form (eq. 26) needs an
estimate of \eqn{\alpha}, the localizing-rate exponent. That is a materially harder
follow-on than the fixed-root form, and we have not attempted it.

Also relevant, but not implemented: Phillips, Magdalinos & Giraitis (2010,
J. Econometrics 158(2), 274-279, "Smoothing local-to-moderate unit root
theory", open as Cowles DP 1659)
shows that the moderate-deviations theory above smooths continuously into
the local-to-unity case as \eqn{\alpha \to 0}. It is background theory and
not a standalone feature. It is the citation to reach for if `root_ci()`
ever needs to handle roots close to the local-to-unity boundary and not
only a clearly mildly explosive \eqn{\rho}.

### Update (2026-08-09, Bundle 1): `root_ci_datestamp()` instead of `summary()`

We first planned to wire `explosive_root()`/`root_ci()` directly into
`radf_obj`'s `summary()` method, and that turned out to be a worse fit than
we had scoped. The `summary()` S3 dispatch for `radf_obj`
(`summary_radf.mc_cv`/`.wb_cv`/`.sb_cv`) is built entirely around
`radf_cv` test-statistic critical values, while root CIs need a
`datestamp()` result. That is a structurally different input with no
natural slot in the dispatch chain. Restructuring shared `summary()`
machinery used by three other critical-value types to accommodate a
different kind of output was a bigger and riskier change than "wire it in"
suggested.

We shipped `root_ci_datestamp(object, ds, level = 0.95, type = "normal")`
instead. It is a standalone function that runs `explosive_root()`/`root_ci()`
on every episode in a `datestamp()` result and returns output in the same
per-series named-list shape that `datestamp()` uses. It is tested end to end
in `test-rootstamp.R` against a real `radf()` → `radf_mc_cv()` →
`datestamp()` pipeline, and checked against calling
`explosive_root()`/`root_ci()` directly on the same episode. Root inference
on very short episodes (duration 1-2) degrades as calling `explosive_root()`
directly on 2-3 points would, and we left that visible and did not filter
such episodes silently. The docs point to `datestamp()`'s existing
`min_duration` argument, and we added no second filtering knob.

Replication script:
[replication/dating-and-root-inference/rootstamp_validation.R](#script-rootstamp_validation).

---

## Confidence sets for bubble dates

**Status: evaluated and not implemented; re-triaged 2026-08-10. The
re-triage found real structure, but the work is still multi-step and we did
not complete it.** We read the full PDF (abstract, intro, model) and wrote
no R/C++ code.

### Source

Kurozumi, E. & Skrobotov, A. (2025). "Confidence Sets for the Emergence,
Collapse, and Recovery Dates of a Bubble." arXiv:2511.16172.

### What it is

A CI layered on top of already-estimated dates. It is conceptually the
dating analogue of [root inference](#root-inference)'s `root_ci()`, and it
is not a new detection or point-estimation method. It does **not** use the
obvious approach, the limiting distribution of the breakpoint estimator
itself, because the paper reports that its preliminary simulations found
that approach performs poorly for bubble dates. Instead it builds confidence
sets by *inverting* hypothesis tests for the break location: a
likelihood-ratio-type test (Eo & Morley 2015) and Elliott-Müller-type
(2007) tests, used individually and combined, with new limiting null and
alternative distributions derived for each and evaluated by Monte Carlo. The
three dates (emergence, collapse, recovery) are estimated separately and not
jointly.

### Cost/feasibility note for exuber

This is not a "wrap an existing point estimate in `± z·se`" CI, like
`root_ci()` for the explosive root (Bundle 1). The paper had to invent a
new route *because* the naive analogue does not work here. A faithful
implementation needs several new test statistics (LR-type and
Elliott-Müller-type), each with its own critical values, a rule for
combining them, and per-date confidence-set construction (inverting a test
statistic, not adding a margin to a point estimate), repeated for three
dates. Its own precondition, a WLS/volatility-corrected `dating_pdc()`
variant, has since shipped (see
[below](#wls-dating-under-time-varying-volatility)), and that prompted a
re-triage.

**Re-triaged (2026-08-10) by re-reading rendered pages 8, 10-11 and 46-47.**
The critical values are cheaper than the original "needs new simulation for
everything" framing suggested. The statistic construction, however, bears
out the original estimate of a multi-step cost on the scale of HLS/HLW, and
there is no hidden quick win.

- **Favourable**: not every statistic needs simulated critical values.
  `LR^e_{a,12}` and `EM^e_{a,12}` (their eq. 15-16, the "12"-direction
  tests) have an **exact closed-form chi-square critical value**,
  `cv^e_{LR12,0.05} = λ1·χ²_{1,0.05}` and
  `cv^e_{EM12,0.05} = sqrt(λ1·χ²_{1,0.05})`, with no simulation at all. The
  "21"-direction tests (`LR^e_{a,21}`, `EM^e_{a,21}`, `EM^e_{b,21}`, eq.
  17-19) have no closed form, but the paper *publishes* a response-surface
  regression for their critical values (the equation above their Table 1,
  `cv = a_{0,ℓ} + a_{-1,ℓ}/λ1* + a_{1,ℓ}·λ1* + a_{2,ℓ}·λ1*² + a_{3,ℓ}·λ1*³`,
  with coefficients transcribed from their Table 1). This is a
  MacKinnon-style formula, and not a table to interpolate or a simulation
  to run. It follows the same "published, so not new" pattern that made
  Kurozumi (2020)'s and HB's own boundaries cheap elsewhere in this
  project.
- **Still bigger**: the statistics themselves are not a thin reuse of what
  has shipped. The paper's recommended test (their "`LE^e` test", which
  combines `LR^e_{b,12}` with `EM^e_{a,21}` and was chosen over the naive
  `LR^e_{a,12}` because that one is "over-sized in finite samples") needs
  two things. (a) `LR^e_{b,12}` (their eq. 11) is a `min` over candidate
  break dates of
  `(y²_{T2} - ρ̂_a · Σ_{t=T1+1}^{T2} y²_{t-1}) / (T·φ̂_a^{2(T2-T1)}·σ̂²/2)`.
  Its numerator follows the familiar prefix-sum-window pattern (`Σy²_{t-1}`
  via cumulative sums, exactly as in `dating_hls()`'s `hls_prefix_sums()`).
  But we did not pin down the nuisance-parameter estimators `ρ̂_a`, `φ̂_a`,
  `σ̂²` or the construction rule for the admissible-break-date set
  `Λ^e_{12}`, because that needs more of Section 2's model setup than we
  read. (b) `EM^e_{a,21}` (their eq. 18) is an *integral* over a continuum
  of candidate break points of an `ADF(λ2*, λ1*)` functional (a ratio of
  Brownian-motion-type functionals, defined in eq. 18). That is new
  estimation machinery with no exuber analogue, and it is not a discrete
  min/max search like `LR^e_{b,12}`.

Net: the work is no longer accurately described as "new estimation *and*
new critical-value simulation". By transcribing the formulas we found
that the critical-value half is cheap (closed form or a published response
surface). The statistic construction still needs new machinery for at least
one of the two components of the paper's recommended test, at a scale
comparable to [HLS](#ssrbic-dating-vs-psy-recursive-dating)/HLW and not a
same-day addition. That is before repeating any of it for the collapse and
recovery dates. We did not complete it. A future pass can start from this
note: `LR^e_{b,12}`'s SSR-style numerator is prefix-sum-ready, so find
`ρ̂_a`/`φ̂_a`/`σ̂²`/`Λ^e_{12}` in Section 2 and then tackle `EM^e_{a,21}`'s
integral separately, without re-deriving what Table 1 means.

---

## Improved retrospective dating

**Status: done. The single-bubble omission fix shipped 2026-08-10, and
Section 3's multi-bubble dynamic programme shipped 2026-09-29 as
`dating_knp(breaks = )`.** We read the full PDF (abstract, intro, model,
Theorems 1-2) and re-verified eq. 1-8, the HLS-equivalence footnote, and
Section 3's DP algorithm description against rendered PDF pages 3-4.

### Source

Kejriwal, M., Nguyen, L. & Perron, P. (2025). "An Improved Procedure for
Retrospectively Dating the Emergence and Collapse of Bubbles." *JTSA*,
46(5), 867-883. `doi:10.1111/jtsa.12810`.

### What it is

A different bias fix for the *same* joint-SSR family as HLS, and not a
relative of PDC/KS. The model uses HLS's **fixed** autoregressive
coefficient framework (`rho` fixed, and not Phillips-Magdalinos's "mildly
explosive" `rho_T -> 1`) with an abrupt collapse and not a stationary
transition. Theorem 1 shows that the standard OLS joint-SSR estimator is
inconsistent. The origination-date estimate converges to the *collapse*
date, and the collapse-date estimate converges to a date *after* the true
collapse, offset by the trimming parameter. Both are biased late, for a
different reason than PSY's own threshold-crossing delay. The fix
(Theorem 2) is a "modified SSR" that **omits the single residual at the
implosion date** from the objective function, which restores consistency
for both dates. A footnote shows that this omission is numerically
*equivalent* to a specific one-time-dummy modification of HLS's Model 4
regression. The paper is therefore best read as a bias fix inside the HLS
estimating equation, and not as a third independent dating family. Section
3 also develops a Bai-Perron/Perron-Qu-style dynamic-programming algorithm,
so that the multi-bubble case avoids HLS's brute-force combinatorial grid
search. It exploits the unit-root restriction on the non-bubble regimes to
skip the iterative initial-values step that ordinary Bai-Perron/Perron-Qu
DP needs.

**A structural finding, from re-deriving the closed forms directly.** KNP's
single-bubble model (their eq. 1-3) is not a generic "intercept+AR(1)-dummy"
structure that needs new estimation code. It has exactly the shape of HLS's
Model 2: an unfitted unit root, then an intercept+slope-fitted explosive
regime, then an unfitted unit root resuming after an instantaneous
collapse. Their eq. 2's `delta_hat` formula is algebraically the standard
bivariate OLS slope (the `(x-xbar)*ybar` cross term vanishes by
construction). Regressing the *level* `y_t` on `y_{t-1}` with an intercept
gives the same residuals and SSR as regressing `Delta y_t` on `y_{t-1}` with
an intercept (a fixed reparameterization, `slope' = delta - 1`), and that
is the regression `hls_segment_ssr()` in `dating_hls.R` already computes for
HLS's Model 2 search. The whole fix reduces to
`SSR_om(T1,T2) = SSR(T1,T2) - (Delta y_{T2+1})^2`: one already-computed
squared term subtracted from an already-computed SSR, with no new
regression at all.

### Implementation

Shipped as `dating_knp(data, trim = 0.05, omit = TRUE)` in
`exuber/R/dating_knp.R`, tested in `exuber/tests/testthat/test-knp.R`. The
internal helper `knp_find_break()` reuses `hls_prefix_sums()` and
`hls_segment_ssr()` from `dating_hls.R` directly, since the structural
finding above means no new closed-form derivation is needed. It searches
`(tau1, tau2)` jointly to minimise the omission-corrected SSR. With
`omit = FALSE` it minimises the plain SSR instead, which is provably
inconsistent. We kept that option so the correction's effect can be
demonstrated and tested directly. Unlike
`hls_model23()`, KNP's candidate set imposes no directional sign constraint
on the fitted "peak".

**Validation, including a direct reproduction of the paper's own theorems
and not only a plausibility check.** The formula check is exact: both the
`omit = FALSE` and `omit = TRUE` searches match an exhaustive brute-force
nested-`lm()` search. The Monte Carlo uses KNP's own DGP (unit root,
no-intercept explosive AR(1), instantaneous collapse back near the
pre-bubble level, fresh unit root; 30 reps, `T1=50`, `T2=90`, `T=200`,
`delta=1.05`):

- **Theorem 1 (naive, `omit = FALSE`) reproduced directly.** The
  origination-date estimate's mean error relative to the *true collapse
  date* (`mean|tau1_hat - T2| = 1.0`) is far smaller than its error relative
  to its *own* true origination date (`mean|tau1_hat - T1| = 39.0`). The
  naive estimator's `tau1_hat` therefore tracks the wrong date and converges
  to `T2` and not `T1`, as Theorem 1 predicts.
- **Theorem 2 (omission-corrected, `omit = TRUE`) reproduced directly.** The
  same bias falls from 39.0 to 13.0 observations. That is a real correction
  but not a full elimination at this finite `T`. Theorem 2 is an asymptotic
  `→p` result, so residual finite-sample bias at `T=200` is expected and is
  not a defect. The collapse date and the explosive coefficient `delta_hat`
  are both close to their true values (`mean|tau2_hat - T2| = 1.0`;
  `delta_hat` mean `0.986` vs. true `1.05`).

Replication script:
[replication/dating-and-root-inference/radf_knp_validation.R](#script-radf_knp_validation).

### Implementation, multi-bubble dynamic programme (2026-09-29)

**We re-triaged this from "new algorithmic machinery".** The original
verdict overstated the effort. Section 3.2's algorithm is a textbook
Bai-Perron segment DP, and what makes it cheap is what KNP emphasise: the
unit-root regimes are *restricted* (`mu = 0`, `rho = 1`), so each segment's
restricted SSR is known in closed form and Perron-Qu's iteration over
initial values is unnecessary. Every segment cost is one of two quantities
that `dating_knp()` already computed through `hls_prefix_sums()`: `sum z^2`
over the segment (minus its first term when it follows a collapse and
`omit = TRUE`), or the intercept+slope OLS SSR
(`hls_segment_ssr(..., fit = TRUE)`). The objective is their eq. 11.
Regimes alternate between unit root and explosive, starting with a unit
root, and every unit-root regime after the first omits its first residual.
The DP is `O(m T^2)` with `O(1)` segment costs, and it returns the exact
global minimiser of the grid search.

Shipped as `dating_knp(data, trim, omit, breaks = 2L)` (`knp_dp()` in
`exuber/R/dating_knp.R`). `breaks` is the paper's `m`: two per bubble, or an
odd number to let the last bubble run to the sample end (its collapse is
`NA`). As in the paper, `m` is taken as given, because KNP leave its
selection open (their Section 4 conditions on the correct number).
`breaks = 2` keeps the original exhaustive single-bubble search. With more
breaks, `origination`, `collapse` and `delta` become one-row-per-bubble
matrices.

**Validated**:

- **Exact**: `knp_dp(y, 2)` returns the same break dates and SSR as the
  single-bubble exhaustive search (`|dSSR| = 0`). With 3 and 4 breaks it
  matches a brute-force search over every admissible partition, with
  segment SSRs from `lm()`, both with and without omission (`n = 28`,
  `|dSSR| = 2.7e-14`).
- **Monte Carlo** on a two-bubble version of KNP's DGP (`T = 200`,
  bubbles over 41-70 and 121-150, `delta = 1.05`, instantaneous collapse
  back near the pre-bubble level, 50 reps). Mean absolute date error per
  break (origination 1, collapse 1, origination 2, collapse 2):
  **11.9 / 5.9 / 10.3 / 3.8** observations with the omission correction
  vs. **33.4 / 18.5 / 29.6 / 13.2** without it. Theorem 1's
  inconsistency carries over to the multi-bubble case, and so does the
  fix.

Tests in `exuber/tests/testthat/test-knp.R`; replication script
[replication/dating-and-root-inference/radf_knp_validation.R](#script-radf_knp_validation)
sections 4-5. Ported to pyexuber (`_knp_dp()`), cross-checked against R
and the same brute force.

---

## WLS dating under time-varying volatility

**Status: done (2026-08-09).** Shipped as `dating_pdc(..., type = "wls")`.

### Source

Kurozumi, E. & Skrobotov, A. (2023). "Improving the accuracy of bubble date
estimators under time-varying volatility." arXiv:2306.02977.

### What it is

A direct two-step generalization of PDC/KS's sequential dating estimator,
by the same authors as the
[KS (2023) 4-regime extension](#3-pdc-2021-journal--ks-2023-journal-sequential-sample-splitting)
that `dating_pdc()` already implements. It builds on that estimator
explicitly ("we estimate these break dates as proposed by PDC and Kurozumi
and Skrobotov (2022) and collect the residuals..."). Step 1 is exactly
`pdc_find_break()`/`dating_pdc()` as already implemented: fit the
homoskedastic no-intercept AR(1) break model, and get consistent break-date
(fraction) estimates and their residuals. Step 2 estimates the time-varying
error variance `sigma_t^2` nonparametrically from those residuals, then
re-estimates each break date by minimizing a **weighted** SSR,
`sum(y_t - a*y_{t-1})^2 / sigma_t^2`, in place of the unweighted sum.
Algebraically this is `pdc_find_break()`'s cumulative-sum trick with every
`y_t`, `y_t*y_{t-1}`, `y_{t-1}^2` term divided by `sigma_t^2` before the
prefix sum. It stays closed-form and `O(T)` per breakpoint, given the
volatility weights.

### Implementation

Shipped as `dating_pdc(data, ..., type = c("ols", "wls"))` in
`exuber/R/dating_pdc.R`. `type = "ols"` is the original PDC/KS estimator,
unchanged. `type = "wls"`:

1. Runs the existing sequential `type = "ols"` fit to get step-1 break
   estimates.
2. `pdc_regime_resid(y, breaks)` (new) computes the fitted no-intercept
   AR(1) residual at every `(y_{t-1}, y_t)` pair, using one OLS `rho` per
   regime implied by the step-1 breaks. This is the paper's "collect the
   residuals of the fitted [regime] model."
3. Those residuals feed `nw_spot_vol()`, the Nadaraya-Watson kernel smoother
   with a leave-one-out cross-validated bandwidth that exuber already had.
   We **extracted** it from `kernel_spot_vol()` in `radf_sbz.R` and did not
   duplicate it, so SBZ (which smooths squared first differences) and this
   estimator (which smooths squared regime residuals) share one
   implementation. `kernel_spot_vol(y)` is now a one-line wrapper,
   `nw_spot_vol(diff(y))`, and its behaviour is unchanged: the full-suite
   regression run shows that `test-sbz.R` still passes after the
   extraction.
4. `pdc_find_break()` gained an optional `weights` argument (`NULL` gives
   the original unweighted behaviour, which a dedicated regression test
   verifies is identical). Every cumulative sum is multiplied by the weight
   vector before the prefix sum, so the search remains the same `O(T)`
   closed-form scan.
5. The full sequential search (collapse, then origination, then recovery)
   is re-run once with `weights = 1 / sigma_t^2`, using the correct
   contiguous slice of the full-sample `sigma_t^2` vector for each
   sub-sample regression.

No new critical-value simulation is needed. As with the OLS version, this
is point estimation and not a threshold-crossing test.

**Independent validation**: two Monte Carlo checks (40 seeds each, and not a
single cherry-picked run), following this project's usual bar of
distinguishing "the estimator works" from "it happened to pass once":

- **Homoskedastic DGP** (no volatility signal to exploit): OLS and WLS
  origination-date MAE are statistically indistinguishable (5.42 vs 5.53).
  WLS costs essentially nothing when there is nothing to gain, as expected
  of a nonparametrically weighted estimator when there is no
  heteroskedasticity to detect.
- **Heteroskedastic DGP with a volatility burst in the first 20% of the
  pre-bubble regime** (the scenario for which the paper's own Monte Carlo
  reports the largest gains): origination-date MAE drops from **13.05 (OLS)
  to 2.33 (WLS)**, an improvement of about 5.6 times. This confirms that
  the mechanism works as claimed. OLS's unweighted objective lets the noisy
  early segment dominate the origination split, and WLS downweights it
  through the estimated spot variance. Collapse-date accuracy is unaffected
  either way. Both are already near-exact, since the explosive-to-collapse
  transition dominates the SSR regardless of earlier noise, consistent with
  PDC's stochastic-order argument that the collapse is identified first.

A regression test in `test-pdc.R` checks a loose 2x version of this margin,
so a future change that erodes the benefit is caught without the test being
brittle to the exact numbers on a different RNG/BLAS.

Replication scripts:
[replication/dating-and-root-inference/radf_pdc_wls_heteroskedastic_mae.R](#script-radf_pdc_wls_heteroskedastic_mae),
[radf_pdc_wls_homoskedastic_mae.R](#script-radf_pdc_wls_homoskedastic_mae).

---

## Reverse-regression recovery dating

**Status: implemented (2026-08-10) and shipped with caveats. `f_r` (the
recovery date) validates well. `f_c` (the crisis-origination date) and the
overall false-detection rate under the null are noisier than we hoped, and
we have not fully resolved them.** Shipped as
`radf_recovery()`/`radf_recovery_cv()`. We read the full PDF (abstract
through Section 3.2, the model, both forward and reverse BSDF definitions,
both limit-theory theorems and the finite-sample Monte Carlo setup). A
second, more careful pass over Section 4.3's real-time monitoring extension
corrected a mis-transcription in our first reading, described below.

### Source

Phillips, P.C.B. & Shi, S. (2014). "Financial Bubble Implosion." Working
paper: Cowles Foundation DP 1967. Published as "Financial Bubble Implosion
and Reverse Regression," *Econometric Theory*.

### What it is

The mechanism is simple. Reverse the series (`X*_t := X_{T+1-t}`), run the
*same* BSDF/BSADF recursion PSY already uses on `X*`, and map the
crossing-time fractions back to the original time index (their eqs. 8-9):

```
f_hat_r = 1 - g_hat_e,  g_hat_e = inf{ g in [g0, 1]    : BSDF_g(g0) > scv }   (recovery date)
f_hat_c = 1 - g_hat_c,  g_hat_c = inf{ g in [g_hat_e,1]: BSDF_g(g0) < scv }   (crisis-origination date)
```

**A correction to our first transcription.** `f_hat_c` is *not* "a
further/later correction after recovery". Re-reading the extracted text
directly against eqs. 8-9 (and not just the surrounding prose) shows that
`g_hat_c` is searched only *after* `g_hat_e`, so `f_hat_c <= f_hat_r`
always. `f_hat_c` is the ORIGINAL series' crisis or collapse-onset date,
re-derived by reverse regression as an alternative to the collapse date
PSY's forward test already dates, and `f_hat_r` (recovery) comes
chronologically *after* it. The paper says so explicitly: "market recovery
(`f_hat_r`) following a crash begins when normal market behavior changes to
exuberance in the reverse series (`g_hat_e`)... market collapse in the
original series begins when exuberance in the reverse series shifts to
collapse at (`g_hat_c`)." Reversing a mildly explosive-then-mildly-integrated
collapse process (model (2)/(7) in the paper) turns the collapse regime into
an explosive regime in reverse time, and vice versa. Detecting the
crisis-origination and market-recovery dates thus becomes detecting
explosiveness in the reversed series, which is the right-tailed test PSY
already runs, applied to `rev(x)`.

**The catch, and why this is not free.** Theorem 1 derives the null limiting
distribution of the reverse statistic, `F_g(W, g_0)`, and it is **not** the
same distribution as the forward statistic's `F_f(W, f_0)`. It has an extra
term, because reversing a random walk makes the reversed "lagged" regressor
correlated with the reversed current error
(`E[X*_{T-j+2} * eps_{T-j+2}] != 0`, stated below Theorem 1). This is a
generic consequence of running a regression on time-reversed data and has
nothing to do with bubbles, so it changes the reverse test's critical values
*even under the null* and not just its finite-sample power. We also checked
this empirically. A paired Monte Carlo (same underlying draws, `n=100`,
`minw=20`, 5,000 reps) compared `radf_mc_cv()`'s standard forward critical
values with the same recursive computation run on the reversed path. The
values differ measurably (mean absolute difference about 0.04, max about
0.11 at the 95% level across positions), as the paper's theorem implies.

**A separate mechanism, Section 4.3's real-time monitoring extension**
(their eqs. 10-11), applies the same reverse-regression machinery
repeatedly on a *growing* sample from the collapse date `T_c` forward,
stopping at the first sample end `K` for which a correction is detected.
This produces the paper's "further correction in January 2004... full return
to normal market conditions in May 2004" figure in their dot-com empirical
application (Section 5). It is a *separate, sequential* application that
starts from the peak, and is not part of the eqs. 8-9 `f_c`/`f_r` pair. We
have **not implemented** it. It would be a natural, cheap follow-on, given
that `radf_recovery()` exists: an outer loop that expands the sample and
re-calls `radf_recovery()`-style logic at each step, structurally similar to
`monitor()`. It was out of scope for this pass.

**Numbers from the paper's own dot-com illustration** (Section 5, in prose,
NASDAQ price-dividend ratio): the eqs. 8-9 pair gives a crash from March to
November 2000 (so `f_c` = March 2000 and `f_r` = November 2000), and the
*separate* Section 4.3 monitoring extension then gives the further
correction in January-May 2004. We have not reproduced these, because we
have no access to the underlying NASDAQ price-dividend series.

### Implementation

Shipped as `radf_recovery()` (main function) and `radf_recovery_cv()`
(reversal-calibrated Monte Carlo critical values). The latter mirrors
`radf_mc_cv()`'s simulate-then-quantile construction, including its
`cummax(badf)`-as-`bsadf`-boundary shortcut, with one added `rev()` before
the recursive computation. We chose it over reusing forward critical values
as an approximation of unknown quality. The implementation is as cheap as we
hoped: `radf()`'s existing `bsadf` recursion runs on the reversed series, is
compared with the reversal-calibrated boundary, and the first
up-crossing/down-crossing pair is mapped back via `f = n + 1 - g`. It needs
no C++ and no new point statistic.

**Validation.** The structural invariant `f_c <= f_r` holds by construction
whenever both dates are identified and uncensored, and it held in all
replications we tried. `f_r`'s bias on synthetic collapse-then-recovery
data is small (a few observations), in the same direction and of roughly the
same size as the paper's own Table 5 finding of about 6 observations early.
`f_c`'s bias is materially larger (mean |bias| approaching the length of the
synthetic collapse window itself in some runs). The false-detection rate
under a pure random-walk null (`n=100`, `minw=20`, 95% level, one stable
reversal-calibrated cv reused across 200 fresh draws) is around 29%. That is
higher than comparable forward-test numbers elsewhere in this project, such
as `monitor()`'s cumulative false-alarm rate of about 10% over a 75-point
horizon. We found and fixed one artifact of the synthetic DGP during this
validation. An abrupt level jump at the expansion-to-collapse regime
boundary produced a spurious, narrow bsadf spike right at the junction and
not across the intended collapse regime. We replaced it with a smooth,
continuous mean-reverting transition. That improved `f_r`'s bias
substantially but did not resolve `f_c`'s. We have a plausible explanation
that is not alarming, but we have not tested it. Eq. 9's `inf` operator
defines the *first* down-crossing with no persistence requirement, so a
transient noise-driven dip below the reversal-calibrated boundary is enough
to trigger a premature `f_c`. If so, this is a property of the paper's own
literal crossing rule under finite-sample noise and not necessarily an
implementation bug, but we have not ruled out a subtler code issue.
Replication script:
[replication/dating-and-root-inference/radf_recovery_validation.R](#script-radf_recovery_validation).

We shipped it anyway, at the user's explicit call after being shown the
tradeoff. The mechanical and structural parts are sound and independently
verified (the invariant, the differing critical values, `f_r`'s accuracy),
and debugging the residual `f_c` and false-detection concern is better
scoped as its own follow-up than as an open-ended extension of this pass.
The roxygen docs of `radf_recovery()` carry the same caveat inline.
