---
title: "dating-and-root-inference"
blurb: "Origination, collapse and recovery dates, plus confidence intervals on the explosive root itself."
order: 2
---
﻿# Dating and root inference

Two related post-detection problems, both operating on an episode
`radf()`/`datestamp()` has already flagged as explosive: **dating** (when
exactly did it start/end/recover?) and **root inference** (how explosive —
what's `rho`, with a confidence interval, and how fast is it doubling?).
Grouped together because they're the same step in exuber's workflow
(`radf()` → `datestamp()` → sub-sample → refine), not because the methods
share statistical machinery. Status legend as in
[volatility-robustness.md](/replication/volatility-robustness).

| Method | Paper | Status |
|---|---|---|
| [SSR/BIC dating — PDC/KS route](#ssrbic-dating-vs-psy-recursive-dating) | Pang/Du/Chong (2021), Kurozumi/Skrobotov (2023) | **done** |
| [SSR/BIC dating — HLS/HLW route](#ssrbic-dating-vs-psy-recursive-dating) | Harvey/Leybourne/Sollis (2017), Harvey/Leybourne/Whitehouse (2020) | **done** (both HLS single-bubble and HLW multi-bubble routes) |
| [Root inference (Cauchy CI + normal-t CI)](#root-inference) | Phillips & Magdalinos (2007), Guo, Sun & Wang (2019) | **done** |
| [Confidence sets for bubble dates](#confidence-sets-for-bubble-dates) | Kurozumi & Skrobotov (2025) | re-triaged 2026-08-10 — critical values cheap (closed-form/published response surface), statistic construction still multi-step, not implemented |
| [Improved retrospective dating](#improved-retrospective-dating) | Kejriwal, Nguyen & Perron (2025) | **single-bubble omission fix done**; multi-bubble Bai-Perron/PQ algorithm not implemented |
| [WLS dating under time-varying volatility](#wls-dating-under-time-varying-volatility) | Kurozumi & Skrobotov (2023) | **done** |
| [Reverse-regression recovery dating](#reverse-regression-recovery-dating) | Phillips & Shi (2014/2019) | **done, shipped with caveats** |

All papers: [references.md](/replication/references#dating-and-root-inference).

---

## SSR/BIC dating vs. PSY recursive dating

**Status: PDC/KS route done (2026-08-09); HLS (2017)'s single-bubble
SSR+BIC route done (2026-08-10); HLW (2020)'s multi-bubble two-step
wrapper around it also done (2026-08-10).** The PDC/KS sequential
sample-splitting estimator is implemented as `dating_pdc()` — see
[Implementation (PDC/KS route)](#implementation-pdcks-route) below. HLS's
four-model grid-search BIC approach is implemented as `dating_hls()` — see
[Implementation (HLS route)](#implementation-hls-route) below, including
a real sign-constraint bug found and fixed shortly after it first
shipped. HLW's two-step extension is implemented as `dating_hlw()` — see
[Implementation (HLW route)](#implementation-hlw-route) below.

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

All four were read as full PDFs, not secondary summaries. HLS and HLW were,
for a long time, the two hardest papers in this whole project to obtain —
both are formally paywalled at Elsevier/ScienceDirect, and their Green-OA
mirror at the University of Nottingham's institutional repository
(`nottingham-repository.worktribe.com`) is behind a Cloudflare bot-check
that blocked every automated route tried, including a real Playwright
browser session over an institutional network (see
[references.md](/replication/references#the-two-that-automation-couldnt-get-hls-and-hlw)
for the full history — a version of this file once thought they'd been
recovered via a direct Wayback Machine snapshot at
`.../preview/829483/bubble_dates.pdf`/`.../preview/4418036/bubble_dates.pdf`,
and the same filename convention is in fact how the working copies now in
the paper library finally arrived — supplied directly rather than by any tool
in this project). PDC's working paper is genuinely open (MPRA 92074); KS
is genuinely open (arXiv:2110.04500).

All formulas quoted below (HLS Models 1-4 and BICj, HLW's two-step
Models/BICj, PDC's Step 1/Step 2 estimators) were re-verified by rendering
the actual PDF pages to PNG (`pymupdf`, `matrix=fitz.Matrix(2.5,2.5)`) and
reading the typeset math directly — `pdftotext -layout` scrambles multi-line
stacked formulas and subscripts badly in all four PDFs (e.g. HLS's DGP
definition and SSR-order-of-magnitude table on p.4-5 render as
out-of-order fragments in plain text; the Model/BIC formulas on p.7-8,
checked by image, are not scrambled and match the OCR text closely once
laid out — but this was confirmed by rendering, not assumed). HLW's Table 1
(BIC model-selection frequencies) was also re-verified by image, since the
raw `pdftotext -layout` row/column output has a distinct offset artefact
(labels `A`, `B`, `C`... appear to be shifted relative to their `j=1`/`j=2`
data rows in the plain-text dump) that a first read could easily
misinterpret; the rendered image confirms the actual row/column pairing is
different from — and correct where the OCR text is ambiguous about — the
naive text-order reading.

### What it is

All three lines of work replace PSY's **threshold-crossing** dating rule
(find where the recursive BSADF/ADF t-statistic sequence first
crosses/re-crosses a critical value) with a **model-based SSR-minimisation +
BIC** dating rule. They differ in how they handle *multiple* bubbles in one
series.

#### 1. HLS (2017) — single-bubble SSR+BIC dating

Defines four candidate regime-structure models for a series `y_t`, each a
piecewise OLS regression of `Δy_t` on regime-indicator dummies and
dummy-interacted `y_{t-1}` (own notation, confirmed by image, p.7):

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
DGP/Model (i.e. consistent, not just consistent break-fraction, for the
exact date under a fixed-magnitude bubble).

Model selection across the 4 candidates uses a BIC with a specific penalty —
the number of fitted coefficients *plus* the number of estimated break
dates, times `ln(T)` (confirmed by image, p.8):

```
BIC1 = T·ln{T⁻¹SSR1(τ̂1,1)}        + (2+1)ln(T)
BIC2 = T·ln{T⁻¹SSR2(τ̂1,τ̂2)}       + (2+2)ln(T)
BIC3 = T·ln{T⁻¹SSR3(τ̂1,τ̂2,1)}     + (4+2)ln(T)
BIC4 = T·ln{T⁻¹SSR4(τ̂1,τ̂2,τ̂3)}    + (4+3)ln(T)
jopt = argmin_j BICj
```

Practical implementation (§5) imposes minimum regime durations
(`τ1 ≥ s`, `τ2−τ1 ≥ s`, `τ3−τ2 ≥ s/2`, with `s = 0.1` in the T=200
simulations, `s = 0.05` in the T=389 empirical application) but is otherwise
a **brute-force grid search** over 1, 2, or 3 breakpoints depending on the
model — no dynamic-programming (Bai-Perron-style) speedup is described in
the paper.

#### 2. HLW (2020) — two-step extension to multiple bubbles

Directly addresses PSY's known "late" end-date bias. Rather than building
an exponentially-growing model set for N bubbles (HLS's 4-model set
generalises combinatorially badly), HLW propose a **two-step** procedure
(confirmed by image, p.9, formulas match `pdftotext` text closely):

- **Step 1**: run PSY's existing GSADF/BSADF detection+dating exactly as-is
  to get preliminary start/end fraction estimates `τ̂P SY_j1, τ̂P SY_j2` for
  each of the `N̂` detected bubbles, and use these to carve the sample into
  `N̂` disjoint sub-sample "date windows" `[sj, ej]` (midpoint-split between
  consecutive PSY-detected regimes, with a rule to nudge the split so a
  window always starts inside a fitted post-explosive/unit-root regime, to
  avoid starting a window mid-bubble).
- **Step 2**: apply HLS's SSR+BIC Model 1-4 procedure *independently within
  each date window* (restricting to Models 2 & 4 for all but the final
  window, since a window boundary by construction is a unit-root point, not
  a sample end).

This is explicitly presented as bolting the HLS refinement onto PSY's own
existing output, not a replacement detection method — "we propose a dating
methodology based on minimum sum of squared residual estimators and BIC
model selection, but using prior information gleaned from the PSY dating
procedure as a means of reducing the dimensionality" (HLW, §1).

#### 3. PDC (2021, journal) / KS (2023, journal) — sequential sample-splitting

A structurally different (and computationally much cheaper) approach to the
*same* SSR-minimisation idea for a single bubble episode (3- or 4-regime
model), extended by HLW's own suggestion to the multiple-bubble case by
running it inside each PSY date window instead of HLS's joint fit.

PDC's model is unit-root → explosive → stationary-collapse (3 regimes, 2
breakpoints `τ1_0 < τ2_0`). Instead of HLS's *joint* minimisation of a
2-breakpoint SSR surface, PDC shows (via a stochastic-order argument, their
"Example 3"/Lemmas A.2-A.4) that under the bubble DGP the **collapse date is
always identified first** — the SSR drop at the collapse breakpoint
dominates (higher stochastic order) the SSR drop at the origination
breakpoint — so the two breaks can be estimated **sequentially** rather than
jointly:

```
Step 1: τ̂2 = argmin_{τ∈(0,1)} RSS2,T(τ),  RSS2,T(τ) = Σ_{t≤⌊τT⌋}(yt − β̂x(τ)yt-1)² + Σ_{t>⌊τT⌋}(yt − β̂3(τ)yt-1)²
        where β̂x(τ) = Σ_{t≤⌊τT⌋} yt·yt-1 / Σ_{t≤⌊τT⌋} yt-1²   (no-intercept AR(1) slope, full-sample split at τ)
Step 2: on the left subsample [1, τ̂2T] only, repeat the same one-break RSS minimisation to get τ̂1.
```

(confirmed by image, PDC p.9). This is a **no-intercept, single AR(1)
coefficient per regime** model (unlike HLS's intercept+AR(1) dummies) —
each `β̂(τ)` is a closed-form ratio of two prefix sums (`Σy_t y_{t-1}`,
`Σy_{t-1}²}`), so the whole `RSS(τ)` curve over all candidate `τ` is
computable in `O(T)` via cumulative sums, with **no joint grid search and no
model-selection BIC step at all** — PDC's algorithm always assumes the
3-regime structure holds and estimates its two breaks one at a time.

KS (2023) extends PDC's 3-regime model to 4 regimes (adding a final
unit-root "recovery" regime after the stationary collapse), reusing PDC's
sequential logic for the extra breakpoint, and *explicitly contrasts its own
cost against HLS's*: "we perform the three SSR minimization with one break
each with O(T) computations, while Harvey et al. (2017) requires minimizing
the three break model over all possible combinations of these breaks" (KS,
§1, p.3) — i.e. KS's own reading of HLS's Model 4 grid search is that it is
an un-sped-up multi-dimensional combinatorial search, corroborating the
"no dynamic programming" reading of HLS §5 above. KS also explicitly notes
their method is meant to slot into HLW's per-window second step: "Harvey et
al. (2020) proposed to initially identify the bubble regimes based on [the]
PSY approach ... one can use our approach in the second step" (KS, §1).

### Exact numbers reproduced

**HLS Table 1** (Nasdaq composite real price index, PWY's own 1973:2-2005:6
series), image-verified against `pdftotext` (both matched exactly — no
scrambling found here unlike HLW's Table 1):

| Sample | PSY test | PSY start | PSY end | BICopt model | BICopt start | BICopt end |
|---|---|---|---|---|---|---|
| 1973:2-2005:6 (full) | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2-2000:9 (pseudo-real-time) | 3.07*** | 1998:11 | 2000:9 | 1 | 2000:1 | 2000:9 |
| 1973:2-2000:10 | 3.07*** | 1998:11 | 2000:10 | 1 | 2000:1 | 2000:10 |
| 1973:2-2000:11 | 3.07*** | 1998:11 | 2000:11 | 1 | 1999:12 | 2000:11 |
| 1973:2-2000:12 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2-2001:1 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |

This is the single concrete quantitative comparison HLS publish between the
two dating rules: both agree on the 1998:11 start, but PSY's end date
(2000:12) is 3 months later than BICopt's (2000:9) — a direct illustration
of PSY's original dating strategy being late specifically on the *end*
date, in this worked example (the start dates agree here).

**What could not be reduced to an exact number**: HLS's own Monte Carlo
dating-accuracy comparison (§6, Figures 1-3) — the actual head-to-head
"BICopt vs PSY, % correct within k observations" evidence — is reported
**only as plots** (frequency-vs-bubble-magnitude curves), not as a table.
The text states the qualitative conclusion ("BICopt out-performs PSY... in
finite samples, particularly with respect to the bubble's end date") but no
specific percentage/RMSE figure is stated in prose either. Same situation in
HLW §4 (Figures 1-6 are histograms of PSY vs. BIC start/end date estimates;
prose describes patterns — "PSY estimates ... typically fall somewhat later
than the true date", "BIC estimated end dates equal the true end date in
almost every replication" — with no numeric table). This mirrors the
[SBZ finding](/replication/volatility-robustness#why-this-cant-be-a-bit-exact-numeric-cross-check-test):
any Monte Carlo comparison reported as a figure rather than a table cannot
be cross-checked bit-for-bit; it is reported here as a qualitative, not
quantitative, claim.

**HLW Table 1** (BIC model-selection frequencies across 6 DGPs A-F,
image-verified — the `pdftotext -layout` dump of this table has row labels
(`A`,`B`,`C`,`D`,`E`,`F`) that read as if offset from their `j=1/j=2/j=3`
data rows; the rendered image resolves this and the values below are read
from the image, not the raw text dump):

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

BIC picks the true model most often in every row; the weakest cases are
DGP A/j=2 and DGP C/j=2 (correct-model rate ≈ 44-45%, both cases where a
post-collapse reversion to unit root close to the end of the fitted window
is easily missed) — matching the paper's own caveat text exactly. This is a
genuine dating-relevant number (frequency of correct model identification),
but note it is **not** the same thing as dating accuracy — the paper is
explicit that incorrect model selection (e.g. Model 3 chosen over Model 4)
still typically gives accurate `τ̂1`/`τ̂2` estimates.

**KS empirical application** (§6, prose, not tabulated as such but
individually stated dates): NASDAQ Composite, Jan 1985-Aug 2013 monthly —
BIC4 (3387.727) beats BIC3 (3404.268) and BIC2 (3409.296), so the 4-regime
model is selected; collapse dated Feb 2000, origination (from the
left-subsample re-split) Aug 1998, recovery (right-subsample) Sep 2001. US
house price index (FHFA, real), Jan 1991-Dec 2012 — BIC4 (-2918.223) again
selected; collapse Nov 2006, origination Sep 1997, recovery May 2011.

**KS Monte Carlo** (§5): reported as histograms (Figures 1-4), not tables;
the only numbers stated directly in the prose are hedged/approximate ("the
frequency of selecting the true break date is not very high (it is around
30% for T=400 and 65% for T=800)... increases to almost 100%..."; "close to
100%"; "approximately 75% and 100% for T=400 and 800, respectively") — these
are quoted verbatim above as the paper's own words, not independently
computed, and should be read as order-of-magnitude, not exact, figures.

### Cost/feasibility note for exuber

exuber's current `datestamp()` (`exuber/R/radf-methods.R`, `datestamp.radf_obj`,
plus helpers `stamp()`/`stamp_to_index()`/`add_peak()` in the same file) is
**pure post-processing on top of statistics `radf()` already computed**: it
takes the already-computed `bsadf`/`badf` recursive-statistic sequences
(from `augment_join(object, cv)`), compares them against the simulated/wild
-bootstrap critical values at a chosen significance level
(`tstat > crit`), and finds contiguous runs of exceedance (`stamp()`,
literally `which(diff(x) != 1)` on the indices where the inequality holds).
No new regression is fit anywhere in this path — it is threshold-crossing
on numbers that already exist from `radf()`.

**This is why a `datestamp(method = "bic")` "contained, high-value
addition" framing undersells the actual work for the HLS/HLW route, but is
roughly right for the PDC/KS route** — these are two different amounts of
work:

1. **HLS-style BIC dating is not a post-processing step — DONE
   (2026-08-10) anyway.** It required, per bubble episode: (a) fitting up
   to 4 new OLS regression *specifications* (regime-dummy-interacted
   AR(1) models with intercept), none of which match exuber's existing
   `radf()` machinery; (b) a grid search over 1-3 breakpoints jointly per
   model; (c) a BIC comparison across the 4 fitted models. The `O(T^m)`
   scaling worry turned out to be avoidable without a Bai-Perron-style
   dynamic-programming rewrite: because the four models' regime dummies
   never overlap, any candidate partition's total SSR decomposes exactly
   into independent per-segment closed-form OLS fits, each a `O(1)`
   lookup into precomputed cumulative sums — the same style of trick
   `dating_pdc()`'s own (differently-specified) breakpoint search already
   uses. This keeps the *grid search itself* fast (`~2` seconds at
   `T=400` for the full 4-model search) even though it is still, in
   absolute terms, an `O(T)`/`O(T^2)`/`O(T^3)` search across the four
   models respectively — see [Implementation (HLS route)](#implementation-hls-route)
   below. HLW's step 1 (splitting into date windows) can legitimately
   reuse exuber's *existing* `radf()`/`datestamp()` output almost as-is
   (it's exactly PSY's own detected start/end dates, already what
   `datestamp()` returns) — step 2 (applying `dating_hls()` inside each
   window) is now a thin wrapper around already-shipped code rather than
   the "100% new estimation code" this note originally worried about, but
   the wrapper itself — window construction, the Models-2-and-4-only
   restriction for non-final windows — is still not built; see
   [Implementation (HLS route)](#implementation-hls-route) for why it was
   left as a documented follow-on.

2. **PDC/KS-style sequential sample-splitting is a meaningfully smaller,
   "contained" addition.** Each breakpoint estimate is a single `O(T)` scan
   of a *closed-form* ratio-of-cumulative-sums statistic (no intercept, no
   joint multi-dimensional grid, no BIC model-selection loop across 4
   specifications) — same complexity class as
   [STADF's `gls_dfstat_grid()`](/replication/volatility-robustness#formulas-implemented-verified-against-rendered-pdf-pages-not-ocr-text),
   i.e. plausibly implementable in R via `cumsum()`/`outer()` without new
   C++. The tradeoff is a narrower model: PDC/KS assume the regime count (3
   or 4) is known/fixed rather than BIC-selected among HLS's 4 alternatives,
   so there is no automatic distinction between "bubble that collapses" vs.
   "bubble that doesn't" vs. "bubble ongoing at sample end" — the practical
   value HLS/HLW's model-selection step actually buys. **This is the
   concrete recommended follow-up if this section gets picked up next.**

3. **Either way, this is a new statistical layer, not a flag on the
   existing rule.** Both routes need: new regime-dummy or no-intercept
   AR(1) fitting code (not reachable by editing `stamp()`/`datestamp()`),
   an explicit minimum-regime-duration trimming parameter (`s` in HLS/HLW,
   the `ψ`-style hyperparameters in PDC/KS) analogous to exuber's existing
   `psy_ds()`/`minw` but for a different purpose, and — for the HLS/HLW
   route specifically — a BIC computation and 1-3 dimensional grid search
   with no existing exuber analogue to build on. A plausible
   `datestamp(..., method = "bic")` is a genuine, self-contained R-only
   feature (PDC/KS variant) or a two-tier R feature reusing existing
   `datestamp()` output for step 1 (HLW variant) — but in neither case is
   it a small diff on top of `stamp()`/`add_peak()`; it is new estimation
   code of a size comparable to (HLS/HLW) or somewhat smaller than (PDC/KS)
   the [SBZ effort](/replication/volatility-robustness#implementation), not
   comparable to STADF's near-zero-cost formula reuse.

### Implementation (HLS route)

Shipped as `dating_hls(data, trim = 0.05)` in `exuber/R/dating_hls.R`, tested
in `exuber/tests/testthat/test-hls.R`. Fits all four of HLS's regime-dummy
models, each by exact SSR minimisation over its candidate breakpoint(s)
(subject to the paper's own minimum-regime-duration trim and
upward-bubble sign constraint, `y_{tau2} > y_{tau1}`), and selects among
them by the paper's own BIC formula (`n*log(SSR/n) + df*log(n)`, `df` in
`{3,4,6,7}` for Models 1-4 respectively). Returns the selected model
number and its origination/collapse/recovery dates (`NA` for whichever
of these the selected model doesn't have), plus every candidate model's
BIC value for inspecting how close the selection was.

**The key implementation move, not just an optimisation**: HLS's own
paper describes a brute-force grid search with no speedup, which the
cost/feasibility note above worried would scale as `O(T^m)` for an
`m`-breakpoint model with a large constant. Because the four models'
regime dummies never overlap, any candidate partition's SSR is exactly
the *sum* of independent per-segment OLS fits — a segment with no active
dummy has zero fitted parameters (`SSR = sum(Delta y_t^2)`); a segment
with an active dummy is a plain intercept+slope fit. Both reduce to a
closed-form ratio of cumulative sums (`Sx, Sxx, Sz, Szz, Sxz`), so every
candidate breakpoint (or breakpoint pair/triple) evaluates in `O(1)`
given precomputed prefix sums — no repeated `lm()` calls, the same style
of trick `dating_pdc()`'s own (differently-specified, no-intercept)
breakpoint search already uses. This keeps the search itself fast:
`dating_hls()` returns in well under 2 seconds at `T=400`, HLW's own
sample size, on plain R with no C++.

**Validation**: formula-exact — `hls_segment_ssr()` matches a
brute-force `lm()` fit's SSR exactly (tolerance `1e-8`) for arbitrary
segments, and, more importantly, the *full joint 3-breakpoint grid
search* for Model 4 matches an exhaustive brute-force nested-`lm()`
search bit-for-bit on a small synthetic series (not just the per-segment
formula in isolation). Monte Carlo, on synthetic regime DGPs built with
a large positive base level (`100 + cumsum(rnorm(...))`) so the
explosive signal is a genuinely large absolute-scale departure from
noise:

- **A true Model 4 DGP** (unit-root, bubble, mean-reverting collapse,
  unit-root recovery): BIC selects Model 3 or 4 (i.e. correctly
  identifies that a distinct collapse regime exists) in 100% of 30
  replications, split 80%/20% between the two — consistent with HLW's
  own reported finding that distinguishing a genuine final recovery
  regime from continued collapse is the hardest case in their own Table
  1 (correct-model rate as low as ~44-45% in their weakest DGP). The
  origination date's mean absolute bias is ~0 observations; the collapse
  date's mean absolute bias (among Model 3/4 selections) is ~1
  observation.
- **A true Model 2 DGP** (bubble fully reverting to a unit root, no
  distinct collapse regime): BIC selects Model 2 in 100% of 30
  replications; origination bias exactly 0 in every replication.
- **A true Model 1 DGP** (bubble ongoing at the sample end): BIC selects
  Model 1 in 100% of 30 replications; origination bias exactly 0 in
  every replication.
- **Pure `H0` (no bubble at all)**: BIC never selects the most complex
  model (0% Model 4 across 30 replications), splitting instead between
  the three simpler ones (30% Model 1, 57% Model 2, 13% Model 3) — no
  runaway complex-model selection under pure noise, though HLS's method
  has no "no bubble" option built in (it is designed to run downstream
  of a PSY-detected episode), so this is a curiosity check, not a core
  requirement.

**A DGP construction pitfall found and fixed during this validation, not
a code bug**: an early version of the synthetic bubble DGP built the
explosive regime as `unit1[n1] * c^i` where `unit1[n1]` is itself a
small, mean-zero random-walk value (potentially even negative) — this
gave a weak, sometimes-backwards "explosive" signal relative to
cumulative noise, and BIC correctly, rationally preferred the simplest
model (Model 1) under that weak signal ~80% of the time regardless of
the true DGP. Diagnosing this (by comparing the selected models' actual
SSR/BIC values, not just the pass/fail selection frequency) traced it to
DGP signal strength, not the search or BIC code — confirmed by rebuilding
the DGP on a large positive base level (matching the successful
construction already used for `quantile_test()`'s Monte Carlo checks
earlier this session).

**A real sign-constraint gap found and fixed shortly after this item was
first shipped, while reading HLW (2020) for the multi-bubble extension**:
HLW's own paper restates HLS's Models 1/3/4 with sign constraints this
file's first pass had abbreviated — Model 1 requires `y_T > y_{tau1}`
(the series must end above where the bubble started), and Models 3/4
require the fitted peak `y_{tau2}` to exceed not just the bubble's
starting level but also wherever the series ends up after the fitted
collapse regime (`y_T` for Model 3, `y_{tau3}` for Model 4) — otherwise
`y_{tau2}` isn't a genuine peak. The initial `dating_hls()` only checked
`y_{tau2} > y_{tau1}`, missing the second half of each constraint. Fixed
in all three affected search functions, with brute-force validation
updated to match; re-running the Monte Carlo checks after the fix showed
a further accuracy improvement, not a regression — Model 1 and Model 2
selection went from 97% to 100% correct with bias dropping to exactly 0
in every replication (numbers above already reflect the fix). Replication
script:
[replication/dating-and-root-inference/radf_hls_validation.R](#script-radf_hls_validation).

### Implementation (HLW route)

Shipped as `dating_hlw(data, cv = NULL, minw = NULL, trim = 0.1, min_duration = NULL, nboot = 199L, seed = NULL)`
in `exuber/R/dating_hlw.R`, tested in `exuber/tests/testthat/test-hlw.R`. A
two-step wrapper around `dating_hls()`, exactly as HLW's paper describes:

1. **Step 1**: run PSY's existing detection and dating
   (`radf()`/`datestamp()`, with `min_duration` defaulting to
   `psy_ds(n)` — HLW's own `ln(T)` minimum-episode-length rule) to get a
   preliminary start/end position for each detected explosive episode.
2. **Step 2**: carve the sample into disjoint date windows using HLW's
   own formula — the end of window `j` is the midpoint between that
   episode's PSY-detected end and the next episode's PSY-detected start
   (`e_j = tau2_psy[j] + floor((tau1_psy[j+1] - tau2_psy[j])/2)`, the
   last window running to the sample end) — then fit each window with
   `hls_fit_series()` (a new shared internal helper factored out of
   `dating_hls()` for exactly this reuse), restricted to Models 2 and 4
   for every window but the last, since a window boundary is by
   construction a unit-root/collapse point, not a genuine sample end.
   After fitting window `j`, HLW's sequential-adjustment rule sets the
   *next* window's start to the first observation of the just-fitted
   post-explosive regime (`s_{j+1} = s_j + tau2_local` for a Model 2 fit,
   `s_j + tau3_local` for Model 4) — this is what prevents a later window
   from starting mid-bubble, a failure mode HLW's own paper discusses
   explicitly.

**Two real bugs found and fixed while building this, not design
choices**:
1. `datestamp()` itself raises a hard error ("Cannot reject H0 at the 5%
   significance level"), not a warning, when no series has any detected
   episode — an initial version of `dating_hlw()` only caught `warning`
   conditions from this call, so a genuine no-bubble series crashed
   instead of returning the intended empty result. Fixed by also
   catching `error`.
2. While reading this paper for the window-construction rule, its own
   restatement of HLS (2017)'s Models 1/3/4 sign constraints turned out
   to be more complete than what the just-shipped `dating_hls()` checked —
   see the sign-constraint fix noted above. Fixed there, not here, since
   it affects `dating_hls()` directly and `dating_hlw()` inherits the fix by
   reuse.

**Validation** (structural, Monte Carlo, and a strong equivalence
check): `hlw_local_to_global()`'s index arithmetic checked directly.
Monte Carlo on a synthetic two-bubble DGP (20 reps): PSY step-1
detection found exactly 2 windows in 13/20 reps (the rest fragmented
into 3 or more spurious sub-windows — HLW's own paper explicitly
discusses this exact failure mode and proposes a run-joining heuristic
for it, not implemented here, see "Not implemented" below); among the
13 clean reps, origination and collapse date bias for *both* bubbles
was exactly 0 in every replication. Windows were always correctly
ordered/non-overlapping across all 20 reps. Under a pure `H0` null (no
bubble at all), `dating_hlw()` never errors and returned 0 detected
windows in all 20 reps tried (no false positives). On a single clean
bubble episode (15 reps), the wrapper's *final* window matched
standalone `dating_hls()` run on the whole series *exactly* (model,
origination, and collapse all identical) in 100% of the reps where a
final window existed — a strong structural correctness check, since
HLW's own paper states the two-step procedure should reduce to plain
HLS when there is only one episode. Replication script:
[replication/dating-and-root-inference/radf_hlw_validation.R](#script-radf_hlw_validation).

**Not implemented**: HLW's own run-joining heuristic for step-1
fragmentation ("if up to 3 non-rejections are surrounded on either side
by an explosive regime of length `ln(T)`, treat them as a single
episode") — `dating_hlw()` uses `datestamp()`'s regimes as detected,
un-joined, so a single true bubble can occasionally surface as multiple
windows in `dating_hlw()`'s output when PSY's own step-1 detection
fragments it (quantified above: 13/20 clean vs. 7/20 fragmented on the
two-bubble DGP, 12/15 vs. 3/15 on the single-bubble DGP). This is a
property of PSY's own step-1 detection noise, not of the window-
construction or per-window fitting logic (both validated exactly
correct above) — scoped out as its own small, well-defined follow-on
rather than folded into this pass.

### Implementation (PDC/KS route)

Shipped as `dating_pdc(data, regimes = 3L, trim = 0.05)` in `exuber/R/dating_pdc.R`,
tested in `exuber/tests/testthat/test-pdc.R`. Internal helper
`pdc_find_break(y, trim)` implements the single-breakpoint no-intercept
AR(1) RSS minimiser exactly as specified above (§"3. PDC/KS") — closed-form
`β̂(τ)` via prefix sums of `Σy_t y_{t-1}` and `Σy_{t-1}²`, so the whole
`RSS(τ)` curve is `O(T)`. `dating_pdc()` calls it sequentially per PDC's
Step 1/Step 2 (collapse first, then origination on the left subsample), and
once more on the right subsample for KS's 4-regime recovery date when
`regimes = 4`.

**Validation methodology** (independent of whoever wrote the original
code — this file's function was found already present, uncommitted, with no
test coverage, and was validated from scratch before being trusted):

1. **Formula-exact check**: `pdc_find_break()`'s closed-form break index and
   RSS match a brute-force scan that fits `lm(y[t] ~ y[t-1] - 1)` separately
   on every candidate left/right split, bit-for-bit (`tolerance = 1e-8`).
   Confirms the cumulative-sum algebra is not just fast but correct.
2. **Consistency in the low-noise/long-series/strong-effect limit**: on a
   synthetic 3-regime series (unit-root → explosive → collapse) and a
   4-regime series (...→ recovery), `dating_pdc()` recovers the true break
   dates to within 1-2 observations. This is the same technique used
   elsewhere in this project (e.g. SBZ, common-bubble) to distinguish "the
   estimator is asymptotically correct" from "the estimator happens to pass
   one finite-sample tolerance check" — a check that only holds at moderate
   sample sizes can't tell those apart, but convergence as noise shrinks
   and sample size grows can.
3. **A DGP-construction pitfall worth recording**: the first attempt at the
   4-regime consistency check used a *deterministic* exponential decay
   (`target + (peak - target) * exp(-k*t)`) for the collapse regime. This
   numerically saturates well before the regime ends, and its noisy "flat"
   tail near `target` is then statistically indistinguishable from the
   following recovery regime's random walk (both look like "small
   increments near a level" locally) — `pdc_find_break()` correctly finds a
   break, just the wrong one, *inside* the collapse regime rather than at
   its true boundary with the recovery regime. This is a test-construction
   bug, not an estimator bug: `pdc_find_break()` fits a genuine no-intercept
   AR(1), and a deterministic trend plus noise simply isn't one. Fixed by
   generating the collapse regime as an actual stationary AR(1) recursion
   (`y_t = ρ·y_{t-1} + ε_t`, `ρ = 0.5`) instead — homogeneous dynamics
   throughout the regime, cleanly distinguishable from the `ρ = 1` regimes
   on either side.
4. **Honest finite-sample characterization, not a rosy one**: a separate
   moderate-`T` test (`T ≈ 130` per series, 20 seeds) only asserts the
   estimates are finite and vary with the data — it does **not** assert
   tight recovery. This matches KS's own Monte Carlo (§5, quoted above):
   they report only ~30% exact-date recovery at `T = 400`, rising to ~65%
   at `T = 800`. A synthetic test that demanded tight accuracy at small `T`
   would either be lucky/cherry-picked or silently contradict what the
   source paper itself reports about its method — recorded here the same
   way the [root-CI coverage finding](#root-inference) below was.

**Result**: 13/13 assertions pass (`devtools::test()`, full suite: 306
passed, 0 failed). No implementation bug found — the only fix needed was to
the test's synthetic data-generating process, not to `dating_pdc.R` itself.

Replication script:
[replication/dating-and-root-inference/radf_pdc_validation.R](#script-radf_pdc_validation).

---

## Root inference

**Status: done.** Guo, Sun & Wang's normal-t CI + doubling time, the
Phillips-Magdalinos Cauchy CI, and per-episode datestamp integration are
all implemented in `exuber/R/rootstamp.R` as a single S3 generic,
`rootstamp()` (`default` method for a single sub-sample, `radf_obj` method
for every `datestamp()` episode at once), tested in
`exuber/tests/testthat/test-rootstamp.R`. **Update (2026-08-18):**
originally shipped as three separate functions (`explosive_root()`,
`root_ci()`, `root_ci_datestamp()`); consolidated into `rootstamp()` before
release once the three-call handoff (and `root_ci_datestamp()`'s
backwards-reading name) turned out to be UX friction worth fixing early,
not a stable public API worth preserving. Historical mentions of the three
old names below are left as-is (they describe what was true at each dated
entry); read them as `rootstamp()` in current terms.

### Source

- Phillips, P. C. B., & Magdalinos, T. (2007). Limit theory for moderate
  deviations from a unit root. Journal of Econometrics, 136(1), 115-130.
  Working paper: Cowles Foundation DP 1471 (July 2004), **open**,
  `cowles.yale.edu/sites/default/files/2022-08/d1471.pdf` — not paywalled,
  just not found in the original search pass; Theorem 4.3 verified by PNG
  page-render 2026-08-09.
- Guo, G., Sun, Y., & Wang, S. (2019). Testing for moderate explosiveness.
  The Econometrics Journal, 22(3), 279-303. Originally not directly
  accessed (paywalled); recovered via institutional access. —
  `root_ci()` was implemented from a secondary restatement before this was
  found (see below); worth cross-checking against the primary source now
  that it's available.

The originally-implemented result is sourced from a secondary review that
restates both papers' results: Skrobotov, A. (2023) — "Testing for
explosive bubbles: a review", open on arXiv (`arxiv.org/pdf/2207.08249`) —
the same review this whole project started from.

### What's implemented, and why only part of it

The review states two distinct results:

1. **Phillips-Magdalinos (2007a)**: for a mildly explosive AR(1)
   \eqn{\rho_T = 1 + c/k_T}, a specific normalization of
   \eqn{(\hat\rho_T - \rho_T)} converges to a **standard Cauchy**
   distribution — remarkable because it holds under non-Gaussian errors
   too. The review gives the resulting CI as
   \eqn{\hat\rho_T \pm \frac{\sqrt{\hat\rho_T^{2n}-1}}{\hat\rho_T^n} C_\xi}.
2. **Guo, Sun & Wang (2019)**: generalizing to allow a drift term, they show
   the **ordinary regression t-statistic** for \eqn{\rho_T} is
   asymptotically **standard normal** under i.i.d. errors (Student's-t/HAR
   under dependence) — a practically simpler result that doesn't need to
   know \eqn{c}, \eqn{k_T}, or the exact rate.

Only (2) was implemented at first. Reason: result (1)'s exact formula came
through `pdftotext` with subscripts/exponents visibly mangled (the review
itself is a secondary source restating a third paper's theorem, so there
was originally no primary-source PDF to render as an image and verify) —
shipping it with an unverified exponent would risk a silently wrong
confidence interval, which is worse than not shipping it. Result (2)
survived as unambiguous prose ("the t-statistic... is asymptotically
standard normal"), which was enough to implement correctly and low-risk:
`root_ci()` is just `rho_hat +/- qnorm(...) * se`, using the standard error
from an ordinary no-intercept OLS fit (`explosive_root()`) -- unusual only
in that the asymptotic justification for treating this as normal, despite
`rho > 1`, is Guo/Sun/Wang's CLT, not the classical stationary-AR one.

### Exact number verified

The review's footnote 17 gives the two-sided Cauchy percentiles used in
result (1): C_0.10 = 6.315 (≈6.314), C_0.05 = 12.7, C_0.01 = 63.65674. The
standard Cauchy distribution is identically Student's *t* with 1 degree of
freedom, a pure mathematical fact independent of whichever paper states it
— so this is checkable exactly, no simulation, no primary-source risk:
`qt(0.95, df = 1) = 6.313752`, `qt(0.975, df = 1) = 12.7062`,
`qt(0.995, df = 1) = 63.65674`, matching the footnote to the precision
given. Tested in `test-rootstamp.R`. (This doesn't validate the *exponent*
in result (1)'s CI formula, which is why that formula wasn't shipped
immediately — only that the percentile table quoted is correct.)

### What `root_ci()` gives you

Given a sub-sample (e.g. an explosive episode already identified by
`datestamp()`, converted to row positions), `explosive_root()` fits the
no-intercept AR(1) OLS regression Phillips-Magdalinos's model specifies,
and `root_ci()` reports:
- `rho`, `rho_ci`: point estimate and Wald interval,
- `doubling_time`, `doubling_time_ci`: `log(2)/log(rho)` (periods for the
  bubble to double at the estimated rate) with its own interval, obtained
  by transforming the `rho` interval's endpoints (doubling time decreases
  in `rho`, so the bounds flip).

### Empirical coverage, honestly reported

A 500-replication simulation at a modest sample size (explosive episode of
149 observations, ρ = 1.03) gave empirical coverage of the nominal 95% CI
around 90%, not 95% — finite-sample undercoverage, not a bug (verified the
formula reduces to the textbook Wald interval correctly; the CLT this
relies on is a T → ∞ result and this specific estimator is known in the
literature to converge slowly). `test-rootstamp.R`'s coverage test checks
against a deliberately loose bound (>80%) reflecting this, rather than
asserting the nominal rate — documented here rather than glossed over.

### Independent validation (2026-08-09)

Different parameters and seed from `test-rootstamp.R` throughout (that suite
uses `rho=1.03, n=150`; here `rho=1.05, n=200`), to check the finding isn't
an artifact of one specific parameterization.

**Point estimate**: on a fresh simulated explosive AR(1) (`seed=8675309`),
`explosive_root()` recovered `rho_hat = 1.0500` against `rho_true = 1.0500`
— matching to 4 decimal places, with a 95% CI so tight it's indistinguishable
from the point estimate at that precision. This is not a bug: it's the
expected "super-consistency" of explosive-root estimation — the regressor
`y_{t-1}` itself grows geometrically, so `sxx = Σy_{t-1}²` (the denominator
of the OLS standard error) explodes at rate `ρ^{2n}`, and by `n=200` with
`ρ=1.05` that denominator is astronomically large, collapsing the standard
error far faster than in the unit-root or stationary case. Consistent with
the theory this CI is built on (Guo, Sun & Wang's CLT for a genuinely
explosive root).

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

Result: **94.6%** coverage of the nominal 95% interval — a marked
improvement over the ~90% the package's own test observes at `rho=1.03,
n=150`. This is the expected direction: a larger `n` and a root further
from 1 (`1.05` vs `1.03`) both push the estimator closer to its asymptotic
(`T→∞`) regime faster, so coverage should improve toward nominal as either
increases — which is exactly what this independent run shows. Reinforces
rather than contradicts the existing "known finite-sample undercoverage"
disclosure: the undercoverage shrinks in the direction the theory predicts,
which is itself a check that the CLT justification is the right one and
not a coincidence.

**Full existing suite**: `test-rootstamp.R` — **8 passed, 0 failed** (1
skipped, CRAN-only).

**Conclusion**: no issues found. Point estimate and CI behave exactly as
the underlying asymptotic theory predicts, including the direction of
convergence as `(n, ρ)` move further from the boundary.

### Update: primary source found, exponent now verified

The Phillips-Magdalinos (2007) primary source turned out to be open after
all — see Source above. `pdftotext` mangles the formula the same way it
mangled SBZ's Table 1 (exponents/fraction bars merge into runs like
"nnn"), so page 14 was rendered to a PNG with PyMuPDF (same method as
STADF/SBZ) and read directly. **Theorem 4.3** (their eq. 26), for the
moderate-deviations model \eqn{\rho_n = 1 + c/n^\alpha}, \eqn{c>0},
\eqn{\alpha \in (0,1)}:

\deqn{\frac{n^\alpha \rho_n^n}{2c}(\hat\rho_n - \rho_n) \Rightarrow C}

with `C` standard Cauchy. Remark (i)/eq. 27 also gives the simpler
**fixed-root exact-explosive case** (White 1958, restated here), which
needs no \eqn{\alpha} or \eqn{c} and is the more directly usable form for a
plug-in CI:

\deqn{\frac{\rho^n}{\rho^2-1}(\hat\rho_n - \rho) \Rightarrow C}

i.e., plugging in \eqn{\hat\rho} for the unknown \eqn{\rho} on the
normalization (standard practice for this kind of self-normalized pivot),
a two-sided CI is \eqn{\hat\rho \pm q_{\alpha/2} \cdot
(\hat\rho^2-1)/\hat\rho^n} where \eqn{q_{\alpha/2}} is a standard-Cauchy
quantile (`qcauchy()`, or equivalently `qt(., df = 1)`, already verified
exactly above). This is a **different, simpler formula than the review's
garbled restatement** quoted earlier (`sqrt(rho^(2n)-1)/rho^n`) — that
restatement should be treated as wrong/superseded now that the primary
source has been read directly, not as an equivalent alternative form.

**Update (2026-08-09, Bundle 1): implemented.** `root_ci(x, type = "cauchy")`
now ships this as a second CI type alongside the default `"normal"` one,
using the fixed-root form (eq. 27) derived above. Roxygen docs note the
Cauchy interval assumes a *fixed* explosive root while the default
normal-t interval (Guo/Sun/Wang) allows drift/dependence and stays the
safer default. Tested in `test-rootstamp.R` (brackets the point estimate;
matches the closed-form eq. 27 formula exactly, not just approximately).
The true moderate-deviations form (eq. 26, needing an estimate of
\eqn{\alpha}, the localizing-rate exponent) is a materially harder
follow-on than the fixed-root form and still isn't attempted.

Also relevant, not implemented: Phillips, Magdalinos & Giraitis (2010,
J. Econometrics 158(2), 274-279, "Smoothing local-to-moderate unit root
theory," open as Cowles DP 1659)
shows the moderate-deviations theory above smooths continuously into the
local-to-unity case as \eqn{\alpha \to 0} — background theory, not a
standalone feature, but the citation to reach for if `root_ci()` ever
needs to handle roots close to the local-to-unity boundary rather than
assuming a clearly mildly-explosive \eqn{\rho}.

### Update (2026-08-09, Bundle 1): `root_ci_datestamp()` instead of `summary()`

Wiring `explosive_root()`/`root_ci()` directly into `radf_obj`'s
`summary()` method turned out to be a worse fit than originally scoped:
`summary.radf_obj()`'s S3 dispatch (`summary_radf.mc_cv`/`.wb_cv`/`.sb_cv`)
is built entirely around `radf_cv` test-statistic critical values — root
CIs need a `datestamp()` result instead, a structurally different kind of
input with no natural slot in that dispatch chain. Restructuring shared
`summary()` machinery used by three other critical-value types to
accommodate a fundamentally different output was a bigger, riskier change
than "wire it in" suggested.

Shipped instead: `root_ci_datestamp(object, ds, level = 0.95, type =
"normal")`, a standalone function that runs `explosive_root()`/`root_ci()`
on every episode in a `datestamp()` result and returns output in the same
per-series-named-list shape `datestamp()` itself uses. Tested end-to-end
in `test-rootstamp.R` against a real `radf()` → `radf_mc_cv()` → `datestamp()`
pipeline, checked against calling `explosive_root()`/`root_ci()` directly
on the same episode. Root inference on very short episodes (duration 1-2)
is left to degrade honestly (matching what calling `explosive_root()`
directly on 2-3 points would do) rather than silently filtered — the
docs point at `datestamp()`'s existing `min_duration` argument instead of
adding a second, redundant filtering knob.

Replication script:
[replication/dating-and-root-inference/rootstamp_validation.R](#script-rootstamp_validation).

---

## Confidence sets for bubble dates

**Status: evaluated, not implemented — re-triaged 2026-08-10, finding
real structure but still a multi-step lift, not completed this pass.**
Full PDF read (abstract, intro, model). No R/C++ code written.

### Source

Kurozumi, E. & Skrobotov, A. (2025). "Confidence Sets for the Emergence,
Collapse, and Recovery Dates of a Bubble." arXiv:2511.16172.

### What it is

A CI layered on top of already-point-estimated dates — conceptually the
dating analogue of [root inference](#root-inference)'s `root_ci()`, not a
new detection or point-estimation method. It explicitly does **not** use
the obvious approach (the limiting distribution of the breakpoint
estimator itself): the paper states its own preliminary simulations found
that performs poorly for bubble dates. Instead it builds confidence sets
by *inverting* hypothesis tests for the break location — a likelihood-
ratio-type test (Eo & Morley 2015) and Elliott-Müller-type (2007) tests,
used individually and combined — deriving new limiting null and
alternative distributions for each, evaluated via Monte Carlo. Three
dates (emergence, collapse, recovery), estimated separately rather than
jointly.

### Cost/feasibility note for exuber

Not a "wrap an existing point estimate in `± z·se`" CI the way
`root_ci()` was for the explosive root (Bundle 1) — the paper had to
invent this whole route *because* the naive analogue of that approach
doesn't work here. A faithful implementation needs: multiple new test
statistics (LR-type, Elliott-Müller-type) each with their own critical
values, a rule for combining them, and per-date confidence-set
construction (invert-a-test-statistic, not point-estimate-plus-margin)
repeated for three separate dates. Its own precondition — a
WLS/volatility-corrected `dating_pdc()` variant — is now shipped (see
[below](#wls-dating-under-time-varying-volatility)), which prompted a
re-triage.

**Re-triaged (2026-08-10), re-reading rendered pages 8, 10-11, 46-47.**
The critical-value machinery is genuinely cheaper than the original
"needs new simulation for everything" framing suggested — but the
statistic construction itself confirms the original "multi-step, HLS
/HLW-scale" cost estimate, not a hidden quick win:

- **Favorable**: not every statistic needs simulated critical values.
  `LR^e_{a,12}`/`EM^e_{a,12}` (their eq. 15-16, the "12"-direction tests)
  have an **exact closed-form chi-square critical value** —
  `cv^e_{LR12,0.05} = λ1·χ²_{1,0.05}`,
  `cv^e_{EM12,0.05} = sqrt(λ1·χ²_{1,0.05})` — no simulation at all. The
  "21"-direction tests (`LR^e_{a,21}`, `EM^e_{a,21}`, `EM^e_{b,21}`, eq.
  17-19) don't have a closed form, but the paper *publishes* a
  response-surface regression for their critical values (their eq.
  above Table 1, `cv = a_{0,ℓ} + a_{-1,ℓ}/λ1* + a_{1,ℓ}·λ1* +
  a_{2,ℓ}·λ1*² + a_{3,ℓ}·λ1*³`, coefficients transcribed from their
  Table 1) — a MacKinnon-style formula, not a table to interpolate or a
  simulation to run, the same "published, not new" pattern that made
  Kurozumi (2020)'s and HB's own boundaries cheap elsewhere in this
  project.
- **Still genuinely bigger**: the statistics themselves are not a thin
  reuse of what's already shipped. The paper's own actually-recommended
  test (their "`LE^e` test," combining `LR^e_{b,12}` with `EM^e_{a,21}`,
  chosen over the naive `LR^e_{a,12}` specifically because that one is
  "over-sized in finite samples") needs: (a) `LR^e_{b,12}` (their eq.
  11) — a `min` over candidate break dates of `(y²_{T2} - ρ̂_a ·
  Σ_{t=T1+1}^{T2} y²_{t-1}) / (T·φ�̂_a^{2(T2-T1)}·σ̂²/2)` — the numerator
  *is* the familiar prefix-sum-window pattern (`Σy²_{t-1}` via
  cumulative sums, exactly `dating_hls()`'s own `hls_prefix_sums()`
  construction), but the nuisance-parameter estimators `ρ̂_a`, `φ̂_a`,
  `σ̂²` and the admissible-break-date set `Λ^e_{12}`'s construction rule
  were not pinned down this pass (need more of Section 2's model setup
  than was read); (b) `EM^e_{a,21}` (their eq. 18) is an *integral* over
  a continuum of candidate break points of an `ADF(λ2*, λ1*)` functional
  (eq. 18's own definition, a ratio of Brownian-motion-type functionals)
  — genuinely new estimation machinery with no exuber analogue, not a
  discrete min/max search the way `LR^e_{b,12}` is.

Net: this is no longer accurately described as "new estimation *and*
new critical-value-simulation work" — the critical-value half is now
known to be cheap (closed-form or published response-surface, confirmed
by direct formula transcription, not assumed). But the statistic-
construction half is confirmed, not merely assumed, to need genuinely
new machinery for at least one of the two components the paper's own
recommended test combines, on a scale comparable to
[HLS](#ssrbic-dating-vs-psy-recursive-dating)/HLW rather than a same-day
addition — and that's before repeating any of it for the collapse and
recovery dates. Not completed this pass; flagged precisely enough that a
future pass can start from "`LR^e_{b,12}`'s SSR-style numerator is
prefix-sum-ready, go find `ρ̂_a`/`φ̂_a`/`σ̂²`/`Λ^e_{12}` in Section 2, then
tackle `EM^e_{a,21}`'s integral separately" rather than re-deriving
Table 1's meaning from scratch.

---

## Improved retrospective dating

**Status: single-bubble omission fix done (2026-08-10); the
multi-bubble Bai-Perron/Perron-Qu dynamic-programming algorithm not
implemented.** Full PDF read (abstract, intro, model, Theorems 1-2 —
re-verified against rendered PDF pages 3-4 for eq. 1-8 — the
HLS-equivalence footnote, and Section 3's DP algorithm description).

### Source

Kejriwal, M., Nguyen, L. & Perron, P. (2025). "An Improved Procedure for
Retrospectively Dating the Emergence and Collapse of Bubbles." *JTSA*,
46(5), 867-883. `doi:10.1111/jtsa.12810`.

### What it is

A different bias fix for the *same* joint-SSR family as HLS, not a
relative of PDC/KS. The model uses HLS's **fixed** autoregressive
coefficient framework (`rho` fixed, not Phillips-Magdalinos "mildly
explosive" `rho_T -> 1`) with an abrupt (not stationary-transition)
collapse. Theorem 1 shows the standard OLS joint-SSR estimator is
inconsistent: the origination-date estimate converges to the *collapse*
date, and the collapse-date estimate converges to a date *after* the true
collapse, offset by the trimming parameter — both biased late, and for a
different reason than PSY's own threshold-crossing delay bias. The fix
(Theorem 2) is a "modified SSR" that **omits the single residual at the
implosion date** from the objective function, which restores consistency
for both dates. A footnote establishes this omission is numerically
*equivalent* to a specific one-time-dummy modification of HLS's own Model
4 regression — i.e., this paper is best read as "a bias-fix inside the
HLS estimating equation," not a third independent dating family. Section
3 also develops a Bai-Perron/Perron-Qu-style dynamic-programming
algorithm so the multi-bubble case avoids HLS's brute-force combinatorial
grid search, exploiting the unit-root restriction on the non-bubble
regimes to skip the iterative-initial-values step that ordinary
Bai-Perron/Perron-Qu DP would need.

**A structural finding, confirmed by re-deriving the closed forms
directly**: KNP's single-bubble model (their eq. 1-3) is not a generic
"intercept+AR(1)-dummy" structure needing new estimation code — it is
*exactly* HLS's own Model 2 shape (unfitted unit root, then an
intercept+slope-fitted explosive regime, then an unfitted unit root
resuming after an instantaneous collapse). Their eq. 2's `delta_hat`
formula is algebraically the standard bivariate OLS slope (the
`(x-xbar)*ybar` cross term vanishes by construction), and regressing the
*level* `y_t` on `y_{t-1}` with an intercept gives identical residuals/SSR
to regressing `Delta y_t` on `y_{t-1}` with an intercept (a fixed
reparameterization, `slope' = delta - 1`) — the same regression
`dating_hls.R`'s `hls_segment_ssr()` already computes for HLS's own Model 2
search. The entire fix reduces to `SSR_om(T1,T2) = SSR(T1,T2) -
(Delta y_{T2+1})^2`: one already-computed squared term subtracted from an
already-computed SSR, no new regression at all.

### Implementation

Shipped as `dating_knp(data, trim = 0.05, omit = TRUE)` in
`exuber/R/dating_knp.R`, tested in `exuber/tests/testthat/test-knp.R`.
Internal helper `knp_find_break()` reuses `dating_hls.R`'s
`hls_prefix_sums()`/`hls_segment_ssr()` directly (no new closed-form
derivation needed, per the structural finding above), jointly searching
`(tau1, tau2)` to minimise the omission-corrected SSR — or, with
`omit = FALSE`, the plain (provably inconsistent) SSR, kept specifically
so the correction's effect can be demonstrated and tested directly
rather than only asserted. Unlike `hls_model23()`, KNP's own candidate
set imposes no directional sign constraint on the fitted "peak".

**Validation, including a direct reproduction of the paper's own
theorems, not just a plausibility check**: formula-exact — both the
`omit = FALSE` and `omit = TRUE` searches match an exhaustive brute-force
nested-`lm()` search exactly. Monte Carlo on KNP's own DGP (unit root,
no-intercept explosive AR(1), instantaneous collapse back near the
pre-bubble level, fresh unit root; 30 reps, `T1=50`, `T2=90`, `T=200`,
`delta=1.05`):

- **Theorem 1 (naive, `omit = FALSE`) reproduced directly**: the
  origination-date estimate's mean mistake relative to the *true
  collapse date* (`mean|tau1_hat - T2| = 1.0`) is dramatically smaller
  than its mistake relative to its *own* true origination date
  (`mean|tau1_hat - T1| = 39.0`) — i.e., the naive estimator's `tau1_hat`
  essentially tracks the wrong date, converging to `T2` rather than `T1`,
  exactly as Theorem 1 predicts, not merely "somewhat biased."
- **Theorem 2 (omission-corrected, `omit = TRUE`) reproduced directly**:
  the same bias falls from 39.0 to 13.0 observations (a genuine
  correction, not a full elimination at this finite `T` — Theorem 2 is
  an asymptotic `→p` result, so residual finite-sample bias at `T=200`
  is expected, not a defect); the collapse date and the explosive
  coefficient `delta_hat` are both close to their true values
  (`mean|tau2_hat - T2| = 1.0`; `delta_hat` mean `0.986` vs. true `1.05`).

Replication script:
[replication/dating-and-root-inference/radf_knp_validation.R](#script-radf_knp_validation).

**Not implemented**: Section 3's Bai-Perron/Perron-Qu-style dynamic
programming algorithm for the multi-bubble case. This is genuinely new
algorithmic machinery exuber has no analogue of (PDC/KS's sequential
`O(T)` scan and HLS's own per-episode grid search are both much simpler
than a DP over an unknown number of breaks) — scoped out as its own
follow-on, the same decision already made for HLW's multi-bubble
fragmentation-joining heuristic elsewhere in this file.

---

## WLS dating under time-varying volatility

**Status: done (2026-08-09).** Shipped as `dating_pdc(..., type = "wls")`.

### Source

Kurozumi, E. & Skrobotov, A. (2023). "Improving the accuracy of bubble date
estimators under time-varying volatility." arXiv:2306.02977.

### What it is

A direct two-step generalization of PDC/KS's own sequential dating
estimator — same authors as the [KS (2023) 4-regime extension](#3-pdc-2021-journal--ks-2023-journal--sequential-sample-splitting)
`dating_pdc()` already implements, explicitly built on top of it ("we
estimate these break dates as proposed by PDC and Kurozumi and Skrobotov
(2022) and collect the residuals..."). Step 1 is exactly
`pdc_find_break()`/`dating_pdc()` as already implemented: fit the
homoskedastic no-intercept AR(1) break model, get consistent break-date
(fraction) estimates and their residuals. Step 2 nonparametrically
estimates the time-varying error variance `sigma_t^2` from those
residuals, then re-estimates each break date by minimizing a **weighted**
SSR — `sum(y_t - a*y_{t-1})^2 / sigma_t^2` instead of the unweighted sum.
Algebraically this is `pdc_find_break()`'s cumulative-sum trick with every
`y_t`, `y_t*y_{t-1}`, `y_{t-1}^2` term divided by `sigma_t^2` before the
prefix sum — still closed-form, still `O(T)` per breakpoint given the
volatility weights.

### Implementation

Shipped as `dating_pdc(data, ..., type = c("ols", "wls"))` in
`exuber/R/dating_pdc.R`. `type = "ols"` is the original PDC/KS estimator,
unchanged. `type = "wls"`:

1. Runs the existing sequential `type = "ols"` fit to get step-1 break
   estimates.
2. `pdc_regime_resid(y, breaks)` (new) computes the fitted no-intercept
   AR(1) residual at every `(y_{t-1}, y_t)` pair, one OLS `rho` per
   regime implied by the step-1 breaks — "collect the residuals of the
   fitted [regime] model," per the paper.
3. Those residuals feed `nw_spot_vol()` — the Nadaraya-Watson kernel
   smoother with leave-one-out-cross-validated bandwidth that exuber
   already had, **extracted** (not duplicated) from
   `radf_sbz.R`'s `kernel_spot_vol()` so both SBZ (smooths squared first
   differences) and this estimator (smooths squared regime residuals)
   share one implementation. `kernel_spot_vol(y) := nw_spot_vol(diff(y))`
   is now a one-line wrapper around it; behavior is unchanged (verified
   by full-suite regression: `test-sbz.R` still green after the
   extraction).
4. `pdc_find_break()` gained an optional `weights` argument (`NULL` =
   original unweighted behavior, verified identical by a dedicated
   regression test) — every cumulative sum is multiplied by the weight
   vector before the prefix sum, so the search stays the same `O(T)`
   closed-form scan.
5. The full sequential (collapse, then origination, then recovery) search
   is re-run once more with `weights = 1 / sigma_t^2`, using the correct
   contiguous slice of the full-sample `sigma_t^2` vector for each
   sub-sample regression.

No new critical-value simulation needed — like the OLS version, this is
point estimation, not a threshold-crossing test.

**Independent validation**: two Monte Carlo checks (40 seeds each, not a
single cherry-picked run), matching this project's usual bar of
distinguishing "the estimator works" from "it happened to pass once":

- **Homoskedastic DGP** (no volatility signal to exploit): OLS and WLS
  origination-date MAE are statistically indistinguishable (5.42 vs
  5.53) — WLS costs essentially nothing when there is nothing to gain,
  as expected from a nonparametrically-weighted estimator with no true
  heteroskedasticity to detect.
- **Heteroskedastic DGP with a volatility burst in the first 20% of the
  pre-bubble regime** (the specific scenario the paper's own Monte Carlo
  reports the largest gains for): origination-date MAE drops from
  **13.05 (OLS) to 2.33 (WLS)**, a ~5.6x improvement, confirming the
  mechanism actually works as claimed — OLS's unweighted objective lets
  the noisy early segment dominate the origination split; WLS correctly
  downweights it via the estimated spot variance. Collapse-date accuracy
  is unaffected either way (both already near-exact, since the
  explosive-to-collapse transition dominates the SSR regardless of
  earlier noise — consistent with PDC's own stochastic-order argument
  for why collapse identifies first).

Regression-tested against the loose 2x version of this margin in
`test-pdc.R` so a future change that erodes the benefit gets caught
without being brittle to the exact numbers on a different RNG/BLAS.

Replication scripts:
[replication/dating-and-root-inference/radf_pdc_wls_heteroskedastic_mae.R](#script-radf_pdc_wls_heteroskedastic_mae),
[radf_pdc_wls_homoskedastic_mae.R](#script-radf_pdc_wls_homoskedastic_mae).

---

## Reverse-regression recovery dating

**Status: implemented (2026-08-10), shipped with honest caveats —
`f_r` (recovery date) validates well, `f_c` (crisis-origination date) and
the overall false-detection rate under the null are noisier than hoped
and not fully resolved.** Shipped as `radf_recovery()`/
`radf_recovery_cv()`. Full PDF read (abstract through Section 3.2 —
model, both forward and reverse BSDF definitions, both limit-theory
theorems, the finite-sample Monte Carlo setup, and — on a second, more
careful pass — Section 4.3's real-time monitoring extension, which
corrected a mis-transcription in an earlier pass of this file, see
below).

### Source

Phillips, P.C.B. & Shi, S. (2014). "Financial Bubble Implosion." Working
paper: Cowles Foundation DP 1967. Published as "Financial Bubble Implosion
and Reverse Regression," *Econometric Theory*.

### What it is

The mechanism is, as hoped, genuinely simple: reverse the series
(`X*_t := X_{T+1-t}`), run the *same* BSDF/BSADF recursion PSY already
uses on `X*`, and map the crossing-time fractions back to the original
time index (their eqs. 8-9):

```
f_hat_r = 1 - g_hat_e,  g_hat_e = inf{ g in [g0, 1]    : BSDF_g(g0) > scv }   (recovery date)
f_hat_c = 1 - g_hat_c,  g_hat_c = inf{ g in [g_hat_e,1]: BSDF_g(g0) < scv }   (crisis-origination date)
```

**Correcting a mis-transcription from an earlier pass of this file**:
`f_hat_c` is *not* "a further/later correction after recovery" — re-
reading the extracted text directly against eqs. 8-9 (not just the
surrounding prose) shows `g_hat_c` is searched only *after* `g_hat_e`, so
`f_hat_c <= f_hat_r` always: `f_hat_c` is the ORIGINAL series' crisis/
collapse-onset date, re-derived via reverse regression as an alternative
to the collapse date PSY's forward test already dates, and `f_hat_r`
(recovery) comes chronologically *after* it. The paper states this
explicitly: "market recovery (`f_hat_r`) following a crash begins when
normal market behavior changes to exuberance in the reverse series
(`g_hat_e`)... market collapse in the original series begins when
exuberance in the reverse series shifts to collapse at (`g_hat_c`)."
Reversing a mildly-explosive-then-mildly-integrated-collapse process
(model (2)/(7) in the paper) turns the collapse regime into an explosive
regime in reverse time and vice versa, so "detect the crisis-origination
and market-recovery dates" becomes "detect explosiveness in the reversed
series" — the same right-tailed test PSY already runs, just on `rev(x)`.

**The catch, and why this isn't free**: Theorem 1 derives the null
limiting distribution of the reverse statistic, `F_g(W, g_0)`, and it is
**not** the same distribution as the forward statistic's `F_f(W, f_0)` —
it has an extra term arising because reversing a random walk makes the
reversed "lagged" regressor correlated with the reversed current error
(`E[X*_{T-j+2} * eps_{T-j+2}] != 0`, stated explicitly below Theorem 1).
This isn't a bubble-specific effect — it's a generic consequence of
running a regression on time-reversed data — so it affects the reverse
test's critical values *even under the null*, not just its finite-sample
power. **Verified empirically, not just asymptotically**: a paired Monte
Carlo (same underlying draws, `n=100`, `minw=20`, 5,000 reps) comparing
`radf_mc_cv()`'s standard forward critical values against the same
recursive computation run on the reversed path gives measurably
different values (mean absolute difference ~0.04, max ~0.11 at the 95%
level across positions) — consistent with, not merely assumed from, the
paper's own theorem.

**A separate mechanism, not this one — Section 4.3's real-time
monitoring extension** (their eqs. 10-11): apply the same reverse-
regression machinery repeatedly on a *growing* sample from the collapse
date `T_c` forward, stopping at the first sample-end `K` for which a
correction is detected. This is what actually produces the paper's
"further correction in January 2004... full return to normal market
conditions in May 2004" figure (their dot-com empirical application,
Section 5) — a *separate, sequential* application starting from the
peak, not part of the eqs. 8-9 `f_c`/`f_r` pair. **Not implemented here**;
a natural, cheap follow-on given `radf_recovery()` already exists (it
would be an outer loop expanding the sample and re-calling
`radf_recovery()`-style logic each step, structurally similar to
`radf_monitor()`), but explicitly out of scope for this pass.

**Concrete empirical numbers from the paper's own dot-com illustration**
(Section 5, prose, NASDAQ price-dividend ratio): the eqs. 8-9 pair gives
crash March-November 2000 (so `f_c` = March 2000, `f_r` = November
2000); the *separate* Section 4.3 monitoring extension then gives the
further correction January-May 2004. Not independently reproduced here
(no access to the underlying NASDAQ price-dividend series).

### Implementation

Shipped as `radf_recovery()` (main function) and `radf_recovery_cv()`
(reversal-calibrated Monte Carlo critical values, mirroring
`radf_mc_cv()`'s own simulate-then-quantile construction — including its
`cummax(badf)`-as-`bsadf`-boundary shortcut — with one added `rev()`
before the recursive computation, rather than reusing forward critical
values as an approximation of unknown quality). Mechanically as cheap as
hoped: `radf()`'s existing `bsadf` recursion run on the reversed series,
compared against the reversal-calibrated boundary, with the first up-
crossing/down-crossing pair mapped back via `f = n + 1 - g`. No C++, no
new point statistic.

**Validation, reported honestly**: the structural invariant `f_c <= f_r`
holds by construction whenever both dates are identified and uncensored
(confirmed across all replications tried). `f_r`'s bias against synthetic
collapse-then-recovery data is small (a few observations, same direction
and rough magnitude as the paper's own Table 5 finding of ~6 observations
early). `f_c`'s bias is materially larger (mean |bias| approaching the
length of the synthetic collapse window itself in some runs), and the
empirical false-detection rate under a pure random-walk null (`n=100`,
`minw=20`, 95% level, one stable reversal-calibrated cv reused across 200
fresh draws) is around 29% — higher than comparable forward-test numbers
found elsewhere in this project (e.g. `radf_monitor()`'s ~10% cumulative
false-alarm rate over a 75-point horizon). One real synthetic-DGP
artifact was found and fixed during this validation (an abrupt level-jump
at the expansion-to-collapse regime boundary was producing a spurious,
narrow bsadf spike concentrated right at the junction rather than
spanning the intended collapse regime — replaced with a smooth,
continuous mean-reverting transition); fixing it improved `f_r`'s bias
substantially but did not resolve `f_c`'s. A plausible, non-alarming
explanation was identified but not confirmed: eq. 9's `inf` operator
defines the *first* down-crossing with no persistence requirement, so a
transient noise-driven dip below the reversal-calibrated boundary is
enough to trigger a premature `f_c` — this would be a genuine property
of the paper's own literal crossing rule under finite-sample noise, not
necessarily an implementation bug, but this has not been fully
distinguished from a subtler code issue. Replication script:
[replication/dating-and-root-inference/radf_recovery_validation.R](#script-radf_recovery_validation).

Shipped anyway (user's explicit call after being presented the tradeoff)
because the mechanical/structural parts are sound and independently
verified (the invariant, the differing-CV finding, `f_r`'s accuracy), and
because further debugging the residual `f_c`/false-detection concern is
better scoped as its own follow-up than as an open-ended extension of
this pass. `radf_recovery()`'s own roxygen docs carry the same caveat
inline.
