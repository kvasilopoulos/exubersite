---
title: "Simulation DGPs"
blurb: "Data-generating processes for the axes the original sim_*() functions do not cover."
order: 6
---
Not a test/statistic taxonomy file like the others — this one catalogues
**data-generating processes** (DGPs) for simulating bubble series, as used
in Monte Carlo sections across the literature this project has already
read, plus a few papers pulled in specifically for this pass. The question
this file answers: what does exuber's own `sim.R` (`sim_psy1`, `sim_psy2`,
`sim_ps1`, `sim_ps2`, `sim_blan`, `sim_evans`, `sim_div`) *not* cover?

All seven existing `sim_*` functions share one structural property: fixed
(homoskedastic), i.i.d. Gaussian innovations. They differ in the mean
equation (PSY-style regime-switching AR(1) vs. Blanchard/Evans rational
bubbles) and in bubble count/collapse shape, but none has time-varying
volatility, GARCH, non-Gaussian innovations, stochastic switch timing, or a
multi-series/factor structure. Every item below is flagged only if it
differs on one of those axes from an actual, genuinely different
mechanism — not a reparameterization of an equation already implemented
(see the exclusions at the end).

Status legend as in [volatility-robustness.md](/replication/volatility-robustness).

**Status: all 15 catalogued DGPs implemented (2026-08-11)**, in
`exuber/R/sim.R`, tested in `exuber/tests/testthat/test-sim-dgp.R` (39
assertions: formula-exact checks against independent brute-force
reimplementations where feasible, published-property checks — e.g.
`sim_tree()`'s price floor, `sim_msbubble()`'s empirical transition-rate
match — moment/reproducibility checks otherwise). Two axes (#1/#2/#3/#8/#9
/#10/#11/#12, stochastic-volatility/innovation-distribution/level-shift/
long-memory/stochastic-coefficient) turned out to compose cleanly as new
optional arguments on the *existing* `sim_psy1()` (`e`, `shifts`,
`coef_noise`) plus five small, independently-reusable generator functions
(`sim_innov()`, `sim_vol_garch()`, `sim_vol_cir()`, `sim_vol_sv()`,
`sim_fi()`) — not thirteen near-duplicate clones of `sim_psy1`'s loop.
`sim_blan()` gained a `type = "rotermann_wilfling"` option (#13). The
remaining six architecturally-distinct DGPs (#4-7, #14, #15) shipped as
standalone functions: `sim_common()`, `sim_coexplosive()`, `sim_tree()`,
`sim_mar()`, `sim_msbubble()`, `sim_falsebubble()`.

**Important scope note, unchanged from the original survey**: this is DGP
machinery only. Building a DGP does not imply the *paired test* is
implemented. Update 2026-08-11: HLTZ (2025)'s level-shift-robustness
result (§3's own motivation for DGP #3) is now implemented — see
[volatility-robustness.md](/replication/volatility-robustness#level-shift-robustness-hltz-2025)
for `radf_sign()`'s level-shift robustness and the new `radf_sign_dm()`.
`sim_falsebubble()`/`sim_msbubble()` still have no dedicated exuber test
at all (they're stress-test/demo series for the *existing* PSY/GSADF
machinery, same role `sim_evans()` already plays).
See each item below for what, if anything, it's meant to validate against.

**A real bug found and fixed during validation, not a design choice**: the
first `sim_fi()` implementation used `stats::filter(..., sides = 1)` for
the truncated MA(∞) convolution, with the innovation vector the same
length as the filter itself — `sides = 1` filtering requires strictly
*more* input than filter taps to produce any non-`NA` output, so this
returned `NA` for the entire series bar one point. Fixed by generating `m`
extra (truncation-lag) innovations beyond the requested length and
convolving by hand (`vapply()`) rather than via `stats::filter()`/
`stats::convolve()` — both of which also turned out to segfault
unconditionally on this project's R installation (Windows, R 4.6.1),
independent of exuber entirely (reproduces on `stats::filter(rnorm(120),
rep(0.01, 21), method = "convolution")` with no package loaded), so the
hand-rolled convolution is also a portability fix, not just a correctness
one. A second real bug: `sim_tree()`'s `p_t = Φ(X_t)` can round to exactly
0 or 1 in a long series' tail (floating-point underflow), which silently
turned `ξ_1t`/`ε_t` into `0/0` (`NaN`) and then propagated `NaN` through
every subsequent `y[t]` — fixed by clipping `p_t` into `[1e-10, 1 - 1e-10]`.

## Already-collected papers with a distinct DGP

| # | DGP | Source | Single/multi | Implementation |
|---|---|---|---|---|
| 1 | CIR-type stochastic volatility | Harvey, Leybourne & Zu (2019) | single | `sim_vol_cir()` |
| 2 | AR(1) lognormal stochastic volatility, near-unit persistence | Sarkar & Wells (2025/2026) | single | `sim_vol_sv()` |
| 3 | Deterministic level shifts (mean jumps) | Harvey, Leybourne, Tatlow & Zu (2025) | single | `sim_psy1(..., shifts = ...)` |
| 4 | Latent common-factor bubble | Chen, Phillips & Shi (2023) | common, multi-series | `sim_common()` |
| 5 | Bivariate co-explosive linkage | Evripidou, Harvey, Leybourne & Sollis (2022) | bivariate | `sim_coexplosive()` |
| 6 | Stochastic branching-tree (random-coefficient RCA) | Gourieroux & Jasiak (2025) | single | `sim_tree()` |
| 7 | Mixed causal-noncausal AR (MAR), heavy-tailed | Blasques, Koopman, Mingoli & Telg (2025) | single, self-terminating | `sim_mar()` |
| 8 | Fractionally-integrated (long-memory) innovations | Lui, Phillips & Yu (2024) | single | `sim_fi()` |
| 9 | PSY equation with heavy-tailed/skewed innovations | Wu, Shi & Wu (2025) | single | `sim_innov(dist = "t"/"skew_t")` |
| 10 | GARCH(1,1) / logistic smooth-transition volatility | Whitehouse, Harvey & Leybourne (2025); Harvey, Leybourne, Taylor & Zu (2024) | single | `sim_vol_garch()` |
| 11 | Stochastic explosive coefficient (persistence itself random) | Kurozumi & Nishi (2025) | single | `sim_psy1(..., coef_noise = ...)` |

### 1. CIR-type stochastic volatility

Harvey, Leybourne & Zu (2019)
robustness design beyond their main deterministic-volatility one:

```
dσ²(r) = 0.03(0.25 − σ²(r))dr + 0.1σ(r)dB(r)
```

A square-root (Cox-Ingersoll-Ross) diffusion for volatility "representative
of Bollerslev and Zhou (2002)," simulated each replication via NIID(0,1)
Brownian increments, independent of the level innovations, feeding the same
PSY-style level process `sim_psy1` already implements. What's new relative
to `sim.R`: continuous-time *stochastic* (not fixed/deterministic) volatility.

### 2. AR(1) lognormal stochastic volatility

Sarkar & Wells (2025, eq. in §2)
reused in Sarkar & Wells (2026, eq. 7):

```
y_t = ρ_n y_{t-1} + u_t,   u_t = σ_t ε_t,   log σ_t² = φ_n log σ_{t-1}² + η_t
```

with both `ρ_n → 1` and `φ_n → 1` — "double local-to-unity" in mean *and*
variance. Discrete-time stochastic volatility with persistent, near-integrated
log-variance, applied to a mildly-explosive alternative — no `sim_*`
function has a persistent (as opposed to i.i.d. or deterministic) variance
process.

### 3. Deterministic level shifts

Harvey, Leybourne, Tatlow & Zu (2025):

```
y_t = y_{t-1} + Σ_{j=1}^{m} δ_j·1(t≥τ_j) + ε_t     (null)
```

An arbitrary number `m` of jump discontinuities of magnitude `δ_j` at
unknown dates `τ_j`, with an explosive regime layered on top for the
alternative (§5.2). What's new: structural breaks in the *mean level* —
nothing in `sim.R` has jump components at all, only smooth regime-switching.

### 4. Latent common-factor bubble

Chen, Phillips & Shi (2023)
eq. 2.3/2.7-2.9:

```
X_t = Λf_t + e_t
```

`N` observed series, one latent factor `f_t` following a PSY-style
unit-root→explosive→collapse sequence, loadings `Λ ~ U[0,2]`, idiosyncratic
noise `σ_e = 0.1`; collapse-splicing construction (paper lines 960-977)
avoids discontinuities at the regime boundary. What's new: a factor-model
DGP where one common bubble drives many series jointly — every existing
`sim_*` function is single-series.

### 5. Bivariate co-explosive linkage

Evripidou, Harvey, Leybourne & Sollis (2022):

`x_t` generated from HLS's own regime-dummy Models 1-4 (individually
PSY-style); a second series is linked to it:

```
y_t = μ_y + φ_x·x_{t-i} + φ_z·z_t + ε_{y,t}
```

lead/lagged copy of the explosive series `x_t` plus a third latent
explosive series `z_t`, with heteroskedasticity timed to the regime
changes. What's new: a two-series lead/lag linkage of explosive
components, not a single univariate path — different from the common-factor
case above (#4), which shares one factor rather than linking two distinct
explosive series.

### 6. Stochastic branching-tree bubble

Gourieroux & Jasiak (2025):
an affine autoregression with a *stochastic* (not fixed) coefficient,
representable as a random-coefficient AR process generated by a binomial
tree with stochastic branching intensity (vs. Cox-Ross-Rubinstein's
deterministic branches). Blanchard & Watson (1982) — already `sim_blan` —
is the special case of *constant* intensity. What's new: the branching/tree
mechanism generates the path directly, no fixed collapse probability.

### 7. Mixed causal-noncausal AR (MAR)

Blasques, Koopman, Mingoli & Telg (2025):

```
(1 − φ1 L)(1 − ψ1 L⁻¹) y_t = ε_t
```

`φ1 = 0.7` causal, `ψ1` noncausal root, innovations Cauchy or Student-t(2).
The noncausal component autonomously generates transient, self-terminating
local bubbles — no scripted regime dates at all. What's new: different on
every axis — lag-polynomial form instead of regime-dummy AR, an implicit
(not scripted) collapse, and heavy-tailed non-Gaussian innovations.

### 8. Fractionally-integrated innovations

Lui, Phillips & Yu (2024):

```
y_t = y_{t-1} + u_t,   u_t = Δ⁻ᵈ ε_t   (d > 0, ε_t iid with finite (2+δ) moments)
```

Innovations themselves are `FI(d)` (long-memory), not i.i.d.; an explosive-
alternative analogue is developed in §4. What's new: non-i.i.d.,
long-range-dependent noise feeding the unit-root/explosive equation.

### 9. PSY equation with heavy-tailed/skewed innovations

Wu, Shi & Wu (2025)
eq. 6: structurally the standard PSY unit-root→explosive→collapse equation,
but innovations drawn from `N(0,1)`, `t(3)`, skewed-`t(3,−0.75)`,
skewed-`t(3,+0.75)`. Flagged only for the innovation distribution — the
mean equation itself is not new — since every existing `sim_*` function
uses fixed Gaussian noise.

### 10. GARCH(1,1) / smooth-transition volatility

Whitehouse, Harvey & Leybourne (2025):

```
z_t = h_t^{1/2}·ε_t,   h_t = 0.1 + 0.1z²_{t-1} + 0.8h_{t-1}
```

same GARCH(1,1) spec reused in Harvey, Leybourne, Taylor & Zu (2024)
alongside their main design, a logistic smooth-transition volatility
function:

```
σ(r) = σ1 + (σ2 − σ1)/(1 + exp{−κ(r−δ)})
```

What's new: proper conditional GARCH heteroskedasticity, and a smooth
(not instant-jump) transition between volatility regimes.

### 11. Stochastic explosive coefficient

Kurozumi & Nishi (2025)
already documented in [volatility-robustness.md](/replication/volatility-robustness#implementation--ssu-done-2026-08-10):
`1 + c1/T + a·u_t/√T` in place of the deterministic `1+c/Tᵅ` — the
persistence parameter itself is random, not just the noise scale. Listed
here for completeness (it meets this file's bar), not new information.

## New papers pulled in for this pass (not previously in the paper library)

Found by searching beyond the 48 already-collected papers specifically for
DGPs; downloaded and text-extracted the same way as the rest of this
project's library (`pdftotext -layout`, PNG-render fallback for garbled
formulas).

### 12. TGARCH(1,1) with leverage effect

Monschang, V. & Wilfling, B. (2021). "Sup-ADF-style bubble-detection
methods under test." *Empirical Economics*, 61, 145-172.
`doi:10.1007/s00181-020-01859-7`. Open access (CQE Working Paper 78/2019):

```
ε_t = s_t·h_t^{1/2},   h_t = ω + α·ε²_{t-1} + β·h_{t-1} + γ·ε²_{t-1}·1(ε_{t-1}<0)
```

calibrated to NASDAQ estimates (`α=.4387, γ_indicator=.1306, β=.9319`).
What's new relative to #10 above: *asymmetric* (sign-dependent)
shock-response (leverage effect), not symmetric/smooth-transition variance.

**Implemented as `sim_vol_garch(omega, alpha, beta, gamma)`** — the same
function as #10, with `gamma` as the shared leverage parameter
(`sim_vol_garch(omega = 0.4387, alpha = 0, beta = 0.9319, gamma = 0.1306)`
reproduces this paper's NASDAQ calibration exactly; `gamma = 0`, the
default, is plain GARCH(1,1) / item #10).

### 13. Lognormal-mixture rational bubble (Rotermann-Wilfling)

Same paper, eq. 4:

```
B_{t+1} = (B_t·u_t/δ)                    w.p. π
        = ((1−πδ)/(1−π))·B_t·u_t          w.p. 1−π,   u_t ~ iid lognormal
```

Produces periodically-recurring, *stochastically deflating* trajectories —
not a single full collapse to a fixed floor like `sim_blan`/`sim_evans`.
What's new: a different collapse mechanism from either existing rational-
bubble DGP (partial, probabilistic deflation vs. total collapse to noise).

**Implemented as `sim_blan(type = "rotermann_wilfling", delta, rw_sigma)`**
— a new branch on the existing function rather than a standalone one,
since it shares `sim_blan`'s two-regime-with-probability-`pi` structure
and only the per-regime update rule differs. Verified to stay strictly
positive (a structural invariant of the multiplicative recursion) in
`test-sim-dgp.R`.

### 14. Markov-switching present-value bubble

Chan, J.C.C. & Santi, C. (2021). "Speculative Bubbles in Present-Value
Models: A Bayesian Markov-Switching State Space Approach." *J. Economic
Dynamics and Control*, 127, 104101. Open access (author's site):

```
b_t = (1/λ_{S_t+1})·b_{t-1} + ε_bt,   S_t ∈ {1,2} first-order Markov, transition probs p11/p22
```

Regime 1 "surviving" (`λ<1`, explosive), regime 2 "collapsing" (`λ>1`,
mean-reverting), embedded in a full present-value state-space model with
time-varying expected returns/dividend growth; paper's own §"Simulated
Datasets" (Table 2) generates artificial data from these parameters. What's
new: switch *timing* is itself stochastic (Markov chain), not a
deterministic break fraction — every PSY-style DGP in `sim.R` and above
scripts `te`/`tf` as fixed fractions of `n`.

**Implemented as `sim_msbubble(p11, p22, lambda1, lambda2, sigma_b)`** —
just the bubble component `b_t`, not the full present-value/state-space
model (the paper's expected-returns/dividend-growth machinery is
orthogonal to what a stress-test DGP needs). One indexing simplification:
the source's own eq. 16 applies the regime coefficient via `S_{t+1}`
(next period's realized regime); this implementation uses the
contemporaneous `S_t` instead, which doesn't change the qualitative
Markov-switching mechanism, only which time step a given draw of `S`
labels. Verified in `test-sim-dgp.R`: the simulated regime path's empirical
self-transition rate matches `p11`/`p22` to within Monte Carlo tolerance.

### 15. Deterministic technology-adoption "false bubble" null

Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). "Technology
Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI
Episodes." arXiv:2604.25826. Open

Embeds a hump-shaped (triangular/Gaussian/Beta/Gamma) deterministic
technology-adoption shock into the Campbell-Shiller present-value
fundamental, making the *fundamental* price locally explosive during
adoption with **no bubble present at all** — engineered to make PSY-style
tests spuriously reject. Appendix E.5 extends this to a Bayesian-updated
stochastic version. What's new: this is a null (no-bubble) DGP with a
smooth deterministic drift, not a bubble-present alternative — none of
`sim.R`'s functions produce a "fundamentals-only, no bubble, still locally
explosive-looking" series.

**Implemented as `sim_falsebubble(t1, t2, kappa, shape, amplitude, mu, r)`**
— `shape = "triangular"` reproduces the paper's own eq. 4 worked example
exactly; `shape = "gaussian"` is one alternative from the paper's own list
of "any hump-shaped specification" (Beta/Gamma noted but not added — no
qualitative difference from triangular/Gaussian for this purpose, per the
paper's own robustness claim). Because the technology shock is
deterministic, its price contribution is an *exact* forward-looking
discounted sum (`T_t = Σ_{s>t} β^{s-t}τ_s`), not a further simulation
approximation. This is a simplified, single-shock reproduction of the
mechanism (deterministic hump → hump-shaped fundamental price, no bubble),
not the paper's full DOLS/multi-functional-form robustness machinery.
Verified in `test-sim-dgp.R`: `amplitude = 0` reduces to `sim_div()`'s own
fundamental-price formula bit-for-bit (a clean special-case check), and
the technology term is exactly zero outside `[t1, t2]`.

### Adjacent, not a price-level DGP

Richter, S., Wang, W. & Wu, W.B. (2023, orig. 2018). "A supreme test for
periodic explosive GARCH." *Econometrics* (MDPI); arXiv:1812.03475. Open
A piecewise/periodic explosive GARCH(1,1) where the explosiveness lives in
the *volatility recursion* (`α`, `β` temporarily driven toward/through the
IGARCH boundary `αΣ+βΣ≥1`) rather than the price level. Kept in the library
as a volatility-bubble reference, not counted above since it targets a
different detection problem (volatility bubbles, not price bubbles) than
anything `radf()`/`sim.R` addresses.

## Excluded as reparameterizations, not new DGPs

- Harvey, Leybourne & Whitehouse (2020)'s 6 Monte Carlo DGPs "A-F" (2- and
  3-bubble sequences) — same regime-dummy AR(1) family as HLS/PSY, only
  regime count and parameters differ from `sim_psy2`.
- Kurozumi & Skrobotov (2025/2026) confidence-sets paper — same
  linear-AR(1)-with-switching-coefficient form as the rest of the
  Kurozumi/Skrobotov papers already in the library, with a named
  "recovery" regime added.
- Giancaterini, Hecq, Jasiak, Manafi Neyazi (2025) noncausal green-bubble
  paper (arXiv:2505.14911) — same mixed causal-noncausal mechanism as #7.

## Checked, not accessible (not verified from primary text, no claims made)

- Breitung & Kruse (2013), "When bubbles burst: econometric tests based on
  structural breaks," *Statistical Papers* 54(4) — Springer paywalled, no
  working-paper mirror found.
- "Testing for explosive bubbles in the presence of non-Gaussian
  conditions," *Economics Letters* 233 (2023) — ScienceDirect paywalled,
  author not even confirmed from search snippets.
- Montanino & De Luca, "The Bubble Crash GARCH model," SSRN 5604452 — SSRN
  gate page, no PDF retrievable.
- Horvath, Trapani & Wang (2024), "Sequential Monitoring for Explosive
  Volatility Regimes," arXiv:2404.17885 — downloaded, but its DGP wasn't
  clearly isolated in the time available and likely overlaps heavily with
  the already-collected Horvath-Trapani (2026) RCA-monitoring paper (same
  authors).
- Lin, Ren & Sornette (2009), LPPLS finite-time-singularity model,
  arXiv:0905.0128 — read; it's a curve-fitting paradigm, not a Monte Carlo
  DGP used to stress-test right-tailed unit-root tests. Out of scope, not
  "inaccessible."

## Candidacy, not just a survey

Unlike the other taxonomy files, none of these DGPs corresponds to a single
missing *statistic* — but that doesn't mean they're all equally far from
being worth shipping. Checked against `exuber/R/` directly (`grep`, no
`radf_ls`/level-shift code exists, confirming the point below):

**Cheap and useful now** — both are a parameter/mechanism swap inside
`sim_psy1`'s existing loop, not a new architecture, and both directly
support heteroskedasticity-robust tests exuber *already ships* (STADF, SBZ,
kernel-purge, sign-based, SSU — see
[volatility-robustness.md](/replication/volatility-robustness)) which currently have
no reusable public DGP to demo/validate against beyond each paper's own
one-off Monte Carlo script:

- **#10 GARCH(1,1)/TGARCH innovations** — replace fixed `sigma` with a
  GARCH(1,1) recursion (`h_t = ω + α·ε²_{t-1} + β·h_{t-1}`, optionally the
  TGARCH leverage term from #12). One conditional branch.
- **#9 non-Gaussian/heavy-tailed innovations** — a `dist =` argument
  swapping `rnorm()` for `rt()`/skew-t. No new mechanism at all.

**Cheap, but blocked on something else** — mechanically just as simple, but
pairs with a test that isn't implemented, so building the DGP first would
be dead weight:

- **#3 deterministic level shifts** — one jump term added to the mean
  equation, but Harvey, Leybourne, Tatlow & Zu (2025)'s level-shift test
  itself has no `radf_*` implementation yet. Worth building once that
  test lands (see [volatility-robustness.md](/replication/volatility-robustness)),
  not before.

**Real value, bigger lift** — worth doing, but not a quick addition:

- **#1/#2 CIR / lognormal-AR stochastic volatility** — proper SDE or
  near-integrated-log-variance simulation, not a parameter swap; would
  stress-test the same already-shipped tests beyond what GARCH covers.

**Genuinely a menu for later, not a queue** — #4-8, #11, #13-15 each need
either a fundamentally different simulation architecture (multi-series
factor model, bivariate linkage, latent state-space with a Markov chain,
branching tree, MAR lag polynomial) or pair with a multivariate/paradigm
test that also isn't implemented — these are correctly a survey rather
than near-term candidates.
