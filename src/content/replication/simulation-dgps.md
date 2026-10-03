---
title: "Simulation DGPs"
blurb: "Data-generating processes for the axes the original sim_*() functions do not cover."
order: 6
---
This file differs from the other family files, because it catalogues
data-generating processes (DGPs) for simulating bubble series and not tests or
statistics. The DGPs come from the Monte Carlo sections of the literature the
project has already read, plus a few papers pulled in specifically for this
pass. The question the file answers is what the `sim.R` of exuber (`sim_psy1`,
`sim_psy2`, `sim_ps1`, `sim_ps2`, `sim_blan`, `sim_evans`, `sim_div`) does not
cover.

All seven existing `sim_*` functions share one structural property, which is
fixed (homoskedastic), i.i.d. Gaussian innovations. They differ in the mean
equation (PSY-style regime-switching AR(1) against Blanchard/Evans rational
bubbles) and in the number and shape of bubbles and collapses. None has
time-varying volatility, GARCH, non-Gaussian innovations, stochastic switch
timing or a multi-series or factor structure. We flag an item below only if it
differs on one of those axes through a genuinely different mechanism, and not
through a reparameterization of an equation that is already implemented (see
the exclusions at the end).

The status legend is the one in [volatility-robustness.md](/replication/volatility-robustness).

Status: all 15 catalogued DGPs were implemented on 2026-08-11 in
`exuber/R/sim.R` and tested in `exuber/tests/testthat/test-sim-dgp.R` (39
assertions). The tests use formula-exact checks against independent
brute-force reimplementations where feasible, checks of published properties
(for example the price floor of `sim_tree()` and the empirical transition-rate
match of `sim_msbubble()`), and moment and reproducibility checks otherwise.
Two groups of DGPs compose cleanly as optional arguments on the existing
`sim_psy1()` (`e`, `shifts`, `coef_noise`), plus five small, independently
reusable generator functions (`sim_innov()`, `sim_vol_garch()`,
`sim_vol_cir()`, `sim_vol_sv()`, `sim_fi()`). The first group is #1/#2/#3, and
the second is #8/#9/#10/#11/#12 (stochastic volatility, innovation
distribution, level shift, long memory and stochastic coefficient). This avoids
thirteen near-duplicate clones of the loop in `sim_psy1`. `sim_blan()` gained a
`type = "rotermann_wilfling"` option (#13). The remaining six DGPs with
distinct architectures (#4-7, #14, #15) shipped as standalone functions:
`sim_common()`, `sim_coexplosive()`, `sim_tree()`, `sim_mar()`,
`sim_msbubble()` and `sim_falsebubble()`.

This is DGP machinery only, as in the original survey. Building a DGP does not
mean that its paired test is implemented. As of 2026-08-11, the level-shift
robustness result of HLTZ (2025), the motivation for DGP #3 in their §3, is
implemented. See
[volatility-robustness.md](/replication/volatility-robustness#level-shift-robustness-hltz-2025)
for the level-shift robustness of `radf_sign()` and the new `radf_sign_dm()`.
`sim_falsebubble()` and `sim_msbubble()` still have no dedicated exuber test.
They are stress-test and demo series for the existing PSY/GSADF machinery, the
same role that `sim_evans()` already plays. Each item below says what, if
anything, it is meant to validate against.

Validation found two real bugs, and neither was a design choice. The first
`sim_fi()` implementation used `stats::filter(..., sides = 1)` for the
truncated MA(∞) convolution, with an innovation vector of the same length as
the filter itself. Filtering with `sides = 1` needs strictly more input than
filter taps to produce any non-`NA` output, so the function returned `NA` for
the entire series except one point. We fixed it by generating `m` extra
(truncation-lag) innovations beyond the requested length and convolving by hand
(`vapply()`), and no longer use `stats::filter()` or `stats::convolve()`. Both
of those also turned out to segfault unconditionally on this project's R
installation (Windows, R 4.6.1), independently of exuber (it reproduces on
`stats::filter(rnorm(120), rep(0.01, 21), method = "convolution")` with no
package loaded). The hand-rolled convolution is therefore a portability fix as
well as a correctness fix. The second bug was in `sim_tree()`. Its
`p_t = Φ(X_t)` can round to exactly 0 or 1 in the tail of a long series
(floating-point underflow). That silently turned `ξ_1t` and `ε_t` into `0/0`
(`NaN`), which then propagated through every subsequent `y[t]`. We fixed it by
clipping `p_t` into `[1e-10, 1 - 1e-10]`.

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
use this as a robustness design beyond their main deterministic-volatility one:

```
dσ²(r) = 0.03(0.25 − σ²(r))dr + 0.1σ(r)dB(r)
```

It is a square-root (Cox-Ingersoll-Ross) diffusion for volatility
"representative of Bollerslev and Zhou (2002)". Each replication simulates it
from NIID(0,1) Brownian increments, independent of the level innovations, and
it feeds the same PSY-style level process that `sim_psy1` already implements.
The new element relative to `sim.R` is continuous-time stochastic volatility,
where `sim.R` has fixed or deterministic volatility.

### 2. AR(1) lognormal stochastic volatility

Sarkar & Wells (2025, eq. in §2)
reused in Sarkar & Wells (2026, eq. 7):

```
y_t = ρ_n y_{t-1} + u_t,   u_t = σ_t ε_t,   log σ_t² = φ_n log σ_{t-1}² + η_t
```

Both `ρ_n → 1` and `φ_n → 1`, which makes the process "double local-to-unity"
in the mean and in the variance. It is a discrete-time stochastic volatility
model with a persistent, near-integrated log-variance, applied to a mildly
explosive alternative. No existing `sim_*` function has a persistent variance
process, as opposed to an i.i.d. or deterministic one.

### 3. Deterministic level shifts

Harvey, Leybourne, Tatlow & Zu (2025):

```
y_t = y_{t-1} + Σ_{j=1}^{m} δ_j·1(t≥τ_j) + ε_t     (null)
```

The process has an arbitrary number `m` of jump discontinuities of magnitude
`δ_j` at unknown dates `τ_j`, with an explosive regime layered on top for the
alternative (§5.2). The new element is structural breaks in the mean level.
Nothing in `sim.R` has jump components at all, and it has only smooth
regime-switching.

### 4. Latent common-factor bubble

Chen, Phillips & Shi (2023)
eq. 2.3/2.7-2.9:

```
X_t = Λf_t + e_t
```

There are `N` observed series and one latent factor `f_t` that follows a
PSY-style unit-root, explosive, collapse sequence. The loadings are
`Λ ~ U[0,2]`, and the idiosyncratic noise has `σ_e = 0.1`. A collapse-splicing
construction (paper lines 960-977) avoids discontinuities at the regime
boundary. The new element is a factor-model DGP in which one common bubble
drives many series jointly. Every existing `sim_*` function is single-series.

### 5. Bivariate co-explosive linkage

Evripidou, Harvey, Leybourne & Sollis (2022):

`x_t` is generated from the regime-dummy Models 1-4 of HLS (each individually
PSY-style), and a second series is linked to it:

```
y_t = μ_y + φ_x·x_{t-i} + φ_z·z_t + ε_{y,t}
```

The second series combines a lead or lagged copy of the explosive series `x_t`
with a third, latent explosive series `z_t`, and the heteroskedasticity is
timed to the regime changes. The new element is a two-series lead/lag linkage
of explosive components, in place of a single univariate path. It differs from
the common-factor case above (#4), which shares one factor and does not link
two distinct explosive series.

### 6. Stochastic branching-tree bubble

Gourieroux & Jasiak (2025)
describe an affine autoregression with a stochastic coefficient. It can be
represented as a random-coefficient AR process generated by a binomial tree
with stochastic branching intensity, in contrast to the deterministic branches
of Cox-Ross-Rubinstein. The Blanchard & Watson (1982) bubble, already
`sim_blan`, is the special case of constant intensity. The new element is that
the branching mechanism generates the path directly, with no fixed collapse
probability.

### 7. Mixed causal-noncausal AR (MAR)

Blasques, Koopman, Mingoli & Telg (2025):

```
(1 − φ1 L)(1 − ψ1 L⁻¹) y_t = ε_t
```

The causal root is `φ1 = 0.7`, `ψ1` is the noncausal root, and the innovations
are Cauchy or Student-t(2). The noncausal component autonomously generates
transient, self-terminating local bubbles, with no scripted regime dates at
all. The process differs on every axis. It has a lag-polynomial form in place of
a regime-dummy AR, an implicit and unscripted collapse, and heavy-tailed
non-Gaussian innovations.

### 8. Fractionally-integrated innovations

Lui, Phillips & Yu (2024):

```
y_t = y_{t-1} + u_t,   u_t = Δ⁻ᵈ ε_t   (d > 0, ε_t iid with finite (2+δ) moments)
```

The innovations themselves are `FI(d)` (long-memory) and not i.i.d. The paper
develops an explosive-alternative analogue in §4. The new element is
non-i.i.d., long-range-dependent noise feeding the unit-root or explosive
equation.

### 9. PSY equation with heavy-tailed/skewed innovations

Wu, Shi & Wu (2025)
eq. 6, use structurally the standard PSY unit-root, explosive, collapse
equation, but draw the innovations from `N(0,1)`, `t(3)`, skewed-`t(3,−0.75)`
and skewed-`t(3,+0.75)`. We flag it only for the innovation distribution. The
mean equation itself is not new, but every existing `sim_*` function uses fixed
Gaussian noise.

### 10. GARCH(1,1) / smooth-transition volatility

Whitehouse, Harvey & Leybourne (2025):

```
z_t = h_t^{1/2}·ε_t,   h_t = 0.1 + 0.1z²_{t-1} + 0.8h_{t-1}
```

Harvey, Leybourne, Taylor & Zu (2024) reuse the same GARCH(1,1) specification
Alongside their main design, they use a logistic smooth-transition volatility
function:

```
σ(r) = σ1 + (σ2 − σ1)/(1 + exp{−κ(r−δ)})
```

The new elements are proper conditional GARCH heteroskedasticity, and a smooth
transition between volatility regimes in place of an instant jump.

### 11. Stochastic explosive coefficient

Kurozumi & Nishi (2025)
is already documented in [volatility-robustness.md](/replication/volatility-robustness#implementation-ssu-done-2026-08-10).
The coefficient `1 + c1/T + a·u_t/√T` replaces the deterministic `1+c/Tᵅ`, so
the persistence parameter itself is random and not only the noise scale. We
list it here for completeness because it meets the bar of this file. It adds no
new information.

## New papers pulled in for this pass (not previously in the paper library)

We searched beyond the 48 papers already collected, specifically for DGPs. We
downloaded and text-extracted the new papers the same way as the rest of the
library of this project (`pdftotext -layout`, with a PNG-render fallback for
garbled formulas).

### 12. TGARCH(1,1) with leverage effect

Monschang, V. & Wilfling, B. (2021). "Sup-ADF-style bubble-detection
methods under test." *Empirical Economics*, 61, 145-172.
`doi:10.1007/s00181-020-01859-7`. Open access (CQE Working Paper 78/2019):

```
ε_t = s_t·h_t^{1/2},   h_t = ω + α·ε²_{t-1} + β·h_{t-1} + γ·ε²_{t-1}·1(ε_{t-1}<0)
```

The process is calibrated to NASDAQ estimates (`α=.4387, γ_indicator=.1306,
β=.9319`). Relative to #10 above, the new element is an asymmetric
(sign-dependent) shock response, the leverage effect, in place of a symmetric
or smooth-transition variance.

It is implemented as `sim_vol_garch(omega, alpha, beta, gamma)`, the same
function as #10, with `gamma` as the shared leverage parameter.
`sim_vol_garch(omega = 0.4387, alpha = 0, beta = 0.9319, gamma = 0.1306)`
reproduces the NASDAQ calibration of this paper exactly, and the default
`gamma = 0` gives plain GARCH(1,1), which is item #10.

### 13. Lognormal-mixture rational bubble (Rotermann-Wilfling)

The same paper, eq. 4:

```
B_{t+1} = (B_t·u_t/δ)                    w.p. π
        = ((1−πδ)/(1−π))·B_t·u_t          w.p. 1−π,   u_t ~ iid lognormal
```

It produces periodically recurring, stochastically deflating trajectories and
no single full collapse to a fixed floor, as `sim_blan` and `sim_evans` do. The
new element is a collapse mechanism that differs from both existing
rational-bubble DGPs, with partial, probabilistic deflation in place of total
collapse to noise.

It is implemented as `sim_blan(type = "rotermann_wilfling", delta,
rw_sigma)`, a new branch of the existing function and not a standalone one,
because it shares the two-regime-with-probability-`pi` structure of `sim_blan`
and only the update rule per regime differs. `test-sim-dgp.R` verifies that it
stays strictly positive, which is a structural invariant of the multiplicative
recursion.

### 14. Markov-switching present-value bubble

Chan, J.C.C. & Santi, C. (2021). "Speculative Bubbles in Present-Value
Models: A Bayesian Markov-Switching State Space Approach." *J. Economic
Dynamics and Control*, 127, 104101. Open access (author's site):

```
b_t = (1/λ_{S_t+1})·b_{t-1} + ε_bt,   S_t ∈ {1,2} first-order Markov, transition probs p11/p22
```

Regime 1 is "surviving" (`λ<1`, explosive) and regime 2 is "collapsing"
(`λ>1`, mean-reverting). The bubble is embedded in a full present-value
state-space model with time-varying expected returns and dividend growth, and
the section "Simulated Datasets" of the paper (Table 2) generates artificial
data from these parameters. The new element is that the timing of the switch is
itself stochastic (a Markov chain) and not a deterministic break fraction.
Every PSY-style DGP in `sim.R` and above scripts `te` and `tf` as fixed
fractions of `n`.

It is implemented as `sim_msbubble(p11, p22, lambda1, lambda2, sigma_b)`. The
function covers only the bubble component `b_t` and not the full
present-value state-space model, because the expected-returns and dividend-growth
machinery of the paper is orthogonal to what a stress-test DGP needs. It has one
indexing simplification. Eq. 16 of the source applies the regime coefficient via
`S_{t+1}` (the realized regime of the next period), whereas this implementation
uses the contemporaneous `S_t`. That does not change the qualitative
Markov-switching mechanism. It changes only which time step a given draw of `S`
labels. `test-sim-dgp.R` verifies that the empirical self-transition rate of
the simulated regime path matches `p11` and `p22` to within Monte Carlo
tolerance.

### 15. Deterministic technology-adoption "false bubble" null

Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). "Technology
Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI
Episodes." arXiv:2604.25826. Open

The paper embeds a hump-shaped (triangular, Gaussian, Beta or Gamma)
deterministic technology-adoption shock into the Campbell-Shiller present-value
fundamental. This makes the fundamental price locally explosive during adoption
with no bubble present at all, and the design aims to make PSY-style tests
reject spuriously. Appendix E.5 extends it to a Bayesian-updated stochastic
version. The new element is a null (no-bubble) DGP with a smooth deterministic
drift, and not a bubble-present alternative. None of the functions in `sim.R`
produces a series of fundamentals only, with no bubble, that still looks
locally explosive.

It is implemented as `sim_falsebubble(t1, t2, kappa, shape, amplitude, mu, r)`.
`shape = "triangular"` reproduces the worked example of eq. 4 in the paper
exactly. `shape = "gaussian"` is one alternative from the paper's own list of
"any hump-shaped specification". We did not add Beta or Gamma, because by the
paper's own robustness claim they make no qualitative difference from
triangular or Gaussian for this purpose. The technology shock is deterministic,
so its price contribution is an exact forward-looking discounted sum
(`T_t = Σ_{s>t} β^{s-t}τ_s`) and no further simulation approximation. This is a
simplified, single-shock reproduction of the mechanism (a deterministic hump
gives a hump-shaped fundamental price, with no bubble), and it leaves out the
full DOLS and multi-functional-form robustness machinery of the paper.
`test-sim-dgp.R` verifies that `amplitude = 0` reduces to the fundamental-price
formula of `sim_div()` bit for bit (a clean special-case check), and that the
technology term is exactly zero outside `[t1, t2]`.

### Adjacent, not a price-level DGP

Richter, S., Wang, W. & Wu, W.B. (2023, orig. 2018). "A supreme test for
periodic explosive GARCH." *Econometrics* (MDPI); arXiv:1812.03475. Open
The paper studies a piecewise or periodic explosive GARCH(1,1) in which the
explosiveness lives in the volatility recursion (`α` and `β` are temporarily
driven toward or through the IGARCH boundary `αΣ+βΣ≥1`) and not in the price
level. We keep it in the library as a volatility-bubble reference. We do not
count it above, because it targets a different detection problem (volatility
bubbles, not price bubbles) than anything that `radf()` or `sim.R` addresses.

## Excluded as reparameterizations, not new DGPs

- The 6 Monte Carlo DGPs "A-F" of Harvey, Leybourne & Whitehouse (2020), with
  2- and 3-bubble sequences. They belong to the same regime-dummy AR(1) family
  as HLS/PSY, and only the regime count and the parameters differ from
  `sim_psy2`.
- The confidence-sets paper of Kurozumi & Skrobotov (2025/2026). It has the
  same linear-AR(1)-with-switching-coefficient form as the other
  Kurozumi/Skrobotov papers already in the library, with a named "recovery"
  regime added.
- The noncausal green-bubble paper of Giancaterini, Hecq, Jasiak and Manafi
  Neyazi (2025) (arXiv:2505.14911), which uses the same mixed causal-noncausal
  mechanism as #7.

## Checked, not accessible (not verified from primary text, no claims made)

- Breitung & Kruse (2013), "When bubbles burst: econometric tests based on
  structural breaks," *Statistical Papers* 54(4). Springer paywalls it, and we
  found no working-paper mirror.
- "Testing for explosive bubbles in the presence of non-Gaussian
  conditions," *Economics Letters* 233 (2023). ScienceDirect paywalls it, and
  the search snippets did not even confirm the author.
- Montanino & De Luca, "The Bubble Crash GARCH model," SSRN 5604452. The SSRN
  gate page gave no retrievable PDF.
- Horvath, Trapani & Wang (2024), "Sequential Monitoring for Explosive
  Volatility Regimes," arXiv:2404.17885. We downloaded it, but we did not
  isolate its DGP in the time available, and it likely overlaps heavily with
  the already-collected Horvath-Trapani (2026) RCA-monitoring paper by the same
  authors.
- Lin, Ren & Sornette (2009), LPPLS finite-time-singularity model,
  arXiv:0905.0128. We read it. It is a curve-fitting paradigm and not a Monte
  Carlo DGP used to stress-test right-tailed unit-root tests. It is out of
  scope, and "inaccessible" does not describe it.

## Candidacy, not just a survey

Unlike the other family files, none of these DGPs corresponds to a single
missing statistic, but they are not all equally far from being worth shipping.
We checked against `exuber/R/` directly with `grep` and found no
`radf_ls` or level-shift code, which confirms the point below. The assessment
that follows dates from the survey, before we built the DGPs.

Cheap and useful at the time of the survey. Both were a parameter or mechanism
swap inside the existing loop of `sim_psy1` and no new architecture. Both
directly support heteroskedasticity-robust tests that exuber already ships
(STADF, SBZ, kernel-purge, sign-based and SSU, see
[volatility-robustness.md](/replication/volatility-robustness)), which had no reusable
public DGP to demo or validate against beyond the one-off Monte Carlo script of
each paper.

- #10 GARCH(1,1)/TGARCH innovations replace the fixed `sigma` with a
  GARCH(1,1) recursion (`h_t = ω + α·ε²_{t-1} + β·h_{t-1}`, optionally with the
  TGARCH leverage term from #12). This takes one conditional branch.
- #9 non-Gaussian and heavy-tailed innovations need a `dist =` argument that
  swaps `rnorm()` for `rt()` or a skew-t. No new mechanism is involved at all.

Cheap but blocked on something else. These are mechanically just as simple, but
they pair with a test that was not implemented at the time, so building the DGP
first would have been dead weight.

- #3 deterministic level shifts add one jump term to the mean equation, but the
  level-shift test of Harvey, Leybourne, Tatlow & Zu (2025) had no `radf_*`
  implementation yet. It was worth building once that test landed (see
  [volatility-robustness.md](/replication/volatility-robustness)) and not before.

Real value, bigger lift. These are worth doing but are no quick addition.

- #1/#2 CIR and lognormal-AR stochastic volatility need a proper SDE or
  near-integrated-log-variance simulation, which is more than a parameter swap.
  They would stress-test the same already-shipped tests beyond what GARCH
  covers.

A menu for later and not a queue. #4-8, #11 and #13-15 each need either a
fundamentally different simulation architecture (a multi-series factor model, a
bivariate linkage, a latent state space with a Markov chain, a branching tree
or a MAR lag polynomial) or a multivariate or paradigm-specific test that was
also not implemented. They are correctly a survey and not near-term candidates.
