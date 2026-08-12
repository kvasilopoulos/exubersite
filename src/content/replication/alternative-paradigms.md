---
title: "Alternative paradigms"
blurb: "Non-ADF-family approaches, principally the quantile-based global test and its recursive monitoring extension."
order: 5
---
Methods that address the same problem (detecting explosive/bubble dynamics)
from outside the ADF/SADF/GSADF/BSADF recursive-regression family exuber is
built around.

| Method | Paper | Fit with exuber |
|---|---|---|
| [Quantile-based detection](#quantile-based-detection) | Pavlidis (2025); Wu, Shi & Wu (2025) | **global test AND `QPWY` monitoring done (2026-08-11)**; `QPSY` (double recursion) not implemented; Pavlidis's `Un`/`QKS` attempted then withdrawn — failed its own bootstrap-calibration validation, see below |
| [Noncausal / local explosive dynamics](#noncausal--local-explosive-dynamics) | Blasques, Koopman, Mingoli & Telg (2025) | evaluated, not implemented — different paradigm, no ADF machinery |
| [Spectral fragility](#spectral-fragility-out-of-scope) | Bhandari (arXiv) | out of scope |
| [Stochastic tree asset pricing](#stochastic-tree-asset-pricing-out-of-scope) | Gourieroux & Jasiak (2025) | out of scope (pricing, not testing) |

All from the *JTSA* 46(5) special issue except Bhandari (arXiv). Papers:
[references.md](/replication/references#alternative-paradigms).

## Quantile-based detection

**Status: Wu/Shi/Wu's "global test" done (2026-08-10); their `QPWY`
recursive monitoring extension also done (2026-08-11, `radf_qpwy()`),
re-triaged from "genuinely more expensive" — that verdict held for
`QPSY` (their double recursion, still not implemented) but not for
`QPWY` specifically, which was originally bundled with `QPSY` under the
same cost verdict without separating their very different profiles.
Pavlidis's quantile-autoregressive `Un`/`QKS` tests attempted
(2026-08-10) then withdrawn — implemented, found and fixed one real bug,
but still failed its own bootstrap-calibration validation against the
paper's own published Table 2, and was removed rather than shipped
broken (see "Deliberately not implemented" below for the full
diagnostic trail).**
Both PDFs read (introduction and statistic definitions; Wu/Shi/Wu's
formulas re-verified against rendered PDF pages 7-8 and 14, not just the
raw text extraction — the raw extraction had scrambled eq. 18 in a way
that would have produced a wrong statistic if implemented directly from
it, see "Implementation" below; Pavlidis's formulas re-verified against
rendered PDF pages 6-7).

### Sources

- Pavlidis, E.G. (2025). "Bubbles and crashes: A tale of quantiles." *JTSA*,
  46(5), 884-907.
- Wu, R., Shi, S. & Wu, J. (2025). "Quantile analysis for financial bubble
  detection and surveillance." *JTSA*, 46(5), 908-931.
  Shi is a PSY co-author, so this specific paper is a quantile-based test
  from inside the same lineage rather than a fully competing school.

### What it is

A genuinely distinct branch: instead of a recursive mean-regression ADF
statistic, these characterize explosive behavior via **quantile
regression (QR)** at chosen points of the conditional distribution.
Wu/Shi/Wu's global test statistic (their eq. 18, confirmed against a
rendered PDF page — the raw-text extraction garbles this equation into
`sqrt(f_hat(b_tau)/(1-tau))`, which is wrong) is the QR analogue of the
DF t-ratio:

```
t_T(tau) = [f_hat(b_tau) / sqrt(tau * (1 - tau))] * (Y_{-1}' P_Z Y_{-1})^{1/2} * (alpha_hat(tau) - 1)
```

where `alpha_hat(tau)` is the **quantile-regression** (not OLS) estimator
of `y_t` on an intercept and `y_{t-1}` at quantile `tau` (their eq. 13),
`Y_{-1}' P_Z Y_{-1}` is the demeaned sum of squares of the lagged level
(`P_Z` a pure demeaning projector, `Z` a column of ones), and
`f_hat(b_tau)` is a **kernel density estimate** of the first-differenced
series' density at its own `tau`-th sample quantile (their eq. 19,
`f_hat(b_tau) = (Th)^{-1} sum_t K((b_hat_tau - u_hat_t)/h)`, `u_hat_t =
y_t - y_{t-1}`) — needed to studentize the QR coefficient the way OLS's
residual-variance estimate studentizes the DF t-ratio. They also define
`QPWY`/`QPSY` monitoring statistics — the same recursive/double-recursive
sup-scan structure as PWY/PSY, but with this QR t-ratio computed at every
window instead of the OLS one (not implemented, see "Implementation"
below). Pavlidis's paper is the same underlying idea from a different
framing: unit-root quantile-autoregressive models, where the largest
autoregressive root is allowed to vary by quantile (below 1 at low
quantiles/crashes, above 1 at high quantiles/expansions) — not
implemented either.

**The critical value (their eq. 22-23), and a structural finding that
made it cheap**: the limiting null distribution of `t_T(tau)` is
`U(tau) = sqrt(1 - delta(tau)^2) * z + delta(tau) * Q`, with `z ~ N(0,1)`
and `delta(tau)` a correlation coefficient estimated directly from the
data (`cor(u_hat_t, tau - 1{u_hat_t < b_hat_tau})`, the correlation
between the innovation and its own quantile-check score). `Q` (their eq.
23, `(int W_bar^2)^{-1/2} int W_bar dW`) is **exactly the standard,
demeaned Dickey-Fuller t-statistic distribution** — verified empirically,
not just recognized from the formula: simulating it via the same
random-walk-plus-OLS-t-stat construction `radf_mc_cv()` already uses for
its own `adf` critical value gives values bit-for-bit identical to
`radf()`'s own single-shot `adf` field computed on the same simulated
series. So simulating `U(tau)`'s quantile needs no new statistical
machinery: one fresh standard normal draw combined with a critical value
this package already knows how to simulate.

The paper's own optimal-quantile-selection criterion (their eq. 33,
`tau* = argmin_{tau} tau(1-tau) / f_hat(b_tau)^2`) is a cheap grid search
over the same `f_hat(b_tau)` computation the test statistic already
needs — not the separate, harder problem the "Cost/feasibility" note
below originally worried it might be.

### Implementation

Shipped as `radf_quantile()` — the *global* test only (their Section
3.1: a single static QR fit at one quantile, no recursion), exactly the
"reasonable minimum-viable first cut" this note originally suggested.
Adds `quantreg` (>= 5.9) as exuber's **first genuinely new estimation
dependency** (`DESCRIPTION`'s `Imports`) — everything else in the
package stayed pure-R-plus-`exubercore`-C++ until this item, as
previously noted. `tau = "optimal"` (default) runs eq. 33's grid search
over `tau_grid` (default `seq(0.2, 0.8, by = 0.05)`, matching the
paper's own recommended practical range); a fixed `tau` can also be
passed directly. No C++ needed — `quantreg::rq()` handles the QR fit,
and the rest (density-at-quantile, the demeaned sum of squares, the
critical-value simulation) is plain R reusing existing patterns from
`radf_mc_cv()`.

**Validation**: empirical size under a pure random-walk null is 5.0%
(`tau = "optimal"`, nominal 5%, 100 reps) and 3.0% (fixed `tau = 0.5`) —
both essentially exact, unlike several other Monte-Carlo-validated items
in this project that ran conservative or inflated. Power under a genuine
explosive alternative is 100%, matching a standard SADF test run on the
same DGP as a rough cross-check. The optimal-`tau` selection lands
sensibly within the search grid across replications (never degenerates
to a boundary). Structural check: the critical value's `Q` component was
confirmed bit-for-bit identical to `radf()`'s own single-shot `adf`
field on the same simulated series, not merely assumed equal from the
formula. Replication script:
[replication/alternative-paradigms/radf_quantile_validation.R](#script-radf_quantile_validation).

### Implementation (QPWY) — done (2026-08-11)

**Re-triaged from "genuinely more expensive."** The recursive/double
-recursive *scanning structure* is familiar — same shape as `radf()`'s
own `sadf`/`gsadf` scan — but the *per-window computation* is not: every
one of PDC/KS, STADF, SBZ, kernel-purge, and the sign-based test reduces
its per-window estimate to a closed-form ratio of cumulative sums
(`O(1)` per window given prefix sums, `O(T)` total per scan); QR has no
such closed form (confirmed by re-reading rendered pages 10-11, their
Corollary 1-2, not assumed), so this part of the original assessment
holds for both `QPWY` and `QPSY`. What the original pass missed by
bundling the two together: `QPWY_r(tau) := t_T^{0,r}(tau)` is a
**single** recursion (window start fixed at `1`, only the end `r`
grows — exactly `radf()`'s own `badf` shape), needing `O(T)` actual QR
fits, the *same cost order* as `badf` itself — tractable, not
"genuinely more expensive" the way `QPSY`'s `O(T^2)` double recursion
(additionally optimizing over the window start too) genuinely is.

A second favorable finding: `QPWY`'s critical-value machinery reuses
`radf_quantile()`'s own already-validated construction, not new theory.
Their Corollary 1 decomposes the limiting distribution of the general
windowed statistic as `U'^{r1,r2}(tau) = sqrt(1-delta(tau)^2)*z +
delta(tau)*Q_{r1,r2}` — the *same* decomposition `radf_quantile()`'s own
critical value already uses — and their Corollary 2 identifies
`Q_{0,r}` (the case relevant to `QPWY`) with exactly `radf()`'s own
`badf` sequence under a simulated null path. One `radf()` call per
Monte Carlo replicate therefore gives the *whole* boundary-relevant path
at once; no new simulation theory, only evaluating the existing
construction along a path instead of at a single endpoint.

**Implementation**: shipped as `radf_qpwy(data, tau = 0.5, minw, nrep,
level, seed)` in `exuber/R/radf_qpwy.R`. `qpwy_stat_path()` is the `O(T)`
loop of actual `quantreg::rq()` fits (mirroring `radf_quantile()`'s own
per-window t-ratio construction exactly, eq. 18, just repeated over a
growing window); `qpwy_boundary_sim()` reuses `radf()` directly (one
call per null replicate) to get the whole `Q_{0,r}` path at once.

**A genuine bug found and fixed by Monte Carlo validation, not caught by
the formula alone**: an initial version used the **per-`r` marginal**
quantile of the simulated `U`-paths as an `r`-varying boundary — this
gave a false-alarm rate of `50%` against a nominal `5%` under `H0`, a
textbook case of the difference between a pointwise quantile and a
supremum/first-crossing-controlled boundary. Fixed by taking each
simulated path's own **maximum** first, then the quantile of those
maxima across replicates — exactly how `radf_mc_cv()`'s own `sadf_cv` is
constructed — giving a single flat boundary that correctly controls
`P(sup_r stat(r) > boundary)`, not the marginal probability at each `r`
separately. After the fix: false-alarm rate under `H0` (`n=150`, 60
reps) is `6.7%` against nominal `5%`; detection power on a genuine
post-training explosive DGP (60 reps) is `50.0%`, comparable to standard
`SADF`'s `55.0%` on the identical DGP — QPWY trades a little power for
its robustness to non-Gaussian innovations (the paper's own stated
motivation for the QR-based approach), a sensible tradeoff rather than a
defect. New `test-qpwy.R` (7 tests, including a check that the
supremum-calibrated boundary stochastically dominates any single-column
marginal quantile — the structural signature of the fix). Replication
script:
[replication/alternative-paradigms/radf_qpwy_validation.R](#script-radf_qpwy_validation).

**Still not implemented**: `QPSY` (the double recursion) — needs
`O(T^2)` genuine QR fits, a materially different computational cost
class, before any critical-value simulation multiplies it further.

**Pavlidis's quantile-autoregressive framing — attempted (2026-08-10),
withdrawn after failing its own validation, not shipped.** Re-triaged on
the theory that, like Wu/Shi/Wu above, "lower priority, different
parameterization" might be underselling it once actually read. It
wasn't: pages 6-7 (rendered PDF, eq. 5-11) give the same ADF regression
form `radf()` itself uses, fit by quantile regression at chosen quantiles
`tau` — `Un(tau) := n*(alpha1_hat(tau) - 1)` (coefficient-based) and
`QKS := sup_{tau in T} Un(tau)` (their eq. 11) — with critical values
from an explicit residual/sieve bootstrap (their page 7, steps 1-5)
structurally very close to `radf_sb_()`'s already-shipped Pedersen-
Schütte bootstrap (fit an `AR(q)` to `diff(y)` under `H0`, resample
centered residuals, recursively regenerate and cumulate). `quantreg` is
already an exuber dependency (added for `radf_quantile()`), so this
looked like a well-scoped, moderate-effort item with a rare bonus: the
paper's own Table 2 gives published Monte Carlo empirical size numbers
(N(0,1)/t3/t2 errors, `n = 100/200/300/400`) to validate directly
against — no other item in this bundle had that.

Implemented as `radf_qar()`, then validated against Table 2's `n=100`,
`N(0,1)` row rather than just trusted. It failed: `Un(tau=0.5)`'s
empirical size matched almost exactly (0.050 vs. published 0.053), but
`Un` at higher `tau` and `QKS` ran visibly oversized (0.075-0.100 vs.
0.052-0.063 published). Chasing this down found a real bug — the
bootstrap DGP was fitting the wrong AR order (`AR(lag+1)` instead of
`AR(lag)`, and an unwanted `AR(1)` at `lag=0` instead of Pavlidis's own
step 1, which reduces to a plain i.i.d. resample when `q=0`) — inherited
from adapting `radf_sb_()`'s specific convention without re-deriving from
Pavlidis's own formula. Fixed. But re-validating after the fix made
things *worse*, not better (0.10-0.20 oversized), which ruled out the
first fix as the whole story.

A decoupled diagnostic isolated the real issue: comparing the *oracle*
finite-sample null distribution of `Un`/`QKS` (1,000 genuine i.i.d.
random walks, no bootstrap at all) against the bootstrap's own implied
critical values (one series, `nboot` up to 1,999) showed the bootstrap
critical value staying substantially and persistently *below* the oracle
at `tau = 0.5/0.8/0.9` even at `nboot = 1999` — ruling out "just needs
more bootstrap replications" (that would shrink toward the oracle, not
plateau below it). Only `tau = 0.95` and `QKS` came close to the oracle
at large `nboot`. This is a genuine, structural calibration problem in
the bootstrap procedure at low/mid quantiles specifically, not sampling
noise and not the AR-order bug (already fixed by that point) — likely
related to `Un(tau) = n*(alpha1_hat(tau) - 1)`'s well-known extreme
sensitivity to small numerical differences in `alpha1_hat` (near 1 for a
unit root, so `n ≈ 99` amplifies third-decimal-place differences by
~100x), but the exact mechanism was not pinned down.

**Not shipped.** Per this project's standing rule against shipping code
that fails its own validation, `radf_qar()` and its tests were removed
rather than committed. Flagged precisely so a future pass starts from
"the AR-order fix was correct but insufficient; the remaining gap is at
low/mid `tau`, `QKS`/`tau=0.95` look calibrated even at large `nboot`" —
not from scratch. If revisited: either investigate the `alpha1_hat`
sensitivity directly (compare bootstrap vs. oracle `alpha1_hat`
*distributions*, not just the `n*(...)`-transformed statistic, at
several `tau`), or restrict to `QKS`-only (the paper's own recommended
headline statistic and the one that validated cleanly) rather than
exposing the individual `Un(tau)` statistics that didn't.

## Noncausal / local explosive dynamics

**Status: evaluated, not implemented (2026-08-09).** PDF read (abstract,
introduction).

### Source

Blasques, F., Koopman, S.J., Mingoli, G. & Telg, S. (2025). "A Novel Test
for the Presence of Local Explosive Dynamics." *JTSA*, 46(5), 966-980.
`doi:10.1111/jtsa.70001`.

### What it is

A different paradigm again, confirmed on reading: the test is built for
**mixed causal-noncausal autoregressive processes** — a model class where
part of the dynamics depends on *future* shocks (anticipative/noncausal
terms), not just past ones. The premise is that bubbles are generated by
an extreme shock acting through the forward-looking (noncausal) component
of the model, rather than emerging from a recursively-estimated explosive
AR root on the past alone. The test statistic's distribution is either
derived analytically or approximated numerically depending on the
assumed error distribution; the empirical application is a monthly oil
price index, framed partly as a Value-at-Risk-style risk-assessment tool
rather than purely a bubble detector.

### Cost/feasibility note for exuber

Confirmed as genuinely out of exuber's architecture, not merely stated as
such from the abstract: mixed causal-noncausal AR models require their
own specialized estimation (there is no closed-form OLS/QR reduction —
noncausal components are typically estimated via approximate/simulated
maximum likelihood under a specified non-Gaussian error distribution,
since noncausal processes are only identifiable when innovations are
non-Gaussian). None of `exubercore`'s recursive-least-squares machinery,
or any transform-then-reuse trick that worked for STADF/sign-based/PDC
above, applies here — this would be a from-scratch statistical framework
with a new estimation dependency (a noncausal-AR fitting routine; no
mainstream R package was identified for this during this pass). Not a
contained addition under any framing. Lowest priority of the three live
items in this file.

## Spectral fragility (out of scope)

Bhandari, A. "Rational Bubbles at the Spectral Edge: An Operator-Spectral
Theory of Fragility, Identification and Finite-Sample Certification."
arXiv:2607.03933.

Factor/co-movement spectral-fragility detection — identifies market
fragility via the strength of a dominant factor extracted from
cross-sectional co-movement (fewer independent factors during crises than
calm periods), not a right-tailed unit-root test on a single series at
all. Detects fragility contemporaneously, not predictively. Not built on
ADF/SADF machinery in any way; would need a from-scratch implementation
with no reuse of exuber's existing code, and arguably answers a different
question (market-wide fragility, not a specific series's explosiveness)
than what `radf()` and friends are for. Noted here so it doesn't get
re-discovered and re-evaluated later, not because it's a live candidate.

## Stochastic tree asset pricing (out of scope)

Gourieroux, C. & Jasiak, J. (2025). "A Stochastic Tree for Bubble Asset
Modelling and Pricing." *JTSA*, 46(5), 932-944.

An asset-*pricing* model (a stochastic tree representation for modelling,
forecasting, and pricing bubbles, with closed-form option-pricing
formulas), not a test for the presence of a bubble. exuber's scope is
testing (`radf()` and friends test for explosiveness); this doesn't fit
regardless of how it's implemented. Noted for awareness, not a candidate.
