---
title: "Alternative paradigms"
blurb: "Non-ADF-family approaches, principally the quantile-based global test and its recursive monitoring extension."
order: 5
---
This file covers methods that address the same problem, detecting
explosive or bubble dynamics, from outside the ADF/SADF/GSADF/BSADF
recursive-regression family around which exuber is built.

| Method | Paper | Fit with exuber |
|---|---|---|
| [Quantile-based detection](#quantile-based-detection) | Pavlidis (2025); Wu, Shi & Wu (2025) | the global test and `QPWY` monitoring are done (2026-08-11); `QPSY` (double recursion) is not implemented; Pavlidis's `Un`/`QKS` was attempted and withdrawn after it failed its own bootstrap-calibration validation, see below |
| [Noncausal / local explosive dynamics](#noncausal--local-explosive-dynamics) | Blasques, Koopman, Mingoli & Telg (2025) | evaluated, not implemented; a different paradigm with no ADF machinery |
| [Spectral fragility](#spectral-fragility-out-of-scope) | Bhandari (arXiv) | out of scope |
| [Stochastic tree asset pricing](#stochastic-tree-asset-pricing-out-of-scope) | Gourieroux & Jasiak (2025) | out of scope (pricing, not testing) |

All papers are from the *JTSA* 46(5) special issue except Bhandari (arXiv).
Papers: [references.md](/replication/references). PDFs:
the paper library.

## Quantile-based detection

Status: the "global test" of Wu, Shi & Wu was done on 2026-08-10, and their
`QPWY` recursive monitoring extension followed on 2026-08-11 as
`monitor_quantile()`. We had first rated both as genuinely more expensive. That
verdict holds for `QPSY` (their double recursion, still not implemented) but not
for `QPWY`, which the original assessment had bundled with `QPSY` under one cost
verdict without separating their very different profiles. We attempted
Pavlidis's quantile-autoregressive `Un`/`QKS` tests on 2026-08-10 and then
withdrew them. The implementation had one real bug, which we found and fixed,
but it still failed its own bootstrap-calibration validation against the
published Table 2 of the paper, so we removed it and did not ship it broken. The subsection on Pavlidis's quantile-autoregressive framing below has the full diagnostic trail.

We read both PDFs (the introduction and the statistic definitions). We
re-verified the formulas of Wu/Shi/Wu against rendered PDF pages 7-8 and 14,
and not only against the raw text extraction, because the raw extraction had
scrambled eq. 18 in a way that would have produced a wrong statistic if
implemented directly from it (see "Implementation" below). We re-verified the
formulas of Pavlidis against rendered PDF pages 6-7.

### Sources

- Pavlidis, E.G. (2025). "Bubbles and crashes: A tale of quantiles." *JTSA*,
  46(5), 884-907.
- Wu, R., Shi, S. & Wu, J. (2025). "Quantile analysis for financial bubble
  detection and surveillance." *JTSA*, 46(5), 908-931.
  Shi is a PSY co-author, so this paper is a quantile-based test from inside
  the same lineage and not a fully competing school.

### What it is

These papers form a distinct branch. They characterize explosive behavior
through quantile regression (QR) at chosen points of the conditional
distribution, in place of a recursive mean-regression ADF statistic. The
global test statistic of Wu/Shi/Wu (their eq. 18, confirmed against a rendered
PDF page, because the raw-text extraction garbles this equation into the wrong
`sqrt(f_hat(b_tau)/(1-tau))`) is the QR analogue of the DF t-ratio:

```
t_T(tau) = [f_hat(b_tau) / sqrt(tau * (1 - tau))] * (Y_{-1}' P_Z Y_{-1})^{1/2} * (alpha_hat(tau) - 1)
```

Here `alpha_hat(tau)` is the quantile-regression (not OLS) estimator of `y_t`
on an intercept and `y_{t-1}` at quantile `tau` (their eq. 13).
`Y_{-1}' P_Z Y_{-1}` is the demeaned sum of squares of the lagged level (`P_Z`
is a pure demeaning projector and `Z` a column of ones). `f_hat(b_tau)` is a
kernel density estimate of the density of the first-differenced series at its
own `tau`-th sample quantile (their eq. 19,
`f_hat(b_tau) = (Th)^{-1} sum_t K((b_hat_tau - u_hat_t)/h)`, with `u_hat_t =
y_t - y_{t-1}`). We need it to studentize the QR coefficient, in the way that
the residual-variance estimate of OLS studentizes the DF t-ratio. The paper
also defines `QPWY` and `QPSY` monitoring statistics. They have the same
recursive and double-recursive sup-scan structure as PWY and PSY, but they
compute this QR t-ratio at every window instead of the OLS one (not
implemented, see "Implementation" below). Pavlidis's paper has the same
underlying idea in a different framing, namely unit-root quantile-autoregressive
models in which the largest autoregressive root may vary by quantile (below 1
at low quantiles and crashes, above 1 at high quantiles and expansions). We have
not implemented that either.

The critical value (their eq. 22-23) rests on a structural finding that made it
cheap. The limiting null distribution of `t_T(tau)` is
`U(tau) = sqrt(1 - delta(tau)^2) * z + delta(tau) * Q`, with `z ~ N(0,1)` and
`delta(tau)` a correlation coefficient estimated directly from the data
(`cor(u_hat_t, tau - 1{u_hat_t < b_hat_tau})`, the correlation between the
innovation and its own quantile-check score). `Q` (their eq. 23,
`(int W_bar^2)^{-1/2} int W_bar dW`) is exactly the standard demeaned
Dickey-Fuller t-statistic distribution. We verified this empirically and did not
merely recognize it from the formula. Simulating `Q` with the same
random-walk-plus-OLS-t-stat construction that `radf_mc_cv()` already uses for
its own `adf` critical value gives values bit for bit identical to the
single-shot `adf` field of `radf()` on the same simulated series. Simulating the
quantile of `U(tau)` therefore needs no new statistical machinery. It takes one
fresh standard normal draw combined with a critical value that this package
already knows how to simulate.

The criterion of the paper for choosing the optimal quantile (their eq. 33,
`tau* = argmin_{tau} tau(1-tau) / f_hat(b_tau)^2`) is a cheap grid search over
the same `f_hat(b_tau)` computation that the test statistic already needs. The
"Cost/feasibility" note below had originally worried that it might be a
separate, harder problem.

### Implementation

The function ships as `quantile_test()`. It implements the global test only
(their Section 3.1, a single static QR fit at one quantile with no recursion),
which is the "reasonable minimum-viable first cut" that this note originally
suggested. It adds `quantreg` (>= 5.9) as the first new estimation dependency of
exuber (`DESCRIPTION`'s `Imports`). Until this item, everything else in the
package stayed pure R plus the `exubercore` C++ code, as noted earlier.
`tau = "optimal"` (the default) runs the grid search of eq. 33 over `tau_grid`
(default `seq(0.2, 0.8, by = 0.05)`, which matches the practical range that the
paper recommends). A fixed `tau` can also be passed directly. No C++ is needed.
`quantreg::rq()` handles the QR fit, and the rest (the density at the quantile,
the demeaned sum of squares and the critical-value simulation) is plain R that
reuses existing patterns from `radf_mc_cv()`.

Validation shows an empirical size under a pure random-walk null of 5.0%
(`tau = "optimal"`, nominal 5%, 100 reps) and 3.0% (fixed `tau = 0.5`). Both are
essentially exact, unlike several other Monte-Carlo-validated items in this
project that ran conservative or inflated. Power under an explosive alternative
is 100%, which matches a standard SADF test on the same DGP as a rough
cross-check. The optimal-`tau` selection lands inside the search grid across
replications and never degenerates to a boundary. As a structural check, we
confirmed that the `Q` component of the critical value is bit for bit identical
to the single-shot `adf` field of `radf()` on the same simulated series, and did
not merely assume it equal from the formula. Replication script:
[replication/alternative-paradigms/radf_quantile_validation.R](#script-radf_quantile_validation).

### Implementation (QPWY), done (2026-08-11)

We had first rated this item as genuinely more expensive, and re-triage changed
that. The recursive and double-recursive scanning structure is familiar, with
the same shape as the `sadf`/`gsadf` scan of `radf()`. The per-window
computation is not familiar. PDC/KS, STADF, SBZ, kernel-purge and the
sign-based test all reduce the per-window estimate to a closed-form ratio of
cumulative sums (`O(1)` per window given prefix sums, `O(T)` total per scan).
QR has no such closed form. We confirmed that by re-reading rendered pages
10-11 (their Corollary 1-2) and did not assume it, so this part of the original
assessment holds for both `QPWY` and `QPSY`. The original pass missed a
difference by bundling the two together. `QPWY_r(tau) := t_T^{0,r}(tau)` is a
single recursion. The window start is fixed at `1` and only the end `r` grows,
which is exactly the `badf` shape of `radf()`. It needs `O(T)` actual QR fits,
the same cost order as `badf` itself. It is tractable, and not "genuinely more
expensive" in the way that the `O(T^2)` double recursion of `QPSY` is, which
also optimizes over the window start.

A second favorable finding is that the critical-value machinery of `QPWY`
reuses the construction in `quantile_test()` that we had already validated, and
needs no new theory. Corollary 1 of the paper decomposes the limiting
distribution of the general windowed statistic as
`U'^{r1,r2}(tau) = sqrt(1-delta(tau)^2)*z + delta(tau)*Q_{r1,r2}`. This is the
same decomposition that the critical value of `quantile_test()` already uses.
Their Corollary 2 identifies `Q_{0,r}` (the case relevant to `QPWY`) with exactly
the `badf` sequence of `radf()` under a simulated null path. One `radf()` call
per Monte Carlo replicate therefore gives the whole boundary-relevant path at
once. We need no new simulation theory, only the existing construction
evaluated along a path and not at a single endpoint.

The function ships as `monitor_quantile(data, tau = 0.5, minw, nrep, level,
seed)` in `exuber/R/monitor_quantile.R`. `qpwy_stat_path()` is the `O(T)` loop of
actual `quantreg::rq()` fits. It mirrors the per-window t-ratio construction of
`quantile_test()` exactly (eq. 18), repeated over a growing window.
`qpwy_boundary_sim()` reuses `radf()` directly, with one call per null
replicate, to get the whole `Q_{0,r}` path at once.

Monte Carlo validation found a real bug that the formula alone did not reveal.
An initial version used the per-`r` marginal quantile of the simulated
`U`-paths as an `r`-varying boundary. It gave a false-alarm rate of `50%`
against a nominal `5%` under `H0`, a textbook case of the difference between a
pointwise quantile and a boundary that controls the supremum or first crossing.
We fixed it by taking each simulated path's own maximum first and then the
quantile of those maxima across replicates, which is how the `sadf_cv` of
`radf_mc_cv()` is constructed. The result is a single flat boundary that
controls `P(sup_r stat(r) > boundary)` correctly, and not the marginal
probability at each `r` separately. After the fix, the false-alarm rate under
`H0` (`n=150`, 60 reps) is `6.7%` against a nominal `5%`. Detection power on a
post-training explosive DGP (60 reps) is `50.0%`, comparable to the `55.0%` of
standard `SADF` on the identical DGP. `QPWY` trades a little power for
robustness to non-Gaussian innovations, which is the paper's own stated
motivation for the QR-based approach. That is a sensible tradeoff and no
defect. `test-qpwy.R` is new (7 tests, including a check that the
supremum-calibrated boundary stochastically dominates any single-column marginal
quantile, which is the structural signature of the fix). Replication script:
[replication/alternative-paradigms/radf_qpwy_validation.R](#script-radf_qpwy_validation).

`QPSY` (the double recursion) is still not implemented. It needs `O(T^2)` actual
QR fits, a materially different computational cost class, before any
critical-value simulation multiplies the cost further.

#### Pavlidis's quantile-autoregressive framing: attempted and withdrawn (2026-08-10)

We attempted it and withdrew it after it failed its own validation, so it is not
shipped. We re-triaged it on the theory that, like Wu/Shi/Wu above, "lower
priority, different parameterization" might undersell the method once we read it
closely. It did not. Pages 6-7 (rendered PDF, eq. 5-11) give the same ADF
regression form that `radf()` itself uses, fitted by quantile regression at
chosen quantiles `tau`. The statistics are `Un(tau) := n*(alpha1_hat(tau) - 1)`
(coefficient-based) and `QKS := sup_{tau in T} Un(tau)` (their eq. 11). The
critical values come from an explicit residual or sieve bootstrap (their page
7, steps 1-5), which is structurally very close to the Pedersen-Schütte
bootstrap of `radf_sb_()` that we had already shipped. It fits an `AR(q)` to
`diff(y)` under `H0`, resamples centered residuals, and regenerates and
cumulates the series recursively. `quantreg` is already an exuber dependency
(added for `quantile_test()`), so this looked like a well-scoped,
moderate-effort item with a rare bonus. Table 2 of the paper gives published
Monte Carlo empirical sizes (N(0,1), t3 and t2 errors, `n = 100/200/300/400`) to
validate against directly, which no other item in this bundle offered.

We implemented it as `radf_qar()` and validated it against the `n=100`, `N(0,1)`
row of Table 2 and did not simply trust it. It failed. The empirical size of
`Un(tau=0.5)` matched almost exactly (0.050 against 0.053 published), but `Un`
at higher `tau` and `QKS` ran visibly oversized (0.075-0.100 against 0.052-0.063
published). Chasing this down turned up a real bug. The bootstrap DGP fitted the
wrong AR order (`AR(lag+1)` instead of `AR(lag)`, and an unwanted `AR(1)` at
`lag=0` instead of step 1 of Pavlidis, which reduces to a plain i.i.d. resample
when `q=0`). The bug came from adapting the specific convention of `radf_sb_()`
without re-deriving it from the formula of Pavlidis. We fixed it, but
re-validating after the fix made things worse and not better (0.10-0.20
oversized), so the first fix was not the whole story.

A decoupled diagnostic isolated the real issue. We compared the oracle
finite-sample null distribution of `Un` and `QKS` (1,000 genuine i.i.d. random
walks, with no bootstrap at all) against the critical values implied by the
bootstrap itself (one series, `nboot` up to 1,999). The bootstrap critical value
stayed substantially and persistently below the oracle at `tau = 0.5/0.8/0.9`,
even at `nboot = 1999`. That rules out "just needs more bootstrap replications",
which would shrink toward the oracle and would not plateau below it. Only
`tau = 0.95` and `QKS` came close to the oracle at large `nboot`. This is a
structural calibration problem in the bootstrap procedure, at low and middle
quantiles specifically. It is neither sampling noise nor the AR-order bug, which
we had already fixed by then. It is likely related to the well-known extreme
sensitivity of `Un(tau) = n*(alpha1_hat(tau) - 1)` to small numerical
differences in `alpha1_hat`. The coefficient is near 1 for a unit root, so
`n ≈ 99` amplifies third-decimal-place differences by about 100 times. We did
not pin down the exact mechanism.

It is not shipped. Our standing rule is not to ship code that fails its own
validation, so we removed `radf_qar()` and its tests and did not commit them. We
record the state precisely so that a future pass starts from here and not from
scratch. The AR-order fix was correct but insufficient, and the remaining gap is
at low and middle `tau`, while `QKS` and `tau=0.95` look calibrated even at large
`nboot`. A revisit could do one of two things. It could investigate the
`alpha1_hat` sensitivity directly, comparing the bootstrap and oracle
distributions of `alpha1_hat` (and not just the `n*(...)`-transformed statistic)
at several `tau`. Or it could restrict itself to `QKS` alone, the headline
statistic recommended by the paper and the one that validated cleanly, instead
of exposing the individual `Un(tau)` statistics that did not.

## Noncausal / local explosive dynamics

Status: evaluated, not implemented (2026-08-09). We read the PDF (abstract and
introduction).

### Source

Blasques, F., Koopman, S.J., Mingoli, G. & Telg, S. (2025). "A Novel Test
for the Presence of Local Explosive Dynamics." *JTSA*, 46(5), 966-980.
`doi:10.1111/jtsa.70001`.

### What it is

This is again a different paradigm. The test is built for mixed
causal-noncausal autoregressive processes, a model class where part of the
dynamics depends on future shocks (anticipative or noncausal terms) as well as
past ones. The premise is that bubbles come from an extreme shock acting
through the forward-looking (noncausal) component of the model, and do not
emerge from a recursively estimated explosive AR root on the past alone. The
distribution of the test statistic is either derived analytically or
approximated numerically, depending on the assumed error distribution. The
empirical application is a monthly oil price index, framed partly as a
Value-at-Risk-style risk-assessment tool and not purely as a bubble detector.

### Cost/feasibility note for exuber

On reading, the method is out of the architecture of exuber, and we did not
conclude this from the abstract alone. Mixed causal-noncausal AR models need
their own specialized estimation. There is no closed-form OLS or QR reduction,
and noncausal components are typically estimated by approximate or simulated
maximum likelihood under a specified non-Gaussian error distribution, because
noncausal processes are identifiable only when the innovations are non-Gaussian.
None of the recursive-least-squares machinery of `exubercore` applies, nor does
any transform-then-reuse trick that worked for STADF, the sign-based test or PDC
above. This would be a from-scratch statistical framework with a new estimation
dependency (a noncausal-AR fitting routine, and we identified no mainstream R
package for it in this pass). It is not a contained addition under any framing,
and it has the lowest priority of the three live items in this file.

## Spectral fragility (out of scope)

Bhandari, A. "Rational Bubbles at the Spectral Edge: An Operator-Spectral
Theory of Fragility, Identification and Finite-Sample Certification."
arXiv:2607.03933.

This is factor and co-movement spectral-fragility detection. It identifies
market fragility through the strength of a dominant factor extracted from
cross-sectional co-movement (fewer independent factors during crises than in
calm periods), and it is not a right-tailed unit-root test on a single series at
all. It detects fragility contemporaneously and not predictively. It is not
built on ADF/SADF machinery in any way, so it would need a from-scratch
implementation with no reuse of the existing code of exuber. It also arguably
answers a different question (market-wide fragility, and not the explosiveness of
a specific series) from the one `radf()` and its relatives address. We note it
here so that nobody rediscovers and re-evaluates it later. It is not a live
candidate.

## Stochastic tree asset pricing (out of scope)

Gourieroux, C. & Jasiak, J. (2025). "A Stochastic Tree for Bubble Asset
Modelling and Pricing." *JTSA*, 46(5), 932-944.

The paper presents an asset-pricing model, a stochastic tree representation for
modelling, forecasting and pricing bubbles, with closed-form option-pricing
formulas. It is not a test for the presence of a bubble. The scope of exuber is
testing (`radf()` and its relatives test for explosiveness), so the paper does
not fit however we implement it. We note it for awareness and it is not a
candidate.
