---
title: "Alternative paradigms"
blurb: "Non-ADF-family approaches, principally the quantile-based global test and its recursive monitoring extension."
order: 5
---
This file covers methods that address the same problem, detecting explosive or bubble dynamics, from outside the ADF/SADF/GSADF/BSADF recursive-regression family around which exuber is built.

| Method | Paper | Status |
|---|---|---|
| [Quantile-based detection](#quantile-based-detection) | Pavlidis (2025); Wu, Shi & Wu (2025) | the global test, `QPWY` and `QPSY` monitoring are done (`QPSY` carries a small-sample caveat away from the median); Pavlidis's `Un`/`QKS` is not implemented |
| [Noncausal / local explosive dynamics](#noncausal--local-explosive-dynamics) | Blasques, Koopman, Mingoli & Telg (2025) | evaluated, not implemented |
| [Spectral fragility](#spectral-fragility-out-of-scope) | Bhandari (arXiv) | out of scope |
| [Stochastic tree asset pricing](#stochastic-tree-asset-pricing-out-of-scope) | Gourieroux & Jasiak (2025) | out of scope (pricing, not testing) |

All papers are from the *JTSA* 46(5) special issue except Bhandari (arXiv). See [references.md](/replication/references).

## Quantile-based detection

Status: the global test is `quantile_test()`, and the recursive monitoring extensions `QPWY` and `QPSY` are `monitor_quantile()`.

### Sources

- Pavlidis, E. G. (2025). Bubbles and crashes: A tale of quantiles. *JTSA*, 46(5), 884–907.
- Wu, R., Shi, S. & Wu, J. (2025). Quantile analysis for financial bubble detection and surveillance. *JTSA*, 46(5), 908–931. Shi is a coauthor of PSY, so this is a quantile-based test from inside the same lineage.

### Idea

These papers characterise explosive behaviour through quantile regression (QR) at chosen points of the conditional distribution, in place of a recursive mean-regression ADF statistic.

The global test of Wu, Shi & Wu (eq. 18) is the QR analogue of the DF $t$-ratio:

$$
t_T(\tau) = \frac{\hat f(b_\tau)}{\sqrt{\tau(1-\tau)}}\, \bigl(Y_{-1}' P_Z Y_{-1}\bigr)^{1/2} \bigl(\hat\alpha(\tau) - 1\bigr).
$$

Here $\hat\alpha(\tau)$ is the quantile-regression estimator of $y_t$ on an intercept and $y_{t-1}$ at quantile $\tau$ (eq. 13). The term $Y_{-1}' P_Z Y_{-1}$ is the demeaned sum of squares of the lagged level ($P_Z$ is a demeaning projector and $Z$ a column of ones). The density $\hat f(b_\tau)$ is a kernel estimate of the density of the first-differenced series at its own $\tau$-th sample quantile (eq. 19, $\hat f(b_\tau) = (Th)^{-1} \sum_t K\bigl((\hat b_\tau - \hat u_t)/h\bigr)$ with $\hat u_t = y_t - y_{t-1}$). It studentises the QR coefficient, as the residual variance does for the OLS DF ratio. The paper's `QPWY` and `QPSY` statistics have the sup-scan structure of PWY and PSY, with this QR $t$-ratio computed in every window.

**Critical values** (eq. 22–23). The limiting null distribution of $t_T(\tau)$ is

$$
U(\tau) = \sqrt{1 - \delta(\tau)^2}\; z + \delta(\tau)\, Q, \qquad z \sim N(0,1),
$$

where $\delta(\tau)$ is a correlation estimated from the data, $\operatorname{cor}\bigl(\hat u_t,\ \tau - \mathbf 1\{\hat u_t < \hat b_\tau\}\bigr)$, between the innovation and its quantile-check score. $Q = (\int \bar W^2)^{-1/2} \int \bar W\, dW$ (eq. 23) is the standard demeaned Dickey–Fuller $t$-distribution. Simulating $Q$ by the random-walk-plus-OLS-$t$ construction of `radf_mc_cv()` gives values identical, bit for bit, to the `adf` field of `radf()` on the same series. The quantiles of $U(\tau)$ therefore need one fresh normal draw and a critical value that the package already simulates.

**Optimal quantile.** The paper's criterion (eq. 33) is $\tau^* = \arg\min_\tau \tau(1-\tau)/\hat f(b_\tau)^2$, a grid search over the same $\hat f(b_\tau)$.

### Implementation: global test

`quantile_test()` implements the global test of Section 3.1, a single static QR fit at one quantile with no recursion. It adds `quantreg` (>= 5.9) as an estimation dependency of exuber. `tau = "optimal"` (the default) runs the grid search of eq. 33 over `tau_grid` (default `seq(0.2, 0.8, by = 0.05)`, the practical range the paper recommends), and a fixed `tau` can be passed. `quantreg::rq()` fits the regression, and the density, the demeaned sum of squares and the critical-value simulation are plain R.

**Validation.** The empirical size under a pure random-walk null is 5.0% for `tau = "optimal"` and 3.0% for fixed `tau = 0.5` (nominal 5%, 100 replications). Power under an explosive alternative is 100%, as for SADF on the same DGP. The selected $\tau$ lies inside the search grid across replications and does not collapse to a boundary. Replication script: [replication/alternative-paradigms/radf_quantile_validation.R](#script-radf_quantile_validation).

### Implementation: QPWY

$\mathrm{QPWY}_r(\tau) := t_T^{0,r}(\tau)$ is a single recursion. The window start is fixed at 1 and only the end $r$ grows, which is the `badf` shape of `radf()`, so it needs $O(T)$ QR fits. QR has no closed form in cumulative sums, so each window needs a genuine fit (Corollaries 1–2, pages 10–11).

Corollary 1 decomposes the limit of the windowed statistic as

$$
U'^{r_1, r_2}(\tau) = \sqrt{1 - \delta(\tau)^2}\; z + \delta(\tau)\, Q_{r_1, r_2},
$$

the same decomposition used for `quantile_test()`. Corollary 2 identifies $Q_{0,r}$ with the `badf` sequence of `radf()` under a simulated null path.

`monitor_quantile(data, tau = 0.5, minw, nrep, level, seed)` is in `exuber/R/monitor_quantile.R`. `qpwy_stat_path()` is the $O(T)$ loop of `quantreg::rq()` fits and repeats the per-window $t$-ratio of `quantile_test()` (eq. 18) over a growing window. The boundary is a single flat value, the quantile across replicates of each simulated path's own maximum, as in the `sadf_cv` of `radf_mc_cv()`. A pointwise marginal quantile at each $r$ would not control the first crossing: it gave a false-alarm rate of 50% against a nominal 5%.

The limit in Theorem 1 is $\int \tilde W\, dB_\psi / \sqrt{\int \tilde W^2}$, where $B_\psi$ is a Brownian motion with correlation $\delta$ to $W$. Writing $B_\psi = \delta W + \sqrt{1 - \delta^2}\, V$ with $V$ independent of $W$ gives

$$
\delta\, Q_{r_1, r_2} + \sqrt{1 - \delta^2}\, Z_{r_1, r_2}, \qquad Z = \frac{\int \tilde W\, dV}{\sqrt{\int \tilde W^2}} .
$$

For a single window $Z$ is exactly $N(0,1)$, which is all that `quantile_test()` needs. It varies across windows, however, and a monitoring boundary is a quantile of path suprema, so a single constant draw of $z$ understates it. At $n = 200$ (4000 replications, Gaussian innovations, $\tau = 0.5$) the single-$z$ 95% boundary and the correct one are:

| $\delta$ | single-$z$ boundary | correct boundary | true size of single-$z$ |
|---|---|---|---|
| 0.8 | 1.534 | 1.665 | 0.066 |
| 0.5 | 1.657 | 2.060 | 0.108 |
| 0.2 | 1.704 | 2.405 | 0.196 |

The distortion is largest where `QPWY` is meant to help, with heavy tails and non-central quantiles, because $\delta$ is small there. `quantile_boundary_sim()` simulates $Q$ and $Z$ jointly for every window from prefix sums of $(e_t, v_t)$, at $O(1)$ per window, takes each path's supremum for the data-estimated $\delta$ and matches a per-window brute force to $2 \times 10^{-15}$.

**Validation.** The false-alarm rate under $H_0$ at a nominal 5% ($n = 150$, 200 replications) is 0.035 for Gaussian innovations and for $t_3$ at $\tau = 0.5$. For $t_3$ at $\tau = 0.2$, $0.8$ and $0.9$ it is 0.075, 0.085 and 0.125. The last is well above nominal, in line with the paper's advice to avoid extreme quantiles in small samples (its Table II reports 0.07 for the global test at $\tau = 0.9$ under $t_3$). Power on a post-training explosive DGP (60 replications) is 50.0%, against 55.0% for SADF on the same DGP. `QPWY` gives up a little power for robustness to non-Gaussian innovations, which is the paper's motivation. `test-qpwy.R` has 7 tests, including a check that the supremum-calibrated boundary stochastically dominates any single-column marginal quantile.

### Implementation: QPSY

$\mathrm{QPSY}_r(\tau, r_0) = \sup_{r_1} t^{r_1, r}(\tau)$ (eq. 26) needs $O(T^2)$ QR fits, about $T^2/2$ per series. In R with `rq.fit(method = "br")` that takes 2.7 s at $n = 100$ and 12 s at $n = 200$. The boundary needs none of them. The $\sup_{r_1} [\delta\, Q_{r_1, r} + \sqrt{1 - \delta^2}\, Z_{r_1, r}]$ of Corollary 2 comes from the same $Q$/$Z$ simulation as `QPWY`, run over the full $(r_1, r)$ grid instead of the $r_1 = 0$ row.

`monitor_quantile(..., type = "qpsy")` uses windows $[r_1, r]$ with at least `minw` regression observations (the same first-window floor as `QPWY`) and the same flat, supremum-calibrated boundary.

**Validation.** `qpsy_stat_path()` matches a `quantreg::rq()` brute force over every window exactly. The grid simulation matches a per-window brute force to $2 \times 10^{-15}$, and the `QPSY` suprema dominate those of `QPWY` on shared draws, as they must. At a nominal 5% ($n = 100$, 80 to 100 replications), Gaussian innovations give size 0.040 at $\tau = 0.5$ and 0.350 at $\tau = 0.9$. With $t_3$ innovations the size is 0.037 at $\tau = 0.5$, 0.212 at $\tau = 0.8$ and 0.440 at $\tau = 0.9$. For power ($n = 100$, explosive from $t = 71$, $\rho = 1.04$, Gaussian innovations, 60 replications), `QPWY` gives 0.400, `QPSY` 0.433 and SADF 0.483. The OLS test is ahead with Gaussian errors, as in Table V of the paper.

**Caveat.** The asymptotic boundary is exact and well sized at the median, and does not hold in the small early windows (about 20 observations) away from the median. A double supremum over thousands of such windows amplifies the finite-sample error, which the single supremum of `QPWY` mostly avoids. The paper uses bootstrap critical values for monitoring (Algorithm 1, used for Table V) and advises against extreme quantiles. With `tau` away from 0.5, `type = "qpsy"` emits a caveat as a message and as `attr(x, "caveat")`, and `?monitor_quantile` gives the numbers. pyexuber has `monitor_quantile(type="qpsy")` with a `UserWarning` for the caveat.

The bootstrap of Algorithm 1 is not implemented. It resamples the centred $\Delta y$ and recomputes the whole statistic path for every replicate, so each replicate repeats the full $O(T^2)$ QR sweep. That takes about 9 minutes per series at $n = 100$ with 199 replicates.

Tests are in `test-qpwy.R`. Replication script: [replication/alternative-paradigms/radf_qpwy_validation.R](#script-radf_qpwy_validation). The sizes are Monte Carlo estimates from 80 to 200 replications, good to about 2 or 3 points.

### Pavlidis's quantile-autoregressive tests: not implemented

Pavlidis characterises bubbles through unit-root quantile-autoregressive models in which the largest autoregressive root may vary by quantile (below 1 at low quantiles and crashes, above 1 at high quantiles and expansions). Pages 6–7 (eq. 5–11) give the same ADF regression form as `radf()`, fitted by quantile regression at chosen $\tau$. The statistics are $U_n(\tau) = n\,(\hat\alpha_1(\tau) - 1)$ (coefficient-based) and $\mathrm{QKS} = \sup_{\tau \in \mathcal T} U_n(\tau)$ (eq. 11). Critical values come from a residual or sieve bootstrap (page 7, steps 1–5) that is close to the Pedersen–Schütte bootstrap of `radf_sb_()`: fit an AR($q$) to $\operatorname{diff}(y)$ under $H_0$, resample the centred residuals and rebuild the series by cumulation. Table 2 of the paper gives empirical sizes (N(0,1), $t_3$ and $t_2$ errors, $n = 100, 200, 300, 400$).

The tests are not implemented because the bootstrap does not reproduce the published sizes. In a prototype, the size of $U_n(\tau = 0.5)$ at $n = 100$ with N(0,1) errors matched Table 2 (0.050 against 0.053), but $U_n$ at higher $\tau$ and $\mathrm{QKS}$ were oversized (0.075–0.100 against 0.052–0.063 published). Comparing the oracle null distribution of $U_n$ and $\mathrm{QKS}$ (1000 i.i.d. random walks, no bootstrap) with the critical values implied by the bootstrap on a single series, the bootstrap critical value stays substantially below the oracle at $\tau = 0.5, 0.8, 0.9$ even at 1999 bootstrap replications, so the gap is not a matter of replication count. Only $\tau = 0.95$ and $\mathrm{QKS}$ came close to the oracle at large `nboot`. The statistic $U_n(\tau) = n(\hat\alpha_1(\tau) - 1)$ is very sensitive to small differences in $\hat\alpha_1$, because the coefficient is near 1 and $n \approx 99$ amplifies third-decimal differences by about 100. A next step would be to compare the bootstrap and oracle distributions of $\hat\alpha_1$ itself at several $\tau$, or to implement $\mathrm{QKS}$ alone, the headline statistic of the paper, which calibrated well.

## Noncausal / local explosive dynamics

Status: evaluated, not implemented.

### Source

Blasques, F., Koopman, S. J., Mingoli, G. & Telg, S. (2025). A Novel Test for the Presence of Local Explosive Dynamics. *JTSA*, 46(5), 966–980, `doi:10.1111/jtsa.70001`.

### Idea

The test is built for mixed causal-noncausal autoregressive processes, a model class in which part of the dynamics depends on future shocks (anticipative or noncausal terms) as well as past ones. The premise is that bubbles come from an extreme shock acting through the forward-looking component of the model, and not from a recursively estimated explosive AR root on the past alone. The distribution of the test statistic is derived analytically or approximated numerically, depending on the assumed error distribution. The application is a monthly oil price index, framed partly as a Value-at-Risk-style risk-assessment tool.

### Fit with exuber

Mixed causal-noncausal AR models need their own estimation. There is no closed-form OLS or QR reduction. Noncausal components are typically estimated by approximate or simulated maximum likelihood under a specified non-Gaussian error distribution, because noncausal processes are identifiable only with non-Gaussian innovations. None of the recursive least squares of `exubercore` applies, nor does any transform-and-reuse approach of the kind used for STADF, the sign-based test or PDC. It would be a separate statistical framework with a new estimation dependency, and no mainstream R package for noncausal AR fitting is available.

## Spectral fragility (out of scope)

Bhandari, A. Rational Bubbles at the Spectral Edge: An Operator-Spectral Theory of Fragility, Identification and Finite-Sample Certification. arXiv:2607.03933.

The paper detects factor and co-movement spectral fragility. It identifies market fragility through the strength of a dominant factor extracted from cross-sectional co-movement (fewer independent factors during crises than in calm periods). It is not a right-tailed unit-root test on a single series, and it detects fragility contemporaneously and not predictively. It is not built on ADF machinery, and it addresses market-wide fragility and not the explosiveness of a specific series, so it falls outside the scope of exuber.

## Stochastic tree asset pricing (out of scope)

Gourieroux, C. & Jasiak, J. (2025). A Stochastic Tree for Bubble Asset Modelling and Pricing. *JTSA*, 46(5), 932–944.

The paper presents a stochastic-tree representation for modelling, forecasting and pricing bubbles, with closed-form option-pricing formulas. It is an asset-pricing model and not a test for the presence of a bubble, and exuber is concerned with testing.
