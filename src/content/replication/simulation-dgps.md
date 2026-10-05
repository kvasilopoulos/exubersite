---
title: "Simulation DGPs"
blurb: "Data-generating processes for the axes the original sim_*() functions do not cover."
order: 6
---
This file differs from the other family files, because it catalogues data-generating processes (DGPs) for simulating bubble series and not tests or statistics. The DGPs come from the Monte Carlo sections of the papers in this project, plus a few papers added for this purpose. The question the file answers is what the `sim.R` of exuber (`sim_psy1`, `sim_psy2`, `sim_ps1`, `sim_ps2`, `sim_blan`, `sim_evans`, `sim_div`) does not cover.

The seven existing `sim_*` functions share one property: fixed (homoskedastic), i.i.d. Gaussian innovations. They differ in the mean equation (PSY-style regime-switching AR(1) against Blanchard and Evans rational bubbles) and in the number and shape of bubbles and collapses. None has time-varying volatility, GARCH, non-Gaussian innovations, stochastic switch timing or a multi-series or factor structure. A DGP is catalogued below only if it differs on one of those axes through a different mechanism, and not through a reparameterisation of an equation that already exists (see the exclusions at the end).

All 15 catalogued DGPs are implemented in `exuber/R/sim.R` and tested in `exuber/tests/testthat/test-sim-dgp.R` (39 assertions). The tests use formula-exact checks against independent brute-force reimplementations where feasible, checks of published properties (for example the price floor of `sim_tree()` and the empirical transition rate of `sim_msbubble()`), and moment and reproducibility checks otherwise.

Items 1, 2 and 3 compose as optional arguments of `sim_psy1()` (`e`, `shifts`, `coef_noise`). Five small, reusable generators cover items 8 to 12 (stochastic volatility, innovation distribution, level shift, long memory and stochastic coefficient): `sim_innov()`, `sim_vol_garch()`, `sim_vol_cir()`, `sim_vol_sv()` and `sim_fi()`. `sim_blan()` has a `type = "rotermann_wilfling"` option (item 13). The six DGPs with their own architecture (items 4 to 7, 14 and 15) are standalone functions: `sim_common()`, `sim_coexplosive()`, `sim_tree()`, `sim_mar()`, `sim_msbubble()` and `sim_falsebubble()`.

A DGP does not imply that its paired test exists. The level-shift robustness result of Harvey, Leybourne, Tatlow & Zu (2025), which motivates item 3, is covered by `radf_sign()` and `radf_sign_dm()` (see [volatility-robustness.md](/replication/volatility-robustness#level-shift-robustness-hltz-2025)). `sim_falsebubble()` and `sim_msbubble()` have no dedicated test in exuber. They are stress-test and demonstration series for the existing PSY/GSADF machinery, the role that `sim_evans()` already plays.

Two implementation notes. `sim_fi()` convolves the truncated MA($\infty$) filter by hand over `m` extra truncation-lag innovations. `stats::filter(..., sides = 1)` needs strictly more input than filter taps to return non-`NA` values, and `stats::filter()` and `stats::convolve()` also crashed on the R installation used for this project (Windows, R 4.6.1) independently of exuber. `sim_tree()` clips $p_t = \Phi(X_t)$ into $[10^{-10}, 1 - 10^{-10}]$, because $p_t$ can round to exactly 0 or 1 in the tail of a long series, which would turn $\xi_{1t}$ and $\varepsilon_t$ into $0/0$.

## DGPs from papers already in the project

| # | DGP | Source | Single/multi | Implementation |
|---|---|---|---|---|
| 1 | CIR-type stochastic volatility | Harvey, Leybourne & Zu (2019) | single | `sim_vol_cir()` |
| 2 | AR(1) lognormal stochastic volatility, near-unit persistence | Sarkar & Wells (2025/2026) | single | `sim_vol_sv()` |
| 3 | Deterministic level shifts (mean jumps) | Harvey, Leybourne, Tatlow & Zu (2025) | single | `sim_psy1(..., shifts = ...)` |
| 4 | Latent common-factor bubble | Chen, Phillips & Shi (2023) | common, multi-series | `sim_common()` |
| 5 | Bivariate co-explosive linkage | Evripidou, Harvey, Leybourne & Sollis (2022) | bivariate | `sim_coexplosive()` |
| 6 | Stochastic branching-tree (random-coefficient RCA) | Gourieroux & Jasiak (2025) | single | `sim_tree()` |
| 7 | Mixed causal-noncausal AR (MAR), heavy-tailed | Blasques, Koopman, Mingoli & Telg (2025) | single, self-terminating | `sim_mar()` |
| 8 | Fractionally integrated (long-memory) innovations | Lui, Phillips & Yu (2024) | single | `sim_fi()` |
| 9 | PSY equation with heavy-tailed or skewed innovations | Wu, Shi & Wu (2025) | single | `sim_innov(dist = "t"/"skew_t")` |
| 10 | GARCH(1,1) and logistic smooth-transition volatility | Whitehouse, Harvey & Leybourne (2025); Harvey, Leybourne, Taylor & Zu (2024) | single | `sim_vol_garch()` |
| 11 | Stochastic explosive coefficient (persistence itself random) | Kurozumi & Nishi (2025) | single | `sim_psy1(..., coef_noise = ...)` |

### 1. CIR-type stochastic volatility

Harvey, Leybourne & Zu (2019) use a square-root (Cox–Ingersoll–Ross) diffusion for volatility, "representative of Bollerslev and Zhou (2002)", as a robustness design beyond their main deterministic-volatility one:

$$
d\sigma^2(r) = 0.03\,\bigl(0.25 - \sigma^2(r)\bigr)\,dr + 0.1\,\sigma(r)\,dB(r).
$$

Each replication simulates it from NIID(0,1) Brownian increments, independent of the level innovations, and feeds the same PSY-style level process that `sim_psy1` implements. The new element is continuous-time stochastic volatility, where `sim.R` has fixed or deterministic volatility.

### 2. AR(1) lognormal stochastic volatility

Sarkar & Wells (2025, Section 2) reused in Sarkar & Wells (2026, eq. 7):

$$
y_t = \rho_n\, y_{t-1} + u_t, \qquad u_t = \sigma_t\, \varepsilon_t, \qquad \log \sigma_t^2 = \phi_n \log \sigma_{t-1}^2 + \eta_t .
$$

Both $\rho_n \to 1$ and $\phi_n \to 1$, which makes the process double local-to-unity in the mean and in the variance. It is a discrete-time stochastic volatility model with a persistent, near-integrated log-variance. No existing `sim_*` function has a persistent variance process.

### 3. Deterministic level shifts

Harvey, Leybourne, Tatlow & Zu (2025):

$$
y_t = y_{t-1} + \sum_{j=1}^{m} \delta_j\, \mathbf 1(t \ge \tau_j) + \varepsilon_t \qquad \text{(null)}.
$$

The process has an arbitrary number $m$ of jumps of size $\delta_j$ at unknown dates $\tau_j$, with an explosive regime layered on top for the alternative (Section 5.2). Nothing in `sim.R` has jump components.

### 4. Latent common-factor bubble

Chen, Phillips & Shi (2023) eq. 2.3 and 2.7–2.9:

$$
X_t = \Lambda f_t + e_t .
$$

There are $N$ observed series and one latent factor $f_t$ that follows a PSY-style unit-root, explosive, collapse sequence. The loadings are $\Lambda \sim U[0,2]$, and the idiosyncratic noise has $\sigma_e = 0.1$. A collapse-splicing construction (paper lines 960–977) avoids discontinuities at the regime boundary. The new element is one common bubble that drives many series jointly. Every existing `sim_*` function is single-series.

### 5. Bivariate co-explosive linkage

Evripidou, Harvey, Leybourne & Sollis (2022) The series $x_t$ is generated from the regime-dummy Models 1–4 of HLS (each individually PSY-style), and a second series is linked to it:

$$
y_t = \mu_y + \phi_x\, x_{t-i} + \phi_z\, z_t + \varepsilon_{y,t} .
$$

The second series combines a lead or lagged copy of the explosive series $x_t$ with a third, latent explosive series $z_t$, and the heteroskedasticity is timed to the regime changes. The new element is a two-series lead and lag linkage of explosive components. It differs from the common-factor case (item 4), which shares one factor and does not link two distinct explosive series.

### 6. Stochastic branching-tree bubble

Gourieroux & Jasiak (2025) describe an affine autoregression with a stochastic coefficient. It can be represented as a random-coefficient AR process generated by a binomial tree with stochastic branching intensity, in contrast to the deterministic branches of Cox–Ross–Rubinstein. The Blanchard & Watson (1982) bubble, already `sim_blan`, is the special case of constant intensity. The branching mechanism generates the path directly, with no fixed collapse probability.

### 7. Mixed causal-noncausal AR (MAR)

Blasques, Koopman, Mingoli & Telg (2025):

$$
(1 - \varphi_1 L)(1 - \psi_1 L^{-1})\, y_t = \varepsilon_t .
$$

The causal root is $\varphi_1 = 0.7$, $\psi_1$ is the noncausal root, and the innovations are Cauchy or Student-$t(2)$. The noncausal component generates transient, self-terminating local bubbles with no scripted regime dates. The process differs on every axis: a lag-polynomial form in place of a regime-dummy AR, an implicit collapse and heavy-tailed non-Gaussian innovations.

### 8. Fractionally integrated innovations

Lui, Phillips & Yu (2024):

$$
y_t = y_{t-1} + u_t, \qquad u_t = \Delta^{-d} \varepsilon_t \quad (d > 0,\ \varepsilon_t \text{ i.i.d. with finite } (2+\delta) \text{ moments}).
$$

The innovations are $\mathrm{FI}(d)$ (long memory) and not i.i.d. The paper develops an explosive-alternative analogue in Section 4. The new element is long-range-dependent noise in the unit-root or explosive equation.

### 9. PSY equation with heavy-tailed or skewed innovations

Wu, Shi & Wu (2025) eq. 6, use the standard PSY unit-root, explosive, collapse equation with innovations from $N(0,1)$, $t(3)$, skewed-$t(3, -0.75)$ and skewed-$t(3, +0.75)$. Only the innovation distribution is new. Every existing `sim_*` function uses fixed Gaussian noise.

### 10. GARCH(1,1) and smooth-transition volatility

Whitehouse, Harvey & Leybourne (2025):

$$
z_t = h_t^{1/2}\, \varepsilon_t, \qquad h_t = 0.1 + 0.1\, z_{t-1}^2 + 0.8\, h_{t-1} .
$$

Harvey, Leybourne, Taylor & Zu (2024) use the same GARCH(1,1) specification. Next to their main design, they use a logistic smooth-transition volatility function:

$$
\sigma(r) = \sigma_1 + \frac{\sigma_2 - \sigma_1}{1 + \exp\{-\kappa (r - \delta)\}} .
$$

The new elements are conditional GARCH heteroskedasticity and a smooth transition between volatility regimes in place of an instant jump.

### 11. Stochastic explosive coefficient

Kurozumi & Nishi (2025) is documented in [volatility-robustness.md](/replication/volatility-robustness#stochastic-explosive-coefficient-test). The coefficient $1 + c_1/T + a\, u_t/\sqrt T$ replaces the deterministic $1 + c/T^\alpha$, so the persistence parameter itself is random and not only the noise scale.

## DGPs from additional papers

### 12. TGARCH(1,1) with leverage effect

Monschang, V. & Wilfling, B. (2021). Sup-ADF-style bubble-detection methods under test. *Empirical Economics*, 61, 145–172, `doi:10.1007/s00181-020-01859-7`. Open access (CQE Working Paper 78/2019):

$$
\varepsilon_t = s_t\, h_t^{1/2}, \qquad h_t = \omega + \alpha\, \varepsilon_{t-1}^2 + \beta\, h_{t-1} + \gamma\, \varepsilon_{t-1}^2\, \mathbf 1(\varepsilon_{t-1} < 0) .
$$

The process is calibrated to NASDAQ estimates ($\alpha = 0.4387$, $\gamma = 0.1306$, $\beta = 0.9319$). Relative to item 10, the new element is an asymmetric (sign-dependent) shock response, the leverage effect.

It is `sim_vol_garch(omega, alpha, beta, gamma)`, the same function as item 10 with `gamma` as the leverage parameter. `sim_vol_garch(omega = 0.4387, alpha = 0, beta = 0.9319, gamma = 0.1306)` reproduces the NASDAQ calibration of the paper, and the default `gamma = 0` gives plain GARCH(1,1).

### 13. Lognormal-mixture rational bubble (Rotermann–Wilfling)

The same paper, eq. 4:

$$
B_{t+1} =
\begin{cases}
B_t\, u_t / \delta & \text{with probability } \pi,\\[2pt]
\dfrac{1 - \pi\delta}{1 - \pi}\, B_t\, u_t & \text{with probability } 1 - \pi,
\end{cases}
\qquad u_t \overset{\text{iid}}{\sim} \text{lognormal}.
$$

It produces recurring, stochastically deflating trajectories and no single full collapse to a fixed floor, unlike `sim_blan` and `sim_evans`. The new element is partial, probabilistic deflation in place of total collapse to noise.

It is `sim_blan(type = "rotermann_wilfling", delta, rw_sigma)`, a branch of the existing function, because it shares the two-regime structure with probability $\pi$ and only the update rule per regime differs. `test-sim-dgp.R` verifies that it stays strictly positive, which is a structural invariant of the multiplicative recursion.

### 14. Markov-switching present-value bubble

Chan, J. C. C. & Santi, C. (2021). Speculative Bubbles in Present-Value Models: A Bayesian Markov-Switching State Space Approach. *Journal of Economic Dynamics and Control*, 127, 104101. Open access (author's site):

$$
b_t = \frac{1}{\lambda_{S_t + 1}}\, b_{t-1} + \varepsilon_{bt}, \qquad S_t \in \{1, 2\} \text{ first-order Markov with transition probabilities } p_{11}, p_{22}.
$$

Regime 1 is "surviving" ($\lambda < 1$, explosive) and regime 2 is "collapsing" ($\lambda > 1$, mean-reverting). The bubble sits inside a full present-value state-space model with time-varying expected returns and dividend growth, and the section "Simulated Datasets" of the paper (Table 2) generates artificial data from these parameters. The new element is that the timing of the switch is itself stochastic (a Markov chain). The PSY-style DGPs in `sim.R` script the dates as fixed fractions of $n$.

It is `sim_msbubble(p11, p22, lambda1, lambda2, sigma_b)`. It covers only the bubble component $b_t$ and not the full present-value state-space model. Eq. 16 of the source applies the regime coefficient through $S_{t+1}$, the realised regime of the next period, while the implementation uses the contemporaneous $S_t$, which only changes which time step a given draw of $S$ labels. `test-sim-dgp.R` verifies that the empirical self-transition rate of the simulated regime path matches $p_{11}$ and $p_{22}$ within Monte Carlo tolerance.

### 15. Deterministic technology-adoption "false bubble" null

Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). Technology Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI Episodes. arXiv:2604.25826. Open

The paper embeds a hump-shaped (triangular, Gaussian, Beta or Gamma) deterministic technology-adoption shock into the Campbell–Shiller present-value fundamental. The fundamental price is then locally explosive during adoption with no bubble present, and PSY-style tests can reject spuriously. Appendix E.5 extends it to a Bayesian-updated stochastic version. The new element is a no-bubble null DGP with a smooth deterministic drift, where the other functions in `sim.R` do not produce a series of fundamentals that looks locally explosive.

It is `sim_falsebubble(t1, t2, kappa, shape, amplitude, mu, r)`. `shape = "triangular"` reproduces the worked example of eq. 4 exactly, and `shape = "gaussian"` is one alternative from the paper's list of hump-shaped specifications. Beta and Gamma are not included, since the paper's robustness claim is that they make no qualitative difference. The shock is deterministic, so its price contribution is an exact forward-looking discounted sum, $T_t = \sum_{s > t} \beta^{s-t} \tau_s$. The function is a single-shock reproduction of the mechanism and leaves out the DOLS and multi-functional-form robustness machinery. `test-sim-dgp.R` verifies that `amplitude = 0` reduces to the fundamental-price formula of `sim_div()` bit for bit, and that the technology term is exactly zero outside $[t_1, t_2]$.

### Adjacent: not a price-level DGP

Richter, S., Wang, W. & Wu, W. B. (2023). A supreme test for periodic explosive GARCH. *Econometrics* (MDPI); arXiv:1812.03475. Open The paper studies a piecewise or periodic explosive GARCH(1,1) in which the explosiveness lives in the volatility recursion ($\alpha$ and $\beta$ are temporarily driven toward or through the IGARCH boundary $\alpha_\Sigma + \beta_\Sigma \ge 1$) and not in the price level. It is a volatility-bubble reference and is not counted above, because it addresses a different detection problem from `radf()` and `sim.R`.

## Excluded as reparameterisations

- The six Monte Carlo DGPs A–F of Harvey, Leybourne & Whitehouse (2020), with two- and three-bubble sequences. They belong to the same regime-dummy AR(1) family as HLS/PSY, and only the regime count and the parameters differ from `sim_psy2`.
- The confidence-sets paper of Kurozumi & Skrobotov (2025). It has the linear AR(1) with switching coefficient of the other Kurozumi and Skrobotov papers, with a named recovery regime added.
- The noncausal green-bubble paper of Giancaterini, Hecq, Jasiak & Manafi Neyazi (2025), arXiv:2505.14911, which uses the same mixed causal-noncausal mechanism as item 7.

## Considered but not used

- Breitung & Kruse (2013), When bubbles burst: econometric tests based on structural breaks, *Statistical Papers*, 54(4). Not freely available.
- Testing for explosive bubbles in the presence of non-Gaussian conditions, *Economics Letters*, 233 (2023). Not freely available.
- Montanino & De Luca, The Bubble Crash GARCH model, SSRN 5604452. No retrievable PDF.
- Horváth, Trapani & Wang (2024), Sequential Monitoring for Explosive Volatility Regimes, arXiv:2404.17885. Likely overlaps with the RCA-monitoring paper of Horváth & Trapani (2026) by the same authors.
- Lin, Ren & Sornette (2009), the LPPLS finite-time-singularity model, arXiv:0905.0128. It is a curve-fitting paradigm and not a Monte Carlo DGP for stress-testing right-tailed unit-root tests.
