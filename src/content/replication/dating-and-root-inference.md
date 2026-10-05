---
title: "Dating and root inference"
blurb: "Origination, collapse and recovery dates, plus confidence intervals on the explosive root itself."
order: 2
---
Once `radf()` or `datestamp()` has flagged an episode as explosive, two questions follow. **Dating** asks when the episode started, ended or recovered. **Root inference** asks how explosive it was: what the autoregressive root $\rho$ is, what confidence interval it carries and how fast the series was doubling. The two sit in one file because they are consecutive steps in the exuber workflow (`radf()`, then `datestamp()`, then a sub-sample, then refinement). The methods share little statistical machinery. Status labels are the ones used in [volatility-robustness.md](/replication/volatility-robustness).

| Method | Paper | Status |
|---|---|---|
| [SSR/BIC dating, PDC/KS route](#ssrbic-dating-vs-psy-recursive-dating) | Pang, Du & Chong (2021); Kurozumi & Skrobotov (2023) | done |
| [SSR/BIC dating, HLS/HLW route](#ssrbic-dating-vs-psy-recursive-dating) | Harvey, Leybourne & Sollis (2017); Harvey, Leybourne & Whitehouse (2020) | done (single-bubble HLS and multi-bubble HLW) |
| [Root inference (Cauchy CI and normal-t CI)](#root-inference) | Phillips & Magdalinos (2007); Guo, Sun & Wang (2019) | done |
| [Confidence sets for bubble dates](#confidence-sets-for-bubble-dates) | Kurozumi & Skrobotov (2025) | evaluated, not implemented |
| [Improved retrospective dating](#improved-retrospective-dating) | Kejriwal, Nguyen & Perron (2025) | done (single bubble, and multi-bubble dynamic programme `dating_knp(breaks = )`) |
| [WLS dating under time-varying volatility](#wls-dating-under-time-varying-volatility) | Kurozumi & Skrobotov (2023) | done |
| [Reverse-regression recovery dating](#reverse-regression-recovery-dating) | Phillips & Shi (2014/2019) | done, with caveats |

All papers are listed in [references.md](/replication/references#dating-and-root-inference).

---

## SSR/BIC dating vs. PSY recursive dating

PSY dates an episode where the recursive BSADF statistic first crosses, and later re-crosses, a critical value. The methods in this section replace that rule with a model-based one that minimises the sum of squared residuals (SSR) over candidate break dates and chooses the regime structure by BIC. They differ in how they handle several bubbles in one series.

Status: the PDC/KS route is `dating_pdc()`, HLS's single-bubble route is `dating_hls()` and HLW's multi-bubble two-step wrapper is `dating_hlw()`.

### Sources

1. Harvey, D. I., Leybourne, S. J. & Sollis, R. (2017). Improving the accuracy of asset price bubble start and end date estimators. *Journal of Empirical Finance*, 40, 121–138, `doi:10.1016/j.jempfin.2016.11.001` ("HLS").
2. Harvey, D. I., Leybourne, S. J. & Whitehouse, E. J. (2020). Date-stamping multiple bubble regimes. *Journal of Empirical Finance*, 58, 226–246, `doi:10.1016/j.jempfin.2020.06.004` ("HLW").
3. Pang, T., Du, L. & Chong, T. T. L. (2021). Estimating multiple breaks in nonstationary autoregressive models. *Journal of Econometrics*, 221(1), 277–311 ("PDC").
4. Kurozumi, E. & Skrobotov, A. (2023). On the asymptotic behavior of bubble date estimators. *Journal of Time Series Analysis*, 44(4), 359–373 ("KS").

### HLS (2017): single-bubble SSR + BIC dating

HLS defines four regime-structure models for a series $y_t$. Each is a piecewise OLS regression of $\Delta y_t$ on regime dummies and dummy-interacted $y_{t-1}$ (p. 7), with $D_t(a,b) = \mathbf 1(\lfloor aT \rfloor < t \le \lfloor bT \rfloor)$:

$$
\begin{aligned}
\text{Model 1:}\quad & \Delta y_t = \mu_1 D_t(\tau_1,1) + \delta_1 D_t(\tau_1,1)\,y_{t-1} + v_{1t} && \text{unit root, then a bubble to the sample end}\\
\text{Model 2:}\quad & \Delta y_t = \mu_1 D_t(\tau_1,\tau_2) + \delta_1 D_t(\tau_1,\tau_2)\,y_{t-1} + v_{2t} && \text{unit root, bubble, unit root}\\
\text{Model 3:}\quad & \Delta y_t = \mu_1 D_t(\tau_1,\tau_2) + \mu_2 D_t(\tau_2,1) + \delta_1 D_t(\tau_1,\tau_2)\,y_{t-1} + \delta_2 D_t(\tau_2,1)\,y_{t-1} + v_{3t} && \text{unit root, bubble, collapse to the end}\\
\text{Model 4:}\quad & \Delta y_t = \mu_1 D_t(\tau_1,\tau_2) + \mu_2 D_t(\tau_2,\tau_3) + \delta_1 D_t(\tau_1,\tau_2)\,y_{t-1} + \delta_2 D_t(\tau_2,\tau_3)\,y_{t-1} + v_{4t} && \text{unit root, bubble, collapse, unit root}
\end{aligned}
$$

For each model the break fractions $(\hat\tau_1, \hat\tau_2, \hat\tau_3)$ jointly minimise the model's residual sum of squares $\mathrm{SSR}_j$ over all candidate dates that satisfy the ordering and sign constraints. The bubble phase must be upward, for example. Theorem 1 shows that $\lfloor \hat\tau_i T \rfloor - \lfloor \tau_{i,0} T \rfloor \to_p 0$ for each correctly paired DGP and model. With a fixed-magnitude bubble the estimator is therefore consistent for the exact date and not only for the break fraction.

The four models are compared by a BIC whose penalty counts the fitted coefficients plus the estimated break dates (p. 8):

$$
\begin{aligned}
\mathrm{BIC}_1 &= T \ln\{T^{-1}\mathrm{SSR}_1(\hat\tau_1, 1)\} + (2+1)\ln T,\\
\mathrm{BIC}_2 &= T \ln\{T^{-1}\mathrm{SSR}_2(\hat\tau_1, \hat\tau_2)\} + (2+2)\ln T,\\
\mathrm{BIC}_3 &= T \ln\{T^{-1}\mathrm{SSR}_3(\hat\tau_1, \hat\tau_2, 1)\} + (4+2)\ln T,\\
\mathrm{BIC}_4 &= T \ln\{T^{-1}\mathrm{SSR}_4(\hat\tau_1, \hat\tau_2, \hat\tau_3)\} + (4+3)\ln T,
\end{aligned}
\qquad j_{\mathrm{opt}} = \arg\min_j \mathrm{BIC}_j .
$$

In practice (Section 5) HLS impose minimum regime durations, $\tau_1 \ge s$, $\tau_2 - \tau_1 \ge s$ and $\tau_3 - \tau_2 \ge s/2$, with $s = 0.1$ in the $T = 200$ simulations and $s = 0.05$ in the $T = 389$ application. Beyond that the method is a brute-force grid search over one, two or three breakpoints, depending on the model.

### HLW (2020): two-step extension to multiple bubbles

HLS's four-model set grows combinatorially with the number of bubbles, and PSY's end dates are known to be biased late. HLW therefore combine the two (p. 9).

- **Step 1.** Run PSY's GSADF/BSADF detection and dating as it stands, which gives preliminary start and end fractions for each of the $\hat N$ detected bubbles. These split the sample into $\hat N$ disjoint date windows $[s_j, e_j]$. Windows are split at the midpoint between consecutive PSY regimes, with a rule that ensures that a window starts inside a post-explosive (unit-root) regime and never in the middle of a bubble.
- **Step 2.** Apply HLS's Model 1–4 procedure independently within each window. For every window but the last, only Models 2 and 4 are allowed, because a window boundary is a unit-root point and not a sample end.

The authors describe this as a refinement of PSY's output: "we propose a dating methodology based on minimum sum of squared residual estimators and BIC model selection, but using prior information gleaned from the PSY dating procedure as a means of reducing the dimensionality" (HLW, Section 1).

### PDC (2021) and KS (2023): sequential sample-splitting

PDC take the same SSR-minimisation idea for a single bubble episode and make it much cheaper. Their model has three regimes (unit root, explosive, stationary collapse) with breaks $\tau_{1,0} < \tau_{2,0}$. A stochastic-order argument (their Example 3 and Lemmas A.2–A.4) shows that under the bubble DGP the collapse date is identified first, because the drop in SSR at the collapse dominates the drop at the origination. The two breaks can therefore be estimated one after the other and not jointly (PDC p. 9):

$$
\begin{aligned}
\text{Step 1:}\quad & \hat\tau_2 = \arg\min_{\tau \in (0,1)} \mathrm{RSS}_{2,T}(\tau),\\
& \mathrm{RSS}_{2,T}(\tau) = \sum_{t \le \lfloor \tau T \rfloor} \bigl(y_t - \hat\beta_x(\tau)\, y_{t-1}\bigr)^2 + \sum_{t > \lfloor \tau T \rfloor} \bigl(y_t - \hat\beta_3(\tau)\, y_{t-1}\bigr)^2,\\
& \hat\beta_x(\tau) = \frac{\sum_{t \le \lfloor \tau T \rfloor} y_t\, y_{t-1}}{\sum_{t \le \lfloor \tau T \rfloor} y_{t-1}^2},\\
\text{Step 2:}\quad & \text{repeat the one-break minimisation on the left sub-sample } [1, \hat\tau_2 T] \text{ to get } \hat\tau_1 .
\end{aligned}
$$

The model has no intercept and one AR(1) coefficient per regime, unlike HLS's intercept and dummies. Each $\hat\beta(\tau)$ is a ratio of two prefix sums, so the whole $\mathrm{RSS}(\tau)$ curve costs $O(T)$. There is no joint grid search and no BIC step. The algorithm always assumes that the three-regime structure holds.

KS add a fourth regime, a unit-root recovery after the stationary collapse, and reuse the sequential logic for the extra break. They contrast their cost with that of HLS: "we perform the three SSR minimization with one break each with O(T) computations, while Harvey et al. (2017) requires minimizing the three break model over all possible combinations of these breaks" (KS, Section 1). They also suggest running their method inside each HLW date window in place of HLS's joint fit.

### Published numbers

**HLS Table 1** (Nasdaq composite real price index, 1973:2–2005:6, the series of PWY):

| Sample | PSY test | PSY start | PSY end | BIC model | BIC start | BIC end |
|---|---|---|---|---|---|---|
| 1973:2–2005:6 (full) | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2–2000:9 (pseudo-real-time) | 3.07*** | 1998:11 | 2000:9 | 1 | 2000:1 | 2000:9 |
| 1973:2–2000:10 | 3.07*** | 1998:11 | 2000:10 | 1 | 2000:1 | 2000:10 |
| 1973:2–2000:11 | 3.07*** | 1998:11 | 2000:11 | 1 | 1999:12 | 2000:11 |
| 1973:2–2000:12 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |
| 1973:2–2001:1 | 3.07*** | 1998:11 | 2000:12 | 3 | 1998:11 | 2000:9 |

This is the one concrete head-to-head between the two dating rules in HLS. Both agree on the 1998:11 start, and PSY's end date (2000:12) is three months later than the BIC end date (2000:9).

HLS's Monte Carlo comparison (Section 6, Figures 1–3) and the HLW simulations (Section 4, Figures 1–6) are reported only as plots of frequencies. The text states that the BIC method outperforms PSY in finite samples, particularly for the end date, and that PSY estimates "typically fall somewhat later than the true date", but gives no percentages. These results are therefore recorded as qualitative claims.

**HLW Table 1** (BIC model-selection frequencies across six DGPs, A–F):

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

BIC picks the true model most often in every row. The weakest cases are DGP A, regime 2, and DGP C, regime 2, where the correct-model rate is about 44–45%. A reversion to a unit root close to the end of the window is easily missed, as the paper's own caveat says. This is the frequency of correct model identification and not a measure of dating accuracy. The paper notes that choosing the wrong model, for example Model 3 in place of Model 4, still typically gives accurate break estimates.

**KS empirical application** (Section 6, in prose). NASDAQ Composite, monthly, January 1985 to August 2013: $\mathrm{BIC}_4 = 3387.727$ beats $\mathrm{BIC}_3 = 3404.268$ and $\mathrm{BIC}_2 = 3409.296$, so the four-regime model is chosen. The collapse is dated February 2000, the origination (from the left sub-sample) August 1998 and the recovery (right sub-sample) September 2001. US real house price index (FHFA), January 1991 to December 2012: $\mathrm{BIC}_4 = -2918.223$ is chosen again, with collapse in November 2006, origination in September 1997 and recovery in May 2011.

**KS Monte Carlo** (Section 5) is shown as histograms. The prose gives only approximate figures: the frequency of selecting the true break date "is around 30% for T=400 and 65% for T=800", and other quantities are "close to 100%" or "approximately 75% and 100% for T=400 and 800". They should be read as orders of magnitude.

### How the methods relate to `datestamp()`

`datestamp()` (`exuber/R/radf-methods.R`, with the helpers `stamp()`, `stamp_to_index()` and `add_peak()`) post-processes statistics that `radf()` has already computed. It compares the recursive `bsadf` or `badf` sequences with critical values and finds contiguous runs of exceedance. No regression is fitted in that path.

The SSR-based routes are different in kind. They need new regime-dummy or no-intercept AR(1) fitting code and a minimum-regime-duration trimming parameter (`s` in HLS/HLW). The HLS/HLW route also needs a BIC comparison and a search over one to three breakpoints. The PDC/KS route is smaller, because each break is a single $O(T)$ scan of a closed-form ratio of cumulative sums, but it takes the number of regimes (three or four) as given. It cannot tell a bubble that collapses from one that is still running at the sample end. That distinction is what the model selection of HLS/HLW provides.

### Implementation: HLS route

`dating_hls(data, trim = 0.05)` is in `exuber/R/dating_hls.R` and is tested in `exuber/tests/testthat/test-hls.R`. It fits all four models, each by exact SSR minimisation over its candidate breakpoints, subject to the minimum-regime trim and the sign constraints. It selects among them by BIC, $n \log(\mathrm{SSR}/n) + \mathrm{df}\,\log n$ with $\mathrm{df} \in \{3, 4, 6, 7\}$ for Models 1–4. The function returns the selected model and its origination, collapse and recovery dates (`NA` for any that the model does not have), plus the BIC of every candidate.

The sign constraints follow HLW's statement of HLS's models. Model 1 requires $y_T > y_{\tau_1}$. Models 3 and 4 require the fitted peak $y_{\tau_2}$ to exceed both the level at the start of the bubble and the level reached after the collapse regime ($y_T$ for Model 3, $y_{\tau_3}$ for Model 4).

**Algorithm.** The regime dummies of the four models never overlap, so the SSR of any candidate partition is the sum of independent per-segment OLS fits. A segment with no active dummy has no fitted parameters ($\mathrm{SSR} = \sum \Delta y_t^2$), and a segment with an active dummy is an intercept-and-slope fit. Both are closed-form ratios of cumulative sums ($S_x, S_{xx}, S_z, S_{zz}, S_{xz}$), so each candidate breakpoint, pair or triple is evaluated in $O(1)$ from precomputed prefix sums with no repeated `lm()` calls. The search takes under two seconds at $T = 400$, in plain R.

**Validation.** `hls_segment_ssr()` matches the SSR of a brute-force `lm()` fit (tolerance $10^{-8}$) for arbitrary segments, and the full joint three-breakpoint search of Model 4 matches an exhaustive nested-`lm()` search bit for bit on a small synthetic series. The Monte Carlo runs use DGPs built on a large positive base level (`100 + cumsum(rnorm(...))`), so that the explosive signal is a large departure from noise on the absolute scale.

- True Model 4 DGP (unit root, bubble, mean-reverting collapse, unit-root recovery): BIC selects Model 3 or 4 in 100% of 30 replications, split 80% and 20%. Telling a final recovery regime from a continued collapse is the hardest case, which is consistent with the correct-model rate of 44–45% in HLW's weakest DGP. The mean absolute bias is about 0 observations for the origination and about 1 for the collapse.
- True Model 2 DGP: Model 2 is selected in 100% of 30 replications, with origination bias exactly 0.
- True Model 1 DGP: Model 1 is selected in 100% of 30 replications, with origination bias exactly 0.
- Pure $H_0$ (no bubble): BIC never selects Model 4 (0% across 30 replications) and splits between the three simpler models (30% Model 1, 57% Model 2, 13% Model 3). The method has no "no bubble" option, since it is designed to run downstream of a PSY-detected episode.

Replication script: [replication/dating-and-root-inference/radf_hls_validation.R](#script-radf_hls_validation).

### Implementation: HLW route

`dating_hlw(data, cv = NULL, minw = NULL, trim = 0.1, min_duration = NULL, nboot = 199L, seed = NULL, join = 3L)` is in `exuber/R/dating_hlw.R` and is tested in `exuber/tests/testthat/test-hlw.R`. It wraps `dating_hls()` in HLW's two steps.

1. Run PSY's detection and dating (`radf()` and `datestamp()`, with `min_duration` defaulting to `psy_ds(n)`, HLW's own $\ln T$ minimum-episode rule) to get a start and end position for each detected episode.
2. Carve the sample into disjoint date windows. The end of window $j$ is the midpoint between that episode's PSY end and the next episode's PSY start, $e_j = \tau_2^{\mathrm{PSY}}[j] + \lfloor (\tau_1^{\mathrm{PSY}}[j+1] - \tau_2^{\mathrm{PSY}}[j]) / 2 \rfloor$, and the last window runs to the sample end. Each window is fitted with `hls_fit_series()`, the shared helper of `dating_hls()`, restricted to Models 2 and 4 for every window but the last. After window $j$ is fitted, HLW's sequential-adjustment rule moves the start of the next window to the first observation of the fitted post-explosive regime ($s_{j+1} = s_j + \tau_{2,\mathrm{local}}$ for a Model 2 fit and $s_j + \tau_{3,\mathrm{local}}$ for Model 4), so that a later window cannot start inside a bubble.

A series with no detected episode returns an empty result and no error.

**Run-joining.** PSY's detection can split one true bubble into several runs. HLW's rule treats up to three non-rejections surrounded on both sides by an explosive regime of length $\ln T$ as a single episode. `join = 3L` applies it to the regimes from `datestamp()` before the windows are built (`hlw_join_runs()` in R, `_join_runs()` in pyexuber, unit-tested on identical cases). `join = 0` switches it off. Gaps wider than three non-rejections are not joined, by design.

**Validation.**

- On a synthetic two-bubble DGP (20 replications), step 1 finds exactly two windows in 14 of 20 replications with the run-joining rule, and in 13 of 20 without it. In the clean replications the origination and collapse bias for both bubbles is exactly 0. Windows are correctly ordered and do not overlap in all 20.
- Under a pure $H_0$ null, `dating_hlw()` returns no windows in all 20 replications.
- On a single clean bubble (15 replications), the final window matches standalone `dating_hls()` on the whole series exactly (model, origination and collapse) in every replication that has a final window (12 of 15). The paper states that the two-step procedure reduces to HLS when there is a single episode.

Replication script: [replication/dating-and-root-inference/radf_hlw_validation.R](#script-radf_hlw_validation).

### Implementation: PDC/KS route

`dating_pdc(data, regimes = 3L, trim = 0.05)` is in `exuber/R/dating_pdc.R` and is tested in `exuber/tests/testthat/test-pdc.R`. The helper `pdc_find_break(y, trim)` minimises the single-break no-intercept AR(1) RSS from prefix sums of $\sum y_t y_{t-1}$ and $\sum y_{t-1}^2$. `dating_pdc()` applies it in PDC's order (collapse first, then origination on the left sub-sample) and once more on the right sub-sample for KS's recovery date when `regimes = 4`.

**Validation.**

1. `pdc_find_break()` matches a brute-force scan that fits `lm(y[t] ~ y[t-1] - 1)` on every candidate split, bit for bit (tolerance $10^{-8}$).
2. On a synthetic three-regime series (unit root, explosive, collapse) and a four-regime series (with a recovery), `dating_pdc()` recovers the true breaks to within one or two observations in the low-noise, long-series limit. The collapse regime in these tests is a stationary AR(1) ($\rho = 0.5$). A deterministic exponential decay would be a poor test DGP, because its flat tail is indistinguishable from the random walk of the recovery regime.
3. At moderate $T$ ($T = 350$, 30 seeds) the exact-date recovery rate is 3.3% for the origination and 0% for the collapse, with mean absolute errors of 5.5 and 1.0. This matches the Monte Carlo of KS (Section 5), which reports about 30% exact-date recovery at $T = 400$ and about 65% at $T = 800$.

The full suite has 13 assertions for this function, all passing. Replication script: [replication/dating-and-root-inference/radf_pdc_validation.R](#script-radf_pdc_validation).

---

## Root inference

Status: done. `rootstamp()` in `exuber/R/rootstamp.R` is one S3 generic. The `default` method handles a single sub-sample and the `radf_obj` method handles every episode of a `datestamp()` result at once. Tests are in `exuber/tests/testthat/test-rootstamp.R`.

### Sources

- Phillips, P. C. B. & Magdalinos, T. (2007). Limit theory for moderate deviations from a unit root. *Journal of Econometrics*, 136(1), 115–130. Working paper: Cowles Foundation DP 1471, open at `cowles.yale.edu/sites/default/files/2022-08/d1471.pdf`.
- Guo, G., Sun, Y. & Wang, S. (2019). Testing for moderate explosiveness. *The Econometrics Journal*, 22(3), 279–303.
- Skrobotov, A. (2023). Testing for explosive bubbles: a review. arXiv:2207.08249, which restates both results.

### Two results

1. **Phillips & Magdalinos (2007), Cauchy limit.** For a mildly explosive AR(1) with $\rho_n = 1 + c/n^\alpha$, $c > 0$ and $\alpha \in (0,1)$, Theorem 4.3 (their eq. 26) gives

   $$
   \frac{n^\alpha \rho_n^n}{2c}\,(\hat\rho_n - \rho_n) \Rightarrow C,
   $$

   with $C$ standard Cauchy, also under non-Gaussian errors. Remark (i) (eq. 27) gives the simpler fixed-root case of White (1958), which needs neither $\alpha$ nor $c$:

   $$
   \frac{\rho^n}{\rho^2 - 1}\,(\hat\rho_n - \rho) \Rightarrow C .
   $$

   Replacing $\rho$ by $\hat\rho$ in the normalisation gives the two-sided interval

   $$
   \hat\rho \pm q_{\alpha/2}\,\frac{\hat\rho^2 - 1}{\hat\rho^n},
   $$

   where $q_{\alpha/2}$ is a standard Cauchy quantile (`qcauchy()`, or equivalently `qt(., df = 1)`). The interval of eq. 26 needs an estimate of $\alpha$, and is not implemented.
2. **Guo, Sun & Wang (2019), normal limit.** With a drift term allowed, the ordinary regression $t$-statistic for $\rho_T$ is asymptotically standard normal under i.i.d. errors (Student's $t$ or HAR under dependence). It does not require $c$, $k_T$ or the rate, so it is the simpler interval and the default.

### Interface

`rootstamp()` fits the no-intercept AR(1) regression of the Phillips–Magdalinos model on a sub-sample, for example an episode already identified by `datestamp()`. It reports

- `rho` and `rho_ci`, the point estimate and the Wald interval ($\hat\rho \pm z\,\mathrm{se}$, with the standard error from the no-intercept OLS fit), or the Cauchy interval with `type = "cauchy"`;
- `doubling_time` and `doubling_time_ci`, $\log 2 / \log\hat\rho$, the number of periods for the series to double at the estimated rate. The interval comes from transforming the endpoints of the $\rho$ interval, and the bounds flip because the doubling time decreases in $\rho$.

The Cauchy interval assumes a fixed explosive root, while the normal-$t$ interval allows drift and dependence. Root inference on very short episodes (duration 1 or 2) is as unreliable as fitting a regression to two or three points. Use the `min_duration` argument of `datestamp()` to filter them.

Also relevant background: Phillips, Magdalinos & Giraitis (2010), *J. Econometrics* 158(2), 274–279, show that moderate-deviation theory joins the local-to-unity case smoothly as $\alpha \to 0$. It would be the reference for roots close to the local-to-unity boundary.

### Validation

- **Cauchy percentiles.** The standard Cauchy distribution is Student's $t$ with one degree of freedom, so the percentiles of Skrobotov's footnote 17 can be checked exactly: `qt(0.95, 1) = 6.313752`, `qt(0.975, 1) = 12.7062` and `qt(0.995, 1) = 63.65674`, against $6.315$, $12.7$ and $63.65674$ in the footnote. They are tested in `test-rootstamp.R`.
- **Cauchy interval.** It brackets the point estimate and matches the closed-form eq. 27 formula exactly.
- **Point estimate.** On an explosive AR(1) with $\rho = 1.05$ and $n = 200$, `rootstamp()` recovers $\hat\rho = 1.0500$. The confidence interval is indistinguishable from the point estimate at that precision. This is the super-consistency of explosive-root estimation: the regressor $y_{t-1}$ grows geometrically, so $\sum y_{t-1}^2$ grows at rate $\rho^{2n}$ and the standard error collapses much faster than in the unit-root or stationary case.
- **Coverage.** For $\rho = 1.03$ and an episode of 149 observations, 500 replications give about 90% coverage for a nominal 95% interval. For $\rho = 1.05$ and $n = 200$, 800 replications give 94.6%. Undercoverage is a finite-sample effect of a $T \to \infty$ result for an estimator that converges slowly, and it shrinks as $n$ or $\rho - 1$ grows. The package test asserts a loose bound (above 80%) and not the nominal rate.

Replication script: [replication/dating-and-root-inference/rootstamp_validation.R](#script-rootstamp_validation).

---

## Confidence sets for bubble dates

Status: evaluated, not implemented.

### Source

Kurozumi, E. & Skrobotov, A. (2025). Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble. arXiv:2511.16172.

### Idea

The paper builds a confidence interval for estimated dates. It is the dating analogue of `rootstamp()`, and not a new detection method. It does not use the limiting distribution of the break-date estimator, which the authors found performs poorly for bubble dates. Instead, it inverts hypothesis tests on the break location: a likelihood-ratio-type test (Eo & Morley 2015) and Elliott–Müller-type tests (2007), used separately and combined. New limiting null and alternative distributions are derived for each and evaluated by Monte Carlo. The emergence, collapse and recovery dates are estimated separately.

### Ingredients

- The "12"-direction tests $\mathrm{LR}^e_{a,12}$ and $\mathrm{EM}^e_{a,12}$ (eq. 15–16) have closed-form critical values, $cv^e_{\mathrm{LR}12, 0.05} = \lambda_1\, \chi^2_{1, 0.05}$ and $cv^e_{\mathrm{EM}12, 0.05} = \sqrt{\lambda_1\, \chi^2_{1,0.05}}$.
- The "21"-direction tests ($\mathrm{LR}^e_{a,21}$, $\mathrm{EM}^e_{a,21}$, $\mathrm{EM}^e_{b,21}$, eq. 17–19) have no closed form, but the paper publishes a response-surface regression for their critical values, $cv = a_{0,\ell} + a_{-1,\ell}/\lambda_1^* + a_{1,\ell}\lambda_1^* + a_{2,\ell}\lambda_1^{*2} + a_{3,\ell}\lambda_1^{*3}$, with coefficients in its Table 1.
- The recommended "$\mathrm{LE}^e$" test combines $\mathrm{LR}^e_{b,12}$ with $\mathrm{EM}^e_{a,21}$, because $\mathrm{LR}^e_{a,12}$ alone is over-sized in finite samples. $\mathrm{LR}^e_{b,12}$ (eq. 11) is a minimum over candidate break dates of

  $$\frac{y^2_{T_2} - \hat\rho_a \sum_{t=T_1+1}^{T_2} y^2_{t-1}}{T\, \hat\phi_a^{2(T_2 - T_1)}\, \hat\sigma^2 / 2},$$

  whose numerator follows the prefix-sum pattern of `hls_prefix_sums()`. $\mathrm{EM}^e_{a,21}$ (eq. 18) is an integral over a continuum of candidate break points of an $\mathrm{ADF}(\lambda_2^*, \lambda_1^*)$ functional, which has no counterpart in exuber and is not a discrete search.

The construction needs the estimators $\hat\rho_a$, $\hat\phi_a$ and $\hat\sigma^2$ and the admissible-date set $\Lambda^e_{12}$ from Section 2 of the paper, and then the integral of $\mathrm{EM}^e_{a,21}$, repeated for three dates. The work is of the scale of HLS/HLW.

---

## Improved retrospective dating

Status: done, as `dating_knp(breaks = )`, covering the single-bubble correction and the multi-bubble dynamic programme of Section 3.

### Source

Kejriwal, M., Nguyen, L. & Perron, P. (2025). An Improved Procedure for Retrospectively Dating the Emergence and Collapse of Bubbles. *Journal of Time Series Analysis*, 46(5), 867–883, `doi:10.1111/jtsa.12810`.

### Idea

KNP fix a bias in the joint-SSR estimator of the HLS family. Their model uses HLS's fixed autoregressive coefficient (not the mildly explosive $\rho_T \to 1$ of Phillips–Magdalinos) with an abrupt collapse. Theorem 1 shows that the standard joint-SSR estimator is inconsistent. The origination estimate converges to the collapse date, and the collapse estimate converges to a date after the true collapse, offset by the trimming parameter. The fix (Theorem 2) is a modified SSR that omits the single residual at the implosion date, which restores consistency of both dates. A footnote shows that the omission is numerically equivalent to a one-time dummy in HLS's Model 4 regression.

Section 3 adds a Bai–Perron/Perron–Qu-style dynamic programme for several bubbles. The unit-root regimes are restricted ($\mu = 0$, $\rho = 1$), so each segment cost is known in closed form and the iteration over initial values that Perron–Qu need is not required.

KNP's single-bubble model has the shape of HLS's Model 2: an unfitted unit root, an explosive regime fitted with intercept and slope, and an unfitted unit root after an instantaneous collapse. Regressing the level $y_t$ on $y_{t-1}$ with an intercept gives the same residuals and SSR as regressing $\Delta y_t$ on $y_{t-1}$ (the slope shifts by one), which is the regression that `hls_segment_ssr()` computes. The whole correction is

$$
\mathrm{SSR}_{\mathrm{om}}(T_1, T_2) = \mathrm{SSR}(T_1, T_2) - (\Delta y_{T_2 + 1})^2,
$$

an already computed SSR minus one squared term.

### Implementation

`dating_knp(data, trim = 0.05, omit = TRUE, breaks = 2L)` is in `exuber/R/dating_knp.R` and is tested in `exuber/tests/testthat/test-knp.R`. `knp_find_break()` reuses `hls_prefix_sums()` and `hls_segment_ssr()` and searches $(\tau_1, \tau_2)$ jointly to minimise the omission-corrected SSR. With `omit = FALSE` it minimises the plain SSR, which is inconsistent, so that the effect of the correction can be shown. Unlike `hls_model23()`, the candidate set has no sign constraint on the fitted peak.

`breaks` is the paper's $m$: two per bubble, or an odd number to let the last bubble run to the sample end (its collapse is `NA`). As in the paper, $m$ is taken as given, because KNP leave its selection open. `breaks = 2` keeps the exhaustive single-bubble search. For more breaks `knp_dp()` runs the dynamic programme with the objective of their eq. 11. Regimes alternate between unit root and explosive, starting with a unit root, and every unit-root regime after the first omits its first residual. Each segment cost is $\sum z^2$ over the segment (minus its first term after a collapse when `omit = TRUE`) or the intercept-and-slope OLS SSR (`hls_segment_ssr(..., fit = TRUE)`). The programme is $O(mT^2)$ with $O(1)$ segment costs and returns the exact global minimiser of the grid search. With more breaks, `origination`, `collapse` and `delta` are matrices with one row per bubble.

### Validation

The formula check is exact. `omit = FALSE` and `omit = TRUE` both match an exhaustive nested-`lm()` search, and `knp_dp(y, 2)` returns the same dates and SSR as the single-bubble search ($|\Delta \mathrm{SSR}| = 0$). With three and four breaks it matches a brute-force search over every admissible partition, with and without omission ($n = 28$, $|\Delta\mathrm{SSR}| = 2.7 \times 10^{-14}$).

The Monte Carlo uses KNP's DGP (unit root, no-intercept explosive AR(1), instantaneous collapse back near the pre-bubble level, fresh unit root).

- **Single bubble** ($T = 200$, $T_1 = 50$, $T_2 = 90$, $\delta = 1.05$, 30 replications). Theorem 1 is reproduced for the naive estimator: the mean error of the origination estimate against the true collapse date is $\operatorname{mean}|\hat\tau_1 - T_2| = 1.0$, far below its error against its own origination, $\operatorname{mean}|\hat\tau_1 - T_1| = 39.0$. With the omission the latter falls from 39.0 to 13.0 observations. The residual bias at $T = 200$ is expected, because Theorem 2 is an asymptotic result. The collapse date and the explosive coefficient are close to their true values ($\operatorname{mean}|\hat\tau_2 - T_2| = 1.0$, mean $\hat\delta = 0.986$ against a true $1.05$).
- **Two bubbles** ($T = 200$, bubbles over 41–70 and 121–150, $\delta = 1.05$, 50 replications). The mean absolute date error per break (origination 1, collapse 1, origination 2, collapse 2) is $11.9 / 5.9 / 10.3 / 3.8$ observations with the omission correction and $33.4 / 18.5 / 29.6 / 13.2$ without it. The inconsistency of Theorem 1 carries over to several bubbles, and so does the fix.

`dating_knp()` is also in pyexuber (`_knp_dp()`), where it agrees with R and with the same brute force. Replication script: [replication/dating-and-root-inference/radf_knp_validation.R](#script-radf_knp_validation).

---

## WLS dating under time-varying volatility

Status: done, as `dating_pdc(..., type = "wls")`.

### Source

Kurozumi, E. & Skrobotov, A. (2023). Improving the accuracy of bubble date estimators under time-varying volatility. arXiv:2306.02977.

### Idea

The estimator is a two-step generalisation of the PDC/KS sequential estimator. Step 1 is `dating_pdc()` as described above: fit the homoskedastic no-intercept AR(1) break model and keep the residuals. Step 2 estimates the time-varying error variance $\sigma_t^2$ nonparametrically from those residuals and re-estimates each break by minimising a weighted SSR,

$$
\sum_t \frac{(y_t - a\, y_{t-1})^2}{\sigma_t^2},
$$

in place of the unweighted sum. This is the cumulative-sum construction of `pdc_find_break()` with each $y_t$, $y_t y_{t-1}$ and $y_{t-1}^2$ divided by $\sigma_t^2$ before the prefix sum, so it remains closed form and $O(T)$ per break.

### Implementation

`dating_pdc(data, ..., type = c("ols", "wls"))` is in `exuber/R/dating_pdc.R`. `type = "ols"` is the original estimator. `type = "wls"` does the following.

1. Run the sequential `type = "ols"` fit to get the first-step breaks.
2. `pdc_regime_resid(y, breaks)` computes the fitted no-intercept AR(1) residual at every pair $(y_{t-1}, y_t)$, with one OLS $\rho$ per regime implied by those breaks.
3. `nw_spot_vol()`, a Nadaraya–Watson kernel smoother with a leave-one-out cross-validated bandwidth, turns the residuals into $\hat\sigma_t^2$. It is the smoother shared with SBZ, which `kernel_spot_vol(y)` calls as `nw_spot_vol(diff(y))`.
4. `pdc_find_break()` takes an optional `weights` argument (`NULL` gives the unweighted search). Each cumulative sum is multiplied by the weights before the prefix sum.
5. The sequential search (collapse, then origination, then recovery) runs again with `weights = 1 / sigma_t^2`, using the matching slice of the full-sample variance vector for each sub-sample.

No critical values are needed, because this is point estimation and not a threshold-crossing test.

### Validation

Two Monte Carlo checks with 40 seeds each.

- **Homoskedastic DGP.** The origination-date mean absolute error is 5.42 for OLS and 5.53 for WLS. Weighting costs essentially nothing when there is no heteroskedasticity to exploit.
- **Volatility burst in the first 20% of the pre-bubble regime**, the scenario with the largest gains in the paper's Monte Carlo. The origination-date mean absolute error falls from 13.05 (OLS) to 2.33 (WLS), about 5.6 times lower. The unweighted objective lets the noisy early segment dominate the origination split, and WLS downweights it. The collapse date is near-exact in both cases, because the explosive-to-collapse transition dominates the SSR whatever happens earlier, as PDC's stochastic-order argument says.

A test in `test-pdc.R` asserts a loose version of this margin (a factor of 2).

Replication scripts: [replication/dating-and-root-inference/radf_pdc_wls_heteroskedastic_mae.R](#script-radf_pdc_wls_heteroskedastic_mae), [radf_pdc_wls_homoskedastic_mae.R](#script-radf_pdc_wls_homoskedastic_mae).

---

## Reverse-regression recovery dating

Status: done, as `radf_recovery()` and `radf_recovery_cv()`, with the caveats listed under Validation. The recovery date $f_r$ behaves well. The crisis-origination date $f_c$ and the false-detection rate under the null are noisier.

### Source

Phillips, P. C. B. & Shi, S. (2014). Financial Bubble Implosion. Cowles Foundation DP 1967, published as Financial Bubble Implosion and Reverse Regression, *Econometric Theory*.

### Idea

Reverse the series, $X^*_t = X_{T+1-t}$, run the same BSDF/BSADF recursion that PSY uses on $X^*$, and map the crossing fractions back to the original time index (their eq. 8–9):

$$
\begin{aligned}
\hat f_r &= 1 - \hat g_e, & \hat g_e &= \inf\{ g \in [g_0, 1] : \mathrm{BSDF}_g(g_0) > \mathrm{scv} \} && \text{(recovery date)},\\
\hat f_c &= 1 - \hat g_c, & \hat g_c &= \inf\{ g \in [\hat g_e, 1] : \mathrm{BSDF}_g(g_0) < \mathrm{scv} \} && \text{(crisis-origination date)}.
\end{aligned}
$$

Because $\hat g_c$ is searched only after $\hat g_e$, $\hat f_c \le \hat f_r$ always. $\hat f_c$ is the collapse-onset date of the original series, derived by reverse regression as an alternative to the collapse date that PSY's forward test already gives, and $\hat f_r$ follows it. Reversing a mildly explosive process followed by a mildly integrated collapse turns the collapse regime into an explosive regime in reverse time, so recovery and crisis dating become the right-tailed test of PSY applied to `rev(x)`.

Theorem 1 shows that the null limit of the reverse statistic, $F_g(W, g_0)$, is not the forward distribution $F_f(W, f_0)$. Reversing a random walk makes the reversed lagged regressor correlated with the reversed current error ($E[X^*_{T-j+2}\,\varepsilon_{T-j+2}] \ne 0$), which changes the critical values even under the null. A paired Monte Carlo ($n = 100$, `minw = 20`, 5000 replications) comparing `radf_mc_cv()` with the same recursion on the reversed path gives critical values that differ by about 0.04 on average and up to about 0.11 at the 95% level.

Section 4.3 of the paper has a separate sequential extension (eq. 10–11). It applies the reverse regression repeatedly on a growing sample from the collapse date $T_c$ forward, and stops at the first sample end for which a correction is detected. It produces the further correction of January 2004 and the full return to normal conditions in May 2004 in their dot-com application. It is not part of the eq. 8–9 pair and is not implemented. It would be an outer loop around `radf_recovery()`, similar in structure to `monitor()`.

In that application (Section 5, NASDAQ price–dividend ratio, reported in prose) the eq. 8–9 pair gives a crash from March to November 2000, so $f_c$ is March 2000 and $f_r$ is November 2000. We have not reproduced this, because the underlying series is not available.

### Implementation

`radf_recovery()` runs the `bsadf` recursion of `radf()` on the reversed series, compares it with a reversal-calibrated boundary and maps the first up-crossing and down-crossing back through $f = n + 1 - g$. `radf_recovery_cv()` produces the boundary by the simulate-then-quantile construction of `radf_mc_cv()`, including its `cummax(badf)` shortcut for the `bsadf` boundary, with one added `rev()` before the recursion. No C++ and no new statistic are needed.

### Validation and caveats

- $f_c \le f_r$ holds by construction whenever both dates are identified and uncensored, and held in every replication.
- The bias of $f_r$ on synthetic collapse-then-recovery data is a few observations, in the same direction and of similar size as the roughly six observations early in Table 5 of the paper.
- The bias of $f_c$ is larger, approaching the length of the synthetic collapse window in some runs.
- Under a pure random-walk null ($n = 100$, `minw = 20`, 95% level, one reversal-calibrated critical value reused across 200 fresh draws) the false-detection rate is about 29%. This is higher than comparable forward-test numbers in this project, such as the cumulative false-alarm rate of about 10% for `monitor()` over a 75-point horizon. The `inf` in eq. 9 defines the first down-crossing with no persistence requirement, so a transient noise-driven dip below the boundary triggers a premature $f_c$. This is a property of the literal crossing rule under finite-sample noise.

The roxygen documentation of `radf_recovery()` carries the same caveat. Replication script: [replication/dating-and-root-inference/radf_recovery_validation.R](#script-radf_recovery_validation).
