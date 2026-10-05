---
title: "Real-time monitoring for bubbles"
blurb: "Sequential and real-time detection: training-vs-monitoring orchestration, CUSUM families, and closed-form boundaries."
order: 3
---
Monitoring asks a different question from the retrospective tests. Instead of testing a complete sample, it fixes a training window that is assumed to be free of bubbles and then checks each new observation as it arrives, so that the first alarm comes as soon as the series turns explosive. The methods differ in the detector they use and in how they control the false-alarm rate. There are two families. Family A is the recursive training-maximum detector of Phillips & Shi (2020), which reuses the BSADF sequence of `radf()`. Family B is the CUSUM family of Homm & Breitung (2012) and its relatives.

Status labels are the ones used in [volatility-robustness.md](/replication/volatility-robustness).

| Method | Paper | Status |
|---|---|---|
| [Recursive monitoring, Family A: `monitor()`](#recursive-monitoring-monitor) | Phillips & Shi (2020) | done |
| [CUSUM and CUSUMV: `monitor_cusum()`](#cusum-monitor_cusum) | Homm & Breitung (2012); Astill et al. (2023) | done, including the finite-sample boundary of Homm & Breitung |
| [FLUC statistic: `monitor(boundary = "fluc")`](#fluc-monitorboundary--fluc) | Homm & Breitung (2012) | done |
| [Closed-form SADF and GSADF boundaries](#kurozumi-2020-2021-sadf-and-gsadf-boundaries) | Kurozumi (2020) | done: `monitor(boundary = "kurozumi", s0 = 0/0.4/0.8)` |
| [LBI test and sequential monitoring](#breitung--diegel-2025-lbi-test-and-sequential-extension) | Breitung & Diegel (2025) | done: `lbi_test()`, `monitor_lbi()` with the constant boundary `mCUSUM`/`wCUSUM` |
| [Robust Chebyshev-type monitoring (RCA)](#horváth--trapani-20232026-rca-monitoring) | Horváth & Trapani (2023/2026) | evaluated, not implemented |
| Delay-time paper | Kurozumi (2021) | evaluated, not implemented |

All papers are listed in [references.md](/replication/references#monitoring).

## Sources

1. Homm, U. & Breitung, J. (2012). Testing for speculative bubbles in stock markets: a comparison of alternative methods. *Journal of Financial Econometrics*, 10(1), 198–231, `doi:10.1093/jjfinec/nbr009`.
2. Astill, S., Harvey, D. I., Leybourne, S. J., Taylor, A. M. R. & Zu, Y. (2023). CUSUM-Based Monitoring for Explosive Episodes in Financial Data in the Presence of Time-Varying Volatility. *Journal of Financial Econometrics*, 21(1), 187–227, `doi:10.1093/jjfinec/nbab009` ("AHLTZ"). The earlier companion, Astill, Harvey, Leybourne, Sollis & Taylor (2018), Real-Time Monitoring for Explosive Financial Bubbles, *JTSA*, 39, 863–891 ("AHLST"), compares each monitoring statistic with the maximum of the training sample and not with a CUSUM boundary.
3. Whitehouse, E. J., Harvey, D. I. & Leybourne, S. J. (2025). Real-time monitoring procedures for early detection of bubbles. *International Journal of Forecasting*, 41(3), 1260–1277, `doi:10.1016/j.ijforecast.2024.12.005`. Open access. It supplies the AHLST decision rule and false-positive-rate formula below.
4. Kurozumi, E. (2020). Asymptotic properties of bubble monitoring tests. *Econometric Reviews*, 39(5), 510–538, `doi:10.1080/07474938.2019.1697086`. It extends SADF and GSADF to a monitoring scheme, studies a CUSUM detector next to them and derives monitoring-period critical values under moderate-deviation and local-to-unity asymptotics. Kurozumi, E. (2021). Asymptotic Behavior of Delay Times of Bubble Monitoring Tests. *JTSA*, 42(3), 314–337, `doi:10.1111/jtsa.12569`, concerns the stochastic order of the detection delay.
5. Horváth, L. & Trapani, L. (2026). Real-time monitoring with RCA models. *Econometric Theory*, 42, 514–547. Working paper arXiv:2312.11710.
6. Astill, S., Taylor, A. M. R. & Zu, Y. (2026, forthcoming). Covariate Augmented CUSUM Bubble Monitoring Procedures. *Econometric Theory*. Essex Finance Centre Working Paper No. 94. Section 3 restates the CUSUM statistic and boundary of Homm & Breitung and the volatility-robust modification of Astill et al., with equation numbers and page references.
7. Breitung, J. & Diegel, M. (2025). Sequential Detector Statistics for Speculative Bubbles. *JTSA*, 46(5).

## The two families

**Family A** compares a training maximum with the monitoring statistic (AHLST, Whitehouse et al., Phillips & Shi). The sample is split into a training period $t = 1, \dots, T^*$, assumed free of bubbles, and a monitoring period $t = T^*+1, \dots, T$. A recursive statistic $A_{e,k}$ is computed over sub-samples of fixed length $k$ in both periods. The training maximum $A^*_{\max} = \max A_{e,k}$ becomes a fixed critical value, and monitoring rejects $H_0$ at the first $e$ with $A_{e,k} > A^*_{\max}$. AHLST prove a closed-form asymptotic false-positive rate (FPR) that depends only on the ratio of training to monitoring length (eq. 6 below). Phillips & Shi use a bootstrap in place of the closed form.

**Family B** uses CUSUM and Page-CUSUM detectors (Homm & Breitung, Astill et al., Kurozumi's CUSUM variant, Horváth & Trapani, Breitung & Diegel). A partial sum of standardised first differences $\Delta y_t$ is accumulated from the end of the training sample and compared with a boundary that grows with $t$, for example $c_t \sqrt t$. The boundary is chosen so that the cumulative false-alarm probability over the whole, possibly infinite, monitoring horizon stays below $\alpha$. This is a running standardised sum, not a recursive ADF regression. The volatility-robust variants replace the standardisation by a kernel spot-variance estimate, which resembles the kernel machinery of [SBZ](/replication/volatility-robustness#sbz-wls--kernel-volatility). Horváth & Trapani extend the idea to random-coefficient autoregressions (RCA) with weighted CUSUM and Page-CUSUM detectors that work for transitions in both directions between stationary and explosive regimes.

Kurozumi (2020) places SADF/GSADF-type (Family A) and CUSUM-type (Family B) monitoring in one asymptotic framework, derives monitoring-period critical values for both and studies the detection delay separately. CUSUM detects an early, short bubble faster, and ADF/BSADF-type detectors detect a middle-to-late bubble faster. A union of rejections that combines BSADF and CUSUM is possible, and is the monitoring counterpart of the [SBZ union statistic](/replication/volatility-robustness#sbz-wls--kernel-volatility).

## Recursive monitoring: `monitor()`

Status: done. `monitor(data, r_star = 0.5, minw, nboot, level, adflag, type, seed)` is in `exuber/R/monitor.R`.

Two facts make it a thin layer over existing code.

1. `radf_wb_ps_cv(..., tb = T*)` computes the training-window wild-bootstrap critical value that monitoring needs and broadcasts it as a constant boundary across the monitoring horizon.
2. The BSADF statistic of `radf()` at calendar time $t$ depends only on data up to $t$ and equals the last BSADF value of a fresh `radf(y[1:t])`. The whole monitoring path therefore comes from one full-sample `radf()` call, which is $O(T)$.

To avoid look-ahead, `monitor()` calls `radf_wb_ps_cv()` on `data[1:T*]` only. Inside `radf_wb_ps_cv()` the null-model fit (`adf_res()` in `radf_wb.R`) uses all the data it receives to estimate the bootstrap residuals and coefficients, and `tb` only truncates the simulated bootstrap sample. Passing the full series would let post-$T^*$, possibly explosive, data influence the training-window calibration.

**Validation.**

1. A basic run returns a well-formed object, with `T_star` and a `bsadf` row count of `n - minw`. An alarm always falls strictly after $T^*$ (checked over 10 seeds with a post-training bubble).
2. Under $H_0$ (no bubble, 40 replications, 75-observation monitoring horizon, 95% per-point threshold) the false-alarm rate is 10%. It exceeds the 5% per-point level because it is a cumulative probability over 75 sequential comparisons against a fixed boundary, and it grows with the horizon. The AHLST formula (eq. 6 below) documents the same property.
3. For a bubble that starts after $T^*$ (15 replications) the detection rate is 86.7%. The alarm delay, the alarm date minus the true origination date, is always positive: minimum 6, median 19 and maximum 33 observations.

Tests are in `test-monitor.R`. Replication script: [replication/monitoring/radf_monitor_validation.R](#script-radf_monitor_validation).

Not implemented for Family A: a closed-form or simulated false-alarm-versus-horizon boundary function (AHLST eq. 4–6), a union of rejections across families and date-stamping methods specific to a monitoring result.

## CUSUM: `monitor_cusum()`

Status: done, as `monitor_cusum(data, r_star, b_alpha, type, boundary)`.

The CUSUM procedure of Homm & Breitung (Section 3, eq. 26–30) is a standardised running sum of first differences compared with a closed-form asymptotic boundary derived from an inequality of Chu, Stinchcombe & White (1996). It needs no bootstrap and no simulation. For a training window ending at $T^*$ and a monitoring point $t > T^*$,

$$
S_t = \frac{y_t - y_{T^*}}{\hat\sigma_t}, \qquad \hat\sigma_t^2 = \frac{1}{t-1} \sum_{j=2}^{t} (\Delta y_j)^2,
$$

where the numerator telescopes the post-training first differences and $\hat\sigma_t^2$ is the recursive variance of all differences up to $t$. The boundary is

$$
c_t = \sqrt{b_\alpha + \log(t / T^*)}, \qquad \text{boundary}_t = c_t \sqrt t ,
$$

and the alarm is the first $t$ with $S_t > \text{boundary}_t$. The constant $b_\alpha = 4.6$ is the one-sided asymptotic calibration for a 5% level. Because $\hat\sigma_t$ uses only data up to the current monitoring point, there is no look-ahead.

Homm & Breitung propose two statistics. CUSUM is `monitor_cusum()`. FLUC is `monitor(..., boundary = "fluc")`, described below. The code is an internal `cusum_stat_path()` and the exported `monitor_cusum()`, which returns a `radf_cusum_obj`. It shares no code with `monitor()`, `radf()` or `exubercore`.

### Finite-sample boundary

`boundary = "finite"` replaces $b_\alpha = 4.6$ with the finite-sample constant of Homm & Breitung's Table 8 ("without drift estimation", which matches the raw-first-difference construction here). It is indexed by training length, significance level and the horizon ratio $k = N/T^*$. The values are transcribed from the table and no simulation is run. The option applies to both `type = "standard"` and `type = "kernel"`, because Corollary 1 of Astill et al. shows that the same boundary function serves both statistics.

### CUSUMV: the volatility-robust variant

Status: done, as `monitor_cusum(..., type = "kernel")`.

Astill et al. (2023) allow for time-varying volatility, which can "heavily inflate the false positive rate (FPR) of the CUSUM-based procedure". Their eq. 6–7 standardise each first difference individually by a one-sided Nadaraya–Watson estimate of the spot variance before cumulating. One-sided means that it uses only current and past lags, as real-time monitoring requires:

$$
SV_t = \sum_{j=T^*+1}^{t} \frac{\Delta y_j}{\hat\sigma_{j,N}}, \qquad \hat\sigma^2_{j,N} = \sum_{s=0}^{N} w_s\, (\Delta y_{j-s})^2, \qquad w_s = \frac{K(s/N)}{\sum_{s=0}^{N} K(s/N)} .
$$

Corollary 1 shows that the same boundary $c_t \sqrt t$ still controls the asymptotic false-alarm rate under time-varying volatility.

`one_sided_kernel_spot_vol()` builds a fixed causal kernel weight vector with `stats::filter(..., sides = 1)`, with $\hat\sigma^2_{j,N} := 1$ for $j \le N$ as in the paper, and `cusum_stat_path_kernel()` forms the statistic. The default bandwidth is $N = 20$, the value that AHLTZ recommend ("setting H = 20 delivered a procedure with the best trade-off" between FPR robustness and power). The data-driven cross-validated bandwidth of their eq. 8–9 is not implemented.

### Validation

1. **Formulas.** `cusum_stat_path()`, `one_sided_kernel_spot_vol()` and `cusum_stat_path_kernel()` match brute-force loop recomputations to floating-point precision.
2. **False alarms under a homoskedastic $H_0$** (pure random walk, 100 replications, 75-observation horizon). The asymptotic CUSUM rate is 0%, as expected of a conservative bound (eq. 28, via Chu et al.). `standard` and `kernel` agree, as in Remark 8 of AHLTZ.
3. **False alarms under a heteroskedastic $H_0$** (volatility jumping from 1 to 8 during the monitoring region, 60 replications). `standard` CUSUM rises to 8.3%, and `kernel` CUSUMV stays at 0%, which is the failure mode that AHLTZ describe.
4. **Power** under the post-training bubble DGP used for `monitor()` (30 replications): 30.0% for `standard` and 36.7% for `kernel`, with a median alarm delay of 27 observations for `standard`. This is well below the 86.7% detection rate and 19-observation median delay of `monitor()` on the same DGP. The DGP starts the bubble about 65% of the way into the sample, the middle-to-late regime in which Kurozumi (2020, 2021) finds CUSUM-type detectors lag ADF-type ones. It illustrates why a union of both families is recommended.
5. **Finite-sample boundary** at $T^* = 75$, $n = 150$ ($k = 2$, with $T^*$ snapped to the tabulated $n = 50$): the false-alarm rate is 9.0% under $H_0$, closer to the nominal 5% than the 0% of the asymptotic bound, and detection is 36.7% against 30.0%.

Tests are in `test-cusum.R`. Replication scripts: [replication/monitoring/radf_cusum_validation.R](#script-radf_cusum_validation), [radf_cusum_finite_boundary_validation.R](#script-radf_cusum_finite_boundary_validation), [radf_cusumv_kernel_validation.R](#script-radf_cusumv_kernel_validation).

## FLUC: `monitor(boundary = "fluc")`

Status: done.

Homm & Breitung's second statistic (eq. 27) is $Z_t = (\hat\rho_t - 1)/\hat\sigma_{\hat\rho_t} = \mathrm{DF}_{t/n}$, the ordinary expanding-window OLS ADF $t$-statistic on $\{y_0, \dots, y_t\}$. This is the `badf` sequence of `radf()`, so no new statistic is needed. The rejection rule (eq. 29/31) is $\mathrm{DF}_{t/n} > \kappa_t$ with $\kappa_t = \sqrt{b_{k,\alpha} + \log(t/n)}$. The constant $b_{k,\alpha}$ comes from simulation in the paper, and the paper publishes it (Table 7, part i, "without detrending", which matches the no-trend default of `radf()`). It is tabulated by training length $n \in \{20, 50, 100\}$, level $\alpha \in \{0.10, 0.05, 0.01\}$ and horizon ratio $k = N/n \in \{2, 3, 4, 5, 6, 8, 10\}$, where $N$ is the total sample including training.

`hb_fluc_table` is the $9 \times 7$ table and `hb_fluc_q(level, n_train, k)` snaps `n_train` to the nearest of $\{20, 50, 100\}$ and $k$ to the nearest of $\{2, \dots, 10\}$, and requires `level` to be one of the three tabulated values. It follows the pattern of `kurozumi_sadf_q()` and is the third `boundary` option of `monitor()`, next to `"bootstrap"` and `"kurozumi"`. The table covers training lengths only up to $n = 100$, so larger $T^*$ are snapped to the $n = 100$ row.

**Validation.** The lookups match Table 7 exactly (6 cells, including a tie-breaking snap). Alarms never fire before $T^*$ (10 of 10 replications). Under $H_0$ ($n = 150$, $T^* = 75$, $k = 2$, the smallest and most conservative horizon ratio, 100 replications) the false-alarm rate is 0%. Detection on the post-training bubble DGP (30 replications) is 56.7%, below the 80% of `boundary = "kurozumi"` and the 90% of `"bootstrap"`. This agrees with Homm & Breitung's finding that FLUC and CUSUM generally have less power than a supDF-style test, although FLUC beats their CUSUM.

Tests extend `test-monitor.R`. Replication script: [replication/monitoring/radf_monitor_fluc_boundary_validation.R](#script-radf_monitor_fluc_boundary_validation).

## Kurozumi (2020, 2021): SADF and GSADF boundaries

Status: done, as `monitor(..., boundary = "kurozumi", s0 = ...)` for $\mathrm{SADF}$ ($s_0 = 0$) and $\mathrm{GSADF}_{s_0}$ ($s_0 = 0.4$ or $0.8$).

Kurozumi's $\mathrm{SADF}(k) := \mathrm{ADF}_1^{m+k}$ is the `badf` sequence of `radf()`, which was confirmed bit for bit against a from-scratch OLS ADF $t$-statistic at three check points (tolerance $10^{-8}$). His

$$
\mathrm{GSADF}_{s_0}(k) := \max_{1 \le k_1 \le \lfloor m s_0 \rfloor} \mathrm{ADF}_{k_1}^{m+k}
$$

differs from `radf()$bsadf`. In `bsadf` the range of the window start grows with the current point $t$ (from 1 to $t - \mathrm{minw}$). Here it is capped at a fixed fraction of the training length $m$ whatever $t$ is, so it needs only a bounded band of start points and no recursion. In both cases a published, table-based threshold replaces the wild-bootstrap boundary of `monitor()`.

### Boundary functions and Table 1

The boundary functions are

$$
\begin{aligned}
\mathrm{SADF}: &\quad g_0^{df}(k/m) = q_0^{df} \quad \text{(constant)},\\
\mathrm{GSADF}: &\quad g_{s_0}^{df}(k/m) = q_{s_0}^{df}\,\bigl(a_{s_0} + b_{s_0} \log(c_{s_0} + k/m)\bigr),\\
&\quad \{a, b, c\} = \{0.76, 0.02, 0.34\} \text{ for } s_0 = 0.4, \quad \{0.73, 0.03, 0.90\} \text{ for } s_0 = 0.8,\\
\mathrm{CS}: &\quad g_\gamma^{cs}(k/m) = q_\gamma^{cs}\,(1 + k/m)^{1-\gamma}\,(k/m)^{\gamma}.
\end{aligned}
$$

Table 1 gives the scaling constants $q$ by significance level $\beta$ and monitoring-horizon ratio $\bar s = \bar k / m$, where monitoring runs $\bar k$ observations past the training length $m$. It covers only $\bar s \in \{1, 3, 5\}$.

| $\bar s$ | $\beta$ | $q_0^{df}$ | $q_{0.4}^{df}$ | $q_{0.8}^{df}$ | $q_{0.25}^{cs}$ | $q_{0.45}^{cs}$ |
|---|---|---|---|---|---|---|
| 1 | 0.10 | 0.6946 | 1.3969 | 1.9369 | 1.5071 | 2.1300 |
| 1 | 0.05 | 1.0381 | 1.8081 | 2.3330 | 1.7646 | 2.3948 |
| 1 | 0.01 | 1.6474 | 2.5927 | 3.0941 | 2.2405 | 2.9265 |
| 3 | 0.10 | 1.0299 | 1.7088 | 2.1315 | 1.6772 | 2.1958 |
| 3 | 0.05 | 1.3330 | 2.0737 | 2.4944 | 1.9619 | 2.4638 |
| 3 | 0.01 | 1.8978 | 2.7677 | 3.2136 | 2.4955 | 3.0163 |
| 5 | 0.10 | 1.1308 | 1.7988 | 2.1794 | 1.7326 | 2.2057 |
| 5 | 0.05 | 1.4255 | 2.1480 | 2.5369 | 2.0182 | 2.4844 |
| 5 | 0.01 | 1.9735 | 2.8276 | 3.2616 | 2.5884 | 3.0476 |

$\mathrm{CS}(k)$ is the CUSUM statistic of Homm & Breitung written in Kurozumi's notation. The $q^{cs}$ columns are finite-sample simulated alternatives to the asymptotic $b_\alpha = 4.6$ of `monitor_cusum()` for $\gamma \in \{0.25, 0.45\}$. Kurozumi obtained them by simulation (50,000 replications, with Brownian motion approximated by normalised i.i.d. sums over increments of $1/1000$).

### SADF case ($s_0 = 0$)

`monitor()` has a `boundary = c("bootstrap", "kurozumi")` argument. With `"kurozumi"`, `kurozumi_sadf_q(level, s_bar)` looks up $q_0^{df}$, snapping $\bar s = (n - T^*)/T^*$ to the nearest of $\{1, 3, 5\}$ and requiring `level` in $\{0.90, 0.95, 0.99\}$. The value is compared with `radf()$badf` over the monitoring window. There is no bootstrap, and `nboot`, `type`, `adflag` and `seed` are ignored. The `stat` field of the returned list holds `bsadf` for `boundary = "bootstrap"` and `badf` for `boundary = "kurozumi"`.

**Validation.** The lookups match Table 1 exactly (6 values, `s_bar` snapping and the error for an invalid `level`). With $n = 150$, $T^* = 75$, $\bar s = 1$ and 100 replications, the false-alarm rate is 4.0% against a nominal 5%, while the wild-bootstrap boundary gives 7.0%. Detection on a post-training bubble (30 replications) is 80% for `kurozumi` and 90% for `bootstrap`. The two boundaries calibrate different statistics (`badf` and `bsadf`). Alarms never fire before $T^*$ (10 of 10 replications). Replication script: [replication/monitoring/radf_monitor_kurozumi_boundary_validation.R](#script-radf_monitor_kurozumi_boundary_validation).

### GSADF case ($s_0 = 0.4, 0.8$)

Because $\lfloor m s_0 \rfloor$ is a small fixed cap, $\mathrm{GSADF}_{s_0}(k)$ needs $\mathrm{ADF}_{k_1}^{t}$ only for $k_1$ in a bounded band. Each is a with-intercept OLS ADF $t$-statistic on a fixed window, computed from cumulative-sum differences as in `hls_segment_ssr()`. `kurozumi_gsadf_stat()` in `exuber/R/monitor.R` computes the band with `outer()`-vectorised differences, with an intercept, unlike the no-intercept `gls_dfstat_grid()`. `s0 = 0.4` or `0.8` (the only values with tabulated $a$, $b$, $c$) switches from the flat SADF boundary to the $k$-varying $g_{s_0}^{df}$ above, with $q_{s_0}^{df}$ from the `q04_df` or `q08_df` column. The default `s0 = 0` reproduces the SADF behaviour exactly.

**Validation.** `kurozumi_gsadf_stat()` matches `radf()$badf` to machine precision at `k1_max = 1`, and matches a brute-force `lm()` search over the restricted start band (`|diff| < 1e-14`) at three monitoring points on a 150-observation series. The `q04_df` and `q08_df` lookups are exact, with the expected tie-breaking snap ($s_0 = 0.6$ goes to the lower value, $0.4$). Alarms never fire before $T^*$ (30 of 30 replications). Under $H_0$ (300 replications, $n = 150$, $T^* = 75$) the false-alarm rate is 4.3%, 5.3% and 4.3% for SADF, $\mathrm{GSADF}_{0.4}$ and $\mathrm{GSADF}_{0.8}$, against a nominal 5%. Detection on a post-training bubble (60 replications) is 70.0%, 73.3% and 66.7%. The modest edge of $s_0 = 0.4$ over SADF echoes Kurozumi's finding that GSADF works better than SADF in many cases. The weaker result for $s_0 = 0.8$ is specific to this DGP, since a wider start range dilutes power against some alternatives. Nine tests extend `test-monitor.R`. Replication script: [replication/monitoring/radf_monitor_gsadf_s0_validation.R](#script-radf_monitor_gsadf_s0_validation).

### Kurozumi (2021)

Kurozumi (2021) studies the stochastic order of the detection delay for the same detector families. It supports the split cited above: early or short bubbles favour CUSUM, and middle or late bubbles favour ADF. `monitor_cusum()` and `monitor()` reproduce that split on the post-training bubble DGP. The paper adds dating and inference on top of detection and is not a new detector, so it is not implemented.

## Breitung & Diegel (2025): LBI test and sequential extension

Status: done. `lbi_test()` is the static test for a bubble window that spans the full sample, and `monitor_lbi()` is the sequential extension.

### Static LBI test

The paper proposes a locally best invariant (LBI) statistic. It is robust to heteroskedasticity by construction, through the invariance result of Cavaliere (2005), so it needs no wild bootstrap, and its limiting null distribution is standard normal. For a bubble that spans the whole sample ($y_0 = 0$), eq. 4 gives the telescoping identity

$$
2 \sum \Delta y_t\, y_{t-1} = y_T^2 - T \tilde\sigma^2, \qquad \tilde\sigma^2 = \frac1T \sum \Delta y_t^2,
$$

and substituting it into the numerator of the naive DF-type statistic gives eq. 5, $\mathrm{LBI}_T^2 = y_T^2/(\tilde\sigma^2 T)$. The statistic is the standardised sample endpoint:

$$
\mathrm{LBI}_T = \frac{y_T - y_1}{\tilde\sigma \sqrt{T-1}} .
$$

It is compared with a standard normal quantile (for example $1.645$ at 5%). The test is one-sided, because the paper targets positive bubbles only, on the grounds that negative bubbles are economically implausible for a risky asset. It needs no regression, no recursion, no table and no boundary function.

`lbi_test(data, level = 0.95)` is in `exuber/R/lbi_test.R` and is tested in `exuber/tests/testthat/test-lbi.R`. The telescoping identity of eq. 4 holds exactly on a simulated random walk. In a Monte Carlo under $H_0$ (500 replications) the mean and standard deviation of the statistic are $0.023$ and $0.976$ (theory: 0 and 1), and the false-alarm rate at the 95% level is $0.050$. A Kolmogorov–Smirnov test against $N(0,1)$ gives $p = 0.849$. Detection under an explosive alternative (60 replications) is 100%, the same as a standard SADF test on the same DGP. Replication script: [replication/monitoring/radf_lbi_validation.R](#script-radf_lbi_validation).

### Sequential extension

The authors report that "the exponentially weighted CUSUM detector with a constant boundary function turns out to be most powerful" (Section 4.1, eq. 12 and 15, Table 1).

- **Statistic.** Normalise the index to the monitoring period, $r = j/T_m$ for $j = 1, \dots, T_m$, where $T_m$ is the monitoring horizon fixed in advance, and form the weighted partial sum
  $$\mathrm{LBI}_{[rT]} = \frac{1}{\tilde\sigma \sqrt{T_m}} \sum_{t=1}^{[r T_m]} w_t\, \Delta y_t \Rightarrow W(r),$$
  a standard Brownian motion under $H_0$. The Chu–Stinchcombe–White boundary of `monitor_cusum()` grows like $\sqrt t$. Normalising by the fixed $T_m$ lets a single constant boundary control the size uniformly over the monitoring window. The paper calls this variant `mCUSUM`. It is more powerful than the classical time-varying-boundary CUSUM of Brown et al. (1975), because under an explosive alternative the detector tends to be largest near the end of the window, which a boundary with shrinking relative tolerance penalises.
- **Weights.** Eq. 12 is $w_r^{\bar c} = \sqrt{2\bar c/T_m}\,/\sqrt{e^{2\bar c} - 1}\; e^{\bar c r}$, where $\bar c \ge 0$ up-weights later, more bubble-like observations. $\bar c = 0$ gives flat weights (`mCUSUM`), and $\bar c > 0$ gives `wCUSUM`. The authors suggest $\bar c \approx 2$.
- **Critical values.** Table 1 (page 7, 1,000,000 replications at $T = 10{,}000$) gives one-sided asymptotic critical values. The running maximum of a time-changed Brownian motion has the same distribution whatever the time change, $\sup_r W^*(\eta(r)) =_d \sup_r W(r)$, so one set of values covers every $\bar c$: $1.64$, $1.95$, $2.24$, $2.57$ and $2.80$ at the 10%, 5%, 2.5%, 1% and 0.5% levels.
- **Variance.** $\tilde\sigma^2$ is estimated from the training window only (Section 4.2). `lbi_test()` uses the full-sample $\tilde\sigma^2$ in the static case.

`monitor_lbi(data, r_star = 0.5, c_bar = 0, level = 0.95)` is in `exuber/R/lbi_test.R`, with tests in `test-lbi.R`.

**Validation.**

- The sum of squares of the flat weight vector ($\bar c = 0$) equals 1 exactly, as in the discrete form of eq. 12. For $\bar c > 0$ it equals 1 up to the expected Riemann-sum error ($< 0.5\%$ at $T_m = 500$).
- The final-point statistic under `mCUSUM` matches a hand-computed telescoped value with the training-window $\tilde\sigma^2$ to machine precision.
- The table lookups match Table 1 exactly, with a clean error for an untabulated level, and alarms never fire before the end of the training window (50 of 50 replications).
- The false-alarm rate under $H_0$ (1000 replications, $n = 200$, $T^* = 100$) is 3.7% for `mCUSUM` and 4.0% for `wCUSUM`, against a nominal 5%.
- Detection on a post-training bubble (60 replications) is 41% for `mCUSUM` and 44% for `wCUSUM`, above the 31% of `monitor_cusum(type = "standard")` on the same DGP. This confirms the paper's claim that the constant-boundary LBI detector is more powerful than the classical CUSUM, and `wCUSUM` is at least as powerful as `mCUSUM`.

Replication script: [replication/monitoring/radf_lbi_monitor_validation.R](#script-radf_lbi_monitor_validation).

Not implemented: the Table 1 row for the classical time-varying-boundary CUSUM of Brown et al. (1975), which the paper uses as a comparison, and the DF-statistic monitoring variant of Section 4.2, which is a `badf`-based analogue with a different boundary and addresses the same problem as `monitor()`.

## Whitehouse, Harvey & Leybourne (2025): AHLST decision rule and FPR

The DGP is $y_t = \mu + u_t$, with $u_t = u_{t-1} + \varepsilon_t$ for $t \le \lfloor \tau T \rfloor$ and $u_t = (1+\delta) u_{t-1} + \varepsilon_t$ afterwards. The statistic (eq. 2) is

$$
A_{e,k} = \frac{B_{e,k}}{\sqrt{C_{e,k}}}, \qquad B_{e,k} = \sum_{t=e-k+1}^{e} (t - e + k)\,\Delta y_t, \qquad C_{e,k} = \sum_{t=e-k+1}^{e} \bigl\{(t - e + k)\,\Delta y_t\bigr\}^2 .
$$

The training-sample maximum $A^*_{\max} = \max_{e \in [k+1, T^*]} A_{e,k}$ is the critical value, and the rule "reject $H_0$ at time $e$ if $A_{e,k} > A^*_{\max}$" defines the $\mathrm{AMAX}(k)$ procedure. Under $H_0$, for a monitoring point $T'$ (eq. 4–5),

$$
\lim_{T \to \infty} P\Bigl(\max_{e \in [T^*+k, T']} A_{e,k} > \max_{e \in [k+1, T^*]} A_{e,k}\Bigr) = \tau = \lim \frac{T' - T^*}{T'} .
$$

The approximate FPR at $T'$ (eq. 6) is

$$
\alpha \approx \frac{T' - T^* - k + 1}{T' - 2k + 1},
$$

so monitoring can run until $T' \approx (T^* + k - 1 - \alpha(2k-1))/(1-\alpha)$ at a chosen FPR $\alpha$. In contrast to CUSUM-based approaches (Homm & Breitung, Astill et al., Horváth & Trapani), this gives an exact, usable FPR with no asymptotic boundary and no conservatism, but the FPR necessarily grows with the monitoring horizon, so it suits short-range monitoring. CUSUM-style methods can hold a fixed FPR (for example 0.05) over an arbitrarily long horizon, at the cost of lower power (a lower true positive rate).

Table 1 ($k = 10$, NIID and GARCH(1,1) errors) reports an empirical FPR for the baseline $\mathrm{AMAX}(k)$ of 0.006 (NIID) at $T' = 200$, rising monotonically to 0.147 at $T' = 230$. The two variance-standardised variants, $A^{AR,\max}(k)$ and $A^{T,\max}(k)$, run a little higher (0.015 and 0.013 at $T' = 200$), so they are less conservative, which is the design goal of Theorem 1 (the same asymptotic FPR with different finite-sample behaviour).

In the empirical application, $A^{AR,\max}(k)$ detects the bubble in the US house price-to-rent ratio that preceded the 2007/08 financial crisis as early as 1999:Q1, against 2000:Q1 for $\mathrm{AMAX}(k)$, an improvement of four quarters (Table 2).

## Horváth & Trapani (2023/2026): RCA monitoring

Status: evaluated, not implemented.

The WLS-residual CUSUM detector (eq. 2.4) over a training window of length $m$ is

$$
Z_m(k) = \sum_{i=m+1}^{m+k} \frac{(y_i - \hat\theta_m y_{i-1})\, y_{i-1}}{1 + y_{i-1}^2}, \qquad k \ge 1 .
$$

For the open-ended or long-horizon case the boundary function (eq. 2.5/2.9) is

$$
g_{m,\gamma}(k) = c_{\gamma,\alpha}\, s\, m^{1/2} \left(1 + \frac{k}{m}\right) \left(\frac{k}{m+k}\right)^{\gamma}, \qquad 0 \le \gamma < \tfrac12,
$$

and a short-horizon variant (eq. 2.10) is $g_{m,\gamma}(k) = c_{\gamma,\alpha}\, s\, m^{1/2 - \gamma} k$, used when the monitoring horizon $m'$ is $o(m)$. The stopping time is $\tau_{m,\gamma} = \inf\{k \ge 1 : Z_m(k) \ge g_{m,\gamma}(k)\}$. The constant $c_{\gamma,\alpha}$ controls size, like $b_\alpha$ of Homm & Breitung. The paper also defines a Page-CUSUM variant for a shorter detection delay.

Table 5.4 (median detection delay, no covariates, $m = 200$): in Case I ($\delta_0 = 0.5$) the standard weighted CUSUM ($\gamma = 0$) has a median delay of 54 and the standardised CUSUM ($c_{\gamma,0.5}$) has 37, a reduction of about 30%, at the price of lower empirical power (a rejection frequency of 0.465 against 0.705). In the application to Los Angeles daily housing prices ($m, m' \in \{100, 200\}$), the ex-post analysis dates the break at 4 February 2009. The real-time procedure with no covariates flags a change point on 2–15 June 2009, depending on the windows, a delay of about four months. Adding covariates (interest-rate proxies, VXO, the Weekly Economic Indicator) moves the flag to 18 May 2009 in the richest specification (Table 6.2).

**Why it is not implemented.** The statistic $Z_m(k)$ is a single `cumsum()` after an OLS coefficient $\hat\theta_m$ from the training window, so it is as cheap as `monitor_cusum()`. The obstacles are the critical values and a nuisance parameter.

- Several boundary regimes need their own critical values: open-ended (eq. 2.5), closed-ended long-horizon (eq. 2.9) and closed-ended short-horizon (eq. 2.10).
- For $\gamma < 1/2$ the critical value solves a Brownian-motion sup-norm probability, $P(\sup |W(u)| < c_\gamma)$ (eq. 3.4). For $\gamma = 1/2$ it follows a Darling–Erdős extreme-value asymptotic (eq. 3.5–3.6), $c_{\alpha,0.5} = [x + b(\log m)]/a(\log m)$ with $a(x) = \sqrt{2 \log x}$, $b(x) = 2 \log x + \tfrac12 \log\log x - \tfrac12 \log \pi$ and $x$ solved from $\exp(-\exp(-x)) = 1 - \alpha$. The authors note that these asymptotic values are "bound to be inaccurate due to the slow convergence to the Extreme Value distribution... leading to low power", and propose a finite-sample correction (eq. 3.7–3.8) that requires solving an implicit equation for $c$ with a tuning parameter $h_m$ (recommended $h_m = \sqrt{\log m}$).
- At $\gamma = 0$, the limiting probabilities of Theorems 3.1 and 3.3 reduce to $P\{\sup_{0 < u \le 1} |W(u)| < c_{\alpha,0}\}$, the classical sup-norm distribution behind the two-sided Kolmogorov–Smirnov statistic, which has a closed-form alternating series and can be inverted with `uniroot()`. The general case $\gamma \in (0, 1/2)$ has no such form, and the paper obtains its critical values "by simulation".
- The normalising constant $\mathfrak s^2$ of the boundary (eq. 2.6) is defined by a case split on the Lyapunov-type exponent $E[\log |\beta_0 + \varepsilon_{0,1}|]$. For a plain unit root ($\beta_0 = 1$, i.i.d. innovations) the exponent is negative (the paper's Case III gives $-0.007$, "the STUR model"), and then $\mathfrak s^2 = a_1\sigma_1^2 + a_2\sigma_2^2$ with $a_1 = E[(\bar y_0^2/(1 + \bar y_0^2))^2]$ and $a_2 = E[(\bar y_0/(1 + \bar y_0^2))^2]$, expectations over the stationary distribution of the RCA(1) process, which has no general closed form. The Monte Carlo of the paper computes critical values from the true DGP parameters (Theorems 3.2/3.6), and the paper gives no estimator that works on data with unknown parameters. A usable implementation needs such an estimator.

## Remaining items

- A closed-form false-alarm-versus-horizon boundary for the wild-bootstrap route of `monitor()`, which still recalibrates by simulation.
- A non-bootstrap (asymptotic or Monte Carlo) training critical value for `monitor()`, which `radf_mc_cv()` could supply.
- The Page-CUSUM family of Horváth & Trapani.
