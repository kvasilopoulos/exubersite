---
title: "Real-time monitoring for bubbles"
blurb: "Sequential and real-time detection: training-vs-monitoring orchestration, CUSUM families, and closed-form boundaries."
order: 3
---
Status: the Phillips & Shi (2020) procedure (Family A, `monitor()`) was
implemented on 2026-08-09. On 2026-08-10 we added the CUSUM procedure of
Family B: Homm & Breitung (2012)'s original statistic and the
volatility-robust "CUSUMV" kernel variant of Astill et al. (2023),
`monitor_cusum(..., type = "standard"/"kernel")`, together with HB's own
finite-sample boundary (`monitor_cusum(..., boundary = "finite")`). The
same day we added Kurozumi (2020)'s closed-form `SADF` boundary and its
`GSADF_{s0}` generalization (`monitor(..., boundary = "kurozumi",
s0 = 0/0.4/0.8)`), HB's second statistic FLUC (`monitor(...,
boundary = "fluc")`), and Breitung & Diegel (2025)'s static LBI test
(`lbi_test()`) with their sequential extension (`monitor_lbi()`, constant
boundary `mCUSUM`/`wCUSUM`). Two items were read but not implemented:
Horváth and Trapani's RCA framework, which still lacks a feasible
nuisance-parameter estimator (see below), and Kurozumi's 2021 delay-time
paper. "Implementation" below describes what shipped, and "Cost/feasibility
note" describes what did not and why.

Monitoring is the largest item in this research programme. It spans two
structurally different detector families (recursive training-maximum
detectors and CUSUM/Page-CUSUM detectors), and building both with their own
S3 infrastructure and size-control theory is still estimated at three to
five times the cost of the [SBZ item](/replication/volatility-robustness#sbz-wls--kernel-volatility).
The one sub-item that changed the estimate was Family A. The earlier
cost note had marked it as the only place where reuse of existing exuber
machinery was real and not hoped for: the Phillips & Shi (2020)
wild-bootstrap training critical value, already implemented as
`radf_wb_ps_cv()`. It needed only a thin orchestration wrapper and no new
statistical theory, and it is now built.

## Implementation

Family A shipped as `monitor(data, r_star = 0.5, minw, nboot, level,
adflag, type, seed)` in a new file, `exuber/R/monitor.R`. It matches the
cost note's prediction that the work would be "a loop and a stopping
condition around machinery that already exists."

Two structural facts made it nearly free once we tried it.

1. `radf_wb_ps_cv(..., tb = T*)` already computes the training-window
   wild-bootstrap critical value that monitoring needs, and it already
   broadcasts that value as a constant boundary across the whole
   monitoring horizon. We checked this by passing the full-length series
   with `tb = T*`: the returned `bsadf_cv` matrix has `nrow = n - minw`,
   not `T* - minw`. The code was built to produce a horizon-length
   boundary and had simply never been wired into a monitoring workflow.
2. `radf()`'s BSADF statistic at calendar time `t` depends only on data up
   to `t`. We verified that it is bit-identical to the last BSADF value of
   a fresh `radf(y[1:t])` call. The whole monitoring path therefore comes
   from a single full-sample `radf()` call, which is `O(T)` because of the
   efficient recursive computation, and we do not need to refit at every
   monitoring point.

One design decision is worth recording because it avoids look-ahead
leakage. `monitor()` calls `radf_wb_ps_cv()` on `data[1:T*]` only, not on
the full series with `tb = T*`. The null-model fit inside
`radf_wb_ps_cv()` (`adf_res()`, in `radf_wb.R`) uses all the data it
receives to estimate the bootstrap DGP's residuals and coefficients. The
argument `tb` only truncates the length of the simulated bootstrap sample,
and it does not change which part of the real data feeds the residual
estimation. Passing the full series would therefore let post-`T*`,
possibly explosive, data leak into the training-window null calibration.
That would have been a correctness error in a naive "just call
`radf_wb_ps_cv(full_data, tb = T*)`" reading of the cost note. Slicing to
`data[1:T*]` first avoids it, at the price of one extra step that
`monitor()` takes on the caller's behalf.

Validation:

1. Structural checks. A basic run returns a well-formed object, and
   `T_star` and the `bsadf` row count (`n - minw`) match expectations. An
   alarm, when raised, always falls strictly after `T*`. We checked this
   over 10 seeds with a post-training bubble, and monitoring never fired
   on training data.
2. False-alarm rate under H0 (no bubble anywhere, 40 replications,
   75-observation monitoring horizon, 95% per-point threshold): 10%. This
   is above the 5% per-point nominal level, as it should be. The 10% is a
   cumulative false-alarm probability over 75 sequential comparisons
   against a fixed boundary, and the false-alarm rate of Family A or
   PSY-style monitoring grows with the monitoring horizon. The
   AHLST/Whitehouse discussion below documents the same property through
   its closed-form eq. 6, which we have not implemented (see item 2 of the
   cost/feasibility note).
3. Detection power under a bubble that starts strictly after `T*` (15
   replications): 86.7%. The alarm delay, meaning the date of the alarm
   minus the true origination date, was always positive and bounded, with
   a minimum of 6, a median of 19 and a maximum of 33 observations. A
   positive delay is expected for any threshold-crossing procedure, since
   the BSADF statistic needs several post-origination observations to
   accumulate enough evidence to cross a fixed boundary.
4. A bug found and fixed while writing the tests. The training boundary
   vector lost its series names when we extracted it from a single-series
   critical-value matrix (`cv$gsadf_cv[, "95%"]` drops to an unnamed scalar
   under R's default `drop = TRUE`). That broke the print method, because
   the `data.frame()` construction failed on a 0-versus-1 row-count
   mismatch. We now re-attach the names explicitly with
   `setNames(cv$gsadf_cv[, lvl_lab], snames)`.

The new tests are in `test-monitor.R`, and the full package suite passes.

Replication script:
[replication/monitoring/radf_monitor_validation.R](#script-radf_monitor_validation).

Still out of scope: the CUSUM/Page-CUSUM detector family beyond the one
item below, any closed-form or simulated false-alarm-rate-versus-horizon
boundary function for Family A (AHLST/Whitehouse's eq. 4-6), a
union-of-rejections combination across families, and date-stamping methods
specific to a monitoring result. See "Cost/feasibility note" below.

## Implementation (CUSUM)

Status: done (2026-08-10). Shipped as `monitor_cusum()`.

Homm & Breitung (2012)'s CUSUM procedure (their Section 3, eq. 26-30) has
no structural relation to Family A. It is not a recursive ADF regression.
It is a standardized running sum of first differences, compared against a
closed-form asymptotic boundary derived from the inequality in Chu,
Stinchcombe & White (1996) (eq. 28). It needs no wild bootstrap, no Monte
Carlo simulation and no new dependency, which makes it the cheapest item
in the monitoring bundle. That matches the original cost note word for
word: "the CUSUM partial-sum statistic itself [is] trivial, one `cumsum()`
in R, no C++ needed."

HB propose two monitoring statistics, CUSUM and FLUC. We confirmed this
by reading the primary source, since the earlier draft relied on a
restatement. Both are implemented. CUSUM is `monitor_cusum()`, described
here, and FLUC is `monitor(..., boundary = "fluc")` (2026-08-10). Its
point statistic turns out to be `radf()`'s own `badf` sequence and not a
new one (see "Implementation (FLUC)" below).

The statistic (eq. 26, 29-30), for a training window ending at `T*` and a
monitoring point `t > T*`:

```
S_t = (y_t - y_{T*}) / sigma_hat_t         (telescoped cumulative sum of post-training first differences)
sigma_hat_t^2 = (t-1)^-1 * sum_{j=2}^{t} (Delta y_j)^2   (recursive sample variance of ALL differences up to t)
c_t = sqrt(b_alpha + log(t / T*))
boundary_t = c_t * sqrt(t)
reject H0 (alarm) at the first t where S_t > boundary_t
```

The constant `b_alpha = 4.6` is HB's one-sided asymptotic calibration for
a 5% significance level. In contrast to `monitor()`'s wild-bootstrap
boundary, `sigma_hat_t` is re-estimated from all data up to the current
monitoring point `t`, and not only from the training window. This causes
no look-ahead problem of the kind we found for `radf_wb_ps_cv()`, because
at each real monitoring instant `t` only data up to `t` is ever used, so
no future information can leak.

### Implementation shape

The code is an internal `cusum_stat_path()` and an exported
`monitor_cusum(data, r_star, b_alpha)`, which returns a `radf_cusum_obj`.
It reuses no code from Family A, from `radf()` or from `exubercore`. We
predicted this and then confirmed it by building the function: the
statistic shares no structure with the recursive-ADF family, which matches
the original cost note's assessment of the CUSUM/Page-CUSUM family as 0%
reusable and all new.

### Independent validation

1. Formula check. `cusum_stat_path()` matches an independent brute-force
   loop recomputation, with separate `diff()` and `sum()` calls at each
   monitoring point instead of the vectorized `cumsum()`, to floating-point
   precision (maximum absolute difference of `0`).
2. False-alarm rate under H0 (pure random walk, no bubble, 100
   replications, 75-observation monitoring horizon): 0%. This agrees with
   HB's eq. 28 being a conservative asymptotic upper bound, via Chu,
   Stinchcombe & White (1996), and not an exact size. A rate well below
   the nominal 5% is the expected signature of a conservative bound and
   does not indicate miscalibration.
3. Detection power under the same post-training bubble DGP that we used
   to validate `monitor()` (30 replications): a detection rate of 30% and
   a median alarm delay of 27 observations. This is substantially lower
   than Family A's 86.7% detection rate and 19-observation median delay on
   the identical DGP. We report it as it came out. It agrees with the
   literature cited elsewhere in this file: Kurozumi (2020, 2021) reports
   that CUSUM-type detectors do better than ADF-type ones for early, short
   bubbles and worse for middle-to-late ones. Our test DGP starts the
   bubble about 65% of the way into the sample, which is the
   middle-to-late regime where the literature predicts CUSUM should lag
   ADF-type detection, and it does. Neither implementation is at fault.
   The result also illustrates why the field recommends a
   union-of-rejections strategy that runs both families at once, which we
   have not implemented.

The new tests are in `test-cusum.R`, and the full package suite passes.

Replication script:
[replication/monitoring/radf_cusum_validation.R](#script-radf_cusum_validation).

`boundary = "finite"` (2026-08-10): HB's finite-sample boundary constant
(their Table 8, "without drift estimation," which matches this
implementation's raw-first-difference construction) replaces the fixed
asymptotic `b_alpha = 4.6`. It is the same table-lookup shortcut we used
for Kurozumi's and FLUC's boundaries. We transcribed the values directly
from HB's published table, which is indexed by training length,
significance level and monitoring-horizon ratio `k = N/T*`, and ran no new
simulation. The option applies to both `type = "standard"` and `type =
"kernel"` (CUSUMV). Corollary 1 of Astill et al., already validated below,
shows that the same boundary function works for both statistics, so
extending the finite-sample table to CUSUMV follows directly from that
result and needs no separate validation.

We checked that the table lookups are exact. At `T*=75` and `n=150`
(`k=2`, with `T*` snapped to the nearest tabulated `n=50`), the
finite-sample boundary gives a false-alarm rate of 9.0% under H0. That is
closer to the nominal 5% level than the 0% of the asymptotic bound, and
it comes with higher detection power on the same post-training bubble DGP
(36.7% against 30.0%). The improvement is real and not cosmetic, because
the constant is now calibrated for finite samples instead of being an
asymptotic upper bound. New tests extend `test-cusum.R`. Replication
script:
[replication/monitoring/radf_cusum_finite_boundary_validation.R](#script-radf_cusum_finite_boundary_validation).

### CUSUMV: the volatility-robust variant (Astill et al. 2023)

Status: done (2026-08-10). Shipped as `monitor_cusum(..., type = "kernel")`.

Astill, Harvey, Leybourne, Taylor & Zu (2023, *JFEC* 21(1), 187-227;
"AHLTZ") generalize HB's CUSUM procedure to allow time-varying volatility.
Their abstract states that "such behavior can heavily inflate the false
positive rate (FPR) of the CUSUM-based procedure." Their fix (eq. 6-7)
replaces HB's single running variance with a one-sided Nadaraya-Watson
kernel-weighted estimate of the spot variance. One-sided means causal: it
uses only current and past lags, as real-time monitoring requires. Each
first difference is standardized individually before the cumulation:

```
SV_t = sum_{j=T*+1}^{t} Delta y_j / sigma_hat_{j,N}     (eq. 6)
sigma_hat_j,N = sum_{s=0}^{N} w_s * Delta y_{j-s}^2,   w_s = K(s/N) / sum_s K(s/N)   (eq. 7)
```

Their Corollary 1 proves that the same boundary function `c_t * sqrt(t)`
(with the same `b_α`) used for HB's original statistic still gives a
controlled asymptotic false-alarm rate for the modified statistic, even
under time-varying volatility. The boundary stays as it was, and only the
construction of the statistic changes.

Implementation: `one_sided_kernel_spot_vol()` computes a fixed causal
kernel weight vector with `stats::filter(..., sides = 1)` and follows
AHLTZ's convention `σ̂²_{j,N} := 1` for `j <= N`. It is paired with
`cusum_stat_path_kernel()`. Both are wired into `monitor_cusum()` as a
`type = "kernel"` option that reuses the existing boundary and
decision-rule logic unchanged, so there is no separate function, since
only the numerator of the statistic differs. The default bandwidth is
`N = 20`, which is the value AHLTZ recommend empirically ("setting H = 20
delivered a procedure with the best trade-off" between FPR robustness and
power). We did not implement their data-driven local cross-validation
bandwidth selection (eq. 8-9). This is a documented simplification, and we
do not claim to match their exact finite-sample procedure.

Independent validation:

1. Formula check. `one_sided_kernel_spot_vol()` and
   `cusum_stat_path_kernel()` both match independent brute-force loop
   recomputations to floating-point precision.
2. Under a homoskedastic H0, `standard` and `kernel` give comparable
   false-alarm rates (here both are 0%, equally conservative). This
   agrees with AHLTZ's Remark 8: "in the case where the innovations are
   homoskedastic... [both] lead to the same limiting null distribution."
3. Under a heteroskedastic H0 (a volatility jump from 1 to 8 partway
   through the monitoring region, 60 replications), the false-alarm rate
   of `standard` CUSUM rises to 8.3%, visibly above its homoskedastic
   rate. This is the failure mode that AHLTZ's abstract describes. The
   `kernel` CUSUMV stays at 0%. The paper's central claim therefore holds
   up in independent simulation and not only in asymptotic theory.
4. Detection power under the same post-training bubble DGP (30
   replications, homoskedastic): 30.0% for `standard` and 36.7% for
   `kernel`. We observed no power sacrifice in this scenario. AHLTZ claim
   that the modification "sacrifices only a small amount of power... when
   the shocks are homoskedastic," and our result is consistent with that,
   though slightly better.

New tests extend `test-cusum.R`, and the full package suite passes.

Replication script:
[replication/monitoring/radf_cusumv_kernel_validation.R](#script-radf_cusumv_kernel_validation).

### Implementation (FLUC)

Status: done (2026-08-10). Shipped as `monitor(..., boundary = "fluc")`.

Homm & Breitung's second monitoring statistic is their eq. 27. We checked
it by rendering the PDF page. It reads `Z_t = (rho_hat_t - 1) /
sigma_hat_{rho_t} = DF_{t/n}`, the ordinary recursive (expanding-window)
OLS ADF t-statistic on the sample `{y_0, ..., y_t}`. That is exactly
`radf()`'s existing `badf` sequence, the same statistic that, as the
Kurozumi subsection below shows, `monitor(..., boundary = "kurozumi")`
already reuses. No new point statistic is needed.

The rejection rule (eq. 29/31) has the same form as CUSUM's boundary,
`DF_{t/n} > kappa_t` with `kappa_t = sqrt(b_{k,alpha} + log(t/n))`, but it
uses a different calibration constant `b_{k,alpha}`. Unlike CUSUM's
closed-form `b_alpha = 4.6`, HB's text says they "determine... by means of
simulation," which puts it in the cost tier that the cost note first
flagged: new critical-value work. It turned out to be cheap anyway because
HB publish the simulated constant (their Table 7, part i, "without
detrending," which matches `radf()`'s no-trend default). It is tabulated by
training length `n` in `{20, 50, 100}`, by significance level `alpha` in
`{0.10, 0.05, 0.01}`, and by monitoring-horizon ratio `k = N/n` in `{2, 3,
4, 5, 6, 8, 10}`, where `N` is the total sample including training. The
lookup is closed-form and needs no simulation on exuber's side, the same
shortcut that Kurozumi's Table 1 gave for `boundary = "kurozumi"`.

Implementation: `hb_fluc_table` is a 9x7 table transcribed from the
rendered PDF page. `hb_fluc_q(level, n_train, k)` snaps `n_train` to the
nearest of `{20, 50, 100}` and `k` to the nearest of `{2,...,10}`, and it
requires `level` to match one of the three tabulated significance levels
exactly. It follows the same lookup-and-snap pattern as
`kurozumi_sadf_q()` and is wired into `monitor()` as a third `boundary`
option alongside `"bootstrap"` and `"kurozumi"`.

Validation: the table lookups match Table 7 exactly (6 checked cells,
including a tie-breaking snap case). A basic run prints without error, and
alarms never fire before `T*` (a structural invariant, checked in 10 of 10
replications). The false-alarm rate under H0 (`n=150`, `T*=75`, `k=2`, the
smallest and most conservative tabulated horizon ratio, 100 replications)
is 0%. HB's table is a conservative finite-sample-simulated bound at this
`k` and is not miscalibrated, which is comparable to the conservative
behavior of CUSUM noted above. Detection power under the same
post-training bubble DGP used elsewhere in this file (30 replications) is
56.7%. That is below the 80% of `boundary = "kurozumi"` and the 90% of
`"bootstrap"` on the identical DGP, and it agrees with HB's own finding
that FLUC and CUSUM monitoring generally have less power than a supDF-style
test, although FLUC beats their CUSUM.

One caveat. HB's table covers training lengths only up to `n=100`, which
is smaller than typical financial series, so larger `T*` values are
snapped to the `n=100` row and not interpolated or extrapolated. This is a
documented approximation and does not give exact finite-sample calibration
at realistic sample sizes. New tests extend `test-monitor.R`, and the full
package suite passes.

Replication script:
[replication/monitoring/radf_monitor_fluc_boundary_validation.R](#script-radf_monitor_fluc_boundary_validation).

## Source

There are seven items. On 2026-08-09, institutional access recovered the
primary-source PDFs for items 1, 2, 4 and 7. When this evaluation was first
written, only items 3, 5 and 6 were open. The rest were either known at
abstract level only (Kurozumi) or read through the restatement in item 6
(Homm & Breitung, Astill 2021/2023). The formulas below that came from the
restatement have not yet been re-verified against the primary PDFs, and
that re-verification remains follow-up work.

1. Homm, U. & Breitung, J. (2012). "Testing for speculative bubbles in
   stock markets: a comparison of alternative methods." *Journal of
   Financial Econometrics*, 10(1), 198–231. `doi:10.1093/jjfinec/nbr009`.
   We now have access through our institution, and the local copy is
   The formulas below (the CUSUM statistic, the boundary function and the
   Table 8 finite-sample `b_α` calibration) still come from the
   restatement in item 6 and have not been checked against this primary
   copy.

2. Astill, S., Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y.
   (2021/2023). "CUSUM-Based Monitoring for Explosive Episodes in
   Financial Data in the Presence of Time-Varying Volatility." *Journal of
   Financial Econometrics*, 21(1), 187–227 (advance access March 2021,
   print issue Winter 2023). `doi:10.1093/jjfinec/nbab009`. The literature
   cites it as "Astill et al. (2021)" or "(2024)", depending on which
   author's own citation list you use, but it is the same paper,
   `nbab009`. We have access through our institution, and the local copy
   is
   The formulas below are still taken from item 6's restatement (labelled
   "AHLTZ") and have not been checked against this primary copy.

   An earlier companion paper by an overlapping author set, Astill, S.,
   Harvey, D.I., Leybourne, S.J., Sollis, R. & Taylor, A.M.R. (2018),
   "Real-Time Monitoring for Explosive Financial Bubbles," *Journal of
   Time Series Analysis*, 39, 863–891 (often abbreviated AHLST), uses a
   different monitoring design. It compares against the maximum of the
   training sample instead of using a CUSUM, and it is essential
   background because item 3 below builds directly on it. We have access
   through our institution, and the local copy is
   We have read it only through its restatement in item 3 and have not
   yet re-read it directly.

3. Whitehouse, E.J., Harvey, D.I. & Leybourne, S.J. (2025). "Real-time
   monitoring procedures for early detection of bubbles." *International
   Journal of Forecasting*, 41(3), 1260–1277.
   `doi:10.1016/j.ijforecast.2024.12.005`. The article is open access
   (CC-BY), and we downloaded it directly from the White Rose Research
   Online repository. The local copy is
   We read it with `pdftotext -layout`. For pages 2–4 (the introduction,
   the AHLST model and decision rule, and equations 2–9) we also rendered
   the pages to PNG at 2.5x (PyMuPDF) and read the typeset math directly,
   because the OCR text mangled subscripts and superscripts (for example
   `A_{e,k}` came out as `Ae,k` merged into the surrounding prose) and
   scrambled the summation limits of the FPR formula across the two-column
   layout. This is the paper that supplied the AHLST decision rule and FPR
   formula used below.

4. Kurozumi, E. Two papers, one on monitoring critical values and one on
   detection delay.
   - Kurozumi, E. (2020). "Asymptotic properties of bubble monitoring
     tests." *Econometric Reviews*, 39(5), 510–538.
     `doi:10.1080/07474938.2019.1697086`. It extends SADF and GSADF to a
     monitoring scheme and studies a CUSUM detector alongside them. It
     derives new monitoring-period critical values and compares ADF-type
     and CUSUM-type detection under moderate-deviation and local-to-unity
     asymptotics. We have access through our institution, and the local
     copy is
   - Kurozumi, E. (2021). "Asymptotic Behavior of Delay Times of Bubble
     Monitoring Tests." *Journal of Time Series Analysis*, 42(3), 314–337.
     `doi:10.1111/jtsa.12569`. It is specifically about detection delay,
     namely the stochastic order of the stopping time, for ADF-type and
     CUSUM-type monitoring statistics. We have access through our
     institution, and the local copy is

   What we report below for Kurozumi (2021) is still abstract-level only,
   cross-checked against a secondary source (Skrobotov 2023). Now that
   both PDFs are on hand, the most useful follow-up in this file is to
   read them directly and to replace that qualitative summary with
   verified formulas and numbers.

5. Horvath, L. & Trapani, L. "Real-time monitoring with RCA models."
   Working paper arXiv:2312.11710 (Dec 2023), published as Horváth, L. &
   Trapani, L. (2026), "Real-time monitoring with RCA models,"
   *Econometric Theory*, 42, 514–547. The paper is open access, and we
   downloaded it directly from arXiv. The local copy is
   We read it with `pdftotext -layout`. We rendered the pages with the
   core detector and boundary-function definitions (section 2, equations
   2.4–2.11) to PNG and checked them there, because the OCR badly mangled
   the summation- and subscript-heavy statistic definitions.

6. Astill, S., Taylor, A.M.R. & Zu, Y. (2026, forthcoming). "Covariate
   Augmented CUSUM Bubble Monitoring Procedures." *Econometric Theory*.
   It is an open-access working paper, Essex Finance Centre Working Paper
   No. 94. The local copy is
   It was the most useful document we obtained in the original pass. Its
   Section 3 ("CUSUM-based Bubble Detection Procedures") restates, with
   full equation numbers and explicit page citations back to the originals
   (for example "HB... Table 8, p221"), both the original Homm & Breitung
   (2012) CUSUM statistic and boundary function (item 1) and the Astill et
   al. (2021/2023) volatility-robust modification (item 2). Two of its
   authors (Taylor and Zu) are co-authors of item 2. We rendered the
   formula pages (pp. 10–11 of the PDF) to PNG and read them directly,
   because the OCR text mangled the subscripts on `S^t_T`, `SV^t_T` and
   `c_t`, and the kernel-weight definitions, too badly to use.

7. Breitung, J. & Diegel, M. (2025). "Sequential Detector Statistics for
   Speculative Bubbles." *JTSA*, 46(5). The local copy is
   We found it in the *JTSA* 46(5) special issue, and it matches the intent
   of this section directly and can be cited. We have not yet read it in
   depth. A correction to an earlier lead: the Horvath and Trapani 2025
   JTSA piece in the same special issue, "Sequential Monitoring for
   Changes in GARCH(1,1) Models Without Assuming Stationarity," is about
   GARCH change-point monitoring in general and not about bubble
   monitoring. It is tangential to this section, and it is a different
   Horváth & Trapani paper from item 5.

## What it is

All seven papers address the same underlying problem, and they fall into
two structurally different families.

Family A is the recursive/BSADF comparison of a training maximum with the
monitoring statistic (AHLST and Whitehouse). We split the sample into a
training period `t = 1,...,T*`, assumed to be free of bubbles, and a
monitoring period `t = T*+1,...,T`. We compute a recursive statistic
`A_{e,k}` over rolling sub-samples of fixed length `k` in both periods.
AHLST's statistic is a sub-sample regression of `Δy_t` on a linear trend,
studentized with a White-type variance estimator. It is structurally close
to exuber's ADF-family statistics but not identical to them. The maximum
over the training sample, `A*_max = max A_{e,k}`, becomes the fixed
critical value, and monitoring rejects `H0` at the first `e` where
`A_{e,k} > A*_max`. AHLST prove a closed-form asymptotic FPR that is a
simple function of the ratio of training length to monitoring length (the
formula is eq. 6 below). This is the "PSY-style" monitoring philosophy to
which Phillips & Shi (2020) also belong, and exuber has already
implemented part of it (see Cost/feasibility), although PSY use a
bootstrap and not a closed-form FPR.

Family B is the family of full-sample CUSUM and Page-CUSUM detectors
(Homm-Breitung, Astill et al. 2021/2026, Kurozumi's CUSUM variant,
Horvath-Trapani and Breitung-Diegel). We compute a partial-sum (CUSUM)
statistic of the standardized first differences `Δy_t` from the end of the
training sample onward. At each monitoring date we compare it with a
boundary function that grows with `t` (for example `c_t·√t`), chosen so
that the cumulative false-alarm probability over the whole, possibly
infinite, monitoring horizon stays below a target `α`. This is a different
kind of statistic. It is a running standardized sum and not a recursive
ADF regression. Its volatility-robust variants (Astill et al.) replace the
standardization with a kernel spot-variance estimate, which resembles the
kernel machinery of [SBZ/STADF](/replication/volatility-robustness#sbz-wls--kernel-volatility)
in spirit although it serves a different target. Horvath-Trapani extend
this further to a Random-Coefficient-Autoregressive (RCA) framework with
weighted CUSUM and Page-CUSUM detectors. These detectors work
"symmetrically" for transitions from stationary to explosive and from
explosive to stationary, and they do not require us to know in advance
which regime we start in.

Kurozumi's contribution spans both families. He puts SADF/GSADF-type
(Family A) and CUSUM-type (Family B) monitoring statistics into one
asymptotic framework, derives new monitoring-period critical values for
both, and separately studies the distribution of the detection delay (the
stopping time). According to the secondary source, CUSUM detects an early,
short bubble faster, while ADF/BSADF-type detectors detect a
middle-to-late bubble faster. A union of rejections that combines BSADF and
CUSUM is possible, and it is the monitoring counterpart of the [SBZ union
statistic](/replication/volatility-robustness#sbz-wls--kernel-volatility).

## Exact numbers/formulas reproduced

### Homm & Breitung (2012) CUSUM statistic and boundary (via item 6, not yet re-verified against the now-available primary copy)

The CUSUM statistic (their eq. 7, with training sample `t=1,...,T*` and
monitoring points `t>T*`) is

```
S^t_{T*} := (1/σ̃_t) · Σ_{j=T*+1}^{t} Δy_j
```

with `σ̃_t² := (t-1)^{-1} Σ_{j=2}^{t} (Δy_j)²`. Under `H0`, `T*^{-1/2}
S^{⌊Tr⌋}_{T*} ⇒ W(r) − W(1)` (eq. 8). By Theorem 3.4 of Chu et al. (1996),
for any `λ>1`:

```
lim_{T→∞} Pr(|S^t_{T*}| > c_t·√t for some t ∈ {T*+1,...,⌊λT*⌋}) ≤ exp(−b_α/2)   (eq. 9)
```

with the boundary function `c_t := √(b_α + log(t/T*))`. We reject `H0` if
`S^t_{T*} > c_t·√t` and flag the first `t` at which this happens. For a
one-sided test of size `α = 0.05`, the asymptotic setting is `b_α = 4.6`,
which gives a two-sided size of at most 0.10 from eq. 9. This asymptotic
setting assumes an infinite monitoring horizon, and the authors of item 6,
citing HB directly, describe it as "extremely conservative in practice."
HB's paper gives finite-sample `b_α` values in its Table 8 (p.221),
calibrated for target FPRs in {0.10, 0.05, 0.01} at specific training and
monitoring lengths. We transcribed and implemented them on 2026-08-10 as
`monitor_cusum(..., boundary = "finite")`, described in "Implementation
(CUSUM)" above.

### Astill et al. (2021/2023) volatility-robust modification (via item 6, not yet re-verified against the now-available primary copy)

The modification replaces `S^t_{T*}` with

```
SV^t_{T*} := Σ_{j=T*+1}^{t} Δy_j / σ̂_{j,N},   t > T*
```

where `σ̂²_{j,N}` is a one-sided kernel smoothing estimator of the spot
variance `σ²_j := σ²(j/T)`:

```
σ̂²_{j,N} := Σ_{s=0}^{N} k_s (Δy_{j-s})²,   k_s := K(s/N) / Σ_{s=0}^N K(s/N)
```

The authors prove that the same boundary function `c_t·√t` still gives a
theoretically controlled FPR when volatility is time-varying, at some cost
in power relative to `S^t_{T*}` under homoskedasticity. Their empirical
application uses Bitcoin prices (this is taken from the OUP abstract, and
we have not yet verified it against the primary copy).

### Whitehouse, Harvey & Leybourne (2025), AHLST decision rule and FPR (verified against rendered PDF pages 2–4)

The DGP is `y_t = μ + u_t`, with `u_t = u_{t-1}+ε_t` for `t ≤ ⌊τT⌋` and
`(1+δ)u_{t-1}+ε_t` afterwards. The statistic (their eq. 2) is

```
A_{e,k} = B_{e,k} / √C_{e,k},   B_{e,k} = Σ_{t=e-k+1}^{e} (t-e+k)Δy_t,   C_{e,k} = Σ_{t=e-k+1}^{e} {(t-e+k)Δy_t}²
```

The training-sample maximum `A*_max = max_{e∈[k+1,T*]} A_{e,k}` is the
critical value for monitoring. The decision rule is "Reject H0 at time e
if `A_{e,k} > A*_max`," which defines the `AMAX(k)` procedure. Under `H0`,
for an arbitrary monitoring point `T'`:

```
lim_{T→∞} P(max_{e∈[T*+k,T']} A_{e,k} > max_{e∈[k+1,T*]} A_{e,k}) = τ = lim(T'-T*)/T'    (eq. 4-5)
```

The approximate FPR at monitoring point `T'` (eq. 6, the formula we can
use directly for calibration) is

```
α ≈ (T' - T* - k + 1) / (T' - 2k + 1)
```

Rearranging gives the point up to which monitoring can run while keeping
the FPR at a chosen level: `T' ≈ (T* + k - 1 - α(2k-1)) / (1-α)`. The
paper contrasts this with CUSUM-based approaches (HB, Astill et al. 2021,
Horvath & Trapani 2026). AHLST's approach gives an exact, usable FPR
formula with no asymptotic boundary and no conservatism, but the FPR
necessarily grows with the monitoring horizon, so it suits short-range
monitoring better. CUSUM-style methods can be tuned to hold a fixed FPR
(for example 0.05) over an arbitrarily long horizon, at the cost of lower
power (a lower true positive rate, TPR).

Empirical Monte Carlo numbers (Table 1, `k=10`, NIID and GARCH(1,1)
errors). We spot-checked a sample of rows with `pdftotext -layout`. The
table is plain numeric with no problematic subscripts, so we did not
verify it separately against a PNG. At `T'=200`, the empirical FPR of the
baseline `AMAX(k)` is 0.006 (NIID), and it rises monotonically to 0.147 by
`T'=230`. The two new variance-standardized variants (`A^{AR,max}(k)` and
`A^{T,max}(k)`) run consistently a little above the baseline (for example
0.015 and 0.013 against 0.006 at `T'=200`). They are slightly less
conservative, which is the paper's stated design goal (Theorem 1: the same
asymptotic FPR with different finite-sample behaviour).

Empirical application: the new `A^{AR,max}(k)` procedure detects the
bubble in the US house price-to-rent ratio that preceded the 2007/08 GFC
as early as 1999:Q1, against 2000:Q1 for the baseline `AMAX(k)`. This is
an improvement of four quarters (Table 2). We verified it in the extracted
text and cross-checked it against two separate passages in the paper that
give the same 1999:Q1 and 2000:Q1 dates.

### Horvath & Trapani (2023/2026) RCA monitoring, evaluated, not implemented (2026-08-10)

Status: evaluated and not implemented. We had already checked the formulas
below against rendered PDF pages. On 2026-08-10 we also read the
asymptotic-theory section (Theorems 3.3-3.5, eq. 3.4-3.8) directly, in
order to assess the feasibility of an implementation and not only to
transcribe the detector formula.

The WLS-residual CUSUM detector (their eq. 2.4) over a training window of
length `m` is

```
Z_m(k) = Σ_{i=m+1}^{m+k} [(y_i - θ̂_m y_{i-1}) y_{i-1}] / (1+y_{i-1}²),   k ≥ 1
```

The boundary function for the open-ended or long-horizon case (eq. 2.5/2.9)
is

```
g_{m,γ}(k) = c_{γ,α}·s·m^{1/2}·(1+k/m)·(k/(m+k))^γ,   0 ≤ γ < 1/2
```

A separate short-horizon variant (eq. 2.10), `g_{m,γ}(k) =
c_{γ,α}·s·(m)^{1/2-γ}·k`, applies when the monitoring horizon `m'` is
`o(m)`. The stopping time is `τ_{m,γ} = inf{k≥1 : Z_m(k) ≥ g_{m,γ}(k)}`.
The constant `c_{γ,α}` controls size and is analogous to HB's `b_α`. It is
calibrated either asymptotically or with the paper's own finite-sample
approximation, which the authors claim is better than the asymptotic
Extreme Value approximation (we have not checked that claim). The paper
also defines a parallel "Page-CUSUM" variant, designed for a shorter
detection delay.

Numeric results (Table 5.4, median detection delay in periods, no
covariates, `m=200`, three DGP cases): in their "Case I" DGP (`δ_0=0.5`),
the standard weighted CUSUM (`γ=0`) has a median delay of 54 and the
standardised CUSUM (`c_{γ,0.5}`) has a median delay of 37. That is a
reduction of roughly 30% from standardizing, at somewhat lower empirical
power (a rejection frequency of 0.465 against 0.705). The result is a
trade-off between delay and power, and the gain in delay is not free.

Empirical application: online monitoring of Los Angeles daily housing
prices, with training and monitoring windows `m,m'∈{100,200}`. The ex-post
companion analysis dates the actual break at Feb 4, 2009. The real-time
procedure with no covariates first flags a changepoint on Jun 2–15, 2009,
depending on window sizes, which is a delay of about four months. Adding
covariates (interest-rate proxies, VXO, the Weekly Economic Indicator)
brings the flag forward to May 18, 2009 in the richest specification. This
supports the authors' claim that covariates meaningfully shorten detection
delay (we checked it in Table 6.2 of the extracted text).

Cost/feasibility note: the detector statistic `Z_m(k)` itself is cheap and
closed-form, with one `cumsum()`, in the same complexity class as
`monitor_cusum()`. It needs an OLS coefficient `θ̂_m` from the training
window and then a weighted running sum. The larger lift, compared with
`monitor_cusum()` and `monitor_cusum(..., type = "kernel")`, is the
critical-value theory, which is considerably more involved than the single
published constant `b_α` of HB and Astill.

1. There are several boundary-function regimes (open-ended, eq. 2.5;
   closed-ended long-horizon, eq. 2.9; closed-ended short-horizon,
   eq. 2.10), and each needs its own critical value.
2. Two different asymptotic theories apply, depending on `γ`. For `γ <
   1/2` (the "weighted CUSUM"), the critical value solves a
   Brownian-motion sup-norm probability (their eq. 3.4, `P(sup|W(u)| <
   c_γ)`). For `γ = 1/2` (the "standardised CUSUM"), it follows a
   Darling-Erdős-style extreme-value asymptotic (eq. 3.5-3.6:
   `c_{α,0.5} = [x + b(log m)] / a(log m)` with `a(x) = sqrt(2 log x)`,
   `b(x) = 2 log x + 0.5 log log x - 0.5 log π`, and `x` solved from the
   target level through `exp(-exp(-x)) = 1-α`). This asymptotic machinery
   is standard of its kind but new, and exuber has nothing like it today.
3. The paper itself says that the asymptotic critical values in (2) are
   inaccurate ("bound to be inaccurate due to the slow convergence to the
   Extreme Value distribution... leading to low power"). It proposes an
   improved finite-sample correction (eq. 3.7-3.8) that requires solving
   an implicit equation for `c`. That equation has no closed form, so it
   would need numerical root-finding such as `uniroot()`, and it has its
   own tuning parameter `h_m` (the recommended default is `h_m =
   sqrt(log m)`).

None of this is prohibitively hard on its own. `uniroot()` is standard base
R and needs no new dependency. Even so, this is a larger package of new
asymptotic theory to port than anything we have shipped in this bundle so
far. Choosing between the closed-ended and open-ended regimes, and getting
the `γ`-dependent critical value right, need careful and dedicated
attention and are not a same-day extension of existing code. We did not
take it up in this pass.

On 2026-08-10 we re-triaged the item by re-reading rendered pages 7-8 and
11-12 (eq. 2.3-2.9, Theorems 3.1 and 3.3), to check the plan above of
implementing the open-ended `γ = 0` case first. We found one favorable
point and one new blocker.

- Favorable: at `ψ = 0`, the limiting probabilities of Theorem 3.1
  (open-ended) and Theorem 3.3 (closed-ended, short-horizon) both collapse
  to `P{sup_{0<u≤1} |W(u)| < c_{α,0}}`, the classical Brownian-motion
  sup-norm distribution that also underlies the two-sided
  Kolmogorov-Smirnov statistic. It has a well-known closed-form
  alternating series (through the reflection principle) and can be
  inverted for `c_{α,0}` with `uniroot()`. This sub-case needs no
  simulation, as the earlier note had guessed. The general case
  `ψ ∈ (0, 1/2)` does not have this classical form, and the paper says its
  critical values are obtained "by simulation," so only `ψ = 0` is free.
- New blocker: the boundary function's normalizing constant `ð²` (eq.
  2.6, a fraktur-s in the original that every text-extraction attempt
  mangled) is defined by a case split on the model's Lyapunov-type
  exponent `E[log|β₀ + ε_{0,1}|]`. exuber's use case (`H0`: a plain unit
  root, `β₀ = 1`, i.i.d. innovations) falls into the `< 0` branch for any
  reasonable innovation variance. The paper's own Case III DGP, with `β₀ =
  1`, gives `E log|β₀+ε_{0,1}| = -0.007`, which is barely negative and
  "corresponds to the STUR model." That branch requires `ð² = a₁σ₁² +
  a₂σ₂²`, where `a₁ = E[(ȳ₀²/(1+ȳ₀²))²]` and `a₂ = E[(ȳ₀/(1+ȳ₀²))²]` are
  expectations over the stationary distribution `ȳᵢ` of the RCA(1)
  process, which has no general closed form. We checked the paper's
  Section 5 (simulations) for a feasible plug-in estimator that would work
  on real data with unknown parameters and found none. Their Monte Carlo
  section says only that critical values are "computed using Theorem
  3.2/3.6." That suffices for their validation because they know the true
  DGP parameters by construction, since the data are synthetic, and not
  because they give a data-driven estimator of `ð²`. A real implementation
  would either have to derive such an estimator, which would be a
  nontrivial extension of the paper's results and beyond what a
  same-source-only implementation can do, or find one in the companion 2023
  working-paper version, which we did not try.

The net result is that critical-value machinery is a smaller obstacle than
we first scoped for the `ψ = 0` case. The harder requirement is the
nuisance-parameter estimator `ð²`, which the detector needs in order to
run on real data with unknown parameters. We did not resolve it and did
not implement the detector. We record the blocker so that a later pass can
go straight to the feasible `ð²` estimator and will not need to re-derive
the whole problem.

### Kurozumi (2020, 2021), SADF and GSADF cases both implemented (2026-08-10)

Status: both the `SADF` (`s0 = 0`) and `GSADF_{s0}` (`s0 = 0.4`/`0.8`) cases
are implemented as `monitor(..., boundary = "kurozumi", s0 = ...)`. We read
the full PDF of Kurozumi (2020), through Theorem 1 and Table 1, and
rendered it to PNG for accurate transcription. We read Kurozumi (2021) at
the abstract and introduction level only (see its own paragraph below).

One structural finding is useful and has been checked empirically.
Kurozumi (2020)'s `SADF(k) := ADF_1^{m+k}` (his eq., Section 3) is exactly
`radf()`'s `badf` sequence. We confirmed this bit for bit against a
from-scratch OLS ADF t-statistic (fixed start at `t=1`, expanding window
end) at three check points, with tolerance `1e-8`. His `GSADF_{s0}(k) :=
max_{1<=k1<=floor(m*s0)} ADF_{k1}^{m+k}`, in contrast, is not the same
construction as `radf()$bsadf`. The search range for the window start in
exuber's `bsadf` grows with the current monitoring point `t` (from `1` to
`t - minw`), whereas Kurozumi's is capped at a fixed fraction of the
training length `m` regardless of `t`. These are different double
recursions, and the match is not just notational. A fixed, small
window-start band nevertheless needed only a bounded closed-form
computation and no new recursion code (see "`GSADF_{s0}` case, DONE"
below), so both cases are now implemented through `monitor(...,
boundary = "kurozumi", s0 = ...)`. The default `s0 = 0` reuses
`radf()$badf` directly and needs no new statistic. The values `s0 =
0.4`/`0.8` use the new `kurozumi_gsadf_stat()`. In either case a new
published, table-based closed-form threshold replaces `monitor()`'s default
wild-bootstrap-calibrated boundary. This is analogous to the way
`radf_tt_cv()` offers a bootstrap-free alternative to `radf_wb_cv()` for
the static (non-monitoring) GSADF case.

The boundary functions (his eq., checked on the rendered PDF page, since
the raw-text extraction scrambled the subscripts and superscripts) are:

```
SADF:  g_0^df(k/m)    := q_0^df                                    (constant)
GSADF: g_{s0}^df(k/m) := q_{s0}^df * (a_{s0} + b_{s0} * log(c_{s0} + k/m))
       {a,b,c} = {0.76, 0.02, 0.34} for s0=0.4;  {0.73, 0.03, 0.90} for s0=0.8
CS:    g_gamma^cs(k/m) := q_gamma^cs * (1 + k/m)^(1-gamma) * (k/m)^gamma
```

Table 1 gives the scaling constants `q`, by significance level `β` and
monitoring-horizon ratio `s̄ = k̄/m`, where monitoring runs `k̄` observations
past the training length `m`. We transcribed it from the rendered page and
not from the raw OCR text. The table covers only `s̄ ∈ {1, 3, 5}`, which
means a monitoring horizon of 1, 3 or 5 times the training window.

| `s̄` | `β` | `q_0^df` | `q_{0.4}^df` | `q_{0.8}^df` | `q_{0.25}^cs` | `q_{0.45}^cs` |
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

`CS(k)` is HB's own CUSUM statistic (eq. 26) written in Kurozumi's
notation. The `q^cs` columns are therefore published finite-sample
simulated alternatives to `monitor_cusum()`'s asymptotic constant `b_α =
4.6`, for the two cases `γ ∈ {0.25, 0.45}` that the table covers.
Kurozumi obtained them by simulation ("50,000 replications... a standard
Brownian motion is approximated by the sum of suitably normalized i.i.d.
pseudo N(0,1) random variates with increments of 1/1000"), and they do not
come from a closed-form asymptotic result like HB's `b_α = 4.6`. They are
still fixed published numbers, and exuber does not need to simulate
anything new.

What we implemented for `SADF` with `s0 = 0`: `monitor()` gained a
`boundary = c("bootstrap", "kurozumi")` parameter. With `boundary =
"kurozumi"`, it looks up `q_0^df` from Table 1 with `kurozumi_sadf_q(level,
s_bar)`. The lookup snaps `s_bar = (n - T*) / T*` to the nearest tabulated
value in `{1, 3, 5}` and requires `level` to be one of `{0.90, 0.95,
0.99}`. It then compares the value with `radf()$badf` over the monitoring
window. There is no bootstrap and no simulation, and `nboot`, `type`,
`adflag` and `seed` are all ignored on this path. We renamed the statistic
field of the returned list from `bsadf` to `stat`. It holds `bsadf` for
`boundary = "bootstrap"` and `badf` for `boundary = "kurozumi"`.

Validation
(`radf_monitor_kurozumi_boundary_validation.R`):
the table lookups match Table 1 exactly (all 6 checked values, plus
`s_bar` snapping and the error path for an invalid `level`). A basic run
prints without error. The false-alarm rate under `H0` (`n=150`, `T*=75`,
`s_bar=1`, 100 replications) is 4.0% against a nominal 5% target, whereas
the existing wild-bootstrap boundary gives 7.0% on the same DGP. Both are
close to nominal, unlike HB's asymptotic CUSUM bound elsewhere in this
file, which is much more conservative. Detection power under a post-training
bubble (30 replications) is 80% for `kurozumi` and 90% for `bootstrap`.
Both are substantial. We expected some gap because the two boundaries
calibrate different statistics (`badf`/SADF and `bsadf`/GSADF). Alarms
never fire before `T*` (a structural invariant, checked in 10 of 10
replications). The tests added to `test-monitor.R` cover exact table
lookups, a structural run, the invalid-level error, a loose Monte Carlo
bound on the false-alarm rate, and the timing invariant for alarms.

`GSADF_{s0}` case, DONE (2026-08-10). We first scoped this case out
because we thought it needed new recursion code, and we then re-triaged
it. In `GSADF_{s0}(k) := max_{1 <= k1 <= floor(m*s0)} ADF_{k1}^{m+k}`, the
parameter `s0` controls how far the window's start point `k1` may range,
as a fixed fraction of the training length `m`. This range does not grow
as the monitoring point `k` advances. `radf()`'s `bsadf` search range, by
contrast, grows with the current time `t` (from `1` to `t - minw`) by
construction (PSY's GSADF). The two are different double recursions and
not a reindexing of the same one. The re-triage showed that "different
double recursion" does not mean "needs new C++." Since `floor(m*s0)` is a
small, fixed cap (independent of the current monitoring point),
`GSADF_{s0}(k)` needs `ADF_{k1}^{t}` only for `k1` in a small bounded
band, not on a growing triangular grid. Each `ADF_{k1}^{t}` is a plain
with-intercept OLS ADF t-statistic on a fixed window. We can compute it
with the same closed-form cumulative-sum-difference construction that the
project already uses elsewhere (`hls_segment_ssr()` in `dating_hls.R` and
`gls_dfstat_grid()` in `radf_tt.R`). It needs no recursion and no C++,
only `O(1)` work per `(k1, t)` cell from prefix sums, restricted to the
bounded band and not to the full grid.

Implementation: `kurozumi_gsadf_stat()` in `exuber/R/monitor.R` computes
this band with `outer()`-vectorized cumulative-sum differences. It
includes an intercept, unlike the no-intercept GLS-demeaned
`gls_dfstat_grid()`. We verified that it matches `radf()$badf` exactly at
`k1_max = 1`, which confirms that the intercept convention is right.
`monitor(..., boundary = "kurozumi", s0 = 0.4)` or `s0 = 0.8` (the only
two values for which the `a/b/c` scaling constants of his boundary are
tabulated) switches from the flat `SADF` boundary to the `k`-varying
`g_{s0}^df(k/m) := q_{s0}^df * (a_{s0} + b_{s0}*log(c_{s0} + k/m))`
transcribed above. The value `q_{s0}^df` comes from the `q04_df` or
`q08_df` column of Table 1. The default `s0 = 0` reproduces the original
`SADF`-only behavior exactly, which we verified bit for bit, so there is no
regression.

Validated: `kurozumi_gsadf_stat()` matches `radf()$badf` to machine
precision at `k1_max = 1`. It also matches a brute-force `lm()` search over
the restricted window-start band exactly (`|diff| < 1e-14`) at three
monitoring points on a 150-observation series. The `q04_df` and `q08_df`
lookups from Table 1 are exact, with the expected tie-breaking snap (`s0 =
0.6` snaps to the first, lower tabulated value, `0.4`). Alarms never fire
before `T*` (30 of 30 replications). The false-alarm rate under `H0` (300
replications, `n=150`, `T*=75`) is `4.3%`/`5.3%`/`4.3%` for
`SADF`/`GSADF_{0.4}`/`GSADF_{0.8}`, against a nominal `5%`, so all three
are close to nominal. Detection power on a post-training-bubble DGP (60
replications) is `70.0%`/`73.3%`/`66.7%` for `SADF`/`GSADF_{0.4}`/
`GSADF_{0.8}`. The modest edge of `s0 = 0.4` over plain `SADF` echoes
Kurozumi's finding that "GSADF works better than SADF in many cases." The
weaker result for `s0 = 0.8` is specific to this DGP: a wider start range
dilutes power on some alternatives, and it is not an improvement in every
case. We report it as it came out. We extended `test-monitor.R` with 9 new
tests. Replication script:
[replication/monitoring/radf_monitor_gsadf_s0_validation.R](#script-radf_monitor_gsadf_s0_validation).

Kurozumi (2021) ("Asymptotic Behavior of Delay Times of Bubble Monitoring
Tests," *JTSA* 42(3), 314-337) was read at the abstract and introduction
level in this pass. It studies the stochastic order of the detection delay
(the stopping time) for the same detector families as the 2020 paper. It
supports the qualitative split, already cited elsewhere in this file, in
which early or short bubbles favor CUSUM and middle or late bubbles favor
ADF. Our own comparison of `monitor_cusum()` with `monitor()` has also
reproduced that split. The paper adds a layer of dating and inference on
top of detection and is not a new detector, so it has lower priority than
validating the boundary functions of the 2020 paper.

Kurozumi (2020)'s formulas and Table 1 boundary constants are transcribed
in full in the subsection of "Exact numbers/formulas reproduced" above,
together with the structural finding that his `SADF`/`GSADF`/`CS`
detectors are `monitor()`'s and `monitor_cusum()`'s existing statistics.
We do not repeat them here.

### Breitung & Diegel (2025), static LBI test AND sequential extension
both done (2026-08-10)

Status: the static LBI test (`lbi_test()`, for a bubble window known to
span the full sample) and the authors' main sequential extension
(`monitor_lbi()`) are both done. We read the abstract and introduction, the
derivation of the core statistic, and Section 4's sequential extension. We
verified them against rendered PDF pages 3-7, because the raw-text
extraction badly scrambles the `σ̃` and summation notation throughout,
which mattered for both pieces.

The paper proposes a different detector that, by its authors' claim, is
more powerful: a locally best invariant (LBI) statistic. It is robust to
heteroskedasticity by construction, citing the invariance result of
Cavaliere (2005), so it does not require a wild bootstrap the way HB's
plain CUSUM does. Its limiting null distribution is standard normal, and
no simulated or bootstrap critical values are needed. Their eq. 4 gives a
telescoping identity for a bubble known to span the whole sample (`y_0 =
0`, their Section 3): `2*sum(Delta y_t * y_{t-1}) = y_T^2 -
T*sigma_tilde^2`, with `sigma_tilde^2 := T^{-1}*sum(Delta y_t^2)`, the
plain sample variance of first differences under `H0`. Substituting it
into the numerator of the naive DF-type statistic gives their eq. 5,
`LBI_T^2 = y_T^2 / (sigma_tilde^2 * T)`. The whole statistic is the
standardized sample endpoint, with no regression and no recursion:

```
LBI_T = (y_T - y_1) / (sigma_tilde * sqrt(T - 1))
```

We compare it with a standard normal quantile (for example `1.645` at 5%).
The comparison is one-sided because the paper explicitly targets positive
bubbles only, on the grounds that negative bubbles are economically
implausible for a risky asset. This is the cheapest statistic in the whole
project. It needs no C++, no bootstrap, no published table and not even a
boundary function, only a mean and a normal quantile.

Implementation: `lbi_test(data, level = 0.95)` in `exuber/R/lbi_test.R`,
tested in `exuber/tests/testthat/test-lbi.R`. The validation includes a
direct check of the paper's claimed null distribution and not only an
approximately correct size. Eq. 4's telescoping identity holds exactly,
and not only approximately, on a simulated random walk. In a Monte Carlo
under `H0` (500 replications), the empirical mean and standard deviation
of the statistic are `0.023` and `0.976` (theory: `0` and `1`), and the
false-alarm rate at the 95% level is `0.050` against a nominal 5%. The
match is almost exact and not merely close. A Kolmogorov-Smirnov test
against `N(0,1)` gives `p = 0.849`, which is strong evidence that the
statistic follows the claimed distribution exactly and not only
approximately. Detection power under an explosive alternative (60
replications) is `100%`, which matches a standard SADF test on the
identical DGP exactly. Replication script:
[replication/monitoring/radf_lbi_validation.R](#script-radf_lbi_validation).

The sequential extension is now done. The abstract claims that "the
exponentially weighted CUSUM detector with a constant boundary function
turns out to be most powerful." That claim was blocked at first by two
unconfirmed details: the exact construction of the exponential weighting
scheme and the calibrated value of the constant boundary. We pinned both
down by re-reading rendered pages 6-7 (Section 4.1, "Sequential Tests
Based on the LBI Detector").

- The statistic (their eq. 15) applies the classical CUSUM idea directly
  to the LBI construction. We normalize the time index to the monitoring
  period, `r = j/T_m` for `j = 1, ..., T_m`, where `T_m` is the fixed
  monitoring horizon chosen in advance, and compute the (optionally
  weighted) partial sum `LBI_[rT] = (1/(σ̃√T_m)) * sum_{t=1}^{[rT_m]} w_t
  Δy_t ⇒ W(r)`, which is a standard Brownian motion under `H0`. The
  Chu-Stinchcombe-White boundary of `monitor_cusum()` grows as `sqrt(t)`.
  Here we normalize by the fixed `T_m` up front, so a single constant
  boundary controls size uniformly across the whole monitoring window.
  The paper calls this constant-boundary variant `mCUSUM`. It is more
  powerful than the classical time-varying-boundary CUSUM (Brown et al.
  1975, also in their Table 1 but not implemented here, since the paper's
  own recommendation supersedes it) because, under an explosive
  alternative, the detector tends to be largest near the end of the
  monitoring window, and a time-varying boundary with shrinking relative
  tolerance penalizes that.
- The exponential weighting (their eq. 12) is `w_r^c̄ =
  sqrt(2c̄/T_m)/sqrt(e^{2c̄}-1) · e^{c̄r}`, where `c̄ ≥ 0` is a single
  tunable parameter that up-weights later, more bubble-like monitoring
  observations. The value `c̄ = 0` recovers flat weights (`mCUSUM`), and
  `c̄ > 0` gives `wCUSUM`. The authors suggest `c̄ ≈ 2` as a default for a
  moderate power boost.
- The critical value comes from their Table 1 (page 7 of the rendered PDF,
  1,000,000 Monte Carlo replications at `T = 10,000`), which gives
  published one-sided asymptotic critical values. By the paper's
  structural argument, `sup_r{W*(η(r))} =_d sup_r{W(r)}`, the running
  maximum of a time-changed Brownian motion has the same distribution
  whatever the time change, so one set of critical values covers every
  `c̄`, for `mCUSUM` and `wCUSUM` alike: `1.64/1.95/2.24/2.57/2.80` at the
  `10%/5%/2.5%/1%/0.5%` levels. We needed no new simulation, which is the
  published-table-lookup shortcut that already worked for Kurozumi (2020)
  and for HB's FLUC and CUSUM boundaries elsewhere in this project.
- The variance `σ̃²` is estimated from the training window only (the
  paper's Section 4.2: "the training set ... is used for estimating
  nuisance parameters such as `σ²`"). This is consistent with `lbi_test()`,
  which in the static case uses the full-sample `σ̃²` and needs no
  training window.

Implementation: `monitor_lbi(data, r_star = 0.5, c_bar = 0, level = 0.95)`
in `exuber/R/lbi_test.R`, with extended tests in
`exuber/tests/testthat/test-lbi.R`. Validation:

- The sum of squares of the flat-weight (`c̄ = 0`) weight vector equals `1`
  exactly, which is eq. 12's discrete closed form and not an
  approximation. For `c̄ > 0` it equals `1` to within the expected
  Riemann-sum approximation error (`< 0.5%` at `T_m = 500`).
- The final monitoring-point statistic under `mCUSUM` matches, to machine
  precision, a hand-computed telescoped value that uses the
  training-window `σ̃²`. It relies on the same telescoping identity as
  `lbi_test()`, which is now cross-checked in the monitoring context too.
- The table lookups match Table 1 exactly, with a clean error for an
  untabulated level, and alarms never fire before the training window ends
  (50 of 50 replications).
- The false-alarm rate under `H0` (1,000 replications, `n=200`, `T*=100`)
  is `3.7%` for `mCUSUM` and `4.0%` for `wCUSUM`, against a nominal `5%`.
  Both are mildly conservative and not oversized, like most of the other
  finite-sample-against-asymptotic boundaries validated in this project.
- Detection power on a post-training-bubble DGP (60 replications, with the
  bubble starting well into the monitoring window) is `41%` for `mCUSUM`
  and `44%` for `wCUSUM`. Both clearly exceed the `31%` of
  `monitor_cusum(type = "standard")` on the identical DGP, which confirms
  the paper's claim that the constant-boundary LBI detector is more
  powerful than HB's classical CSW-boundary CUSUM. As expected from the
  added up-weighting, `wCUSUM ≥ mCUSUM`.

Replication script:
[replication/monitoring/radf_lbi_monitor_validation.R](#script-radf_lbi_monitor_validation).

Not implemented are two things. One is the row of Table 1 for the classical
time-varying-boundary CUSUM (Brown et al. 1975). The paper uses it as a
point of comparison and not as its contribution, and shows it to be less
powerful than `mCUSUM` and `wCUSUM`. The other is Section 4.2's separate
monitoring variant based on the DF statistic. It is a `badf`-based
analogue with a different boundary, based on the supremum of a ratio. It
addresses the same monitoring problem through a structurally different
statistic that `monitor()` already covers in spirit.

## Cost/feasibility note for exuber

Three things are needed, and exuber had substantial existing
infrastructure for exactly one of them.

1. Family-wise size control through a training-sample critical value, DONE
   (2026-08-09), via `monitor()`. `radf_wb_ps_cv()` in
   `exuber/R/radf_wb.R` already implemented the Phillips & Shi (2020) wild
   bootstrap (its roxygen docs cite "Phillips, P. C., & Shi, S. (2020).
   Real time monitoring of asset markets: Bubbles and crises"). It
   already had a `tb` parameter that truncates the bootstrap DGP to a
   training sub-sample of length `tb`, computes `sadf_crit` and
   `gsadf_crit` on that sub-sample only, and broadcasts those
   training-sample quantiles as a constant boundary across the full
   `bsadf_crit` pointer range. `monitor()` is the orchestration layer that
   this note predicted was missing. It (a) fixes "now" at `T*`, (b) reuses
   the `bsadf` sequence of a single full-sample `radf()` call and does not
   walk forward and refit, which is cheap because the computation is an
   efficient `O(T)` recursion and, as we verified, `bsadf[t]` depends only
   on data up to `t`, (c) compares each monitoring-region value with the
   fixed training boundary, and (d) returns the first breach. See
   "Implementation" above for the full account, including the look-ahead
   leakage pitfall we found and avoided: the null-model fit of
   `radf_wb_ps_cv()` has no internal truncation to `tb`, so the training
   window must be sliced before the call, and `tb` alone cannot enforce
   it. `radf_mc_cv()` could play the same role for the non-bootstrap
   (asymptotic or Monte Carlo) case, but we have not implemented that, and
   `monitor()` supports only the wild-bootstrap route for now.

2. Accounting for the FPR as a function of the monitoring horizon, DONE
   (2026-08-10) for both the `SADF` and `GSADF_{s0}` ADF-family detectors
   through Kurozumi (2020). exuber's `radf_wb_ps_cv(tb=...)` gives a
   single bootstrap critical value. Kurozumi's Table 1 supplies the
   published, simulation-calibrated closed-form analogue for the `SADF`
   (`s0=0`) and `GSADF_{s0}` (`s0=0.4`/`0.8`) cases through `monitor(...,
   boundary = "kurozumi", s0 = ...)`, indexed by the monitoring-horizon
   ratio `s̄ = k̄/m`. AHLST/Whitehouse's eq. (4)–(6) provide the same kind of
   accounting for their detector. Still missing is an equivalent for the
   default wild-bootstrap boundary of `monitor()`, which still recalibrates
   by simulation and not from a closed form.

3. The CUSUM/Page-CUSUM detector family is 0% reusable and all new. HB's
   basic CUSUM and the volatility-robust variant of Astill et al. are both
   DONE (2026-08-10), which confirms this assessment. Items 1, 2, 4 (the
   CUSUM variant), 5, 6 and 7 use a structurally different statistic, a
   standardized running sum of `Δy_t` and not a recursive ADF regression.
   None of the RLS machinery in `exubercore/src` (`rls_gsadf.cpp` and
   `radf.hpp`, the recursive OLS based on the matrix inversion lemma)
   applies, and we confirmed this empirically and did not just predict it:
   `monitor_cusum()` shares no code with `monitor()` or `radf()`. What is
   built and what is still missing:
   (a) The CUSUM partial-sum statistic itself is done as `monitor_cusum()`.
   It is closed-form and needs no C++, for the same reasons as in
   [STADF's closed-form statistic](/replication/volatility-robustness#why-not-exubercore-c).
   (b) The boundary-function and critical-value module is partly done. HB's
   asymptotic `c_t√t` constant (`b_α = 4.6` at the 5% level) is
   implemented. By Corollary 1 of Astill et al., the same boundary works
   unchanged for the volatility-robust variant, and we checked this. HB's
   other statistic, FLUC, is `radf()`'s existing `badf` sequence compared
   with a published finite-sample-simulated boundary (their Table 7,
   transcribed directly), so exuber does not have to simulate the boundary
   itself. It is done too, through `monitor(..., boundary = "fluc")`. HB's
   finite-sample `b_{k,α}` values for CUSUM itself (their Table 8, a
   different table from FLUC's Table 7) are also done (2026-08-10), through
   `monitor_cusum(..., boundary = "finite")`.
   (c) The one-sided kernel spot-variance estimator of the
   volatility-robust variant is done as `monitor_cusum(..., type =
   "kernel")`. It is close in structure to the two-sided profile kernel
   estimator already built for
   [`radf_tt.R`](/replication/volatility-robustness#time-transformed-test-stadf--gstadf)
   but is a separate function. The estimator in `radf_tt.R` estimates a
   cumulative variance profile for the time transformation, using both past
   and future data within a bandwidth window. The new one is deliberately
   one-sided and causal and uses only current and past lags, because
   real-time monitoring at time `t` can see only data up to `t`. This
   confirms the original reading that the engineering pattern could be
   reused but the code could not. The pattern (kernel weights with a
   `filter()`-style rolling computation) carried over, and the code did
   not.

Net assessment (updated 2026-08-09). Item 1, the Phillips & Shi
(2020)-style orchestration of training and monitoring, turned out cheaper
than the original estimate of "plausibly 1-2 weeks." It took a single
focused pass (design, empirical checks of the two reuse assumptions,
about 150 lines of R, tests and docs) and not a multi-week project. The
original estimate priced in uncertainty about whether the existing `tb`
and `bsadf` machinery would compose as predicted, and it did compose
cleanly once we checked. The remaining work is essentially as scoped at
the start. The CUSUM/Page-CUSUM family (item 3) is a structurally
different statistic family and 0% reusable. A real FPR-versus-horizon
boundary function (item 2) needs new asymptotic theory ported in, and not
just new R code. Neither has a shortcut of the kind item 1 had. Building
the CUSUM family well, with its own size-control theory and a coherent
monitoring API, is still not a single-pass item. A realistic estimate is
still several weeks, and it is better described as "the CUSUM/Page-CUSUM
family plus FPR-versus-horizon theory" and no longer as "everything,"
since the one tractable piece is done.
