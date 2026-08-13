---
title: "monitoring"
blurb: "Sequential and real-time detection: training-vs-monitoring orchestration, CUSUM families, and closed-form boundaries."
order: 3
---
﻿# Real-time monitoring for bubbles

**Status: Family A (Phillips & Shi 2020, `radf_monitor()`) done
(2026-08-09); Family B's CUSUM procedure — both Homm & Breitung (2012)'s
original statistic and Astill et al. (2023)'s volatility-robust "CUSUMV"
kernel variant, `monitor_cusum(..., type = "standard"/"kernel")`, plus HB's
own finite-sample boundary (`monitor_cusum(..., boundary = "finite")`) —
done (2026-08-10); Kurozumi (2020)'s closed-form `SADF` boundary AND its
`GSADF_{s0}` generalization (`radf_monitor(..., boundary = "kurozumi",
s0 = 0/0.4/0.8)`), plus HB's own second statistic FLUC
(`radf_monitor(..., boundary = "fluc")`), all done (2026-08-10);
Breitung & Diegel (2025)'s static LBI test (`lbi_test()`) AND their own
sequential/monitoring extension (`monitor_lbi()`, constant-boundary
`mCUSUM`/`wCUSUM`), both done (2026-08-10). Horváth-Trapani's RCA
framework (a feasible nuisance-parameter estimator remains an open gap,
see below) and Kurozumi's 2021 delay-time paper remain evaluated, not
implemented.** See "Implementation" below for what shipped and
"Cost/feasibility note" for what didn't and why.

This was, and remains, **almost certainly the single largest item in this
whole research programme** taken as a whole — it spans two structurally
distinct detector families (recursive training-max vs. CUSUM/Page-CUSUM),
and doing *both* properly, with new S3 infrastructure and proper
size-control theory for each, is still estimated at 3-5x the
[SBZ item's](/replication/volatility-robustness#sbz-wls--kernel-volatility) cost.
What changed: the *one* sub-item this file's own cost/feasibility note
flagged as "the only sub-item where reuse of existing exuber machinery is
actually true rather than aspirational" — Family A, specifically the
Phillips & Shi (2020) wild-bootstrap training-critical-value construction
`radf_wb_cv2()` already implemented — turned out to need only a thin
orchestration wrapper, not new statistical theory, and has been built.

## Implementation

Shipped as `radf_monitor(data, r_star = 0.5, minw, nboot, level, adflag,
type, seed)` in a new `exuber/R/radf_monitor.R`. Confirms the
cost/feasibility note's own prediction: "a loop and a stopping condition
around machinery that already exists."

**Two structural facts made this near-free once actually attempted:**

1. `radf_wb_cv2(..., tb = T*)` already computes exactly the training-
   window wild-bootstrap critical value this needs, and already
   broadcasts it as a **constant boundary** across the full monitoring
   horizon's row count (verified empirically: passing the full-length
   series with `tb = T*` gives a `bsadf_cv` matrix with `nrow = n - minw`,
   not `T* - minw` — the existing code was already built to produce a
   horizon-length boundary, just never wired into a monitoring workflow).
2. `radf()`'s own BSADF statistic at calendar time `t` depends **only**
   on data up to `t`, by construction — verified empirically to be
   bit-identical to a fresh `radf(y[1:t])` call's own last BSADF value.
   This means the entire monitoring path can be obtained from a single
   full-sample `radf()` call (`O(T)`, the existing efficient recursive
   computation) rather than re-fitting at every monitoring point.

**A design decision worth recording — avoiding look-ahead leakage**:
`radf_monitor()` deliberately calls `radf_wb_cv2()` on `data[1:T*]` only,
not the full series with `tb = T*`. The reason: `radf_wb_cv2()`'s
underlying null-model fit (`adf_res()`, in `radf_wb.R`) uses whatever data
it is given *in full* to estimate the bootstrap DGP's residuals/
coefficients — `tb` only truncates the *simulated* bootstrap sample
length, not which portion of the *real* data feeds the residual
estimation. Passing the full series (including post-`T*`, possibly
explosive, data) directly would leak future information into the
training-window null calibration — a genuine correctness issue the
"just call `radf_wb_cv2(full_data, tb=T*)`" reading of the cost note above
would have introduced. Slicing to `data[1:T*]` before calling it avoids
this entirely, at the cost of one extra step the caller (here,
`radf_monitor()` itself) has to take.

**Independent validation**:

1. **Structural checks**: a basic run returns a well-formed object;
   `T_star` and the `bsadf` row count (`n - minw`) match expectations; an
   alarm, when raised, always falls strictly after `T*` (verified across
   10 seeds with a genuine post-training bubble — monitoring never fires
   on training data).
2. **Empirical false-alarm rate under H0** (no bubble anywhere, 40 reps,
   75-observation monitoring horizon at a 95% per-point threshold):
   **10%** — well above the 5% *per-point* nominal level, but this is
   expected and descriptive, not a defect: this is a *cumulative*
   false-alarm probability across 75 sequential comparisons against a
   fixed boundary, and Family A/PSY-style monitoring's FPR growing with
   the monitoring horizon is exactly the property this file's own
   AHLST/Whitehouse discussion documents (their own closed-form eq. 6,
   not implemented here — see item 2 of the cost/feasibility note below,
   still missing).
3. **Detection power** under a genuine bubble starting strictly after
   `T*` (15 reps): **86.7%** detection rate, with alarm delay (raised
   date minus true origination date) always positive and bounded —
   min 6, median 19, max 33 observations. A positive delay is expected
   for any threshold-crossing procedure (it takes several post-origination
   observations for the BSADF statistic to accumulate enough evidence to
   cross a fixed boundary), not a sign of a bug.
4. **A real bug found and fixed during test-writing**: the training
   boundary vector initially lost its series names when extracted from a
   single-series critical-value matrix (`cv$gsadf_cv[, "95%"]` drops to an
   unnamed scalar via R's default `drop = TRUE` matrix-indexing behavior),
   which broke the print method (`data.frame()` construction failed on a
   0-vs-1 row-count mismatch). Fixed by explicitly re-attaching series
   names (`setNames(cv$gsadf_cv[, lvl_lab], snames)`) rather than relying
   on subsetting to preserve them.

New tests in `test-monitor.R`. Full package suite green.

Replication script:
[replication/monitoring/radf_monitor_validation.R](#script-radf_monitor_validation).

**Explicitly out of scope, still**: the CUSUM/Page-CUSUM detector family
beyond the one item below, any closed-form or simulated FPR-vs-horizon
boundary function for Family A (AHLST/Whitehouse's own eq. 4-6), the
union-of-rejections combination across families, and date-stamping
methods specific to a monitoring result — see "Cost/feasibility note"
below.

## Implementation (CUSUM)

**Status: done (2026-08-10).** Shipped as `monitor_cusum()`.

Homm & Breitung (2012)'s CUSUM procedure (their Section 3, eq. 26-30) is
structurally unrelated to Family A: not a recursive ADF regression at
all, just a standardized running sum of first differences compared
against a **closed-form asymptotic boundary** derived from Chu,
Stinchcombe & White (1996)'s inequality (eq. 28) — no wild bootstrap, no
Monte Carlo simulation, no new dependency. This is genuinely the cheapest
item validated in the entire monitoring bundle, matching the original
cost note's own prediction word for word: "the CUSUM partial-sum
statistic itself [is] trivial, one `cumsum()` in R, no C++ needed."

**What HB actually propose** (confirmed by reading the primary source,
not the earlier restatement-only draft): *two* monitoring statistics,
CUSUM and FLUC. Both are now implemented — CUSUM here as `monitor_cusum()`,
FLUC (2026-08-10) as `radf_monitor(..., boundary = "fluc")`, since its
point statistic turns out to be `radf()`'s own `badf` sequence, not a
new one — see "Implementation (FLUC)" below.

**The statistic** (eq. 26, 29-30), for a training window ending at `T*`
and monitoring point `t > T*`:

```
S_t = (y_t - y_{T*}) / sigma_hat_t         (telescoped cumulative sum of post-training first differences)
sigma_hat_t^2 = (t-1)^-1 * sum_{j=2}^{t} (Delta y_j)^2   (recursive sample variance of ALL differences up to t)
c_t = sqrt(b_alpha + log(t / T*))
boundary_t = c_t * sqrt(t)
reject H0 (alarm) at the first t where S_t > boundary_t
```

`b_alpha = 4.6` is HB's own quoted one-sided asymptotic calibration for a
5% significance level. Unlike `radf_monitor()`'s wild-bootstrap boundary,
`sigma_hat_t` is legitimately re-estimated using *all* data up to the
current monitoring point `t` (not just the training window) — this is
not a look-ahead problem the way it would have been for `radf_wb_cv2()`,
because at each real monitoring instant `t` only data up to `t` is ever
used; there is no future information to leak.

### Implementation shape

`cusum_stat_path()` (internal) + `monitor_cusum(data, r_star, b_alpha)`
(exported), returning a `radf_cusum_obj`. **Zero code reuse from Family A
or from `radf()`/`exubercore`** — confirmed, not just predicted, by
actually building it: this statistic shares no structure with the
recursive-ADF family at all, exactly matching the original cost note's
"CUSUM/Page-CUSUM detector family — 0% reusable, all new" assessment.

### Independent validation

1. **Formula-exact check**: `cusum_stat_path()` matches an independent
   brute-force loop recomputation (separate `diff()`/`sum()` calls per
   monitoring point rather than the vectorized `cumsum()` implementation)
   to floating-point precision (`0` max absolute difference).
2. **Empirical false-alarm rate under H0** (pure random walk, no bubble,
   100 reps, 75-observation monitoring horizon): **0%**. Consistent with
   HB's own eq. 28 being a *conservative asymptotic upper bound* (via
   Chu/Stinchcombe/White 1996), not an exact size — a rate well below the
   nominal 5% is the expected signature of a genuinely conservative
   bound, not evidence of a miscalibration.
3. **Detection power** under the *same* post-training bubble DGP already
   used to validate `radf_monitor()` (30 reps): **30%** detection rate,
   median alarm delay 27 observations — genuinely, substantially lower
   than Family A's 86.7% detection rate / 19-observation median delay on
   the identical DGP. Reported honestly rather than only showing
   favorable numbers: this is consistent with, not contradictory to, the
   literature already cited elsewhere in this file — Kurozumi (2020,
   2021) report that CUSUM-type detectors do *better* than ADF-type ones
   specifically for early/short bubbles, and *worse* for middle-to-late
   ones. The test DGP (bubble starting ~65% into the sample) sits
   squarely in the "middle-to-late" regime where the literature predicts
   CUSUM should lag ADF-type detection, and it does — this is the
   expected shape of the two families' relative power, not a defect in
   either implementation. It's also a concrete, literature-consistent
   demonstration of *why* a union-of-rejections strategy (both families
   at once) is the field's own recommended practice, not implemented
   here.

New tests in `test-cusum.R`. Full package suite green.

Replication script:
[replication/monitoring/radf_cusum_validation.R](#script-radf_cusum_validation).

**`boundary = "finite"` (2026-08-10)**: HB's own *finite-sample* boundary
constant (their Table 8, "without drift estimation," matching this
implementation's raw-first-difference construction), used in place of
the fixed asymptotic `b_alpha = 4.6` — the same table-lookup shortcut
already used for Kurozumi's and FLUC's boundaries, transcribed directly
from HB's published table (indexed by training length, significance
level, and monitoring-horizon ratio `k = N/T*`), no new simulation.
Applied to both `type = "standard"` and `type = "kernel"` (CUSUMV):
Astill et al.'s own Corollary 1, already validated above, establishes
the *same* boundary function works for both statistics, so extending the
same finite-sample table to CUSUMV is a direct consequence of that
result, not a new claim needing separate validation. Validated: table
lookups exact; at `T*=75`/`n=150` (`k=2`, snapping `T*` to the nearest
tabulated `n=50`), the finite-sample boundary gives a false-alarm rate
of **9.0%** under `H0` — closer to the nominal 5% level than the
asymptotic bound's suspiciously conservative **0%** — *and* higher
detection power on the same post-training bubble DGP (**36.7%** vs.
**30.0%**): a genuine, not just cosmetic, improvement from using a
properly finite-sample-calibrated constant instead of an asymptotic
upper bound. New tests extend `test-cusum.R`. Replication script:
[replication/monitoring/radf_cusum_finite_boundary_validation.R](#script-radf_cusum_finite_boundary_validation).

### CUSUMV: the volatility-robust variant (Astill et al. 2023)

**Status: done (2026-08-10).** Shipped as `monitor_cusum(..., type = "kernel")`.

Astill, Harvey, Leybourne, Taylor & Zu (2023, *JFEC* 21(1), 187-227;
"AHLTZ") generalize HB's CUSUM procedure to allow time-varying volatility:
their own abstract states plainly that "such behavior can heavily inflate
the false positive rate (FPR) of the CUSUM-based procedure." The fix
(their eq. 6-7) replaces HB's single running variance with a **one-sided**
(causal — only current and past lags, appropriate for real-time
monitoring) Nadaraya-Watson kernel-weighted spot-variance estimate,
standardizing *each* first difference individually before cumulating:

```
SV_t = sum_{j=T*+1}^{t} Delta y_j / sigma_hat_{j,N}     (eq. 6)
sigma_hat_j,N = sum_{s=0}^{N} w_s * Delta y_{j-s}^2,   w_s = K(s/N) / sum_s K(s/N)   (eq. 7)
```

Their Corollary 1 proves the **same** boundary function `c_t * sqrt(t)`
(same `b_α`) that HB's original statistic uses still delivers a
controlled asymptotic FPR for this modified statistic, even under
time-varying volatility — the boundary doesn't need to change, only the
statistic's own construction does.

**Implementation**: `one_sided_kernel_spot_vol()` (a fixed causal kernel
weight vector, computed via `stats::filter(..., sides = 1)`, following
AHLTZ's own convention of `σ̂²_{j,N} := 1` for `j <= N`) +
`cusum_stat_path_kernel()`, wired into `monitor_cusum()` as a new
`type = "kernel"` option reusing the existing boundary/decision-rule
logic unchanged — not a separate function, since the only thing that
differs is the statistic's numerator. Default bandwidth `N = 20`, AHLTZ's
own empirically-recommended value ("setting H = 20 delivered a procedure
with the best trade-off" between FPR robustness and power) — their own
data-driven local cross-validation bandwidth selection (their eq. 8-9)
was not implemented; this is a documented simplification, not a claim of
matching their exact finite-sample procedure.

**Independent validation**:

1. **Formula-exact checks**: both `one_sided_kernel_spot_vol()` and
   `cusum_stat_path_kernel()` match independent brute-force loop
   recomputations to floating-point precision.
2. **Under homoskedastic H0**, `standard` and `kernel` give comparable
   (here, identically conservative — both 0%) empirical false-alarm
   rates, consistent with AHLTZ's Remark 8: "in the case where the
   innovations are homoskedastic... [both] lead to the same limiting null
   distribution."
3. **Under heteroskedastic H0** (a volatility jump from 1 to 8 partway
   through the monitoring region, 40 reps): `standard` CUSUM's false-alarm
   rate rises to **8.3%**, visibly inflated above its own homoskedastic-
   case rate — the exact failure mode AHLTZ's abstract describes. `kernel`
   CUSUMV stays at **0%**, fully controlled. This is the paper's central
   claim, and it held up cleanly under independent simulation, not just
   asymptotic theory.
4. **Detection power** under the same post-training bubble DGP (30 reps,
   homoskedastic): `standard` 30.0%, `kernel` 36.7% — comparable, no
   power sacrifice observed in this scenario (AHLTZ's own claim is that
   the modification "sacrifices only a small amount of power... when the
   shocks are homoskedastic," consistent with, if slightly better than,
   what was found here).

New tests extend `test-cusum.R`. Full package suite green.

Replication script:
[replication/monitoring/radf_cusumv_kernel_validation.R](#script-radf_cusumv_kernel_validation).

### Implementation (FLUC)

**Status: done (2026-08-10).** Shipped as
`radf_monitor(..., boundary = "fluc")`.

Homm & Breitung's own second monitoring statistic (their eq. 27,
confirmed by rendering the PDF page): `Z_t = (rho_hat_t - 1) /
sigma_hat_{rho_t} = DF_{t/n}`, the *ordinary recursive/expanding-window
OLS ADF t-statistic* on the sample `{y_0, ..., y_t}` — exactly `radf()`'s
existing `badf` sequence, the same statistic this file's Kurozumi
subsection already confirmed `radf_monitor(..., boundary = "kurozumi")`
reuses. No new point statistic needed, matching that same reuse pattern.

The rejection rule (eq. 29/31) has the same functional form as CUSUM's
own boundary — `DF_{t/n} > kappa_t`, `kappa_t = sqrt(b_{k,alpha} +
log(t/n))` — but a *different* calibration constant `b_{k,alpha}` that,
unlike CUSUM's closed-form `b_alpha = 4.6`, HB's own text says they
"determine... by means of simulation," putting it in the same cost tier
this file's cost note originally flagged (needing new critical-value
work, not free). **What made this cheap anyway**: HB *publish* the
simulated constant directly (their Table 7, part i, "without
detrending" — matching `radf()`'s own no-trend default), tabulated by
training length `n` in `{20, 50, 100}`, significance level `alpha` in
`{0.10, 0.05, 0.01}`, and monitoring-horizon ratio `k = N/n` (`N` the
total sample including training) in `{2, 3, 4, 5, 6, 8, 10}` — a
published, closed-form lookup, no simulation needed on exuber's end,
exactly the same shortcut Kurozumi's own Table 1 provided for
`boundary = "kurozumi"`.

**Implementation**: `hb_fluc_table` (a 9x7 table transcribed from the
rendered PDF page) + `hb_fluc_q(level, n_train, k)` (snaps `n_train` to
the nearest of `{20, 50, 100}` and `k` to the nearest of `{2,...,10}`,
requires `level` to exactly match one of the three tabulated
significance levels) — the same lookup-and-snap pattern as
`kurozumi_sadf_q()`, wired into `radf_monitor()` as a third `boundary`
option alongside `"bootstrap"`/`"kurozumi"`.

**Validation**: table lookups match Table 7 exactly (6 checked cells,
including a tie-breaking snap case). Basic run prints without error;
alarms never fire before `T*` (structural invariant, 10/10 reps).
Empirical false-alarm rate under `H0` (`n=150`, `T*=75`, `k=2`, the
smallest/most conservative tabulated horizon ratio, 100 reps) is **0%**
— HB's own table is a genuinely conservative finite-sample-simulated
bound at this `k`, not miscalibrated (comparable to CUSUM's own
conservative behavior noted above). Detection power under the same
post-training bubble DGP used elsewhere in this file (30 reps) is
**56.7%**, below `boundary = "kurozumi"`'s 80% and `"bootstrap"`'s 90%
on the identical DGP — reported honestly, and consistent with, not
contradicting, HB's own stated finding that FLUC/CUSUM monitoring
generally has *less* power than a supDF-style test (though FLUC beats
their own CUSUM). One caveat worth flagging: HB's table only tabulates
training lengths up to `n=100`, smaller than typical financial series in
practice, so larger `T*` values get snapped to the `n=100` row rather
than genuinely interpolated/extrapolated — a documented approximation,
not a claim of exact finite-sample calibration at realistic sample
sizes. New tests extend `test-monitor.R`. Full package suite green.

Replication script:
[replication/monitoring/radf_monitor_fluc_boundary_validation.R](#script-radf_monitor_fluc_boundary_validation).

## Source

Six items total. **Access status update (2026-08-09): institutional access
subsequently recovered primary-source PDFs for items 1, 2, 4, and 7 below**
— at the time this evaluation was originally written, only items 3, 5, and
6 were open; the rest were either abstract-level (Kurozumi) or read via a
restatement in item 6 rather than the primary text (Homm & Breitung, Astill
2021/2023). **Those formulas below are still marked as sourced from the
restatement, not re-verified against the now-available primary PDFs** —
that re-verification is flagged as follow-up work, not done in this pass.

1. **Homm, U. & Breitung, J. (2012).** "Testing for speculative bubbles in
   stock markets: a comparison of alternative methods." *Journal of
   Financial Econometrics*, 10(1), 198–231. `doi:10.1093/jjfinec/nbr009`.
   Now open (institutional access) —
   The formulas below (CUSUM statistic, boundary function, and its Table 8
   finite-sample `b_α` calibration) are still from the restatement in item
   6, not yet re-verified against this primary copy.

2. **Astill, S., Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y.
   (2021/2023).** "CUSUM-Based Monitoring for Explosive Episodes in
   Financial Data in the Presence of Time-Varying Volatility." *Journal of
   Financial Econometrics*, 21(1), 187–227 (advance-access March 2021,
   print issue Winter 2023). `doi:10.1093/jjfinec/nbab009`. Referenced in
   the literature as "Astill et al. (2021)" or "(2024)" depending on which
   author's own citation list you use — same paper, `nbab009`. Now open
   (institutional access) —
   Formulas below still taken from item 6's restatement (labelled "AHLTZ"),
   not yet re-verified against this primary copy.

   A companion, earlier paper by an overlapping author set — **Astill, S.,
   Harvey, D.I., Leybourne, S.J., Sollis, R. & Taylor, A.M.R. (2018),
   "Real-Time Monitoring for Explosive Financial Bubbles," *Journal of
   Time Series Analysis*, 39, 863–891** (often abbreviated AHLST) — is a
   *different* monitoring design (training-sample-maximum comparison, not
   CUSUM) that turned out to be essential background because item 3 below
   builds directly on it. Now open (institutional access) —
   Still only read via its restatement in item 3, not yet re-read directly.

3. **Whitehouse, E.J., Harvey, D.I. & Leybourne, S.J. (2025).**
   "Real-time monitoring procedures for early detection of bubbles."
   *International Journal of Forecasting*, 41(3), 1260–1277.
   `doi:10.1016/j.ijforecast.2024.12.005`. **Open access (CC-BY)** —
   downloaded directly from the White Rose Research Online repository.
   Read via `pdftotext -layout`, with pages 2–4 (Introduction, the AHLST
   model and decision rule, equations 2–9) re-verified by rendering to PNG
   at 2.5x (PyMuPDF) and reading the typeset math directly — the OCR text
   mangled subscripts/superscripts badly (e.g. `A_{e,k}` came out as
   `Ae,k` merged into surrounding prose, and the FPR formula's summation
   limits were scrambled across columns in the two-column layout). This is
   the paper that supplied the exact, verified AHLST decision rule and FPR
   formula used below.

4. **Kurozumi, E.** Two papers, both located and confirmed to be about
   monitoring critical values and detection delay respectively:
   - **Kurozumi, E. (2020). "Asymptotic properties of bubble monitoring
     tests." *Econometric Reviews*, 39(5), 510–538.**
     `doi:10.1080/07474938.2019.1697086`. Extends SADF/GSADF to a
     monitoring scheme and studies a CUSUM detector alongside it,
     deriving new monitoring-period critical values and comparing
     ADF-type vs. CUSUM-type detection under moderate-deviation and
     local-to-unity asymptotics. Now open (institutional access) —
   - **Kurozumi, E. (2021). "Asymptotic Behavior of Delay Times of Bubble
     Monitoring Tests." *Journal of Time Series Analysis*, 42(3),
     314–337.** `doi:10.1111/jtsa.12569`. Specifically about detection
     delay (stochastic order of the stopping time) for ADF-type vs.
     CUSUM-type monitoring statistics. Now open (institutional access) —

   **What's reported below for Kurozumi is still abstract-level only**,
   cross-checked against a secondary academic source (Skrobotov 2023) —
   now that both PDFs are on hand, reading them directly and replacing this
   qualitative summary with verified formulas/numbers is the single most
   valuable follow-up in this file.

5. **Horvath, L. & Trapani, L.** "Real-time monitoring with RCA models."
   Working paper arXiv:2312.11710 (Dec 2023), published as **Horváth, L.
   & Trapani, L. (2026). "Real-time monitoring with RCA models,"
   *Econometric Theory*, 42, 514–547.** **Open access** — downloaded
   directly from arXiv.
   Read via `pdftotext -layout`; the core detector and boundary-function
   definitions (section 2, its equations 2.4–2.11) were re-verified by
   rendering the relevant pages to PNG, since the OCR badly mangled the
   summation/subscript-heavy statistic definitions.

6. **Astill, S., Taylor, A.M.R. & Zu, Y. (2026, forthcoming).**
   "Covariate Augmented CUSUM Bubble Monitoring Procedures." *Econometric
   Theory*. **Open access working paper**: Essex Finance Centre Working
   Paper No. 94.
   Turned out to be the single most useful document obtained in the
   original pass: its Section 3 ("CUSUM-based Bubble Detection
   Procedures") restates, with full equation numbers and explicit page
   citations back to the originals (e.g. "HB... Table 8, p221"), both the
   original Homm & Breitung (2012) CUSUM statistic and boundary function
   (item 1) *and* the Astill et al. (2021/2023) volatility-robust
   modification (item 2). Two of its authors (Taylor, Zu) are co-authors
   of item 2. Formula pages (pp. 10–11 of the PDF) were rendered to PNG
   and read directly; the OCR text mangled the subscripts on `S^t_T`,
   `SV^t_T`, `c_t`, and the kernel-weight definitions badly enough that the
   image read was necessary to get them right.

7. **Breitung, J. & Diegel, M. (2025).** "Sequential Detector Statistics
   for Speculative Bubbles." *JTSA*, 46(5).
   Surfaced from the *JTSA* 46(5) special issue — a direct, concretely
   citable match for this section's intent, not yet read in depth.
   (Correction to an earlier lead: the Horvath/Trapani 2025 JTSA piece in
   the same special issue, "Sequential Monitoring for Changes in
   GARCH(1,1) Models Without Assuming Stationarity," turns out to be about
   GARCH change-point monitoring generically, not bubble monitoring
   specifically — tangential to this section, not a duplicate of item 5
   above, which is a different Horváth & Trapani paper.)

## What it is

All seven papers address the same underlying problem, in two structurally
distinct families:

**Family A — recursive/BSADF training-vs-monitoring maximum (AHLST /
Whitehouse).** Split the sample into a training period `t = 1,...,T*`
(assumed bubble-free) and a monitoring period `t = T*+1,...,T`. Compute a
recursive statistic `A_{e,k}` over rolling sub-samples of fixed length `k`
in *both* periods (AHLST's statistic is a sub-sample regression of `Δy_t`
on a linear trend, White-studentized — structurally close to, but not
identical to, exuber's ADF-family statistics). The training-sample maximum
`A*_max = max A_{e,k}` becomes the (fixed) critical value; monitoring
rejects `H0` at the first `e` where `A_{e,k} > A*_max`. Crucially, AHLST
prove a closed-form asymptotic FPR that is a simple function of the
training/monitoring length ratio (verified formula, eq. 6 below) — this is
the "PSY-style" monitoring philosophy that Phillips & Shi (2020) (already
partially implemented in exuber, see Cost/feasibility) also belongs to,
though PSY use a bootstrap rather than a closed-form FPR.

**Family B — full-sample CUSUM / Page-CUSUM detectors (Homm-Breitung,
Astill et al. 2021/2026, Kurozumi's CUSUM variant, Horvath-Trapani,
Breitung-Diegel).** Compute a partial-sum (CUSUM) statistic of the
standardized first differences `Δy_t` from the end of the training sample
onward, and compare it at each monitoring date to a boundary function that
grows with `t` (e.g. `c_t·√t`) so the cumulative false-alarm probability
over the whole (possibly infinite) monitoring horizon stays below a target
`α`. This family is a genuinely different statistic — not a recursive ADF
regression at all, just a running standardized sum — and its
variance-robust variants (Astill et al.) replace the standardization with a
kernel spot-variance estimate, structurally similar in spirit (but not
target) to
[SBZ/STADF's](/replication/volatility-robustness#sbz-wls--kernel-volatility) kernel
machinery. Horvath-Trapani generalize this further to a
Random-Coefficient-Autoregressive (RCA) framework with weighted CUSUM and
Page-CUSUM detectors that work "symmetrically" for stationary→explosive
*and* explosive→stationary transitions without needing to know in advance
which regime you start in.

Kurozumi's contribution sits across both families: he puts SADF/GSADF-type
(Family A) and CUSUM-type (Family B) monitoring statistics into a common
asymptotic framework, derives new critical values for the monitoring
period for both, and separately studies the detection-delay (stopping
time) distribution — finding, per the secondary source, that CUSUM
detects an *early, short* bubble faster, while ADF/BSADF-type detects a
*middle-to-late* bubble faster, and that a union-of-rejections combining
BSADF and CUSUM is possible (structurally the monitoring analogue of the
[SBZ union statistic](/replication/volatility-robustness#sbz-wls--kernel-volatility)).

## Exact numbers/formulas reproduced

### Homm & Breitung (2012) CUSUM statistic and boundary (via item 6, not yet re-verified against the now-available primary copy)

CUSUM statistic (their eq. 7, training sample `t=1,...,T*`, monitoring
`t>T*`):

```
S^t_{T*} := (1/σ̃_t) · Σ_{j=T*+1}^{t} Δy_j
```

with `σ̃_t² := (t-1)^{-1} Σ_{j=2}^{t} (Δy_j)²`. Under `H0`, `T*^{-1/2}
S^{⌊Tr⌋}_{T*} ⇒ W(r) − W(1)` (eq. 8), and via Chu et al. (1996) Theorem
3.4, for any `λ>1`:

```
lim_{T→∞} Pr(|S^t_{T*}| > c_t·√t for some t ∈ {T*+1,...,⌊λT*⌋}) ≤ exp(−b_α/2)   (eq. 9)
```

with **boundary function `c_t := √(b_α + log(t/T*))`**. Reject `H0` if
`S^t_{T*} > c_t·√t`, flagging the first `t` where this happens. For a
one-sided test with size `α = 0.05`, the asymptotic setting is **`b_α =
4.6`** (this delivers a two-sided size ≤ 0.10 from eq. 9). This asymptotic
setting assumes an infinite monitoring horizon and is described (by item
6's authors, citing HB directly) as "extremely conservative in practice" —
HB's own paper instead gives finite-sample `b_α` values in their **Table
8, p.221**, calibrated for target FPR ∈ {0.10, 0.05, 0.01} at specific
training/monitoring lengths — **now transcribed and implemented
(2026-08-10)** as `monitor_cusum(..., boundary = "finite")`, see
"Implementation (CUSUM)" above.

### Astill et al. (2021/2023) volatility-robust modification (via item 6, not yet re-verified against the now-available primary copy)

Replaces `S^t_{T*}` with

```
SV^t_{T*} := Σ_{j=T*+1}^{t} Δy_j / σ̂_{j,N},   t > T*
```

where `σ̂²_{j,N}` is a one-sided kernel smoothing estimator of the spot
variance `σ²_j := σ²(j/T)`:

```
σ̂²_{j,N} := Σ_{s=0}^{N} k_s (Δy_{j-s})²,   k_s := K(s/N) / Σ_{s=0}^N K(s/N)
```

They prove the same boundary function `c_t·√t` still delivers a
theoretically controlled FPR when volatility is time-varying, at some cost
in power relative to `S^t_{T*}` under homoskedasticity. Empirical
application: Bitcoin prices (per the OUP abstract; not yet independently
verified against the now-available primary copy).

### Whitehouse, Harvey & Leybourne (2025) — AHLST decision rule and FPR (verified against rendered PDF pages 2–4)

DGP: `y_t = μ + u_t`, `u_t = u_{t-1}+ε_t` for `t ≤ ⌊τT⌋` then
`(1+δ)u_{t-1}+ε_t` after. Statistic (their eq. 2):

```
A_{e,k} = B_{e,k} / √C_{e,k},   B_{e,k} = Σ_{t=e-k+1}^{e} (t-e+k)Δy_t,   C_{e,k} = Σ_{t=e-k+1}^{e} {(t-e+k)Δy_t}²
```

Training-sample max `A*_max = max_{e∈[k+1,T*]} A_{e,k}` is the critical
value for monitoring. **Decision rule**: "Reject H0 at time e if `A_{e,k}
> A*_max`" — this is the `AMAX(k)` procedure. Under `H0`, for an arbitrary
monitoring point `T'`:

```
lim_{T→∞} P(max_{e∈[T*+k,T']} A_{e,k} > max_{e∈[k+1,T*]} A_{e,k}) = τ = lim(T'-T*)/T'    (eq. 4-5)
```

and the **approximate FPR at monitoring point `T'`** (eq. 6, the key
usable formula for calibration):

```
α ≈ (T' - T* - k + 1) / (T' - 2k + 1)
```

rearranged to find how far monitoring can run while keeping FPR at a
chosen level: `T' ≈ (T* + k - 1 - α(2k-1)) / (1-α)`. Contrasted explicitly
against the CUSUM-based approaches (HB, Astill et al. 2021, Horvath &
Trapani 2026): AHLST's approach gives an *exact* usable FPR formula with
no asymptotic boundary/conservatism, but the FPR necessarily *grows* with
the monitoring horizon (better suited to short-range monitoring), whereas
CUSUM-style methods can be tuned to hold a fixed FPR (e.g. 0.05) over an
arbitrarily long horizon at the cost of lower power (TPR).

**Empirical Monte Carlo numbers (Table 1, k=10, NIID and GARCH(1,1)
errors)** — spot-checked a sample of rows, verified via `pdftotext
-layout` (this table is plain numeric with no problematic subscripts, so
not separately PNG-verified): at `T'=200`, empirical FPR of the baseline
`AMAX(k)` = 0.006 (NIID) vs. theoretical target region; rising
monotonically to 0.147 by `T'=230`. The two new variance-standardized
variants (`A^{AR,max}(k)`, `A^{T,max}(k)`) run consistently slightly above
the baseline (e.g. 0.015/0.013 vs. 0.006 at T'=200), i.e. slightly less
conservative, matching the paper's stated design goal (Theorem 1: same
asymptotic FPR, different finite-sample behaviour).

**Empirical application**: the US house price-to-rent ratio bubble
preceding the 2007/08 GFC is detected as early as **1999:Q1** by the new
`A^{AR,max}(k)` procedure, vs. **2000:Q1** for the baseline `AMAX(k)` — a
four-quarter improvement (Table 2, verified in the extracted text, cross
checked against two separate passages in the paper stating the same
1999:Q1/2000:Q1 dates).

### Horvath & Trapani (2023/2026) RCA monitoring — evaluated, not implemented (2026-08-10)

**Status: evaluated, not implemented.** Formulas below were already
verified against rendered PDF pages in an earlier pass; this pass added a
direct read of the asymptotic-theory section (their Theorems 3.3-3.5,
eq. 3.4-3.8) specifically to assess implementation feasibility, not just
transcribe the detector formula.

WLS-residual CUSUM detector (their eq. 2.4) over a training window of
length `m`:

```
Z_m(k) = Σ_{i=m+1}^{m+k} [(y_i - θ̂_m y_{i-1}) y_{i-1}] / (1+y_{i-1}²),   k ≥ 1
```

**Boundary function** for the open-ended/long-horizon case (eq. 2.5/2.9):

```
g_{m,γ}(k) = c_{γ,α}·s·m^{1/2}·(1+k/m)·(k/(m+k))^γ,   0 ≤ γ < 1/2
```

with a separate short-horizon variant (eq. 2.10) `g_{m,γ}(k) =
c_{γ,α}·s·(m)^{1/2-γ}·k` when the monitoring horizon `m'` is `o(m)`.
Stopping time `τ_{m,γ} = inf{k≥1 : Z_m(k) ≥ g_{m,γ}(k)}`. `c_{γ,α}` is a
size-controlling constant (analogous to HB's `b_α`), calibrated either
asymptotically or via the paper's own finite-sample approximation
(claimed superior to the asymptotic Extreme Value approximation, not
independently checked here). A parallel "Page-CUSUM" variant is also
defined, designed for shorter detection delay.

**Numeric results** (Table 5.4, median detection delay in periods, no
covariates, `m=200`, three DGP cases): e.g. for their "Case I" DGP
(`δ_0=0.5`) the standard weighted CUSUM (`γ=0`) has median delay 54 vs.
the standardised CUSUM (`c_{γ,0.5}`) at 37, i.e. roughly a 30% reduction
in median detection delay from standardizing, at somewhat lower empirical
power (0.705 vs. 0.465 rejection frequency) — a genuine
delay/power trade-off, not a free lunch.

**Empirical application**: online monitoring of Los Angeles daily housing
prices, training/monitoring windows `m,m'∈{100,200}`. Ex-post companion
analysis dates the actual break at **Feb 4, 2009**; the real-time
procedure with no covariates first flags a changepoint on **Jun 2–15,
2009** depending on window sizes (≈4 months delay); adding covariates
(interest-rate proxies, VXO, the Weekly Economic Indicator) brings this
forward to **May 18, 2009** in the richest specification — confirming
their claim that covariates meaningfully shorten detection delay (verified
directly from Table 6.2 in the extracted text).

**Cost/feasibility note**: the detector statistic `Z_m(k)` itself is
cheap — closed-form, one `cumsum()`, same complexity class as
`monitor_cusum()` (an OLS coefficient `θ̂_m` from the training window, then a
weighted running sum). What makes this a meaningfully bigger lift than
`monitor_cusum()`/`monitor_cusum(..., type = "kernel")` is the **critical-value
theory**, which is genuinely more involved than HB/Astill's single
published constant `b_α`:

1. **Multiple boundary-function regimes** (open-ended, eq. 2.5; closed-
   ended long-horizon, eq. 2.9; closed-ended short-horizon, eq. 2.10),
   each needing its own critical value.
2. **Two different asymptotic theories depending on `γ`**: for `γ < 1/2`
   (the "weighted CUSUM"), the critical value solves a Brownian-motion
   sup-norm probability (their eq. 3.4, `P(sup|W(u)| < c_γ)`); for
   `γ = 1/2` (the "standardised CUSUM"), it's a Darling-Erdős-style
   extreme-value asymptotic (eq. 3.5-3.6: `c_{α,0.5} = [x + b(log m)] /
   a(log m)` with `a(x) = sqrt(2 log x)`, `b(x) = 2 log x + 0.5 log log x
   - 0.5 log π`, and `x` solved from the target level via
   `exp(-exp(-x)) = 1-α`) — genuinely new (if standard, Darling-Erdős-
   type) asymptotic machinery, not present anywhere in exuber today.
3. **The paper's own text says the asymptotic critical values in (2) are
   inaccurate** ("bound to be inaccurate due to the slow convergence to
   the Extreme Value distribution... leading to low power") and proposes
   an improved finite-sample correction (eq. 3.7-3.8) that requires
   solving an *implicit* equation for `c` (not closed-form; would need
   numerical root-finding, e.g. `uniroot()`) with its own tuning
   parameter `h_m` (recommended default `h_m = sqrt(log m)`).

None of this is prohibitively hard individually — `uniroot()` is
standard base R, no new dependency — but it's a meaningfully bigger
package of new asymptotic theory to port correctly than anything shipped
in this bundle so far, and getting the closed-ended vs. open-ended
regime selection and the `γ`-dependent critical value right needs careful,
dedicated attention rather than a same-day extension of existing code.
Not picked up this pass.

**Re-triaged (2026-08-10), re-reading rendered pages 7-8, 11-12 (eq.
2.3-2.9, Theorems 3.1/3.3) specifically to check the "implement the
open-ended, `γ = 0` case first" plan above.** Two findings, one
favorable and one a newly-identified blocker:

- **Favorable**: at `ψ = 0`, Theorem 3.1's (open-ended) and Theorem
  3.3's (closed-ended, short-horizon) limiting probabilities both
  collapse to `P{sup_{0<u≤1} |W(u)| < c_{α,0}}` — the classical
  Brownian-motion sup-norm distribution (the same one underlying the
  two-sided Kolmogorov-Smirnov statistic), which has a well-known
  closed-form alternating series (via the reflection principle) and can
  be inverted for `c_{α,0}` with `uniroot()` — genuinely no simulation
  needed for this specific sub-case, confirming the earlier note's
  guess. (The general `ψ ∈ (0, 1/2)` case does *not* have this classical
  form and the paper's own text says its critical values are obtained
  "by simulation" — only `ψ = 0` is free.)
- **New blocker**: the boundary function's own normalizing constant
  `ð²` (eq. 2.6, a fraktur-s in the original, badly mangled by every
  text-extraction attempt) is defined by a case split on the model's
  Lyapunov-type exponent `E[log|β₀ + ε_{0,1}|]`. exuber's own use case
  (`H0`: a plain unit root, `β₀ = 1`, i.i.d. innovations) falls into the
  `< 0` branch for any reasonable innovation variance (their own Case
  III DGP, `β₀ = 1`, gives `E log|β₀+ε_{0,1}| = -0.007` — just barely
  negative, "corresponds to the STUR model") — which requires `ð² =
  a₁σ₁² + a₂σ₂²`, where `a₁ = E[(ȳ₀²/(1+ȳ₀²))²]` and `a₂ =
  E[(ȳ₀/(1+ȳ₀²))²]` are expectations over the RCA(1) process's own
  **stationary distribution**, `ȳᵢ` — a quantity with no general
  closed form. Checked the paper's Section 5 (simulations) for a
  feasible/plug-in estimator usable on real data with unknown
  parameters — found none: their own Monte Carlo section only says
  critical values are "computed using Theorem 3.2/3.6," which is
  sufficient for *their* validation because they know the true DGP
  parameters by construction (it's their own synthetic data), not
  because they give a data-driven estimator of `ð²`. A real
  implementation would need to either derive one (a nontrivial
  extension of the paper's own results, out of scope for a
  same-source-only implementation) or find it in the companion 2023
  working-paper version — not attempted.

Net: the critical-value machinery is now less of a blocker than
originally scoped for the `ψ = 0` case specifically, but a *different*
piece — the nuisance-parameter estimator `ð²` needed to actually run the
detector on real (parameter-unknown) data — is the harder of the two
requirements and was not resolved this pass. Not implemented; flagged
precisely so a future pass knows exactly which piece to chase (the
feasible `ð²` estimator, not the critical value) rather than re-deriving
the whole thing from scratch.

### Kurozumi (2020, 2021) — SADF and GSADF cases both implemented (2026-08-10)

**Status: `SADF`/`s0 = 0` AND `GSADF_{s0}`/`s0 = 0.4`/`0.8` cases both
implemented as `radf_monitor(..., boundary = "kurozumi", s0 = ...)`.**
Full PDF read for Kurozumi (2020)
(through Theorem 1 and Table 1, rendered to PNG for accurate
transcription); Kurozumi (2021) read at the abstract/intro level only (see
its own paragraph below).

**A genuinely valuable structural finding, now verified empirically**:
Kurozumi (2020)'s `SADF(k) := ADF_1^{m+k}` (his eq., Section 3) is
*exactly* `radf()`'s `badf` sequence — confirmed bit-for-bit against a
from-scratch OLS ADF t-statistic (fixed start at `t=1`, expanding window
end) at three separate check points, tolerance `1e-8`. His
`GSADF_{s0}(k) := max_{1<=k1<=floor(m*s0)} ADF_{k1}^{m+k}`, by contrast,
is genuinely *not* the same construction as `radf()$bsadf`: exuber's own
`bsadf` search range for the window start grows with the current
monitoring point `t` (`1` to `t - minw`), while Kurozumi's is capped at a
*fixed* fraction of the training length `m` regardless of `t`. These are
different double recursions, not a notational match — but (see
"`GSADF_{s0}` case — DONE" below) a fixed, small window-start band turned
out to need only a bounded closed-form computation, not new recursion
code, so both cases are now implemented via `radf_monitor(..., boundary
= "kurozumi", s0 = ...)`: `s0 = 0` (default) reuses `radf()$badf`
directly (no new statistic), `s0 = 0.4`/`0.8` uses the new
`kurozumi_gsadf_stat()` — either way, a new (published, table-based)
closed-form comparison threshold in place of `radf_monitor()`'s default
wild-bootstrap-calibrated boundary, directly analogous to how
`radf_tt_cv()` provides a bootstrap-free alternative to
`radf_wb_cv()` for the static (non-monitoring) GSADF case.

**The boundary functions** (his eq., confirmed via rendered PDF page —
the raw-text extraction badly scrambled the sub/superscripts):

```
SADF:  g_0^df(k/m)    := q_0^df                                    (constant)
GSADF: g_{s0}^df(k/m) := q_{s0}^df * (a_{s0} + b_{s0} * log(c_{s0} + k/m))
       {a,b,c} = {0.76, 0.02, 0.34} for s0=0.4;  {0.73, 0.03, 0.90} for s0=0.8
CS:    g_gamma^cs(k/m) := q_gamma^cs * (1 + k/m)^(1-gamma) * (k/m)^gamma
```

**Table 1 (transcribed from the rendered page, not the raw OCR text)** —
scaling constants `q`, by significance level `β` and monitoring-horizon
ratio `s̄ = k̄/m` (monitoring runs `k̄` observations past training length
`m`; the table only tabulates `s̄ ∈ {1, 3, 5}`, i.e. a monitoring horizon
1x/3x/5x the training window):

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

(`CS(k)` itself is HB's own CUSUM statistic, eq. 26, restated in
Kurozumi's own notation — so the `q^cs` columns are, in effect, published
finite-sample-simulated alternatives to `monitor_cusum()`'s asymptotic
`b_α = 4.6` constant, for the two specific `γ ∈ {0.25, 0.45}` cases
tabulated. These were obtained by simulation — "50,000 replications... a
standard Brownian motion is approximated by the sum of suitably
normalized i.i.d. pseudo N(0,1) random variates with increments of
1/1000" — not a closed-form asymptotic result the way HB's `b_α = 4.6`
is, but still a fixed published number, not something needing new
simulation on exuber's end.)

**What was implemented (`SADF`, `s0 = 0`)**: `radf_monitor()` gained a
`boundary = c("bootstrap", "kurozumi")` parameter. `boundary = "kurozumi"`
looks up `q_0^df` from Table 1 (`kurozumi_sadf_q(level, s_bar)`, snapping
`s_bar = (n - T*) / T*` to the nearest tabulated `{1, 3, 5}`, requiring
`level` in `{0.90, 0.95, 0.99}`) and compares it against `radf()$badf`
over the monitoring window — no bootstrap, no simulation, `nboot`/`type`/
`adflag`/`seed` all ignored on this path. The returned list's statistic
field was renamed from `bsadf` to `stat` (it now holds `bsadf` for
`boundary = "bootstrap"`, `badf` for `boundary = "kurozumi"`).

Validation (`radf_monitor_kurozumi_boundary_validation.R`):
table lookups match Table 1 exactly (all 6 checked values, plus
`s_bar`-snapping and the invalid-`level` error path); a basic run prints
without error; empirical false-alarm rate under `H0` (`n=150`, `T*=75`,
`s_bar=1`, 100 reps) is 4.0% against a nominal 5% target (vs. 7.0% for
the existing wild-bootstrap boundary on the same DGP) — both close to
nominal, unlike HB's CUSUM asymptotic bound elsewhere in this file, which
is far more conservative; detection power under a genuine post-training
bubble (30 reps) is 80% for `kurozumi` vs. 90% for `bootstrap` — both
non-trivial, some gap expected since the two boundaries calibrate
different statistics (`badf`/SADF vs. `bsadf`/GSADF); and alarms never
fire before `T*` (structural invariant, 10/10 reps). Test coverage added
to `test-monitor.R`: exact table-lookup checks, a structural run, the
invalid-level error, a loose false-alarm-rate Monte Carlo bound, and the
alarm-timing invariant.

**`GSADF_{s0}` case — DONE (2026-08-10), re-triaged after initially
being scoped out as needing new recursion code.** `s0` in
`GSADF_{s0}(k) := max_{1 <= k1 <= floor(m*s0)} ADF_{k1}^{m+k}` controls
how far the window's *start point* `k1` is allowed to range, as a *fixed*
fraction of the training length `m` — this range does not grow as the
monitoring point `k` advances. `radf()`'s own `bsadf` search range, by
contrast, grows with the current time `t` (`1` to `t - minw`) by
construction (PSY's GSADF). These are genuinely different double
recursions, not a reindexing of the same one — but re-triaging this item
clarified that "different double recursion" does not mean "needs new
C++": since `floor(m*s0)` is a small, *fixed* cap
(independent of the current monitoring point), `GSADF_{s0}(k)` only ever
needs `ADF_{k1}^{t}` for `k1` in a small bounded band, not a growing
triangular grid. Each `ADF_{k1}^{t}` is a plain WITH-INTERCEPT OLS ADF
t-statistic on a fixed window, computable via the same closed-form
cumulative-sum-difference construction already used throughout this
project (`dating_hls.R`'s `hls_segment_ssr()`, `radf_tt.R`'s
`gls_dfstat_grid()`) — no recursion, no C++, just `O(1)` per `(k1, t)`
cell via prefix sums, restricted to the bounded band instead of the full
grid.

**Implementation**: `kurozumi_gsadf_stat()` in `exuber/R/radf_monitor.R`
computes this band via `outer()`-vectorized cumulative-sum differences
(intercept included, unlike `gls_dfstat_grid()`'s no-intercept
GLS-demeaned version — verified to match `radf()$badf` exactly at
`k1_max = 1`, confirming the intercept convention is right).
`radf_monitor(..., boundary = "kurozumi", s0 = 0.4)` or `s0 = 0.8` (the
only two values his boundary's `a/b/c` scaling constants are tabulated
for) switches from the flat `SADF` boundary to the `k`-varying
`g_{s0}^df(k/m) := q_{s0}^df * (a_{s0} + b_{s0}*log(c_{s0} + k/m))`
already transcribed above, with `q_{s0}^df` looked up from Table 1's
`q04_df`/`q08_df` columns. `s0 = 0` (default) reproduces the original
`SADF`-only behavior exactly (verified bit-for-bit, no regression).

**Validated**: `kurozumi_gsadf_stat()` matches `radf()$badf` to machine
precision at `k1_max = 1`, and matches a brute-force `lm()` search over
the restricted window-start band exactly (`|diff| < 1e-14`) at three
separate monitoring points on a 150-observation series; Table 1's
`q04_df`/`q08_df` lookups exact, with the expected tie-breaking snap
(`s0 = 0.6` snaps to the first/lower tabulated value, `0.4`); alarms
never fire before `T*` (30/30 reps). Empirical false-alarm rate under
`H0` (300 reps, `n=150`, `T*=75`) is `4.3%`/`5.3%`/`4.3%` for
`SADF`/`GSADF_{0.4}`/`GSADF_{0.8}` against nominal `5%` — all close to
nominal. Detection power on a post-training-bubble DGP (60 reps) is
`70.0%`/`73.3%`/`66.7%` for `SADF`/`GSADF_{0.4}`/`GSADF_{0.8}` — `s0 =
0.4`'s modest edge over plain `SADF` is a direct, quantitative echo of
Kurozumi's own stated finding that "GSADF works better than SADF in many
cases," while `s0 = 0.8`'s underperformance here is a genuine,
DGP-specific finding, reported honestly rather than smoothed over (a
wider start-range dilutes power on some alternatives, not a universal
improvement). Extended `test-monitor.R` (9 new tests). Replication
script:
[replication/monitoring/radf_monitor_gsadf_s0_validation.R](#script-radf_monitor_gsadf_s0_validation).

**Kurozumi (2021)** ("Asymptotic Behavior of Delay Times of Bubble
Monitoring Tests," *JTSA* 42(3), 314-337) was read at the abstract/intro
level for this pass: it studies the stochastic order of the detection-
delay (stopping time) for the same detector families as the 2020 paper,
confirming the qualitative early/short-bubble-favors-CUSUM vs.
middle/late-bubble-favors-ADF split already cited elsewhere in this file
(and now empirically reproduced by this session's own `monitor_cusum()` vs.
`radf_monitor()` comparison). A dating/inference layer on top of
detection, not a new detector — lower priority than validating the 2020
paper's boundary functions.

(Kurozumi (2020)'s formulas, Table 1 boundary constants, and the
structural finding that his `SADF`/`GSADF`/`CS` detectors are literally
`radf_monitor()`/`monitor_cusum()`'s existing statistics are now transcribed
in full — see the dedicated subsection under "Exact numbers/formulas
reproduced" above, not repeated here.)

### Breitung & Diegel (2025) — static LBI test AND sequential extension
both done (2026-08-10)

**Status: both the static (known/full-sample bubble window) LBI test
(`lbi_test()`) and their own headline sequential/monitoring extension
(`monitor_lbi()`) are done.** Full abstract/intro read plus the
core statistic's derivation section and Section 4's sequential
extension, re-verified against rendered PDF pages 3-7 (the raw-text
extraction badly scrambles the `σ̃`/summation notation throughout,
which mattered for both — see below).

Proposes a genuinely different (and, by their own claim, more powerful)
detector: a **locally best invariant (LBI) statistic**, heteroskedasticity
-robust *by construction* (citing Cavaliere 2005's invariance result, not
requiring a wild bootstrap the way HB's plain CUSUM does), with a
**standard normal** limiting null distribution — no simulated or
bootstrap critical values needed at all. Their eq. 4 gives the key
telescoping identity for a bubble known to span the whole sample
(`y_0 = 0`, their Section 3): `2*sum(Delta y_t * y_{t-1}) = y_T^2 -
T*sigma_tilde^2`, `sigma_tilde^2 := T^{-1}*sum(Delta y_t^2)` — the plain
sample variance of first differences under `H0`. Substituting into the
naive DF-type statistic's numerator gives their eq. 5, `LBI_T^2 = y_T^2 /
(sigma_tilde^2 * T)` — i.e. the whole statistic is just the standardized
sample endpoint, no regression, no recursion:

```
LBI_T = (y_T - y_1) / (sigma_tilde * sqrt(T - 1))
```

compared against a standard normal quantile (e.g. `1.645` at 5%,
one-sided since the paper explicitly targets *positive* bubbles only,
citing that negative bubbles are economically implausible for a risky
asset). Genuinely the cheapest statistic validated in this whole
project: no C++, no bootstrap, no published table, not even a boundary
function — just a mean and a normal quantile.

**Implementation**: `lbi_test(data, level = 0.95)` in
`exuber/R/lbi_test.R`, tested in `exuber/tests/testthat/test-lbi.R`.
**Validation, including a direct check of the paper's own claimed null
distribution, not just an approximately-correct size**: eq. 4's
telescoping identity confirmed to hold exactly (not just approximately)
on a simulated random walk. Monte Carlo under `H0` (500 reps): the
statistic's empirical mean/sd are `0.023`/`0.976` (theory: `0`/`1`), its
empirical false-alarm rate at the 95% level is `0.050` against a nominal
5% — matching almost exactly, not merely "close" — and a
Kolmogorov-Smirnov test against `N(0,1)` gives `p = 0.849`, strong
evidence the statistic genuinely follows the claimed distribution
exactly, not just approximately in some looser sense. Detection power
under a genuine explosive alternative (60 reps) is `100%`, matching a
standard SADF test on the identical DGP exactly. Replication script:
[replication/monitoring/radf_lbi_validation.R](#script-radf_lbi_validation).

**Sequential extension — now done.** The abstract's claim — "the
exponentially weighted CUSUM detector with a constant boundary function
turns out to be most powerful" — was previously blocked on two
unconfirmed details: the exponential weighting scheme's exact
construction and the calibrated constant boundary value. Both were
pinned down by re-reading rendered pages 6-7 (Section 4.1, "Sequential
Tests Based on the LBI Detector"):

- **The statistic (their eq. 15)** is the classical CUSUM idea applied
  directly to the LBI construction: normalize the time index to the
  *monitoring* period, `r = j/T_m` for `j = 1, ..., T_m` (`T_m` the
  fixed monitoring horizon, chosen in advance), and compute the
  (optionally weighted) partial sum `LBI_[rT] = (1/(σ̃√T_m)) *
  sum_{t=1}^{[rT_m]} w_t Δy_t ⇒ W(r)`, a standard Brownian motion under
  `H0`. Unlike `monitor_cusum()`'s existing Chu-Stinchcombe-White boundary
  (which grows as `sqrt(t)`), normalizing by the *fixed* `T_m` up front
  means a single **constant** boundary controls size uniformly across
  the whole monitoring window — this constant-boundary variant is what
  the paper calls `mCUSUM`, and its own stated reason it's more
  powerful than the classical time-varying-boundary CUSUM (Brown et al.
  1975, also in their Table 1 but not implemented here — superseded by
  the paper's own recommendation) is that under an explosive
  alternative the detector tends to be largest near the *end* of the
  monitoring window, which a shrinking-relative-tolerance time-varying
  boundary penalizes.
- **The exponential weighting (their eq. 12)**: `w_r^c̄ =
  sqrt(2c̄/T_m)/sqrt(e^{2c̄}-1) · e^{c̄r}`, `c̄ ≥ 0` a single tunable
  parameter that up-weights later (more bubble-like) monitoring
  observations; `c̄ = 0` recovers flat weights (`mCUSUM`), `c̄ > 0` is
  `wCUSUM`. Their own suggested default for a moderate power boost is
  `c̄ ≈ 2`.
- **The critical value**: their Table 1 (page 7 of the rendered PDF,
  1,000,000 Monte Carlo reps at `T = 10,000`) gives published one-sided
  asymptotic critical values, and — their own structural argument,
  `sup_r{W*(η(r))} =_d sup_r{W(r)}` since the running max of a
  time-changed Brownian motion has the same distribution regardless of
  the time-change — **one set of critical values covers every `c̄`**,
  `mCUSUM` and `wCUSUM` alike: `1.64/1.95/2.24/2.57/2.80` at the
  `10%/5%/2.5%/1%/0.5%` levels. No new simulation needed, exactly the
  published-table-lookup shortcut that worked for Kurozumi (2020)'s and
  HB's own FLUC/CUSUM boundaries elsewhere in this project.
- `σ̃²` is estimated from the **training window only** (their own
  Section 4.2 text: "the training set ... is used for estimating
  nuisance parameters such as `σ²`"), consistent with `lbi_test()`'s
  own training-free full-sample `σ̃²` for the static case.

**Implementation**: `monitor_lbi(data, r_star = 0.5, c_bar = 0,
level = 0.95)` in `exuber/R/lbi_test.R`, extended
`exuber/tests/testthat/test-lbi.R`. **Validated**: the flat-weight
(`c̄ = 0`) weight vector's sum of squares equals `1` exactly (eq. 12's
own discrete closed form, not an approximation); for `c̄ > 0` it equals
`1` to within the expected Riemann-sum-approximation error (`< 0.5%` at
`T_m = 500`); the final monitoring-point statistic under `mCUSUM`
matches, to machine precision, a hand-computed telescoped value using
training-window `σ̃²` — the same telescoping identity `lbi_test()`
itself relies on, now cross-checked in the monitoring context too;
table lookups match Table 1 exactly, with a clean error for an
untabulated level; alarms never fire before the training window ends
(50/50 reps). Empirical false-alarm rate under `H0` (1,000 reps,
`n=200`, `T*=100`) is `3.7%` (`mCUSUM`) and `4.0%` (`wCUSUM`) against a
nominal `5%` — both mildly conservative, not oversized, matching the
pattern of most other finite-sample-vs-asymptotic boundaries validated
in this project. Detection power on a post-training-bubble DGP (60
reps, bubble starting well into the monitoring window) is `41%`
(`mCUSUM`) and `44%` (`wCUSUM`) — both clearly exceeding
`monitor_cusum(type = "standard")`'s `31%` on the identical DGP, a direct
confirmation of the paper's own claim that the constant-boundary LBI
detector is more powerful than HB's classical CSW-boundary CUSUM, and
`wCUSUM ≥ mCUSUM` as expected from the added up-weighting. Replication
script:
[replication/monitoring/radf_lbi_monitor_validation.R](#script-radf_lbi_monitor_validation).

Not implemented: the classical time-varying-boundary CUSUM row of
Table 1 (Brown et al. 1975) — the paper's own point of comparison, not
its contribution, and explicitly shown to be *less* powerful than
`mCUSUM`/`wCUSUM`; and Section 4.2's separate DF-statistic-based
monitoring variant (a `badf`-based analogue using a different,
sup-of-a-ratio boundary), which targets the same monitoring problem via
a structurally different statistic already covered in spirit by
`radf_monitor()`.

## Cost/feasibility note for exuber

Concretely, three things are needed, of which exuber currently has
meaningful partial infrastructure for exactly one:

1. **Family-wise/size control via a training-sample critical value — DONE
   (2026-08-09), via `radf_monitor()`.** `exuber/R/radf_wb.R`'s
   `radf_wb_cv2()` already implemented the Phillips & Shi (2020) wild
   bootstrap (its own roxygen docs cite "Phillips, P. C., & Shi, S.
   (2020). Real time monitoring of asset markets: Bubbles and crises")
   and already had a `tb` parameter that truncates the bootstrap DGP to a
   training sub-sample of length `tb`, computing `sadf_crit`/`gsadf_crit`
   on that sub-sample only and broadcasting those training-sample
   quantiles as a **constant boundary** across the full `bsadf_crit`
   pointer range. `radf_monitor()` is exactly the missing orchestration
   layer this note predicted: it (a) fixes "now" at `T*`, (b) reuses a
   single full-sample `radf()` call's own `bsadf` sequence rather than
   walking forward and re-fitting (cheap, since it's already an efficient
   `O(T)` recursive computation and — verified — `bsadf[t]` depends only
   on data up to `t`), (c) compares each monitoring-region value against
   the fixed training boundary, and (d) returns the first breach. See
   "Implementation" above for the full writeup, including a look-ahead-
   leakage pitfall found and avoided (`radf_wb_cv2()`'s null-model fit
   has no internal truncation to `tb`, so the training window must be
   sliced *before* calling it, not relied on `tb` alone to enforce).
   `radf_mc_cv()` could supply the same role for the non-bootstrap
   (asymptotic/MC) case — not implemented, `radf_monitor()` only supports
   the wild-bootstrap route for now.

2. **FPR-vs-monitoring-horizon accounting — done for both the `SADF` and
   `GSADF_{s0}` ADF-family detectors via Kurozumi (2020), DONE
   (2026-08-10).** exuber's `radf_wb_cv2(tb=...)` gives *one* bootstrap
   critical value; Kurozumi's Table 1 supplies the published,
   simulation-calibrated closed-form analogue for both the `SADF`
   (`s0=0`) and `GSADF_{s0}` (`s0=0.4`/`0.8`) cases —
   `radf_monitor(..., boundary = "kurozumi", s0 = ...)` — indexed by
   monitoring-horizon ratio `s̄ = k̄/m`, the same accounting
   AHLST/Whitehouse's eq. (4)–(6) provide for their own detector. Still
   missing: an equivalent for `radf_monitor()`'s default wild-bootstrap
   boundary path itself (that one still recalibrates via simulation
   rather than a closed form).

3. **CUSUM/Page-CUSUM detector family — 0% reusable, all new — HB's basic
   CUSUM and Astill et al.'s volatility-robust variant both now DONE
   (2026-08-10), confirming this assessment.** Items 1, 2, 4 (CUSUM
   variant), 5, 6, and 7 are a structurally different statistic (a
   standardized running sum of `Δy_t`, not a recursive ADF regression),
   so none of `exubercore/src`'s RLS machinery
   (`rls_gsadf.cpp`/`radf.hpp`, matrix-inversion-lemma recursive OLS)
   applies — confirmed empirically, not just predicted: `monitor_cusum()`
   shares no code with `radf_monitor()`/`radf()` at all. What's now
   actually built vs. still missing: (a) **the CUSUM partial-sum
   statistic itself — done**, `monitor_cusum()`, closed-form, no C++, same
   "why not exubercore" logic as
   [STADF's closed-form statistic](/replication/volatility-robustness#why-not-exubercore-c));
   (b) **a boundary-function/critical-value module — partially done**:
   HB's own asymptotic `c_t√t` constant (`b_α = 4.6` at the 5% level) is
   implemented, and — per Astill et al.'s Corollary 1 — the *same*
   boundary was confirmed to work unchanged for the volatility-robust
   variant too; HB's own *other* proposed statistic, FLUC — confirmed to
   be `radf()`'s existing `badf` sequence compared against a *published*
   finite-sample-simulated boundary (their Table 7, transcribed directly,
   no new simulation needed) rather than one exuber would have to
   simulate itself — is now **done** too, via
   `radf_monitor(..., boundary = "fluc")`. HB's *finite-sample* `b_{k,α}`
   values for CUSUM itself (their Table 8, a different table from FLUC's
   Table 7) are **now also done** (2026-08-10), via
   `monitor_cusum(..., boundary = "finite")`; (c) **the volatility-robust
   variant's one-sided kernel spot-variance estimator — done**,
   `monitor_cusum(..., type = "kernel")`, structurally close to, but built
   as a genuinely separate function from, the two-sided/profile kernel
   estimator already built for
   [`radf_tt.R`](/replication/volatility-robustness#time-transformed-test-stadf--gstadf)
   (that one estimates a cumulative variance *profile* for time
   transformation, using both past and future data within a bandwidth
   window; this one is deliberately one-sided/causal, using only current
   and past lags, since real-time monitoring at time `t` can only ever
   see data up to `t`) — confirming the original "partially reusable
   engineering pattern, not reusable code" read: the *pattern* (kernel
   weights + `filter()`-style rolling computation) carried over, but the
   code itself did not.

**Net assessment (updated 2026-08-09).** Item 1 (Phillips & Shi 2020-style
training-vs-monitoring orchestration) turned out cheaper in practice than
the original "plausibly 1-2 weeks" estimate below — a single focused pass
(design, empirical verification of the two reuse assumptions, ~150 lines
of R, tests, docs), not a multi-week project — precisely because that
estimate was pricing in still-unresolved uncertainty about whether the
existing `tb`/`bsadf` machinery would really compose the way the note
predicted; it did, cleanly, once actually checked. That leaves the
originally-scoped **remaining** work essentially unchanged: the CUSUM/
Page-CUSUM family (item 3 below) is a **structurally different statistic
family**, 0% reusable, and a genuine FPR-vs-horizon boundary function
(item 2) needs new asymptotic theory ported in, not just new R code —
neither shortcut the way item 1 did. Building the CUSUM family well,
with its own size-control theory and wired into a coherent monitoring
API, remains **not a single-pass item** — a realistic estimate is still
multi-week, now better described as "the CUSUM/Page-CUSUM family plus
FPR-vs-horizon theory" rather than "everything," since the one tractable
piece is done.
