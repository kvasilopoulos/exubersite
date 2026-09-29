"""Validation of cobubble_test() -- Evripidou, Harvey, Leybourne & Sollis
(2022)'s co-explosive behaviour test. Same folder/base name as the R
script it cross-checks: radf_cobubble_validation.R. See
docs/multivariate.md, "Co-bubble test", for the full write-up.

No RNG-bit-exact comparison against R is attempted (numpy's Generator vs
R's RNG -- see sim.py's module docstring for why that's a project-wide
non-goal); this reproduces the SAME five checks the R script runs, each
against its own independently-generated data, and (like the R script)
mostly cares about direction/magnitude (empirical size near nominal,
power high, lag recovered) rather than exact numbers. Run standalone from
the exuber-project/ root:

  uv run --project pyexuber python docs/replication/multivariate/radf_cobubble_validation.py

Not imported by pyexuber's own pytest suite for the same repo-boundary
reason as radf_common_validation.py (see that file's docstring); these
checks are re-implemented directly in pyexuber/tests/test_multivariate.py.
"""

import numpy as np

from exuber.cobubble_test import _coexplosive_select_lag, _coexplosive_stat, cobubble_test


def check_formula_exact(seed: int = 1) -> None:
    """coexplosive_stat() vs. an independently written brute-force
    computation (separate lstsq call, manual loop-based cumulative sum
    instead of vectorized cumsum) -- mirrors the R script's check 1."""
    rng = np.random.default_rng(seed)
    tn = 100
    x = rng.normal(size=tn)
    y = 2 + 0.5 * x + rng.normal(size=tn)
    lag = 2

    S, _resid, _sigma2, _n = _coexplosive_stat(y, x, lag)

    lo, hi = max(lag, 0), tn + min(lag, 0)
    yy, xx = y[lo:hi], x[lo - lag : hi - lag]
    design = np.column_stack([np.ones(len(yy)), xx])
    beta, *_ = np.linalg.lstsq(design, yy, rcond=None)
    e_brute = yy - design @ beta
    n_brute = len(e_brute)
    sigma2_brute = np.sum(e_brute**2) / n_brute
    running = 0.0
    s_manual_sum = 0.0
    for e in e_brute:
        running += e
        s_manual_sum += running**2
    s_brute = s_manual_sum / (sigma2_brute * n_brute**2)

    diff = abs(S - s_brute)
    print(f"package S={S:.8f}  brute S={s_brute:.8f}  diff={diff:.2e}")
    assert diff < 1e-8
    print("OK: exact match\n")


def _build_dgp(seed: int, hetero: bool) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    tn, te = 150, 90
    ex = np.cumsum(rng.normal(size=te))
    expl = ex[-1] * 1.05 ** np.arange(1, tn - te + 1) + np.cumsum(
        rng.normal(scale=0.3, size=tn - te)
    )
    x = np.concatenate([ex, expl])
    if hetero:
        sd_pattern = np.concatenate([np.ones(tn // 2), np.full(tn - tn // 2, 4.0)])
        y = 1 + 0.8 * x + rng.normal(size=tn) * sd_pattern
    else:
        y = 1 + 0.8 * x + rng.normal(size=tn)
    return y, x


def check_size(seed_base: int, hetero: bool, reps: int = 100) -> float:
    rejections = []
    for i in range(reps):
        y, x = _build_dgp(seed_base + i, hetero)
        out = cobubble_test(y, x, lag=0, nboot=199, seed=1)
        rejections.append(out.reject)
    size = float(np.mean(rejections))
    label = "heteroskedastic" if hetero else "homoskedastic"
    print(f"Empirical size ({label}, {reps} reps): {size:.3f} (nominal 0.05)")
    return size


def check_power(seed_base: int = 3000, reps: int = 60) -> float:
    rejections = []
    for i in range(reps):
        rng = np.random.default_rng(seed_base + i)
        tn, te = 150, 90
        ex = np.cumsum(rng.normal(size=te))
        expl_x = ex[-1] * 1.05 ** np.arange(1, tn - te + 1) + np.cumsum(
            rng.normal(scale=0.3, size=tn - te)
        )
        x = np.concatenate([ex, expl_x])
        ey = np.cumsum(rng.normal(size=te))
        expl_y = ey[-1] * 1.05 ** np.arange(1, tn - te + 1) + np.cumsum(
            rng.normal(scale=0.3, size=tn - te)
        )
        y = np.concatenate([ey, expl_y])
        out = cobubble_test(y, x, lag=0, nboot=199, seed=1)
        rejections.append(out.reject)
    power = float(np.mean(rejections))
    print(f"Empirical power (H1, independent explosive episodes, {reps} reps): {power:.3f}")
    return power


def check_lag_recovery(true_lag: int = 3, seeds: int = 20) -> list[int]:
    estimated = []
    for seed in range(1, seeds + 1):
        rng = np.random.default_rng(seed)
        tn, te = 200, 120
        ex = np.cumsum(rng.normal(size=te))
        expl = ex[-1] * 1.06 ** np.arange(1, tn - te + 1) + np.cumsum(
            rng.normal(scale=0.3, size=tn - te)
        )
        x = np.concatenate([ex, expl])
        y = np.full(tn, np.nan)
        for t in range(true_lag, tn):
            y[t] = 1 + 0.8 * x[t - true_lag] + rng.normal(scale=0.5)
        y[:true_lag] = x[:true_lag] + rng.normal(scale=0.5, size=true_lag)
        estimated.append(_coexplosive_select_lag(y, x, range(-6, 7)))
    print(f"Estimated lags across {seeds} seeds (true={true_lag}): {estimated}")
    return estimated


if __name__ == "__main__":
    print("=== 1. Formula check ===")
    check_formula_exact()
    print("=== 2. Size under H0, homoskedastic errors ===")
    check_size(seed_base=1, hetero=False)
    print("\n=== 3. Size under H0, heteroskedastic errors ===")
    check_size(seed_base=1, hetero=True)
    print("\n=== 4. Power under H1 ===")
    check_power()
    print("\n=== 5. Lag recovery (true lag = 3) ===")
    check_lag_recovery()
