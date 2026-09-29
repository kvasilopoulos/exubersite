"""Validation of contagion_reg() -- Greenaway-McGrevy & Phillips (2016)'s
bubble contagion regression, minimum-viable subset (fixed-window AR(1)
sequence, single-delay Nadaraya-Watson regression, LOOCV bandwidth). Same
folder/base name as the R script it cross-checks: radf_contagion_validation.R.
See docs/multivariate.md, "Contagion regression", for the full write-up
(including two real bugs the R implementation found and fixed: a
window-width off-by-one, and a matrix-orientation bug in the LOOCV SSE
helper -- both are checked directly here, not just inherited by copying
the fixed formulas).

No published numeric table exists to validate against -- the source
paper's own results are Figures 7-8, not tabulated numbers (same
situation the R script's own header documents). Validated instead via
brute-force cross-checks of each closed-form piece (to near machine
precision, no RNG involved -- these are NOT stochastic checks, unlike
radf_common_validation.py/radf_cobubble_validation.py), plus a
directional sensible-behavior check. Run standalone from the
exuber-project/ root:

  uv run --project pyexuber python docs/replication/multivariate/radf_contagion_validation.py

Not imported by pyexuber's own pytest suite for the same repo-boundary
reason as radf_common_validation.py (see that file's docstring); these
checks are re-implemented directly in pyexuber/tests/test_multivariate.py.
"""

import numpy as np

from exuber.contagion_reg import (
    _contagion_bandwidth_cv,
    _contagion_fixed_window_beta,
    _contagion_loocv_sse,
    _contagion_nw_delta2,
    contagion_reg,
)


def check_fixed_window_beta(seed: int = 1) -> None:
    """eq. 1 vs. brute-force lstsq at several window-end dates."""
    rng = np.random.default_rng(seed)
    n, S = 150, 50
    core = np.cumsum(rng.normal(size=n))
    t_core, beta_core = _contagion_fixed_window_beta(core, S)

    pos = {int(t): i for i, t in enumerate(t_core)}
    for t_check in (60, 80, 100, 130, 150):
        win = core[t_check - S : t_check]  # S levels, python 0-indexed
        design = np.column_stack([np.ones(S - 1), win[:-1]])
        beta_lstsq, *_ = np.linalg.lstsq(design, win[1:], rcond=None)
        cf = beta_core[pos[t_check]]
        diff = abs(cf - beta_lstsq[1])
        print(f"t={t_check} closed-form={cf:.8f} lstsq={beta_lstsq[1]:.8f} |diff|={diff:.2e}")
        assert diff < 1e-8


def check_nw_ratio(seed: int = 1) -> None:
    """eq. 6 vs. a manual weighted-least-squares ratio."""
    rng = np.random.default_rng(seed)
    n, S = 150, 50
    core = np.cumsum(rng.normal(size=n))
    y = 0.5 * core + np.cumsum(rng.normal(scale=0.5, size=n))
    t_core, beta_core = _contagion_fixed_window_beta(core, S)
    t_j, beta_j = _contagion_fixed_window_beta(y, S)
    d, r_test, h_test = 2, np.array([0.5]), 0.2

    fast = _contagion_nw_delta2(t_core, beta_core, t_j, beta_j, n, r_test, h_test, d)[0]

    bcore_c = beta_core - beta_core.mean()
    bj_c = beta_j - beta_j.mean()
    pos = {int(t): i for i, t in enumerate(t_core)}
    idx = np.array([pos.get(int(t) - d, -1) for t in t_j])
    valid = idx >= 0
    s2 = t_j[valid]
    bjc = bj_c[valid]
    csh = bcore_c[idx[valid]]
    w = np.exp(-0.5 * ((s2 / n - r_test[0]) / h_test) ** 2) / np.sqrt(2 * np.pi) / h_test
    manual = np.sum(w * bjc * csh) / np.sum(w * csh**2)

    diff = abs(fast - manual)
    print(f"fast={fast:.10f} manual={manual:.10f} |diff|={diff:.2e}")
    assert diff < 1e-10
    return t_core, beta_core, t_j, beta_j, bjc, csh, s2, n, d, h_test


def check_loocv_sse(t_core, beta_core, t_j, beta_j, bjc, csh, s2, n, d, h_test) -> None:
    """eq. 7 vs. a manual leave-one-out double loop (catches the exact
    matrix-orientation bug the R implementation found: K.T @ v, not K @ v,
    since the kernel weight matrix isn't symmetric)."""
    fast_sse = _contagion_loocv_sse(h_test, t_core, beta_core, t_j, beta_j, n, d)

    m = len(s2)
    manual_sse = 0.0
    for i in range(m):
        r_i = s2[i] / m
        num = den = 0.0
        for p in range(m):
            if p == i:
                continue
            wp = np.exp(-0.5 * ((s2[p] / n - r_i) / h_test) ** 2) / np.sqrt(2 * np.pi) / h_test
            num += wp * bjc[p] * csh[p]
            den += wp * csh[p] ** 2
        pred = (num / den) * csh[i]
        manual_sse += (bjc[i] - pred) ** 2

    diff = abs(fast_sse - manual_sse)
    print(f"fast={fast_sse:.10f} manual={manual_sse:.10f} |diff|={diff:.2e}")
    assert diff < 1e-8


def check_bandwidth_cv(t_core, beta_core, t_j, beta_j, n, d) -> None:
    """Interior optimum, SSE no worse than either H_T endpoint (eq. 7)."""
    h_opt = _contagion_bandwidth_cv(t_core, beta_core, t_j, beta_j, n, d)
    m = len(beta_j)
    H_T = (m ** (-1 / 2), m ** (-1 / 10))
    sse_opt = _contagion_loocv_sse(h_opt, t_core, beta_core, t_j, beta_j, n, d)
    sse_lo = _contagion_loocv_sse(H_T[0], t_core, beta_core, t_j, beta_j, n, d)
    sse_hi = _contagion_loocv_sse(H_T[1], t_core, beta_core, t_j, beta_j, n, d)
    print(f"H_T: {H_T}  h_opt: {h_opt}")
    print(f"SSE at h_opt: {sse_opt}  SSE at H_T[0]: {sse_lo}  SSE at H_T[1]: {sse_hi}")
    assert sse_opt <= sse_lo + 1e-8
    assert sse_opt <= sse_hi + 1e-8


def check_sensible_behavior(seed_base: int = 1000, nrep: int = 15) -> None:
    """Directional check (no ground truth to match): a satellite series
    whose local persistence genuinely tracks the core's own should show a
    visibly wider range of estimated delta_2(r) than an independent series."""
    n, S = 150, 50
    planted_range = np.empty(nrep)
    indep_range = np.empty(nrep)
    for i in range(nrep):
        rng = np.random.default_rng(seed_base + i)
        core_i = np.cumsum(rng.normal(size=n))
        y_planted = np.empty(n)
        y_planted[:10] = rng.normal(size=10)
        for t in range(10, n):
            local_rho = 0.5 + 0.4 * np.tanh(
                (core_i[max(t - 3, 0)] - core_i[max(t - 13, 0)]) / 5
            )
            y_planted[t] = local_rho * y_planted[t - 1] + rng.normal()
        y_indep = np.cumsum(rng.normal(size=n))

        out_planted = contagion_reg(y_planted, core_i, S=S, d=3, h=0.3)
        out_indep = contagion_reg(y_indep, core_i, S=S, d=3, h=0.3)
        planted_range[i] = np.ptp(out_planted.delta2)
        indep_range[i] = np.ptp(out_indep.delta2)

    print(f"mean range(delta2), planted contagion: {planted_range.mean():.4f}")
    print(f"mean range(delta2), independent series: {indep_range.mean():.4f}")


if __name__ == "__main__":
    print("=== 1. Fixed-window AR(1) coefficient sequence (eq. 1) ===")
    check_fixed_window_beta()
    print("\n=== 2. Nadaraya-Watson ratio (eq. 6) ===")
    ctx = check_nw_ratio()
    print("\n=== 3. LOOCV SSE (eq. 7) ===")
    check_loocv_sse(*ctx)
    print("\n=== 4. Bandwidth CV ===")
    check_bandwidth_cv(ctx[0], ctx[1], ctx[2], ctx[3], ctx[7], ctx[8])
    print("\n=== 5. Directional sensible-behavior check ===")
    check_sensible_behavior()
    print("\ndone")
