"""Python cross-check of exuber's monitor(..., boundary = "kurozumi",
s0 = 0.4/0.8) (Kurozumi 2020's GSADF_{s0} generalization), mirroring
radf_monitor_gsadf_s0_validation.R's table-lookup, formula-exact and
basic-run checks. The closed-form GSADF_{s0}(k) statistic is entirely
self-contained (no radf()/C++ extension call, unlike s0=0/"fluc"), so
every check below runs fully offline.

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_monitor_gsadf_s0_validation.py

Reference values:
  - Table 1's q04_df/q08_df columns, transcribed into both
    exuber/R/monitor.R and pyexuber/src/exuber/monitor.py.
  - The formula-exact brute-force check reimplements the with-intercept
    OLS ADF t-statistic independently (np.linalg.lstsq per window) and
    compares against the vectorized cumulative-sum construction -- same
    check exuber/R/monitor.R's own validation performed against lm().
  - The Y42/boundary/stat tail values: a direct R run
    (`monitor(y, r_star=0.5, minw=15, boundary="kurozumi", s0=0.4,
    sig_lvl=95)`), same series/command family as
    radf_monitor_kurozumi_boundary_validation.py, run from the
    exuber-project/ root on 2026-09-16.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor import _kurozumi_gsadf_q, _kurozumi_gsadf_stat, monitor  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def _brute_force_gsadf_stat(y: np.ndarray, t_star: int, s0: float) -> np.ndarray:
    n = len(y)
    dy = np.diff(y)
    ylag = y[: n - 1]
    k1_max = max(int(np.floor(t_star * s0)), 1)

    out = []
    for t in range(t_star - 1, n - 1):
        best = -np.inf
        for k1 in range(k1_max):
            yy, dd = ylag[k1 : t + 1], dy[k1 : t + 1]
            x_mat = np.column_stack([np.ones(len(yy)), yy])
            beta, _, _, _ = np.linalg.lstsq(x_mat, dd, rcond=None)
            resid = dd - x_mat @ beta
            sigma2 = np.sum(resid**2) / (len(yy) - 2)
            se = np.sqrt(sigma2 * np.linalg.inv(x_mat.T @ x_mat)[1, 1])
            best = max(best, beta[1] / se)
        out.append(best)
    return np.array(out)


def main() -> None:
    print("=== 1. Table lookup sanity checks (q04_df/q08_df) ===")
    checks = [(95, 1, 0.4, 1.8081), (95, 1, 0.8, 2.3330), (95, 1, 0.6, 1.8081)]
    for sig_lvl, s_bar, s0, expected in checks:
        q = _kurozumi_gsadf_q(sig_lvl, s_bar, s0)
        print(f"sig_lvl={sig_lvl}, s_bar={s_bar}, s0={s0} (expect {expected}): {q}")
        assert q == expected

    print("\n=== 2. Formula-exact check vs brute-force OLS ===")
    for s0 in (0.4, 0.8):
        expected = _brute_force_gsadf_stat(Y42, 40, s0)
        actual = _kurozumi_gsadf_stat(Y42, 40, s0)
        np.testing.assert_allclose(actual, expected, atol=1e-8)
        print(f"s0={s0}: vectorized construction matches brute-force lstsq scan (max |diff| = "
              f"{np.max(np.abs(actual - expected)):.2e})")

    print("\n=== 3. Basic run, boundary='kurozumi', s0=0.4 -- cross-check vs R ===")
    mon = monitor(Y42, r_star=0.5, minw=15, boundary="kurozumi", s0=0.4, sig_lvl=95)
    assert mon.t_star == 40
    expected_boundary_tail = np.array(
        [1.3819348577, 1.3826566781, 1.3833643719, 1.3840584815, 1.3847395186]
    )
    expected_stat_tail = np.array(
        [-1.5282745861, -1.5179824638, -1.5031546411, -1.5618190616, -1.5957597584]
    )
    np.testing.assert_allclose(mon.boundary[-5:], expected_boundary_tail, atol=1e-8)
    np.testing.assert_allclose(mon.stat[-5:, 0], expected_stat_tail, atol=1e-8)
    assert np.isnan(mon.alarm[0])
    print(f"T_star={mon.t_star}, boundary[-5:]={mon.boundary[-5:]}")

    print("\n=== 4. Structural check: alarm never before T_star (s0=0.4 and 0.8) ===")
    rng = np.random.default_rng(0)
    for s0 in (0.4, 0.8):
        for _ in range(15):
            y = np.cumsum(rng.normal(size=150))
            out = monitor(y, r_star=0.5, minw=20, boundary="kurozumi", s0=s0)
            if not np.isnan(out.alarm[0]):
                assert out.alarm[0] >= out.t_star
    print("all alarms (if any) fired at/after T_star, for s0 in {0.4, 0.8}: OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
