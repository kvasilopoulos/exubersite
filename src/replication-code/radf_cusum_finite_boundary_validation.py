"""Python cross-check of exuber's monitor_cusum(..., boundary = "finite")
(Homm & Breitung 2012's finite-sample CUSUM boundary, their Table 8),
mirroring radf_cusum_finite_boundary_validation.R's table-lookup and
basic-run checks.

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_cusum_finite_boundary_validation.py

Reference values: Homm & Breitung (2012) Table 8(i), already transcribed
into both exuber/R/monitor_cusum.R's hb_cusum_finite_table and
pyexuber/src/exuber/monitor.py's _HB_CUSUM_FINITE_TABLE. b_alpha/stat
tail from the same R run cited in radf_cusum_validation.py, with
`boundary="finite"` substituted.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor_cusum import _hb_cusum_finite_q, monitor_cusum  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Table lookup sanity checks (Table 8(i)) ===")
    checks = [(95, 40, 2, 1.43), (95, 100, 2, 1.51), (90, 20, 10, 2.02), (99, 100, 2, 2.86)]
    for sig_lvl, n_train, k, expected in checks:
        q = _hb_cusum_finite_q(sig_lvl, n_train, k)
        print(f"sig_lvl={sig_lvl}, n_train={n_train}, k={k} (expect {expected}): {q}")
        assert q == expected

    print("\n=== 2. Basic run, boundary='finite' -- cross-check vs R ===")
    res = monitor_cusum(Y42, r_star=0.5, boundary="finite", sig_lvl=95)
    assert res.t_star == 40
    assert res.b_alpha == 1.43
    expected_s_tail = np.array(
        [3.6727107022, 4.4016203054, 4.8602235379, 4.0370328378, 3.0016493394]
    )
    np.testing.assert_allclose(res.stat[-5:, 0], expected_s_tail, atol=1e-8)
    print(f"T_star={res.t_star}, b_alpha={res.b_alpha}, S[-5:]={res.stat[-5:, 0]}")

    # The finite-sample boundary is a stricter (smaller) constant than the
    # conservative asymptotic default at this (n, k) -- same qualitative
    # finding as docs/monitoring.md's own R-side validation.
    res_asy = monitor_cusum(Y42, r_star=0.5, b_alpha=4.6, boundary="asymptotic")
    assert res.b_alpha < res_asy.b_alpha
    print(f"finite b_alpha ({res.b_alpha}) < asymptotic b_alpha ({res_asy.b_alpha}): OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
