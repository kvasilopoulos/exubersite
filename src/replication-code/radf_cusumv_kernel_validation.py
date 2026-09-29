"""Python cross-check of exuber's monitor_cusum(..., type = "kernel")
(Astill, Harvey, Leybourne, Taylor & Zu (2023)'s volatility-robust
"CUSUMV" variant), mirroring radf_cusumv_kernel_validation.R's
formula-exact and basic-run checks.

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_cusumv_kernel_validation.py

Reference values: stat tail from the same R run cited in
radf_cusum_validation.py, with `type="kernel", h=20,
kernel="gaussian"` substituted.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor_cusum import _one_sided_kernel_spot_vol, monitor_cusum  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Kernel spot-variance: causal + starts-at-1 checks ===")
    rng = np.random.default_rng(0)
    dy = rng.normal(size=50)
    sigma2 = _one_sided_kernel_spot_vol(dy, h=20)
    np.testing.assert_allclose(sigma2[:20], 1.0)
    print("sigma2_j == 1 for j <= h: OK (AHLTZ's own convention)")

    dy2 = dy.copy()
    dy2[30] += 100.0
    sigma2_2 = _one_sided_kernel_spot_vol(dy2, h=20)
    np.testing.assert_allclose(sigma2[:30], sigma2_2[:30])
    print("perturbing a future observation leaves earlier spot-variance estimates unchanged: OK "
          "(one-sided/causal, matching real-time monitoring's information set)")

    print("\n=== 2. Basic run, type='kernel' -- cross-check vs R ===")
    res = monitor_cusum(Y42, r_star=0.5, type="kernel", h=20, kernel="gaussian")
    expected_tail = np.array(
        [4.2409795835, 5.0516103362, 5.5494862292, 4.607070108, 3.2258025215]
    )
    np.testing.assert_allclose(res.stat[-5:, 0], expected_tail, atol=1e-8)
    print(f"S[-5:]={res.stat[-5:, 0]}")
    assert np.isnan(res.alarm[0])

    print("\n=== 3. Corollary 1: standard and kernel share the same boundary formula ===")
    res_std = monitor_cusum(Y42, r_star=0.5)
    np.testing.assert_allclose(res.boundary, res_std.boundary)
    print("boundary paths identical for type='standard' and type='kernel': OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
