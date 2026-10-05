"""Python cross-check of exuber's monitor_cusum() (Homm & Breitung 2012's
CUSUM real-time monitoring, both the "standard" statistic and Astill,
Harvey, Leybourne, Taylor & Zu (2023)'s volatility-robust "kernel"
variant). Fully deterministic (no bootstrap, no simulation at all), so
every check here is a bit-for-bit cross-check against R, not just a
Monte Carlo size/power reproduction -- unlike radf_cusum_validation.R's
own size/power sections (different RNGs make those non-reproducible
across languages; docs/monitoring.md already records the R-side numbers).

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_cusum_validation.py

Reference values: a direct R run from the exuber-project/ root,
    Rscript -e 'Sys.setenv(NOT_CRAN="true");
    devtools::load_all("exuber", quiet=TRUE); set.seed(42);
    y <- cumsum(rnorm(80));
    mc <- monitor_cusum(y, r_star=0.5, b_alpha=4.6, boundary="asymptotic");
    dput(mc$T_star); dput(round(tail(as.numeric(mc$S), 5), 10));
    dput(round(tail(as.numeric(mc$boundary), 5), 10))'
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor_cusum import monitor_cusum  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Formula-exact check vs brute-force recomputation ===")
    t_star = 40
    dy = np.diff(Y42)
    expected_s, expected_b = [], []
    for t in range(t_star, len(Y42)):
        sigma2 = np.sum(dy[:t] ** 2) / t
        s_t = (Y42[t] - Y42[t_star - 1]) / np.sqrt(sigma2)
        c_t = np.sqrt(4.6 + np.log((t + 1) / t_star))
        expected_s.append(s_t)
        expected_b.append(c_t * np.sqrt(t + 1))
    res = monitor_cusum(Y42, r_star=0.5, b_alpha=4.6, boundary="asymptotic")
    np.testing.assert_allclose(res.stat[:, 0], expected_s, atol=1e-10)
    np.testing.assert_allclose(res.boundary[:, 0], expected_b, atol=1e-10)
    print(f"max |diff| S: {np.max(np.abs(res.stat[:, 0] - expected_s)):.2e}, "
          f"boundary: {np.max(np.abs(res.boundary[:, 0] - expected_b)):.2e}")

    print("\n=== 2. Basic run -- cross-check vs R ===")
    assert res.t_star == 40
    expected_s_tail = np.array(
        [3.6727107022, 4.4016203054, 4.8602235379, 4.0370328378, 3.0016493394]
    )
    expected_b_tail = np.array(
        [19.9594813397, 20.1153995614, 20.2704388473, 20.4246151364, 20.5779438828]
    )
    np.testing.assert_allclose(res.stat[-5:, 0], expected_s_tail, atol=1e-8)
    np.testing.assert_allclose(res.boundary[-5:, 0], expected_b_tail, atol=1e-8)
    print(f"T_star={res.t_star}, S[-5:]={res.stat[-5:, 0]}")
    assert np.isnan(res.alarm[0])

    print("\n=== 3. Structural check: alarm never before T_star ===")
    rng = np.random.default_rng(1)
    for _ in range(10):
        y = np.cumsum(rng.normal(size=150))
        out = monitor_cusum(y, r_star=0.5)
        if not np.isnan(out.alarm[0]):
            assert out.alarm[0] >= out.t_star
    print("all alarms (if any) fired at/after T_star: OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
