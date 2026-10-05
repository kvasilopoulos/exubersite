"""Python cross-check of exuber's monitor_lbi() (Breitung & Diegel 2025's
sequential mCUSUM/wCUSUM extension), mirroring
radf_lbi_monitor_validation.R's weight-normalization, table-lookup and
basic-run checks.

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_lbi_monitor_validation.py

Reference values: a direct R run from the exuber-project/ root,
    Rscript -e '...; ml0 <- monitor_lbi(y, r_star=0.5, c_bar=0, sig_lvl=95);
    ml2 <- monitor_lbi(y, r_star=0.5, c_bar=2, sig_lvl=95);
    dput(ml0$T_star); dput(ml0$boundary);
    dput(round(tail(as.numeric(ml0$stat), 5), 10));
    dput(round(tail(as.numeric(ml2$stat), 5), 10))'
same Y42 series as radf_monitor_kurozumi_boundary_validation.py.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.lbi_test import _bd_cusum_q, _bd_cusum_weights, monitor_lbi  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Table 1 lookup checks ===")
    checks = [(90, 1.64), (95, 1.95), (97.5, 2.24), (99, 2.57), (99.5, 2.80)]
    for sig_lvl, expected in checks:
        q = _bd_cusum_q(sig_lvl)
        print(f"sig_lvl={sig_lvl} (expect {expected}): {q}")
        assert q == expected

    print("\n=== 2. Weight normalization (eq. 12) ===")
    w0 = _bd_cusum_weights(500, 0)
    print(f"c_bar=0: sum(w^2)={np.sum(w0**2):.10f} (expect exactly 1)")
    assert abs(np.sum(w0**2) - 1.0) < 1e-12
    w2 = _bd_cusum_weights(500, 2)
    print(f"c_bar=2: sum(w^2)={np.sum(w2**2):.6f} (expect ~1, Riemann-sum approx error)")
    assert abs(np.sum(w2**2) - 1.0) < 5e-3

    print("\n=== 3. Basic run, c_bar=0 (mCUSUM) -- cross-check vs R ===")
    ml0 = monitor_lbi(Y42, r_star=0.5, c_bar=0, sig_lvl=95)
    assert ml0.t_star == 40
    assert ml0.boundary == 1.95
    expected_tail0 = np.array(
        [0.5187423335, 0.6196912485, 0.6806364851, 0.5642336848, 0.4197078303]
    )
    np.testing.assert_allclose(ml0.stat[-5:, 0], expected_tail0, atol=1e-8)
    print(f"T_star={ml0.t_star}, boundary={ml0.boundary}, stat[-5:]={ml0.stat[-5:, 0]}")

    print("\n=== 4. Basic run, c_bar=2 (wCUSUM) -- cross-check vs R ===")
    ml2 = monitor_lbi(Y42, r_star=0.5, c_bar=2, sig_lvl=95)
    expected_tail2 = np.array(
        [0.4699810166, 0.6453696898, 0.7566848661, 0.5331770251, 0.2414413064]
    )
    np.testing.assert_allclose(ml2.stat[-5:, 0], expected_tail2, atol=1e-8)
    print(f"stat[-5:]={ml2.stat[-5:, 0]}")

    print("\n=== 5. Structural check: alarm never before T_star ===")
    rng = np.random.default_rng(4)
    for _ in range(10):
        y = np.cumsum(rng.normal(size=150))
        out = monitor_lbi(y, r_star=0.5)
        if not np.isnan(out.alarm[0]):
            assert out.alarm[0] >= out.t_star
    print("all alarms (if any) fired at/after T_star: OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
