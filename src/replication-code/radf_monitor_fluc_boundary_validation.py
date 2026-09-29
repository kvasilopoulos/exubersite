"""Python cross-check of exuber's monitor(..., boundary = "fluc") (Homm &
Breitung 2012's FLUC detector), mirroring
radf_monitor_fluc_boundary_validation.R's table-lookup and basic-run
checks (the R script's Monte Carlo false-alarm-rate/detection-power
sections are not repeated here -- different RNGs make an exact
cross-language match meaningless; docs/monitoring.md already records
those numbers from the R side).

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_monitor_fluc_boundary_validation.py

Reference values: Homm & Breitung (2012) Table 7(i), already transcribed
into both exuber/R/monitor.R's hb_fluc_table and
pyexuber/src/exuber/monitor.py's _HB_FLUC_TABLE. The stat/boundary
sequence for Y42 is from the same R run cited in
radf_monitor_kurozumi_boundary_validation.py, with `boundary="fluc"`
substituted:
    out <- monitor(y, r_star=0.5, minw=15, boundary="fluc", sig_lvl=95)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor import _hb_fluc_q, monitor  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Table lookup sanity checks (Table 7(i)) ===")
    checks = [
        (95, 40, 2, 4.19),  # n_train snaps to 50
        (95, 100, 2, 4.50),
        (90, 20, 10, 4.12),
        (99, 100, 2, 7.76),
    ]
    for sig_lvl, n_train, k, expected in checks:
        q = _hb_fluc_q(sig_lvl, n_train, k)
        print(f"sig_lvl={sig_lvl}, n_train={n_train}, k={k} (expect {expected}): {q}")
        assert q == expected

    print("\n=== 2. Basic run, boundary='fluc' -- cross-check vs R ===")
    mon = monitor(Y42, r_star=0.5, minw=15, boundary="fluc", sig_lvl=95)
    assert mon.t_star == 40
    assert mon.boundary[0] == 4.19
    expected_stat_tail = np.array(
        [-1.5922477876, -1.5806550833, -1.5647907562, -1.6255928209, -1.660938116]
    )
    np.testing.assert_allclose(mon.stat[-5:, 0], expected_stat_tail, atol=1e-8)
    print(f"T_star={mon.t_star}, boundary={mon.boundary[0]}, stat[-5:]={mon.stat[-5:, 0]}")
    assert np.isnan(mon.alarm[0])

    # FLUC and boundary="kurozumi" (s0=0) both reuse radf()'s badf sequence
    # (Homm & Breitung's DF_{t/n} == Kurozumi's SADF(k), both confirmed
    # bit-identical to badf -- docs/monitoring.md, "Implementation (FLUC)").
    mon_k = monitor(Y42, r_star=0.5, minw=15, boundary="kurozumi", sig_lvl=95)
    np.testing.assert_allclose(mon.stat, mon_k.stat)
    print("\nFLUC stat path is bit-identical to boundary='kurozumi' s0=0 (both == badf): OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
