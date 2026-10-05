"""Python cross-check of exuber's monitor(..., boundary = "kurozumi", s0 = 0)
(Kurozumi 2020's closed-form SADF boundary), mirroring
radf_monitor_kurozumi_boundary_validation.R's checks 1 and 5 (checks 2-4 of
the R script exercise `boundary = "bootstrap"`, which pyexuber does not
implement -- see docs/parity.md and exuber.monitor's module docstring).

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_monitor_kurozumi_boundary_validation.py

Reference values below (Table 1 lookups, and the deterministic
badf/boundary sequence for a fixed series) come from two sources:
  - Kurozumi (2020) Table 1 itself, already transcribed into both
    exuber/R/monitor.R and pyexuber/src/exuber/monitor.py -- checked here
    for transcription-consistency between the two ports, not re-derived.
  - A direct R run, for a bit-for-bit cross-check of the statistic path:
    `Rscript -e 'Sys.setenv(NOT_CRAN="true"); devtools::load_all("exuber",
    quiet=TRUE); set.seed(42); y <- cumsum(rnorm(80));
    out <- monitor(y, r_star=0.5, minw=15, boundary="kurozumi", sig_lvl=95);
    dput(out$T_star); dput(unname(out$boundary));
    dput(round(tail(as.numeric(out$stat), 5), 10))'`
    run from the exuber-project/ root.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor import _kurozumi_sadf_q, monitor  # noqa: E402

# set.seed(42); y <- cumsum(rnorm(80))
Y42 = np.array(
    [
        1.3709584471, 0.8062602758, 1.1693886871, 1.802251292, 2.2065196152,
        2.1003950991, 3.6119170965, 3.5172580581, 5.535681772, 5.4729676729,
        6.7778373272, 9.0644827199, 7.6756220188, 7.3968332519, 7.2635119156,
        7.8994623136, 7.6152093922, 4.9587539713, 2.5182870427, 3.8384003885,
        3.5317617944, 1.7504533604, 1.5785360046, 2.7932107038, 4.6884041651,
        4.2579350335, 4.0006656507, 2.2375025655, 2.6975999203, 2.0576050444,
        2.5130551676, 3.2178925048, 4.2529960268, 3.6440696514, 4.1490247747,
        2.4320160956, 1.6475570873, 0.7966494931, -1.6175581569, -1.58143555,
        -1.3754369498, -1.7364942483, -0.9783310126, -1.7050358397, -3.0733168841,
        -2.6404988582, -3.4518920344, -2.0077907727, -2.4392369753, -1.7835890919,
        -1.4616638267, -2.2455027676, -0.6697752478, -0.0268759421, 0.0628847045,
        0.3394354518, 1.0187242679, 1.1085571544, -1.8845329287, -1.5996499752,
        -1.9668846179, -1.7816540531, -1.1998303257, 0.1999065016, -0.5273855579,
        0.7751570742, 1.1110051939, 2.1495112926, 3.0702398609, 3.7911180238,
        2.7479990852, 2.6578126986, 3.2813308606, 2.3278075028, 1.7849786883,
        2.3659751859, 3.1341539238, 3.5979215123, 2.7121452149, 1.6123643163,
    ]
)


def main() -> None:
    print("=== 1. Table lookup sanity checks ===")
    checks = [
        (95, 1, 1.0381), (95, 1.2, 1.0381), (95, 3, 1.3330),
        (95, 5, 1.4255), (90, 1, 0.6946), (99, 1, 1.6474),
    ]
    for sig_lvl, s_bar, expected in checks:
        q = _kurozumi_sadf_q(sig_lvl, s_bar)
        print(f"sig_lvl={sig_lvl}, s_bar={s_bar} (expect {expected}): {q}")
        assert q == expected

    print("\n=== 2. Basic run, boundary='kurozumi', s0=0 -- cross-check vs R ===")
    mon = monitor(Y42, r_star=0.5, minw=15, boundary="kurozumi", sig_lvl=95)
    assert mon.t_star == 40
    assert mon.boundary[0] == 1.0381
    expected_stat_tail = np.array(
        [-1.5922477876, -1.5806550833, -1.5647907562, -1.6255928209, -1.660938116]
    )
    np.testing.assert_allclose(mon.stat[-5:, 0], expected_stat_tail, atol=1e-8)
    print(f"T_star={mon.t_star}, boundary={mon.boundary[0]}, stat[-5:]={mon.stat[-5:, 0]}")
    assert np.isnan(mon.alarm[0])

    print("\n=== 3. Structural check: alarm never before T_star ===")
    rng = np.random.default_rng(0)
    for _ in range(10):
        y = np.cumsum(rng.normal(size=150))
        out = monitor(y, r_star=0.5, minw=20, boundary="kurozumi")
        if not np.isnan(out.alarm[0]):
            assert out.alarm[0] >= out.t_star
    print("all alarms (if any) fired at/after T_star: OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
