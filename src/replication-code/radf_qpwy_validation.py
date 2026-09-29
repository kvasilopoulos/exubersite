"""Python cross-check of exuber's monitor_quantile() (Wu, Shi & Wu 2025's
QPWY recursive quantile monitoring extension), mirroring
radf_qpwy_validation.R's formula-exact check of the per-window QR
t-ratio and the boundary-vs-marginal-quantile structural check that
caught a real bug in the R implementation (see docs/alternative-
paradigms.md, "Implementation (QPWY)").

Run standalone: uv run --project pyexuber python
docs/replication/alternative-paradigms/radf_qpwy_validation.py

_qpwy_stat_path() needs no RNG and no radf()/C++ extension, so its check
runs fully offline and matches R bit-for-bit (within the IRLS-vs-simplex
QR-solver tolerance documented in monitor.py's module docstring). Since
2026-09-29 the boundary simulation needs no radf() either (Q and Z from
prefix sums), so sections 4-5 check it directly: the independent-BM term
Z must be a process over windows (the single-z bug), and QPSY's grid
contains QPWY's path.

Reference values: a direct R run from the exuber-project/ root,
2026-09-16:
    Rscript -e 'Sys.setenv(NOT_CRAN="true");
    devtools::load_all("exuber", quiet=TRUE); set.seed(42);
    y <- cumsum(rnorm(80)); minw <- 15;
    r_idx <- (minw + 1L):length(y);
    stat <- exuber:::qpwy_stat_path(y, 0.5, r_idx);
    dput(round(tail(stat, 5), 10)); dput(length(stat))'
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.monitor_quantile import (  # noqa: E402
    _qpsy_stat_path,
    _qpwy_stat_path,
    _quantile_boundary_sim,
)
from exuber.quantile_test import quantile_test  # noqa: E402

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
    print("=== 1. Formula-exact check: qpwy_stat_path -- cross-check vs R ===")
    minw = 15
    r_idx = np.arange(minw + 1, len(Y42) + 1)
    stat = _qpwy_stat_path(Y42, 0.5, r_idx)
    assert len(stat) == 65
    expected_tail = np.array(
        [-1.128241465, -1.2210267401, -1.2384527008, -1.2255226955, -1.1087521367]
    )
    np.testing.assert_allclose(stat[-5:], expected_tail, atol=1e-4)
    print(f"len={len(stat)}, tail={stat[-5:]}")

    print("\n=== 2. Structural check: Q_{0,r} at r=n equals quantile_test()'s tstat ===")
    # Corollary 2's own claim: QPWY_n(tau) (the full-window recursion
    # endpoint) is exactly quantile_test()'s own single-shot tstat.
    qt = quantile_test(Y42, tau=0.5, nrep=10, sig_lvl=95, seed=1)
    assert abs(stat[-1] - qt.tstat[0]) < 1e-6
    print(f"stat[-1]={stat[-1]}, quantile_test tstat={qt.tstat[0]}: match")

    print(
        "\n=== 3. Boundary-vs-marginal-quantile sanity (the bug the R port found "
        "and fixed) ==="
    )
    # Reproduce the shape of the bug from docs/alternative-paradigms.md: a
    # per-r MARGINAL quantile of simulated paths badly inflates the
    # false-alarm rate relative to a SUPREMUM-calibrated one. Demonstrated
    # here directly on simulated standard-normal paths.
    rng = np.random.default_rng(0)
    nrep, n_mon = 500, 40
    paths = rng.normal(size=(nrep, n_mon))  # stand-in for the Q_{0,r} paths
    marginal_boundary = np.quantile(paths, 0.95, axis=0)  # WRONG (per-r)
    sup_boundary = np.quantile(paths.max(axis=1), 0.95)  # RIGHT (supremum)

    fpr_marginal = np.mean(np.any(paths > marginal_boundary, axis=1))
    fpr_sup = np.mean(paths.max(axis=1) > sup_boundary)
    print(f"per-r marginal boundary false-alarm rate: {fpr_marginal:.3f} (badly inflated)")
    print(f"supremum-calibrated boundary false-alarm rate: {fpr_sup:.3f} (~nominal 0.05)")
    assert fpr_marginal > 3 * fpr_sup
    assert abs(fpr_sup - 0.05) < 0.03

    print("\n=== 4. The 2026-09-29 bug: Z is a process over windows, not one z per path ===")
    # delta = 0 leaves only Z: one shared z per path would make sup_r Z
    # exactly N(0,1), 95% quantile 1.645
    sim = _quantile_boundary_sim(150, 20, 400, np.array([0.0]), False, np.random.default_rng(3))
    q95 = np.quantile(sim[:, 0], 0.95)
    print(f"95% quantile of sup_r Z_r: {q95:.3f} (single-z construction: 1.645)")
    assert q95 > 2

    print("\n=== 5. QPSY: its grid contains QPWY's path ===")
    sy = _qpsy_stat_path(Y42, 0.5, r_idx, minw)
    assert np.all(sy >= stat - 1e-12) and abs(sy[0] - stat[0]) < 1e-12
    d = np.array([0.2, 0.8])
    wy = _quantile_boundary_sim(60, 12, 20, d, False, np.random.default_rng(5))
    sb = _quantile_boundary_sim(60, 12, 20, d, True, np.random.default_rng(5))
    assert np.all(sb >= wy - 1e-12)
    print(f"QPSY path >= QPWY path everywhere; max QPSY = {sy.max():.4f}")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
