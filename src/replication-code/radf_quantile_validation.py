"""Python cross-check of exuber's quantile_test() (Wu, Shi & Wu 2025's
quantile-based global test), mirroring radf_quantile_validation.R's
structural check that Q matches radf()'s own adf statistic, plus a
direct bit-for-bit cross-check of the deterministic point statistics
(tstat/tau/delta) against R.

Run standalone: uv run --project pyexuber python
docs/replication/alternative-paradigms/radf_quantile_validation.py

Reference values: a direct R run from the exuber-project/ root,
2026-09-16:
    Rscript -e 'Sys.setenv(NOT_CRAN="true");
    devtools::load_all("exuber", quiet=TRUE); set.seed(42);
    y <- cumsum(rnorm(80));
    dput(exuber:::quantile_adf_tstat(y));
    qcd <- exuber:::quantile_check_density(diff(y), 0.5);
    dput(qcd$b_tau); dput(qcd$f_hat);
    qt <- quantile_test(y, tau=0.5, nrep=50, sig_lvl=95, seed=7);
    dput(unname(qt$tstat)); dput(round(unname(qt$delta), 10));
    qt_opt <- quantile_test(y, tau="optimal", nrep=50, sig_lvl=95, seed=7);
    dput(unname(qt_opt$tau)); dput(unname(qt_opt$tstat))'

Note on tolerance: pyexuber's QR fit uses an IRLS solver (see monitor.py's
module docstring) rather than R's quantreg::rq() simplex/interior-point
solver, so tstat matches to ~1e-6, not to machine precision the way the
fully closed-form pieces (badf/cusum/lbi) do elsewhere in this family.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.quantile_test import (  # noqa: E402
    _quantile_adf_tstat,
    _quantile_check_density,
    quantile_test,
)

# set.seed(42); y <- cumsum(rnorm(80)) -- same series as
# docs/replication/monitoring/radf_monitor_kurozumi_boundary_validation.py.
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
    print("=== 1. quantile_adf_tstat / quantile_check_density -- cross-check vs R ===")
    tstat = _quantile_adf_tstat(Y42)
    np.testing.assert_allclose(tstat, -1.66093811600216, atol=1e-6)
    print(f"quantile_adf_tstat(y) = {tstat} (expect -1.66093811600216)")

    b_tau, f_hat = _quantile_check_density(np.diff(Y42), 0.5)
    np.testing.assert_allclose(b_tau, 0.0898328865790818, atol=1e-8)
    np.testing.assert_allclose(f_hat, 0.369113058891683, atol=1e-6)
    print(f"b_tau={b_tau}, f_hat={f_hat}")

    print("\n=== 2. quantile_test(tau=0.5) -- cross-check vs R ===")
    qt = quantile_test(Y42, tau=0.5, nrep=50, sig_lvl=95, seed=7)
    np.testing.assert_allclose(qt.tstat[0], -1.1087521367351, atol=1e-4)
    np.testing.assert_allclose(qt.delta[0], 0.7832341878, atol=1e-6)
    print(f"tstat={qt.tstat[0]} (expect -1.1087521367351), delta={qt.delta[0]}")

    print("\n=== 3. quantile_test(tau='optimal') -- cross-check vs R ===")
    qt_opt = quantile_test(Y42, tau="optimal", nrep=50, sig_lvl=95, seed=7)
    assert qt_opt.tau[0] == 0.8
    np.testing.assert_allclose(qt_opt.tstat[0], -0.440485749709493, atol=1e-4)
    print(f"tau={qt_opt.tau[0]} (expect 0.8), tstat={qt_opt.tstat[0]}")

    print("\n=== 4. Empirical size under H0 (100 reps, tau=0.5) ===")
    rng = np.random.default_rng(11)
    rejections = []
    for _ in range(100):
        y_null = np.cumsum(rng.normal(size=100))
        res = quantile_test(y_null, tau=0.5, nrep=200, sig_lvl=95, seed=1)
        rejections.append(bool(res.detected[0]))
    fpr = np.mean(rejections)
    print(f"empirical size: {fpr:.3f} (nominal 0.05)")

    print("\n=== 5. Detection power under a genuine explosive alternative ===")
    rng = np.random.default_rng(12)
    detect = []
    for _ in range(30):
        n_bubble = 100
        e = rng.normal(size=n_bubble - 1)
        y_bubble = np.concatenate(([0.0], np.cumsum(1.03 ** np.arange(n_bubble - 1) + e)))
        res = quantile_test(y_bubble, tau=0.5, nrep=200, sig_lvl=95, seed=1)
        detect.append(bool(res.detected[0]))
    print(f"detection rate: {np.mean(detect):.3f}")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
