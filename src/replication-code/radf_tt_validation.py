"""Python counterpart of radf_tt_validation.R -- cross-checks pyexuber's
port of radf_tt()/radf_tt_cv() (STADF/GSTADF, Kurozumi, Skrobotov & Tsarev
2024, R/radf_tt.R) in exuber.volatility.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_tt_validation.py

Part 1 checks gls_dfstat_grid()/variance_profile()/radf_tt() bit-for-bit
against R, fed the SAME deterministic input series (no RNG involved at the
formula level, so R and numpy agree exactly) -- produced with:

    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    set.seed(7); y <- round(cumsum(rnorm(40)), 8)
    exuber:::gls_dfstat_grid(y, 10)
    exuber:::variance_profile(y, kernel = "uniform")
    radf_tt(y, minw = 10, kernel = "uniform")

Part 2 checks radf_tt_cv() (Monte Carlo, numpy's own RNG -- see sim.py's
module docstring for why this can't be bit-exact against R) against
Whitehouse (2019)'s published STADF triple, quoted in the paper's footnote
4 for minw/n = 0.1: (2.319, 2.626, 3.223). Since this is the fixed
asymptotic target of a pivotal simulation (not an R-RNG-specific value),
a wide-tolerance comparison is a real, RNG-agnostic cross-check, not just
a shape/sanity check.
"""

import numpy as np

from exuber._gls_dfstat import gls_dfstat_grid
from exuber.radf_tt import radf_tt, radf_tt_cv, variance_profile

MINW = 10

Y_VEC = np.array(
    [
        2.28724716, 1.09047548, 0.39618297, -0.01610998, -0.98678332, -1.93406327,
        -1.18592393, -1.30287915, -1.15022153, 1.03975658, 1.39674281, 4.11349459,
        6.39494652, 6.71896706, 8.61503413, 9.08271464, 8.18891391, 7.88158561,
        7.87676319, 8.86492734, 9.7046777, 10.41001953, 11.71598425, 10.32798804,
        11.6009049, 11.78509767, 12.53737757, 13.12912262, 12.14607002, 11.87000607,
        10.99915505, 11.7178656, 11.82851848, 11.75005171, 11.32956125, 10.76743537,
        11.76494882, 10.65981876, 10.51753093, 10.83252583,
    ]
)

R_GLS_BADF = np.array(
    [
        -0.4869476194, -0.5993680238, -0.2024272131, -0.0927192182, 0.4725944183,
        0.5988636710, 0.2121796002, 0.1113111528, 0.1065742445, 0.3591482854,
        0.5876420656, 0.7788759548, 1.1488183477, 0.5616782653, 0.8780267990,
        0.8944064816, 1.0772802859, 1.2107378397, 0.8360708791, 0.7377989079,
        0.4970657730, 0.6456036924, 0.6603799133, 0.6309975992, 0.5245109634,
        0.3967069003, 0.5846610107, 0.3332265155, 0.3047853676, 0.3604560044,
    ]
)
R_GLS_SADF = 1.2107378397
R_GLS_GSADF = 1.9637889135

R_VP_ETAGRID = np.array(
    [
        0.0, 0.0446045107, 0.0644273035, 0.0713730614, 0.1025117707, 0.1199565554,
        0.1522699431, 0.1523421936, 0.1534211678, 0.3066249818, 0.3113364130,
        0.3113364130, 0.4621388947, 0.4623333174, 0.5629110394, 0.5650758544,
        0.6055428527, 0.6150329465, 0.6170647429, 0.6362261844, 0.6488411132,
        0.6577037608, 0.6964175182, 0.7756513090, 0.8175841271, 0.8179259903,
        0.8317696672, 0.8375285024, 0.8734079342, 0.8764756454, 0.9011078597,
        0.9177485559, 0.9189042776, 0.9191068751, 0.9229030275, 0.9299740903,
        0.9702629790, 0.9942297001, 0.9942346351, 1.0,
    ]
)
R_VP_OMEGA2 = 0.8233400111

R_TT_ADF = 0.3401194521
R_TT_SADF = 1.1571296746
R_TT_GSADF = 1.8910439255


def check_gls_dfstat_grid_matches_r() -> None:
    res = gls_dfstat_grid(Y_VEC, MINW)
    np.testing.assert_allclose(res["badf"], R_GLS_BADF, atol=1e-6)
    np.testing.assert_allclose(res["sadf"], R_GLS_SADF, atol=1e-6)
    np.testing.assert_allclose(res["gsadf"], R_GLS_GSADF, atol=1e-6)
    print("gls_dfstat_grid(): matches R's exuber:::gls_dfstat_grid() to 1e-6.")


def check_variance_profile_matches_r() -> None:
    eta_grid, omega2 = variance_profile(Y_VEC, kernel="uniform")
    np.testing.assert_allclose(eta_grid, R_VP_ETAGRID, atol=1e-6)
    np.testing.assert_allclose(omega2, R_VP_OMEGA2, atol=1e-6)
    print("variance_profile(): matches R's exuber:::variance_profile() to 1e-6.")


def check_radf_tt_matches_r() -> None:
    res = radf_tt(Y_VEC, minw=MINW, kernel="uniform")
    np.testing.assert_allclose(res.adf[0], R_TT_ADF, atol=1e-6)
    np.testing.assert_allclose(res.sadf[0], R_TT_SADF, atol=1e-6)
    np.testing.assert_allclose(res.gsadf[0], R_TT_GSADF, atol=1e-6)
    print("radf_tt(): matches R's radf_tt() to 1e-6.")


def check_radf_tt_cv_against_published_whitehouse() -> None:
    """radf_tt_cv()'s sadf_cv is a Monte Carlo estimate of a pivotal
    (RNG-agnostic) asymptotic target -- comparable directly to the
    published triple even though numpy's RNG differs from R's."""
    cv = radf_tt_cv(n=300, minw=30, nrep=1500, seed=555)
    published = np.array([2.319, 2.626, 3.223])
    diff = np.abs(cv.sadf_cv - published)
    assert np.all(diff < 0.35), f"sadf_cv={cv.sadf_cv} too far from published {published}"
    print(f"radf_tt_cv() sadf_cv={cv.sadf_cv} vs published {published} (max diff {diff.max():.3f}).")


def check_radf_tt_cv_badf_cv_identity() -> None:
    """badf_cv's last row must be bit-identical to adf_cv -- adf is
    literally badf's last point per replicate (see radf_tt.R's own note)."""
    cv = radf_tt_cv(n=60, minw=15, nrep=200, seed=1)
    np.testing.assert_allclose(cv.badf_cv[-1], cv.adf_cv)
    print("radf_tt_cv(): badf_cv's last row matches adf_cv exactly.")


if __name__ == "__main__":
    check_gls_dfstat_grid_matches_r()
    check_variance_profile_matches_r()
    check_radf_tt_matches_r()
    check_radf_tt_cv_against_published_whitehouse()
    check_radf_tt_cv_badf_cv_identity()
    print("All checks passed.")
