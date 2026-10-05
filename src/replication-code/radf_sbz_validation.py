"""Python counterpart of radf_sbz_validation.R -- cross-checks pyexuber's
port of radf_sbz()/radf_sbz_cv()/radf_sbz_union() (SBZ: WLS + kernel
volatility, Harvey, Leybourne & Zu 2019, R/radf_sbz.R) in
exuber.volatility.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_sbz_validation.py

Part 1 checks wls_dfstat_grid()/radf_sbz() bit-for-bit against R, fed the
SAME deterministic input series used by radf_tt_validation.py (no RNG
involved at the formula level):

    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    set.seed(7); y <- round(cumsum(rnorm(40)), 8)
    vol <- exuber:::kernel_spot_vol(y, kernel = "gaussian")
    exuber:::wls_dfstat_grid(y, vol$sigma2, 10)
    radf_sbz(y, minw = 10)

Part 2 reproduces R's own H0 empirical-size check (see
docs/volatility-robustness.md, "SBZ"): a
pure random walk should reject at roughly the nominal rate, not grossly
more. Unlike radf_wb_cv(), radf_sbz()/radf_sbz_cv()'s own wild-bootstrap
DGP is pure Python (no exuber._core call), so this part runs everywhere.
radf_sbz_union() itself needs exuber._core (supDF via the compiled
radf_stat()), so its own check is CI-only.
"""

import numpy as np

from exuber._kernel_vol import kernel_spot_vol
from exuber.radf_sbz import radf_sbz, radf_sbz_cv, radf_sbz_union, wls_dfstat_grid

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

R_WLS_SADF = 1.6435722442
R_WLS_GSADF = 1.7667532509
R_SBZ_ADF = 0.0608511749
R_SBZ_SADF = 1.6435722442
R_SBZ_GSADF = 1.7667532509


def check_wls_dfstat_grid_matches_r() -> None:
    sigma2, _ = kernel_spot_vol(Y_VEC, kernel="gaussian")
    res = wls_dfstat_grid(Y_VEC, sigma2, MINW)
    np.testing.assert_allclose(res["sadf"], R_WLS_SADF, atol=1e-6)
    np.testing.assert_allclose(res["gsadf"], R_WLS_GSADF, atol=1e-6)
    print("wls_dfstat_grid(): matches R's exuber:::wls_dfstat_grid() to 1e-6.")


def check_radf_sbz_matches_r() -> None:
    res = radf_sbz(Y_VEC, minw=MINW)
    np.testing.assert_allclose(res.adf[0], R_SBZ_ADF, atol=1e-6)
    np.testing.assert_allclose(res.sadf[0], R_SBZ_SADF, atol=1e-6)
    np.testing.assert_allclose(res.gsadf[0], R_SBZ_GSADF, atol=1e-6)
    print("radf_sbz(): matches R's radf_sbz() to 1e-6.")


def check_empirical_size_under_h0() -> None:
    """Empirical rejection rate under H0 (pure random walk) should be close
    to nominal 0.05, not grossly oversized."""
    rng = np.random.default_rng(13579)
    n, nrep, nboot = 100, 60, 150
    rej_sbz = 0
    for _ in range(nrep):
        y = np.cumsum(rng.normal(size=n))
        res = radf_sbz(y, minw=15)
        cv = radf_sbz_cv(y, minw=15, nboot=nboot, seed=int(rng.integers(1_000_000_000)))
        if res.sadf[0] > cv.sadf_cv[0, 1]:
            rej_sbz += 1
    rate = rej_sbz / nrep
    print(f"supBZ empirical rejection rate under H0 (n={n}, nrep={nrep}, nboot={nboot}): {rate:.3f}")
    print("Target ~0.05 nominal.")
    assert rate < 0.25, f"rate={rate} is grossly oversized"


def check_radf_sbz_union_runs() -> None:
    """Needs exuber._core (supDF via the compiled radf_stat())."""
    rng = np.random.default_rng(1)
    y = np.cumsum(rng.normal(size=100))
    res = radf_sbz_union(y, minw=15, nboot=150, seed=1)
    print(f"supDF={res.supDF[0]:.4f} supBZ={res.supBZ[0]:.4f} U={res.U[0]:.4f} "
          f"p_supDF={res.p_supDF[0]:.3f} p_supBZ={res.p_supBZ[0]:.3f} p_U={res.p_U[0]:.3f}")
    assert res.U[0] >= res.supDF[0]  # U := max(supDF, ratio*supBZ)


if __name__ == "__main__":
    check_wls_dfstat_grid_matches_r()
    check_radf_sbz_matches_r()
    check_empirical_size_under_h0()
    check_radf_sbz_union_runs()
    print("All checks passed.")
