"""Python counterpart of sign_based_finite_T_crosscheck.R -- cross-checks
pyexuber's port of radf_sign()/radf_sign_cv() (sign-based sGSADF, Harvey,
Leybourne & Zu 2020, R/radf_sign.R) in exuber.volatility.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/sign_based_finite_T_crosscheck.py

Part 1 checks gls_dfstat_grid(sign_transform(y))/gls_dfstat_grid(
sign_demean_transform(y)) bit-for-bit against R on the SAME deterministic
input series used by radf_tt_validation.py (no RNG at the formula level):

    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    set.seed(7); y <- round(cumsum(rnorm(40)), 8)
    exuber:::gls_dfstat_grid(exuber:::sign_transform(y), 10)
    exuber:::gls_dfstat_grid(exuber:::sign_demean_transform(y), 10)

Part 2 checks radf_sign_cv() against the paper's own Table 1 at the
EXACT finite T = 200 (not the T = Inf asymptotic row): the paper's own
text documents sPSY's finite-sample critical values converging to the
asymptotic limit much more slowly than sPWY's, so comparing at the exact
matching finite T avoids that convergence ambiguity (same reasoning
sign_based_finite_T_crosscheck.R uses). Like radf_tt_validation.py's
Whitehouse check, this is a real RNG-agnostic cross-check against a fixed
published target, not just a shape/sanity check.
"""

import numpy as np

from exuber._gls_dfstat import gls_dfstat_grid
from exuber.radf_sign import (
    radf_sign,
    radf_sign_cv,
    radf_sign_dm,
    sign_demean_transform,
    sign_transform,
)

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

R_SIGN_SADF = 0.6739661980
R_SIGN_GSADF = 2.0040132787
R_SIGNDM_SADF = 3.6784167432
R_SIGNDM_GSADF = 3.6784167432


def check_radf_sign_matches_r() -> None:
    res = radf_sign(Y_VEC, minw=MINW)
    np.testing.assert_allclose(res.sadf[0], R_SIGN_SADF, atol=1e-6)
    np.testing.assert_allclose(res.gsadf[0], R_SIGN_GSADF, atol=1e-6)
    print("radf_sign(): matches R's radf_sign() to 1e-6.")


def check_radf_sign_dm_matches_r() -> None:
    res = radf_sign_dm(Y_VEC, minw=MINW)
    np.testing.assert_allclose(res.sadf[0], R_SIGNDM_SADF, atol=1e-6)
    np.testing.assert_allclose(res.gsadf[0], R_SIGNDM_GSADF, atol=1e-6)
    print("radf_sign_dm(): matches R's radf_sign_dm() to 1e-6.")


def check_exact_invariance_to_heteroskedasticity() -> None:
    """The paper's central claim: identical sadf/gsadf under wildly
    different volatility scaling of the SAME sign pattern."""
    rng = np.random.default_rng(7)
    n2, te = 150, 90
    base_incr = rng.normal(size=te)
    expl_incr = np.concatenate([[rng.normal(loc=3)], rng.normal(size=n2 - te - 1)])
    raw_dy = np.concatenate([base_incr, expl_incr])
    y_homo = np.cumsum(raw_dy)
    vol_pattern = np.concatenate([np.full(40, 0.1), np.full(60, 10.0), np.full(n2 - 100, 1.0)])
    y_hetero = np.cumsum(raw_dy * vol_pattern)

    r_homo = radf_sign(y_homo, minw=20)
    r_hetero = radf_sign(y_hetero, minw=20)
    print(f"homoskedastic:   sadf={r_homo.sadf[0]:.6f} gsadf={r_homo.gsadf[0]:.6f}")
    print(f"heteroskedastic: sadf={r_hetero.sadf[0]:.6f} gsadf={r_hetero.gsadf[0]:.6f}")
    np.testing.assert_allclose(r_homo.sadf, r_hetero.sadf)
    np.testing.assert_allclose(r_homo.gsadf, r_hetero.gsadf)
    print("Exact invariance confirmed: bit-identical sadf/gsadf under a wildly different volatility path.")


def check_radf_sign_cv_against_published_table1_t200() -> None:
    cv = radf_sign_cv(n=200, minw=20, nrep=1500, seed=1)
    published_sadf = np.array([2.405, 2.735, 3.434])  # sPWY, T=200
    published_gsadf = np.array([3.469, 3.901, 4.957])  # sPSY, T=200
    print(f"Simulated sadf_cv (sPWY):  {cv.sadf_cv} vs published {published_sadf}")
    print(f"Simulated gsadf_cv (sPSY): {cv.gsadf_cv} vs published {published_gsadf}")
    assert np.all(np.abs(cv.sadf_cv - published_sadf) < 0.4)
    assert np.all(np.abs(cv.gsadf_cv - published_gsadf) < 0.6)


def check_power() -> None:
    """Empirical power on a clear mildly explosive alternative, using
    radf_sign_cv()'s own simulated critical value."""
    rng = np.random.default_rng(2)
    cv = radf_sign_cv(n=150, minw=20, nrep=500, seed=2)

    def run_once() -> bool:
        te = 90
        normal_part = np.cumsum(rng.normal(size=te))
        expl_len = 150 - te
        expl_part = normal_part[-1] * 1.05 ** np.arange(1, expl_len + 1) + np.cumsum(
            rng.normal(scale=0.5, size=expl_len)
        )
        y = np.concatenate([normal_part, expl_part])
        return radf_sign(y, minw=20).gsadf[0] > cv.gsadf_cv[1]

    power = np.mean([run_once() for _ in range(30)])
    print(f"Empirical power (30 reps, gsadf vs simulated 95% cv at n=150): {power:.3f}")


if __name__ == "__main__":
    check_radf_sign_matches_r()
    check_radf_sign_dm_matches_r()
    check_exact_invariance_to_heteroskedasticity()
    check_radf_sign_cv_against_published_table1_t200()
    check_power()
    print("All checks passed.")
