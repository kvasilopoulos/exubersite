"""Python counterpart of radf_ssu_validation.R -- cross-checks pyexuber's
port of ssu_test() (stochastic explosive-coefficient SSU test, Kurozumi &
Nishi 2025, minimum-viable subset, R/ssu_test.R) in exuber.volatility.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_ssu_validation.py

Part 1 checks ssu_stat_path()/ssu_test() bit-for-bit against R, fed the
SAME deterministic input series used by radf_tt_validation.py (no RNG
involved at the formula level):

    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    set.seed(7); y <- round(cumsum(rnorm(40)), 8)
    ps <- exuber:::ssu_prefix_sums(y)
    exuber:::ssu_stat_path(ps, 10:39)

Part 2 also cross-checks ssu_stat_path()'s bilinear cross-moment
expansion against a from-scratch brute-force computation (two separately
fitted OLS regressions plus a manual residual cross-moment) at several
window sizes, independent of both R and the closed-form implementation --
the same style of check test-ssu.R's own R validation uses.

Part 3 reproduces the Table I lookup checks and a detection-power check
on Kurozumi & Nishi's own eq. 2 style stochastic-explosive-coefficient DGP
(their own alternative, own RNG).
"""

import numpy as np

from exuber.ssu_test import ssu_prefix_sums, ssu_q, ssu_stat_path, ssu_test

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

R_SSU_STAT = np.array(
    [
        0.3507221612, 0.3706956340, 0.8058323864, -0.6875656648, -0.5776884647,
        -1.2064421069, -0.5760735436, -0.7290714153, -0.9597233521, -1.2301899913,
        -1.5466598868, -1.8660594834, -2.0114222483, -0.6075219402, -0.7668093101,
        -1.0384324273, -1.3022337320, -1.5645114773, -1.1043525912, -1.1980506482,
        -1.0578407671, -1.2097931568, -1.3608122755, -1.4783438471, -1.5161782631,
        -1.5230627989, -1.6201706488, -1.4294408532, -1.4956945921, -1.5827171047,
    ]
)


def check_ssu_stat_path_matches_r() -> None:
    ps = ssu_prefix_sums(Y_VEC)
    hi_idx = np.arange(MINW, len(Y_VEC))
    stat = ssu_stat_path(ps, hi_idx)
    np.testing.assert_allclose(stat, R_SSU_STAT, atol=1e-6)
    print("ssu_stat_path(): matches R's exuber:::ssu_stat_path() to 1e-6.")


def check_ssu_test_matches_r() -> None:
    res = ssu_test(Y_VEC, minw=MINW, sig_lvl=95)
    np.testing.assert_allclose(res.sadf[0], R_SSU_STAT.max(), atol=1e-6)
    print(f"ssu_test(): sadf={res.sadf[0]:.10f} matches max(R_SSU_STAT)={R_SSU_STAT.max():.10f}.")


def check_formula_vs_brute_force() -> None:
    """Independent of R: ssu_stat_path()'s bilinear expansion against a
    from-scratch computation (two OLS fits + manual cross-moment)."""
    rng = np.random.default_rng(2)
    n = 150
    y = np.cumsum(rng.normal(size=n))
    ps = ssu_prefix_sums(y)

    def brute_force_stat(hi: int) -> float:
        win = np.arange(hi)
        x1 = y[win]
        d1 = y[win + 1] - x1
        x2 = x1**2
        d2 = d1**2

        a1 = np.column_stack([np.ones(hi), x1])
        beta1, *_ = np.linalg.lstsq(a1, d1, rcond=None)
        eps_hat = d1 - a1 @ beta1

        a2 = np.column_stack([np.ones(hi), x2])
        beta2, *_ = np.linalg.lstsq(a2, d2, rcond=None)
        eta_hat = d2 - a2 @ beta2

        sigma2_eps = np.sum(eps_hat**2) / (hi - 2)
        sigma2_eta = np.sum(eta_hat**2) / (hi - 2)
        sigma2_epseta = np.sum(eps_hat * eta_hat) / (hi - 1)
        sigma_eps, sigma_eta = np.sqrt(sigma2_eps), np.sqrt(sigma2_eta)
        psi_hat = sigma2_epseta / (sigma_eps * sigma_eta)

        omega_hat = beta2[1]
        sxx2_c = np.sum((x2 - x2.mean()) ** 2)
        t_omega = omega_hat / np.sqrt(sigma2_eta / sxx2_c)

        num_corr = np.sum((x2 - x2.mean()) * d1)
        den_corr = np.sqrt(sxx2_c)
        correction = (psi_hat / sigma_eps) * num_corr / den_corr
        return (t_omega - correction) / np.sqrt(1 - psi_hat**2)

    for hi in (50, 80, 120, 149):
        fast = ssu_stat_path(ps, np.array([hi]))[0]
        manual = brute_force_stat(hi)
        print(f"hi={hi} fast={fast:.8f} manual={manual:.8f} |diff|={abs(fast - manual):.2e}")
        assert abs(fast - manual) < 1e-6


def check_table_lookup() -> None:
    for lv, expect in ((90, 2.90), (95, 3.30), (99, 4.20)):
        assert ssu_q(lv) == expect
        print(f"sig_lvl={lv} -> crit={ssu_q(lv):.2f}")
    try:
        ssu_q(80)
        print("FAIL: expected an error")
    except ValueError as e:
        print(f"untabulated sig_lvl=80: OK, errored: {e}")


def check_power_on_stochastic_coefficient_dgp() -> None:
    """Detection power on Kurozumi & Nishi's own eq. 2 style
    stochastic-explosive-coefficient alternative."""

    def make_stochastic_bubble(rng: np.random.Generator, n: int, te_frac: float = 0.5,
                                c1: float = 3.0, a: float = 4.0) -> np.ndarray:
        y = np.empty(n)
        y[0] = rng.normal()
        te = round(te_frac * n)
        for t in range(1, n):
            if t < te:
                y[t] = y[t - 1] + rng.normal()
            else:
                rho_t = 1 + c1 / n + a * rng.normal() / np.sqrt(n)
                y[t] = rho_t * y[t - 1] + rng.normal()
        return y

    rng = np.random.default_rng(2)
    n, nrep = 200, 60
    detected = sum(
        ssu_test(make_stochastic_bubble(rng, n), sig_lvl=95).detected[0] for _ in range(nrep)
    )
    power = detected / nrep
    print(f"SSU power on stochastic-coefficient bubble DGP: {power:.3f}")


if __name__ == "__main__":
    check_ssu_stat_path_matches_r()
    check_ssu_test_matches_r()
    check_formula_vs_brute_force()
    check_table_lookup()
    check_power_on_stochastic_coefficient_dgp()
    print("All checks passed.")
