"""Python replication script for dating_knp() (Kejriwal, Nguyen & Perron
2025 bias-corrected single-bubble dating), cross-checking pyexuber's port
against radf_knp_validation.R in this same folder.

RNG note: numpy's Generator, not R's RNG -- independent seeds, same
qualitative checks (see rootstamp_validation.py's module docstring for
the same convention elsewhere in this port). dating_knp() is pure numpy
(no C++ extension needed), so this script WAS run directly on the dev
machine that wrote the port -- all numbers below are real local runs.

R's own script (docs/dating-and-root-inference.md, "Improved
retrospective dating"), reports for context (not
asserted bit-for-bit, different RNG): naive (omit=FALSE) mean|tau1-T1|
=39.0 vs mean|tau1-T2|=1.0 (Theorem 1: tau1_hat tracks the COLLAPSE
date, not origination); omission-corrected mean|tau1-T1|=13.0 (Theorem
2: a genuine correction, not full elimination at this finite T);
delta_hat mean 0.986 vs true 1.05.

Run standalone:
uv run --project pyexuber python
docs/replication/dating-and-root-inference/radf_knp_validation.py
"""

import math

import numpy as np

from exuber.dating_knp import _knp_dp, _knp_find_break, dating_knp


def _ols_ssr(xseg: np.ndarray, zseg: np.ndarray) -> float:
    a = np.vstack([xseg, np.ones_like(xseg)]).T
    coef, *_ = np.linalg.lstsq(a, zseg, rcond=None)
    return float(np.sum((zseg - a @ coef) ** 2))


def check_formula_exact() -> None:
    print("=== 1. Formula-exact: _knp_find_break() vs brute-force nested OLS ===")
    rng = np.random.default_rng(3)
    n = 26
    y = np.cumsum(rng.normal(size=n))
    n1 = n - 1
    x, z = y[:n1], np.diff(y)
    k_min = max(2, math.ceil(0.1 * n1))

    for omit in (False, True):
        tau1, tau2, ssr = _knp_find_break(y, trim=0.1, omit=omit)
        best_ssr, best = math.inf, None
        for t1 in range(k_min, n1 - 2 * k_min + 1):
            for t2 in range(t1 + k_min, n1 - k_min + 1):
                s = np.sum(z[:t1] ** 2) + _ols_ssr(x[t1:t2], z[t1:t2]) + np.sum(z[t2:n1] ** 2)
                if omit:
                    s -= z[t2] ** 2
                if s < best_ssr:
                    best_ssr, best = s, (t1, t2)
        print(f"  omit={omit}: vectorized=({tau1},{tau2},{ssr:.4f})  brute={best},{best_ssr:.4f}")
        assert (tau1, tau2) == best
        assert abs(ssr - best_ssr) < 1e-6


def _sim_knp(seed, t1=50, t2=90, t=200, delta=1.05):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    for i in range(1, t1):
        y[i] = y[i - 1] + rng.normal()
    for i in range(t1, t2):
        y[i] = delta * y[i - 1] + rng.normal()
    y[t2] = y[t1 - 1] + rng.normal()
    for i in range(t2 + 1, t):
        y[i] = y[i - 1] + rng.normal()
    return y, t1, t2


def check_theorem_1_and_2() -> None:
    print("\n=== 2. Reproducing Theorem 1 (naive inconsistency) vs Theorem 2")
    print("     (omission-corrected consistency) ===")
    print("(KNP's own DGP: unit root -> no-intercept explosive AR(1) -> an")
    print("instantaneous collapse back near the pre-bubble level -> fresh unit")
    print("root. Theorem 1 proves plain OLS's origination-date estimate")
    print("converges to the TRUE COLLAPSE date, not the true origination date;")
    print("Theorem 2 proves the single-residual omission fixes this.)\n")

    def run(seed, omit):
        y, t1, t2 = _sim_knp(seed)
        tau1, tau2, _ssr = _knp_find_break(y, trim=0.05, omit=omit)
        return tau1, tau2, t1, t2

    res_naive = [run(s, False) for s in range(30)]
    res_om = [run(s, True) for s in range(30)]

    bias_naive_t1 = np.mean([abs(tau1 - t1) for tau1, _tau2, t1, _t2 in res_naive])
    bias_naive_t2 = np.mean([abs(tau1 - t2) for tau1, _tau2, _t1, t2 in res_naive])
    bias_om_t1 = np.mean([abs(tau1 - t1) for tau1, _tau2, t1, _t2 in res_om])
    bias_om_t2 = np.mean([abs(tau2 - t2) for _tau1, tau2, _t1, t2 in res_om])

    print(
        f"  Naive (omit=False): mean|tau1-T1|={bias_naive_t1:.1f}  "
        f"mean|tau1-T2|={bias_naive_t2:.1f}  (tau1 should track T2, not T1)"
    )
    print(
        f"  Omission-corrected: mean|tau1-T1|={bias_om_t1:.2f}  "
        f"mean|tau2-T2|={bias_om_t2:.2f}"
    )
    assert bias_naive_t2 < bias_naive_t1
    assert bias_om_t1 < bias_naive_t1 / 2


def check_delta_accuracy() -> None:
    print("\n=== 3. delta_hat accuracy under omission correction (true delta=1.05) ===")
    deltas = []
    for s in range(30):
        y, _t1, _t2 = _sim_knp(s)
        out = dating_knp(y, trim=0.05, omit=True)
        deltas.append(out.delta[0])
    deltas = np.array(deltas)
    print(
        f"  mean delta_hat = {deltas.mean():.3f} (true = 1.05), "
        f"mean|bias| = {np.mean(np.abs(deltas - 1.05)):.3f}"
    )
    assert np.mean(np.abs(deltas - 1.05)) < 0.3


# set.seed(7); y <- round(cumsum(rnorm(40)), 8) -- shared with the SSU script
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


def check_multi_bubble_dp() -> None:
    """Section 3.2 DP: breaks=2 reproduces the single-bubble search; 3 and 4
    breaks reproduce R's dates and coefficients (Rscript: dating_knp(y,
    trim = 0.1, breaks = b) on Y_VEC)."""
    print("\n=== Multi-bubble dynamic programme ===")
    y = np.cumsum(np.random.default_rng(11).normal(size=80))
    tau, ssr = _knp_dp(y, 2, 0.05, True)
    t1, t2, fssr = _knp_find_break(y, 0.05, True)
    assert tau == [t1, t2] and abs(ssr - fssr) < 1e-10
    print(f"  breaks=2: DP {tau} = single-bubble search ({t1}, {t2})")
    r3 = dating_knp(Y_VEC, trim=0.1, breaks=3)
    r4 = dating_knp(Y_VEC, trim=0.1, breaks=4)
    assert list(r3.origination[:, 0]) == [9, 19] and np.isnan(r3.collapse[1, 0])
    assert list(r4.collapse[:, 0]) == [14, 23]
    np.testing.assert_allclose(r4.delta[:, 0], [0.8410590945, 1.0818531838], atol=1e-8)
    print(f"  breaks=4: origination {r4.origination[:, 0]}, collapse {r4.collapse[:, 0]} (R: 9 19 / 14 23)")


def main() -> None:
    check_formula_exact()
    check_theorem_1_and_2()
    check_delta_accuracy()
    check_multi_bubble_dp()
    print("\nAll dating_knp() checks passed.")


if __name__ == "__main__":
    main()
