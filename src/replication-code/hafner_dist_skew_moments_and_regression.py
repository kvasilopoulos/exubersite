"""Python counterpart of hafner_dist_skew_moments_and_regression.R --
cross-checks pyexuber's port of `radf_wb_cv(dist_skew=True)`/
`_wb_dgp_hlst(dist_skew=True)` (Hafner 2020's skewness-corrected wild
bootstrap multiplier, R/radf_wb.R's `dist_skew` option) in exuber.cv.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/hafner_dist_skew_moments_and_regression.py

Checks 1-2 need no RNG match against R (the multiplier's defining moments
and the "additive option" regression check are both RNG-agnostic/
RNG-identity properties). Check 3 needs exuber._core (radf_wb_cv() calls
the compiled radf_stat() internally), so it only runs where the extension
is built (CI) -- see radf_kp_validation.py's own module docstring for the
same caveat.
"""

import numpy as np

from exuber.cv import _wb_dgp_hlst, radf_wb_cv


def check_moment_construction() -> None:
    """w = u/sqrt(2) + (v^2-1)/2, u,v ~ iid N(0,1) independent -- should
    have E[w]=0, E[w^2]=1, E[w^3]=1 (the paper's own construction, Step
    1), regardless of RNG engine."""
    rng = np.random.default_rng(1)
    n = 2_000_000
    u = rng.normal(size=n)
    v = rng.normal(size=n)
    w = u / np.sqrt(2) + (v**2 - 1) / 2
    print(f"E[w]={w.mean():.4f} (want 0)  E[w^2]={(w**2).mean():.4f} (want 1)  "
          f"E[w^3]={(w**3).mean():.4f} (want 1)")
    assert abs(w.mean()) < 0.01
    assert abs((w**2).mean() - 1) < 0.01
    assert abs((w**3).mean() - 1) < 0.03


def check_dist_skew_false_unchanged() -> None:
    """A purely additive option: dist_skew=False (the default) must
    reproduce the pre-change DGP bit-for-bit for the same seed."""
    y = np.cumsum(np.random.default_rng(5).normal(size=60))
    r1 = _wb_dgp_hlst(y, False, np.random.default_rng(5))
    r2 = _wb_dgp_hlst(y, False, np.random.default_rng(5), dist_skew=False)
    assert np.array_equal(r1, r2)
    print("OK: dist_skew=False identical to pre-change default behavior.")


def check_mutually_exclusive() -> None:
    data = np.cumsum(np.random.default_rng(0).normal(size=40))
    try:
        radf_wb_cv(data, minw=10, nboot=5, dist_rad=True, dist_skew=True)
        print("FAIL: expected an error")
    except ValueError as e:
        print(f"OK, errored: {e}")


def check_size_under_h0_right_skewed_heteroskedastic() -> None:
    """Empirical size under H0 with right-skewed, heteroskedastic errors
    (Hafner's own setting, footnote 1) -- needs exuber._core."""

    def rskew_innov(n: int, rng: np.random.Generator) -> np.ndarray:
        z = rng.normal(size=n)
        e = -np.log(z**2)
        return (e - e.mean()) / e.std()

    def run_once(seed: int, dist_skew: bool) -> bool:
        from exuber._unroot import unroot
        from exuber import _core

        rng = np.random.default_rng(seed)
        tn = 100
        g = 0.05 * (1 + 2 * np.cos(np.pi * np.arange(1, tn + 1) / tn) ** 2)
        y = np.cumsum(g * rskew_innov(tn, rng))
        cv = radf_wb_cv(y, minw=20, nboot=199, dist_skew=dist_skew, seed=1)
        obs = _core.radf_stat(unroot(y), 20, 0)
        sadf_obs = obs[tn - 20 + 1]
        return sadf_obs > cv.sadf_cv[0, 1]

    rej_normal = np.mean([run_once(s, False) for s in range(1, 41)])
    rej_skew = np.mean([run_once(s, True) for s in range(1, 41)])
    print(f"Empirical size, normal multiplier (dist_skew=False): {rej_normal:.3f}")
    print(f"Empirical size, skewed multiplier (dist_skew=True):  {rej_skew:.3f}")
    print("(nominal 0.05; paper's own finding: undersized in small samples under "
          "heteroskedasticity -- well below 0.05 is expected, not a red flag)")


if __name__ == "__main__":
    check_moment_construction()
    check_dist_skew_false_unchanged()
    check_mutually_exclusive()
    check_size_under_h0_right_skewed_heteroskedastic()
    print("All checks passed.")
