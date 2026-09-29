"""Python counterpart of hafner_dist_skew_power_and_size.R -- power and
size checks for pyexuber's port of `radf_wb_cv(dist_skew=True)` (Hafner
2020's skewness-corrected wild bootstrap, exuber.cv). See
hafner_dist_skew_moments_and_regression.py for the moment-construction
and additive-option checks.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/hafner_dist_skew_power_and_size.py

Both checks below need exuber._core (radf() and radf_wb_cv() call the
compiled radf_stat() internally), so they only run where the extension
is built (CI).
"""

import numpy as np

from exuber.cv import radf_wb_cv
from exuber.radf import radf


def check_power_ordinary_normal_innovations() -> None:
    """Power check with dist_skew=True bootstrap, ordinary (non-skewed)
    normal innovations -- confirms the skewed multiplier doesn't harm
    basic detection ability."""

    def run_once(seed: int) -> bool:
        rng = np.random.default_rng(seed)
        tn, te = 100, 60
        normal_part = np.cumsum(rng.normal(size=te))
        expl_len = tn - te
        expl_part = normal_part[-1] * 1.06 ** np.arange(1, expl_len + 1) + np.cumsum(
            rng.normal(scale=0.3, size=expl_len)
        )
        y = np.concatenate([normal_part, expl_part])
        obs = radf(y, minw=20).sadf[0]
        cv = radf_wb_cv(y, minw=20, nboot=199, dist_skew=True, seed=1)
        return obs > cv.sadf_cv[0, 1]

    power = np.mean([run_once(s) for s in range(1, 31)])
    print(f"Empirical power (30 reps): {power:.3f}")


def check_size_under_h0_paper_own_distribution() -> None:
    """Size under H0 with the paper's own right-skewed innovation
    distribution (negative log-chi-square(1)), no heteroskedasticity
    added on top."""

    def rskew_innov(n: int, rng: np.random.Generator) -> np.ndarray:
        z = rng.normal(size=n)
        e = -np.log(z**2)
        return (e - e.mean()) / e.std()

    def run_size(seed: int, dist_skew: bool) -> bool:
        rng = np.random.default_rng(seed)
        tn = 150
        y = np.cumsum(rskew_innov(tn, rng))
        obs = radf(y, minw=20).sadf[0]
        cv = radf_wb_cv(y, minw=20, nboot=199, dist_skew=dist_skew, seed=1)
        return obs > cv.sadf_cv[0, 1]

    rej_normal = np.mean([run_size(s, False) for s in range(1, 61)])
    rej_skew = np.mean([run_size(s, True) for s in range(1, 61)])
    print(f"Empirical size, normal multiplier: {rej_normal:.3f}")
    print(f"Empirical size, skewed multiplier: {rej_skew:.3f}")
    print("(nominal 0.05)")


if __name__ == "__main__":
    check_power_ordinary_normal_innovations()
    check_size_under_h0_paper_own_distribution()
    print("All checks passed.")
