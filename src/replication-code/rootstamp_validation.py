"""Python replication script for rootstamp() (root inference: Guo, Sun &
Wang 2019 normal-t CI + Phillips-Magdalinos 2007 Cauchy CI), cross-checking
pyexuber's port against rootstamp_validation.R in this same folder.

RNG note: pyexuber uses numpy's Generator, not R's RNG, so this does not
reproduce the R script's draws bit-for-bit (see cv.py's module docstring
for the same convention elsewhere in this port) -- it independently
re-derives the same qualitative findings with its own seeds. The
R-derived reference numbers hardcoded below come from:

    cd exuber-project && Rscript docs/replication/dating-and-root-inference/rootstamp_validation.R

R 4.6.1:
  - qt(0.95,1)=6.313752, qt(0.975,1)=12.7062, qt(0.995,1)=63.65674
  - coverage (rho=1.05, n=200, 800 reps, seed=24601): 0.94625

Run standalone: uv run --project pyexuber python
docs/replication/dating-and-root-inference/rootstamp_validation.py
"""

import math

import numpy as np

from exuber.rootstamp import rootstamp


def check_cauchy_percentiles() -> None:
    print("=== 1. Cauchy percentiles (Skrobotov 2023 review's footnote 17) ===")
    published = {"10%": 6.313752, "5%": 12.7062, "1%": 63.65674}
    for label, p in [("10%", 0.95), ("5%", 0.975), ("1%", 0.995)]:
        q = math.tan(math.pi * (p - 0.5))  # standard Cauchy quantile, == qt(p, df=1)
        print(f"  {label}: computed={q:.6f}  R's qt()={published[label]:.6f}")
        assert abs(q - published[label]) < 1e-4


def check_point_estimate() -> None:
    print("\n=== 2. Point estimate: super-consistency at rho=1.05, n=200 ===")
    rng = np.random.default_rng(8675309)
    n = 200
    y = np.zeros(n)
    e = rng.normal(size=n)
    for t in range(1, n):
        y[t] = 1.05 * y[t - 1] + e[t]
    fit = rootstamp(y)
    print(f"  rho_hat={fit.rho:.4f} vs rho_true=1.0500")
    print(f"  95% CI: [{fit.rho_ci[0]:.6f}, {fit.rho_ci[1]:.6f}]")
    assert abs(fit.rho - 1.05) < 0.001


def check_coverage() -> None:
    print("\n=== 3. Coverage (rho=1.05, n=200, 800 reps, seed=24601) ===")
    print("R's own run (rootstamp_validation.R, different RNG): 0.94625")
    rng = np.random.default_rng(24601)
    n = 200
    covered = 0
    for _ in range(800):
        y = np.zeros(n)
        e = rng.normal(size=n)
        for t in range(1, n):
            y[t] = 1.05 * y[t - 1] + e[t]
        ci = rootstamp(y)
        if ci.rho_ci[0] <= 1.05 <= ci.rho_ci[1]:
            covered += 1
    coverage = covered / 800
    print(f"  Python coverage: {coverage}")
    # loose band, not an exact match -- different RNG, same qualitative
    # finite-sample-undercoverage finding documented in
    # docs/dating-and-root-inference.md's "Root inference" section.
    assert 0.85 < coverage < 1.0


def check_cauchy_formula_exact() -> None:
    print("\n=== 4. Cauchy-type CI (eq. 27, Phillips-Magdalinos 2007) formula-exact check ===")
    rng = np.random.default_rng(11)
    n = 150
    y = np.zeros(n)
    e = rng.normal(size=n)
    for t in range(1, n):
        y[t] = 1.04 * y[t - 1] + e[t]
    ci = rootstamp(y, type="cauchy", sig_lvl=95)
    q = math.tan(math.pi * (0.975 - 0.5))
    half_width_formula = q * (ci.rho**2 - 1) / ci.rho**ci.n
    half_width_fn = ci.rho_ci[1] - ci.rho
    print(f"  rootstamp(type='cauchy') half-width: {half_width_fn}")
    print(f"  eq. 27 formula half-width:           {half_width_formula}")
    assert abs(half_width_fn - half_width_formula) < 1e-10


def main() -> None:
    check_cauchy_percentiles()
    check_point_estimate()
    check_coverage()
    check_cauchy_formula_exact()
    print("\nAll rootstamp() checks passed.")


if __name__ == "__main__":
    main()
