"""Python cross-check of exuber's lbi_test() (Breitung & Diegel 2025's
static locally best invariant test), mirroring
radf_lbi_validation.R's telescoping-identity, null-distribution, and
basic-run checks.

Run standalone: uv run --project pyexuber python
docs/replication/monitoring/radf_lbi_validation.py

Reference values: a direct R run from the exuber-project/ root,
    Rscript -e 'Sys.setenv(NOT_CRAN="true");
    devtools::load_all("exuber", quiet=TRUE); set.seed(42);
    y <- cumsum(rnorm(80)); res <- lbi_test(y, sig_lvl=95);
    dput(unname(res$stat)); dput(res$crit)'
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "pyexuber" / "src"))

from exuber.lbi_test import lbi_test  # noqa: E402

from radf_monitor_kurozumi_boundary_validation import Y42  # noqa: E402


def main() -> None:
    print("=== 1. Telescoping identity (eq. 4) ===")
    rng = np.random.default_rng(2)
    y = np.cumsum(rng.normal(size=60))
    n = len(y)
    dy = np.diff(y)
    ylag = y[: n - 1]
    lhs = 2 * np.sum(dy * ylag)
    sigma2_tilde = np.mean(dy**2)
    rhs = y[-1] ** 2 - y[0] ** 2 - (n - 1) * sigma2_tilde
    np.testing.assert_allclose(lhs, rhs, rtol=1e-10)
    print(f"2*sum(dy*ylag)={lhs:.6f}, y_T^2 - y_0^2 - (n-1)*sigma_tilde^2={rhs:.6f}: OK")

    print("\n=== 2. Basic run -- cross-check vs R ===")
    res = lbi_test(Y42, sig_lvl=95)
    np.testing.assert_allclose(res.stat[0], 0.025526221649211, atol=1e-8)
    np.testing.assert_allclose(res.crit, 1.64485362695147, atol=1e-10)
    print(f"stat={res.stat[0]}, crit={res.crit}, detected={res.detected[0]}")
    assert not res.detected[0]

    print("\n=== 3. Empirical null distribution is standard normal (500 reps) ===")
    rng = np.random.default_rng(3)
    stats = np.array([lbi_test(np.cumsum(rng.normal(size=100))).stat[0] for _ in range(500)])
    print(f"mean={stats.mean():.4f} (theory 0), sd={stats.std():.4f} (theory 1)")
    assert abs(stats.mean()) < 0.15
    assert abs(stats.std() - 1) < 0.15
    fpr = np.mean(stats > 1.64485362695147)
    print(f"empirical false-alarm rate at 95%: {fpr:.3f} (nominal 0.05)")

    print("\n=== 4. Detection power under a genuine explosive alternative ===")
    rng = np.random.default_rng(5)
    detect = []
    for _ in range(60):
        n_bubble = 60
        e = rng.normal(size=n_bubble - 1)
        y_bubble = np.concatenate(([0.0], np.cumsum(1.05 ** np.arange(n_bubble - 1) + e)))
        detect.append(lbi_test(y_bubble, sig_lvl=95).detected[0])
    print(f"detection rate: {np.mean(detect):.3f}")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
