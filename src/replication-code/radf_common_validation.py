"""Validation of radf_common()/radf_common_cv() -- Chen, Phillips & Shi
(2023)'s common-bubble detection via PCA + PSY. Same folder/base name as
the R script it cross-checks: radf_common_validation.R. See
docs/multivariate.md, "Common-bubble detection via PCA + PSY", for the
full narrative: Theorem
4.3's claim that the PSY-on-PC1 statistic's null is asymptotically
identical to plain univariate GSADF's does NOT hold at practical panel
widths N -- the true null quantile *grows* with N. That is why
radf_common_cv() (not radf_mc_cv()) must be used for critical values, and
why this script reproduces the SAME directional finding rather than
hiding it, per this project's instruction to keep the honest caveat in
the Python port too.

Reference numbers in the comments below (not asserted bit-for-bit --
cv.py's own module docstring already establishes numpy's Generator vs R's
RNG can't be matched bit-for-bit for any stochastic simulation in this
project) were produced by:

  Rscript -e '
    Sys.setenv(NOT_CRAN = "true")
    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    cv4  <- radf_common_cv(n = 100, N = 4,  nrep = 300, seed = 42)
    cv20 <- radf_common_cv(n = 100, N = 20, nrep = 300, seed = 42)
    cat("N=4  gsadf_cv:", cv4$gsadf_cv, "\n")
    cat("N=20 gsadf_cv:", cv20$gsadf_cv, "\n")
  '
  (run from exuber-project/ root)

  N=4  gsadf_cv: 2.184956 2.393796 3.348111   (90/95/99%)
  N=20 gsadf_cv: 3.281803 3.478511 3.968891   (90/95/99%)

  -- confirms the N-dependence: N=20's 95% (3.48) is well above N=4's
  (2.39), the same direction as docs/multivariate.md's own published
  table (N=6: 2.66, N=20: 3.45-3.48, N=50: 4.33-4.37, N=100: 5.20).

Two checks, run standalone from the exuber-project/ root:
  uv run --project pyexuber python docs/replication/multivariate/radf_common_validation.py
  1. PCA formula-exact check (pure numpy, no _core/radf() needed): the
     module's internal _pca() helper against an independent
     eigendecomposition of the sample covariance matrix.
  2. N-dependence sanity check (needs the compiled _core extension --
     see pyexuber/CLAUDE.md for why it can't be built on every machine):
     radf_common_cv() at two panel widths, checking the 95% gsadf_cv
     grows with N, matching the R finding above.

Not imported by pyexuber's own pytest suite: pyexuber is developed and
CI-tested in its own repo (kvasilopoulos/pyexuber), which doesn't check
out this umbrella repo's docs/ -- so pyexuber/tests/test_multivariate.py
re-implements these same two checks directly instead of importing this
file across that repo boundary. Keep the two in sync by hand if either
changes.
"""

import numpy as np

from exuber.radf_common import _pca, radf_common_cv


def check_pca_formula_exact(seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    n, nc, r = 60, 5, 2
    x = np.cumsum(rng.normal(size=(n, nc)), axis=0)

    loadings, scores, explained = _pca(x, r)

    # Independent recomputation via eigendecomposition of the sample
    # covariance matrix -- a different numerical path than SVD.
    xc = x - x.mean(axis=0)
    cov = (xc.T @ xc) / (n - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Eigenvector sign is arbitrary -- align each column's sign before
    # comparing, then require agreement to near machine precision.
    for j in range(r):
        sign = np.sign(np.dot(loadings[:, j], eigvecs[:, j])) or 1.0
        diff = float(np.max(np.abs(loadings[:, j] - sign * eigvecs[:, j])))
        assert diff < 1e-10, f"loadings column {j}: diff={diff:.2e}"

    expected_ratio = eigvals[:r] / np.sum(eigvals)
    assert np.allclose(explained, expected_ratio, atol=1e-10)
    assert np.allclose(scores, xc @ loadings, atol=1e-10)

    print(f"PCA formula-exact check OK (max loadings diff < 1e-10, r={r})")


def check_cv_grows_with_n(seed: int = 1) -> None:
    cv_small = radf_common_cv(n=60, N=4, nrep=150, seed=seed)
    cv_large = radf_common_cv(n=60, N=20, nrep=150, seed=seed)

    small_95 = cv_small.gsadf_cv[1]
    large_95 = cv_large.gsadf_cv[1]
    print(f"N=4  95% gsadf_cv: {small_95:.3f}")
    print(f"N=20 95% gsadf_cv: {large_95:.3f}")
    assert large_95 > small_95, (
        "expected N=20's null quantile to exceed N=4's (see R reference numbers "
        "in this file's module docstring) -- Theorem 4.3's N-independence claim "
        "does not hold at practical panel widths"
    )
    print("N-dependence check OK: N=20's 95% gsadf_cv > N=4's")


if __name__ == "__main__":
    check_pca_formula_exact()
    check_cv_grows_with_n()
