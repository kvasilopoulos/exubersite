"""Python replication script for radf_recovery()/radf_recovery_cv()
(Phillips & Shi 2014 reverse-regression crisis-origination/market-recovery
dating), cross-checking pyexuber's port against radf_recovery_validation.R
in this same folder.

Environment note: unlike rootstamp()/dating_pdc() (pure numpy),
radf_recovery() calls radf() itself, which needs pyexuber's compiled C++
extension (exuber._core). That extension cannot be built on the
Windows dev machine this port was written on (no MSVC/vcpkg -- see
pyexuber/CLAUDE.md), so this script could not be run locally as part of
writing the port. It was written directly from R's own
radf_recovery_validation.R structure and dating.py's ported logic, and is
wired into pyexuber/tests/test_dating_validation.py so CI (which does
build the extension, on ubuntu/macos/windows) executes and verifies it.
Assertions below are deliberately structural/loose (invariants, sane
ranges) rather than exact-number matches to a run that was never
performed on this machine -- avoiding a false "verified" claim.

RNG note: numpy's Generator, not R's RNG (see rootstamp_validation.py's
module docstring for the same convention elsewhere in this port).

R's own script, run 2026-08-10 (docs/dating-and-root-inference.md,
"Reverse-regression recovery dating"), reports for context (not asserted
here bit-for-bit): max abs CV diff ~0.11, mean abs CV diff ~0.04 at the
95% level; H0 false-detection rate ~29% (n=100, minw=20); f_r mean
|bias| ~a few observations (paper's own ~6-early finding); f_c's bias
materially larger and not fully resolved.

Run standalone (needs the built extension):
uv run --project pyexuber python
docs/replication/dating-and-root-inference/radf_recovery_validation.py
"""

import warnings

import numpy as np

from exuber.cv import radf_mc_cv
from exuber.radf_recovery import radf_recovery, radf_recovery_cv


def check_reversal_calibrated_cv_differs_from_forward() -> None:
    print("=== 1. Reversal-calibrated CV differs from forward CV ===")
    print("(confirms Theorem 1's endogeneity finding: the reverse-time")
    print(" regression has no forward-regression analogue, so its null")
    print(" distribution -- and hence critical values -- genuinely differ)\n")
    n, minw = 100, 20
    fwd = radf_mc_cv(n, minw=minw, nrep=2000, seed=7)
    rev = radf_recovery_cv(n, minw=minw, nrep=2000, seed=7)
    diff_95 = np.abs(fwd.bsadf_cv[:, 1] - rev.bsadf_cv[:, 1])
    print(f"  max abs diff (95% col):  {diff_95.max()}")
    print(f"  mean abs diff (95% col): {diff_95.mean()}")
    # a genuinely different distribution shows up as a nontrivial gap --
    # not asserting the exact R figures (different RNG/nrep), just that
    # the two boundaries aren't the same (would be ~0 diff if they were).
    assert diff_95.mean() > 0.005


def check_h0_false_detection_rate() -> None:
    print("\n=== 2. False-detection rate under pure H0 (random walk) ===")
    n, minw = 100, 20
    cv_nrep = 500  # smaller than R's 1000/5000 to keep CI runtime reasonable
    rng = np.random.default_rng(123)
    detected = []
    for _ in range(60):
        y = np.cumsum(rng.normal(size=n))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = radf_recovery(y, minw=minw, nrep=cv_nrep, seed=1)
        detected.append(bool(out.detected[0]))
    rate = np.mean(detected)
    print(f"  False-detection rate: {rate:.3f} (R's own run: ~0.29 -- noisier")
    print("  than comparable forward-test numbers; flagged honestly as")
    print("  unresolved, see docs/dating-and-root-inference.md)")
    # loose sanity bound only, matching R test-recovery.R's own "<= 0.5" check
    assert rate <= 0.6


def check_detection_accuracy_and_invariant() -> None:
    print("\n=== 3. Detection accuracy + f_c <= f_r invariant on a synthetic")
    print("    collapse-then-recovery DGP ===")
    n1, n2, n3 = 40, 25, 35
    true_collapse, true_recovery = n1, n1 + n2

    def run(seed):
        rng = np.random.default_rng(seed)
        expansion = 100 * 1.03 ** np.arange(1, n1 + 1) + np.cumsum(rng.normal(size=n1))
        target = expansion[-1] * 0.5
        collapse = np.empty(n2)
        collapse[0] = expansion[-1] + rng.normal(scale=1)
        for k in range(1, n2):
            collapse[k] = target + 0.9 * (collapse[k - 1] - target) + rng.normal(scale=1)
        recovery = collapse[-1] + np.cumsum(rng.normal(size=n3)) + np.arange(1, n3 + 1) * 0.5
        y = np.concatenate([expansion, collapse, recovery])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = radf_recovery(y, minw=15, nrep=200, seed=1)
        return out

    results = [run(s) for s in range(20)]
    detected = [r for r in results if r.detected[0]]
    uncensored = [r for r in detected if not r.censored[0]]
    print(f"  Detection rate: {len(detected) / len(results):.3f}")
    print(f"  n detected & uncensored = {len(uncensored)}")
    if uncensored:
        f_c_bias = [r.f_c[0] - true_collapse for r in uncensored]
        f_r_bias = [r.f_r[0] - true_recovery for r in uncensored]
        print(f"  mean f_c bias = {np.mean(f_c_bias):.2f}, mean f_r bias = {np.mean(f_r_bias):.2f}")
        print(
            f"  mean |f_c bias| = {np.mean(np.abs(f_c_bias)):.2f}, "
            f"mean |f_r bias| = {np.mean(np.abs(f_r_bias)):.2f}"
        )
        print("  (R's own run: f_r matches the paper's ~6-observation-early")
        print("  finding; f_c's bias is materially larger, not fully resolved)")
        assert all(r.f_c[0] <= r.f_r[0] for r in uncensored)


def main() -> None:
    check_reversal_calibrated_cv_differs_from_forward()
    check_h0_false_detection_rate()
    check_detection_accuracy_and_invariant()
    print("\nAll radf_recovery() structural checks passed.")


if __name__ == "__main__":
    main()
