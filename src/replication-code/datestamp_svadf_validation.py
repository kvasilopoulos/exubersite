"""Python counterpart of datestamp_svadf_validation.R -- cross-checks
pyexuber's port of `datestamp(option="svadf")` (Sarkar & Wells 2026's
SV-ADF asymmetric-threshold dating, R/svadf.R + R/radf-methods.R's
datestamp_svadf()) in exuber.datestamp.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/datestamp_svadf_validation.py

This is a preprint (not peer-reviewed) -- flagged explicitly, same bar as
the R implementation and the port's own runtime warning.

Checks 1-2 (threshold formulas, structural collapse-never-before
-origination) need no exuber._core. Checks 3-4 (dating accuracy on a
synthetic bubble, false-alarm rate under H0) call radf() to build the
badf sequence datestamp() dates against, so they only run where the
compiled extension is available (CI) -- see radf_kp_validation.py's own
module docstring for the same caveat.
"""

import numpy as np

from exuber.datestamp import datestamp, svadf_threshold
from exuber.radf import psy_ds, radf


def check_threshold_formulas() -> None:
    t = np.array([100.0])
    orig = svadf_threshold(t, "origination")
    coll = svadf_threshold(t, "collapse")
    print(f"svadf_threshold(100, 'origination') = {orig[0]:.6f}  expect log(100)/10 = {np.log(100) / 10:.6f}")
    print(f"svadf_threshold(100, 'collapse')     = {coll[0]:.6f}  expect log(100)/2  = {np.log(100) / 2:.6f}")
    np.testing.assert_allclose(orig, np.log(100) / 10)
    np.testing.assert_allclose(coll, np.log(100) / 2)


def check_collapse_never_before_origination() -> None:
    """Needs exuber._core (radf() call). Structural check: since the
    collapse search only starts after the origination row, End can never
    precede Start when both are found -- across 20 reps."""
    rng = np.random.default_rng(3)
    ok = True
    for _ in range(20):
        yy = np.cumsum(rng.normal(size=150))
        rr = radf(yy, lag=0)
        out = datestamp(rr, option="svadf", min_duration=psy_ds(150))
        if out:
            ep = next(iter(out.values()))[0]
            if not ep.ongoing and ep.end is not None and ep.end <= ep.start:
                ok = False
    print(f"collapse always after origination when detected: {ok}")
    assert ok


def check_dating_accuracy_on_synthetic_bubble() -> None:
    """Needs exuber._core. Same synthetic bubble+collapse construction as
    the R validation script (large base bubble DGP)."""
    orig_err, coll_err = [], []
    detected = 0
    nrep = 20
    for i in range(nrep):
        r = np.random.default_rng(100 + i)
        n1 = 60
        yy1 = 100 + np.cumsum(r.normal(size=n1))
        n2 = 40
        bubble = yy1[-1] * 1.04 ** np.arange(1, n2 + 1) + np.cumsum(r.normal(scale=1, size=n2))
        n3 = 40
        coll = bubble[-1] - np.cumsum(np.abs(r.normal(loc=3, scale=1, size=n3)))
        yy = np.concatenate([yy1, bubble, coll])
        rr = radf(yy, lag=0)
        out = datestamp(rr, option="svadf", min_duration=psy_ds(len(yy)))
        if out:
            detected += 1
            ep = next(iter(out.values()))[0]
            orig_err.append(abs(ep.start - n1))
            if not ep.ongoing and ep.end is not None:
                coll_err.append(abs(ep.end - (n1 + n2)))
    print(f"detection rate: {detected / nrep:.3f}")
    if orig_err:
        print(f"mean |origination error|: {np.mean(orig_err):.3f}")
    if coll_err:
        print(f"mean |collapse error|: {np.mean(coll_err):.3f}")


def check_false_alarm_rate_under_h0() -> None:
    """Needs exuber._core. Pure random walk, 60 reps."""
    fa = 0
    nrep = 60
    for i in range(nrep):
        r = np.random.default_rng(1000 + i)
        yy = np.cumsum(r.normal(size=150))
        rr = radf(yy, lag=0)
        out = datestamp(rr, option="svadf", min_duration=psy_ds(150))
        if out:
            fa += 1
    print(f"false origination-alarm rate: {fa / nrep:.3f}")


if __name__ == "__main__":
    check_threshold_formulas()
    check_collapse_never_before_origination()
    check_dating_accuracy_on_synthetic_bubble()
    check_false_alarm_rate_under_h0()
    print("done")
