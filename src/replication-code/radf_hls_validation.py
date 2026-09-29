"""Python replication script for dating_hls() (Harvey, Leybourne & Sollis
2017 SSR+BIC dating), cross-checking pyexuber's port against
radf_hls_validation.R in this same folder.

RNG note: numpy's Generator, not R's RNG -- independent seeds, same
qualitative checks (see rootstamp_validation.py's module docstring for
the same convention elsewhere in this port). dating_hls() is pure numpy
(no C++ extension needed), so -- unlike radf_recovery_validation.py --
this script WAS run directly on the dev machine that wrote the port; all
numbers below are real local runs, not placeholders.

R's own script (docs/dating-and-root-inference.md, "Implementation (HLS
route)"), run 2026-08-10, reports for context (not asserted bit-for-bit,
different RNG): Model 4 DGP -- BIC selects Model 3/4 in 100% of 30 reps,
split ~80/20; Model 2 DGP -- 100% Model 2, origination bias exactly 0;
Model 1 DGP -- 100% Model 1, bias exactly 0; pure H0 -- 0% Model 4,
split ~30/57/13 across Models 1/2/3.

Run standalone:
uv run --project pyexuber python
docs/replication/dating-and-root-inference/radf_hls_validation.py
"""

import math

import numpy as np

from exuber._hls_common import _hls_model4, _hls_prefix_sums, _hls_segment_ssr
from exuber.dating_hls import dating_hls


def _ols_ssr(xseg: np.ndarray, zseg: np.ndarray) -> float:
    a = np.vstack([xseg, np.ones_like(xseg)]).T
    coef, *_ = np.linalg.lstsq(a, zseg, rcond=None)
    return float(np.sum((zseg - a @ coef) ** 2))


def check_segment_ssr_formula_exact() -> None:
    print("=== 1. Formula-exact: _hls_segment_ssr() vs brute-force OLS ===")
    rng = np.random.default_rng(1)
    y = np.cumsum(rng.normal(size=40))
    ps = _hls_prefix_sums(y)
    n1 = len(y) - 1
    x_all, z_all = y[:n1], np.diff(y)
    for lo, hi in [(0, 10), (10, 25), (4, 39)]:
        manual = _hls_segment_ssr(ps, lo, hi, True)
        brute = _ols_ssr(x_all[lo:hi], z_all[lo:hi])
        print(f"  segment ({lo},{hi}]: manual={manual:.6f} brute={brute:.6f}")
        assert abs(manual - brute) < 1e-8


def check_model4_formula_exact() -> None:
    print("\n=== 2. Full Model 4 (3-breakpoint) joint grid search vs brute force ===")
    rng = np.random.default_rng(5)
    n = 24
    y = np.cumsum(rng.normal(size=n))
    ps = _hls_prefix_sums(y)
    tau1, tau2, tau3, ssr = _hls_model4(y, ps, trim=0.1)

    n1 = n - 1
    x, z = y[:n1], np.diff(y)
    k_min = max(2, math.ceil(0.1 * n1))
    best_ssr, best = math.inf, None
    for t1 in range(k_min, n1 - 3 * k_min + 1):
        for t2 in range(t1 + k_min, n1 - 2 * k_min + 1):
            if y[t2] <= y[t1]:
                continue
            for t3 in range(t2 + k_min, n1 - k_min + 1):
                if y[t2] <= y[t3]:
                    continue
                s = (
                    np.sum(z[:t1] ** 2)
                    + _ols_ssr(x[t1:t2], z[t1:t2])
                    + _ols_ssr(x[t2:t3], z[t2:t3])
                    + np.sum(z[t3:n1] ** 2)
                )
                if s < best_ssr:
                    best_ssr, best = s, (t1, t2, t3)
    print(f"  vectorized: {(tau1, tau2, tau3)}, ssr={ssr}")
    print(f"  brute:      {best}, ssr={best_ssr}")
    assert (tau1, tau2, tau3) == best
    assert abs(ssr - best_ssr) < 1e-6


def check_performance() -> None:
    print("\n=== 3. Performance at realistic sample sizes ===")
    import time

    for n in (100, 200, 400):
        rng = np.random.default_rng(1)
        y = np.cumsum(rng.normal(size=n))
        t0 = time.time()
        dating_hls(y, trim=0.05)
        print(f"  n={n}: {time.time() - t0:.2f} sec")


def _sim_model4(seed, n1=60, n2=25, n3=25, n4=40, base=100.0, c_bubble=1.05):
    rng = np.random.default_rng(seed)
    unit1 = base + np.cumsum(rng.normal(size=n1))
    bubble = unit1[-1] * c_bubble ** np.arange(1, n2 + 1) + np.cumsum(rng.normal(size=n2))
    target = bubble[-1] * 0.5
    collapse = np.empty(n3)
    collapse[0] = bubble[-1] + rng.normal()
    for k in range(1, n3):
        collapse[k] = target + 0.85 * (collapse[k - 1] - target) + rng.normal()
    recovery = collapse[-1] + np.cumsum(rng.normal(size=n4))
    return np.concatenate([unit1, bubble, collapse, recovery]), n1, n1 + n2


def _sim_model2(seed, n1=60, n2=30, n3=60, base=100.0, c_bubble=1.05):
    rng = np.random.default_rng(seed)
    unit1 = base + np.cumsum(rng.normal(size=n1))
    bubble = unit1[-1] * c_bubble ** np.arange(1, n2 + 1) + np.cumsum(rng.normal(size=n2))
    unit2 = bubble[-1] + np.cumsum(rng.normal(size=n3))
    return np.concatenate([unit1, bubble, unit2]), n1


def _sim_model1(seed, n1=80, n2=60, base=100.0, c_bubble=1.05):
    rng = np.random.default_rng(seed)
    unit1 = base + np.cumsum(rng.normal(size=n1))
    bubble = unit1[-1] * c_bubble ** np.arange(1, n2 + 1) + np.cumsum(rng.normal(size=n2))
    return np.concatenate([unit1, bubble]), n1


def check_monte_carlo_model_selection() -> None:
    print("\n=== 4. Monte Carlo: model-selection accuracy and breakpoint bias by DGP ===")
    print("(bubble/collapse regimes on a large positive base (100), matching R's")
    print(" own fix for a weak-signal DGP that biased selection toward Model 1)\n")

    res4 = []
    for s in range(30):
        y, t1, t2 = _sim_model4(s)
        out = dating_hls(y, trim=0.05)
        coll = out.collapse[0] - t2 if not math.isnan(out.collapse[0]) else None
        res4.append((out.model[0], out.origination[0] - t1, coll))
    models4 = [m for m, _, _ in res4]
    freq4 = {m: models4.count(m) / 30 for m in (1, 2, 3, 4)}
    print(f"  Model 4 DGP -- selection freq: {freq4}")
    orig_bias4 = np.mean([abs(b) for _, b, _ in res4])
    coll_biases = [c for m, _, c in res4 if m in (3, 4) and c is not None]
    coll_bias4 = np.mean([abs(c) for c in coll_biases]) if coll_biases else float("nan")
    print(f"  origination mean|bias|={orig_bias4:.2f}; collapse mean|bias| (M3/4)={coll_bias4:.2f}")
    assert freq4[3] + freq4[4] > 0.8  # BIC should favor a distinct-collapse model

    res2 = []
    for s in range(30):
        y, t1 = _sim_model2(s)
        out = dating_hls(y, trim=0.05)
        res2.append((out.model[0], out.origination[0] - t1))
    models2 = [m for m, _ in res2]
    freq2 = {m: models2.count(m) / 30 for m in (1, 2, 3, 4)}
    print(f"  Model 2 DGP -- selection freq: {freq2}; origination mean|bias|="
          f"{np.mean([abs(b) for _, b in res2]):.2f}")
    assert freq2[2] > 0.8

    res1 = []
    for s in range(30):
        y, t1 = _sim_model1(s)
        out = dating_hls(y, trim=0.05)
        res1.append((out.model[0], out.origination[0] - t1))
    models1 = [m for m, _ in res1]
    freq1 = {m: models1.count(m) / 30 for m in (1, 2, 3, 4)}
    print(f"  Model 1 DGP -- selection freq: {freq1}; origination mean|bias|="
          f"{np.mean([abs(b) for _, b in res1]):.2f}")
    assert freq1[1] > 0.8

    models0 = []
    for s in range(30):
        rng = np.random.default_rng(s)
        y0 = 100 + np.cumsum(rng.normal(size=150))
        out = dating_hls(y0, trim=0.05)
        models0.append(out.model[0])
    freq0 = {m: models0.count(m) / 30 for m in (1, 2, 3, 4)}
    print(f"  Pure H0 (no bubble) -- selection freq: {freq0}")
    assert freq0[4] < 0.3  # no runaway complex-model selection under pure noise


def main() -> None:
    check_segment_ssr_formula_exact()
    check_model4_formula_exact()
    check_performance()
    check_monte_carlo_model_selection()
    print("\nAll dating_hls() checks passed.")


if __name__ == "__main__":
    main()
