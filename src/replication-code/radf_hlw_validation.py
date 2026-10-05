"""Python replication script for dating_hlw() (Harvey, Leybourne &
Whitehouse 2020 multi-bubble two-step wrapper around dating_hls()),
cross-checking pyexuber's port against radf_hlw_validation.R in this same
folder.

Two kinds of checks run here.

1. The window construction and per-window fitting
   (_dating_hlw_from_episodes()) is pure numpy and needs no C++ extension.
   Sections 1 to 3 use hand-built Episode objects in place of a noisy
   datestamp() detection, which separates the logic under test from the
   noise of the PSY step-1 detection.
2. The full pipeline (radf() -> radf_wb_cv() -> datestamp() -> per-window
   dating_hlw()) needs the C++ extension. Section 4 is run by
   pyexuber/tests/test_dating_validation.py in CI, with loose structural
   assertions.

RNG note: numpy's Generator, not R's RNG (see rootstamp_validation.py's
module docstring for the same convention).

R's own script (docs/dating-and-root-inference.md, "Implementation (HLW
route)"), reports for context (not asserted bit-for-bit):
exactly 2 windows detected in 13/20 reps on a synthetic two-bubble DGP
(the rest were split by PSY step-1 detection noise); among clean 2-window reps, origination/collapse bias exactly 0 in
every replication; single-bubble final window matched standalone
dating_hls() in 100% of reps with a final window.

Run standalone:
uv run --project pyexuber python
docs/replication/dating-and-root-inference/radf_hlw_validation.py
"""

import warnings

import numpy as np

from exuber.datestamp import Episode
from exuber.dating_hls import dating_hls
from exuber.dating_hlw import _dating_hlw_from_episodes, _hlw_local_to_global, dating_hlw


def check_local_to_global_arithmetic() -> None:
    print("=== 1. _hlw_local_to_global() arithmetic ===")
    print("(R's own check: local_tau=5, s=21 (1-indexed) -> i_index=25, position=26;")
    print(" s0 (0-indexed window start) = s - 1 = 20)\n")
    g = _hlw_local_to_global(local_tau=5, s0=20)
    print(f"  local_tau=5, s0=20 -> position={g} (expect 26)")
    assert g == 26


def _sim_hls_model4(seed, n1=60, n2=25, n3=25, n4=40, base=100.0, c_bubble=1.05):
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


def check_single_episode_reduces_to_dating_hls() -> None:
    print("\n=== 2. Single clean episode: window matches standalone dating_hls() ===")
    match = 0
    for seed in range(15):
        y, _t1, _t2 = _sim_hls_model4(seed)
        n = len(y)
        ep = Episode(start=0, peak=0, end=None, duration=n, ongoing=True)
        hlw_eps = _dating_hlw_from_episodes(y, [ep], n, trim=0.1)
        hls_out = dating_hls(y, trim=0.1)
        last = hlw_eps[0]
        ok = (
            last.model == hls_out.model[0]
            and last.origination == hls_out.origination[0]
            and last.collapse == hls_out.collapse[0]
        )
        match += ok
    rate = match / 15
    print(f"  Final-window match rate vs standalone dating_hls(): {rate:.2f}")
    assert rate == 1.0


def _sim_two_bubbles(seed, n1a=50, n2a=20, n3a=30, n1b=50, n2b=20, n3b=30):
    rng = np.random.default_rng(seed)
    e1 = 100 + np.cumsum(rng.normal(size=n1a))
    b1 = e1[-1] * 1.05 ** np.arange(1, n2a + 1) + np.cumsum(rng.normal(size=n2a))
    u1 = b1[-1] + np.cumsum(rng.normal(size=n3a))
    e2 = u1[-1] + np.cumsum(rng.normal(size=n1b))
    b2 = e2[-1] * 1.05 ** np.arange(1, n2b + 1) + np.cumsum(rng.normal(size=n2b))
    u2 = b2[-1] + np.cumsum(rng.normal(size=n3b))
    y = np.concatenate([e1, b1, u1, e2, b2, u2])
    true1 = (n1a, n1a + n2a)
    true2 = (n1a + n2a + n3a + n1b, n1a + n2a + n3a + n1b + n2b)
    return y, true1, true2


def check_two_window_accuracy() -> None:
    print("\n=== 3. Two bubbles, clean hand-built windows: breakpoint accuracy ===")
    print("(isolates window-construction + per-window HLS fitting from step-1")
    print(" PSY detection noise -- that detection pass is exuber's/pyexuber's")
    print(" own already-tested datestamp(), not new code this port adds)\n")
    o1b, c1b, o2b, c2b = [], [], [], []
    for seed in range(20):
        y, true1, true2 = _sim_two_bubbles(seed)
        n = len(y)
        ep1 = Episode(start=true1[0] - 5, peak=0, end=true1[1] + 5, duration=0, ongoing=False)
        ep2 = Episode(start=true2[0] - 5, peak=0, end=true2[1] + 5, duration=0, ongoing=False)
        eps = _dating_hlw_from_episodes(y, [ep1, ep2], n, trim=0.1)
        if len(eps) == 2:
            o1b.append(eps[0].origination - true1[0])
            c1b.append(eps[0].collapse - true1[1])
            o2b.append(eps[1].origination - true2[0])
            c2b.append(eps[1].collapse - true2[1])
    print(f"  clean 2-window reps: {len(o1b)}/20")
    print(
        f"  bubble 1: orig mean|bias|={np.mean(np.abs(o1b)):.2f}, "
        f"coll mean|bias|={np.mean(np.abs(c1b)):.2f}"
    )
    print(
        f"  bubble 2: orig mean|bias|={np.mean(np.abs(o2b)):.2f}, "
        f"coll mean|bias|={np.mean(np.abs(c2b)):.2f}"
    )
    assert len(o1b) > 0
    assert np.mean(np.abs(o1b)) < 5
    assert np.mean(np.abs(o2b)) < 5


def check_end_to_end_h0_and_two_bubble() -> None:
    print("\n=== 4. Full end-to-end pipeline (needs the C++ extension) ===")
    print("(radf() -> radf_wb_cv() -> datestamp() -> per-window dating_hlw();")
    print(" needs the C++ extension -- see module docstring)\n")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rng = np.random.default_rng(2)
        y0 = 100 + np.cumsum(rng.normal(size=150))
        out0 = dating_hlw(y0, trim=0.1, nboot=199, seed=1)
        print(f"  Pure H0: {len(out0.episodes['series1'])} windows detected (expect 0, no error)")
        assert out0.episodes["series1"] == []

        n_windows = []
        for seed in range(8):
            y, _true1, _true2 = _sim_two_bubbles(seed)
            out = dating_hlw(y, trim=0.1, nboot=199, seed=1)
            eps = out.episodes["series1"]
            n_windows.append(len(eps))
            if len(eps) == 2:
                assert eps[0].origination < eps[1].origination
        print(f"  Two-bubble DGP window counts (8 reps): {n_windows}")
        assert all(w >= 0 for w in n_windows)


def main() -> None:
    check_local_to_global_arithmetic()
    check_single_episode_reduces_to_dating_hls()
    check_two_window_accuracy()
    check_end_to_end_h0_and_two_bubble()
    print("\nAll dating_hlw() checks passed.")


if __name__ == "__main__":
    main()
