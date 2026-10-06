"""Timing of exuber.radf() in Python: median elapsed time per call on a random
walk of length n (median of five block medians), with minw=30 and lag=1, for n = 100, 200, ..., 1000.

This is the same grid and the same settings as speed_r.R, which times the R
package. speed_plot.R draws both series. Run it from the project root with a
pyexuber build that includes the compiled core (``pip install pyexuber``)::

    python docs/replication/speed/speed_py.py
"""

from __future__ import annotations

import csv
import statistics
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np

import exuber

MINW = 30
LAG = 1
SAMPLE_SIZE = range(100, 1001, 100)
TIMES = 60  # calls per block
BLOCKS = 5  # the reported time is the median of the block medians


def median_ms(n: int) -> float:
    """Median time of ``radf()`` in milliseconds on a random walk of length n."""
    rw = np.cumsum(np.random.default_rng(123).standard_normal(n))
    exuber.radf(rw, minw=MINW, lag=LAG)  # warm-up call
    blocks = []
    for _ in range(BLOCKS):
        runs = []
        for _ in range(TIMES):
            t0 = time.perf_counter()
            exuber.radf(rw, minw=MINW, lag=LAG)
            runs.append((time.perf_counter() - t0) * 1e3)
        blocks.append(statistics.median(runs))
    return statistics.median(blocks)


def main() -> None:
    rows = []
    for n in SAMPLE_SIZE:
        t_ms = median_ms(n)
        print(f"n = {n:4d}   pyexuber {version('pyexuber')}   {t_ms:9.3f} ms")
        rows.append((f"pyexuber {version('pyexuber')} (Python)", n, t_ms))

    out = Path("docs/replication/speed/speed-py.csv")
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["software", "n", "time_ms"])
        w.writerows(rows)
    print("written:", out)


if __name__ == "__main__":
    main()
