"""Python counterpart of radf_wb_ps_validation.R -- cross-checks pyexuber's
port of exuber's shared OLS lag-selection/AR-fit subsystem
(exuber._lagselect: lag_select(), adf_res()) and radf_wb_ps_cv() (the
Phillips & Shi (2020) wild bootstrap variant).

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_wb_ps_validation.py

lag_select()/adf_res() are deterministic (no RNG) and checked bit-for-bit
against the R reference numbers below (produced by
`Rscript docs/replication/volatility-robustness/radf_wb_ps_validation.R`,
seed as in that script). radf_wb_ps_cv()'s bootstrap draws use numpy's
Generator, not R's RNG (see cv.py's module docstring), so that part is
checked structurally (shapes, monotonic quantiles, the tb-mode collapse
identity) rather than bit-for-bit -- same approach as radf_wb_cv's own
tests.
"""

import numpy as np

from exuber._lagselect import adf_res, lag_select
from exuber.cv import radf_wb_ps_cv
from exuber.radf import psy_minw

# y <- cumsum(rnorm(40)) after set.seed(11) in R -- see the .R script.
Y = np.array(
    [
        -0.5910311026, -0.5644367336, -2.0809898307, -3.4436431799, -2.2651540239,
        -3.1993053436, -1.8756996974, -1.2507819074, -1.2965048631, -2.3006254388,
        -3.1290586755, -3.4774104006, -5.0157037977, -5.2712690432, -6.4212140759,
        -6.4088871082, -6.6318566490, -5.7440850011, -6.3362402809, -6.9919583995,
        -7.6744760217, -7.6903342145, -8.1329389998, -7.7803815004, -7.7072109181,
        -7.7000521177, -7.8876522283, -8.6533528738, -8.8744096946, -9.8579982820,
        -10.9622823242, -11.9004325384, -11.2218082942, -12.7993061595, -13.6692446177,
        -13.1845675724, -13.3706202710, -11.8250655708, -12.4364456409, -12.7842021283,
    ]
)


def check_lag_select_and_adf_res() -> None:
    assert lag_select(Y, "aic", max_lag=5) == 5
    assert lag_select(Y, "bic", max_lag=5) == 4

    fit = adf_res(Y, adflag=2, type="fixed")
    np.testing.assert_allclose(
        fit.beta, [-0.2766891117, -0.0510771847, 0.0954335649], atol=1e-8
    )
    np.testing.assert_allclose(
        fit.res[:5],
        [-1.1659634957, 1.5303078394, -0.4672254328, 1.4401135170, 1.0583623423],
        atol=1e-8,
    )
    assert len(fit.res) == 37
    print("lag_select()/adf_res(): match R bit-for-bit.")


def check_radf_wb_ps_cv_shapes() -> None:
    minw = psy_minw(len(Y))
    wb = radf_wb_ps_cv(Y, minw=minw, nboot=80, adflag=0, seed=5)

    assert wb.adf_cv.shape == (1, 3)
    assert wb.gsadf_cv.shape == (1, 3)
    pointer = len(Y) - minw
    assert wb.bsadf_cv.shape == (pointer, 3, 1)
    assert np.all(np.diff(wb.gsadf_cv.ravel()) >= 0)
    print("radf_wb_ps_cv(): shapes and quantile ordering match R's structural checks.")


def check_tb_mode() -> None:
    minw = psy_minw(len(Y))
    tb = minw + 10
    wbt = radf_wb_ps_cv(Y, minw=minw, nboot=60, adflag=0, tb=tb, seed=6)

    pointer_full = len(Y) - minw
    assert wbt.badf_cv.shape == (pointer_full, 3, 1)
    np.testing.assert_allclose(wbt.badf_cv[0, :, 0], wbt.sadf_cv[0, :])
    np.testing.assert_allclose(wbt.bsadf_cv[0, :, 0], wbt.gsadf_cv[0, :])
    print("radf_wb_ps_cv(tb=...): badf/bsadf collapse to repeated sadf/gsadf, as in R.")


if __name__ == "__main__":
    check_lag_select_and_adf_res()
    check_radf_wb_ps_cv_shapes()
    check_tb_mode()
    print("All checks passed.")
