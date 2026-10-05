"""Python counterpart of radf_sb_cv_aic_bic_validation.R -- cross-checks
pyexuber's port of radf_sb_cv()/radf_sb_distr() (Pavlidis et al. 2016
sieve bootstrap, R/radf_sb.R) and its Pedersen & Schuette (2020) AIC/BIC
automatic lag selection, which reuses the shared lag_select() from
exuber._lagselect (the same subsystem radf_wb_ps_validation.py checks).

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_sb_cv_aic_bic_validation.py

lag_select() is deterministic and checked bit-for-bit against R (two
independent series below, reusing radf_wb_ps_validation.R's own
seed-11 series for the nonzero-lag case so the reference numbers aren't
duplicated). The bootstrap-DGP part of radf_sb_cv() itself uses numpy's
Generator, not R's RNG, so it's checked structurally, as in
radf_wb_ps_validation.py.
"""

import numpy as np

from exuber._lagselect import lag_select
from exuber.cv import radf_sb_cv
from exuber.radf import psy_minw

# y <- cumsum(rnorm(40)) after set.seed(11) -- same series as
# radf_wb_ps_validation.R/.py; lag_select(y, 'bic', max_lag=5) == 4 there.
Y_NONZERO_LAG = np.array(
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

# y <- cumsum(rnorm(120)) after set.seed(901) (the s=1 draw of the existing
# R script's "modal lag 0 on pure random-walk data" check, reproduced with:
# set.seed(1 + 900); y <- cumsum(rnorm(120))) --
# lag_select(y, 'bic', max_lag=6) == 0 in R.
Y_ZERO_LAG = np.array(
    [
        0.8280005254, 0.9037403631, 1.0297638781, 3.2953730540, 2.4571820742,
        1.9881313218, 2.0036552915, 2.9649966337, 2.6837130519, 2.7534546993,
        2.3854981840, 2.1791112564, 2.7034419784, 1.1850348963, 0.8810668612,
        1.2207721547, 0.0086998789, -0.0403166404, -0.3933561298, 0.7179381554,
        0.5767380450, 2.8443491868, 2.3306774450, 1.8950008302, 0.5124493609,
        0.7907321379, 1.4771141254, 2.1673821143, 3.3046109976, 2.2997063424,
        5.4637789994, 6.7368314769, 6.8170202655, 6.7320646259, 6.2367888184,
        6.5863425516, 6.5152545756, 7.1525906165, 7.3899634326, 8.2405404171,
        6.7936506516, 4.8922572610, 3.3091199930, 4.4928717647, 5.8216340213,
        5.1727613201, 6.4639183405, 7.6068791535, 6.8849483108, 5.3793020469,
        4.3203202966, 3.2886314359, 5.2465118704, 5.4659803669, 6.2022668447,
        6.4388892920, 6.5529224693, 7.7872883608, 8.1538456891, 8.1673367230,
        9.7791618427, 10.3803305345, 9.7183931202, 9.9031786321, 10.7545557098,
        10.6928446347, 9.5782195803, 9.9318977725, 9.3332014070, 9.3562971436,
        10.4522518578, 11.9305100522, 12.5635943914, 12.1295041898, 10.8993394960,
        10.9233483727, 10.2472321968, 11.9068619299, 11.1611669364, 10.7446409638,
        10.7294083363, 12.4878449174, 12.1559236311, 11.0333089876, 9.4305997938,
        9.8618063784, 9.9668011181, 9.1291369253, 10.8708990612, 11.2428864427,
        10.1246813256, 10.1516123071, 10.0709268466, 10.8126420199, 10.4041225606,
        8.7313835197, 8.6688233086, 9.4735560265, 8.4051664325, 5.9722189285,
        7.9492864599, 6.4858189000, 6.7738346470, 7.1397413025, 6.7550253996,
        8.0640710638, 9.3552276060, 9.1812074033, 9.7607815192, 10.6397064963,
        9.8385103895, 8.1932508953, 8.3612089825, 8.2865638915, 8.8064160287,
        10.5658044873, 10.1533012384, 8.5183272011, 8.4386235415, 7.3939110514,
    ]
)


def check_lag_select_bit_for_bit() -> None:
    assert lag_select(Y_NONZERO_LAG, "bic", max_lag=5) == 4
    assert lag_select(Y_NONZERO_LAG, "aic", max_lag=5) == 5
    assert lag_select(Y_ZERO_LAG, "bic", max_lag=6) == 0
    print("lag_select(): matches R bit-for-bit (nonzero-lag bic=4/aic=5; zero-lag bic=0).")


def check_radf_sb_cv_fixed_vs_default() -> None:
    """type='fixed' (the default) is unaffected by the new type/max_lag
    arguments -- mirrors the R script's check 1."""
    minw = psy_minw(len(Y_NONZERO_LAG))
    a = radf_sb_cv(Y_NONZERO_LAG, minw=minw, lag=2, nboot=40, seed=22)
    b = radf_sb_cv(Y_NONZERO_LAG, minw=minw, lag=2, type="fixed", nboot=40, seed=22)
    np.testing.assert_array_equal(a.gsadf_panel_cv, b.gsadf_panel_cv)
    np.testing.assert_array_equal(a.bsadf_panel_cv, b.bsadf_panel_cv)
    print("radf_sb_cv(): default type='fixed' identical to explicit type='fixed'.")


def check_radf_sb_cv_shapes_lag_gt_0() -> None:
    """For lag > 0, bsadf_panel_cv has `nr - minw - lag` rows."""
    minw = psy_minw(len(Y_NONZERO_LAG))
    lag = 2
    sb = radf_sb_cv(Y_NONZERO_LAG, minw=minw, lag=lag, nboot=30, seed=1)
    expected_pointer = len(Y_NONZERO_LAG) - minw - lag
    assert sb.bsadf_panel_cv.shape == (expected_pointer, 3), (
        f"got {sb.bsadf_panel_cv.shape}, expected ({expected_pointer}, 3)"
    )
    assert sb.lag == lag
    print(f"radf_sb_cv(lag={lag}): bsadf_panel_cv shape is the full (nr-minw-lag, 3).")


if __name__ == "__main__":
    check_lag_select_bit_for_bit()
    check_radf_sb_cv_fixed_vs_default()
    check_radf_sb_cv_shapes_lag_gt_0()
    print("All checks passed.")
