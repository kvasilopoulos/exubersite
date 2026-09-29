"""Python counterpart of radf_kp_validation.R -- cross-checks pyexuber's
port of radf_kp() (kernel-purge test, Harvey, Leybourne, Taylor & Zu 2024,
R/radf_kp.R) in exuber.volatility.

Run standalone: uv run --project pyexuber python
docs/replication/volatility-robustness/radf_kp_validation.py

Checks kernel_spot_vol()/kernel_purge()/radf_kp() bit-for-bit against R,
fed the SAME deterministic input series used by radf_tt_validation.py and
sign_based_finite_T_crosscheck.py (no RNG involved at the formula level):

    options(exuber.parallel = FALSE, exuber.show_progress = FALSE)
    devtools::load_all("exuber", quiet = TRUE)
    set.seed(7); y <- round(cumsum(rnorm(40)), 8)
    exuber:::kernel_spot_vol(y, kernel = "gaussian")
    exuber:::kernel_purge(y, kernel = "gaussian")
    radf_kp(y, minw = 10)

radf_kp() itself needs exuber._core (it calls radf() internally), so the
last check only runs where the compiled extension is available (CI).
"""

import numpy as np

from exuber._kernel_vol import kernel_spot_vol
from exuber.radf_kp import kernel_purge, radf_kp

MINW = 10

Y_VEC = np.array(
    [
        2.28724716, 1.09047548, 0.39618297, -0.01610998, -0.98678332, -1.93406327,
        -1.18592393, -1.30287915, -1.15022153, 1.03975658, 1.39674281, 4.11349459,
        6.39494652, 6.71896706, 8.61503413, 9.08271464, 8.18891391, 7.88158561,
        7.87676319, 8.86492734, 9.7046777, 10.41001953, 11.71598425, 10.32798804,
        11.6009049, 11.78509767, 12.53737757, 13.12912262, 12.14607002, 11.87000607,
        10.99915505, 11.7178656, 11.82851848, 11.75005171, 11.32956125, 10.76743537,
        11.76494882, 10.65981876, 10.51753093, 10.83252583,
    ]
)

R_KSV_SIGMA2 = np.array(
    [
        0.4959746477, 0.5312260867, 0.5716433342, 0.6169370612, 0.6874127672,
        0.8412310930, 1.1577525057, 1.6723606201, 2.3107950220, 2.8881768207,
        3.1903140074, 3.0969314316, 2.6568967882, 2.0450515634, 1.4477027141,
        0.9864934574, 0.7115504863, 0.6200852851, 0.6744228299, 0.8193921797,
        0.9912734704, 1.1188219752, 1.1404876822, 1.0422514405, 0.8736254828,
        0.7103030842, 0.5977433306, 0.5317693774, 0.4823425959, 0.4254286845,
        0.3600372811, 0.3091871957, 0.3044556076, 0.3582846686, 0.4454091397,
        0.5159624129, 0.5322132895, 0.4917189956, 0.4177513067,
    ]
)
R_KSV_H = 0.0533043701

R_PURGED = np.array(
    [
        -1.7873025622, -2.7815277011, -3.3347316743, -4.5671929336, -5.7535320421,
        -4.8328468384, -4.9568680837, -4.8312994075, -3.3989648256, -3.1994688263,
        -1.7867645834, -0.5703660507, -0.3766526498, 0.9732799329, 1.3858598183,
        0.3777725907, -0.0501407948, -0.0571310824, 1.2267423867, 2.1748043745,
        2.8717529025, 4.0401985161, 2.8275287339, 4.0251706963, 4.2266895665,
        5.1809673847, 5.9841391137, 4.6283810407, 4.2390555386, 2.9405290986,
        4.1438219706, 4.3633073113, 4.1957073655, 3.4324194005, 2.6141599989,
        3.8927306251, 2.4765274358, 2.2724582562, 2.8167023625,
    ]
)

R_KP_ADF = -0.9285832329
R_KP_SADF = 0.6113330426
R_KP_GSADF = 0.8288627249


def check_kernel_spot_vol_matches_r() -> None:
    sigma2, h = kernel_spot_vol(Y_VEC, kernel="gaussian")
    np.testing.assert_allclose(sigma2, R_KSV_SIGMA2, atol=1e-6)
    np.testing.assert_allclose(h, R_KSV_H, atol=1e-6)
    print("kernel_spot_vol(): matches R's exuber:::kernel_spot_vol() to 1e-6.")


def check_kernel_purge_matches_r() -> None:
    purged = kernel_purge(Y_VEC, kernel="gaussian")
    np.testing.assert_allclose(purged, R_PURGED, atol=1e-6)
    print("kernel_purge(): matches R's exuber:::kernel_purge() to 1e-6.")


def check_radf_kp_matches_r() -> None:
    res = radf_kp(Y_VEC, minw=MINW)
    np.testing.assert_allclose(res.adf[0], R_KP_ADF, atol=1e-6)
    np.testing.assert_allclose(res.sadf[0], R_KP_SADF, atol=1e-6)
    np.testing.assert_allclose(res.gsadf[0], R_KP_GSADF, atol=1e-6)
    print("radf_kp(): matches R's radf_kp() to 1e-6.")


if __name__ == "__main__":
    check_kernel_spot_vol_matches_r()
    check_kernel_purge_matches_r()
    check_radf_kp_matches_r()
    print("All checks passed.")
