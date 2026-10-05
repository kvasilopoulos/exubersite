"""Python replication script for dating_pdc() (Pang, Du & Chong 2021 /
Kurozumi & Skrobotov 2023 sequential sample-splitting bubble dating, both
the base OLS route and the WLS volatility correction), cross-checking
pyexuber's port against radf_pdc_validation.R and the two
radf_pdc_wls_*_mae.R scripts in this same folder.

RNG note: numpy's Generator, not R's RNG -- independent seeds, same
qualitative checks (see rootstamp_validation.py's module docstring for the
same convention). R-derived reference numbers come from:

    cd exuber-project && Rscript docs/replication/dating-and-root-inference/radf_pdc_validation.R
    cd exuber-project && Rscript docs/replication/dating-and-root-inference/radf_pdc_wls_homoskedastic_mae.R
    cd exuber-project && Rscript docs/replication/dating-and-root-inference/radf_pdc_wls_heteroskedastic_mae.R

R 4.6.1:
  - formula check: closed-form break_idx == brute-force break_idx, abs
    diff in rss = 1.28e-13
  - 3-regime low-noise limit: origination |err|=1, collapse |err|=1 (of 400/600)
  - 4-regime low-noise limit: all three breaks |err|=1 (of 300/450/650)
  - homoskedastic MAE: OLS origination=5.42, WLS=5.53 (WLS costs ~nothing)
  - heteroskedastic (volatility burst) MAE: OLS origination=13.05, WLS=2.33
    (WLS materially better)

Run standalone: uv run --project pyexuber python
docs/replication/dating-and-root-inference/radf_pdc_validation.py
"""

import numpy as np

from exuber.dating_pdc import _pdc_find_break, dating_pdc


def check_formula_exact() -> None:
    print("=== 1. Formula check: _pdc_find_break() vs brute-force RSS scan ===")
    rng = np.random.default_rng(9001)
    y = np.cumsum(rng.normal(size=120))
    trim = 0.05
    break_idx, rss = _pdc_find_break(y, trim)

    n1 = len(y) - 1
    ylag, ycur = y[:n1], y[1 : n1 + 1]
    k_min = max(2, int(np.ceil(trim * n1)))
    k_max = n1 - k_min

    def rss_ols(x, yv):
        beta = np.sum(x * yv) / np.sum(x * x)
        return np.sum((yv - beta * x) ** 2)

    ks = list(range(k_min, k_max + 1))
    rss_brute = [rss_ols(ylag[:k], ycur[:k]) + rss_ols(ylag[k:n1], ycur[k:n1]) for k in ks]
    brute_break = ks[int(np.argmin(rss_brute))]

    print(f"  closed-form break_idx: {break_idx}  brute-force break_idx: {brute_break}")
    print(f"  abs diff in rss: {abs(rss - min(rss_brute))}")
    assert break_idx == brute_break
    assert abs(rss - min(rss_brute)) < 1e-8


def _three_regime(rng, n1_len, n2_len, n3_len, scales=(1, 0.1, 0.5), rho3=0.5, c=1.08):
    s1, s2, s3 = scales
    regime1 = np.cumsum(rng.normal(size=n1_len, scale=s1))
    regime2 = regime1[-1] * c ** np.arange(1, n2_len + 1) + np.cumsum(rng.normal(size=n2_len, scale=s2))
    peak = regime2[-1]
    regime3 = np.zeros(n3_len)
    regime3[0] = rho3 * peak + rng.normal(scale=s3)
    for t in range(1, n3_len):
        regime3[t] = rho3 * regime3[t - 1] + rng.normal(scale=s3)
    return regime1, regime2, regime3


def check_3regime_consistency() -> None:
    print("\n=== 2. 3-regime consistency in the low-noise/long-series/strong-effect limit ===")
    rng = np.random.default_rng(4001)
    n1_len, n2_len, n3_len = 400, 200, 250
    regime1, regime2, regime3 = _three_regime(rng, n1_len, n2_len, n3_len)
    y = np.concatenate([regime1, regime2, regime3])

    out = dating_pdc(y, regimes=3, trim=0.05)
    true_o, true_c = n1_len, n1_len + n2_len
    print(f"  origination: true={true_o}  est={out.origination[0]}  |err|={abs(out.origination[0]-true_o)}")
    print(f"  collapse:    true={true_c}  est={out.collapse[0]}  |err|={abs(out.collapse[0]-true_c)}")
    assert abs(out.origination[0] - true_o) <= 2
    assert abs(out.collapse[0] - true_c) <= 2


def check_4regime_consistency() -> None:
    print("\n=== 3. 4-regime (KS extension) consistency in the same low-noise limit ===")
    rng = np.random.default_rng(4)
    n1_len, n2_len, n3_len, n4_len = 400, 200, 250, 250
    regime1, regime2, regime3 = _three_regime(rng, n1_len, n2_len, n3_len, scales=(0.3, 0.05, 0.5))
    regime4 = regime3[-1] + np.cumsum(rng.normal(size=n4_len, scale=0.3))
    y = np.concatenate([regime1, regime2, regime3, regime4])

    out = dating_pdc(y, regimes=4, trim=0.05)
    true_o, true_c, true_r = n1_len, n1_len + n2_len, n1_len + n2_len + n3_len
    print(f"  origination: true={true_o}  est={out.origination[0]}  |err|={abs(out.origination[0]-true_o)}")
    print(f"  collapse:    true={true_c}  est={out.collapse[0]}  |err|={abs(out.collapse[0]-true_c)}")
    print(f"  recovery:    true={true_r}  est={out.recovery[0]}  |err|={abs(out.recovery[0]-true_r)}")
    assert abs(out.origination[0] - true_o) <= 3
    assert abs(out.collapse[0] - true_c) <= 3
    assert abs(out.recovery[0] - true_r) <= 5


def check_wls_homoskedastic_mae() -> None:
    print("\n=== 4. WLS vs OLS MAE under homoskedasticity (should be close) ===")
    n1_len, n2_len, n3_len = 150, 80, 100

    def run(seed):
        rng = np.random.default_rng(seed)
        regime1, regime2, regime3 = _three_regime(
            rng, n1_len, n2_len, n3_len, scales=(0.5, 0.15, 0.5), c=1.07
        )
        y = np.concatenate([regime1, regime2, regime3])
        true_o, true_c = n1_len, n1_len + n2_len
        out_ols = dating_pdc(y, regimes=3, trim=0.05, type="ols")
        out_wls = dating_pdc(y, regimes=3, trim=0.05, type="wls")
        return (
            abs(out_ols.origination[0] - true_o), abs(out_ols.collapse[0] - true_c),
            abs(out_wls.origination[0] - true_o), abs(out_wls.collapse[0] - true_c),
        )

    res = np.array([run(s) for s in range(40)])
    mae = res.mean(axis=0)
    print(f"  Origination MAE: OLS={mae[0]:.2f}  WLS={mae[2]:.2f}")
    print(f"  Collapse MAE:     OLS={mae[1]:.2f}  WLS={mae[3]:.2f}")
    # WLS shouldn't be dramatically worse when there's no volatility signal
    assert mae[2] < 2 * mae[0] + 1


def check_wls_heteroskedastic_mae() -> None:
    print("\n=== 5. WLS vs OLS MAE under a volatility burst (WLS should win) ===")
    n1_len, n2_len, n3_len = 150, 80, 100
    burst_len = round(0.2 * n1_len)

    def run(seed):
        rng = np.random.default_rng(seed)
        e1 = np.concatenate(
            [rng.normal(size=burst_len, scale=4), rng.normal(size=n1_len - burst_len, scale=0.3)]
        )
        regime1 = np.cumsum(e1)
        regime2 = regime1[-1] * 1.07 ** np.arange(1, n2_len + 1) + np.cumsum(
            rng.normal(size=n2_len, scale=0.15)
        )
        peak = regime2[-1]
        rho3 = 0.5
        regime3 = np.zeros(n3_len)
        regime3[0] = rho3 * peak + rng.normal(scale=0.5)
        for t in range(1, n3_len):
            regime3[t] = rho3 * regime3[t - 1] + rng.normal(scale=0.5)
        y = np.concatenate([regime1, regime2, regime3])
        true_o = n1_len
        out_ols = dating_pdc(y, regimes=3, trim=0.05, type="ols")
        out_wls = dating_pdc(y, regimes=3, trim=0.05, type="wls")
        return abs(out_ols.origination[0] - true_o), abs(out_wls.origination[0] - true_o)

    res = np.array([run(s) for s in range(40)])
    mae_ols, mae_wls = res.mean(axis=0)
    print(f"  Origination MAE: OLS={mae_ols:.2f}  WLS={mae_wls:.2f}  (WLS better if smaller)")
    assert mae_wls < mae_ols


def main() -> None:
    check_formula_exact()
    check_3regime_consistency()
    check_4regime_consistency()
    check_wls_homoskedastic_mae()
    check_wls_heteroskedastic_mae()
    print("\nAll dating_pdc() checks passed.")


if __name__ == "__main__":
    main()
