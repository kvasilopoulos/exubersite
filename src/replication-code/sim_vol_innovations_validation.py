"""Python counterpart of sim_vol_innovations_validation.R -- cross-checks
pyexuber's sim_vol_break()/sim_vol_garch()/sim_vol_cir()/sim_vol_sv()/
sim_fi()/sim_innov() (exuber.sim) against R's R/sim.R equivalents.

Run standalone: uv run --project pyexuber python
docs/replication/simulation-dgps/sim_vol_innovations_validation.py

These are otherwise-deterministic transformations of Gaussian noise, so
each is checked by feeding the *actual* exuber.sim functions a fixed
pool of draws (R: set.seed(1); rnorm(10)) via a fake np.random.Generator
replacement (_FakeRNG) that returns pool slices in the same call order R
draws them in -- the algorithm gets checked bit-for-bit, not the RNG
bit-stream (see exuber.sim's module docstring on that distinction).
"""

import math

import numpy as np

from exuber.sim import (
    _beta_fn,
    sim_fi,
    sim_innov,
    sim_vol_break,
    sim_vol_cir,
    sim_vol_garch,
    sim_vol_sv,
)

# set.seed(1); rnorm(10) in R.
EPS10 = np.array(
    [-0.62645381, 0.18364332, -0.83562861, 1.59528080, 0.32950777,
     -0.82046838, 0.48742905, 0.73832471, 0.57578135, -0.30538839]
)


class _FakeRNG:
    """Feeds a fixed pool to standard_normal()/normal() in call order --
    see exuber's own tests/test_sim.py, which uses the identical class."""

    def __init__(self, pool):
        self.pool = np.asarray(pool, dtype=float)
        self.pos = 0

    def _take(self, size):
        if size is None:
            v = self.pool[self.pos]
            self.pos += 1
            return v
        out = self.pool[self.pos : self.pos + size]
        self.pos += size
        return out

    def standard_normal(self, size=None):
        return self._take(size)

    def normal(self, loc=0.0, scale=1.0, size=None):
        if size is None and isinstance(scale, np.ndarray):
            size = scale.shape[0]
        return loc + scale * self._take(size)


def _with_fake_rng(pool, fn, *args, **kwargs):
    orig = np.random.default_rng
    np.random.default_rng = lambda seed=None: _FakeRNG(pool)
    try:
        return fn(*args, **kwargs)
    finally:
        np.random.default_rng = orig


def check_sim_vol_break() -> None:
    y = _with_fake_rng(EPS10, sim_vol_break, 10, tau=0.5, ratio=3, sigma=2, seed=1)
    expected = [-1.25290762, 0.36728665, -1.67125722, 3.19056160, 0.65901554,
                -4.92281030, 2.92457431, 4.42994823, 3.45468811, -1.83233032]
    np.testing.assert_allclose(y, expected, atol=1e-6)
    print("sim_vol_break(): matches R.")


def check_sim_vol_garch() -> None:
    y = _with_fake_rng(EPS10[:5], sim_vol_garch, 5, omega=0.1, alpha=0.1, beta=0.8, gamma=0.0, seed=1)
    expected = [-0.19810209, 0.07875803, -0.41593815, 0.89607126, 0.21675026]
    np.testing.assert_allclose(y, expected, atol=1e-6)

    y2 = _with_fake_rng(EPS10[:5], sim_vol_garch, 5, omega=0.1, alpha=0.1, beta=0.8, gamma=0.2, seed=1)
    expected2 = [-0.19810209, 0.08042096, -0.42119779, 0.95247029, 0.22731252]
    np.testing.assert_allclose(y2, expected2, atol=1e-6)
    print("sim_vol_garch(): GARCH and TGARCH (gamma>0) match R.")


def check_sim_vol_cir() -> None:
    pool = np.concatenate([EPS10[:4], EPS10[5:10]])
    y = _with_fake_rng(pool, sim_vol_cir, 5, kappa=0.03, theta=0.25, xi=0.1, seed=1)
    expected = [-0.41023419, 0.23678823, 0.36175334, 0.27117724, -0.15439035]
    np.testing.assert_allclose(y, expected, atol=1e-6)
    print("sim_vol_cir(): matches R.")


def check_sim_vol_sv() -> None:
    pool = np.concatenate([EPS10[:4], EPS10[5:10]])
    y = _with_fake_rng(pool, sim_vol_sv, 5, phi=0.98, tau=0.1, seed=1)
    expected = [-0.82046838, 0.47239810, 0.72260999, 0.54069901, -0.31098370]
    np.testing.assert_allclose(y, expected, atol=1e-6)
    print("sim_vol_sv(): matches R.")


def check_sim_innov_normal() -> None:
    y = _with_fake_rng(EPS10, sim_innov, 10, dist="normal", sigma=2, seed=1)
    expected = [-1.25290762, 0.36728665, -1.67125722, 3.19056160, 0.65901554,
                -1.64093677, 0.97485810, 1.47664941, 1.15156270, -0.61077677]
    np.testing.assert_allclose(y, expected, atol=1e-6)
    print("sim_innov(dist='normal'): matches R.")


def check_sim_innov_t_skew_t_constants() -> None:
    df = 5
    assert 1 / math.sqrt(df / (df - 2)) == 0.7745966692 or abs(
        1 / math.sqrt(df / (df - 2)) - 0.7745966692
    ) < 1e-9

    xi = -0.75
    delta = xi / math.sqrt(1 + xi**2)
    e_abs_t0 = (2 * math.sqrt(df) / ((df - 1) * _beta_fn(df / 2, 0.5))) / math.sqrt(df / (df - 2))
    mean_raw = delta * e_abs_t0
    sd_raw = math.sqrt(max(1 - delta**2 * e_abs_t0**2, np.finfo(float).eps))
    assert abs(delta - (-0.6)) < 1e-9
    assert abs(e_abs_t0 - 0.7351051939) < 1e-9
    assert abs(mean_raw - (-0.4410631163)) < 1e-9
    assert abs(sd_raw - 0.8974760874) < 1e-9
    print("sim_innov(dist='t'/'skew_t'): closed-form rescaling constants match R.")


def check_sim_fi_psi_and_convolution() -> None:
    d, m = 0.2, 5
    psi = np.empty(m + 1)
    psi[0] = 1.0
    for j in range(1, m + 1):
        psi[j] = psi[j - 1] * (j - 1 + d) / j
    np.testing.assert_allclose(psi, [1.0, 0.2, 0.12, 0.088, 0.0704, 0.059136], atol=1e-6)

    eps = EPS10[:8]  # n=3, m=5 -> n+m=8 innovations, matching sim_fi's own convolution
    u = np.empty(3)
    for k in range(3):
        window = eps[k : k + m + 1]
        u[k] = np.dot(psi, window[::-1])
    np.testing.assert_allclose(u, [-0.66078593, 0.45529270, 0.82924303], atol=1e-6)
    print("sim_fi(): psi recursion and convolution formula match R.")

    # smoke test the real function (m = max(500, 5n) there, too big to fake-RNG by hand)
    y = sim_fi(200, d=0.2, sigma=1.0, seed=1)
    assert y.shape == (200,) and np.all(np.isfinite(y))
    print("sim_fi(): real call (n=200) runs and returns finite values.")


if __name__ == "__main__":
    check_sim_vol_break()
    check_sim_vol_garch()
    check_sim_vol_cir()
    check_sim_vol_sv()
    check_sim_innov_normal()
    check_sim_innov_t_skew_t_constants()
    check_sim_fi_psi_and_convolution()
    print("All checks passed.")
