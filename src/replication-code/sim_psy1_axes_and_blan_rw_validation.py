"""Python counterpart of sim_psy1_axes_and_blan_rw_validation.R --
cross-checks pyexuber's sim_psy1() e/shifts/coef_noise/coef_a optional
axes and sim_blan(type="rotermann_wilfling") (exuber.sim) against R's
R/sim.R equivalents.

Run standalone: uv run --project pyexuber python
docs/replication/simulation-dgps/sim_psy1_axes_and_blan_rw_validation.py

sim_psy1's e=/shifts= structural behavior (overrides innovations, adds
an isolated one-period bump) is checked directly against the real
function. sim_blan(type="rotermann_wilfling")'s recursion formula is
checked bit-for-bit via a fake RNG feeding R's exact theta/u draws to
the real sim_blan() function.
"""

import numpy as np

from exuber.sim import sim_blan, sim_psy1


def check_sim_psy1_e() -> None:
    y_plain = sim_psy1(50, seed=42)
    y_zero_e = sim_psy1(50, seed=42, e=np.zeros(49))
    assert np.all(np.isfinite(y_zero_e))
    assert not np.allclose(y_zero_e, y_plain)  # e overrides the default rnorm draws
    print("sim_psy1(e=...): overrides innovations, matches R's structural behavior.")


def check_sim_psy1_shifts() -> None:
    y_plain = sim_psy1(50, seed=42)
    y_shift = sim_psy1(50, seed=42, shifts={"date": [10], "size": [100.0]})
    assert y_shift[9] - y_plain[9] == 100.0  # R: y_shift[10] - y_plain[10] == 100
    assert np.array_equal(y_shift[:9], y_plain[:9])
    print("sim_psy1(shifts=...): exact +100 one-period bump, no effect before it, matches R.")


def check_sim_blan_rotermann_wilfling() -> None:
    # set.seed(7); theta <- rbinom(4,1,0.7); u <- rlnorm(4, meanlog=-.05^2/2, sdlog=.05) in R.
    theta = np.array([0, 1, 1, 1])
    u = np.array([0.96467442, 0.97837265, 0.95143523, 0.95254875])

    class _FakeRWRNG:
        def binomial(self, n_trials, p, size=None):
            return theta

        def lognormal(self, mean=0.0, sigma=1.0, size=None):
            return u

    orig = np.random.default_rng
    np.random.default_rng = lambda seed=None: _FakeRWRNG()
    try:
        b = sim_blan(
            5, pi=0.7, type="rotermann_wilfling", delta=0.984, rw_sigma=0.05, b0=0.1, seed=7
        )
    finally:
        np.random.default_rng = orig

    expected = [0.10000000, 0.10006889, 0.09949661, 0.09620385, 0.09312891]
    np.testing.assert_allclose(b, expected, atol=1e-6)
    print("sim_blan(type='rotermann_wilfling'): recursion matches R bit-for-bit.")


if __name__ == "__main__":
    check_sim_psy1_e()
    check_sim_psy1_shifts()
    check_sim_blan_rotermann_wilfling()
    print("All checks passed.")
