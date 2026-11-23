"""Verification tests for the computational-discovery material."""

import math
import sys
from fractions import Fraction

sys.path.insert(0, "src")

from euler.discovery import (  # noqa: E402
    bernoulli_numbers,
    brachistochrone_parameter,
    brachistochrone_time,
    continued_fraction,
    continued_fraction_convergents,
    elliptic_e,
    elliptic_k,
    euler_sum_harmonic_over_square,
    faulhaber_sum,
    integrate_euler_top,
    integrate_two_body,
    rigid_body_energy,
    two_body_angular_momentum,
    two_body_energy,
)


def test_bernoulli_and_faulhaber_values():
    assert bernoulli_numbers(6) == [
        Fraction(1), Fraction(-1, 2), Fraction(1, 6), Fraction(0),
        Fraction(-1, 30), Fraction(0), Fraction(1, 42)
    ]
    assert [faulhaber_sum(10, k) for k in range(1, 5)] == [55, 385, 3025, 25333]


def test_continued_fraction_convergents():
    coefficients = continued_fraction(Fraction(355, 113))
    assert coefficients == [3, 7, 16]
    assert continued_fraction_convergents(coefficients)[-1] == Fraction(355, 113)


def test_elliptic_integrals_at_zero():
    assert abs(elliptic_k(0.0, 200) - math.pi / 2) < 1e-10
    assert abs(elliptic_e(0.0, 200) - math.pi / 2) < 1e-10


def test_brachistochrone_is_finite():
    theta = brachistochrone_parameter(1.0, 1.0)
    assert 0 < theta < 2 * math.pi
    assert 0 < brachistochrone_time(1.0, 1.0) < 1.0


def test_euler_sum_converges_to_known_identity():
    # The slowly convergent direct sum is intentionally tested loosely.
    value = euler_sum_harmonic_over_square(100_000)
    assert abs(value - 2.4041138063191885) < 2e-4


def test_euler_top_energy_is_numerically_conserved():
    moments = (1.0, 2.0, 3.0)
    trajectory = integrate_euler_top((1.0, 0.2, 0.1), moments,
                                     step=0.002, steps=1000)
    energies = [rigid_body_energy(state, moments) for state in trajectory]
    assert max(energies) - min(energies) < 1e-10


def test_circular_two_body_motion_preserves_invariants():
    orbit = integrate_two_body((1.0, 0.0, 0.0, 1.0), step=0.002, steps=3141)
    initial_energy = two_body_energy(orbit[0])
    final_energy = two_body_energy(orbit[-1])
    initial_momentum = two_body_angular_momentum(orbit[0])
    final_momentum = two_body_angular_momentum(orbit[-1])
    assert abs(initial_energy - final_energy) < 1e-9
    assert abs(initial_momentum - final_momentum) < 1e-9
    # The chosen step count is close to, but not exactly, one period; check the
    # circular radius rather than requiring a phase error that the grid cannot provide.
    assert abs(math.hypot(orbit[-1][0], orbit[-1][1]) - 1.0) < 1e-9
