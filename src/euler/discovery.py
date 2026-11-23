"""Computational experiments for Euler-style mathematical discovery.

The functions in this module deliberately expose intermediate objects (terms,
convergents, trajectories, and errors) so that a student can observe a pattern
before being shown its proof.  Exact arithmetic is used for identities whenever
possible; numerical routines return ordinary floats and document their error
controls.
"""

from fractions import Fraction
from math import comb, cos, pi, sin, sqrt
from typing import Callable, Iterable, List, Sequence, Tuple


def bernoulli_numbers(n: int) -> List[Fraction]:
    """Return ``B_0, ..., B_n`` with the convention ``B_1 = -1/2``."""
    if n < 0:
        raise ValueError("n must be non-negative")
    result = [Fraction(1)]
    for m in range(1, n + 1):
        value = -sum(comb(m + 1, j) * result[j] for j in range(m))
        result.append(value / (m + 1))
    return result


def faulhaber_sum(n: int, power: int) -> int:
    """Return ``1**power + ... + n**power`` using Bernoulli numbers."""
    if n < 0 or power < 0:
        raise ValueError("n and power must be non-negative")
    b = bernoulli_numbers(power)
    value = sum(
        Fraction((-1) ** j * comb(power + 1, j), power + 1)
        * b[j] * n ** (power + 1 - j)
        for j in range(power + 1)
    )
    if value.denominator != 1:
        raise ArithmeticError("Faulhaber result was not integral")
    return value.numerator


def continued_fraction(value, max_terms: int = 32) -> List[int]:
    """Return the first continued-fraction coefficients of a number.

    ``Fraction`` inputs are exact.  Floating-point inputs are iterated only up
    to ``max_terms`` because their final coefficients are not reliable.
    """
    if max_terms <= 0:
        return []
    exact = isinstance(value, Fraction)
    x = value
    result = []
    for _ in range(max_terms):
        a = x.numerator // x.denominator if exact else int(x)
        result.append(a)
        remainder = x - a
        if remainder == 0:
            break
        x = 1 / remainder
    return result


def continued_fraction_convergents(coefficients: Sequence[int]) -> List[Fraction]:
    """Return all rational convergents associated with coefficients."""
    p_minus_2, p_minus_1 = 0, 1
    q_minus_2, q_minus_1 = 1, 0
    result = []
    for a in coefficients:
        p = a * p_minus_1 + p_minus_2
        q = a * q_minus_1 + q_minus_2
        result.append(Fraction(p, q))
        p_minus_2, p_minus_1 = p_minus_1, p
        q_minus_2, q_minus_1 = q_minus_1, q
    return result


def euler_sum_harmonic_over_square(n_terms: int = 100_000) -> float:
    """Approximate ``sum H_n/n**2``; the limit is ``2*zeta(3)``."""
    harmonic = 0.0
    total = 0.0
    for n in range(1, n_terms + 1):
        harmonic += 1.0 / n
        total += harmonic / (n * n)
    return total


def simpson_integral(f: Callable[[float], float], a: float, b: float,
                     intervals: int = 10_000) -> float:
    """Composite Simpson integration with an even number of intervals."""
    if intervals <= 0 or intervals % 2:
        raise ValueError("intervals must be a positive even number")
    h = (b - a) / intervals
    total = f(a) + f(b)
    total += 4 * sum(f(a + i * h) for i in range(1, intervals, 2))
    total += 2 * sum(f(a + i * h) for i in range(2, intervals, 2))
    return total * h / 3


def elliptic_k(modulus: float, intervals: int = 10_000) -> float:
    """Numerically evaluate the complete elliptic integral ``K(k)``."""
    if not 0 <= modulus < 1:
        raise ValueError("the modulus must satisfy 0 <= k < 1")
    return simpson_integral(
        lambda theta: 1.0 / sqrt(1.0 - modulus * modulus * sin(theta) ** 2),
        0.0, pi / 2, intervals)


def elliptic_e(modulus: float, intervals: int = 10_000) -> float:
    """Numerically evaluate the complete elliptic integral ``E(k)``."""
    if not 0 <= modulus <= 1:
        raise ValueError("the modulus must satisfy 0 <= k <= 1")
    return simpson_integral(
        lambda theta: sqrt(1.0 - modulus * modulus * sin(theta) ** 2),
        0.0, pi / 2, intervals)


def brachistochrone_parameter(x: float, y: float, tolerance: float = 1e-12) -> float:
    """Solve ``y/x = (1-cos(theta))/(theta-sin(theta))`` for the cycloid.

    Coordinates are measured downward from the starting point.  The returned
    parameter is in ``(0, 2*pi)`` for positive ``x`` and ``y``.
    """
    if x <= 0 or y <= 0:
        raise ValueError("x and y must be positive")
    ratio = y / x

    def f(theta: float) -> float:
        return (1 - cos(theta)) / (theta - sin(theta)) - ratio

    lo, hi = 1e-4, 2 * pi - 1e-10
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tolerance:
            break
    return (lo + hi) / 2


def brachistochrone_time(x: float, y: float, gravity: float = 9.80665) -> float:
    """Return the ideal cycloid travel time ``theta*sqrt(a/g)``."""
    theta = brachistochrone_parameter(x, y)
    a = y / (1 - cos(theta))
    return theta * sqrt(a / gravity)


def euler_top_rhs(state: Sequence[float], moments: Sequence[float]) -> Tuple[float, float, float]:
    """Right-hand side for angular-momentum coordinates of a free rigid body."""
    l1, l2, l3 = state
    i1, i2, i3 = moments
    return (
        (i2 - i3) * l2 * l3 / (i2 * i3),
        (i3 - i1) * l3 * l1 / (i3 * i1),
        (i1 - i2) * l1 * l2 / (i1 * i2),
    )


def integrate_euler_top(initial: Sequence[float], moments: Sequence[float],
                        step: float = 0.001, steps: int = 1000
                        ) -> List[Tuple[float, float, float]]:
    """Integrate the torque-free Euler top with classical RK4."""
    if len(initial) != 3 or len(moments) != 3:
        raise ValueError("initial and moments must each have length 3")
    state = [float(v) for v in initial]
    trajectory = [tuple(state)]
    for _ in range(steps):
        f = lambda s: euler_top_rhs(s, moments)
        k1 = f(state)
        s2 = [state[i] + step * k1[i] / 2 for i in range(3)]
        k2 = f(s2)
        s3 = [state[i] + step * k2[i] / 2 for i in range(3)]
        k3 = f(s3)
        s4 = [state[i] + step * k3[i] for i in range(3)]
        k4 = f(s4)
        state = [state[i] + step * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6
                 for i in range(3)]
        trajectory.append(tuple(state))
    return trajectory


def rigid_body_energy(angular_momentum: Sequence[float], moments: Sequence[float]) -> float:
    """Rotational kinetic energy for principal-axis coordinates."""
    return sum(l * l / (2 * i) for l, i in zip(angular_momentum, moments))


def integrate_two_body(initial: Sequence[float], step: float = 0.001,
                       steps: int = 1000, gravitational_parameter: float = 1.0
                       ) -> List[Tuple[float, float, float, float]]:
    """Integrate planar two-body motion ``(x,y,vx,vy)`` with RK4."""
    if len(initial) != 4:
        raise ValueError("initial state must have length 4")

    def rhs(s):
        x, y, vx, vy = s
        r3 = (x * x + y * y) ** 1.5
        return (vx, vy, -gravitational_parameter * x / r3,
                -gravitational_parameter * y / r3)

    state = [float(v) for v in initial]
    orbit = [tuple(state)]
    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs([state[i] + step * k1[i] / 2 for i in range(4)])
        k3 = rhs([state[i] + step * k2[i] / 2 for i in range(4)])
        k4 = rhs([state[i] + step * k3[i] for i in range(4)])
        state = [state[i] + step * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6
                 for i in range(4)]
        orbit.append(tuple(state))
    return orbit


def two_body_energy(state: Sequence[float], gravitational_parameter: float = 1.0) -> float:
    """Specific mechanical energy of a planar two-body state."""
    x, y, vx, vy = state
    return 0.5 * (vx * vx + vy * vy) - gravitational_parameter / sqrt(x * x + y * y)


def two_body_angular_momentum(state: Sequence[float]) -> float:
    """Specific angular momentum ``x*vy - y*vx``."""
    x, y, vx, vy = state
    return x * vy - y * vx
