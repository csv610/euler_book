"""
Euler's Infinite Series

This module implements Euler's work on infinite series including:
- Harmonic series and Euler's constant γ
- Basel problem: ζ(2) = π²/6
- Alternating series
- Zeta function at even integers
- Leibniz formula for π

References:
- Euler, L. (1735). "De summis serierum reciprocarum"
- Euler, L. (1748). "Introductio in analysin infinitorum"
"""

from math import pi, log, factorial
from typing import List, Tuple


def harmonic_series(n: int) -> float:
    """Sum of first n terms of harmonic series H_n = 1 + 1/2 + ... + 1/n."""
    return sum(1.0 / k for k in range(1, n + 1))


def harmonic_series_partial(n: int) -> List[float]:
    """Return list of partial sums H_1, H_2, ..., H_n."""
    partials = []
    s = 0.0
    for k in range(1, n + 1):
        s += 1.0 / k
        partials.append(s)
    return partials


def basel_problem(n_terms: int = None) -> float:
    """
    Basel problem: ζ(2) = ∑_{n=1}^∞ 1/n² = π²/6
    If n_terms is None, return the exact value π²/6.
    Otherwise return partial sum.
    """
    if n_terms is None:
        return pi ** 2 / 6
    return sum(1.0 / (k * k) for k in range(1, n_terms + 1))


def basel_partial_sum(n: int) -> List[float]:
    """Return list of partial sums for Basel problem."""
    partials = []
    s = 0.0
    for k in range(1, n + 1):
        s += 1.0 / (k * k)
        partials.append(s)
    return partials


def euler_constant(n: int) -> float:
    """
    Approximate Euler-Mascheroni constant γ using H_n - ln(n).
    γ = lim_{n→∞} (H_n - ln n) ≈ 0.5772156649...
    """
    return harmonic_series(n) - log(n)


def euler_constant_partial(n: int) -> List[float]:
    """Return sequence of γ approximations H_k - ln(k) for k=1..n."""
    return [harmonic_series(k) - log(k) for k in range(1, n + 1)]


def euler_mascheroni_gamma() -> float:
    """Return the known value of Euler-Mascheroni constant."""
    return 0.57721566490153286060651209008240243104215933593992


def alternating_harmonic(n: int) -> float:
    """Alternating harmonic series: ∑_{k=1}^n (-1)^{k+1}/k → ln(2)."""
    return sum((1.0 if k % 2 == 1 else -1.0) / k for k in range(1, n + 1))


def leibniz_pi(n: int) -> float:
    """
    Leibniz formula for π: π/4 = ∑_{k=0}^{∞} (-1)^k/(2k+1)
    Returns approximation using n terms.
    """
    return 4 * sum((1.0 if k % 2 == 0 else -1.0) / (2 * k + 1) for k in range(n))


def zeta_function(s: float, n_terms: int = 100000) -> float:
    """
    Riemann zeta function ζ(s) = ∑_{n=1}^∞ 1/n^s for Re(s) > 1.
    For even integers, Euler found exact values: ζ(2n) = (-1)^{n+1} B_{2n} (2π)^{2n} / (2(2n)!)
    """
    if s <= 1:
        raise ValueError("Zeta function series converges only for s > 1")
    return sum(1.0 / (k ** s) for k in range(1, n_terms + 1))


def zeta_even(n: int) -> float:
    """
    Euler's formula for ζ(2n):
    ζ(2n) = (-1)^{n+1} * B_{2n} * (2π)^{2n} / (2 * (2n)!)
    
    Using known Bernoulli numbers for small n:
    B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, B_8 = -1/30, B_10 = 5/66, ...
    """
    # Bernoulli numbers B_2, B_4, B_6, B_8, B_10, B_12, B_14, B_16
    bernoulli = {
        1: 1/6,
        2: -1/30,
        3: 1/42,
        4: -1/30,
        5: 5/66,
        6: -691/2730,
        7: 7/6,
        8: -3617/510,
    }
    
    if n not in bernoulli:
        raise ValueError(f"Bernoulli number B_{2*n} not available for n={n}")
    
    B = bernoulli[n]
    sign = 1 if (n + 1) % 2 == 0 else -1  # (-1)^{n+1}
    return sign * B * (2 * pi) ** (2 * n) / (2 * factorial(2 * n))


def euler_acceleration_alternating(series_fn, n: int) -> float:
    """
    Euler's transformation for accelerating alternating series.
    If S = ∑(-1)^k a_k, then S = ∑_{k=0}^∞ (-1)^k Δ^k a_0 / 2^{k+1}
    where Δ is the forward difference operator.
    """
    # Compute first n+1 terms
    a = [series_fn(k) for k in range(n + 1)]
    
    # Compute forward differences
    result = 0.0
    for k in range(n + 1):
        # Δ^k a_0
        diff = a[0]
        for _ in range(k):
            diff = a[_ + 1] - a[_]
        result += ((-1) ** k) * diff / (2 ** (k + 1))
    return result


def series_for_e(n: int) -> float:
    """Euler's series for e: e = ∑_{k=0}^∞ 1/k!"""
    return sum(1.0 / factorial(k) for k in range(n))


def euler_constant_series(n: int) -> float:
    """
    Faster converging series for γ using Gregory coefficients:
    γ = ∑_{k=1}^∞ G_k / k where G_k are Gregory coefficients.
    """
    # Gregory coefficients: G_1 = 1/2, G_2 = -1/12, G_3 = 1/24, G_4 = -19/720, ...
    gregory = [0, 1/2, -1/12, 1/24, -19/720, 3/160, -863/60480]
    return sum(gregory[k] / k for k in range(1, min(n + 1, len(gregory))))


def sum_of_reciprocal_squares(n: int) -> float:
    """∑_{k=1}^n 1/k² - partial sum for Basel problem."""
    return sum(1.0 / (k * k) for k in range(1, n + 1))


def sum_of_reciprocal_powers(n: int, p: int) -> float:
    """∑_{k=1}^n 1/k^p"""
    return sum(1.0 / (k ** p) for k in range(1, n + 1))


# Example usage and demonstrations
if __name__ == "__main__":
    print("=== Euler's Infinite Series Demonstrations ===\n")
    
    # Harmonic series
    print("Harmonic Series H_n:")
    for n in [1, 2, 3, 4, 5, 10, 100]:
        print(f"  H_{n} = {harmonic_series(n):.6f}")
    
    print("\nEuler-Mascheroni constant γ = lim(H_n - ln n):")
    for n in [10, 100, 1000, 10000]:
        print(f"  n={n}: γ ≈ {euler_constant(n):.10f}")
    print(f"  Known value: {euler_mascheroni_gamma():.10f}")
    
    print("\nBasel Problem ζ(2) = π²/6:")
    for n in [1, 2, 5, 10, 100, 1000]:
        approx = basel_problem(n)
        exact = basel_problem()
        print(f"  n={n}: {approx:.10f} (error: {abs(exact - approx):.2e})")
    print(f"  Exact: π²/6 = {exact:.10f}")
    
    print("\nZeta at even integers (Euler's formula):")
    for n in range(1, 6):
        exact = zeta_even(n)
        approx = zeta_function(2 * n, 100000)
        print(f"  ζ({2*n}) = {exact:.10f} (series: {approx:.10f})")
    
    print("\nLeibniz formula for π:")
    for n in [1, 10, 100, 1000, 10000]:
        print(f"  n={n}: π ≈ {leibniz_pi(n):.10f}")
    
    print("\nAlternating harmonic series → ln(2):")
    for n in [10, 100, 1000, 10000]:
        print(f"  n={n}: {alternating_harmonic(n):.10f} (ln(2)={log(2):.10f})")
    
    print("\nEuler's series for e:")
    for n in [5, 10, 15, 20]:
        print(f"  n={n}: e ≈ {series_for_e(n):.10f}")