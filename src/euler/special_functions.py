"""
Euler's Special Functions

This module implements Euler's contributions to special functions:
- Gamma function Γ(z)
- Beta function B(x, y)
- Euler's reflection formula
- Legendre's duplication formula
- Digamma function

References:
- Euler, L. (1729). "De progressionibus transcendentibus"
- Euler, L. (1768). "Institutiones calculi integralis"
"""

from math import gamma, pi, sin, cos, log, exp, sqrt, factorial
from typing import Tuple
import cmath


def gamma_function(z: float) -> float:
    """
    Gamma function Γ(z) for real z > 0.
    Γ(z) = ∫_0^∞ t^{z-1} e^{-t} dt
    """
    return gamma(z)


def gamma_function_lanczos(z: float, g: int = 7) -> float:
    """
    Lanczos approximation for Gamma function (works for complex z too).
    More accurate than built-in for large values.
    """
    # Lanczos coefficients for g=7
    p = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ]
    
    if z < 0.5:
        # Reflection formula
        return pi / (sin(pi * z) * gamma_function_lanczos(1 - z, g))
    
    z -= 1
    x = p[0]
    for i in range(1, g + 2):
        x += p[i] / (z + i)
    
    t = z + g + 0.5
    return sqrt(2 * pi) * t ** (z + 0.5) * exp(-t) * x


def beta_function(x: float, y: float) -> float:
    """
    Beta function B(x, y) = Γ(x)Γ(y)/Γ(x+y)
    B(x, y) = ∫_0^1 t^{x-1} (1-t)^{y-1} dt
    """
    return gamma(x) * gamma(y) / gamma(x + y)


def beta_function_integral(x: float, y: float, n: int = 10000) -> float:
    """
    Compute Beta function via numerical integration:
    B(x, y) = ∫_0^1 t^{x-1} (1-t)^{y-1} dt
    """
    h = 1.0 / n
    total = 0.0
    for i in range(n):
        t = (i + 0.5) * h
        if t == 0 or t == 1:
            continue
        total += t ** (x - 1) * (1 - t) ** (y - 1)
    return total * h


def digamma_function(z: float) -> float:
    """
    Digamma function ψ(z) = Γ'(z)/Γ(z).
    For large z: ψ(z) ≈ log(z) - 1/(2z) - 1/(12z²) + 1/(120z⁴) - ...
    """
    # Use recurrence to shift to z > 8
    result = 0.0
    while z < 8:
        result -= 1 / z
        z += 1
    
    # Asymptotic expansion
    inv_z = 1 / z
    inv_z2 = inv_z * inv_z
    result += log(z) - inv_z / 2 - inv_z2 / 12 + inv_z2 * inv_z2 / 120 - inv_z2 ** 3 / 252
    
    return result


def euler_reflection_formula(z: float) -> Tuple[float, float]:
    """
    Euler's reflection formula for Gamma function:
    Γ(z) Γ(1-z) = π / sin(πz)
    
    Returns (LHS, RHS) for verification.
    """
    lhs = gamma(z) * gamma(1 - z)
    rhs = pi / sin(pi * z)
    return (lhs, rhs)


def euler_reflection_demo(z_values: list = None) -> list:
    """Demonstrate Euler's reflection formula."""
    if z_values is None:
        z_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    
    results = []
    for z in z_values:
        lhs, rhs = euler_reflection_formula(z)
        results.append({
            'z': z,
            'LHS': lhs,
            'RHS': rhs,
            'difference': abs(lhs - rhs),
            'relative_error': abs(lhs - rhs) / abs(lhs) if lhs != 0 else float('inf')
        })
    return results


def legendre_duplication(z: float) -> Tuple[float, float]:
    """
    Legendre's duplication formula:
    Γ(z) Γ(z + 1/2) = 2^{1-2z} √π Γ(2z)
    
    Returns (LHS, RHS) for verification.
    """
    lhs = gamma(z) * gamma(z + 0.5)
    rhs = 2 ** (1 - 2 * z) * sqrt(pi) * gamma(2 * z)
    return (lhs, rhs)


def legendre_duplication_demo(z_values: list = None) -> list:
    """Demonstrate Legendre's duplication formula."""
    if z_values is None:
        z_values = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
    
    results = []
    for z in z_values:
        lhs, rhs = legendre_duplication(z)
        results.append({
            'z': z,
            'LHS': lhs,
            'RHS': rhs,
            'difference': abs(lhs - rhs),
            'relative_error': abs(lhs - rhs) / abs(lhs) if lhs != 0 else float('inf')
        })
    return results


def euler_beta_integral_demo() -> dict:
    """
    Demonstrate Beta function integral representation.
    """
    test_cases = [
        (0.5, 0.5),  # B(1/2, 1/2) = π
        (1, 1),      # B(1, 1) = 1
        (2, 3),      # B(2, 3) = 1/12
        (3, 3),      # B(3, 3) = 1/30
        (0.5, 1.5),  # B(1/2, 3/2) = π/2
    ]
    
    results = []
    for x, y in test_cases:
        exact = beta_function(x, y)
        approx = beta_function_integral(x, y, 50000)
        results.append({
            'x': x,
            'y': y,
            'exact': exact,
            'integral': approx,
            'error': abs(exact - approx)
        })
    return results


def gamma_function_values() -> dict:
    """Return some notable values of the Gamma function."""
    return {
        'Γ(1)': gamma(1),
        'Γ(1/2)': gamma(0.5),
        'Γ(3/2)': gamma(1.5),
        'Γ(2)': gamma(2),
        'Γ(5/2)': gamma(2.5),
        'Γ(3)': gamma(3),
        'Γ(4)': gamma(4),
        'Γ(1/3)': gamma(1/3),
        'Γ(2/3)': gamma(2/3),
        'Γ(1/4)': gamma(1/4),
        'Γ(3/4)': gamma(3/4),
    }


def stirling_approximation(n: int) -> float:
    """
    Stirling's approximation for n! (derived from Gamma function):
    n! ≈ √(2πn) (n/e)^n
    
    More terms:
    n! ≈ √(2πn) (n/e)^n (1 + 1/(12n) + 1/(288n²) - ...)
    """
    if n == 0:
        return 1.0
    return sqrt(2 * pi * n) * (n / exp(1)) ** n


def stirling_approximation_demo() -> dict:
    """Show Stirling's approximation accuracy."""
    results = []
    for n in [1, 2, 3, 5, 10, 20, 50, 100]:
        exact = factorial(n)
        approx = stirling_approximation(n)
        results.append({
            'n': n,
            'exact': exact,
            'approx': approx,
            'relative_error': abs(exact - approx) / exact
        })
    return results


def polygamma_function(n: int, z: float) -> float:
    """
    Polygamma function ψ^{(n)}(z) = d^{n+1}/dz^{n+1} log Γ(z).
    n=0: digamma, n=1: trigamma, etc.
    """
    if n == 0:
        return digamma_function(z)
    
    # For n > 0, use recurrence
    # ψ^{(n)}(z+1) = ψ^{(n)}(z) + (-1)^n n! / z^{n+1}
    # Shift to large z
    result = 0.0
    while z < 10:
        result += (-1) ** n * factorial(n) / (z ** (n + 1))
        z += 1
    
    # Asymptotic expansion for polygamma
    inv_z = 1 / z
    # ψ^{(n)}(z) ≈ (-1)^{n+1} (n-1)! / z^n + ...
    if n == 1:
        result += 1 / z + 1 / (2 * z * z) + 1 / (6 * z ** 3) - 1 / (30 * z ** 5)
    elif n == 2:
        result -= 1 / (z ** 2) + 1 / (z ** 3) + 1 / (2 * z ** 4) - 1 / (6 * z ** 6)
    
    return (-1) ** (n + 1) * result


def euler_infinite_product_gamma(z: float, n: int) -> float:
    """
    Euler's infinite product for Gamma function:
    Γ(z) = lim_{n→∞} n^z n! / (z(z+1)...(z+n))
    """
    result = 1.0
    for k in range(1, n + 1):
        result *= (1 + 1/k) ** z / (1 + z/k)
    return result * (n ** z) / z


def euler_gamma_demo():
    """Demonstrate Euler's Gamma function properties."""
    print("=== Euler's Special Functions ===\n")
    
    print("Gamma Function Values:")
    vals = gamma_function_values()
    for k, v in vals.items():
        print(f"  {k} = {v:.10f}")
    
    print("\nEuler's Reflection Formula Γ(z)Γ(1-z) = π/sin(πz):")
    results = euler_reflection_demo()
    for r in results:
        print(f"  z={r['z']:.1f}: LHS={r['LHS']:.8f}, RHS={r['RHS']:.8f}, rel.err={r['relative_error']:.2e}")
    
    print("\nLegendre's Duplication Formula Γ(z)Γ(z+½) = 2^{1-2z}√π Γ(2z):")
    results = legendre_duplication_demo()
    for r in results:
        print(f"  z={r['z']:.1f}: LHS={r['LHS']:.8f}, RHS={r['RHS']:.8f}, rel.err={r['relative_error']:.2e}")
    
    print("\nBeta Function B(x,y) = Γ(x)Γ(y)/Γ(x+y):")
    results = euler_beta_integral_demo()
    for r in results:
        print(f"  B({r['x']}, {r['y']}) = {r['exact']:.10f} (integral: {r['integral']:.10f}, err={r['error']:.2e})")
    
    print("\nStirling's Approximation n! ≈ √(2πn)(n/e)^n:")
    results = stirling_approximation_demo()
    for r in results:
        print(f"  {r['n']:3d}!: exact={r['exact']}, approx={r['approx']:.2f}, rel.err={r['relative_error']:.2e}")
    
    print("\nEuler's Infinite Product for Γ(0.5):")
    for n in [10, 50, 100, 500, 1000]:
        approx = euler_infinite_product_gamma(0.5, n)
        exact = gamma(0.5)
        print(f"  n={n:4d}: {approx:.10f} (exact={exact:.10f}, err={abs(exact-approx):.2e})")


if __name__ == "__main__":
    euler_gamma_demo()