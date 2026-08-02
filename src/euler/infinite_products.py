"""
Euler's Infinite Products

This module implements Euler's infinite product formulas:
- Euler product formula for the Riemann zeta function
- Wallis product for π/2
- Sinc product formula
- Euler's product formula for sine function
- Product formula for cosine function

References:
- Euler, L. (1737). "Variae observationes circa series infinitas"
- Euler, L. (1748). "Introductio in analysin infinitorum"
"""

from math import pi, sin, cos, prod
from typing import List, Tuple


def prime_sieve(limit: int) -> List[int]:
    """Generate primes up to limit using Sieve of Eratosthenes."""
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i, is_p in enumerate(sieve) if is_p]


def euler_product_formula(s: float, n_primes: int = 1000) -> float:
    """
    Euler product formula for the Riemann zeta function:
    ζ(s) = ∏_{p prime} 1/(1 - p^{-s})  for Re(s) > 1
    
    For s=2: ζ(2) = ∏_p 1/(1 - 1/p²) = π²/6
    
    Args:
        s: Exponent (must be > 1 for convergence)
        n_primes: Number of primes to use
    
    Returns:
        Approximation of ζ(s)
    """
    if s <= 1:
        raise ValueError("Euler product converges only for s > 1")
    
    primes = prime_sieve(10000)[:n_primes]
    product = 1.0
    for p in primes:
        product *= 1 / (1 - p ** (-s))
    return product


def euler_product_partial(s: float, n: int) -> List[float]:
    """Return partial products of Euler product formula."""
    primes = prime_sieve(10000)
    partials = []
    product = 1.0
    for p in primes[:n]:
        product *= 1 / (1 - p ** (-s))
        partials.append(product)
    return partials


def wallis_product(n: int) -> float:
    """
    Wallis product for π/2:
    π/2 = ∏_{n=1}^{∞} (2n * 2n) / ((2n-1)(2n+1))
         = 2/1 * 2/3 * 4/3 * 4/5 * 6/5 * 6/7 * ...
    
    Args:
        n: Number of terms (each term contributes 2 factors)
    
    Returns:
        Approximation of π/2
    """
    product = 1.0
    for k in range(1, n + 1):
        product *= (2 * k) ** 2 / ((2 * k - 1) * (2 * k + 1))
    return product


def wallis_product_partial(n: int) -> List[float]:
    """Return list of partial products for Wallis formula."""
    partials = []
    product = 1.0
    for k in range(1, n + 1):
        product *= (2 * k) ** 2 / ((2 * k - 1) * (2 * k + 1))
        partials.append(product)
    return partials


def sinc_product(x: float, n: int) -> float:
    """
    Sinc product formula:
    sinc(x) = sin(x)/x = ∏_{k=1}^{∞} (1 - x²/(k²π²))
    
    This is Euler's product formula for sin(x)/x.
    """
    if x == 0:
        return 1.0
    product = 1.0
    for k in range(1, n + 1):
        product *= 1 - x * x / (k * k * pi * pi)
    return product


def sinc_product_partial(x: float, n: int) -> List[float]:
    """Return partial products for sinc product."""
    if x == 0:
        return [1.0] * n
    partials = []
    product = 1.0
    for k in range(1, n + 1):
        product *= 1 - x * x / (k * k * pi * pi)
        partials.append(product)
    return partials


def euler_product_sine(x: float, n: int) -> float:
    """
    Euler's product formula for sine:
    sin(x) = x ∏_{k=1}^{∞} (1 - x²/(k²π²))
    
    This was Euler's key insight for solving the Basel problem.
    """
    if x == 0:
        return 0.0
    product = x
    for k in range(1, n + 1):
        product *= 1 - x * x / (k * k * pi * pi)
    return product


def euler_product_sine_partial(x: float, n: int) -> List[float]:
    """Return partial products for Euler's sine product."""
    if x == 0:
        return [0.0] * n
    partials = []
    product = x
    for k in range(1, n + 1):
        product *= 1 - x * x / (k * k * pi * pi)
        partials.append(product)
    return partials


def euler_cosine_product(x: float, n: int) -> float:
    """
    Euler's infinite product for cosine:
    cos(x) = ∏_{k=1}^{∞} (1 - 4x²/((2k-1)²π²))
    """
    product = 1.0
    for k in range(1, n + 1):
        product *= 1 - 4 * x * x / ((2 * k - 1) ** 2 * pi * pi)
    return product


def euler_cosine_product_partial(x: float, n: int) -> List[float]:
    """Return partial products for Euler's cosine product."""
    partials = []
    product = 1.0
    for k in range(1, n + 1):
        product *= 1 - 4 * x * x / ((2 * k - 1) ** 2 * pi * pi)
        partials.append(product)
    return partials


def euler_product_zeta_even(n: int, n_primes: int = 1000) -> float:
    """
    Euler product for ζ(2n):
    ζ(2n) = ∏_p 1/(1 - p^{-2n})
    
    Compare with Euler's exact formula:
    ζ(2n) = (-1)^{n+1} B_{2n} (2π)^{2n} / (2(2n)!)
    """
    primes = prime_sieve(10000)[:n_primes]
    product = 1.0
    for p in primes:
        product *= 1 / (1 - p ** (-2 * n))
    return product


def viete_formula(n: int) -> float:
    """
    Viète's formula for π (discovered before Euler, but related):
    2/π = ∏_{k=1}^{∞} cos(π/2^{k+1})
        = cos(π/4) * cos(π/8) * cos(π/16) * ...
    """
    product = 1.0
    for k in range(1, n + 1):
        product *= cos(pi / (2 ** (k + 1)))
    return 2 / product


def leibniz_product(n: int) -> float:
    """
    Product form of Leibniz series:
    π/4 = ∏_{k=1}^{∞} (1 - 1/(2k+1)²)^{-1} ... not quite.
    Actually this is related to the wallis product.
    """
    return wallis_product(n) / 2


def euler_product_demonstration():
    """Demonstrate convergence of various infinite products."""
    print("=== Euler's Infinite Products ===\n")
    
    # Wallis product
    print("Wallis Product for π/2:")
    for n in [1, 2, 5, 10, 50, 100, 1000]:
        approx = wallis_product(n)
        exact = pi / 2
        print(f"  n={n:4d}: {approx:.10f} (π/2={exact:.10f}, error={abs(exact-approx):.2e})")
    
    print("\nEuler Product for ζ(2) = π²/6:")
    for n in [10, 50, 100, 500, 1000]:
        approx = euler_product_formula(2, n)
        exact = pi * pi / 6
        print(f"  {n:4d} primes: {approx:.10f} (exact={exact:.10f}, error={abs(exact-approx):.2e})")
    
    print("\nEuler's Sine Product sin(x) = x∏(1 - x²/(k²π²)):")
    for x in [0.5, 1.0, pi/4, pi/2]:
        for n in [5, 10, 50]:
            approx = euler_product_sine(x, n)
            exact = sin(x)
            print(f"  x={x:.3f}, n={n:2d}: {approx:.8f} (exact={exact:.8f}, error={abs(exact-approx):.2e})")
            break
        print()
    
    print("\nEuler's Cosine Product cos(x) = ∏(1 - 4x²/((2k-1)²π²)):")
    for x in [0.5, pi/6, pi/4]:
        for n in [5, 10, 50]:
            approx = euler_cosine_product(x, n)
            exact = cos(x)
            print(f"  x={x:.3f}, n={n:2d}: {approx:.8f} (exact={exact:.8f}, error={abs(exact-approx):.2e})")
            break
        print()
    
    print("\nViète's Formula 2/π = ∏cos(π/2^{k+1}):")
    for n in [5, 10, 20, 50]:
        approx = viete_formula(n)
        exact = pi
        print(f"  n={n:2d}: π ≈ {approx:.10f} (error={abs(exact-approx):.2e})")


if __name__ == "__main__":
    euler_product_demonstration()