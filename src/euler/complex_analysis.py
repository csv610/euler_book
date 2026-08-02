"""
Euler's Complex Analysis

This module implements Euler's contributions to complex analysis:
- Euler's formula: e^{ix} = cos(x) + i sin(x)
- Euler's identity: e^{iπ} + 1 = 0
- Complex exponential, logarithm, and powers
- De Moivre's formula
- Roots of unity
- Series expansions

References:
- Euler, L. (1748). "Introductio in analysin infinitorum"
- Euler, L. (1755). "Institutiones calculi differentialis"
"""

from cmath import exp, log, sin, cos, tan, pi, e, sqrt, phase, polar, rect
from math import factorial, sqrt as math_sqrt, cos as math_cos, sin as math_sin, pi as math_pi, e as math_e
from typing import List, Tuple
import cmath

def euler_formula(x: float) -> complex:
    """
    Euler's formula: e^{ix} = cos(x) + i sin(x)
    
    This is Euler's most famous formula, connecting exponential
    and trigonometric functions.
    """
    return exp(1j * x)


def euler_identity() -> complex:
    """
    Euler's identity: e^{iπ} + 1 = 0
    
    Often called the most beautiful equation in mathematics,
    connecting e, i, π, 1, and 0.
    """
    return exp(1j * pi) + 1


def complex_exponential(z: complex) -> complex:
    """Complex exponential function e^z."""
    return exp(z)


def complex_logarithm(z: complex, branch: int = 0) -> complex:
    """
    Complex logarithm (multi-valued).
    Principal branch: Log(z) = ln|z| + i Arg(z) where Arg(z) ∈ (-π, π]
    General branch: log(z) = Log(z) + 2πik for k ∈ ℤ
    """
    return log(z) + 2j * pi * branch


def complex_power(base: complex, exponent: complex) -> complex:
    """
    Complex exponentiation: base^exponent = exp(exponent * log(base)).
    Multi-valued in general.
    """
    return exp(exponent * log(base))


def de_moivre_formula(n: int, theta: float) -> complex:
    """
    De Moivre's formula: (cos θ + i sin θ)^n = cos(nθ) + i sin(nθ)
    
    Derived from Euler's formula: (e^{iθ})^n = e^{inθ}
    """
    return exp(1j * n * theta)


def roots_of_unity(n: int) -> List[complex]:
    """
    The n-th roots of unity: solutions to z^n = 1.
    
    z_k = e^{2πik/n} = cos(2πk/n) + i sin(2πk/n) for k = 0, 1, ..., n-1
    
    These form a regular n-gon on the unit circle.
    """
    return [exp(2j * pi * k / n) for k in range(n)]


def euler_formula_series(x: float, n_terms: int) -> complex:
    """
    Taylor series expansion of e^{ix}:
    e^{ix} = Σ_{k=0}^{∞} (ix)^k / k!
          = 1 + ix - x²/2! - ix³/3! + x⁴/4! + ...
    
    Demonstrates the connection between exponential and trigonometric series.
    """
    result = 0 + 0j
    for k in range(n_terms):
        result += (1j * x) ** k / factorial(k)
    return result


def cos_from_euler(x: float) -> float:
    """cos(x) = (e^{ix} + e^{-ix}) / 2"""
    return (exp(1j * x) + exp(-1j * x)).real / 2


def sin_from_euler(x: float) -> float:
    """sin(x) = (e^{ix} - e^{-ix}) / (2i)"""
    return (exp(1j * x) - exp(-1j * x)).imag / 2


def tan_from_euler(x: float) -> float:
    """tan(x) = sin(x)/cos(x) using Euler's formula"""
    e_ix = exp(1j * x)
    e_minus_ix = exp(-1j * x)
    sin_x = (e_ix - e_minus_ix) / (2j)
    cos_x = (e_ix + e_minus_ix) / 2
    return (sin_x / cos_x).real


def euler_formula_proof_terms(x: float, n_terms: int) -> List[complex]:
    """
    Show the terms in Euler's formula series expansion:
    e^{ix} = Σ (ix)^k/k! = Σ_{even} (-1)^{k/2} x^k/k! + i Σ_{odd} (-1)^{(k-1)/2} x^k/k!
    """
    terms = []
    for k in range(n_terms):
        terms.append((1j * x) ** k / factorial(k))
    return terms


def separate_real_imaginary(terms: List[complex]) -> Tuple[List[float], List[float]]:
    """Separate real and imaginary parts of series terms."""
    real_parts = [t.real for t in terms]
    imag_parts = [t.imag for t in terms]
    return real_parts, imag_parts


def verify_euler_identity() -> dict:
    """
    Verify Euler's identity e^{iπ} + 1 = 0 with increasing precision.
    """
    results = []
    for n in [5, 10, 20, 50, 100]:
        approx = euler_formula_series(pi, n)
        results.append({
            'terms': n,
            'approximation': approx,
            'identity_value': approx + 1,
            'error': abs(approx + 1)
        })
    return results


def roots_of_unity_visualization(n: int) -> List[Tuple[float, float]]:
    """Return (x, y) coordinates of n-th roots of unity for plotting."""
    roots = roots_of_unity(n)
    return [(r.real, r.imag) for r in roots]


def complex_trig_identities() -> dict:
    """
    Derive trigonometric identities using Euler's formula.
    """
    return {
        'cos_sum': 'cos(a+b) = cos(a)cos(b) - sin(a)sin(b)',
        'sin_sum': 'sin(a+b) = sin(a)cos(b) + cos(a)sin(b)',
        'cos_double': 'cos(2x) = 2cos²(x) - 1 = 1 - 2sin²(x)',
        'sin_double': 'sin(2x) = 2sin(x)cos(x)',
        'cos_triple': 'cos(3x) = 4cos³(x) - 3cos(x)',
        'sin_triple': 'sin(3x) = 3sin(x) - 4sin³(x)',
    }


def euler_cosine_product(x: float, n: int) -> float:
    """
    Euler's infinite product for cosine:
    cos(x) = ∏_{k=1}^{∞} (1 - 4x²/((2k-1)²π²))
    """
    product = 1.0
    for k in range(1, n + 1):
        product *= 1 - 4 * x * x / ((2 * k - 1) ** 2 * math_pi ** 2)
    return product


def euler_sine_product(x: float, n: int) -> float:
    """
    Euler's infinite product for sine:
    sin(x) = x ∏_{k=1}^{∞} (1 - x²/(k²π²))
    """
    if x == 0:
        return 0.0
    product = x
    for k in range(1, n + 1):
        product *= 1 - x * x / (k * k * math_pi ** 2)
    return product


def complex_arctan(z: complex) -> complex:
    """Complex arctangent using log: arctan(z) = i/2 * log((i+z)/(i-z))"""
    return 1j / 2 * log((1j + z) / (1j - z))


def complex_arcsin(z: complex) -> complex:
    """Complex arcsin: arcsin(z) = -i log(iz + sqrt(1-z²))"""
    return -1j * log(1j * z + sqrt(1 - z * z))


def complex_arccos(z: complex) -> complex:
    """Complex arccos: arccos(z) = -i log(z + i sqrt(1-z²))"""
    return -1j * log(z + 1j * sqrt(1 - z * z))


# Example usage
if __name__ == "__main__":
    print("=== Euler's Complex Analysis ===\n")
    
    print("Euler's Formula e^{ix} = cos(x) + i sin(x):")
    for x in [0, math_pi/6, math_pi/4, math_pi/3, math_pi/2, math_pi]:
        e_ix = euler_formula(x)
        print(f"  x={x:.4f}: e^{{ix}} = {e_ix:.6f}, cos+isin = {math_cos(x):.6f}+{math_sin(x):.6f}i")
    
    print("\nEuler's Identity e^{iπ} + 1 = 0:")
    e_id = euler_identity()
    print(f"  e^{{iπ}} + 1 = {e_id:.10e}")
    
    print("\nVerification with series:")
    for r in verify_euler_identity():
        print(f"  n={r['terms']:3d}: e^{{iπ}}+1 = {r['identity_value']:.2e}")
    
    print("\nDe Moivre's Formula (cos θ + i sin θ)^n = cos(nθ) + i sin(nθ):")
    for n in [2, 3, 4, 5]:
        for theta in [math_pi/6, math_pi/4, math_pi/3]:
            lhs = de_moivre_formula(n, theta)
            rhs = math_cos(n * theta) + 1j * math_sin(n * theta)
            print(f"  n={n}, θ={theta:.3f}: LHS={lhs:.6f}, RHS={rhs:.6f}, match={abs(lhs-rhs)<1e-10}")
    
    print("\nRoots of Unity:")
    for n in [3, 4, 5, 6, 8]:
        roots = roots_of_unity(n)
        print(f"  {n}-th roots: {[f'{r.real:.3f}+{r.imag:.3f}i' for r in roots]}")
        # Verify they sum to 0 (for n > 1)
        if n > 1:
            print(f"    Sum = {sum(roots):.2e} (should be 0)")
    
    print("\nEuler's Sine Product:")
    for x in [0.5, 1.0, math_pi/2]:
        for n in [5, 10, 50]:
            approx = euler_sine_product(x, n)
            exact = math_sin(x)
            print(f"  x={x:.3f}, n={n}: {approx:.8f} (exact: {exact:.8f})")
            break
        print()
    
    print("\nEuler's Cosine Product:")
    for x in [0.5, math_pi/4, math_pi/3]:
        for n in [5, 10, 50]:
            approx = euler_cosine_product(x, n)
            exact = math_cos(x)
            print(f"  x={x:.3f}, n={n}: {approx:.8f} (exact: {exact:.8f})")
            break
        print()