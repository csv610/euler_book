"""
Euler's Calculus and Numerical Methods

This module implements Euler's contributions to calculus and numerical analysis:
- Euler's method for ODEs
- Improved Euler (Heun's) method
- Runge-Kutta 4th order method
- Euler-Maclaurin summation formula
- Euler's constant γ

References:
- Euler, L. (1768). "Institutiones calculi integralis"
- Euler, L. (1755). "Institutiones calculi differentialis"
"""

from math import factorial, log, exp, pi, e, isclose
from typing import Callable, List, Tuple, Any, Optional


def euler_method(
    f: Callable[[float, float], float],
    x0: float,
    y0: float,
    h: float,
    n: int
) -> List[Tuple[float, float]]:
    """
    Euler's method for solving dy/dx = f(x, y) with y(x0) = y0.
    
    y_{n+1} = y_n + h * f(x_n, y_n)
    x_{n+1} = x_n + h
    
    Args:
        f: Function f(x, y) = dy/dx
        x0: Initial x value
        y0: Initial y value
        h: Step size
        n: Number of steps
    
    Returns:
        List of (x, y) pairs
    """
    results = [(x0, y0)]
    x, y = x0, y0
    
    for _ in range(n):
        y = y + h * f(x, y)
        x = x + h
        results.append((x, y))
    
    return results


def euler_method_system(
    f: Callable[[float, List[float]], List[float]],
    x0: float,
    y0: List[float],
    h: float,
    n: int
) -> List[Tuple[float, List[float]]]:
    """
    Euler's method for systems of ODEs.
    
    dy/dx = f(x, y) where y is a vector.
    """
    results = [(x0, y0[:])]
    x = x0
    y = y0[:]
    
    for _ in range(n):
        dy = f(x, y)
        y = [y[i] + h * dy[i] for i in range(len(y))]
        x = x + h
        results.append((x, y[:]))
    
    return results


def improved_euler_method(
    f: Callable[[float, float], float],
    x0: float,
    y0: float,
    h: float,
    n: int
) -> List[Tuple[float, float]]:
    """
    Improved Euler method (Heun's method) - 2nd order Runge-Kutta.
    
    Predictor: y* = y_n + h * f(x_n, y_n)
    Corrector: y_{n+1} = y_n + h/2 * (f(x_n, y_n) + f(x_{n+1}, y*))
    """
    results = [(x0, y0)]
    x, y = x0, y0
    
    for _ in range(n):
        k1 = f(x, y)
        y_star = y + h * k1
        k2 = f(x + h, y_star)
        y = y + h * (k1 + k2) / 2
        x = x + h
        results.append((x, y))
    
    return results


def runge_kutta_4(
    f: Callable[[float, float], float],
    x0: float,
    y0: float,
    h: float,
    n: int
) -> List[Tuple[float, float]]:
    """
    Classical 4th order Runge-Kutta method.
    
    k1 = h * f(x_n, y_n)
    k2 = h * f(x_n + h/2, y_n + k1/2)
    k3 = h * f(x_n + h/2, y_n + k2/2)
    k4 = h * f(x_n + h, y_n + k3)
    y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4) / 6
    """
    results = [(x0, y0)]
    x, y = x0, y0
    
    for _ in range(n):
        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        x = x + h
        results.append((x, y))
    
    return results


def euler_maclaurin_summation(
    f: Callable[[int], float],
    a: int,
    b: int,
    k: int = 10
) -> float:
    """
    Euler-Maclaurin summation formula:
    ∑_{n=a}^b f(n) = ∫_a^b f(x) dx + (f(a) + f(b))/2 
                    + ∑_{j=1}^k B_{2j}/(2j)! * (f^{(2j-1)}(b) - f^{(2j-1)}(a))
                    + R_k
    
    Args:
        f: Function to sum
        a: Lower bound (inclusive)
        b: Upper bound (inclusive)
        k: Number of Bernoulli correction terms
    
    Returns:
        Approximation of the sum
    """
    # Bernoulli numbers B_2, B_4, B_6, B_8, B_10, B_12, B_14, B_16
    bernoulli = [
        1/6,      # B_2
        -1/30,    # B_4
        1/42,     # B_6
        -1/30,    # B_8
        5/66,     # B_10
        -691/2730,# B_12
        7/6,      # B_14
        -3617/510,# B_16
    ]
    
    # Direct sum for small ranges
    if b - a < 100:
        return sum(f(i) for i in range(a, b + 1))
    
    # Numerical integration using trapezoidal rule
    def integral_approx():
        n_segments = 1000
        dx = (b - a) / n_segments
        total = (f(a) + f(b)) / 2
        for i in range(1, n_segments):
            total += f(a + i * dx)
        return total * dx
    
    integral = integral_approx()
    
    # Boundary correction
    boundary = (f(a) + f(b)) / 2
    
    # Bernoulli corrections (simplified - need derivatives)
    # For demonstration, we just return integral + boundary
    # Full implementation requires symbolic or numerical derivatives
    
    return integral + boundary


def euler_maclaurin_integral(
    f: Callable[[float], float],
    a: float,
    b: float,
    n: int = 1000
) -> float:
    """
    Euler-Maclaurin formula for approximating integrals:
    ∫_a^b f(x) dx ≈ (h/2)(f(a) + f(b)) + h ∑_{k=1}^{n-1} f(a + kh)
    where h = (b-a)/n
    
    This is the composite trapezoidal rule.
    """
    h = (b - a) / n
    total = (f(a) + f(b)) / 2
    for k in range(1, n):
        total += f(a + k * h)
    return total * h


def euler_constant_limit(n: int) -> float:
    """
    Euler's constant γ as a limit:
    γ = lim_{n→∞} (H_n - ln n)
    where H_n is the n-th harmonic number.
    """
    harmonic = sum(1.0 / k for k in range(1, n + 1))
    return harmonic - log(n)


def euler_constant_series(n: int) -> float:
    """
    Faster converging series for γ using Gregory coefficients:
    γ = ∑_{k=1}^∞ G_k / k
    where G_k are Gregory coefficients.
    """
    # Gregory coefficients
    gregory = [
        0,
        1/2,
        -1/12,
        1/24,
        -19/720,
        3/160,
        -863/60480,
        275/24192,
        -33953/3628800,
    ]
    
    return sum(gregory[k] / k for k in range(1, min(n + 1, len(gregory))))


def euler_constant_integral(n: int = 1000) -> float:
    """
    Integral representation of γ:
    γ = ∫_0^1 (1 - e^{-t})/t dt - ∫_1^∞ e^{-t}/t dt
    """
    # First integral from 0 to 1
    h = 1.0 / n
    int1 = 0.0
    for i in range(n):
        t = (i + 0.5) * h
        if t == 0:
            val = 1.0  # limit as t->0
        else:
            val = (1 - exp(-t)) / t
        int1 += val * h
    
    # Second integral from 1 to infinity (approximate to 10)
    h2 = 9.0 / n
    int2 = 0.0
    for i in range(n):
        t = 1 + (i + 0.5) * h2
        val = exp(-t) / t
        int2 += val * h2
    
    return int1 - int2


def harmonic_number(n: int) -> float:
    """H_n = 1 + 1/2 + ... + 1/n"""
    return sum(1.0 / k for k in range(1, n + 1))


def digamma_function(x: float, n: int = 1000) -> float:
    """
    Digamma function ψ(x) = Γ'(x)/Γ(x)
    For large x: ψ(x) ≈ ln(x) - 1/(2x) - 1/(12x²) + ...
    """
    # Use recurrence to shift to large x
    while x < 8:
        x += 1
    
    # Asymptotic expansion
    inv_x = 1 / x
    result = log(x) - inv_x / 2 - inv_x * inv_x / 12 + inv_x**4 / 120 - inv_x**6 / 252
    
    # Add back the recurrence terms
    while x > 8:
        x -= 1
        result -= 1 / x
    
    return result


def bernoulli_numbers(n: int) -> List[float]:
    """
    Compute Bernoulli numbers B_0, B_1, ..., B_n using Akiyama-Tanigawa algorithm.
    """
    A = [0.0] * (n + 1)
    for m in range(n + 1):
        A[m] = 1.0 / (m + 1)
        for j in range(m, 0, -1):
            A[j-1] = j * (A[j-1] - A[j])
    return A


def euler_maclaurin_summation_full(
    f: Callable[[float], float],
    f_derivatives: List[Callable[[float], float]],
    a: int,
    b: int,
    k: int = 8
) -> float:
    """
    Full Euler-Maclaurin summation with Bernoulli corrections.
    
    f_derivatives[j] should be the (2j-1)-th derivative of f.
    """
    B = bernoulli_numbers(2 * k)
    
    # Direct sum for small range
    if b - a < 50:
        return sum(f(i) for i in range(a, b + 1))
    
    # Integral term (trapezoidal rule approximation)
    integral = euler_maclaurin_integral(f, a, b, max(1000, (b - a) * 10))
    
    # Boundary term
    boundary = (f(a) + f(b)) / 2
    
    # Bernoulli corrections
    bernoulli_sum = 0.0
    for j in range(1, k + 1):
        if j < len(f_derivatives):
            deriv_a = f_derivatives[j-1](a)
            deriv_b = f_derivatives[j-1](b)
            bernoulli_sum += B[2*j] / factorial(2*j) * (deriv_b - deriv_a)
    
    return integral + boundary + bernoulli_sum


# Example usage
if __name__ == "__main__":
    print("=== Euler's Calculus and Numerical Methods ===\n")
    
    # Test ODE: dy/dx = y, y(0) = 1, exact: y = e^x
    def f(x, y):
        return y
    
    exact = lambda x: exp(x)
    
    print("Euler's Method for dy/dx = y, y(0)=1:")
    for h in [0.5, 0.2, 0.1, 0.05]:
        sol = euler_method(f, 0, 1, h, int(1/h))
        x, y = sol[-1]
        print(f"  h={h:.2f}: y(1)≈{y:.6f}, exact={exact(1):.6f}, error={abs(y-exact(1)):.2e}")
    
    print("\nImproved Euler (Heun) Method:")
    for h in [0.5, 0.2, 0.1, 0.05]:
        sol = improved_euler_method(f, 0, 1, h, int(1/h))
        x, y = sol[-1]
        print(f"  h={h:.2f}: y(1)≈{y:.6f}, exact={exact(1):.6f}, error={abs(y-exact(1)):.2e}")
    
    print("\nRunge-Kutta 4th Order:")
    for h in [0.5, 0.2, 0.1, 0.05]:
        sol = runge_kutta_4(f, 0, 1, h, int(1/h))
        x, y = sol[-1]
        print(f"  h={h:.2f}: y(1)≈{y:.10f}, exact={exact(1):.10f}, error={abs(y-exact(1)):.2e}")
    
    print("\nEuler's Constant γ:")
    for n in [10, 100, 1000, 10000, 100000]:
        gamma = euler_constant_limit(n)
        print(f"  n={n:6d}: γ ≈ {gamma:.10f}")
    print(f"  Series (10 terms): {euler_constant_series(10):.10f}")
    print(f"  Integral: {euler_constant_integral(10000):.10f}")
    print(f"  Known value: 0.5772156649")
    
    print("\nEuler-Maclaurin for ∑ 1/k² (Basel problem):")
    # Sum 1/k² from 1 to N
    def f(k):
        return 1.0 / (k * k)
    
    for N in [10, 100, 1000]:
        direct = sum(f(k) for k in range(1, N + 1))
        approx = euler_maclaurin_summation(f, 1, N, 5)
        print(f"  N={N}: direct={direct:.10f}, EM={approx:.10f}, exact(π²/6)={pi*pi/6:.10f}")