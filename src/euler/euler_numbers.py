"""
Euler Numbers and Tangent/Secant Expansions

This module implements Euler's work on tangent and secant number expansions,
including the Euler zigzag numbers (up/down permutations) and the classical
Euler numbers used in Taylor series of tan(x) and sec(x).

Euler discovered that the Taylor series of tan(x) and sec(x) have integer
coefficients related to counting alternating permutations. These numbers
appear throughout combinatorics and analysis.

References:
    - Euler, L. (1755). Institutiones calculi differentialis.
    - Knuth, T. & Buckholz, P. (1969). "Computation of Tangent, Euler, and Secant Numbers."
"""

from math import factorial, sin, cos, pi
from typing import List, Tuple


def euler_number_zigzag(n: int) -> int:
    """
    Compute the n-th Euler zigzag number (A_n).

    These count the number of alternating permutations of {1, 2, ..., n}.
    An alternating permutation satisfies a_1 < a_2 > a_3 < a_4 > ...

    A_0 = 1, A_1 = 1, A_2 = 1, A_3 = 2, A_4 = 5, A_5 = 16, A_6 = 61,
    A_7 = 272, A_8 = 1385, A_9 = 7936, A_10 = 50521

    Uses the Entringer triangle DP method.

    Args:
        n: non-negative integer

    Returns:
        The n-th Euler zigzag number

    Raises:
        ValueError: if n < 0
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    if n <= 1:
        return 1

    # 2*A_{n} = sum_{k=0}^{n-1} binomial(n-1,k) * A_k * A_{n-1-k}
    from math import comb

    A = [1, 1]  # A_0=1, A_1=1
    for nn in range(1, n):
        total = 0
        for k in range(nn + 1):
            total += comb(nn, k) * A[k] * A[nn - k]
        A.append(total // 2)

    return A[n]


def euler_number_classical(n: int) -> int:
    """
    Compute the n-th classical Euler number E_n (secant number).

    These are the coefficients in the Taylor expansion of sec(x):
    sec(x) = sum_{n=0}^{inf} E_{2n} * x^{2n} / (2n)!

    E_0 = 1, E_2 = 1, E_4 = 5, E_6 = 61, E_8 = 1385, ...
    Odd-indexed Euler numbers are 0.

    Args:
        n: non-negative integer (even index for non-zero values)

    Returns:
        The n-th Euler number

    Raises:
        ValueError: if n < 0
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    if n % 2 == 1:
        return 0

    return euler_number_zigzag(n)


def tangent_taylor(x: float, terms: int = 10) -> float:
    """
    Compute tan(x) using its Taylor series with Euler numbers.

    tan(x) = x + x^3/3 + 2*x^5/15 + 17*x^7/315 + 62*x^9/2835 + ...

    The coefficients relate to Euler zigzag numbers:
    tan(x) + sec(x) = sum_{n=0}^{inf} A_n * x^n / n!

    Args:
        x: angle in radians
        terms: number of terms to use

    Returns:
        Approximation of tan(x)
    """
    result = 0.0
    for n in range(1, terms + 1):
        e_idx = 2 * n - 1  # odd indices for tan
        e_val = euler_number_zigzag(e_idx)
        coeff = e_val / factorial(e_idx)
        result += coeff * (x ** e_idx)

    return result


def secant_taylor(x: float, terms: int = 10) -> float:
    """
    Compute sec(x) using its Taylor series with Euler numbers.

    sec(x) = 1 + x^2/2 + 5*x^4/24 + 61*x^6/720 + 1385*x^8/40320 + ...

    Args:
        x: angle in radians
        terms: number of terms to use

    Returns:
        Approximation of sec(x)
    """
    result = 0.0
    for n in range(terms):
        e_idx = 2 * n  # even indices for secant
        e_val = euler_number_zigzag(e_idx)
        coeff = e_val / factorial(e_idx)
        result += coeff * (x ** e_idx)

    return result


def alternating_permutations_count(n: int) -> int:
    """
    Count the number of alternating (up-down) permutations of {1, 2, ..., n}.

    Also known as Zigzag numbers or Andre permutations.
    A_0 = 1, A_1 = 1, A_2 = 1, A_3 = 2, A_4 = 5, A_5 = 16, A_6 = 61, ...

    Args:
        n: positive integer

    Returns:
        Number of alternating permutations of length n
    """
    return euler_number_zigzag(n)


def entringer_triangle(rows: int) -> List[List[int]]:
    """
    Construct the Entringer triangle (Seidel-Entringer-Arnold triangle).

    The triangle is built as follows:
    Row 0: [1]
    Row n: E[n][k] = E[n][k-1] + E[n-1][n-k] for k = 1, ..., n
           E[n][0] = E[n-1][0]

    The leftmost entry of each row gives the Euler zigzag number.

    Args:
        rows: number of rows to compute

    Returns:
        List of rows, where each row is a list of integers
    """
    if rows <= 0:
        return []

    triangle = [[1]]  # Row 0

    for n in range(1, rows):
        prev = triangle[-1]
        curr = [0] * (n + 1)
        curr[0] = prev[0]
        for k in range(1, n + 1):
            idx = n - k
            if idx < len(prev):
                curr[k] = curr[k - 1] + prev[idx]
            else:
                curr[k] = curr[k - 1]
        triangle.append(curr)

    return triangle


def euler_secant_numbers(n_max: int) -> List[int]:
    """
    Compute Euler secant numbers E_0, E_2, E_4, ..., E_{2n_max}.

    These are the non-zero Euler numbers (even indices).
    E_0=1, E_2=1, E_4=5, E_6=61, E_8=1385, E_10=50521, ...

    Args:
        n_max: maximum half-index (returns E_0 through E_{2*n_max})

    Returns:
        List of Euler secant numbers
    """
    return [euler_number_classical(2 * n) for n in range(n_max + 1)]


def verify_taylor_accuracy(tol: float = 1e-10, x_max: float = 1.0) -> dict:
    """
    Verify the accuracy of Taylor series approximations for tan and sec.

    Args:
        tol: tolerance for comparison
        x_max: maximum x value to test

    Returns:
        Dictionary with accuracy information
    """
    import math

    results = {
        'tan_tests': [],
        'sec_tests': [],
        'all_passed': True
    }

    test_points = [0.1, 0.3, 0.5, 0.7, 0.9]

    for x in test_points:
        if x > x_max:
            break

        # Test tan
        tan_exact = math.tan(x)
        for terms in [3, 5, 7, 10]:
            tan_approx = tangent_taylor(x, terms)
            error = abs(tan_approx - tan_exact)
            passed = error < tol
            if not passed:
                results['all_passed'] = False
            results['tan_tests'].append({
                'x': x, 'terms': terms, 'exact': tan_exact,
                'approx': tan_approx, 'error': error, 'passed': passed
            })

        # Test sec
        sec_exact = 1.0 / math.cos(x)
        for terms in [3, 5, 7, 10]:
            sec_approx = secant_taylor(x, terms)
            error = abs(sec_approx - sec_exact)
            passed = error < tol
            if not passed:
                results['all_passed'] = False
            results['sec_tests'].append({
                'x': x, 'terms': terms, 'exact': sec_exact,
                'approx': sec_approx, 'error': error, 'passed': passed
            })

    return results


def euler_number_demo():
    """Demonstrate Euler zigzag and secant numbers."""
    print("=== Euler Zigzag Numbers (Alternating Permutations) ===\n")

    print("n  | A_n  | Description")
    print("-" * 40)
    for n in range(12):
        e = euler_number_zigzag(n)
        print(f"{n:2d} | {e:4d}  |")

    print(f"\nFirst 12 Euler zigzag numbers: ", end="")
    print([euler_number_zigzag(n) for n in range(12)])

    print("\n=== Entringer Triangle ===\n")
    triangle = entringer_triangle(8)
    for i, row in enumerate(triangle):
        print(f"Row {i}: {row}")

    print("\n=== Euler Secant Numbers (Even Indices) ===\n")
    secs = euler_secant_numbers(8)
    print(f"E_0, E_2, ..., E_16: {secs}")

    print("\n=== Taylor Series Verification ===\n")
    import math
    x = 0.5
    print(f"x = {x}")
    print(f"tan({x}) = {math.tan(x):.12f}")
    for terms in [3, 5, 7, 10]:
        approx = tangent_taylor(x, terms)
        print(f"  tan({x}) [{terms} terms] = {approx:.12f}  (error: {abs(math.tan(x) - approx):.2e})")

    print(f"\nsec({x}) = {1/math.cos(x):.12f}")
    for terms in [3, 5, 7, 10]:
        approx = secant_taylor(x, terms)
        print(f"  sec({x}) [{terms} terms] = {approx:.12f}  (error: {abs(1/math.cos(x) - approx):.2e})")


if __name__ == "__main__":
    euler_number_demo()
