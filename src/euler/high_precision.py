"""
High-Precision Computation for Euler's Constants

This module implements arbitrary-precision arithmetic using pure Python integers,
with no external dependencies (no mpmath, no gmpy2). It provides:
- Arbitrary-precision integer arithmetic scaled to fixed decimal places
- High-precision computation of pi, e, gamma
- Machin-like formulas for pi
- Continued fraction and series methods for e
- Euler-Maclaurin acceleration for gamma

All algorithms are based on Euler's original series and products,
accelerated with modern convergence techniques.

References:
    - Brent, R. (1976). "Multiple-precision zero-finding methods."
    - Borwein, J. (1989). "Computer Experimentation in Number Theory."
    - Euler, L. (1734). "De progressionibus transcendentibus."
"""

from math import factorial, comb, floor, log10, sqrt as math_sqrt
from typing import List, Tuple, Optional


# ======================================================================
# Core: Scaled Integer Arithmetic
# ======================================================================

class ScaledInt:
    """
    Arbitrary-precision decimal using scaled integers.

    Value = integer / scale_factor, where scale_factor = 10^precision.
    All arithmetic is exact integer arithmetic.
    """

    def __init__(self, value, precision: int = 50):
        self.precision = precision
        self.scale = 10 ** precision

        if isinstance(value, ScaledInt):
            # Convert to common precision
            if value.precision < precision:
                self.integer = value.integer * (10 ** (precision - value.precision))
            elif value.precision > precision:
                self.integer = value.integer // (10 ** (value.precision - precision))
            else:
                self.integer = value.integer
            self.precision = precision
            self.scale = 10 ** precision
            return

        if isinstance(value, int):
            self.integer = value * self.scale
            return

        if isinstance(value, float):
            # Convert via string to avoid float imprecision
            s = f"{value:.{precision + 5}g}"
            self._from_string(s)
            return

        if isinstance(value, str):
            self._from_string(value)
            return

        raise ValueError(f"Cannot create ScaledInt from {type(value)}")

    def _from_string(self, s: str):
        """Parse a string representation."""
        s = s.strip()
        negative = s.startswith('-')
        if negative:
            s = s[1:]

        parts = s.split('.')
        int_part = parts[0] if parts[0] else '0'
        frac_part = parts[1] if len(parts) > 1 else ''

        # Combine into integer
        combined = int_part.lstrip('0') + frac_part
        if not combined:
            combined = '0'

        self.integer = int(combined)
        if negative:
            self.integer = -self.integer

        # Adjust for precision
        frac_len = len(frac_part)
        if frac_len > self.precision:
            # Truncate
            self.integer //= (10 ** (frac_len - self.precision))
        elif frac_len < self.precision:
            # Pad
            self.integer *= (10 ** (self.precision - frac_len))

    def __str__(self):
        """String representation with decimal point."""
        neg = self.integer < 0
        abs_int = abs(self.integer)

        # Pad with leading zeros
        s = str(abs_int).zfill(self.precision + 1)

        # Insert decimal point
        int_part = s[:-self.precision] if self.precision > 0 else s
        frac_part = s[-self.precision:] if self.precision > 0 else ''

        # Remove trailing zeros from frac
        frac_part = frac_part.rstrip('0')

        result = int_part + ('.' + frac_part if frac_part else '')
        if neg and result != '0':
            result = '-' + result
        return result

    def __repr__(self):
        return f"ScaledInt({str(self)}, precision={self.precision})"

    def __add__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)

        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale

        # Convert to common precision
        if other.precision < self.precision:
            other_scaled = other.integer * (10 ** (self.precision - other.precision))
        elif other.precision > self.precision:
            self_scaled = self.integer * (10 ** (other.precision - self.precision))
            other_scaled = other.integer
            result.integer = self_scaled + other_scaled
            result.precision = other.precision
            result.scale = 10 ** result.precision
            return result
        else:
            other_scaled = other.integer
            self_scaled = self.integer

        result.integer = self_scaled + other_scaled
        return result

    def __sub__(self, other):
        return self + (-ScaledInt(other, self.precision))

    def __neg__(self):
        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale
        result.integer = -self.integer
        return result

    def __mul__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)

        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale

        # Multiply then divide by scale
        product = self.integer * other.integer
        result.integer = product // self.scale
        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)

        if other.integer == 0:
            raise ZeroDivisionError("division by zero")

        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale

        # Multiply numerator by scale, then divide
        dividend = self.integer * self.scale
        result.integer = dividend // other.integer
        return result

    def __floordiv__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)
        if other.integer == 0:
            raise ZeroDivisionError("division by zero")
        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale
        result.integer = (self.integer * self.scale) // (other.integer * self.scale)
        # Simplify
        result.integer = self.integer // other.integer
        return result

    def __lt__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)
        return self.integer < other.integer

    def __le__(self, other):
        return self.integer <= other.integer

    def __gt__(self, other):
        return self.integer > other.integer

    def __ge__(self, other):
        return self.integer >= other.integer

    def __eq__(self, other):
        if not isinstance(other, ScaledInt):
            other = ScaledInt(other, self.precision)
        return self.integer == other.integer

    def __abs__(self):
        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale
        result.integer = abs(self.integer)
        return result

    def sqrt(self) -> 'ScaledInt':
        """Compute square root using integer Newton's method."""
        if self.integer < 0:
            raise ValueError("Cannot take sqrt of negative number")
        if self.integer == 0:
            return ScaledInt(0, self.precision)

        # Newton's method: x_{n+1} = (x_n + n/x_n) / 2
        # Start with integer sqrt estimate
        x = self.integer
        # Initial guess
        guess = int(math_sqrt(abs(x / self.scale))) * self.scale
        if guess == 0:
            guess = self.scale  # at least 1

        for _ in range(50):
            new_guess = (guess + (x * self.scale) // guess) // 2
            if new_guess == guess:
                break
            guess = new_guess

        result = ScaledInt.__new__(ScaledInt)
        result.precision = self.precision
        result.scale = self.scale
        result.integer = guess
        return result

    def power(self, n: int) -> 'ScaledInt':
        """Raise to integer power."""
        if n < 0:
            return ScaledInt(1, self.precision) / self ** (-n)
        if n == 0:
            return ScaledInt(1, self.precision)
        if n == 1:
            return ScaledInt(self, self.precision)

        result = ScaledInt(1, self.precision)
        base = ScaledInt(self, self.precision)
        while n > 0:
            if n % 2 == 1:
                result = result * base
            base = base * base
            n //= 2
        return result

    def __pow__(self, n: int):
        return self.power(n)

    def reciprocal(self) -> 'ScaledInt':
        """Compute 1/self using Newton's method."""
        if self.integer == 0:
            raise ZeroDivisionError("division by zero")

        # Initial guess
        float_val = float(str(self))
        x = ScaledInt(1.0 / float_val, self.precision)

        # Newton iterations: x_{n+1} = x_n * (2 - self * x_n)
        two = ScaledInt(2, self.precision)
        for _ in range(10):
            x = x * (two - self * x)

        return x


# ======================================================================
# High-Precision Pi
# ======================================================================

def pi_machin(digits: int = 50) -> str:
    """
    Compute pi using Machin's formula: pi/4 = 4*arctan(1/5) - arctan(1/239).

    arctan(x) = x - x^3/3 + x^5/5 - x^7/7 + ...

    Args:
        digits: number of decimal digits

    Returns:
        String representation of pi
    """
    prec = digits + 10
    hp = ScaledInt

    def arctan_inv(n: int, max_terms: int) -> ScaledInt:
        """Compute arctan(1/n) using Taylor series."""
        result = hp(0, prec)
        x = hp(1, prec) / hp(n, prec)
        x_sq = x * x
        power = x

        for k in range(max_terms):
            sign = 1 if k % 2 == 0 else -1
            term = power / hp(2 * k + 1, prec)
            if sign > 0:
                result = result + term
            else:
                result = result - term
            power = power * x_sq

        return result

    # Number of terms: each term reduces error by n^2
    # For 1/5: error ~ (1/25)^k, need k ~ digits/log10(25)
    terms_5 = int(digits * 2.5) + 10
    terms_239 = int(digits * 1.5) + 10

    atan1_5 = arctan_inv(5, terms_5)
    atan1_239 = arctan_inv(239, terms_239)

    four = hp(4, prec)
    sixteen = hp(16, prec)
    pi_val = sixteen * atan1_5 - four * atan1_239
    return str(pi_val)


def pi_gauss_agm(digits: int = 50) -> str:
    """
    Compute pi using the Brent-Salamin (Gauss-Legendre) AGM algorithm.

    This converges quadratically: each iteration doubles the number of
    correct digits. Much faster than Machin's formula.

    Algorithm:
    a_0 = 1, b_0 = 1/sqrt(2), t_0 = 1/4, p_0 = 1
    a_{n+1} = (a_n + b_n)/2
    b_{n+1} = sqrt(a_n * b_n)
    t_{n+1} = t_n - p_n * (a_n - a_{n+1})^2
    p_{n+1} = 2 * p_n
    pi ≈ (a_n + b_n)^2 / (4 * t_n)

    Args:
        digits: number of decimal digits

    Returns:
        String representation of pi
    """
    prec = digits + 10
    hp = ScaledInt

    # a = 1
    a = hp(1, prec)
    # b = 1/sqrt(2) = sqrt(2)/2
    sqrt2 = hp(2, prec).sqrt()
    b = sqrt2 / hp(2, prec)
    # t = 1/4
    t = hp(1, prec) / hp(4, prec)
    # p = 1
    p = hp(1, prec)

    # Iterate until convergence (quadratic: ~log2(digits) iterations)
    n_iter = int(digits / 3) + 5

    for _ in range(n_iter):
        a_new = (a + b) / hp(2, prec)
        b_new = hp.sqrt(a * b)
        diff = a - a_new
        t_new = t - p * diff * diff
        p_new = p * hp(2, prec)
        a, b, t, p = a_new, b_new, t_new, p_new

    # pi ≈ (a + b)^2 / (4t)
    ab = a + b
    pi_val = ab * ab / (hp(4, prec) * t)
    return str(pi_val)


# ======================================================================
# High-Precision e
# ======================================================================

def e_series(digits: int = 50) -> str:
    """
    Compute e using the series: e = sum_{n=0}^{inf} 1/n!

    Converges rapidly: each term adds ~1 digit.

    Args:
        digits: number of decimal digits

    Returns:
        String representation of e
    """
    prec = digits + 10
    hp = ScaledInt

    result = hp(0, prec)
    factorial_n = hp(1, prec)

    for n in range(prec + 20):
        if n > 0:
            factorial_n = factorial_n * hp(n, prec)
        term = hp(1, prec) / factorial_n
        result = result + term

    return str(result)


def e_continued_fraction(digits: int = 50) -> str:
    """
    Compute e using its continued fraction expansion.

    e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, ...]
    Pattern: 1, 2k, 1, 1, 2(k+1), 1, ...

    Converges much faster than the series.

    Args:
        digits: number of decimal digits

    Returns:
        String representation of e
    """
    prec = digits + 10
    hp = ScaledInt

    # Build continued fraction coefficients
    n_terms = digits * 2 + 10

    coeffs = [2]
    k = 1
    while len(coeffs) < n_terms:
        coeffs.extend([1, 2 * k, 1])
        k += 1

    # Evaluate from bottom
    result = hp(coeffs[-1], prec)
    for i in range(len(coeffs) - 2, -1, -1):
        result = hp(coeffs[i], prec) + hp(1, prec) / result

    return str(result)


# ======================================================================
# High-Precision Euler's Constant gamma
# ======================================================================

def euler_gamma(digits: int = 50) -> str:
    """
    Compute Euler's constant gamma using the accelerated series.

    gamma = sum_{k=1}^{inf} (1/k - ln(1+1/k))
    = sum_{k=1}^{N} (1/k - ln(1+1/k)) + 1/(2N) - 1/(12N^2) + 1/(120N^4) - ...

    Uses ln(1+1/k) = 2*atanh(1/(2k+1)) for fast convergence of each term.
    """
    prec = digits + 10
    hp = ScaledInt

    # Choose N based on desired digits
    N = 10 ** (digits // 2 + 3)
    N = min(N, 500000)  # cap for performance

    # Compute partial sum
    result = hp(0, prec)
    for k in range(1, N + 1):
        inv_k = hp(1, prec) / hp(k, prec)
        # ln(1 + 1/k) = 2 * atanh(1/(2k+1))
        denom = hp(2 * k + 1, prec)
        x = hp(1, prec) / denom
        x_sq = x * x
        ln_val = hp(0, prec)
        power = x
        for j in range(1, 80):
            term = power / hp(2 * j - 1, prec)
            ln_val = ln_val + term
            power = power * x_sq
            if power.integer == 0:
                break
        ln_val = hp(2, prec) * ln_val

        result = result + (inv_k - ln_val)

    # Euler-Maclaurin correction:
    # sum_{k=1}^{N} (1/k - ln(1+1/k)) = H_N - ln(N+1)
    # H_N - ln(N+1) = gamma - 1/(2N) + 1/(12N^2) - 1/(120N^4) + ...
    # So: gamma = sum + 1/(2N) - 1/(12N^2) + 1/(120N^4) - ...
    inv_N = hp(1, prec) / hp(N, prec)
    correction = hp(0, prec)
    correction = correction + inv_N / hp(2, prec)
    correction = correction - inv_N ** 2 / hp(12, prec)
    correction = correction + inv_N ** 4 / hp(120, prec)
    correction = correction - inv_N ** 6 / hp(252, prec)
    result = result + correction

    return str(result)


# ======================================================================
# High-Precision Zeta Function
# ======================================================================

def zeta_high_precision(s: int, digits: int = 50) -> str:
    """
    Compute zeta(s) for positive integer s.

    For even s, uses Euler's formula with pi.
    For odd s, uses the accelerated series.

    Args:
        s: positive integer
        digits: number of decimal digits

    Returns:
        String representation of zeta(s)
    """
    prec = digits + 10
    hp = ScaledInt

    if s == 2:
        pi_str = pi_machin(prec)
        pi_hp = hp(pi_str, prec)
        return str(pi_hp * pi_hp / hp(6, prec))

    if s == 4:
        pi_str = pi_machin(prec)
        pi_hp = hp(pi_str, prec)
        return str(pi_hp ** 4 / hp(90, prec))

    if s == 6:
        pi_str = pi_machin(prec)
        pi_hp = hp(pi_str, prec)
        return str(pi_hp ** 6 / hp(945, prec))

    # General case: accelerated series
    # zeta(s) = sum_{k=1}^{inf} 1/k^s
    n_terms = digits * 3 + 50
    result = hp(0, prec)

    for k in range(1, n_terms + 1):
        result = result + hp(1, prec) / hp(k, prec) ** s

    return str(result)


# ======================================================================
# Demonstration
# ======================================================================

def high_precision_demo():
    """Demonstrate high-precision computations."""
    print("=== High-Precision Computation Demo ===\n")

    # Test basic arithmetic
    print("--- ScaledInt Type ---")
    a = ScaledInt("1.23456789012345678901234567890", 30)
    b = ScaledInt("0.98765432109876543210987654321", 30)
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a + b = {a + b}")
    print(f"a * b = {a * b}")
    print(f"a - b = {a - b}")
    print(f"a / b = {a / b}")
    print()

    # Test sqrt
    print("--- Square Root ---")
    two = ScaledInt(2, 50)
    sqrt2 = two.sqrt()
    print(f"sqrt(2) = {sqrt2}")
    print()

    # Compute pi
    print("--- Pi Computation ---")
    print(f"pi (Machin, 50 digits): {pi_machin(50)}")
    print(f"pi (AGM, 50 digits): {pi_gauss_agm(50)}")
    print()

    # Compute e
    print("--- e Computation ---")
    print(f"e (series, 50 digits): {e_series(50)}")
    print(f"e (continued fraction, 50 digits): {e_continued_fraction(50)}")
    print()

    # Compute gamma
    print("--- Euler's Constant gamma ---")
    print(f"gamma (accelerated series, 30 digits): {euler_gamma(30)}")
    print()

    # Compute zeta values
    print("--- Zeta Function ---")
    print(f"zeta(2) = {zeta_high_precision(2, 50)}")
    print(f"zeta(4) = {zeta_high_precision(4, 50)}")
    print()


if __name__ == "__main__":
    high_precision_demo()
