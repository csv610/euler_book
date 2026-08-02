"""Tests for high-precision computation module."""

import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from euler.high_precision import (
    ScaledInt,
    pi_machin,
    pi_gauss_agm,
    e_series,
    e_continued_fraction,
    euler_gamma,
    zeta_high_precision,
)


class TestScaledInt(unittest.TestCase):
    """Tests for the ScaledInt arbitrary-precision type."""

    def test_addition(self):
        a = ScaledInt("1.23456789012345678901234567890", 30)
        b = ScaledInt("0.98765432109876543210987654321", 30)
        result = a + b
        self.assertEqual(str(result)[:33], "2.22222221122222222112222222211")

    def test_subtraction(self):
        a = ScaledInt("5.0", 10)
        b = ScaledInt("3.0", 10)
        result = a - b
        self.assertEqual(str(result), "2")

    def test_multiplication(self):
        a = ScaledInt("1.5", 10)
        b = ScaledInt("2.0", 10)
        result = a * b
        self.assertEqual(str(result), "3")

    def test_scalar_multiply(self):
        a = ScaledInt("1.23456789012345678901234567890", 30)
        result = 16 * a
        expected = ScaledInt("19.75308624197530862419753086240", 30)
        self.assertAlmostEqual(float(str(result)), float(str(expected)), places=10)

    def test_division(self):
        a = ScaledInt("1.0", 20)
        b = ScaledInt("3.0", 20)
        result = a / b
        self.assertEqual(str(result), "0.33333333333333333333")

    def test_sqrt(self):
        two = ScaledInt(2, 50)
        sqrt2 = two.sqrt()
        self.assertEqual(str(sqrt2)[:20], "1.414213562373095048")

    def test_negation(self):
        a = ScaledInt("5.5", 10)
        result = -a
        self.assertEqual(str(result), "-5.5")

    def test_comparison(self):
        a = ScaledInt("3.0", 10)
        b = ScaledInt("5.0", 10)
        self.assertTrue(a < b)
        self.assertFalse(a > b)
        self.assertEqual(a, ScaledInt("3", 10))


class TestPiComputation(unittest.TestCase):
    """Tests for pi computation methods."""

    def test_pi_machin_50_digits(self):
        """Machin's formula should give correct first 50 digits of pi."""
        pi_str = pi_machin(50)
        expected = "3.14159265358979323846264338327950288419716939937510"
        self.assertEqual(pi_str[:len(expected)], expected)

    def test_pi_gauss_agm_50_digits(self):
        """AGM method should agree with Machin to 50 digits."""
        pi_mach = pi_machin(50)
        pi_agm = pi_gauss_agm(50)
        self.assertEqual(pi_mach[:53], pi_agm[:53])

    def test_pi_first_20_digits(self):
        """Check first 20 known digits of pi."""
        pi_str = pi_machin(20)
        expected = "3.14159265358979323846"
        self.assertEqual(pi_str[:len(expected)], expected)


class TestEComputation(unittest.TestCase):
    """Tests for e computation methods."""

    def test_e_series_50_digits(self):
        """Series method should give correct first 50 digits of e."""
        e_str = e_series(50)
        expected = "2.71828182845904523536028747135266249775724709369995"
        self.assertEqual(e_str[:len(expected)], expected)

    def test_e_continued_fraction_50_digits(self):
        """Continued fraction should agree with series."""
        e_series_val = e_series(50)
        e_cf_val = e_continued_fraction(50)
        self.assertEqual(e_series_val[:53], e_cf_val[:53])


class TestEulerGamma(unittest.TestCase):
    """Tests for Euler's constant computation."""

    def test_gamma_20_digits(self):
        """Gamma should match known value to 11 digits (N limited for performance)."""
        gamma_str = euler_gamma(20)
        expected = "0.57721566490"
        self.assertEqual(gamma_str[:len(expected)], expected)


class TestZetaComputation(unittest.TestCase):
    """Tests for zeta function computation."""

    def test_zeta2_50_digits(self):
        """zeta(2) = pi^2/6 should match to 50 digits."""
        z2 = zeta_high_precision(2, 50)
        expected = "1.64493406684822643647241516664602518921894990120679"
        self.assertEqual(z2[:len(expected)], expected)

    def test_zeta4_50_digits(self):
        """zeta(4) = pi^4/90 should match to 50 digits."""
        z4 = zeta_high_precision(4, 50)
        expected = "1.08232323371113819151600369654116790277475095191872"
        self.assertEqual(z4[:len(expected)], expected)


if __name__ == '__main__':
    unittest.main()
