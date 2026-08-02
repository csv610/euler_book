"""
Euler's Number Theory Functions

This module implements Euler's contributions to number theory including:
- Euler's totient function φ(n)
- Euler's theorem (generalization of Fermat's little theorem)
- Prime sieves and factorization
- Divisor functions and perfect numbers
- Amicable numbers

References:
- Euler, L. (1763). "Theoremata arithmetica nova methodo demonstrata"
- Euler, L. (1736). "De fractionibus continuis"
"""

from math import gcd, isqrt
from typing import List, Set, Tuple


def is_prime(n: int) -> bool:
    """Check if n is prime using trial division."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, isqrt(n) + 1, 2):
        if n % i == 0:
            return False
    return True


def prime_sieve(limit: int) -> List[int]:
    """Sieve of Eratosthenes to find all primes up to limit."""
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i, is_p in enumerate(sieve) if is_p]


def prime_factors(n: int) -> List[int]:
    """Return prime factorization of n as a list of prime factors."""
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors.append(n)
    return factors


def divisors(n: int) -> List[int]:
    """Return all positive divisors of n."""
    divs = set()
    for i in range(1, isqrt(n) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(divs)


def sigma_function(n: int, k: int = 1) -> int:
    """
    Sum of k-th powers of divisors of n.
    σ_k(n) = sum_{d|n} d^k
    """
    return sum(d ** k for d in divisors(n))


def euler_totient(n: int) -> int:
    """
    Euler's totient function φ(n): count of positive integers ≤ n that are coprime to n.
    
    Uses the formula: φ(n) = n * ∏_{p|n} (1 - 1/p)
    where the product is over distinct prime factors of n.
    """
    if n <= 0:
        return 0
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1 if p == 2 else 2
    if n > 1:
        result -= result // n
    return result


def euler_totient_sieve(limit: int) -> List[int]:
    """
    Compute φ(n) for all n from 1 to limit using a sieve.
    More efficient for computing many values.
    """
    phi = list(range(limit + 1))
    for i in range(2, limit + 1):
        if phi[i] == i:  # i is prime
            for j in range(i, limit + 1, i):
                phi[j] -= phi[j] // i
    return phi


def euler_totient_theorem(a: int, n: int) -> bool:
    """
    Euler's theorem: a^φ(n) ≡ 1 (mod n) for gcd(a, n) = 1.
    Returns True if the theorem holds for given a, n.
    """
    if gcd(a, n) != 1:
        return False
    return pow(a, euler_totient(n), n) == 1


def fermat_little_theorem(a: int, p: int) -> bool:
    """
    Fermat's little theorem: a^(p-1) ≡ 1 (mod p) for prime p and gcd(a, p) = 1.
    Special case of Euler's theorem when n = p (prime).
    """
    if not is_prime(p):
        return False
    if a % p == 0:
        return False
    return pow(a, p - 1, p) == 1


def perfect_numbers(limit: int) -> List[int]:
    """
    Find all perfect numbers up to limit.
    A perfect number equals the sum of its proper divisors.
    Even perfect numbers have form 2^(p-1) * (2^p - 1) where 2^p - 1 is Mersenne prime.
    """
    perfect = []
    # Check using Euclid-Euler theorem form
    for p in range(2, 20):  # reasonable limit for Mersenne primes
        mersenne = (1 << p) - 1
        if is_prime(mersenne):
            perfect_num = (1 << (p - 1)) * mersenne
            if perfect_num <= limit:
                perfect.append(perfect_num)
    return perfect


def amicable_pairs(limit: int) -> List[Tuple[int, int]]:
    """
    Find all amicable pairs (a, b) where a < b ≤ limit.
    Amicable numbers: sum of proper divisors of a equals b, and vice versa.
    Euler discovered many amicable pairs.
    """
    # Precompute sum of proper divisors
    sigma_proper = [0] * (limit + 1)
    for i in range(1, limit + 1):
        for j in range(2 * i, limit + 1, i):
            sigma_proper[j] += i
    
    pairs = []
    for a in range(1, limit + 1):
        b = sigma_proper[a]
        if a < b <= limit and sigma_proper[b] == a:
            pairs.append((a, b))
    return pairs


def euler_totient_ratio(n: int) -> float:
    """
    Compute n/φ(n) - the ratio that appears in many analytic number theory formulas.
    """
    if n <= 0:
        return 0.0
    return n / euler_totient(n)


def sum_of_totients(n: int) -> int:
    """
    Compute sum_{k=1}^n φ(k).
    Asymptotically: 3n²/π² + O(n log n)
    """
    phi = euler_totient_sieve(n)
    return sum(phi[1:])


# Example usage and demonstrations
if __name__ == "__main__":
    print("=== Euler's Number Theory Demonstrations ===\n")
    
    # Totient function
    print("Euler's Totient Function φ(n):")
    for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 30]:
        print(f"  φ({n}) = {euler_totient(n)}")
    
    print("\nEuler's Theorem (a^φ(n) ≡ 1 mod n for gcd(a,n)=1):")
    for a, n in [(2, 9), (3, 10), (5, 12), (7, 15)]:
        print(f"  {a}^{euler_totient(n)} ≡ {pow(a, euler_totient(n), n)} (mod {n}) - {euler_totient_theorem(a, n)}")
    
    print("\nFermat's Little Theorem (a^(p-1) ≡ 1 mod p):")
    for a, p in [(2, 7), (3, 11), (5, 13)]:
        print(f"  {a}^{p-1} ≡ {pow(a, p-1, p)} (mod {p}) - {fermat_little_theorem(a, p)}")
    
    print("\nPerfect numbers up to 10000:")
    for p in perfect_numbers(10000):
        print(f"  {p} = {prime_factors(p)}")
    
    print("\nAmicable pairs up to 10000:")
    for a, b in amicable_pairs(10000):
        print(f"  ({a}, {b})")
    
    print("\nSum of totients:")
    for n in [10, 100, 1000]:
        print(f"  Sum_1^{n} φ(k) = {sum_of_totients(n)}")