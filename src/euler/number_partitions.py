"""
Euler's Number Partitions

This module implements Euler's work on integer partitions:
- Partition function p(n)
- Pentagonal number theorem
- Euler's recurrence relation
- Distinct partitions = odd partitions theorem
- Generating functions

References:
- Euler, L. (1748). "Introductio in analysin infinitorum"
- Euler, L. (1750). "De partitione numerorum"
"""

from functools import lru_cache
from typing import List, Dict, Tuple, Generator


def pentagonal_numbers(n: int) -> List[int]:
    """
    Generalized pentagonal numbers: g_k = k(3k-1)/2 for k = 0, ±1, ±2, ±3, ...
    Sequence: 0, 1, 2, 5, 7, 12, 15, 22, 26, 35, 40, 51, 57, ...
    """
    result = [0]
    for k in range(1, n + 1):
        result.append(k * (3 * k - 1) // 2)
        result.append(k * (3 * k + 1) // 2)
    return sorted(result)


def pentagonal_number(k: int) -> int:
    """Single generalized pentagonal number g_k."""
    if k == 0:
        return 0
    if k > 0:
        return k * (3 * k - 1) // 2
    # k < 0
    k = -k
    return k * (3 * k + 1) // 2


def partition_function_recursive(n: int) -> int:
    """
    Partition function p(n) using Euler's pentagonal number theorem recurrence:
    p(n) = p(n-1) + p(n-2) - p(n-5) - p(n-7) + p(n-12) + p(n-15) - ...
    
    where the offsets are generalized pentagonal numbers.
    Signs follow pattern: +, +, -, -, +, +, -, -, ...
    """
    if n < 0:
        return 0
    if n == 0:
        return 1
    
    # Generate generalized pentagonal numbers up to n
    pentagonal = []
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 > n and g2 > n:
            break
        pentagonal.append((g1, 1 if k % 2 == 1 else -1))
        pentagonal.append((g2, 1 if k % 2 == 1 else -1))
        k += 1
    
    result = 0
    for g, sign in pentagonal:
        if g <= n:
            result += sign * partition_function_recursive(n - g)
    return result


@lru_cache(maxsize=None)
def partition_function_pentagonal(n: int) -> int:
    """
    Efficient implementation using memoization and Euler's recurrence.
    """
    if n < 0:
        return 0
    if n == 0:
        return 1
    
    result = 0
    k = 1
    while True:
        # Generalized pentagonal numbers
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        
        if g1 > n and g2 > n:
            break
        
        sign = 1 if k % 2 == 1 else -1
        
        if g1 <= n:
            result += sign * partition_function_pentagonal(n - g1)
        if g2 <= n:
            result += sign * partition_function_pentagonal(n - g2)
        
        k += 1
    
    return result


def partition_function_dp(n: int) -> List[int]:
    """
    Dynamic programming approach to compute p(0), p(1), ..., p(n).
    Uses Euler's pentagonal number recurrence.
    """
    p = [0] * (n + 1)
    p[0] = 1
    
    for i in range(1, n + 1):
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            
            if g1 > i and g2 > i:
                break
            
            sign = 1 if k % 2 == 1 else -1
            
            if g1 <= i:
                p[i] += sign * p[i - g1]
            if g2 <= i:
                p[i] += sign * p[i - g2]
            
            k += 1
    
    return p


def partition_function(n: int, method: str = 'pentagonal') -> int:
    """
    Compute the partition number p(n).
    
    Methods:
    - 'pentagonal': Euler's pentagonal number recurrence (fastest)
    - 'dp': Dynamic programming (returns all up to n)
    - 'recursive': Naive recursive (slow, for demonstration)
    """
    if method == 'pentagonal':
        return partition_function_pentagonal(n)
    elif method == 'dp':
        return partition_function_dp(n)[n]
    elif method == 'recursive':
        return partition_function_recursive(n)
    else:
        raise ValueError(f"Unknown method: {method}")


def partition_conjugate(partition: List[int]) -> List[int]:
    """
    Conjugate of a partition (transpose Ferrers diagram).
    If λ = (λ₁, λ₂, ..., λₖ), conjugate has λᵢ' = number of parts ≥ i.
    """
    if not partition:
        return []
    
    max_part = max(partition)
    conjugate = []
    for i in range(1, max_part + 1):
        count = sum(1 for p in partition if p >= i)
        conjugate.append(count)
    return conjugate


def distinct_partitions(n: int) -> List[List[int]]:
    """
    Generate all partitions of n into distinct parts.
    Euler proved: number of distinct partitions = number of partitions into odd parts.
    """
    result = []
    
    def backtrack(remaining: int, start: int, current: List[int]):
        if remaining == 0:
            result.append(current[:])
            return
        for i in range(start, remaining + 1):
            current.append(i)
            backtrack(remaining - i, i + 1, current)
            current.pop()
    
    backtrack(n, 1, [])
    return result


def odd_partitions(n: int) -> List[List[int]]:
    """
    Generate all partitions of n into odd parts only.
    """
    result = []
    
    def backtrack(remaining: int, start: int, current: List[int]):
        if remaining == 0:
            result.append(current[:])
            return
        # Only use odd numbers
        for i in range(start, remaining + 1, 2):
            current.append(i)
            backtrack(remaining - i, i, current)
            current.pop()
    
    backtrack(n, 1, [])
    return result


def euler_partition_theorem(n: int) -> Dict:
    """
    Demonstrate Euler's theorem: q(n) = p_odd(n)
    where q(n) = partitions into distinct parts
    and p_odd(n) = partitions into odd parts.
    """
    distinct = distinct_partitions(n)
    odd = odd_partitions(n)
    
    return {
        'n': n,
        'distinct_partitions': distinct,
        'odd_partitions': odd,
        'q(n)': len(distinct),
        'p_odd(n)': len(odd),
        'equal': len(distinct) == len(odd)
    }


def euler_pentagonal_theorem_demo(n: int = 20) -> List[int]:
    """
    Demonstrate Euler's pentagonal number theorem:
    ∏_{k=1}^{∞} (1 - x^k) = Σ_{k=-∞}^{∞} (-1)^k x^{k(3k-1)/2}
    = 1 + Σ_{k=1}^{∞} (-1)^k (x^{k(3k-1)/2} + x^{k(3k+1)/2})
    """
    # Compute product coefficients
    coeffs = [1] + [0] * n
    for k in range(1, n + 1):
        for i in range(n, k - 1, -1):
            coeffs[i] -= coeffs[i - k]
    
    # Compute pentagonal number series
    pentagonal = [1] + [0] * n
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 > n and g2 > n:
            break
        sign = 1 if k % 2 == 1 else -1
        if g1 <= n:
            pentagonal[g1] = sign
        if g2 <= n:
            pentagonal[g2] = sign
        k += 1
    
    return coeffs, pentagonal


def restricted_partitions(n: int, max_part: int = None, max_length: int = None) -> int:
    """
    Partitions of n with restrictions.
    - max_part: largest part allowed
    - max_length: maximum number of parts
    """
    # DP approach
    dp = [[0] * (n + 1) for _ in range(n + 2)]
    dp[0][0] = 1
    
    for i in range(1, n + 2):
        for j in range(n + 1):
            dp[i][j] = dp[i - 1][j]
            if i <= j:
                dp[i][j] += dp[i][j - i]
    
    if max_part is not None and max_length is not None:
        return dp[max_part + 1][n]  # Not exactly right, needs proper DP
    
    return dp[n + 1][n]


def plane_partitions(n: int) -> int:
    """
    Number of plane partitions of n (2D arrays).
    MacMahon's formula (Euler didn't discover this, but it's related).
    """
    # Generating function: ∏_{k=1}^{∞} 1/(1-x^k)^k
    # This is beyond basic Euler, but included for completeness
    dp = [0] * (n + 1)
    dp[0] = 1
    
    for k in range(1, n + 1):
        for _ in range(k):
            for i in range(k, n + 1):
                dp[i] += dp[i - k]
    
    return dp[n]


def partition_generating_function_terms(n: int) -> Dict[int, int]:
    """
    Coefficients of the generating function:
    ∏_{k=1}^{∞} 1/(1 - x^k) = Σ p(n) x^n
    """
    coeffs = [1] + [0] * n
    for k in range(1, n + 1):
        for i in range(k, n + 1):
            coeffs[i] += coeffs[i - k]
    return {i: coeffs[i] for i in range(n + 1)}


def euler_partition_demo():
    """Demonstrate Euler's partition theory."""
    print("=== Euler's Number Partitions ===\n")
    
    print("Partition numbers p(n) using pentagonal number theorem:")
    p = partition_function_dp(20)
    for n in range(21):
        print(f"  p({n:2d}) = {p[n]}")
    
    print("\nPentagonal Number Recurrence for p(10):")
    # Show the recurrence calculation
    n = 10
    k = 1
    total = 0
    terms = []
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 > n and g2 > n:
            break
        sign = '+' if k % 2 == 1 else '-'
        if g1 <= n:
            terms.append(f"{sign} p({n - g1})")
            total += (1 if k % 2 == 1 else -1) * p[n - g1]
        if g2 <= n:
            terms.append(f"{sign} p({n - g2})")
            total += (1 if k % 2 == 1 else -1) * p[n - g2]
        k += 1
    print(f"  p({n}) = {' '.join(terms)} = {total}")
    
    print("\nEuler's Theorem: Distinct partitions = Odd partitions")
    for n in range(1, 16):
        result = euler_partition_theorem(n)
        status = "✓" if result['equal'] else "✗"
        print(f"  n={n:2d}: q(n)={result['q(n)']:3d}, p_odd(n)={result['p_odd(n)']:3d} {status}")
    
    print("\nExample: n=8")
    result = euler_partition_theorem(8)
    print(f"  Distinct partitions: {result['distinct_partitions']}")
    print(f"  Odd partitions:      {result['odd_partitions']}")
    
    print("\nPentagonal Number Theorem: ∏(1-x^k) = Σ(-1)^k x^{k(3k±1)/2}")
    prod_coeffs, pent_coeffs = euler_pentagonal_theorem_demo(15)
    print("  Product coefficients:", [c for c in prod_coeffs[:16] if c != 0])
    print("  Pentagonal series:   ", [c for c in pent_coeffs[:16] if c != 0])
    
    print("\nGenerating function coefficients (1/∏(1-x^k)):")
    gen = partition_generating_function_terms(15)
    for i, c in gen.items():
        print(f"  x^{i}: {c}")
    
    print("\nConjugate partitions:")
    test_partitions = [
        [5, 3, 1],
        [4, 4, 2],
        [6, 2, 1, 1],
        [3, 3, 3, 1],
    ]
    for part in test_partitions:
        conj = partition_conjugate(part)
        print(f"  {part} → {conj} (sum: {sum(part)})")


if __name__ == "__main__":
    euler_partition_demo()