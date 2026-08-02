# Euler's Mathematics

A comprehensive reference book and companion Python library exploring Leonhard Euler's foundational contributions to mathematics.

## Overview

This project contains:

- **Euler's Mathematics** (53-page reference book) — a scholarly text covering Euler's work in analysis, number theory, graph theory, and more
- **Euler Student Guide** (235-page) — extended edition with exercises, solutions, and computational explorations
- **Python library** (`src/euler/`) — 14 modules implementing Euler's algorithms and formulas
- **98 tests** — comprehensive test suite verifying correctness

## Book Contents

| Chapter | Topic | Lines |
|---------|-------|-------|
| 1 | Euler's Totient Function | 491 |
| 2 | Prime Numbers | 665 |
| 3 | Perfect and Amicable Numbers | 369 |
| 4 | Harmonic Series | 409 |
| 5 | The Basel Problem | 356 |
| 6 | Infinite Products | 464 |
| 7 | The Königsberg Bridges | 405 |
| 8 | Polyhedra and Euler's Formula | 413 |
| 9 | Euler's Formula (Reference) | 378 |
| 10 | The Gamma Function | 479 |
| 12 | Euler's Formula (Student) | 390 |
| 18 | Streamlines | 587 |
| 19 | Partial Differential Equations | 384 |
| 20 | Euler Numbers | 473 |
| 21 | High-Precision Computation | 555 |

**Total**: 6,818 lines across 15 chapters (avg 455 lines/chapter)

## Quick Start

### Build PDFs

```bash
make pdf
```

### Run Tests

```bash
make test
```

### Install Package

```bash
make install
# or
pip install -e .
```

## Project Structure

```
Euler/
├── Euler_Book.tex                    # Main reference book
├── Euler_Student_Guide.tex           # Student guide with exercises
├── Makefile                          # Build automation
├── pyproject.toml                    # Python package config
├── .gitignore                        # Excludes build artifacts
├── appendix/
│   ├── appendix_constants.tex
│   ├── appendix_glossary.tex
│   └── appendix_solutions.tex
├── chapters/                         # 15 chapter source files
│   ├── ch01_totient.tex
│   ├── ch02_primes.tex
│   └── ...
├── figures/                          # Diagrams and illustrations
├── src/euler/                        # Python library
│   ├── __init__.py
│   ├── calculus.py
│   ├── complex_analysis.py
│   ├── euler_numbers.py
│   ├── fluid_dynamics.py
│   ├── graph_theory.py
│   ├── high_precision.py
│   ├── infinite_products.py
│   ├── infinite_series.py
│   ├── number_partitions.py
│   ├── number_theory.py
│   ├── pde_solvers.py
│   ├── polyhedron.py
│   └── special_functions.py
└── tests/
    ├── test_euler.py
    └── test_high_precision.py
```

## Python Library

The `euler` package implements algorithms from the book:

```python
from euler import (
    euler_totient,           # Euler's totient function
    primorial,               # Product of primes up to n
    euler_acceleration,      # Series acceleration
    catalan_constant,        # Catalan's G
    euler_mascheroni,        # Euler-Mascheroni γ
)

from euler.high_precision import (
    ScaledInt,               # Arbitrary precision integers
    compute_pi_gauss_agm,    # AGM algorithm for π
    compute_e_continued,     # Continued fraction for e
)
```

## Quality Metrics

| Metric | Value |
|--------|-------|
| Chapter average | 455 lines |
| Minimum chapter | 356 lines |
| Tests passing | 98/98 |
| Euler_Book PDF | 53 pages |
| Student Guide PDF | 235 pages |
| LaTeX errors | 0 |

## Requirements

- Python 3.9+
- TeX Live (for PDF generation)
- pytest (for tests)

No external Python dependencies required — the library is pure Python.

## Building

```bash
# Build both PDFs
make pdf

# Run all tests
make test

# Build PDF and run tests
make all

# Clean build artifacts
make clean
```

## License

MIT License

## References

- Euler, L. (1707-1783) — Original works and correspondence
- The book references standard mathematical texts and primary sources
