"""
Euler's Mathematics - A Python library for exploring Euler's mathematical discoveries.

This package provides pure Python implementations of Euler's major mathematical
contributions, including number theory, infinite series, infinite products,
graph theory, complex analysis, calculus, number partitions, polyhedron geometry,
and special functions.

All implementations are in pure Python with zero external dependencies.

References:
- Euler, L. (1735). "De summis serierum reciprocarum"
- Euler, L. (1736). "Solutio problematis ad geometriam situs pertinentis"
- Euler, L. (1748). "Introductio in analysin infinitorum"
- Euler, L. (1752). "Elementa doctrinae solidorum"
- Euler, L. (1768). "Institutiones calculi integralis"
"""

__version__ = "1.0.0"
__author__ = "Euler Student Guide Project"

from .number_theory import (
    euler_totient,
    euler_totient_sieve,
    euler_totient_theorem,
    fermat_little_theorem,
    is_prime,
    prime_sieve,
    prime_factors,
    divisors,
    sigma_function,
    perfect_numbers,
    amicable_pairs,
    sum_of_totients,
)

from .infinite_series import (
    harmonic_series,
    harmonic_series_partial,
    basel_problem,
    basel_partial_sum,
    euler_constant,
    euler_constant_partial,
    alternating_harmonic,
    leibniz_pi,
    zeta_function,
    zeta_even,
    euler_mascheroni_gamma,
)

from .infinite_products import (
    euler_product_formula,
    wallis_product,
    sinc_product,
    euler_product_sine,
    wallis_product_partial,
    sinc_product_partial,
    viete_formula,
)

from .graph_theory import (
    euler_characteristic,
    euler_characteristic_polyhedron,
    euler_formula_planar,
    konigsberg_bridges,
    eulerian_path_exists,
    eulerian_circuit_exists,
    is_eulerian_graph,
    is_semi_eulerian_graph,
    fleury_algorithm,
    hierholzer_algorithm,
)

from .complex_analysis import (
    euler_formula,
    euler_identity,
    complex_exponential,
    complex_logarithm,
    complex_power,
    de_moivre_formula,
    roots_of_unity,
    euler_formula_series,
)

from .calculus import (
    euler_method,
    euler_method_system,
    improved_euler_method,
    runge_kutta_4,
    euler_maclaurin_summation,
    euler_maclaurin_integral,
    euler_constant_limit,
)

from .number_partitions import (
    partition_function,
    partition_function_pentagonal,
    partition_function_recursive,
    partition_function_dp,
    pentagonal_numbers,
    euler_pentagonal_theorem_demo,
    partition_conjugate,
    distinct_partitions,
    odd_partitions,
    euler_partition_theorem,
)

from .polyhedron import (
    euler_polyhedron_formula,
    platonic_solids,
    verify_platonic_solids,
    euler_characteristic_surface,
    polyhedron_dual,
    euler_characteristic_genus,
    archimedean_solids,
    verify_archimedean_solids,
    platonic_solid_volume,
    platonic_solid_surface_area,
)

from .special_functions import (
    gamma_function,
    gamma_function_lanczos,
    beta_function,
    beta_function_integral,
    digamma_function,
    euler_reflection_formula,
    euler_beta_integral_demo,
    legendre_duplication,
)

from .fluid_dynamics import (
    streamline_integrate,
    uniform_flow,
    source_flow,
    sink_flow,
    vortex_flow,
    doublet_flow,
    flow_past_cylinder,
    stream_function,
    velocity_potential,
    bernoulli_pressure,
    stagnation_points,
    compute_streamline_data,
    euler_equations_steady,
)

from .pde_solvers import (
    heat_equation_explicit,
    heat_equation_crank_nicolson,
    wave_equation_explicit,
    laplace_equation_jacobi,
    poisson_equation_gauss_seidel,
    von_neumann_stability,
    solve_heat_1d,
    solve_wave_1d,
)

from .euler_numbers import (
    euler_number_zigzag,
    euler_number_classical,
    tangent_taylor,
    secant_taylor,
    alternating_permutations_count,
    entringer_triangle,
    euler_secant_numbers,
)

from .high_precision import (
    ScaledInt,
    pi_machin,
    pi_gauss_agm,
    e_series,
    e_continued_fraction,
    euler_gamma,
    zeta_high_precision,
)

__all__ = [
    # number_theory
    "euler_totient",
    "euler_totient_sieve",
    "euler_totient_theorem",
    "fermat_little_theorem",
    "is_prime",
    "prime_sieve",
    "prime_factors",
    "divisors",
    "sigma_function",
    "perfect_numbers",
    "amicable_pairs",
    "sum_of_totients",
    # infinite_series
    "harmonic_series",
    "harmonic_series_partial",
    "basel_problem",
    "basel_partial_sum",
    "euler_constant",
    "euler_constant_partial",
    "alternating_harmonic",
    "leibniz_pi",
    "zeta_function",
    "zeta_even",
    "euler_mascheroni_gamma",
    # infinite_products
    "euler_product_formula",
    "wallis_product",
    "sinc_product",
    "euler_product_sine",
    "wallis_product_partial",
    "sinc_product_partial",
    "viete_formula",
    # graph_theory
    "euler_characteristic",
    "euler_characteristic_polyhedron",
    "euler_formula_planar",
    "konigsberg_bridges",
    "eulerian_path_exists",
    "eulerian_circuit_exists",
    "is_eulerian_graph",
    "is_semi_eulerian_graph",
    "fleury_algorithm",
    "hierholzer_algorithm",
    # complex_analysis
    "euler_formula",
    "euler_identity",
    "complex_exponential",
    "complex_logarithm",
    "complex_power",
    "de_moivre_formula",
    "roots_of_unity",
    "euler_formula_series",
    # calculus
    "euler_method",
    "euler_method_system",
    "improved_euler_method",
    "runge_kutta_4",
    "euler_maclaurin_summation",
    "euler_maclaurin_integral",
    "euler_constant_limit",
    # number_partitions
    "partition_function",
    "partition_function_pentagonal",
    "partition_function_recursive",
    "partition_function_dp",
    "pentagonal_numbers",
    "euler_pentagonal_theorem_demo",
    "partition_conjugate",
    "distinct_partitions",
    "odd_partitions",
    "euler_partition_theorem",
    # polyhedron
    "euler_polyhedron_formula",
    "platonic_solids",
    "verify_platonic_solids",
    "euler_characteristic_surface",
    "polyhedron_dual",
    "euler_characteristic_genus",
    "archimedean_solids",
    "verify_archimedean_solids",
    "platonic_solid_volume",
    "platonic_solid_surface_area",
    # special_functions
    "gamma_function",
    "gamma_function_lanczos",
    "beta_function",
    "beta_function_integral",
    "digamma_function",
    "euler_reflection_formula",
    "euler_beta_integral_demo",
    "legendre_duplication",
    # fluid_dynamics
    "streamline_integrate",
    "uniform_flow",
    "source_flow",
    "sink_flow",
    "vortex_flow",
    "doublet_flow",
    "flow_past_cylinder",
    "stream_function",
    "velocity_potential",
    "bernoulli_pressure",
    "stagnation_points",
    "compute_streamline_data",
    "euler_equations_steady",
    # pde_solvers
    "heat_equation_explicit",
    "heat_equation_crank_nicolson",
    "wave_equation_explicit",
    "laplace_equation_jacobi",
    "poisson_equation_gauss_seidel",
    "von_neumann_stability",
    "solve_heat_1d",
    "solve_wave_1d",
    # euler_numbers
    "euler_number_zigzag",
    "euler_number_classical",
    "tangent_taylor",
    "secant_taylor",
    "alternating_permutations_count",
    "entringer_triangle",
    "euler_secant_numbers",
    # high_precision
    "pi_machin",
    "pi_gauss_agm",
    "e_series",
    "e_continued_fraction",
    "euler_gamma",
    "zeta_high_precision",
    "ScaledInt",
]