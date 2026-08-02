"""
Euler's Numerical PDE Solvers

This module implements numerical methods for solving partial differential
equations, following Euler's pioneering work in calculus of variations
and the calculus of differences.

- Heat equation (diffusion)
- Wave equation
- Laplace/Poisson equation
- Finite difference methods (explicit, implicit, Crank-Nicolson)
- Stability analysis (von Neumann)

References:
- Euler, L. (1755). "Lezioni di statica dei fluidi"
- Euler, L. (1759). "Recherches sur la nature et la propagation du bruit"
- Fourier, J. (1822). "Theorie analytique de la chaleur"
"""

from math import sin, cos, pi, exp, sqrt
from typing import Callable, List, Tuple, Optional


def heat_equation_explicit(
    u0: List[List[float]],
    alpha: float,
    dx: float,
    dt: float,
    nx: int,
    ny: int,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
    bc_top: float = 0.0,
    bc_bottom: float = 0.0,
    n_steps: int = 100,
) -> List[List[List[float]]]:
    """
    Solve 2D heat equation: du/dt = alpha * (d2u/dx2 + d2u/dy2)
    Using explicit finite differences (FTCS).

    Stability requires: alpha * dt / dx^2 <= 0.25 (for 2D).

    Args:
        u0: initial condition (ny x nx array)
        alpha: thermal diffusivity
        dx: spatial step size
        dt: time step size
        nx: number of x points
        ny: number of y points
        bc_left: boundary value at x=0
        bc_right: boundary value at x=L
        bc_top: boundary value at y=H
        bc_bottom: boundary value at y=0
        n_steps: number of time steps

    Returns:
        List of solution arrays at each time step
    """
    u = [list(row) for row in u0]  # deep copy
    solutions = [list(list(row) for row in u)]  # single initial snapshot

    r = alpha * dt / (dx * dx)

    # Stability check
    if r > 0.25:
        raise ValueError(
            f"Explicit scheme unstable: r={r:.4f} > 0.25. "
            f"Reduce dt or increase dx. Try dt <= {0.25 * dx * dx / alpha}"
        )

    for _ in range(n_steps):
        u_new = [list(row) for row in u]

        # Apply boundary conditions
        u_new[0] = [bc_bottom] * nx
        u_new[-1] = [bc_top] * nx
        for i in range(ny):
            u_new[i][0] = bc_left
            u_new[i][-1] = bc_right

        # Interior points: FTCS stencil
        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                u_new[i][j] = (u[i][j]
                               + r * (u[i+1][j] - 2*u[i][j] + u[i-1][j])
                               + r * (u[i][j+1] - 2*u[i][j] + u[i][j-1]))

        u = u_new
        solutions.append([list(row) for row in u])

    return solutions


def heat_equation_crank_nicolson(
    u0: List[List[float]],
    alpha: float,
    dx: float,
    dt: float,
    nx: int,
    ny: int,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
    bc_top: float = 0.0,
    bc_bottom: float = 0.0,
    n_steps: int = 100,
) -> List[List[List[float]]]:
    """
    Solve 2D heat equation using Crank-Nicolson method.
    Unconditionally stable, second-order accurate in time and space.

    Uses ADI (Alternating Direction Implicit) splitting for efficiency.
    """
    u = [list(row) for row in u0]
    solutions = [list(list(row) for row in u)]
    r = alpha * dt / (dx * dx)

    for _ in range(n_steps):
        # Sweep x-direction (implicit)
        for i in range(1, ny - 1):
            n = nx - 2
            if n <= 0:
                continue

            rhs = [r * u[i][j-1] + (1 - 2*r) * u[i][j] + r * u[i][j+1] for j in range(1, nx-1)]

            # Thomas algorithm
            c = [-r / (1 + 2*r)] * (n - 1)
            d = [rhs[0] / (1 + 2*r)]
            for k in range(1, n):
                d.append((rhs[k] + r * d[k-1]) / (1 + 2*r - r * c[k-1]))

            x_solve = [0.0] * n
            x_solve[-1] = d[-1]
            for k in range(n-2, -1, -1):
                x_solve[k] = d[k] + c[k] * x_solve[k+1]

            for j in range(1, nx-1):
                u[i][j] = x_solve[j-1]

        # Sweep y-direction (implicit)
        for j in range(1, nx - 1):
            n = ny - 2
            if n <= 0:
                continue

            rhs = [r * u[i-1][j] + (1 - 2*r) * u[i][j] + r * u[i+1][j] for i in range(1, ny-1)]

            c = [-r / (1 + 2*r)] * (n - 1)
            d = [rhs[0] / (1 + 2*r)]
            for k in range(1, n):
                d.append((rhs[k] + r * d[k-1]) / (1 + 2*r - r * c[k-1]))

            x_solve = [0.0] * n
            x_solve[-1] = d[-1]
            for k in range(n-2, -1, -1):
                x_solve[k] = d[k] + c[k] * x_solve[k+1]

            for i in range(1, ny-1):
                u[i][j] = x_solve[i-1]

        # Apply boundary conditions
        u[0] = [bc_bottom] * nx
        u[-1] = [bc_top] * nx
        for i in range(ny):
            u[i][0] = bc_left
            u[i][-1] = bc_right

        solutions.append([list(row) for row in u])

    return solutions


def wave_equation_explicit(
    u0: List[List[float]],
    v0: List[List[float]],
    c: float,
    dx: float,
    dt: float,
    nx: int,
    ny: int,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
    bc_top: float = 0.0,
    bc_bottom: float = 0.0,
    n_steps: int = 100,
) -> List[List[List[float]]]:
    """
    Solve 2D wave equation: d2u/dt2 = c^2 * (d2u/dx2 + d2u/dy2)
    Using central differences in time and space.

    CFL condition: c * dt / dx <= 1/sqrt(2) for 2D.

    Args:
        u0: initial displacement
        v0: initial velocity
        c: wave speed
        dx: spatial step
        dt: time step
        nx, ny: grid dimensions
        bc_*: boundary conditions
        n_steps: number of time steps

    Returns:
        List of solution arrays
    """
    u = [list(row) for row in u0]
    u_prev = [list(row) for row in u0]  # u at t = -dt (from v0)
    u_next = [[0.0]*nx for _ in range(ny)]

    co = c * dt / dx
    if co > 1.0 / sqrt(2):
        raise ValueError(
            f"Wave equation unstable: Courant number {co:.4f} > "
            f"{1.0/sqrt(2):.4f}. Reduce dt or dx."
        )

    # Initialize u_prev using v0
    for i in range(ny):
        for j in range(nx):
            u_prev[i][j] = u[i][j] - dt * v0[i][j]

    solutions = [list(list(row) for row in u)]

    for _ in range(n_steps):
        # Apply boundary conditions
        u_next[0] = [bc_bottom] * nx
        u_next[-1] = [bc_top] * nx
        for i in range(ny):
            u_next[i][0] = bc_left
            u_next[i][-1] = bc_right

        # Interior: central difference stencil
        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                u_next[i][j] = (2 * u[i][j] - u_prev[i][j]
                                + co * co * (
                                    u[i+1][j] - 2*u[i][j] + u[i-1][j]
                                    + u[i][j+1] - 2*u[i][j] + u[i][j-1]
                                ))

        u_prev = [list(row) for row in u]
        u = [list(row) for row in u_next]
        solutions.append(list(list(row) for row in u))

    return solutions


def laplace_equation_jacobi(
    u0: List[List[float]],
    nx: int,
    ny: int,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
    bc_top: float = 0.0,
    bc_bottom: float = 0.0,
    max_iter: int = 10000,
    tol: float = 1e-8,
) -> Tuple[List[List[float]], int]:
    """
    Solve Laplace equation: d2u/dx2 + d2u/dy2 = 0
    Using Jacobi iteration.

    Args:
        u0: initial guess (ny x nx array)
        bc_*: boundary values
        max_iter: maximum iterations
        tol: convergence tolerance

    Returns:
        (solution, num_iterations)
    """
    u = [list(row) for row in u0]
    u_new = [[0.0]*nx for _ in range(ny)]
    converged = False

    for it in range(max_iter):
        # Boundary conditions
        u_new[0] = [bc_bottom] * nx
        u_new[-1] = [bc_top] * nx
        for i in range(ny):
            u_new[i][0] = bc_left
            u_new[i][-1] = bc_right

        # Interior: average of neighbors
        max_diff = 0.0
        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                u_new[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1])
                diff = abs(u_new[i][j] - u[i][j])
                if diff > max_diff:
                    max_diff = diff

        u = [list(row) for row in u_new]

        if max_diff < tol:
            converged = True
            return u, it + 1

    return u, max_iter


def poisson_equation_gauss_seidel(
    f_source: List[List[float]],
    nx: int,
    ny: int,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
    bc_top: float = 0.0,
    bc_bottom: float = 0.0,
    max_iter: int = 10000,
    tol: float = 1e-8,
) -> Tuple[List[List[float]], int]:
    """
    Solve Poisson equation: d2u/dx2 + d2u/dy2 = f(x,y)
    Using Gauss-Seidel iteration with successive over-relaxation (SOR).

    Args:
        f_source: source term f(x,y) on grid
        bc_*: boundary values
        max_iter: maximum iterations
        tol: convergence tolerance

    Returns:
        (solution, num_iterations)
    """
    omega = 1.9  # relaxation parameter (optimal ~2 for 2D)
    u = [[0.0]*nx for _ in range(ny)]

    for it in range(max_iter):
        max_diff = 0.0

        # Boundary conditions
        u[0] = [bc_bottom] * nx
        u[-1] = [bc_top] * nx
        for i in range(ny):
            u[i][0] = bc_left
            u[i][-1] = bc_right

        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                u_old = u[i][j]
                u[i][j] = (1 - omega) * u[i][j] + \
                          omega * 0.25 * (u[i+1][j] + u[i-1][j] +
                                          u[i][j+1] + u[i][j-1] -
                                          f_source[i][j])
                diff = abs(u[i][j] - u_old)
                if diff > max_diff:
                    max_diff = diff

        if max_diff < tol:
            return u, it + 1

    return u, max_iter


def von_neumann_stability(
    scheme: str,
    alpha: float = 1.0,
    c: float = 1.0,
    dx: float = 0.1,
    dt: float = 0.01,
) -> dict:
    """
    Von Neumann stability analysis for common PDE schemes.

    Args:
        scheme: 'heat_explicit', 'heat_cn', 'wave_explicit', 'advection'
        alpha: diffusivity (for heat equation)
        c: wave/advection speed
        dx: spatial step
        dt: time step

    Returns:
        Dictionary with stability info
    """
    info = {
        'scheme': scheme,
        'dx': dx,
        'dt': dt,
        'stable': True,
        'constraints': [],
    }

    if scheme == 'heat_explicit':
        r = alpha * dt / (dx * dx)
        info['r'] = r
        if r > 0.25:
            info['stable'] = False
            info['constraints'].append(
                f"alpha*dt/dx^2 = {r:.4f} must be <= 0.25 for 2D stability"
            )
        else:
            info['constraints'].append(
                f"alpha*dt/dx^2 = {r:.4f} <= 0.25 (stable)"
            )

    elif scheme == 'heat_cn':
        info['stable'] = True  # Crank-Nicolson is unconditionally stable
        info['constraints'].append("Crank-Nicolson is unconditionally stable")

    elif scheme == 'wave_explicit':
        co = c * dt / dx
        info['courant'] = co
        if co > 1.0 / sqrt(2):
            info['stable'] = False
            info['constraints'].append(
                f"c*dt/dx = {co:.4f} must be <= {1.0/sqrt(2):.4f} for 2D stability"
            )
        else:
            info['constraints'].append(
                f"c*dt/dx = {co:.4f} <= {1.0/sqrt(2):.4f} (stable)"
            )

    elif scheme == 'advection':
        co = c * dt / dx
        info['courant'] = co
        if co > 1.0:
            info['stable'] = False
            info['constraints'].append(
                f"c*dt/dx = {co:.4f} must be <= 1.0 (CFL condition)"
            )
        else:
            info['constraints'].append(
                f"c*dt/dx = {co:.4f} <= 1.0 (stable)"
            )

    return info


def solve_heat_1d(
    u0: Callable[[float], float],
    alpha: float,
    L: float,
    T: float,
    nx: int = 100,
    nt: int = 1000,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Solve 1D heat equation: du/dt = alpha * d2u/dx2

    Args:
        u0: initial condition function u0(x)
        alpha: thermal diffusivity
        L: domain length
        T: final time
        nx: number of spatial points
        nt: number of time steps
        bc_left, bc_right: boundary values

    Returns:
        (x, t, u_final) arrays
    """
    dx = L / (nx - 1)
    dt = T / nt
    r = alpha * dt / (dx * dx)

    if r > 0.5:
        raise ValueError(f"1D scheme unstable: r={r:.4f} > 0.5")

    x = [i * dx for i in range(nx)]
    u = [u0(xi) for xi in x]

    for _ in range(nt):
        u_new = list(u)
        u_new[0] = bc_left
        u_new[-1] = bc_right
        for i in range(1, nx - 1):
            u_new[i] = u[i] + r * (u[i+1] - 2*u[i] + u[i-1])
        u = u_new

    t = [j * dt for j in range(nt)]
    return x, t, u


def solve_wave_1d(
    u0: Callable[[float], float],
    v0: Callable[[float], float],
    c: float,
    L: float,
    T: float,
    nx: int = 100,
    nt: int = 1000,
    bc_left: float = 0.0,
    bc_right: float = 0.0,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Solve 1D wave equation: d2u/dt2 = c^2 * d2u/dx2

    Args:
        u0: initial displacement
        v0: initial velocity
        c: wave speed
        L: domain length
        T: final time
        nx: spatial points
        nt: time steps
        bc_left, bc_right: boundary values

    Returns:
        (x, t, u_final) arrays
    """
    dx = L / (nx - 1)
    dt = T / nt
    co = c * dt / dx

    if co > 1.0:
        raise ValueError(f"1D wave scheme unstable: co={co:.4f} > 1.0")

    x = [i * dx for i in range(nx)]
    u = [u0(xi) for xi in x]
    u_prev = [u[i] - dt * v0(xi) for i, xi in enumerate(x)]

    for _ in range(nt):
        u_new = [0.0] * nx
        u_new[0] = bc_left
        u_new[-1] = bc_right
        for i in range(1, nx - 1):
            u_new[i] = (2*u[i] - u_prev[i]
                        + co*co * (u[i+1] - 2*u[i] + u[i-1]))
        u_prev = list(u)
        u = u_new

    t = [j * dt for j in range(nt)]
    return x, t, u
