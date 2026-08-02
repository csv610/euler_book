"""
Euler's Fluid Dynamics and Streamline Theory

This module implements Euler's contributions to fluid dynamics:
- Inviscid flow equations (Euler equations)
- Streamline computation and visualization data
- Stream function and velocity potential
- Bernoulli's principle along streamlines
- Potential flow around cylinders
- Vortex and source/sink flows

References:
- Euler, L. (1755). "Lezioni di statica dei fluidi"
- Euler, L. (1757). "Principes generaux du mouvement des fluides"
- D'Alembert, J. (1768). "Essai d'une nouvelle theorie de la resistance des fluides"
"""

from math import sin, cos, pi, sqrt, log, atan2, exp
from typing import List, Tuple, Optional

# Pure Python 2D array utilities (no numpy dependency)
def make_grid(nx: int, ny: int, val: float = 0.0) -> List[List[float]]:
    return [[val for _ in range(nx)] for _ in range(ny)]

def grid_get(grid: List[List[float]], i: int, j: int) -> float:
    return grid[i][j]

def grid_set(grid: List[List[float]], i: int, j: int, val: float):
    grid[i][j] = val


def euler_equations_steady(
    vx: callable, vy: callable,
    rho: float = 1.0, gx: float = 0.0, gy: float = 0.0
) -> Tuple[callable, callable]:
    """
    Compute pressure gradient from steady Euler equations.

    For steady inviscid flow with body forces:
        rho * (u * du/dx + v * du/dy) = -dp/dx + rho * gx
        rho * (u * dv/dx + v * dv/dy) = -dp/dy + rho * gy

    Args:
        vx: velocity component u(x,y) as callable
        vy: velocity component v(x,y) as callable
        rho: fluid density
        gx: body force in x-direction
        gy: body force in y-direction

    Returns:
        (dpdx, dpdy) pressure gradient functions
    """
    def dpdx(x, y, dx=1e-5):
        ux = vx(x, y)
        uy = vx(x + dx, y)
        vx_val = vy(x, y)
        vx_dy = vy(x + dx, y)
        du_dx = (uy - ux) / dx
        du_dy = (vy(x, y + dx) - vx(x, y)) / dx
        return -rho * (ux * du_dx + vx_val * du_dy) + rho * gx

    def dpdy(x, y, dx=1e-5):
        ux = vx(x, y)
        vx_val = vy(x, y)
        dv_dx = (vy(x + dx, y) - vx_val) / dx
        dv_dy = (vy(x, y + dx) - vy(x, y)) / dx
        return -rho * (ux * dv_dx + vx_val * dv_dy) + rho * gy

    return dpdx, dpdy


def streamline_integrate(
    vx: callable, vy: callable,
    x0: float, y0: float,
    ds: float = 0.01,
    n_steps: int = 1000,
    max_speed: float = 1e-10
) -> List[Tuple[float, float]]:
    """
    Integrate a streamline from (x0, y0) using RK4.

    A streamline satisfies dx/u = dy/v, so the tangent at each point
    is (u(x,y), v(x,y)).

    Args:
        vx: x-velocity component u(x,y)
        vy: y-velocity component v(x,y)
        x0: starting x coordinate
        y0: starting y coordinate
        ds: step size
        n_steps: maximum number of steps
        max_speed: minimum speed below which integration stops

    Returns:
        List of (x, y) points along the streamline
    """
    points = [(x0, y0)]
    x, y = x0, y0

    for _ in range(n_steps):
        u = vx(x, y)
        v = vy(x, y)
        speed = sqrt(u * u + v * v)

        if speed < max_speed:
            break

        # RK4 integration
        k1u, k1v = u, v

        k2u = vx(x + 0.5 * ds * k1u, y + 0.5 * ds * k1v)
        k2v = vy(x + 0.5 * ds * k1u, y + 0.5 * ds * k1v)

        k3u = vx(x + 0.5 * ds * k2u, y + 0.5 * ds * k2v)
        k3v = vy(x + 0.5 * ds * k2u, y + 0.5 * ds * k2v)

        k4u = vx(x + ds * k3u, y + ds * k3v)
        k4v = vy(x + ds * k3u, y + ds * k3v)

        x += ds * (k1u + 2 * k2u + 2 * k3u + k4u) / 6.0
        y += ds * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0

        points.append((x, y))

    return points


def uniform_flow(vx_val: float = 1.0, vy_val: float = 0.0):
    """Uniform flow with constant velocity."""
    def vx(x, y):
        return vx_val
    def vy(x, y):
        return vy_val
    return vx, vy


def source_flow(strength: float = 1.0):
    """Source flow at origin. Velocity decays as 1/r."""
    def vx(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return strength * x / r2
    def vy(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return strength * y / r2
    return vx, vy


def sink_flow(strength: float = 1.0):
    """Sink flow at origin (negative source)."""
    return source_flow(-strength)


def vortex_flow(circulation: float = 1.0):
    """Free vortex at origin. Tangential velocity ~ 1/r."""
    def vx(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return -circulation * y / (2 * pi * r2)
    def vy(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return circulation * x / (2 * pi * r2)
    return vx, vy


def doublet_flow(strength: float = 1.0):
    """Doublet flow (source-sink pair at origin)."""
    def vx(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return strength * (x * x - y * y) / (r2 * r2)
    def vy(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        return -strength * 2 * x * y / (r2 * r2)
    return vx, vy


def flow_past_cylinder(radius: float = 1.0, freestream: float = 1.0):
    """
    Potential flow around a circular cylinder (d'Alembert's paradox).

    Superposition of uniform flow + doublet.
    Velocity potential: phi = U*(r + a^2/r)*cos(theta)
    Stream function: psi = U*(r - a^2/r)*sin(theta)
    """
    def vx(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        r = sqrt(r2)
        if r < radius * 0.1:
            return 0.0
        ratio = radius * radius
        return freestream * (1.0 - ratio * (x*x - y*y) / (r2*r2))
    def vy(x, y):
        r2 = x * x + y * y
        if r2 < 1e-15:
            return 0.0
        r = sqrt(r2)
        if r < radius * 0.1:
            return 0.0
        ratio = radius * radius
        return -freestream * 2.0 * ratio * x * y / (r2*r2)
    return vx, vy


def stream_function(vx_func, vy_func, x_vals: List[float], y_vals: List[float]):
    """
    Compute stream function psi by integrating velocity field.

    psi(x,y) satisfies: dpsi/dy = vx, dpsi/dx = -vy
    Uses numerical integration from a reference point.
    """
    nx, ny = len(x_vals), len(y_vals)
    psi = make_grid(nx, ny, 0.0)

    # Integrate along x first
    for j in range(1, ny):
        for i in range(1, nx):
            dx = x_vals[i] - x_vals[i-1]
            psi[i][j] = psi[i-1][j] - vy_func(x_vals[i], y_vals[j]) * dx

    # Then correct along y
    psi_y = make_grid(nx, ny, 0.0)
    for j in range(1, ny):
        for i in range(nx):
            dy = y_vals[j] - y_vals[j-1]
            psi_y[i][j] = psi_y[i][j-1] + vx_func(x_vals[i], y_vals[j]) * dy

    # Average the two paths
    for i in range(nx):
        for j in range(ny):
            psi[i][j] = (psi[i][j] + psi_y[i][j]) / 2.0

    return psi


def velocity_potential(vx_func, vy_func, x_vals: List[float], y_vals: List[float]):
    """
    Compute velocity potential phi where vx = dphi/dx, vy = dphi/dy.
    """
    nx, ny = len(x_vals), len(y_vals)
    phi = make_grid(nx, ny, 0.0)

    # Integrate along x first
    for j in range(ny):
        for i in range(1, nx):
            dx = x_vals[i] - x_vals[i-1]
            phi[i][j] = phi[i-1][j] + vx_func(x_vals[i], y_vals[j]) * dx

    # Then correct along y
    phi_y = make_grid(nx, ny, 0.0)
    for j in range(1, ny):
        for i in range(nx):
            dy = y_vals[j] - y_vals[j-1]
            phi_y[i][j] = phi_y[i][j-1] + vy_func(x_vals[i], y_vals[j]) * dy

    phi = phi_y
    for i in range(nx):
        for j in range(ny):
            phi[i][j] = (phi[i][j] + phi_y[i][j]) / 2.0

    return phi


def bernoulli_pressure(vx_func, vy_func, rho: float = 1.0,
                       p_inf: float = 1.0, v_inf: float = 1.0):
    """
    Compute pressure field from Bernoulli's equation.

    For steady, inviscid, incompressible flow:
        p + 0.5*rho*|v|^2 = p_inf + 0.5*rho*v_inf^2 = constant

    Args:
        vx_func: x-velocity function
        vy_func: y-velocity function
        rho: fluid density
        p_inf: far-field pressure
        v_inf: far-field velocity magnitude

    Returns:
        Pressure function p(x,y)
    """
    const = p_inf + 0.5 * rho * v_inf * v_inf

    def p(x, y):
        u = vx_func(x, y)
        v = vy_func(x, y)
        speed_sq = u * u + v * v
        return const - 0.5 * rho * speed_sq

    return p


def stagnation_points(vx_func, vy_func, x_range=(-5, 5), y_range=(-5, 5),
                      nx=200, ny=200, tol=1e-6):
    """
    Find stagnation points where velocity = 0.

    Uses Newton's method from multiple initial guesses.
    """
    points = []
    for i in range(nx):
        for j in range(ny):
            x0 = x_range[0] + i * (x_range[1] - x_range[0]) / nx
            y0 = y_range[0] + j * (y_range[1] - y_range[0]) / ny

            u0 = vx_func(x0, y0)
            v0 = vy_func(x0, y0)
            if sqrt(u0*u0 + v0*v0) < tol * 10:
                # Newton's method
                x, y = x0, y0
                for _ in range(50):
                    u = vx_func(x, y)
                    v = vy_func(x, y)
                    if sqrt(u*u + v*v) < tol:
                        # Check not already found
                        is_new = True
                        for px, py in points:
                            if abs(x - px) < 0.1 and abs(y - py) < 0.1:
                                is_new = False
                                break
                        if is_new:
                            points.append((x, y))
                        break

                    # Jacobian
                    h = 1e-5
                    dudx = (vx_func(x+h, y) - u) / h
                    dudy = (vx_func(x, y+h) - u) / h
                    dvdx = (vy_func(x+h, y) - v) / h
                    dvdy = (vy_func(x, y+h) - v) / h

                    det = dudx * dvdy - dudy * dvdx
                    if abs(det) < 1e-15:
                        break

                    dx = -(u * dvdy - v * dudy) / det
                    dy = -(dudx * v - u * dvdx) / det
                    x += dx
                    y += dy
                    if abs(dx) < tol and abs(dy) < tol:
                        break

    return points


def compute_streamline_data(
    vx_func, vy_func,
    origins: List[Tuple[float, float]],
    x_range: Tuple[float, float] = (-5.0, 5.0),
    y_range: Tuple[float, float] = (-5.0, 5.0),
    n_points: int = 200,
    n_steps: int = 500,
    ds: float = 0.05
) -> dict:
    """
    Compute full streamline data for visualization.

    Args:
        vx_func: x-velocity function
        vy_func: y-velocity function
        origins: list of (x0, y0) starting points
        x_range: (xmin, xmax) for grid
        y_range: (ymin, ymax) for grid
        n_points: grid resolution
        n_steps: steps per streamline
        ds: step size

    Returns:
        Dictionary with 'streamlines', 'origins', 'x_range', 'y_range',
        'stagnation_points', and 'speed_field'
    """
    streamlines = []
    for ox, oy in origins:
        pts = streamline_integrate(vx_func, vy_func, ox, oy, ds=ds, n_steps=n_steps)
        if len(pts) > 2:
            streamlines.append(pts)

    # Compute speed field on grid
    xs = [x_range[0] + i * (x_range[1] - x_range[0]) / (n_points - 1) for i in range(n_points)]
    ys = [y_range[0] + j * (y_range[1] - y_range[0]) / (n_points - 1) for j in range(n_points)]
    speed = make_grid(n_points, n_points, 0.0)
    for i in range(n_points):
        for j in range(n_points):
            u = vx_func(xs[i], ys[j])
            v = vy_func(xs[i], ys[j])
            speed[i][j] = sqrt(u*u + v*v)

    stag_pts = stagnation_points(vx_func, vy_func, x_range, y_range, nx=100, ny=100)

    return {
        'streamlines': streamlines,
        'origins': origins,
        'x_range': x_range,
        'y_range': y_range,
        'stagnation_points': stag_pts,
        'speed_field': speed.tolist(),
        'xs': xs.tolist(),
        'ys': ys.tolist(),
    }
