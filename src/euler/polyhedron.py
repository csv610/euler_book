"""
Euler's Polyhedron Formula and Geometry

This module implements Euler's work on polyhedra and topology:
- Euler's formula for convex polyhedra: V - E + F = 2
- Platonic solids
- Dual polyhedra
- Euler characteristic for surfaces

References:
- Euler, L. (1752). "Elementa doctrinae solidorum"
- Euler, L. (1758). "Demonstratio nonnullarum insignium proprietatum"
"""

from math import sqrt
from typing import Dict, List, Tuple, Optional


def euler_polyhedron_formula(V: int, E: int, F: int) -> int:
    """
    Euler's polyhedron formula: V - E + F = 2
    For any convex polyhedron.
    """
    return V - E + F


def platonic_solids() -> Dict[str, Dict]:
    """
    Return all five Platonic solids with their properties.
    Each has congruent regular polygon faces, same number of faces at each vertex.
    """
    return {
        'tetrahedron': {
            'name': 'Tetrahedron',
            'faces': 4,
            'edges': 6,
            'vertices': 4,
            'face_type': 'triangle',
            'faces_per_vertex': 3,
            'schläfli': '{3,3}',
            'dual': 'tetrahedron',  # self-dual
            'symmetry_group': 'Td (order 24)',
            'dihedral_angle': 'arccos(1/3) ≈ 70.53°',
        },
        'cube': {
            'name': 'Cube (Hexahedron)',
            'faces': 6,
            'edges': 12,
            'vertices': 8,
            'face_type': 'square',
            'faces_per_vertex': 3,
            'schläfli': '{4,3}',
            'dual': 'octahedron',
            'symmetry_group': 'Oh (order 48)',
            'dihedral_angle': '90°',
        },
        'octahedron': {
            'name': 'Octahedron',
            'faces': 8,
            'edges': 12,
            'vertices': 6,
            'face_type': 'triangle',
            'faces_per_vertex': 4,
            'schläfli': '{3,4}',
            'dual': 'cube',
            'symmetry_group': 'Oh (order 48)',
            'dihedral_angle': 'arccos(-1/3) ≈ 109.47°',
        },
        'dodecahedron': {
            'name': 'Dodecahedron',
            'faces': 12,
            'edges': 30,
            'vertices': 20,
            'face_type': 'pentagon',
            'faces_per_vertex': 3,
            'schläfli': '{5,3}',
            'dual': 'icosahedron',
            'symmetry_group': 'Ih (order 120)',
            'dihedral_angle': 'arccos(-1/√5) ≈ 116.57°',
        },
        'icosahedron': {
            'name': 'Icosahedron',
            'faces': 20,
            'edges': 30,
            'vertices': 12,
            'face_type': 'triangle',
            'faces_per_vertex': 5,
            'schläfli': '{3,5}',
            'dual': 'dodecahedron',
            'symmetry_group': 'Ih (order 120)',
            'dihedral_angle': 'arccos(-√5/3) ≈ 138.19°',
        },
    }


def verify_platonic_solids() -> Dict[str, Dict]:
    """Verify Euler's formula V - E + F = 2 for all Platonic solids."""
    solids = platonic_solids()
    results = {}
    for key, solid in solids.items():
        V = solid['vertices']
        E = solid['edges']
        F = solid['faces']
        chi = V - E + F
        results[key] = {
            **solid,
            'euler_characteristic': chi,
            'satisfies_euler': chi == 2,
        }
    return results


def euler_characteristic_surface(V: int, E: int, F: int, genus: int = 0) -> int:
    """
    Euler characteristic for surfaces: χ = V - E + F = 2 - 2g
    where g is the genus (number of holes).
    For sphere (g=0): χ = 2
    For torus (g=1): χ = 0
    For double torus (g=2): χ = -2
    """
    return V - E + F


def euler_characteristic_genus(V: int, E: int, F: int) -> int:
    """Compute genus from Euler characteristic: g = 1 - χ/2 = 1 - (V-E+F)/2"""
    chi = V - E + F
    return (2 - chi) // 2


def polyhedron_dual(V: int, E: int, F: int) -> Tuple[int, int, int]:
    """
    Dual polyhedron swaps vertices and faces:
    V* = F, F* = V, E* = E
    """
    return (F, E, V)


def regular_polygon_area(n_sides: int, side_length: float = 1.0) -> float:
    """Area of a regular n-gon with given side length."""
    from math import tan, pi
    return n_sides * side_length * side_length / (4 * tan(pi / n_sides))


def platonic_solid_volume(name: str, edge_length: float = 1.0) -> float:
    """Volume of Platonic solid with given edge length."""
    from math import sqrt
    volumes = {
        'tetrahedron': sqrt(2) / 12,
        'cube': 1.0,
        'octahedron': sqrt(2) / 3,
        'dodecahedron': (15 + 7 * sqrt(5)) / 4,
        'icosahedron': 5 * (3 + sqrt(5)) / 12,
    }
    return volumes.get(name.lower(), 0) * edge_length ** 3


def platonic_solid_surface_area(name: str, edge_length: float = 1.0) -> float:
    """Surface area of Platonic solid with given edge length."""
    from math import sqrt
    areas = {
        'tetrahedron': sqrt(3),
        'cube': 6.0,
        'octahedron': 2 * sqrt(3),
        'dodecahedron': 3 * sqrt(25 + 10 * sqrt(5)),
        'icosahedron': 5 * sqrt(3),
    }
    return areas.get(name.lower(), 0) * edge_length ** 2


def euler_characteristic_polyhedron(V: int, E: int, F: int) -> int:
    """Euler's formula for convex polyhedra."""
    return V - E + F


def euler_formula_planar(V: int, E: int, F: int) -> int:
    """Euler's formula for connected planar graphs: V - E + F = 2 (including outer face)."""
    return V - E + F


def is_planar_graph(V: int, E: int) -> bool:
    """Check if a graph could be planar using E ≤ 3V - 6 (for V ≥ 3)."""
    if V < 3:
        return True
    return E <= 3 * V - 6


def archimedean_solids() -> Dict[str, Dict]:
    """
    The 13 Archimedean solids (semi-regular convex polyhedra).
    All vertices identical, faces are regular polygons (not all same).
    """
    return {
        'truncated_tetrahedron': {
            'name': 'Truncated Tetrahedron',
            'faces': 8, 'edges': 18, 'vertices': 12,
            'face_types': ['triangle', 'hexagon'],
            'vertex_config': '3.6.6',
        },
        'cuboctahedron': {
            'name': 'Cuboctahedron',
            'faces': 14, 'edges': 24, 'vertices': 12,
            'face_types': ['triangle', 'square'],
            'vertex_config': '3.4.3.4',
        },
        'truncated_cube': {
            'name': 'Truncated Cube',
            'faces': 14, 'edges': 36, 'vertices': 24,
            'face_types': ['triangle', 'octagon'],
            'vertex_config': '3.8.8',
        },
        'truncated_octahedron': {
            'name': 'Truncated Octahedron',
            'faces': 14, 'edges': 36, 'vertices': 24,
            'face_types': ['square', 'hexagon'],
            'vertex_config': '4.6.6',
        },
        'rhombicuboctahedron': {
            'name': 'Rhombicuboctahedron',
            'faces': 26, 'edges': 48, 'vertices': 24,
            'face_types': ['triangle', 'square'],
            'vertex_config': '3.4.4.4',
        },
        'truncated_cuboctahedron': {
            'name': 'Truncated Cuboctahedron',
            'faces': 26, 'edges': 72, 'vertices': 48,
            'face_types': ['square', 'hexagon', 'octagon'],
            'vertex_config': '4.6.8',
        },
        'snub_cube': {
            'name': 'Snub Cube',
            'faces': 38, 'edges': 60, 'vertices': 24,
            'face_types': ['triangle', 'square'],
            'vertex_config': '3.3.3.3.4',
        },
        'icosidodecahedron': {
            'name': 'Icosidodecahedron',
            'faces': 32, 'edges': 60, 'vertices': 30,
            'face_types': ['triangle', 'pentagon'],
            'vertex_config': '3.5.3.5',
        },
        'truncated_dodecahedron': {
            'name': 'Truncated Dodecahedron',
            'faces': 32, 'edges': 90, 'vertices': 60,
            'face_types': ['triangle', 'decagon'],
            'vertex_config': '3.10.10',
        },
        'truncated_icosahedron': {
            'name': 'Truncated Icosahedron (Buckyball)',
            'faces': 32, 'edges': 90, 'vertices': 60,
            'face_types': ['pentagon', 'hexagon'],
            'vertex_config': '5.6.6',
        },
        'rhombicosidodecahedron': {
            'name': 'Rhombicosidodecahedron',
            'faces': 62, 'edges': 120, 'vertices': 60,
            'face_types': ['triangle', 'square', 'pentagon'],
            'vertex_config': '3.4.5.4',
        },
        'truncated_icosidodecahedron': {
            'name': 'Truncated Icosidodecahedron',
            'faces': 62, 'edges': 180, 'vertices': 120,
            'face_types': ['square', 'hexagon', 'decagon'],
            'vertex_config': '4.6.10',
        },
        'snub_dodecahedron': {
            'name': 'Snub Dodecahedron',
            'faces': 92, 'edges': 150, 'vertices': 60,
            'face_types': ['triangle', 'pentagon'],
            'vertex_config': '3.3.3.3.5',
        },
    }


def verify_archimedean_solids() -> Dict[str, Dict]:
    """Verify Euler's formula for Archimedean solids."""
    solids = archimedean_solids()
    results = {}
    for key, solid in solids.items():
        V = solid['vertices']
        E = solid['edges']
        F = solid['faces']
        chi = V - E + F
        results[key] = {
            **solid,
            'euler_characteristic': chi,
            'satisfies_euler': chi == 2,
        }
    return results


def dual_polyhedron_info(name: str) -> Dict:
    """Information about the dual of a Platonic solid."""
    solids = platonic_solids()
    if name.lower() not in solids:
        return {'error': 'Not a Platonic solid'}
    
    solid = solids[name.lower()]
    dual_name = solid['dual']
    dual = solids.get(dual_name.lower(), {})
    
    return {
        'original': solid['name'],
        'dual': dual.get('name', dual_name),
        'vertex_face_swap': f"V={solid['vertices']} ↔ F={dual.get('faces', '?')}",
        'edges_preserved': f"E={solid['edges']} (same)",
    }


def euler_characteristic_demo():
    """Demonstrate Euler characteristic for various surfaces."""
    print("Euler Characteristic χ = V - E + F = 2 - 2g")
    print("=" * 50)
    
    # Convex polyhedra (sphere topology, g=0)
    print("\nConvex Polyhedra (genus 0):")
    for key, solid in platonic_solids().items():
        V, E, F = solid['vertices'], solid['edges'], solid['faces']
        print(f"  {solid['name']:20s}: V={V:2d}, E={E:2d}, F={F:2d}, χ={V-E+F} (g={euler_characteristic_genus(V, E, F)})")
    
    for key, solid in archimedean_solids().items():
        V, E, F = solid['vertices'], solid['edges'], solid['faces']
        print(f"  {solid['name']:30s}: V={V:3d}, E={E:3d}, F={F:2d}, χ={V-E+F} (g={euler_characteristic_genus(V, E, F)})")
    
    # Torus (genus 1)
    print("\nTorus (genus 1): χ = 0")
    # Minimal triangulation of torus: V=7, E=21, F=14
    print(f"  Minimal triangulation: V=7, E=21, F=14, χ={7-21+14} (g={euler_characteristic_genus(7, 21, 14)})")
    
    # Double torus (genus 2)
    print("\nDouble Torus (genus 2): χ = -2")
    print(f"  Example: V=10, E=30, F=20, χ={10-30+20} (g={euler_characteristic_genus(10, 30, 20)})")
    
    # General formula
    print("\nGeneral formula: χ = 2 - 2g")
    for g in range(6):
        print(f"  Genus {g}: χ = {2 - 2*g}")


if __name__ == "__main__":
    print("=== Euler's Polyhedron Formula ===\n")
    
    # Verify Platonic solids
    print("Platonic Solids (V - E + F = 2):")
    results = verify_platonic_solids()
    for name, data in results.items():
        print(f"  {data['name']:20s}: V={data['vertices']:2d}, E={data['edges']:2d}, F={data['faces']:2d}, χ={data['euler_characteristic']} ✓")
    
    print("\nArchimedean Solids:")
    results = verify_archimedean_solids()
    for name, data in results.items():
        print(f"  {data['name']:30s}: V={data['vertices']:3d}, E={data['edges']:3d}, F={data['faces']:2d}, χ={data['euler_characteristic']} ✓")
    
    print("\nDual Pairs:")
    for name in ['tetrahedron', 'cube', 'octahedron', 'dodecahedron', 'icosahedron']:
        info = dual_polyhedron_info(name)
        print(f"  {info['original']} ↔ {info['dual']}: {info['vertex_face_swap']}, {info['edges_preserved']}")
    
    print("\nVolumes and Surface Areas (edge length = 1):")
    for name in ['tetrahedron', 'cube', 'octahedron', 'dodecahedron', 'icosahedron']:
        vol = platonic_solid_volume(name)
        area = platonic_solid_surface_area(name)
        print(f"  {name:15s}: Volume = {vol:.6f}, Surface Area = {area:.6f}")
    
    print("\n" + "=" * 50)
    euler_characteristic_demo()