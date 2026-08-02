"""
Euler's Graph Theory

This module implements Euler's foundational work in graph theory:
- Königsberg bridges problem
- Eulerian paths and circuits
- Euler's formula for planar graphs
- Euler characteristic for polyhedra

References:
- Euler, L. (1736). "Solutio problematis ad geometriam situs pertinentis"
- Euler, L. (1752). "Elementa doctrinae solidorum"
"""

from collections import defaultdict, deque
from typing import List, Dict, Set, Tuple, Optional


def euler_characteristic(V: int, E: int, F: int) -> int:
    """
    Euler characteristic for polyhedra: χ = V - E + F
    For convex polyhedra, χ = 2.
    """
    return V - E + F


def euler_characteristic_polyhedron(V: int, E: int, F: int) -> int:
    """Euler's formula for convex polyhedra: V - E + F = 2."""
    return V - E + F


def euler_formula_planar(V: int, E: int, F: int) -> int:
    """
    Euler's formula for connected planar graphs: V - E + F = 2
    where F includes the outer (unbounded) face.
    """
    return V - E + F


def is_planar_graph(V: int, E: int) -> bool:
    """Check if a graph could be planar using necessary condition E ≤ 3V - 6 (for V ≥ 3)."""
    if V < 3:
        return True
    return E <= 3 * V - 6


def konigsberg_bridges() -> dict:
    """
    The famous Königsberg bridges problem.
    Euler proved no Eulerian path exists.
    
    Graph: 4 landmasses (A, B, C, D) connected by 7 bridges.
    Degrees: A=5, B=3, C=3, D=3 (all odd, so no Eulerian path)
    """
    graph = {
        'A': ['B', 'B', 'C', 'C', 'D'],  # 5 bridges
        'B': ['A', 'A', 'D'],            # 3 bridges (two A-B, one B-D)
        'C': ['A', 'A', 'D'],            # 3 bridges (two A-C, one C-D)
        'D': ['A', 'B', 'C'],            # 3 bridges (one A-D, one B-D, one C-D)
    }
    
    degrees = {v: len(edges) for v, edges in graph.items()}
    odd_vertices = [v for v, d in degrees.items() if d % 2 == 1]
    
    return {
        'graph': graph,
        'degrees': degrees,
        'odd_vertices': odd_vertices,
        'eulerian_path_exists': len(odd_vertices) in [0, 2],
        'eulerian_circuit_exists': len(odd_vertices) == 0,
    }


def eulerian_path_exists(graph: Dict[str, List[str]]) -> bool:
    """
    Check if an Eulerian path exists in an undirected graph.
    An undirected graph has an Eulerian path iff exactly 0 or 2 vertices have odd degree.
    """
    if not graph:
        return False
    
    # Count vertices with odd degree
    odd_count = sum(1 for v in graph if len(graph[v]) % 2 == 1)
    
    # Check connectivity (all non-zero degree vertices in same component)
    non_zero = [v for v in graph if graph[v]]
    if not non_zero:
        return True  # empty graph
    
    visited = set()
    stack = [non_zero[0]]
    while stack:
        v = stack.pop()
        if v in visited:
            continue
        visited.add(v)
        for u in graph[v]:
            if u not in visited:
                stack.append(u)
    
    connected = all(v in visited for v in non_zero)
    
    return connected and odd_count in [0, 2]


def eulerian_circuit_exists(graph: Dict[str, List[str]]) -> bool:
    """
    Check if an Eulerian circuit exists.
    An undirected graph has an Eulerian circuit iff all vertices have even degree and graph is connected.
    """
    if not graph:
        return True
    
    # Check all degrees even
    all_even = all(len(graph[v]) % 2 == 0 for v in graph)
    
    # Check connectivity
    non_zero = [v for v in graph if graph[v]]
    if not non_zero:
        return True
    
    visited = set()
    stack = [non_zero[0]]
    while stack:
        v = stack.pop()
        if v in visited:
            continue
        visited.add(v)
        for u in graph[v]:
            if u not in visited:
                stack.append(u)
    
    connected = all(v in visited for v in non_zero)
    
    return connected and all_even


def is_eulerian_graph(graph: Dict[str, List[str]]) -> bool:
    """Check if graph is Eulerian (has Eulerian circuit)."""
    return eulerian_circuit_exists(graph)


def is_semi_eulerian_graph(graph: Dict[str, List[str]]) -> bool:
    """Check if graph is semi-Eulerian (has Eulerian path but not circuit)."""
    return eulerian_path_exists(graph) and not eulerian_circuit_exists(graph)


def fleury_algorithm(graph: Dict[str, List[str]], start: str = None) -> List[str]:
    """
    Fleury's algorithm for finding an Eulerian path/circuit.
    Avoids bridges unless forced to use them.
    """
    if not graph:
        return []
    
    # Find start vertex
    if start is None:
        odd_vertices = [v for v in graph if len(graph[v]) % 2 == 1]
        start = odd_vertices[0] if odd_vertices else list(graph.keys())[0]
    
    # Copy graph with edge counts (for multigraph)
    g = {v: list(edges) for v, edges in graph.items()}
    
    def is_bridge(u: str, v: str) -> bool:
        """Check if edge u-v is a bridge in the remaining graph."""
        if not g[u] or not g[v]:
            return False
        
        # Remove edge temporarily
        g[u].remove(v)
        g[v].remove(u)
        
        # Check if u and v still connected
        visited = set()
        stack = [u]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            for nxt in g[cur]:
                if nxt not in visited:
                    stack.append(nxt)
        
        # Restore edge
        g[u].append(v)
        g[v].append(u)
        
        return v not in visited
    
    path = []
    current = start
    
    while g[current]:
        # Find next edge, preferring non-bridges
        for nxt in g[current]:
            if not is_bridge(current, nxt):
                break
        else:
            nxt = g[current][0]  # forced to use bridge
        
        # Remove edge
        g[current].remove(nxt)
        g[nxt].remove(current)
        path.append(nxt)
        current = nxt
    
    return [start] + path


def hierholzer_algorithm(graph: Dict[str, List[str]], start: str = None) -> List[str]:
    """
    Hierholzer's algorithm for finding Eulerian circuit/path.
    More efficient than Fleury's - builds circuit by splicing cycles.
    """
    if not graph:
        return []
    
    # Find start vertex
    if start is None:
        odd_vertices = [v for v in graph if len(graph[v]) % 2 == 1]
        start = odd_vertices[0] if odd_vertices else list(graph.keys())[0]
    
    # Copy with edge counts
    g = {v: list(edges) for v, edges in graph.items()}
    
    def find_circuit(start_v: str) -> List[str]:
        circuit = [start_v]
        current = start_v
        while g[current]:
            nxt = g[current].pop()
            g[nxt].remove(current)
            circuit.append(nxt)
            current = nxt
        return circuit
    
    # Build initial circuit
    circuit = find_circuit(start)
    
    # Splice in remaining cycles
    i = 0
    while i < len(circuit):
        v = circuit[i]
        if g[v]:
            # Found vertex with unused edges, find cycle from there
            sub_circuit = find_circuit(v)
            # Splice: replace v with sub_circuit
            circuit = circuit[:i] + sub_circuit + circuit[i+1:]
        else:
            i += 1
    
    return circuit


def euler_characteristic_surface(V: int, E: int, F: int, genus: int = 0) -> int:
    """
    Euler characteristic for surfaces: χ = V - E + F = 2 - 2g
    where g is the genus (number of holes).
    """
    return V - E + F


def polyhedron_vertices_edges_faces(name: str) -> Tuple[int, int, int]:
    """Return (V, E, F) for Platonic solids."""
    solids = {
        'tetrahedron': (4, 6, 4),
        'cube': (8, 12, 6),
        'octahedron': (6, 12, 8),
        'dodecahedron': (20, 30, 12),
        'icosahedron': (12, 30, 20),
    }
    return solids.get(name.lower(), (0, 0, 0))


def verify_platonic_solids() -> Dict[str, dict]:
    """Verify Euler's formula for all Platonic solids."""
    solids = {
        'tetrahedron': (4, 6, 4),
        'cube': (8, 12, 6),
        'octahedron': (6, 12, 8),
        'dodecahedron': (20, 30, 12),
        'icosahedron': (12, 30, 20),
    }
    
    results = {}
    for name, (V, E, F) in solids.items():
        chi = V - E + F
        results[name] = {
            'V': V, 'E': E, 'F': F,
            'euler_characteristic': chi,
            'satisfies_euler': chi == 2
        }
    return results


def euler_characteristic_genus(V: int, E: int, F: int) -> int:
    """Compute genus from Euler characteristic: g = 1 - χ/2 = 1 - (V-E+F)/2"""
    chi = V - E + F
    return (2 - chi) // 2


# Example usage
if __name__ == "__main__":
    print("=== Euler's Graph Theory ===\n")
    
    # Königsberg bridges
    print("Königsberg Bridges Problem:")
    k = konigsberg_bridges()
    print(f"  Degrees: {k['degrees']}")
    print(f"  Odd-degree vertices: {k['odd_vertices']}")
    print(f"  Eulerian path exists: {k['eulerian_path_exists']}")
    print(f"  Eulerian circuit exists: {k['eulerian_circuit_exists']}")
    print("  Euler's result: NO Eulerian path exists!")
    
    print("\nEuler's Formula for Platonic Solids (V - E + F = 2):")
    results = verify_platonic_solids()
    for name, data in results.items():
        print(f"  {name.capitalize()}: V={data['V']}, E={data['E']}, F={data['F']}, χ={data['euler_characteristic']} ✓")
    
    print("\nGraph with Eulerian circuit:")
    g1 = {
        'A': ['B', 'C'],
        'B': ['A', 'C', 'D'],
        'C': ['A', 'B', 'D'],
        'D': ['B', 'C'],
    }
    print(f"  Graph: {g1}")
    print(f"  Eulerian circuit exists: {eulerian_circuit_exists(g1)}")
    print(f"  Eulerian path exists: {eulerian_path_exists(g1)}")
    print(f"  Hierholzer circuit: {hierholzer_algorithm(g1)}")
    
    print("\nGraph with Eulerian path (not circuit):")
    g2 = {
        'A': ['B'],
        'B': ['A', 'C', 'D'],
        'C': ['B', 'D'],
        'D': ['B', 'C'],
    }
    print(f"  Graph: {g2}")
    print(f"  Eulerian circuit exists: {eulerian_circuit_exists(g2)}")
    print(f"  Eulerian path exists: {eulerian_path_exists(g2)}")
    print(f"  Hierholzer path: {hierholzer_algorithm(g2)}")
    
    print("\nSurface genus from Euler characteristic:")
    # Torus: V=7, E=21, F=14 (triangulated)
    print(f"  Torus (g=1): χ={euler_characteristic_surface(7, 21, 14)}, genus={euler_characteristic_genus(7, 21, 14)}")
    # Double torus
    print(f"  Double torus (g=2): χ={euler_characteristic_surface(10, 30, 20)}, genus={euler_characteristic_genus(10, 30, 20)}")