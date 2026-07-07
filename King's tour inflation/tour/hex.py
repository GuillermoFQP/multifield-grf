#!/usr/bin/env python3
"""
Random DFS + backtracking Hamiltonian path on a HEXAGONAL board (hex grid),
with a no-self-intersection constraint on the drawn trajectory.

Board = hexagon of radius R in axial coordinates (q, r):
    max(|q|, |r|, |q+r|) <= R
Number of cells = 1 + 3R(R+1)

Moves = 6 adjacent hex neighbors.

Trajectory = straight segments between hex centers; we disallow any segment
intersection with earlier segments (except sharing the most recent endpoint,
which we skip by construction).
"""

from __future__ import annotations
import json, random
from sys import argv
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from scipy.interpolate import BSpline

Hex = Tuple[int, int]      # axial (q, r)
Pt = Tuple[float, float]   # (x, y) in plot coordinates


# 6 axial neighbor directions
HEX_DIRS = [(1, 0), (1, -1), (0, -1),
            (-1, 0), (-1, 1), (0, 1)]


def hex_cells(radius: int) -> List[Hex]:
    """All axial coords in a hexagon of given radius."""
    cells = []
    for q in range(-radius, radius + 1):
        for r in range(-radius, radius + 1):
            s = -q - r
            if max(abs(q), abs(r), abs(s)) <= radius:
                cells.append((q, r))
    return cells


def neighbors(h: Hex, cellset: Set[Hex]) -> List[Hex]:
    q, r = h
    out = []
    for dq, dr in HEX_DIRS:
        hh = (q + dq, r + dr)
        if hh in cellset:
            out.append(hh)
    return out


def axial_to_xy(h: Hex, size: float = 1.0) -> Pt:
    """
    Pointy-top axial -> 2D.
    x = size * sqrt(3) * (q + r/2)
    y = size * 3/2 * r
    """
    q, r = h
    x = size * np.sqrt(3.0) * (q + 0.5 * r)
    y = size * 1.5 * r
    return (float(x), float(y))


# --- Geometry (segment intersection) ---

def orient(a: Pt, b: Pt, c: Pt) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a: Pt, b: Pt, c: Pt) -> bool:
    return (min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and
            min(a[1], b[1]) <= c[1] <= max(a[1], b[1]))


def segments_intersect(p1: Pt, p2: Pt, q1: Pt, q2: Pt) -> bool:
    """
    Returns True if segments p1p2 and q1q2 intersect in any way
    (including collinear touch/overlap).
    """
    o1 = orient(p1, p2, q1)
    o2 = orient(p1, p2, q2)
    o3 = orient(q1, q2, p1)
    o4 = orient(q1, q2, p2)

    # Collinear/touch cases
    if (o1 == 0 and on_segment(p1, p2, q1)) or (o2 == 0 and on_segment(p1, p2, q2)) \
       or (o3 == 0 and on_segment(q1, q2, p1)) or (o4 == 0 and on_segment(q1, q2, p2)):
        return True

    # Proper intersection
    return (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0)


def would_self_intersect(path: List[Hex], new_hex: Hex, pos2xy: Dict[Hex, Pt]) -> bool:
    """
    Check whether adding segment from path[-1] -> new_hex intersects any earlier segment.
    Skip the most recent existing segment (it shares endpoint path[-1]).
    """
    if len(path) < 2:
        return False

    p_last = pos2xy[path[-1]]
    p_new = pos2xy[new_hex]

    # existing segments: path[i] -> path[i+1], i=0..len(path)-2
    # skip i=len(path)-2 (the last segment)
    for i in range(0, len(path) - 2):
        a = pos2xy[path[i]]
        b = pos2xy[path[i + 1]]
        if segments_intersect(a, b, p_last, p_new):
            return True
    return False


# --- Solver ---

class HexTourDFSNoCross:
    def __init__(self, cells: List[Hex], pos2xy: Dict[Hex, Pt],
                 seed: Optional[int] = None, warnsdorff_noise: float = 0.35):
        if seed is not None:
            random.seed(seed)

        self.cells = cells
        self.cellset = set(cells)
        self.pos2xy = pos2xy
        self.N = len(cells)
        self.warnsdorff_noise = warnsdorff_noise

        self.adj = {h: neighbors(h, self.cellset) for h in self.cells}
        self.calls = 0

    def onward_count(self, h: Hex, visited: Set[Hex]) -> int:
        return sum((n not in visited) for n in self.adj[h])

    def ordered_moves(self, h: Hex, visited: Set[Hex]) -> List[Hex]:
        cand = []
        for n in self.adj[h]:
            if n in visited:
                continue
            score = self.onward_count(n, visited) + random.random() * self.warnsdorff_noise
            cand.append((score, n))
        cand.sort(key=lambda t: t[0])

        # extra randomness among the top few
        k = min(3, len(cand))
        head = cand[:k]
        random.shuffle(head)
        cand = head + cand[k:]
        return [n for _, n in cand]

    def dfs(self, h: Hex, visited: Set[Hex], path: List[Hex]) -> bool:
        self.calls += 1
        if len(path) == self.N:
            return True

        for n in self.ordered_moves(h, visited):
            if would_self_intersect(path, n, self.pos2xy):
                continue
            visited.add(n)
            path.append(n)
            if self.dfs(n, visited, path):
                return True
            path.pop()
            visited.remove(n)
        return False

    def solve(self, start: Optional[Hex] = None, max_restarts: int = 500) -> Optional[List[Hex]]:
        for _ in range(max_restarts):
            self.calls = 0
            s = start if start is not None else random.choice(self.cells)
            visited = {s}
            path = [s]
            if self.dfs(s, visited, path):
                return path
        return None


# --- Plotting ---

def plot_hex_tour(cells: List[Hex], pos2xy: Dict[Hex, Pt], path: List[Hex],
                  hex_size: float = 1.0, annotate_every: int = 4) -> None:
    fig, ax = plt.subplots(figsize=(5,5), frameon=False)

    # draw hex tiles (light checker-ish by (q+r) parity)
    for h in cells:
        x, y = pos2xy[h]
        parity = (h[0] + h[1]) & 1
        face = "white" if parity == 0 else "#e8e8e8"
        hex_patch = RegularPolygon(
            (x, y),
            numVertices=6,
            radius=hex_size,
            orientation=np.radians(30),   # pointy-top
            facecolor=face,
            edgecolor="black",
            linewidth=0.8,
        )
        ax.add_patch(hex_patch)

    # trajectory
    xs = [pos2xy[h][0] for h in path]
    ys = [pos2xy[h][1] for h in path]
    ax.plot(xs, ys, linewidth=20)
    #ax.scatter(xs[0], ys[0], s=90, label="start")
    #ax.scatter(xs[-1], ys[-1], s=90, marker="s", label="end")

    # annotate steps
    for i, h in enumerate(path):
        if i % annotate_every == 0 or i in (0, len(path) - 1):
            x, y = pos2xy[h]
            ax.text(x, y, str(i), ha="center", va="center", fontsize=10)

    '''
    # B-spline representation
    degree = 3; pts = len(path)
    knots = np.zeros(pts+degree+1); knots[-degree:] = 1.0
    knots[degree:-degree] = np.linspace(0.0, 1.0, pts-degree+1)
    spline_x = BSpline(knots, xs, degree)
    spline_y = BSpline(knots, ys, degree)

    t = np.linspace(0.0, 1.0, 1024)
    plt.plot(spline_x(t), spline_y(t), '-', color='tab:orange', linewidth=5)
    '''
    
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    #ax.set_title(f"Hex board tour (cells={len(cells)}) with no self-intersections")
    #ax.legend()
    plt.tight_layout()
    plt.show()


def main():
    # Hexagon radius: R=1 -> 7 cells, R=2 -> 19, R=3 -> 37, R=4 -> 61 ...
    R = 2

    cells = hex_cells(R)
    pos2xy = {h: axial_to_xy(h, size=1.0) for h in cells}

    idx = None if len(argv) < 2 else int(argv[1])

    if idx is not None:
        with open('hex.json', 'r') as catalog:
            path = [(a-b,b) for a,b in json.load(catalog)[idx]]
    else:
        solver = HexTourDFSNoCross(cells, pos2xy, seed=None, warnsdorff_noise=0.35)
        path = solver.solve(start=None, max_restarts=1000)
        print([[a+b,b] for a,b in path])

    if path is None:
        raise RuntimeError(
            "No tour found under constraints. Try smaller R, more restarts, or increase randomness."
        )

    #print(f"Found path length {len(path)} in {solver.calls} DFS calls (last attempt).")
    plot_hex_tour(cells, pos2xy, path, hex_size=1.0, annotate_every=1)


if __name__ == "__main__":
    main()
