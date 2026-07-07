#!/usr/bin/env python3
"""
Random DFS + backtracking to find a king-move Hamiltonian path on a 4x4 board
with the extra constraint: the polyline trajectory must NOT self-intersect.

Intersections are checked between straight segments connecting cell centers.
"""

import json, random
from sys import argv
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

Pos = Tuple[int, int]
Pt = Tuple[float, float]

KING_DIRS = [(-1, -1), (-1, 0), (-1, 1),
             (0, -1),          (0, 1),
             (1, -1),  (1, 0), (1, 1)]


def in_bounds(r: int, c: int, n: int) -> bool:
    return 0 <= r < n and 0 <= c < n


def king_neighbors(p: Pos, n: int) -> List[Pos]:
    r, c = p
    out = []
    for dr, dc in KING_DIRS:
        rr, cc = r + dr, c + dc
        if in_bounds(rr, cc, n):
            out.append((rr, cc))
    return out


def center(p: Pos) -> Pt:
    # Use (x,y) = (col,row) so plotting is natural
    r, c = p
    return (float(c), float(r))


def orient(a: Pt, b: Pt, c: Pt) -> float:
    # Cross product (b-a) x (c-a)
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a: Pt, b: Pt, c: Pt) -> bool:
    # c collinear with ab and within bounding box
    return (min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and
            min(a[1], b[1]) <= c[1] <= max(a[1], b[1]))


def segments_intersect(p1: Pt, p2: Pt, q1: Pt, q2: Pt) -> bool:
    """
    Proper segment intersection test (including collinear touching/overlap).
    For this problem, any intersection is disallowed.
    """
    o1 = orient(p1, p2, q1)
    o2 = orient(p1, p2, q2)
    o3 = orient(q1, q2, p1)
    o4 = orient(q1, q2, p2)

    # General case
    if (o1 == 0 and on_segment(p1, p2, q1)) or (o2 == 0 and on_segment(p1, p2, q2)) \
       or (o3 == 0 and on_segment(q1, q2, p1)) or (o4 == 0 and on_segment(q1, q2, p2)):
        return True

    return (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0)


def would_self_intersect(path: List[Pos], new_pos: Pos) -> bool:
    """
    Check if adding segment from last(path) -> new_pos intersects any prior segment.
    We skip the immediately previous segment (which shares the last endpoint).
    """
    if len(path) < 2:
        return False

    p_last = center(path[-1])
    p_new = center(new_pos)

    # Existing segments: (path[i] -> path[i+1]) for i=0..len(path)-2
    # Skip i=len(path)-2 (the last existing segment), since it shares endpoint path[-1]
    for i in range(0, len(path) - 2):
        a = center(path[i])
        b = center(path[i + 1])
        if segments_intersect(a, b, p_last, p_new):
            return True
    return False


class KingTourDFSNoCross:
    def __init__(self, n=4, seed=None, warnsdorff_noise=0.35):
        self.n = n
        self.N = n * n
        self.warnsdorff_noise = warnsdorff_noise
        if seed is not None:
            random.seed(seed)

        self.adj = {(r, c): king_neighbors((r, c), n)
                    for r in range(n) for c in range(n)}
        self.calls = 0

    def onward_count(self, p: Pos, visited: np.ndarray) -> int:
        return sum(not visited[q] for q in self.adj[p])

    def ordered_moves(self, p: Pos, visited: np.ndarray) -> List[Pos]:
        moves = []
        for q in self.adj[p]:
            if not visited[q]:
                score = self.onward_count(q, visited) + random.random() * self.warnsdorff_noise
                moves.append((score, q))
        moves.sort(key=lambda t: t[0])
        # extra randomness among top few
        k = min(3, len(moves))
        head = moves[:k]
        random.shuffle(head)
        moves = head + moves[k:]
        return [q for _, q in moves]

    def dfs(self, p: Pos, visited: np.ndarray, path: List[Pos]) -> bool:
        self.calls += 1
        if len(path) == self.N:
            return True

        for q in self.ordered_moves(p, visited):
            if would_self_intersect(path, q):
                continue

            visited[q] = True
            path.append(q)
            if self.dfs(q, visited, path):
                return True
            path.pop()
            visited[q] = False
        return False

    def solve(self, start: Optional[Pos] = None, max_restarts: int = 200) -> Optional[List[Pos]]:
        for _ in range(max_restarts):
            self.calls = 0
            s = start if start is not None else (random.randrange(self.n), random.randrange(self.n))

            visited = np.zeros((self.n, self.n), dtype=bool)
            visited[s] = True
            path = [s]

            if self.dfs(s, visited, path):
                return path
        return None


def plot_tour(path: List[Pos], n=4):
    board = np.fromfunction(lambda r, c: (r + c) % 2, (n, n))
    fig, ax = plt.subplots(figsize=(5,5), frameon=False)
    ax.imshow(board, cmap="Greys", extent=(-0.5, n - 0.5, n - 0.5, -0.5))

    xs = [c for r, c in path]
    ys = [r for r, c in path]

    ax.plot(xs, ys, linewidth=20)
    #ax.scatter(xs[0], ys[0], s=90, label="start")
    #ax.scatter(xs[-1], ys[-1], s=90, marker="s", label="end")

    for i, (r, c) in enumerate(path):
        ax.text(c, r, str(i), ha="center", va="center", fontsize=10)

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

    #ax.set_title("4×4 king-move tour (no self-intersections)")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.grid(True)
    #ax.legend()
    plt.tight_layout()
    plt.show()


def main():
    idx = None if len(argv) < 2 else int(argv[1])

    if idx is not None:
        with open('king.json', 'r') as catalog:
            path = [(b,a) for a,b in json.load(catalog)[idx]]
    else:
        solver = KingTourDFSNoCross(n=4, seed=None, warnsdorff_noise=0.35)
        path = solver.solve(start=None, max_restarts=200)
        print([[b,a] for a,b in path])

    if path is None:
        raise RuntimeError("No non-self-intersecting tour found within restart limit. Try increasing max_restarts.")

    #print(f"Found path of length {len(path)} in {solver.calls} DFS calls (last attempt).")
    plot_tour(path, n=4)


if __name__ == "__main__":
    main()
