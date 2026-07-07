#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon


# -------------------------
# Gosper (Peano–Hilbert) L-system
# -------------------------

def gosper_lsystem(order: int) -> str:
    s = "A"
    for _ in range(order):
        out = []
        for ch in s:
            if ch == "A":
                out.append("A-B--B+A++AA+B-")
            elif ch == "B":
                out.append("+A-BB--B-A++A+B")
            else:
                out.append(ch)
        s = "".join(out)
    return s


# -------------------------
# Hex geometry (same as your script)
# -------------------------

def axial_to_xy(q, r, size=1.0):
    """
    Pointy-top axial coordinates → 2D.
    Matches your hex.py exactly.
    """
    x = size * np.sqrt(3) * (q + 0.5 * r)
    y = size * 1.5 * r
    return x, y


# 6 hex directions in axial coords (pointy-top)
HEX_DIRS = [
    (1, 0), (1, -1), (0, -1),
    (-1, 0), (-1, 1), (0, 1)
]


def gosper_hex_curve(order: int):
    instr = gosper_lsystem(order)

    # turtle state
    q, r = 0, 0
    heading = 0

    path = [(q, r)]

    for ch in instr:
        if ch == "+":
            heading = (heading + 1) % 6
        elif ch == "-":
            heading = (heading - 1) % 6
        elif ch in ("A", "B"):
            dq, dr = HEX_DIRS[heading]
            q += dq
            r += dr
            path.append((q, r))

    return path


# -------------------------
# Plotting (style-matched)
# -------------------------

def plot_gosper_with_hexes(path, hex_size=1.0, annotate_every: int = 4):
    fig, ax = plt.subplots(figsize=(8,8), frameon=False)

    cells = set(path)
    pos2xy = {h: axial_to_xy(h[0], h[1], hex_size) for h in cells}

    # draw hex tiles (checker-ish parity like your script)
    for q, r in cells:
        x, y = pos2xy[(q, r)]
        parity = (q + r) & 1
        face = "white" if parity == 0 else "#e8e8e8"

        hex_patch = RegularPolygon(
            (x, y),
            numVertices=6,
            radius=hex_size,
            orientation=np.radians(30),  # pointy-top
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

    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    #ax.set_title(f"Hexagonal Hilbert (Gosper) curve (cells={len(cells)})")
    #ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


# -------------------------
# Main
# -------------------------

if __name__ == "__main__":
    order = 2   # try 1–6; grows fast
    path = gosper_hex_curve(order)
    plot_gosper_with_hexes(path, hex_size=1.0, annotate_every=1)
