#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


def hilbert_points(order: int):
    """Hilbert curve lattice points on an n=2^order grid."""
    dx = [1, 0, -1, 0]  # right, up, left, down
    dy = [0, 1, 0, -1]

    x = y = 0
    direction = 0
    pts = [(x, y)]

    def forward():
        nonlocal x, y
        x += dx[direction]
        y += dy[direction]
        pts.append((x, y))

    def turn_left():
        nonlocal direction
        direction = (direction + 1) % 4

    def turn_right():
        nonlocal direction
        direction = (direction - 1) % 4

    def rec(level, angle):
        if level == 0:
            return

        (turn_left if angle == 90 else turn_right)()
        rec(level - 1, -angle)
        forward()

        (turn_right if angle == 90 else turn_left)()
        rec(level - 1, angle)
        forward()

        rec(level - 1, angle)

        (turn_right if angle == 90 else turn_left)()
        forward()
        rec(level - 1, -angle)

        (turn_left if angle == 90 else turn_right)()

    rec(order, 90)
    return pts


def plot_hilbert(order: int):
    n = 2 ** order

    # checkerboard background (same style)
    board = np.fromfunction(lambda r, c: (r + c) % 2, (n, n))
    fig, ax = plt.subplots(figsize=(5,5), frameon=False)
    ax.imshow(board, cmap="Greys", extent=(-0.5, n - 0.5, n - 0.5, -0.5))

    pts = hilbert_points(order)
    xs = [x for x, y in pts]
    ys = [y for x, y in pts]

    # curve
    ax.plot(xs, ys, linewidth=20)
    #ax.scatter(xs[0], ys[0], s=90, label="start")
    #ax.scatter(xs[-1], ys[-1], s=90, marker="s", label="end")

    # visit index at EVERY lattice point
    for i, (x, y) in enumerate(pts):
        ax.text(x, y, str(i), ha="center", va="center", fontsize=10)

    #ax.set_title(f"Hilbert curve (order {order}, {n}×{n})")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)  # flipped y like your script
    ax.grid(True)
    #ax.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    order = 2  # recommended ≤ 4 for legibility
    plot_hilbert(order)
