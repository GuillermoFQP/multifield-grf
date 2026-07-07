import numpy as np
import matplotlib.pyplot as plt


# Peano–Gosper (a.k.a. Gosper “flowsnake”) curve:
# a classic “hexagonal Hilbert-like” space-filling curve on a hex lattice.
#
# L-system rules (standard):
#   A -> A-B--B+A++AA+B-
#   B -> +A-BB--B-A++A+B
# with turtle:
#   + : turn left 60°
#   - : turn right 60°
#   A,B : draw forward one step


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


def gosper_points(order: int, step: float = 1.0):
    instr = gosper_lsystem(order)

    # 6 directions on a hex grid (pointy-top axial projected to 2D)
    # We'll keep a heading index 0..5. Turning left increments, right decrements.
    # Direction vectors for angles 0, 60, 120, 180, 240, 300 degrees.
    angles = np.deg2rad([0, 60, 120, 180, 240, 300])
    dx = np.cos(angles) * step
    dy = np.sin(angles) * step

    heading = 0
    x, y = 0.0, 0.0
    xs = [x]
    ys = [y]

    for ch in instr:
        if ch == "+":
            heading = (heading + 1) % 6
        elif ch == "-":
            heading = (heading - 1) % 6
        elif ch == "A" or ch == "B":
            x += dx[heading]
            y += dy[heading]
            xs.append(x)
            ys.append(y)
        else:
            # ignore anything else (not expected)
            pass

    return np.array(xs), np.array(ys)


if __name__ == "__main__":
    order = 4      # try 1..6 (grows fast)
    step = 1.0

    xs, ys = gosper_points(order, step=step)

    fig = plt.figure(figsize=(7,7), frameon=False)
    plt.plot(xs, ys, linewidth=3)
    plt.axis("equal")
    plt.axis("off")
    #plt.title(f"Hexagonal 'Hilbert-like' curve (Peano–Gosper), order {order}")
    plt.tight_layout()
    plt.show()
