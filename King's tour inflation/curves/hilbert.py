import matplotlib.pyplot as plt

def hilbert_points(order):
    """
    Generate integer grid points for a Hilbert curve of given order.
    Returns a list of (x, y) points on a (2^order x 2^order) grid.
    """
    # directions: 0 = right, 1 = up, 2 = left, 3 = down
    dx = [1, 0, -1, 0]
    dy = [0, 1, 0, -1]

    x, y = 0, 0          # starting point
    direction = 0        # initial direction: right
    points = [(x, y)]    # list of points along the curve

    def forward():
        nonlocal x, y
        x += dx[direction]
        y += dy[direction]
        points.append((x, y))

    def hilbert(level, angle):
        """
        Recursive construction using a turtle-like algorithm.
        angle = +90 or -90 (degrees), encoded via direction changes.
        """
        nonlocal direction
        if level == 0:
            return

        # turn(angle)
        if angle == 90:
            direction = (direction + 1) % 4
        else:  # angle == -90
            direction = (direction - 1) % 4

        hilbert(level - 1, -angle)
        forward()

        # turn(-angle)
        if angle == 90:
            direction = (direction - 1) % 4
        else:
            direction = (direction + 1) % 4

        hilbert(level - 1, angle)
        forward()

        hilbert(level - 1, angle)

        # turn(-angle)
        if angle == 90:
            direction = (direction - 1) % 4
        else:
            direction = (direction + 1) % 4

        forward()
        hilbert(level - 1, -angle)

        # turn(angle)
        if angle == 90:
            direction = (direction + 1) % 4
        else:
            direction = (direction - 1) % 4

    hilbert(order, 90)
    return points


def hilbert_curve(order):
    """
    Return x, y coordinates of a Hilbert curve of given order
    scaled to the unit square [0, 1] x [0, 1].
    """
    pts = hilbert_points(order)
    n = 2**order - 1  # max coordinate in grid
    scale = 1.0 / n

    xs = [x * scale for x, y in pts]
    ys = [y * scale for x, y in pts]
    return xs, ys


if __name__ == "__main__":
    order = 6  # try changing this to 1..7; higher = more detailed curve
    xs, ys = hilbert_curve(order)

    fig = plt.figure(figsize=(6,6), frameon=False)
    plt.plot(xs, ys, linewidth=3)
    plt.axis("equal")
    plt.axis("off")
    #plt.title(f"Hilbert Curve (order {order})")
    plt.tight_layout()
    plt.show()
