import numpy as np
from matplotlib.axes import Axes


def bezier_curve(P0, P1, P2):
    """
    Define a second-order Bézier curve.
    :param P0: A point with known coordinates.
    :param P1: The control points.
    :param P2: Another point with known coordinates.
    :return: The coordinates of the 50 points on the Bézier curve.
    """
    P = lambda t: (1 - t) ** 2 * P0 + 2 * t * (1 - t) * P1 + t ** 2 * P2  # Define the Bezier curve
    points = np.array([P(t) for t in np.linspace(0, 1, 50)])  # Verify Bézier curves at 50 points in the range [0, 1]
    x, y = points[:, 0], points[:, 1]  # Get the x and y coordinates of the points respectively
    return x, y


def draw_bezier_curve(axes: Axes, x1: float, y1: float, x2: float, y2: float, x0: float = None, y0: float = None,
                      line_width: float = 0.8, line_color: str = 'grey', label: str = None, alpha: float = 0.3):
    """
    Draw Bézier curve.
    :param axes: Polar coordinates object (type=matplotlib.axes.Axes)
    :param x1: The x coordinates of point P0. (type=float)
    :param y1: The y coordinates of point P0. (type=float)
    :param x2: The x coordinates of point P2. (type=float)
    :param y2: The y coordinates of point P2. (type=float)
    :param x0: The x coordinates of control point. (type=float, default=(x1 + x2) / 2)
    :param y0: The y coordinates of control point. (type=float, default=(y1 + y2) / 2)
    :param line_width: Line width. (type=float, default=0.08)
    :param line_color: Line color. (type=str, default='grey')
    :param label: Line2D label. (type=str, default=None)
    :param alpha: Line alpha. (type=float, default=0.5)
    :return: tuple(ax, line2D)
    """
    if not x0:
        x0 = (x1 + x2) / 2
    if not y0:
        y0 = (y1 + y2) / 2
    P0, P1, P2 = np.array([[x1, y1], [x0, y0], [x2, y2]])
    # Get the coordinates of the points on the Bézier curve.
    x, y = bezier_curve(P0, P1, P2)
    # Draw the bezier curve.
    if label:
        line2D, = axes.plot(x, y, transform=axes.transData._b, linewidth=line_width, color=line_color, label=label, alpha=alpha)
    else:
        line2D, = axes.plot(x, y, transform=axes.transData._b, linewidth=line_width, color=line_color, alpha=alpha)
    return axes, line2D
