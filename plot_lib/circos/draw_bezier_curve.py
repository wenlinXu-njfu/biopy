"""
File: draw_bezier_curve.py
Description: Draw bezier curve.
Date: 2022/10/11
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import math
import numpy as np
import matplotlib.axes
import matplotlib.pyplot as plt
from plot_lib.plot_base import _bezier_curve
from plot_lib.circos.draw_chromosome import draw_chr
from plot_lib.circos.draw_outer_bar import draw_outer_bar


def coordinates_conversion(theta: float, r: float):
    """
    Convert polar coordinates to rectangular coordinates.
    """
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    return [x, y]


def draw_bezier_curve(axes: matplotlib.axes.Axes, x1: float, y1: float, x2: float, y2: float,
                      line_width: float = 0.8, line_color: str = 'grey', label: str = None, alpha: float = 0.3):
    """
    Draw bezier curve.
    :param axes: Polar coordinates object (type=matplotlib.axes.Axes)
    :param x1: The x-polar coordinates of point P0. (type=float)
    :param y1: The y-polar coordinates of point P0. (type=float)
    :param x2: The x-polar coordinates of point P2. (type=float)
    :param y2: The y-polar coordinates of point P2. (type=float)
    :param line_width: Line width. (type=float, default=0.08)
    :param line_color: Line color. (type=str, default='grey')
    :param label: Line2D label. (type=str, default=None)
    :param alpha: Line alpha. (type=float, default=0.5)
    :return: tuple(ax, line2D)
    """
    #  Convert polar coordinates to rectangular coordinates.
    P0, P1, P2 = np.array([coordinates_conversion(x1, y1),
                           coordinates_conversion((x1 + x2) / 2, 0),
                           coordinates_conversion(x2, y2)])
    # Get the coordinates of the points on the bezier curve.
    x, y = _bezier_curve(P0, P1, P2)
    # Draw the bezier curve.
    if label:
        line2D, = axes.plot(x, y, transform=axes.transData._b, linewidth=line_width, color=line_color, label=label, alpha=alpha)
    else:
        line2D, = axes.plot(x, y, transform=axes.transData._b, linewidth=line_width, color=line_color, alpha=alpha)
    return axes, line2D


# Test
if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42

    # Draw chromosome bar.
    len_dict, theta_dict, width_dict, ax, bar_container = draw_chr('test_data/Ptc_chr_len.txt')

    # Draw outer bar.
    ax, artist_list = draw_outer_bar(theta_dict, width_dict, ax, 'test_data/stat.txt', 7.8)

    # Draw bezier curve.
    used_color = []
    for line in open('test_data/link.txt'):
        split = line.strip().split('\t')
        chr1, start1, end1 = split[0], int(split[1]), int(split[2])
        chr2, start2, end2 = split[3], int(split[4]), int(split[5])
        color, desc = split[-2], split[-1]
        angle1 = theta_dict[chr1] - width_dict[chr1] / 2
        angle2 = theta_dict[chr2] - width_dict[chr2] / 2
        x1 = (end1 + start1) / 2 / len_dict[chr1] * width_dict[chr1] + angle1
        x2 = (end2 + start2) / 2 / len_dict[chr2] * width_dict[chr2] + angle2
        if color not in used_color:
            used_color.append(color)
            _, line_2D = draw_bezier_curve(ax, x1, 7.7, x2, 7.7, 1, color, desc)
            artist_list.append(line_2D)
        else:
            draw_bezier_curve(ax, x1, 7.7, x2, 7.7, 1, color)

    # Set legend.
    artist_list.insert(0, bar_container)
    plt.legend(artist_list, [i.get_label() for i in artist_list], fontsize=5, loc=(0.9, 0.9), ncol=2)

    # save or show picture
    # plt.savefig('E:/Desktop/1.pdf', bbox_inches='tight')
    plt.show()
