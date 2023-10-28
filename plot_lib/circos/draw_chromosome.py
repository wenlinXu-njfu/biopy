"""
File: draw_chromosome.py
Description: Draw chromosome bar.
CreateDate: 2022/10/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
import numpy as np
import matplotlib.pyplot as plt
from plot_lib.circos.circos_base import get_chr_theta_width


def draw_chr(chr_len_file, height: float = 1, bottom: Union[float, list] = 9.0,
             facecolor: str = 'lightgrey', edgecolor: str = 'black',
             linewidth: float = 0.4, font_size: float = 6, figsize: tuple = (6.4, 4.8)):
    """
    Draw chromosome bar.
    :param chr_len_file: Chromosome length file. (Chr_num\\tLength\\n)
    :param height: Bar height of each chromosome.
    :param bottom: Y-axis of bar bottom for each chromosome.
    :param facecolor: Face color of each chromosome bar.
    :param edgecolor: Edge color of each chromosome bar.
    :param linewidth: Line width of each chromosome bar.
    :param font_size: Font size of each chromosome bar.
    :param figsize: Figure size.
    :return: tuple(theta: dict, width: dict, ax: matplotlib.axes.Axes, BarContainer: matplotlib.container.BarContainer)
             length -> raw chromosome length dict. {Chr_num: length, ...}
             theta -> chromosome theta dict. {Chr_num: theta, ...}
             width -> chromosome width dict. {Chr_num: width, ...}
             ax -> Axes object in matplotlib.
             BarContainer -> Chromosome bar container.
    """
    length, theta, width = get_chr_theta_width(chr_len_file)
    fig = plt.figure(figsize=figsize, dpi=200)
    ax = fig.add_subplot(111, polar=True, frame_on=False)
    BarContainer = ax.bar(theta.values(), height, width.values(), bottom,
                          facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, label='chromosome')
    i = 0
    for chr_num, chr_theta in theta.items():
        if isinstance(bottom, float):
            if 0 <= chr_theta < np.pi:
                plt.text(chr_theta, bottom + height / 2, chr_num, rotation=180 / np.pi * chr_theta - 90,
                         horizontalalignment='center', verticalalignment='center', fontsize=font_size)
            else:
                plt.text(chr_theta, bottom + height / 2, chr_num, rotation=180 * chr_theta / np.pi - 270,
                         horizontalalignment='center', verticalalignment='center', fontsize=font_size)
        elif isinstance(bottom, list):
            if 0 <= chr_theta < np.pi:
                plt.text(chr_theta, bottom[i] + height / 2, chr_num, rotation=180 / np.pi * chr_theta - 90,
                         horizontalalignment='center', verticalalignment='center', fontsize=font_size)
            else:
                plt.text(chr_theta, bottom[i] + height / 2, chr_num, rotation=180 * chr_theta / np.pi - 270,
                         horizontalalignment='center', verticalalignment='center', fontsize=font_size)
            i += 1
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return length, theta, width, ax, BarContainer


# Test
if __name__ == '__main__':
    draw_chr('test_data/Ptc_chr_len.txt', bottom=[9, 9, 9, 9, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10])
    # plt.legend()
    plt.show()
