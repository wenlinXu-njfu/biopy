"""
File: draw_gene_density2.py
Description: 
CreateDate: 2022/10/13
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
import numpy as np
import matplotlib.axes
import matplotlib.pyplot as plt
from plot_lib.circos.draw_chromosome import draw_chr
from plot_lib.circos.draw_outer_bar import draw_outer_bar
from plot_lib.circos.draw_bezier_curve import draw_bezier_curve


def draw_gene_density(gene_density_file: str, chr_theta_dict: dict, chr_width_dict: dict, axes: matplotlib.axes.Axes,
                      bottom: Union[float, list] = 7.7, color: str = None, label: str = 'gene density'):
    """
    Draw gene density.
    :param gene_density_file: Gene density file. (Chr_num\\tStart\\tEnd\\tCount\\n)
    :param chr_theta_dict: Theta angle for each chromosome bar.
    :param chr_width_dict: Width for each chromosome bar.
    :param axes: Axes object of matplotlib.axes.Axes.
    :param bottom: Y-axis coordinate of bottom for gene density plot.
    :param color: Gene density plot color.
    :param label: Gene density plot label.
    :return: tuple(ax, Line2D_obj)
    """
    current_chr = None
    y = []
    total_line = sum(1 for _ in open(gene_density_file))
    line_count = 0
    j = 0
    for line in open(gene_density_file):
        line_count += 1
        split = line.strip().split('\t')
        if line_count != total_line:
            if not current_chr:
                current_chr = split[0]
                y.append(float(split[-1]))
            else:
                if current_chr != split[0]:
                    if isinstance(bottom, list):
                        current_bottom = bottom[j]
                        j += 1
                    else:
                        current_bottom = bottom
                    x = np.linspace(chr_theta_dict[current_chr] - chr_width_dict[current_chr] / 2,
                                    chr_theta_dict[current_chr] + chr_width_dict[current_chr] / 2, len(y))
                    y = [current_bottom + 0.1 + i / max(y) for i in y]
                    line2D, = axes.plot(x, y, linewidth=0.3, color=color, label=label)
                    axes.bar(chr_theta_dict[current_chr], 1.1, chr_width_dict[current_chr], current_bottom,
                             color='w', edgecolor='k', linewidth=0.5)
                    current_chr = split[0]
                    y = [float(line.strip().split('\t')[-1])]
                else:
                    y.append(float(line.strip().split('\t')[-1]))
        else:
            if isinstance(bottom, list):
                current_bottom = bottom[j]
                j += 1
            else:
                current_bottom = bottom
            x = np.linspace(chr_theta_dict[current_chr] - chr_width_dict[current_chr] / 2,
                            chr_theta_dict[current_chr] + chr_width_dict[current_chr] / 2, len(y))
            y = [current_bottom + 0.1 + i / max(y) for i in y]
            axes.plot(x, y, linewidth=0.3, color=color)
            axes.bar(chr_theta_dict[current_chr], 1.1, chr_width_dict[current_chr], current_bottom,
                     color='w', edgecolor='k', linewidth=0.5)
    return axes, line2D


# Test
if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42

    # Draw chromosome bar.
    # length, theta, width, ax, bar_container = draw_chr('test_data/Ptc_chr_len.txt')
    length, theta, width, ax, bar_container = draw_chr('test_data/Ptc_chr_len.txt',
                                                       bottom=[9, 9, 9, 9, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10])

    # Draw outer bar.
    # ax, artist_list = draw_outer_bar(theta, width, ax, 'test_data/stat.txt')
    ax, artist_list = draw_outer_bar(theta, width, ax, 'test_data/stat.txt',
                                     bottom=[10.05, 10.05, 10.05, 10.05, 10.05, 11.05, 10.05, 10.05, 10.05, 10.05,
                                             10.05, 10.05, 10.05, 10.05, 10.05, 10.05, 10.05, 10.05, 11.05])

    # Draw gene density.
    draw_gene_density('test_data/gene_density.txt', theta, width, ax,
                      [7.7, 7.7, 7.7, 7.7, 7.7, 8.7, 7.7, 7.7, 7.7, 7.7,
                       7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 8.7])

    # Draw bezier curve.
    used_color = []
    for line2 in open('test_data/link.txt'):
        split2 = line2.strip().split('\t')
        chr1, start1, end1 = split2[0], int(split2[1]), int(split2[2])
        chr2, start2, end2 = split2[3], int(split2[4]), int(split2[5])
        angle1 = theta[chr1] - width[chr1] / 2
        angle2 = theta[chr2] - width[chr2] / 2
        x1 = (end1 + start1) / 2 / length[chr1] * width[chr1] + angle1
        x2 = (end2 + start2) / 2 / length[chr2] * width[chr2] + angle2
        if chr1 == 'Chr06' or chr1 == 'Chr19':
            y1 = 8.6
        else:
            y1 = 7.6
        if chr2 == 'Chr06' or chr2 == 'Chr19':
            y2 = 8.6
        else:
            y2 = 7.6
        if split2[-1] not in used_color:
            used_color.append(split2[-1])
            _, line_2D = draw_bezier_curve(ax, x1, y1, x2, y2, 0.8, split2[-2], split2[-1])
            artist_list.insert(2, line_2D)
        else:
            draw_bezier_curve(ax, x1, y1, x2, y2, 0.8, split2[-2])

    # Set legend.
    artist_list.insert(0, bar_container)
    plt.legend(artist_list, [i.get_label() for i in artist_list], fontsize=6, loc=(0.9, 0.9), ncol=2)

    # plt.savefig('E:/Desktop/1.pdf', bbox_inches='tight')
    plt.show()
