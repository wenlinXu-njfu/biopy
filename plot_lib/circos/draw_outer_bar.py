"""
File: draw_outer_bar.py
Description: Draw outer bar.
Date: 2022/10/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import math
from typing import Union
import matplotlib.axes
import matplotlib.pyplot as plt
from plot_lib.circos.draw_chromosome import draw_chr


def draw_outer_bar(chr_theta_dict: dict, chr_width_dict: dict, axes: matplotlib.axes.Axes, stat_file,
                   bottom: Union[float, list] = 10.05):
    """
    Draw the statistical bar chart of the outermost circle of the circos graph.
    :param chr_theta_dict: Chromosome theta angle dict.
    :param chr_width_dict:Chromosome width angle dict.
    :param axes: Axes object of matplotlib.axes.Axes.
    :param stat_file: Feature statistics file. (Chr_num\\tFeature_type\\tCount\\tColor\\n)
    :param bottom: Y-axis coordinate of each bar bottom.
    :return: tuple(axes, bar_list)
             axes -> Axes object of matplotlib.axes.Axes.
             bar_list -> Bar container.
    """
    stat_dict = {}  # {Chr01: {'type': [str, ...], 'count': [int, ...], 'color': [str, ...]}, ...}
    all_count = []
    for line in open(stat_file):
        split = line.strip().split('\t')
        chr_num, _type, count, color = split[0], split[1], math.log2(int(split[2]) + 1), split[3]
        all_count.append(count)
        if chr_num in stat_dict:
            stat_dict[chr_num]['type'].append(_type)
            stat_dict[chr_num]['count'].append(count)
            stat_dict[chr_num]['color'].append(color)
        else:
            stat_dict[chr_num] = {'type': [_type], 'count': [count], 'color': [color]}
    j = 0
    for chr_num, d in stat_dict.items():
        bar_list = []
        total_bar_num = len(d['type']) + 2
        chr_theta, chr_width = chr_theta_dict[chr_num], chr_width_dict[chr_num]
        bar_width = chr_width / total_bar_num
        Xi = chr_theta - chr_width / 2 + 3 / 2 * bar_width
        x = [Xi + i * bar_width for i in range(len(d['type']))]
        height = [d['count'][i] / max(all_count) for i in range(len(d['type']))]
        if isinstance(bottom, float):
            k = 0
            while k < len(d['type']):
                height = d['count'][k] / max(all_count)
                bar_list.append(axes.bar(Xi, height, bar_width, bottom, color=d['color'][k], label=d['type'][k]))
                Xi += bar_width
                k += 1
        else:
            rets = axes.bar(x, height, bar_width, [bottom[j] for _ in range(len(d['type']))],
                            color=d['color'], label=d['type'])
            n = 0
            for ret in rets:
                ret.set_label(d['type'][n])
                bar_list.append(ret)
                n += 1
        j += 1
    return axes, bar_list


# Test
if __name__ == '__main__':
    length, theta, width, ax, bar_container = draw_chr('test_data/Ptc_chr_len.txt')
    ax, bars = draw_outer_bar(theta, width, ax, 'test_data/stat.txt')
    bars.insert(0, bar_container)
    plt.legend(bars, ['chromosome', 'exonic', 'intronic', 'intergenic', 'antisense'],
               fontsize=10, loc=10, bbox_to_anchor=(0.5, 0.5))
    plt.show()
