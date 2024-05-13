#!/usr/bin/env python
"""
File: gff.py
Description: Plot mRNA structure based on GFF annotation file
CreateDate: 2022/4/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from click import Choice
import matplotlib.pyplot as plt
from pybioinformatic.gff import Gff


def plot_mRNA_structure(gff_file: str,
                        utr_color: str,
                        cds_color: str,
                        edge_color: str = None,
                        figure_width: float = 20,
                        figure_height: float = 10,
                        out_file: str = 'mRNA_structure.pdf',
                        utr_hatch: Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None,
                        cds_hatch: Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None):
    """Plot mRNA structure based on GFF annotation file."""
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    d = Gff(gff_file).get_mRNA_dict()  # {mRNA_id: [{feature_type: str, start: int, end: int, strand: str}, ...], ...}
    plt.figure(figsize=(figure_width, figure_height))
    utr = cds = intron = None
    i = 0  # mRNA number
    for mRNA_id, l in d.items():
        l.sort(key=lambda item: -item['start'])
        right = 0
        for feature_dict in l:
            width = feature_dict['end'] - feature_dict['start'] + 1
            if feature_dict['feature_type'] == 'three_prime_UTR':
                if feature_dict == l[0]:
                    utr = plt.barh(
                        y=i + 1,
                        width=width,
                        left=feature_dict['start'],
                        color=utr_color,
                        edgecolor=edge_color,
                        hatch=utr_hatch,
                        linewidth=0.3,
                        zorder=1
                    )
                    right = feature_dict['start'] - 1
                else:
                    left = feature_dict['end'] - 1
                    if left + 1 != right:
                        intron, = plt.plot([right, left], [i + 1, i + 1], color='k', linewidth=0.4, zorder=0)
                    plt.barh(
                        y=i + 1,
                        width=width,
                        left=feature_dict['start'],
                        color=utr_color,
                        edgecolor=edge_color,
                        hatch=utr_hatch,
                        linewidth=0.3,
                        zorder=1
                    )
                    right = feature_dict['start'] - 1
            elif feature_dict['feature_type'] == 'five_prime_UTR':
                left = feature_dict['end'] - 1
                if left + 1 != right:
                    intron, = plt.plot([left, right], [i + 1, i + 1], color='k', linewidth=0.4, zorder=0)
                plt.barh(
                    y=i + 1,
                    width=width,
                    left=feature_dict['start'],
                    color=utr_color,
                    edgecolor=edge_color,
                    hatch=utr_hatch,
                    linewidth=0.3,
                    zorder=1
                )
                right = feature_dict['start'] - 1
            elif feature_dict['feature_type'] == 'CDS':
                if feature_dict == l[0]:
                    cds = plt.barh(
                        y=i + 1,
                        width=width,
                        left=feature_dict['start'],
                        color=cds_color,
                        edgecolor=edge_color,
                        hatch=cds_hatch,
                        linewidth=0.3,
                        zorder=1
                    )
                    right = feature_dict['start'] - 1
                else:
                    left = feature_dict['end'] - 1
                    if left + 1 != right:
                        intron, = plt.plot([left, right], [i + 1, i + 1], color='k', linewidth=0.4, zorder=0)
                    cds = plt.barh(
                        y=i + 1,
                        width=width,
                        left=feature_dict['start'],
                        color=cds_color,
                        edgecolor=edge_color,
                        hatch=cds_hatch,
                        linewidth=0.3,
                        zorder=1
                    )
                    right = feature_dict['start'] - 1
        i += 1
    # Set ticks and ticks label of y
    label_list = list(d.keys())
    plt.yticks(range(1, len(label_list) + 1), label_list, size=8)
    # Hide y-axis and three frames
    plt.tick_params('y', length=0)
    ax = plt.gcf().gca()
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    ax.spines['left'].set_color(None)
    ax.spines['bottom'].set_position(('data', 0))
    # Set ylim
    if i <= 3:
        ax.set_ylim([0, len(label_list) * 10])
    # Set legend
    plt.legend([utr, cds, intron], ['UTR', 'CDS', 'intron'], labelcolor='g', shadow=True, loc=3, bbox_to_anchor=(0.9, 0.9))
    # Save figure
    plt.savefig(out_file, bbox_inches='tight', dpi=300)
