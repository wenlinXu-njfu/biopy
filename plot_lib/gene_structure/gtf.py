#!/usr/bin/env python
"""
File: gtf.py
Description: Plot gene structure based on GTF annotation file
CreateDate: 2022/4/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from click import Choice
import matplotlib.pyplot as plt
from pybioinformatic.gtf import Gtf


def plot_gene_structure(gtf_file: str, exon_color: str, intron_color: str = 'black', edge_color: str = None,
                        figure_width: float = 20, figure_height: float = 10,
                        out_path: str = './', out_prefix: str = 'gene_structure', out_suffix: str = 'pdf',
                        exon_hatch: Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None):
    """Plot gene structure based on GTF annotation file."""
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    d = Gtf(gtf_file).get_gene_dict()  # {gene_id: [{start: int, end: int, strand: str}, {}, ...], ...}
    plt.figure(figsize=(figure_width, figure_height))
    exon = intron = None
    i = 0
    for gene_id, l in d.items():
        left = 0
        for exon_dict in l:
            width = exon_dict['end'] - exon_dict['start'] + 1
            exon = plt.barh(i + 1, width, left=exon_dict['start'], color=exon_color, edgecolor=edge_color, hatch=exon_hatch, linewidth=0.3)
            if left:
                right = exon_dict['start'] - 1
                intron, = plt.plot([left, right], [i + 1, i + 1], color=intron_color, linewidth=0.3)
                left = exon_dict['end'] + 1
            else:
                left = exon_dict['end'] + 1
        i += 1
    # Set ticks and ticks label of y
    label_list = list(d.keys())
    plt.yticks(range(1, len(label_list) + 1), label_list, size=8)
    # Hide y-axis and three frames
    plt.tick_params('y', color='w')
    ax = plt.gcf().gca()
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    ax.spines['left'].set_color(None)
    ax.spines['bottom'].set_position(('data', 0))
    # Set ylim
    if i <= 3:
        ax.set_ylim([0, len(label_list) * 10])
    # Set legend
    plt.legend([exon, intron], ['exon', 'intron'], labelcolor='g', shadow=True, loc=3, bbox_to_anchor=(0.99, 0.99))
    # Save figure
    plt.savefig(f'{out_path}/{out_prefix}.{out_suffix}', bbox_inches='tight', dpi=300)
