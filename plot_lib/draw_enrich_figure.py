#!/usr/bin/env python
"""
File: draw_enrich_figure.py
Description: Visualize cluster profiler enrich results.
CreateDate: 2024/10/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union
from re import findall
from pandas import read_table
import matplotlib.pyplot as plt
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(enrich_results_file: Union[str, TextIOWrapper],
         color_map: str,
         output_file: str):
    df = read_table(enrich_results_file, index_col=0)
    df['Rich_Factor'] = df.apply(
        lambda row: int(row['GeneRatio'].split('/')[0]) / int(row['BgRatio'].split('/')[0]),
        axis=1
    )
    # adjust dot size
    if 0 <= df.Count.max() - df.Count.min() <= 11:
        dot_size_fold = 20
        df['Size'] = df.Count * dot_size_fold
    elif 11 < df.Count.max() - df.Count.min() < 41:
        dot_size_fold = 5
        df['Size'] = df.Count * dot_size_fold
    elif 41 <= df.Count.max() - df.Count.min() < 51:
        dot_size_fold = 3
        df['Size'] = df.Count * dot_size_fold
    elif 51 <= df.Count.max() - df.Count.min() < 61:
        dot_size_fold = 2
        df['Size'] = df.Count * dot_size_fold
    else:
        dot_size_fold = 1
        df['Size'] = df.Count
    df.sort_values(by='p.adjust', inplace=True)
    df = df.head(n=20)

    plt.style.use('ggplot')
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(111)
    scatter = ax.scatter(
        x='Rich_Factor',
        y=range(len(df)),
        c='p.adjust',
        s='Size',
        cmap=color_map,
        data=df
    )

    cbar = plt.colorbar(scatter)  # add color bar
    cbar.set_label(label='P value', loc='center')
    if len(df['Size'].unique()) > 5:
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, markersize=100, num=5)
    else:
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, markersize=100)
    if dot_size_fold != 1:
        labels = ['%.0f' % (int(findall(r'\d+', i)[0]) / dot_size_fold) for i in labels]
    ax.legend(handles, labels, loc="upper right", title="Count")

    plt.xlabel('Rich Factor')
    plt.yticks(ticks=range(len(df)), labels=df.Description)
    plt.savefig(output_file, bbox_inches='tight', dpi=300)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--enrich-results', 'enrich_results',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input cluster profiler enrich results file.'
                   r'(header=ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount)')
@click.option('-c', '--color-map', 'color_map',
              metavar='<str>', default='RdYlBu', show_default=True,
              help='Color map. '
                   'More reference see \033[1;34mhttps://matplotlib.org/stable/users/explain/colors/colormaps.html\033[0m')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file>', default='enrich_plot.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(enrich_results, color_map, output_file):
    """Visualize cluster profiler enrich results."""
    main(enrich_results, color_map, output_file)


if __name__ == '__main__':
    run()
