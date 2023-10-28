#!/usr/bin/env python
"""
File: plot_cluster_heatmap.py
Description: Plot gene expression cluster heatmap
CreateDate: 2022/4/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gene_exp_file: str, color_map: str, out_file: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    df = pd.read_table(gene_exp_file, index_col=0)
    sns.clustermap(df, cmap=color_map, standard_scale=0)
    plt.savefig(out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file',
              metavar='<exp file>', required=True,
              help='Input gene expression file. (TAB split)')
@click.option('-c', '--color_map', 'color_map',
              metavar='<str>', default='vlag', show_default=True,
              help='Color map.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<figure file>', default='cluster_heatmap.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gene_exp_file, color_map, outfile):
    """Plot gene expression cluster heatmap."""
    main(gene_exp_file, color_map, outfile)


if __name__ == '__main__':
    run()
