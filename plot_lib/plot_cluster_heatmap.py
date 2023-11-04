#!/usr/bin/env python
"""
File: plot_cluster_heatmap.py
Description: Plot gene expression cluster heatmap.
CreateDate: 2022/4/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from numpy import log2
from matplotlib.pyplot import rcParams, savefig
from seaborn import clustermap
from pybioinformatic import read_in_gene_expression_as_dataframe, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gene_exp_file: str,
         color_map: str,
         row_cluster: bool,
         col_cluster: bool,
         out_file: str):
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'
    df = read_in_gene_expression_as_dataframe(gene_exp_file)
    df = log2(df + 1)
    df.index.name = None
    if len(df) <= 150:
        tree_kws = {'linewidth': 0.3}
    elif 150 < len(df) <= 500:
        tree_kws = {'linewidth': 0.1}
    else:
        tree_kws = {'linewidth': 0.01, 'alpha': 0.3}
    clustermap(df, cmap=color_map, yticklabels=False, tree_kws=tree_kws,
               row_cluster=row_cluster, col_cluster=col_cluster,
               cbar_pos=(1.02, 0.5, 0.025, 0.25))
    savefig(out_file, bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file',
              metavar='<text file>', required=True,
              help='Input gene expression file (support format: txt, csv, xls, and xlsx).')
@click.option('-c', '--color_map', 'color_map',
              metavar='<str>', default='PiYG', show_default=True,
              help='Color map (eg. PiYG, vlag, OrRd, YlOrRd, Spectral).')
@click.option('--row_cluster', 'row_cluster',
              type=bool, is_flag=True, flag_value=True,
              help='If specified, cluster the rows.')
@click.option('--col_cluster', 'col_cluster',
              type=bool, is_flag=True, flag_value=True,
              help='If specified, cluster the columns.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<figure file>', default='cluster_heatmap.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gene_exp_file, color_map, row_cluster, col_cluster, outfile):
    """Plot gene expression cluster heatmap."""
    main(gene_exp_file, color_map, row_cluster, col_cluster, outfile)


if __name__ == '__main__':
    run()
