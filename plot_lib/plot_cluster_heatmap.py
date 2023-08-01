#!/usr/bin/env python
"""
File: plot_cluster_heatmap.py
Description: Plot gene expression cluster heatmap
Date: 2022/4/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(gene_exp_file: str, color_map: str, out_file: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    df = pd.read_table(gene_exp_file, index_col=0)
    sns.clustermap(df, cmap=color_map, standard_scale=0)
    plt.savefig(out_file)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--gene_exp', 'gene_exp_file', help='Input gene expression file. (TAB split)')
@click.option('-c', '--color_map', 'color_map', default='vlag', help='[optional] Color map. {default: vlag}')
@click.option('-o', '--output_file', 'outfile', default='cluster_heatmap.pdf',
              help='[optional] Output file. {default: cluster_heatmap.pdf}')
def run(gene_exp_file, color_map, outfile):
    """Plot gene expression cluster heatmap."""
    main(gene_exp_file, color_map, outfile)


if __name__ == '__main__':
    run()
