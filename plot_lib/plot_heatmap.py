#!/usr/bin/env python
"""
File: plot_heatmap.py
Description: Plot gene expression heatmap
Date: 2022/4/15
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import click
from Biolib.statistics import read_in_gene_expression_as_dataframe
from Biolib.show_info import Displayer


def main(gene_exp_file: str, out_file: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    df = read_in_gene_expression_as_dataframe(gene_exp_file)
    df = df + 1
    df = np.log2(df)
    sns.clustermap(df, cmap='Pastel1', yticklabels=False, col_cluster=False, z_score=0, figsize=(10, 20))
    plt.savefig(out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file', help='Input gene expression file. (support format: txt, xls, xlsx, csv)')
@click.option('-o', '--output_file', 'output_file', default='heatmap.pdf', show_default=True, help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gene_exp_file, output_file):
    """Plot gene expression heatmap."""
    main(gene_exp_file, output_file)


if __name__ == '__main__':
    run()
