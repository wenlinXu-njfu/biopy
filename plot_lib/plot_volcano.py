#!/usr/bin/env python
"""
File: plot_volcano.py
Description: Plot the volcano of differentially expressed genes
Date: 2022/2/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, List, Union
from Biolib.show_info import Displayer


def main(gene_exp_file: str, log2fc: float, padj: float, figure_size: Union[Tuple, List], out_prefix: str, fmt: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    df = pd.read_table(gene_exp_file, index_col='gene_id')
    df['-lg_p'] = -np.log10(df['padj'])  # calculate -log10(padj) for each gene
    # mark gene
    df['sig'] = 'normal'
    df.loc[(df.log2FoldChange > log2fc) & (df.padj < padj), 'sig'] = 'up'
    df.loc[(df.log2FoldChange < -log2fc) & (df.padj < padj), 'sig'] = 'down'
    # plot_figure scatter figure
    fig = plt.figure(figsize=figure_size)
    plt.scatter(df[df['sig'] == 'up']['log2FoldChange'], df[df['sig'] == 'up']['-lg_p'],
                color='salmon', marker='o', label='up', s=5)
    plt.scatter(df[df['sig'] == 'normal']['log2FoldChange'], df[df['sig'] == 'normal']['-lg_p'],
                color='lightgrey', marker='o', label='normal', s=5)
    plt.scatter(df[df['sig'] == 'down']['log2FoldChange'], df[df['sig'] == 'down']['-lg_p'],
                color='lightgreen', marker='o', label='down', s=5)
    plt.xlabel(r'$log_2Fold$' + ' Change')
    plt.ylabel(r'$-log_{10}P$' + ' adjust')
    plt.xlim([-max(abs(df['log2FoldChange'])), max(abs(df['log2FoldChange']))])
    plt.legend(loc='upper center')
    plt.savefig(f'{out_prefix}.{fmt}', bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file',
              help="Expression matrix file. (header must including 'gene_id', 'padj' and 'log2FoldChange')")
@click.option('-l', '--log2fc', 'log2fc', type=float, default=1.5, show_default=True,
              help='Set log2(fold change) as the threshold.')
@click.option('-p', '--p_adjust', 'p_adjust', type=float, default=0.05, show_default=True, help='Set padj as the threshold.')
@click.option('-s', '--figure_size', 'figure_size', default='6.0x8.0', show_default=True, help='Specify figure size.')
@click.option('-o', '--output_prefix', 'output_prefix', default='volcano', show_default=True,
              help='Prefix of output file.')
@click.option('-O', '--output_format', 'output_format',
              type=click.Choice(['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']),
              default='pdf', show_default=True, help='Specify the format of output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gene_exp_file, log2fc, p_adjust, figure_size, output_prefix, output_format):
    """Plot the volcano of differentially expressed genes."""
    figure_size = [float(i) for i in figure_size.split(',')]
    main(gene_exp_file, log2fc, p_adjust, figure_size, output_prefix, output_format)


if __name__ == '__main__':
    run()
