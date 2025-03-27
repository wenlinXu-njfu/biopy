#!/usr/bin/env python
"""
File: plot_volcano.py
Description: Plot the volcano of differentially expressed genes
CreateDate: 2022/2/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, List, Union
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(gene_exp_file: str,
         log2fc: float,
         padj: float,
         figure_size: Union[Tuple, List],
         output_file: str):
    # set figure attr
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    plt.figure(figsize=figure_size)
    # read in data
    df = pd.read_table(gene_exp_file, index_col='gene_id')
    df['-lg_p'] = -np.log10(df['padj'])  # calculate -log10(padj) for each gene
    # mark significance
    df['sig'] = 'normal'
    df.loc[(df.log2FoldChange > log2fc) & (df.padj < padj), 'sig'] = 'up'
    df.loc[(df.log2FoldChange < -log2fc) & (df.padj < padj), 'sig'] = 'down'
    # plot scatter figure
    plt.scatter(df[df['sig'] == 'up']['log2FoldChange'],
                df[df['sig'] == 'up']['-lg_p'],
                color='salmon', marker='o', label='up', s=3, zorder=1)
    plt.scatter(df[df['sig'] == 'normal']['log2FoldChange'],
                df[df['sig'] == 'normal']['-lg_p'],
                color='lightgrey', marker='o', label='normal', s=3, zorder=0)
    plt.scatter(df[df['sig'] == 'down']['log2FoldChange'],
                df[df['sig'] == 'down']['-lg_p'],
                color='lightgreen', marker='o', label='down', s=3, zorder=1)
    # plot line figure
    max_y = float(df['-lg_p'].replace([np.inf, -np.inf], np.nan).max()) + 5
    max_x = float(df['log2FoldChange'].max()) + 1
    min_x = float(df['log2FoldChange'].min()) - 1
    max_abs_x = max(abs(max_x), abs(min_x))
    plt.plot([log2fc, log2fc], [0, max_y], '--', alpha=0.8, color='k', linewidth=0.8, zorder=2)
    plt.plot([-log2fc, -log2fc], [0, max_y], '--', alpha=0.8, color='k', linewidth=0.8, zorder=2)
    plt.plot([-max_abs_x, max_abs_x], [-np.log10(padj), -np.log10(padj)], '--', alpha=0.8, color='k', linewidth=0.8, zorder=2)
    # set axes label
    plt.xlabel(r'$log_2Fold$' + ' Change')
    plt.ylabel(r'$-log_{10}P$' + ' adjust')
    # set axes lim
    plt.xlim([-max_abs_x, max_abs_x])
    # show legend and save figure
    plt.legend()
    plt.savefig(output_file, bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file',
              metavar='<file|stdin>', required=True, type=click.File('r'),
              help="Expression matrix file. (header must including 'gene_id', 'padj' and 'log2FoldChange')")
@click.option('-l', '--log2fc', 'log2fc', type=float, default=1.5, show_default=True,
              help='Set log2(fold change) as the threshold.')
@click.option('-p', '--p_adjust', 'p_adjust',
              metavar='<float>', type=float, default=0.05, show_default=True,
              help='Set padj as the threshold.')
@click.option('-s', '--figure_size', 'figure_size', default='6.0x8.0', show_default=True, help='Specify figure size.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', default='volcano.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gene_exp_file, log2fc, p_adjust, figure_size, output_file):
    """Plot the volcano of differentially expressed genes."""
    figure_size = [float(i) for i in figure_size.split('x')]
    main(
        gene_exp_file=gene_exp_file,
        log2fc=log2fc,
        padj=p_adjust,
        figure_size=figure_size,
        output_file=output_file,
    )


if __name__ == '__main__':
    run()
