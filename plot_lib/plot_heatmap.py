#!/usr/bin/env python
"""
File: plot_heatmap.py
Description: Plot gene expression heatmap
CreateDate: 2022/4/15
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from numpy import log2
from matplotlib.pyplot import rcParams, savefig
from seaborn import heatmap
import click
from pybioinformatic import read_in_gene_expression_as_dataframe, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gene_exp_file: str, color_map: str, out_file: str):
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'
    df = read_in_gene_expression_as_dataframe(gene_exp_file)
    df = log2(df + 1)
    df.index.name = None
    heatmap(df, cmap=color_map, yticklabels=False)
    savefig(out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gene_exp', 'gene_exp_file',
              metavar='<exp file>', required=True,
              help='Input gene expression file. (support format: txt, xls, xlsx, and csv)')
@click.option('-c', '--color_map', 'color_map',
              metavar='<str>', default='PiYG', show_default=True,
              help='Color map (eg. PiYG, vlag, OrRd, YlOrRd, Spectral).')
@click.option('-o', '--output_file', 'output_file',
              metavar='<figure file>', default='heatmap.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gene_exp_file, color_map, output_file):
    """Plot gene expression heatmap."""
    main(gene_exp_file, color_map, output_file)


if __name__ == '__main__':
    run()
