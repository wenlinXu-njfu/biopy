#!/usr/bin/env python
"""
File: plot_venn.py
Description: Draw the venn plot
Date: 2022/2/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import listdir
from venn import venn
import matplotlib.pyplot as plt
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(in_dir: str, figure_size: tuple, out_prefix: str, fmt: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    files = listdir(in_dir)
    data_set_dict = {}
    for file in files:
        data_set_dict[file] = set(line.strip() for line in open(f'{in_dir}/{file}') if line.strip())
    venn(data_set_dict, legend_loc=(0.8, 0.8), figsize=figure_size)
    plt.savefig(f'{out_prefix}.{fmt}')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input_dir', 'input_dir', help='Input directory where data files are stored.')
@click.option('-s', '--figure_size', default='8.0x8.0', help='[optional] Specify figure size. {default: 8.0x8.0}')
@click.option('-o', '--output_prefix', 'output_prefix', default='venn',
              help='Prefix of output file. {default: venn}')
@click.option('-O', '--output_format', 'output_format', default='pdf',
              type=click.Choice(['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']),
              help='[optional] Specify the format of output file. {default: pdf}')
def run(input_dir, figure_size, out, output_format):
    """Draw the venn plot."""
    figure_size = tuple(float(i) for i in figure_size.split('x'))
    main(input_dir, figure_size, out, output_format)


if __name__ == '__main__':
    run()
