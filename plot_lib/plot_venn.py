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
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(in_dir: str, figure_size: tuple, out_prefix: str, fmt: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    files = listdir(in_dir)
    data_set_dict = {}
    for file in files:
        data_set_dict[file] = set(line.strip() for line in open(f'{in_dir}/{file}') if line.strip())
    venn(data_set_dict, legend_loc=(0.8, 0.8), figsize=figure_size)
    plt.savefig(f'{out_prefix}.{fmt}')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_dir', 'input_dir', help='Input directory where data files are stored.')
@click.option('-s', '--figure_size', default='8.0x8.0', show_default=True, help='Figure size.')
@click.option('-o', '--output_prefix', 'output_prefix', default='venn', show_default=True, help='Prefix of output file.')
@click.option('-O', '--output_format', 'output_format',
              type=click.Choice(['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']),
              default='pdf', show_default=True, help='The format of output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_dir, figure_size, out, output_format):
    """Draw the venn plot."""
    figure_size = tuple(float(i) for i in figure_size.split('x'))
    main(input_dir, figure_size, out, output_format)


if __name__ == '__main__':
    run()
