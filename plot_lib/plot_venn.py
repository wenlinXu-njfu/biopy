#!/usr/bin/env python
"""
File: plot_venn.py
Description: Draw the venn plot.
CreateDate: 2022/2/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Tuple
from venn import venn
import matplotlib.pyplot as plt
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(input_files: Tuple[TextIOWrapper],
         group_name: Tuple[str],
         figure_size: Tuple[float, float],
         output_file: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    data_set_dict = {}
    for file, name in zip(input_files, group_name):
        data_set_dict[name] = set(line.strip() for line in file if line.strip())
    venn(data_set_dict, legend_loc=(0.8, 0.8), figsize=figure_size)
    plt.savefig(output_file, bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('input_files', metavar='<files|stdin>', type=click.File('r'), nargs=-1, required=True)
@click.option('-g', '--group-name', 'group_name',
              metavar='<str>', required=True,
              help='Separated by commas group name of different sets. '
                   'The number of names must be the same as the number of input files.')
@click.option('-f', '--figure-size', 'figure_size',
              metavar='<str>', default='6.4x4.8', show_default=True,
              help='Figure size.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', default='venn.pdf', show_default=True,
              help='Output file (support eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, and tiff).')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_files, group_name, figure_size, output_file):
    """Draw the venn plot."""
    figure_size = tuple(float(i) for i in figure_size.split('x'))
    group_names = tuple(group_name.split(','))
    if len(figure_size) != 2:
        click.echo('\033[31mError: Invalid value of "-s, --figure_size" option.\033[0m', err=True)
    elif len(input_files) != len(group_names):
        click.echo('\033[31mError: The number of input files is inconsistent with the number of collection names.\033[0m', err=True)
    else:
        main(input_files=input_files,
             group_name=group_names,
             figure_size=figure_size,
             output_file=output_file)


if __name__ == '__main__':
    run()
