#!/usr/bin/env python
"""
File: draw_kdeplot.py
Description: Plot univariate or bivariate distributions using kernel density estimation.
CreateDate: 2025/5/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from pandas import read_table
import matplotlib.pyplot as plt
from seaborn import kdeplot
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(data_file: TextIOWrapper,
         output_file: str,
         x: str = None,
         y: str = None,
         hue: str = None):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    df = read_table(data_file, sep='\t')
    kdeplot(data=df, x=x, y=y, hue=hue, multiple="stack")
    plt.savefig(output_file, bbox_inches='tight', dpi=300)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('data_file', nargs=1, metavar='<file|stdin>', type=click.File('r'), required=True)
@click.option('-x', '--x-var', 'x', metavar='<str>',
              help='Variables that specify positions on the x axis.')
@click.option('-y', '--y-var', 'y', metavar='<str>',
              help='Variables that specify positions on the y axis.')
@click.option('--hue', 'hue', metavar='<str>',
              help='Semantic variable that is mapped to determine the color of plot elements.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', default='kdeplot.pdf', show_default=True,
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(data_file, x, y, hue, output_file):
    """Plot univariate or bivariate distributions using kernel density estimation."""
    main(
        data_file=data_file,
        x=x,
        y=y,
        hue=hue,
        output_file=output_file
    )


if __name__ == '__main__':
    run()
