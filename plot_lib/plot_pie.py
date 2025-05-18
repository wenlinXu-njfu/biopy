#!/usr/bin/env python
"""
File: plot_pie.py
Description: Draw pie plot.
CreateDate: 2025/5/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from collections import Counter
import matplotlib.pyplot as plt
import click
from pybioinformatic import generate_unique_colors, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def make_autopct(values, total):
    def my_autopct(pct):
        val = int(round(pct * total / 100))
        return f'{pct:.1f}%\n({val})'
    return my_autopct


def main(
    data: TextIOWrapper,
    colors: str = None,
    output_file: TextIOWrapper = None
):
    plt.rcParams['pdf.fonttype'] = 42
    count = Counter([line.strip() for line in data])
    total = sum(count.values())
    if colors is None:
        colors = generate_unique_colors(num_colors=len(count))
    else:
        colors = colors.split(',')
    plt.pie(
        x=count.values(),
        labels=count.keys(),
        colors=colors,
        autopct=make_autopct(count.values(), total),
        shadow=True
    )
    plt.savefig(output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input-file', 'input_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input data file.')
@click.option('-c', '--colors', 'colors', metavar='<str>',
              help='A color sequence separated by commas which the pie chart will cycle (eg. red,yellow,blue). '
                   'If None, a randomly generated color will be used.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', default='pie.pdf', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, colors, output_file):
    """Draw pie plot."""
    main(input_file, colors, output_file)


if __name__ == '__main__':
    run()
