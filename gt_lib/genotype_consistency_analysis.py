#!/usr/bin/env python
"""
File: genotype_consistency_analysis.py
Description: Genotype consistency analysis.
CreateDate: 2023/10/20
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from os import makedirs, getcwd
import matplotlib.pyplot as plt
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.1')


@Timer('Genotype consistency analysis underway.')
def main(gt_file1: Union[str, TextIOWrapper],
         gt_file2: Union[str, TextIOWrapper],
         database_compare: bool,
         color_map: str,
         reverse_cmap: bool,
         output_path: str):
    makedirs(output_path, exist_ok=True)
    if reverse_cmap:
        color_map = plt.get_cmap(color_map).reversed()
    with GenoType(gt_file1) as gt1, GenoType(gt_file2) as gt2:
        gt1.compare(gt2, output_path=output_path, cmap=color_map) \
            if database_compare else \
            gt1.self_compare(gt2, output_path=output_path)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--database-gt', 'database_gt',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input database sample gt file.')
@click.option('-I', '--test-gt', 'test_gt',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input test sample gt file.')
@click.option('--database-compare/--self-compare', 'database_compare',
              default=True, show_default=True,
              help='Specify comparison mode. Database compare means test sample compare with database sample. '
                   'Self compare means each sample compares with its different test batch.')
@click.option('-c', '--color-map', 'color_map',
              metavar='<str>', default='RdYlGn_r', show_default=True,
              help='Color map of colorbar. See url "https://matplotlib.org/stable/users/explain/colors/colormaps.html".')
@click.option('-r', '--reverse-cmap', 'reverse_cmap',
              is_flag=True, flag_value=True,
              help='Reverse color map.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output file path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(database_gt, test_gt, database_compare, color_map, reverse_cmap, output_path):
    """Genotype consistency analysis."""
    main(database_gt, test_gt, database_compare, color_map, reverse_cmap, output_path)


if __name__ == '__main__':
    run()
