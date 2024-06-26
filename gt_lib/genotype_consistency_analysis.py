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
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.1')


def main(gt_file1: Union[str, TextIOWrapper],
         gt_file2: Union[str, TextIOWrapper],
         database_compare: bool,
         font_name: str,
         output_path: str):
    makedirs(output_path, exist_ok=True)
    with GenoType(gt_file1) as gt1, GenoType(gt_file2) as gt2:
        gt1.compare(gt2, output_path=output_path, font_name=font_name) \
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
@click.option('-f', '--font_name', 'font_name',
              metavar='<str>', default='Arial', show_default=True,
              help='Set figure font, if "--database-compare" specified.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output file path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Genotype consistency analysis underway.')
def run(database_gt, test_gt, database_compare, font_name, output_path):
    """Genotype consistency analysis."""
    main(database_gt, test_gt, database_compare, font_name, output_path)


if __name__ == '__main__':
    run()
