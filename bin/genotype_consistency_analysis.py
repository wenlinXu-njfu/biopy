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
from os import mkdir
from os.path import exists
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gt_file1: Union[str, TextIOWrapper],
         gt_file2: Union[str, TextIOWrapper],
         database_compare: bool,
         output_path: str):
    if not exists(output_path):
        mkdir(output_path)
    with GenoType(gt_file1) as gt1:
        with GenoType(gt_file2) as gt2:
            if database_compare:
                gt1.compare(gt2, output_path=output_path)
            else:
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
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default='./', show_default=True,
              help='Output file path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Genotype consistency analysis underway.')
def run(database_gt, test_gt, database_compare, output_path):
    """Genotype consistency analysis."""
    main(database_gt, test_gt, database_compare, output_path)


if __name__ == '__main__':
    run()
