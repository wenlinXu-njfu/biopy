#!/usr/bin/env python
"""
File: genotype_consistency_analysis.py
Description: Genotype consistency analysis.
Date: 2023/10/20
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1])


def main(gt_file1: Union[str, TextIOWrapper],
         gt_file2: Union[str, TextIOWrapper],
         output_prefix: str):
    gt1 = GenoType(gt_file1)
    gt2 = GenoType(gt_file2)
    with open(f'{output_prefix}.consistency.xls', 'w') as o:
        for result in gt1.compare(gt2, output_prefix=output_prefix):
            click.echo(result, o)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--database-gt', 'database_gt',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input database sample gt file.')
@click.option('-I', '--test-gt', 'test_gt',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input test sample gt file.')
@click.option('-o', '--output-prefix', 'output_prefix',
              metavar='<str>', default='TestSample', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Genotype consistency analysis underway.')
def run(database_gt, test_gt, output_prefix):
    """Genotype consistency analysis."""
    main(database_gt, test_gt, output_prefix)


if __name__ == '__main__':
    run()
