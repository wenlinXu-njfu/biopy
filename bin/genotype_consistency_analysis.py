#!/usr/bin/env python
"""
File: genotype_consistency_analysis.py
Description: Genotype consistency analysis.
Date: 2023/10/20
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from pybioinformatic import GenoType, Displayer
displayer = Displayer(__file__.split('/')[-1])


def main(gt_file1, gt_file2, output_prefix):
    gt1 = GenoType(gt_file1)
    gt2 = GenoType(gt_file2)
    for result in gt1.compare(gt2):
        click.echo(result, open(f'{output_prefix}.xls', 'a'))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gt-file1', 'gt_file1',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input gt file.')
@click.option('-I', '--gt-file2', 'gt_file2',
              metavar='<GT file>', type=click.File('r'), required=True,
              help='Input another gt file.')
@click.option('-o', '--output-prefix', 'output_prefix',
              metavar='<str>', default='consistency', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gt_file1, gt_file2, output_prefix):
    """Genotype consistency analysis."""
    main(gt_file1, gt_file2, output_prefix)


if __name__ == '__main__':
    run()
