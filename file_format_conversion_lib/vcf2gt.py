#!/usr/bin/env python
"""
File: vcf2gt.py
Description: Convert the file format from VCF to Genotype.
CreateDate: 2023/10/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import VCF, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(vcf_file: Union[str, TextIOWrapper],
         num_processing: int,
         output_file: Union[str, TextIOWrapper]):
    with VCF(vcf_file) as vcf:
        gt = vcf.to_genotype(num_processing)
        gt.to_csv(output_file, sep='\t')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--vcf_file', 'vcf_file',
              metavar='<vcf file>', type=click.File('r'), required=True,
              help='Input VCF file.')
@click.option('-n', '--num_processing', 'num_processing',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file>', type=click.File('w'), default='GT.xls', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Converting the vcf file to gt format.')
def run(vcf_file, num_processing, output_file):
    """Convert the file format from VCF to Genotype."""
    main(vcf_file, num_processing, output_file)


if __name__ == '__main__':
    run()
