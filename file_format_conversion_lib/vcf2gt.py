#!/usr/bin/env python
"""
File: vcf2gt.py
Description: Convert the file format from VCF to Genotype.
CreateDate: 2023/10/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper, StringIO
from re import sub
from natsort import natsort_key
from pandas import read_table, concat
import click
from pybioinformatic import VCF, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(vcf_files: Union[str, TextIOWrapper],
         output_file: Union[str, TextIOWrapper]):
    dfs = []
    for vcf_file in vcf_files:
        with VCF(vcf_file) as vcf:
            gt = '\n'.join([line for line in vcf.to_genotype()])
            df = read_table(StringIO(gt), index_col=[0, 1, 2, 3])
            dfs.append(df)
    merge = concat(dfs, axis=1)
    merge.reset_index(inplace=True)
    merge.sort_values([merge.columns[1], merge.columns[2]], key=natsort_key, inplace=True)
    merge = merge.to_string(index=False, na_rep='NA').strip()
    merge = sub(r'\n +', '\n', merge)
    merge = sub(r' +', '\t', merge)
    click.echo(merge, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('vcf_files', nargs=-1, metavar='<vcf files|stdin>', type=click.File('r'), required=True)
@click.option('-o', '--output_file', 'output_file',
              metavar='<gt file|stdout>', type=click.File('w'),
              help='Output file, print results to terminal by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Converting the vcf file to gt format.')
def run(vcf_files, output_file):
    """Convert the file format from VCF to Genotype."""
    main(vcf_files, output_file)


if __name__ == '__main__':
    run()
