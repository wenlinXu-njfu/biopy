#!/usr/bin/env python
"""
File: gff2bed.py
Description: Convert the file format from GFF to BED.
CreateDate: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Gff, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gff_file: Union[str, TextIOWrapper],
         feature_type: str,
         bed_file: TextIOWrapper = None):
    with Gff(gff_file) as gff:
        if feature_type:
            feature_type = feature_type.split(',')
        for line in gff.to_bed(feature_type):
            click.echo(line, bed_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gff_file', 'gff_file',
              metavar='<gff file|stdin>', type=click.File('r'), required=True,
              help='Input GFF file.')
@click.option('-t', '--feature_type', 'feature_type', metavar='<str>',
              help='Specify the type of feature to convert to BED format, multiple types are separated by commas. [default: all]')
@click.option('-o', '--output_bed_file', 'output_bed_file',
              metavar='<bed file|stdout>', type=click.File('w'),
              help='Output BED file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gff_file, feature_type, output_bed_file):
    """Convert the file format from GFF to BED."""
    main(gff_file, feature_type, output_bed_file)


if __name__ == '__main__':
    run()
