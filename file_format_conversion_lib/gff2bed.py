#!/usr/bin/env python
"""
File: gff2bed.py
Description: Convert the file format from GFF to BED
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from _io import TextIOWrapper
import click
from Biolib.gff import Gff
from Biolib.show_info import Displayer


def main(gff_file: TextIOWrapper,
         feature_type: str,
         bed_file: TextIOWrapper):
    gff_file_obj = Gff(gff_file)
    if feature_type:
        feature_type = feature_type.split(',')
    content = gff_file_obj.gff_to_bed(feature_type)
    if bed_file:
        with bed_file as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gff_file', 'gff_file', type=click.File('r'), required=True, help='Input GFF file.')
@click.option('-t', '--feature_type', 'feature_type',
              help='Specify the type of feature to convert to BED format, multiple types are separated by commas. [default: all]')
@click.option('-o', '--output_bed_file', 'output_bed_file', type=click.File('w'),
              help='Output BED file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gff_file, feature_type, output_bed_file):
    """Convert the file format from GFF to BED."""
    main(gff_file, feature_type, output_bed_file)


if __name__ == '__main__':
    run()
