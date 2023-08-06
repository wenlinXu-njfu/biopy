#!/usr/bin/env python
"""
File: gff2gtf.py
Description: Convert the file format from GFF to GTF
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gff import Gff
from Biolib.show_info import Displayer


def main(gff_file, gtf_file=None):
    gff_file_obj = Gff(gff_file)
    content = gff_file_obj.gff_to_gtf()
    if gtf_file:
        with gtf_file as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gff_file', 'gff_file', type=click.File('r'), help='Input GFF file.')
@click.option('-o', '--gtf_file', 'gtf_file', type=click.File('w'),
              help='Output GTF file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gff_file, gtf_file):
    """Convert the file format from GFF to GTF."""
    main(gff_file, gtf_file)


if __name__ == '__main__':
    run()
