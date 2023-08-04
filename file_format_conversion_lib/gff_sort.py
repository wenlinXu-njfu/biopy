#!/usr/bin/env python
"""
File: gff_sort.py
Description: Sort the GFF file by sequence ID
Date: 2022/3/31
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gff import Gff
from Biolib.show_info import Displayer


def main(gff_file, out_file):
    gff_file_obj = Gff(gff_file)
    content = gff_file_obj.gff_sort()
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gff_file', 'gff_file', help='Input unsorted GFF file.')
@click.option('-o', '--output_file', 'output_file',
              help='Output sorted GFF file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gff_file, output_file):
    """Sort the GFF file by sequence ID."""
    main(gff_file, output_file)


if __name__ == '__main__':
    run()
