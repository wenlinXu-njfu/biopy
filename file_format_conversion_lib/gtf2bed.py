#!/usr/bin/env python
"""
File: gtf2bed.py
Description: Convert the file format from GTF to BED
Date: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from _io import TextIOWrapper
import click
from Biolib.gtf import Gtf
from Biolib.show_info import Displayer


def main(gtf_file: TextIOWrapper, out_file: TextIOWrapper):
    content = []
    for line in Gtf(gtf_file).gtf_to_bed():
        if out_file:
            content.append(line)
        else:
            print(line.strip())
    if out_file:
        with out_file as o:
            o.write(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gtf_file', 'gtf_file', type=click.File('r'), required=True, help='Input GTF file.')
@click.option('-o', '--output_bed_file', 'output_bed_file', type=click.File('w'),
              help='Output BED file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(gtf_file, output_bed_file):
    """Convert the file format from GTF to BED."""
    main(gtf_file, output_bed_file)


if __name__ == '__main__':
    run()
