#!/usr/bin/env python
"""
File: gtf2bed.py
Description: Convert the file format from GTF to BED.
CreateDate: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Gtf, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gtf_file: Union[str, TextIOWrapper],
         out_file: TextIOWrapper = None):
    with Gtf(gtf_file) as gtf:
        for line in gtf.to_bed():
            click.echo(line, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gtf_file', 'gtf_file',
              metavar='<gtf file|stdin>', type=click.File('r'), required=True,
              help='Input GTF file.')
@click.option('-o', '--output_bed_file', 'output_bed_file',
              metavar='<bed file|stdout>', type=click.File('w'),
              help='Output BED file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gtf_file, output_bed_file):
    """Convert the file format from GTF to BED."""
    main(gtf_file, output_bed_file)


if __name__ == '__main__':
    run()
