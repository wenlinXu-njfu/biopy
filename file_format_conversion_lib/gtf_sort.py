#!/usr/bin/env python
"""
File: gtf_sort.py
Description: Sort naturally in sequence according to chromosome, gene id, transcript id, start, and end.
CreateDate: 2025/4/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Gtf, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gtf_file: TextIOWrapper,
         output_file: TextIOWrapper = None):
    with Gtf(gtf_file) as gtf:
        click.echo(gtf.sort(), output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('gtf_file', nargs=1, metavar='<gtf file|stdin>', type=click.File('r'), required=True)
@click.option('-o', '--output-file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gtf_file, output_file):
    """Sort naturally in sequence according to chromosome, gene id, transcript id, start, and end."""
    main(gtf_file, output_file)


if __name__ == '__main__':
    run()
