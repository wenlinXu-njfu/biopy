#!/usr/bin/env python
"""
File: blast2bed.py
Description: Convert blast format 6 result to BED format.
CreateDate: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Blast, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(blast_file: Union[str, TextIOWrapper],
         ref_seq: click.Choice(['query', 'sbject']),
         out_file: TextIOWrapper = None):
    query_is_ref = True if ref_seq == 'query' else False
    with Blast(blast_file) as blast:
        content = blast.to_bed(query_is_ref)
        click.echo(content, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--blast_result', 'blast_result_file',
              metavar='<blast file>', type=click.File('r'), required=True,
              help='Input blast format 6 result file.')
@click.option('-r', '--ref_seq', 'ref_seq',
              metavar='<query|sbject>', type=click.Choice(['query', 'sbject']),
              default='sbject', show_default=True,
              help='Specify query sequence or sbject sequence as the reference sequence.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<bed file>', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(blast_result_file, ref_seq, output_file):
    """Convert blast format 6 result to BED format."""
    main(blast_result_file, ref_seq, output_file)


if __name__ == '__main__':
    run()
