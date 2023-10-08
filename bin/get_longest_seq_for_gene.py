#!/usr/bin/env python
"""
File: get_longest_seq_for_gene.py
Description: Get the longest transcript of each gene locus.
Date: 2022/3/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(fasta_file: TextIOWrapper,
         regula_exp: str,
         inplace_id: bool,
         out_file: TextIOWrapper = None):
    for seq_obj in Fasta(fasta_file).get_longest_seq(regula_exp, inplace_id):
        click.echo(seq_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-r', '--regular_expression', 'regular_expression',
              metavar='<str>', required=True,
              help='The name of a gene locus represented by a regular expression.')
@click.option('-I', '--inplace_id', 'inplace_id',
              is_flag=True, flag_value=True,
              help='Replace the longest sequence ID with unique ID.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file>', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, regular_expression, inplace_id, outfile):
    """Get the longest transcript of each gene locus."""
    main(fasta_file, regular_expression, inplace_id, outfile)


if __name__ == '__main__':
    run()
