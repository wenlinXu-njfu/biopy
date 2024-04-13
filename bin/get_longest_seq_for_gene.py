#!/usr/bin/env python
"""
File: get_longest_seq_for_gene.py
Description: Get the longest transcript of each gene locus.
CreateDate: 2022/3/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_file: Union[str, TextIOWrapper],
         regula_exp: str,
         inplace_id: bool = False,
         out_file: TextIOWrapper = None):
    with Fasta(fasta_file) as fa:
        for seq_obj in fa.get_longest_seq(regula_exp, inplace_id):
            click.echo(seq_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-r', '--regular_expression', 'regular_expression',
              metavar='<str>', required=True,
              help='The name of a gene locus represented by a regular expression.')
@click.option('-I', '--inplace_id', 'inplace_id',
              is_flag=True, flag_value=True,
              help='Replace the longest sequence ID with unique ID.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, regular_expression, inplace_id, outfile):
    """Get the longest transcript of each gene locus."""
    main(fasta_file, regular_expression, inplace_id, outfile)


if __name__ == '__main__':
    run()
