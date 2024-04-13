#!/usr/bin/env python
"""
File: get_reverse_complementary_seq.py
Description: Get reverse complementary sequence.
CreateDate: 2022/6/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_file: Union[str, TextIOWrapper],
         parse_seqids: bool = True,
         out_file: TextIOWrapper = None):
    with Fasta(fasta_file) as fa:
        for nucl_obj in fa.parse(parse_seqids):
            rev_com_seq = -nucl_obj
            click.echo(rev_com_seq, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence IDs.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, outfile):
    """Get reverse complementary sequence."""
    main(fasta_file, parse_seqids, outfile)


if __name__ == '__main__':
    run()
