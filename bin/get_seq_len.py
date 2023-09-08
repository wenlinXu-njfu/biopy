#!/usr/bin/env python
"""
File: get_seq_len.py
Description: Get each sequence length of FASTA file.
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from Biolib import Fasta, Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(fasta_file: TextIOWrapper, parse_seqids: bool, out_file: TextIOWrapper):
    for nucl_obj in Fasta(fasta_file).parse(parse_seqids):
        click.echo(nucl_obj.get_seq_len_info(), out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', type=click.File('r'), required=True, help='Input FASTA file.')
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence id.')
@click.option('-o', '--output_file', 'outfile', type=click.File('w'),
              help='Output file (Seq_id\\tSeq_len\\n), if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, outfile):
    """Get each sequence length of FASTA file."""
    main(fasta_file, parse_seqids, outfile)


if __name__ == '__main__':
    run()
