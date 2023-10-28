#!/usr/bin/env python
"""
File: get_seq_len.py
Description: Get each sequence length of FASTA file.
CreateDate: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Tuple
from os.path import abspath
import click
from pybioinformatic import Fasta, Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(fasta_files: Tuple[TextIOWrapper], parse_seqids: bool, out_file: TextIOWrapper):
    for fasta_file in fasta_files:
        for nucl_obj in Fasta(fasta_file).parse(parse_seqids):
            click.echo(f'{nucl_obj.get_seq_len_info()}\t{abspath(fasta_file.name)}', out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, metavar='<fasta files>', type=click.File('r'), required=True)
@click.option('-P', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file>', type=click.File('w'),
              help='Output file (Seq_id\\tSeq_len\\n), if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_files, parse_seqids, outfile):
    """Get each sequence length of FASTA file."""
    main(fasta_files, parse_seqids, outfile)


if __name__ == '__main__':
    run()
