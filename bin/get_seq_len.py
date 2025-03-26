#!/usr/bin/env python
"""
File: get_seq_len.py
Description: Get each sequence length of FASTA file.
CreateDate: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Tuple, Union
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_files: Tuple[Union[str, TextIOWrapper]],
         parse_seqids: bool = True,
         out_file: TextIOWrapper = None):
    for fasta_file in fasta_files:
        with Fasta(fasta_file) as fa:
            for nucl_obj in fa.parse(parse_seqids):
                click.echo(f'{nucl_obj.get_seq_len_info()}\t{fa.name}', out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, metavar='<fasta files|stdin>', type=click.File('r'), required=True)
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file|stdout>', type=click.File('w'),
              help=r'Output file (Seq_id\tSeq_len\tFile_path), stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_files, parse_seqids, outfile):
    """Get each sequence length of FASTA file."""
    main(fasta_files, parse_seqids, outfile)


if __name__ == '__main__':
    run()
