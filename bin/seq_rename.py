#!/usr/bin/env python
"""
File: seq_rename.py
Description: Rename sequence id.
CreateDate: 2024/12/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_file: TextIOWrapper,
         map_file: TextIOWrapper,
         parse_seqids: bool,
         remain: bool,
         output_file: TextIOWrapper = None):
    map_dict = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in map_file}
    with Fasta(fasta_file) as fa:
        for ret in fa.rename(map_dict=map_dict, parse_seq_id=parse_seqids, remain=remain):
            click.echo(ret, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta', 'fasta_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input fasta file.')
@click.option('-m', '--map', 'map_file',
              metavar='<file|stdin>', required=True, type=click.File('r'),
              help=r'Sequence id map file. (Old_id\tNew_id\tEtc)')
@click.option('-p', '--parse-seqids', 'parse_seqids', is_flag=True, flag_value=True,
              help='Parse sequence id in fasta file.')
@click.option('-r', '--remain', 'remain',
              is_flag=True, flag_value=True,
              help='Remain the sequence which is not included in map file.')
@click.option('-o', '--output', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, map_file, parse_seqids, remain, output_file):
    """Rename sequence id."""
    main(fasta_file, map_file, parse_seqids, remain, output_file)


if __name__ == '__main__':
    run()
