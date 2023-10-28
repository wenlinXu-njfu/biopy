#!/usr/bin/env python
"""
File: fa2tab.py
Description: 
CreateDate: 2023/10/14
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_file: TextIOWrapper,
         parse_seqids: bool,
         output_file: TextIOWrapper):
    fa = Fasta(fasta_file)
    for ret in fa.fa2tab(parse_seqids):
        click.echo(ret, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input file.')
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file name>', type=click.File('w'),
              help='Output tab delimited text file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, output_file):
    """Description."""
    main(fasta_file, parse_seqids, output_file)


if __name__ == '__main__':
    run()
