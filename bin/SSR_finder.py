#!/usr/bin/env python
"""
File: SSR_finder.py
Description: Find simple sequence repeat (SSR) in the DNA sequences.
Date: 2023/9/9
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from Biolib import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1])


def main(fasta_file: TextIOWrapper, parse_seqids: bool, output_file: TextIOWrapper = None):
    click.echo('# Seq_id\tStart\tEnd\tSSR_unit\tSSR_seq', output_file)
    for nucl in Fasta(fasta_file).parse(parse_seqids):
        resultts = nucl.find_SSR()
        for result in resultts:
            if result:
                click.echo(result, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', type=click.File('r'), required=True, help='Input FASTA file.')
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence IDs.')
@click.option('-o', '--output_file', 'output_file', type=click.File('w'), help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, output_file):
    """Find simple sequence repeat (SSR) in the DNA sequences."""
    main(fasta_file, parse_seqids, output_file)


if __name__ == '__main__':
    run()
