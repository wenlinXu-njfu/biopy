#!/usr/bin/env python
"""
File: SSR_finder.py
Description: Find simple sequence repeat (SSR) in the DNA sequences.
Date: 2023/9/9
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union
import click
from Biolib import Fasta, Nucleotide, TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1])


def sub_processing(nucl: Nucleotide):
    return nucl.find_SSR()


def main(fasta_file: Union[str, TextIOWrapper],
         parse_seqids: bool,
         quiet: bool,
         processing_num: int,
         output_file: TextIOWrapper = None):
    nucl_obj_generator = ((nucl,) for nucl in Fasta(fasta_file).parse(parse_seqids))
    tkm = TaskManager(processing_num=processing_num, params=nucl_obj_generator)
    results = tkm.parallel_run_func(sub_processing)
    results = [result.get() for result in results]
    click.echo('# Seq_id\tStart\tEnd\tSSR_unit\tSSR_seq', output_file)
    for result in results:
        result = result.strip()
        if quiet and 'not found' not in result:
            click.echo(result, output_file)
        elif not quiet:
            click.echo(result, output_file) if 'not found' not in result else click.echo(result, err=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-P', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence IDs.')
@click.option('-q', '--quiet', 'quiet',
              is_flag=True, flag_value=True,
              help='Do not report sequence that not found SSR motif.')
@click.option('-n', '--processing_num', 'processing_num',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file>', type=click.File('w'),
              help='Output file path and name.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, quiet, processing_num, output_file):
    """Find simple sequence repeat (SSR) in the DNA sequences."""
    main(fasta_file, parse_seqids, quiet, processing_num, output_file)


if __name__ == '__main__':
    run()
