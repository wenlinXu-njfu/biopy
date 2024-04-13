#!/usr/bin/env python
"""
File: SSR_finder.py
Description: Find simple sequence repeat (SSR) in the DNA sequences.
CreateDate: 2023/9/9
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union
import click
from pybioinformatic import Fasta, Nucleotide, TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def sub_processing(nucl: Nucleotide):
    return nucl.find_SSR()


def main(fasta_file: Union[str, TextIOWrapper],
         parse_seqids: bool,
         quiet: bool,
         num_processing: int,
         output_file: TextIOWrapper = None):
    click.echo('# Seq_id\tStart\tEnd\tSSR_unit\tSSR_seq', output_file)
    with Fasta(fasta_file) as fa:
        params = ((nucl,) for nucl in fa.parse(parse_seqids))
        tkm = TaskManager(num_processing=num_processing, params=params)
        tkm.parallel_run_func(sub_processing,
                              lambda i: click.echo(i.strip(), err=True) if not quiet and 'not found' in i else
                              click.echo(i.strip(), output_file))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence IDs.')
@click.option('-q', '--quiet', 'quiet',
              is_flag=True, flag_value=True,
              help='Do not report sequence that not found SSR motif.')
@click.option('-n', '--num_processing', 'num_processing',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, quiet, num_processing, output_file):
    """Find simple sequence repeat (SSR) in the DNA sequences."""
    main(fasta_file, parse_seqids, quiet, num_processing, output_file)


if __name__ == '__main__':
    run()
