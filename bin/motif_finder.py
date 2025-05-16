#!/usr/bin/env python
"""
File: motif_finder.py
Description: Find the motif in the sequence.
CreateDate: 2022/3/29
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Tuple, Union
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.3.0')


def main(fasta_files: Tuple[Union[str, TextIOWrapper]],
         motif: str,
         only_forward: bool = False,
         quiet: bool = False,
         log_file: TextIOWrapper = None,
         output_file: TextIOWrapper = None):
    click.echo('Seq_id\tStart\tEnd\tStrand\tMotif\tSequence', output_file)
    for fasta_file in fasta_files:
        with Fasta(fasta_file) as fa:
            for seq_obj in fa.parse():
                ret = seq_obj.find_motif(motif=motif, only_forward=only_forward)
                if 'not found' not in ret:
                    click.echo(ret, output_file)
                else:
                    if not quiet:
                        click.echo(f"\033[33m{ret}\033[0m", err=True, file=log_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, metavar='<fasta files|stdin>', type=click.File('r'), required=True)
@click.option('-m', '--motif', 'motif',
              metavar='<str>', required=True,
              help='Specify motif sequence, support for regular expressions.')
@click.option('-F', '--only-forward', 'only_forward', is_flag=True, flag_value=True,
              help='Only search forward chain.')
@click.option('-q', '--quiet', 'quiet',
              is_flag=True, flag_value=True,
              help='Do not report sequence that not found motif. This conflicts with the "-l --log_file" option and '
                   'takes precedence over the "-l --log_file" option.')
@click.option('-log', '--log_file', 'log_file',
              metavar='<file|stderr>', type=click.File('w'),
              help='Write the sequence that not found motif to logfile, stderr by default. '
                   'This conflicts with the "-q --quiet" option and has a lower priority than the "-q --quiet" option.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file|stdout>', type=click.File('w'),
              help=r'Output file (Seq_id\tStart\tEnd\tStrand\tMotif\tSequence), stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_files, motif, only_forward, quiet, log_file, outfile):
    """Find the motif in the sequence."""
    main(
        fasta_files=fasta_files,
        motif=motif,
        only_forward=only_forward,
        quiet=quiet,
        log_file=log_file,
        output_file=outfile
    )


if __name__ == '__main__':
    run()
