#!/usr/bin/env python
"""
File: ORF_finder.py
Description: Search for open reading frames (ORFs) in the nucleotide sequence you enter, and return maximum length ORF of each nucleotide sequence.
CreateDate: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from tqdm import tqdm
from io import TextIOWrapper
from typing import Tuple, Union
from os import name
from natsort import natsort_key
import click
from pybioinformatic import Fasta, Nucleotide, TaskManager, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def sub_processing(nucl_obj: Nucleotide, min_len: int, complete: bool, only_plus: bool):
    ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
    return ORF


@Timer('ORF predicting.')
def main(fasta_files: Tuple[Union[str, TextIOWrapper]],
         parse_seqids: bool = True,
         min_len: int = 30,
         complete: bool = True,
         only_plus: bool = False,
         log_file: TextIOWrapper = None,
         output_path: str = None,
         num_processes: int = 1):
    tkm = TaskManager(num_processing=num_processes)
    for fasta_file in fasta_files:
        with Fasta(fasta_file) as fa:
            tkm.params = ((nucl_obj, min_len, complete, only_plus) for nucl_obj in fa.parse(parse_seqids))
            # Set output prefix
            if name == 'posix':  # linux
                output_prefix = fa.name.split('/')[-1].replace('.gz', '').replace('<', '').replace('>', '')
            else:  # windows
                output_prefix = fa.name.split('\\')[-1].replace('.gz', '').replace('<', '').replace('>', '')
            output_prefix = '.'.join(output_prefix.split('.')[:-1])
            # Show the progress bar on the command line.
            if (log_file is not None) and (output_path is not None) and ('stdin' not in fa.name):
                with open(f'{output_path}/{output_prefix}_pep.fa', 'w') as output_file:
                    with tqdm(total=fa.seq_num, unit=' sequence', desc=f'[Processing {fa.name}]') as pbar:
                        results = tkm.parallel_run_func(sub_processing, lambda _: pbar.update(1))
                    results = [i.get() for i in results]
                    log = '\n'.join([i for i in results if isinstance(i, str)])
                    click.echo(log, log_file, err=True)
                    ORFs = [ORF for ORF in results if not isinstance(ORF, str)]
                    ORFs.sort(key=lambda i: natsort_key(i.id))
                    for ORF in ORFs:
                        click.echo(ORF, output_file)
            # Do not show the progress bar on the command line.
            else:
                output_file = open(f'{output_path}/{output_prefix}_pep.fa', 'w') if output_path else None
                tkm.parallel_run_func(sub_processing,
                                      lambda ORF: click.echo(ORF, log_file, err=True) if isinstance(ORF, str) else
                                      click.echo(ORF, output_file))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, metavar='<fasta files|stdin>', type=click.File('r'), required=True)
@click.option('-l', '--min-len', 'min_len',
              metavar='<int>', type=int, default=30, show_default=True,
              help='Minimal ORF length.')
@click.option('-p', '--parse-seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-c', '--completed', 'completed',
              is_flag=True, flag_value=True,
              help='Remain completed ORF.')
@click.option('-F', '--forward', 'forward',
              is_flag=True, flag_value=True,
              help='Only predict forward chain.')
@click.option('-log', '--logfile', 'log_file',
              metavar='<file|stderr>', type=click.File('w'),
              help='Write the sequence that not found ORF to logfile, stderr by default.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path|stdout>',
              help='Output path, stdout by default.')
@click.option('-n', '--num_processes', 'num_processes',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processes.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_files, parse_seqids, min_len, completed, forward, log_file, output_path, num_processes):
    """Search for open reading frames (ORFs) in the nucleotide sequence you enter, and return maximum length ORF of each nucleotide sequence."""
    main(fasta_files, parse_seqids, min_len, completed, forward, log_file, output_path, num_processes)


if __name__ == '__main__':
    run()
