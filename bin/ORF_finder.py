#!/usr/bin/env python
"""
File: ORF_finder.py
Description: ORF prediction.
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from tqdm import tqdm
from _io import TextIOWrapper
from typing import Tuple
import click
from gzip import GzipFile
from Biolib.fasta import Fasta
from Biolib.show_info import Displayer


def main(fasta_files: Tuple[TextIOWrapper],
         parse_seqids: bool,
         min_len: int,
         complete: bool,
         only_plus: bool,
         quiet: bool = False,
         log_file: str = None,
         to_file: bool = None):
    for fasta_file in fasta_files:
        content = []
        if fasta_file.name != '<stdin>':
            out_prefix = '.'.join(fasta_file.name.split('/')[-1].split('.')[:-1]).replace('.fa', '') + '_'
        else:
            out_prefix = ''
        # Show the progress bar on the command line.
        if (to_file and log_file) or (to_file and quiet):
            try:
                seq_num = sum(1 for _ in open(fasta_file.name) if _.startswith('>'))
            except UnicodeDecodeError:
                seq_num = sum(1 for _ in GzipFile(fasta_file.name) if str(_, 'utf8').startswith('>'))
            with tqdm(total=seq_num, unit=' sequence') as process_bar:
                for nucl_obj in Fasta(fasta_file.name).parse(parse_seqids):
                    ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
                    if isinstance(ORF, str):
                        if not quiet:
                            click.echo(ORF, err=True, file=open(log_file, 'a')) if log_file else \
                                click.echo(f"\033[33m{ORF}\033[0m", err=True)
                    else:
                        if to_file:
                            content.append(f">{ORF.id}\n{ORF.seq}\n")
                        else:
                            print(ORF)
                    process_bar.update(1)
        # Do not show the progress bar on the command line.
        else:
            for nucl_obj in Fasta(fasta_file).parse(parse_seqids):
                ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
                if isinstance(ORF, str):
                    if not quiet:
                        click.echo(ORF, err=True, file=open(log_file, 'a')) if log_file else \
                            click.echo(f"\033[33m{ORF}\033[0m", err=True)
                else:
                    if to_file:
                        content.append(f">{ORF.id}\n{ORF.seq}\n")
                    else:
                        print(ORF)
        if to_file and content:
            with open(f'./{out_prefix}pep.fa', 'w') as o:
                o.write(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, type=click.File('r'), required=True)
@click.option('-l', '--min_len', 'min_len', type=int, default=30, show_default=True, help='Minimal ORF length.')
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence id.')
@click.option('-c', '--completed', 'completed', is_flag=True, flag_value=True, help='Remain completed ORF.')
@click.option('-p', '--only_plus', 'only_plus', is_flag=True, flag_value=True, help='Only predict plus chain.')
@click.option('-q', '--quiet', 'quiet', is_flag=True, flag_value=True,
              help='Do not report nucleotide sequence that not found ORF. '
                   'This conflicts with the "-L --log_file" option and takes precedence over the "-L --log_file" option.')
@click.option('-log', '--log_file', 'log_file',
              help='Write the sequence that not found ORF to logfile. '
                   'This conflicts with the "-q --quiet" option and has a lower priority than the "-q --quiet" option.')
@click.option('-o', '--to_file', 'to_file', is_flag=True, flag_value=True, show_default=True,
              help='Write the results to file rather than print to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_files, parse_seqids, min_len, completed, only_plus, quiet, log_file, to_file):
    """ORF prediction."""
    main(fasta_files, parse_seqids, min_len, completed, only_plus, quiet, log_file, to_file)


if __name__ == '__main__':
    run()
