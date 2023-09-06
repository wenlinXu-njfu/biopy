#!/usr/bin/env python
"""
File: ORF_finder.py
Description: ORF prediction.
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os.path import abspath
from tqdm import tqdm
from multiprocessing import Pool
from io import TextIOWrapper
from typing import Tuple, Iterable
import click
from gzip import GzipFile
from Biolib import Fasta, Nucleotide, Displayer
displayer = Displayer(__file__.split('/')[-1], version='1.1.0')


def sub_processing(nucl_obj: Nucleotide, min_len: int, complete: bool, only_plus: bool):
    ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
    return ORF


def write_pep_to_file(content: str, output_file: str):
    with open(output_file, 'a') as o:
        o.write(content)


def show_process_bar(fasta_file: TextIOWrapper,
                     params: Iterable[Tuple[Nucleotide, int, bool, bool]],
                     log_file: TextIOWrapper,
                     processes_num: int):
    pool = Pool(processes=processes_num)
    results = []
    try:
        seq_num = sum(1 for _ in fasta_file if _.startswith('>'))
    except UnicodeDecodeError:
        seq_num = sum(1 for _ in GzipFile(fasta_file.name) if str(_, 'utf8').startswith('>'))
    with tqdm(total=seq_num, unit=' sequence',
              desc=f'\033[36m[Processing {abspath(fasta_file.name)}]\033[0m') as pbar:
        for param in params:
            ORF = pool.apply_async(sub_processing, args=param, callback=lambda _: pbar.update(1))
            results.append(ORF)
        pool.close()
        pool.join()
    results = [i.get() for i in results]
    content = [f'>{ORF.id}\n{ORF.seq}\n' for ORF in results if not isinstance(ORF, str)]
    log = '\n'.join([i for i in results if isinstance(i, str)]) + '\n'
    output_prefix = fasta_file.name.replace('.gz', '').split('/')[-1]
    write_pep_to_file(''.join(content), f'./{output_prefix}.orf')
    log_file.write(log)


def main(fasta_files: Tuple[TextIOWrapper],
         parse_seqids: bool,
         min_len: int,
         complete: bool,
         only_plus: bool,
         log_file: TextIOWrapper = None,
         to_file: bool = None,
         processes_num: int = 1):
    for fasta_file in fasta_files:
        params = ((nucl_obj, min_len, complete, only_plus) for nucl_obj in Fasta(fasta_file).parse(parse_seqids))
        # Show the progress bar on the command line.
        if log_file and to_file and fasta_file.name != '<stdin>':
            show_process_bar(fasta_file, params, log_file, processes_num)
        # Do not show the progress bar on the command line.
        else:
            pool = Pool(processes=processes_num)
            results = []
            for param in params:
                ORF = pool.apply_async(sub_processing, args=param)
                results.append(ORF)
            pool.close()
            pool.join()
            results = [i.get() for i in results]
            content = [f'>{ORF.id}\n{ORF.seq}\n' for ORF in results if not isinstance(ORF, str)]
            log = '\n'.join([i for i in results if isinstance(i, str)]) + '\n'
            log_file.write(log) if log_file else click.echo(f'\033[33m{log}\033[0m', err=True)
            output_prefix = fasta_file.name.replace('.gz', '').split('/')[-1] if fasta_file.name != '<stdin>' else 'stdin'
            if to_file:
                with open(f'./{output_prefix}.orf', 'w') as o:
                    o.write(''.join(content))
            else:
                click.echo(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_files', nargs=-1, type=click.File('r'), required=True)
@click.option('-l', '--min_len', 'min_len', type=int, default=30, show_default=True, help='Minimal ORF length.')
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence id.')
@click.option('-c', '--completed', 'completed', is_flag=True, flag_value=True, help='Remain completed ORF.')
@click.option('-p', '--only_plus', 'only_plus', is_flag=True, flag_value=True, help='Only predict plus chain.')
@click.option('-log', '--log_file', 'log_file', type=click.File('a'),
              help='Write the sequence that not found ORF to logfile.')
@click.option('-o', '--to_file', 'to_file', is_flag=True, flag_value=True, show_default=True,
              help='Write the results to file rather than print to terminal as stdout.')
@click.option('-n', '--processes_num', 'processes_num', type=int, default=1, show_default=True, help='Number of processes.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_files, parse_seqids, min_len, completed, only_plus, log_file, to_file, processes_num):
    """ORF prediction."""
    main(fasta_files, parse_seqids, min_len, completed, only_plus, log_file, to_file, processes_num)


if __name__ == '__main__':
    run()
