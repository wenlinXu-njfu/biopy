#!/usr/bin/env python
"""
File: ORF_finder.py
Description: ORF prediction
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from tqdm import tqdm
import click
from gzip import GzipFile
from Biolib.fasta import Fasta
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(fa_file, parse_seqids: click.Choice(['yes', 'no']), min_len, complete,
         only_plus: click.Choice(['yes', 'no']), log, out_file):
    only_plus = {'yes': True, 'no': False}[only_plus]
    parse_seqids = {'yes': True, 'no': False}[parse_seqids]
    content = ''
    if out_file and log:
        try:
            seq_num = sum(1 for _ in open(fa_file) if _.startswith('>'))
        except UnicodeDecodeError:
            seq_num = sum(1 for _ in GzipFile(fa_file) if str(_, 'utf8').startswith('>'))
        with tqdm(total=seq_num, unit=' sequence') as process_bar:
            for nucl_obj in Fasta(fa_file).parse(parse_seqids):
                ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
                if isinstance(ORF, str):
                    click.echo(ORF, err=True, file=open(log, 'a')) if log else \
                        click.echo(f"\033[33m{ORF}\033[0m", err=True)
                else:
                    if out_file:
                        content += f">{ORF.id}\n{ORF.seq}\n"
                    else:
                        print(ORF)
                process_bar.update(1)
    else:
        for nucl_obj in Fasta(fa_file).parse(parse_seqids):
            ORF = nucl_obj.ORF_prediction(min_len, complete, only_plus)
            if isinstance(ORF, str):
                click.echo(ORF, err=True, file=open(log, 'a')) if log else \
                    click.echo(f"\033[33m{ORF}\033[0m", err=True)
            else:
                if out_file:
                    content += f">{ORF.id}\n{ORF.seq}\n"
                else:
                    print(ORF)
    if out_file and content:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--fasta_file', 'fasta_file', help='Input nucleotide sequence file. (format=FASTA)')
@click.option('-p', '-parse_seqids', 'parse_seqids', type=click.Choice(['yes', 'no']), default='yes',
              help='[optional] Specify whether parse sequences IDs. {default=yes}')
@click.option('-L', '--min_len', 'min_len', type=int, default=30, help='[optional] Minimal ORF length. {default=30nt}')
@click.option('-complete/-incomplete', default=True,
              help='[optional] Whether ORF integrity is considered. {default=complete}')
@click.option('-P', '-only_plus', 'only_plus', type=click.Choice(['yes', 'no']), default='no',
              help='[optional] Whether only consider plus chain. {default=no}')
@click.option('-l', '--log_file', 'log_file',
              help='[optional] Output log file, if not specified, the log will print to terminal as stderr.')
@click.option('-o', '--output_file', 'outfile',
              help='[optional] Output peptide chain file, if not specified, print results to terminal as stdout.')
def run(fasta_file, parse_seqids, min_len, complete, only_plus, log_file, outfile):
    """ORF prediction."""
    main(fasta_file, parse_seqids, min_len, complete, only_plus, log_file, outfile)


if __name__ == '__main__':
    run()
