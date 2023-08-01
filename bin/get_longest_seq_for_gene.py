#!/usr/bin/env python
"""
File: get_longest_seq_for_gene.py
Description: Get the longest transcript of each gene
Date: 2022/3/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.fasta import Fasta
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(fasta_file, regula_exp, inplace_id: click.Choice(['yes', 'no']), out_file):
    fa_file_obj = Fasta(fasta_file)
    d = {'yes': True, 'no': False}
    inplace_id = d[inplace_id]
    seq_dict = fa_file_obj.get_longest_seq(regula_exp, inplace_id)
    content = []
    for seq_id, seq in seq_dict.items():
        content.append(f">{seq_id}\n{seq}")
    content = '\n'.join(content)
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--fasta_file', 'fasta_file', help='Input FASTA file.')
@click.option('-r', '--regular_expression', 'regular_expression',
              help='The name of a gene locus represented by a regular expression.')
@click.option('-I', '--inplace_id', 'inplace_id', type=click.Choice(['yes', 'no']), default='no',
              help='[optional] Replace the longest sequence ID with unique ID {default=no}')
@click.option('-o', '--output_file', 'outfile',
              help='[optional] Output file, if not specified, print results to terminal as stdout.')
def run(fasta_file, regular_expression, inplace_id, outfile):
    """Get the longest transcript of each gene."""
    main(fasta_file, regular_expression, inplace_id, outfile)


if __name__ == '__main__':
    run()
