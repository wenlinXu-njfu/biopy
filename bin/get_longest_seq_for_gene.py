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
from Biolib.show_info import Displayer


def main(fasta_file, regula_exp, inplace_id: bool, out_file):
    seq_dict = Fasta(fasta_file).get_longest_seq(regula_exp, inplace_id)
    content = []
    for seq_id, seq in seq_dict.items():
        content.append(f">{seq_id}\n{seq}\n")
    content = ''.join(content)
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', help='Input FASTA file.')
@click.option('-r', '--regular_expression', 'regular_expression',
              help='The name of a gene locus represented by a regular expression.')
@click.option('-I', '--inplace_id', 'inplace_id', is_flag=True, flag_value=True,
              help='[optional] Replace the longest sequence ID with unique ID.')
@click.option('-o', '--output_file', 'outfile',
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_file, regular_expression, inplace_id, outfile):
    """Get the longest transcript of each gene."""
    main(fasta_file, regular_expression, inplace_id, outfile)


if __name__ == '__main__':
    run()
