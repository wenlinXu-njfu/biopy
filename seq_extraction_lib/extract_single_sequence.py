#!/usr/bin/env python
"""
File: extract_single_sequence.py
Description: Extract one sub-sequence from reference sequence file
Date: 2022/5/3
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.fasta import Fasta
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(fasta_file, chr_num, start, end, strand, out_file):
    for nucl in Fasta(fasta_file).parse():
        if nucl.id == chr_num:
            sub_seq = nucl[start - 1:end]
            if strand == '-':
                sub_seq = sub_seq.get_reverse_complementary_seq()
            sub_seq.id = f'{chr_num}:{start}-{end}({strand})'
            if out_file:
                with open(out_file, 'w') as o:
                    o.write(f">{sub_seq.id}\n{sub_seq.seq}\n")
            else:
                print(sub_seq)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--ref_fasta', 'ref_fasta_file', help='Input reference sequence file. (format: FASTA)')
@click.option('-c', '--chr_name', 'chr_name', help='Specify chromosome name. (eg:Chr01)')
@click.option('-s', '--start_site', 'start_site', type=int, help='Specify start site on chromosome.')
@click.option('-e', '--end_site', 'end_site', type=int, help='Specify end site on chromosome.')
@click.option('-S', '--strand', 'strand', type=click.Choice(['+', '-']), help='Specify the direction of the chain.')
@click.option('-o', '--output_file', 'output_file',
              help='[optional] Output file, if not specified, print result to terminal as stdout.')
def run(ref_fasta_file, chr_name, start_site, end_site, strand, output_file):
    """Extract one sub sequence from reference sequence file."""
    main(ref_fasta_file, chr_name, start_site, end_site, strand, output_file)


if __name__ == '__main__':
    run()
