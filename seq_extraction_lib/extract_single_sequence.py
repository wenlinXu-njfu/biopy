#!/usr/bin/env python
"""
File: extract_single_sequence.py
Description: Extract one sub-sequence from reference sequence file.
CreateDate: 2022/5/3
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.1')


def main(fasta_file: TextIOWrapper,
         chr_num: str,
         start: int,
         end: int,
         strand: click.Choice(['+', '-']),
         output_file: TextIOWrapper = None):
    with Fasta(fasta_file) as fa:
        for nucl in fa.parse():
            if nucl.id == chr_num:
                sub_seq = nucl[start - 1:end]
                if strand == '-':
                    sub_seq = sub_seq.get_reverse_complementary_seq()
                if start > len(nucl):
                    click.echo(f'\033[31mError: The interval "{chr_num}:{start}-{end}" is out of {chr_num} sequence range.\033[0m',
                               err=True)
                    exit()
                elif start < len(nucl) < end:
                    end = len(nucl)
                sub_seq.id = f'{chr_num}:{start}-{end}({strand})'
                click.echo(sub_seq, output_file)
                break


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-r', '--ref_fasta', 'ref_fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input reference sequenceFASTA file.')
@click.option('-c', '--chr_name', 'chr_name',
              metavar='<str>', required=True,
              help='Chromosome name. (eg:Chr01)')
@click.option('-s', '--start_site', 'start_site',
              metavar='<int>', type=int, required=True,
              help='Start site on chromosome. (based on 1)')
@click.option('-e', '--end_site', 'end_site',
              metavar='<int>', type=int, required=True,
              help='End site on chromosome. (based on 1)')
@click.option('-S', '--strand', 'strand',
              metavar='<+|->', type=click.Choice(['+', '-']), required=True,
              help='Direction of the chain.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file>', type=click.File('w'),
              help='Output file, if not specified, print result to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(ref_fasta_file, chr_name, start_site, end_site, strand, output_file):
    """Extract one sub sequence from reference sequence file."""
    main(ref_fasta_file, chr_name, start_site, end_site, strand, output_file)


if __name__ == '__main__':
    run()
