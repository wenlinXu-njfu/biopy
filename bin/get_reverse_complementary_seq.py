#!/usr/bin/env python
"""
File: get_reverse_complementary_seq.py
Description: Get reverse complementary sequence.
Date: 2022/6/8
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from _io import TextIOWrapper
import click
from Biolib.fasta import Fasta
from Biolib.show_info import Displayer


def main(fasta_file: TextIOWrapper, out_file: TextIOWrapper):
    content = ''
    for nucl_obj in Fasta(fasta_file).parse(False):
        rev_com_seq = -nucl_obj
        if out_file:
            content += f">{rev_com_seq.id}\n{rev_com_seq.seq}\n"
        else:
            print(rev_com_seq)
    if out_file:
        with out_file as o:
            o.write(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', type=click.File('r'), required=True, help='Input FASTA file.')
@click.option('-o', '--output_file', 'outfile', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_file, outfile):
    """Get reverse complementary sequence."""
    main(fasta_file, outfile)


if __name__ == '__main__':
    run()
