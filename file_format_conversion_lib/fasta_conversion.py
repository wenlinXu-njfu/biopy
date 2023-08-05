#!/usr/bin/env python
"""
File: fasta_conversion.py
Description: Make each sequence to be displayed on a single line or in multiple lines.
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
import click
from click._compat import _NonClosingTextIOWrapper
from click.utils import KeepOpenFile
from Biolib.fasta import Fasta
from Biolib.show_info import Displayer


def main(fasta_file: Union[str, KeepOpenFile],
         char_num: int,
         out_file: str):
    content = []
    if char_num == 0:
        for seq in Fasta(fasta_file).merge_sequence():
            if not out_file:
                print(seq.strip())
            else:
                content.append(seq)
    else:
        for line in Fasta(fasta_file).split_sequence(char_num):
            if out_file:
                content.append(line)
            else:
                print(line.strip())
    if out_file and content:
        with open(out_file, 'w') as o:
            o.write(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', type=click.File('r'), help='Input FASTA file.')
@click.option('-n', '--char_num', 'char_num', type=int, default=60, show_default=True,
              help='Specify how many character show in per line, 0 presents one line show per sequence.')
@click.option('-o', '--output_file', 'outfile',
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_file, char_num, outfile):
    """Make each sequence to be displayed on a single line or in multiple lines."""
    if isinstance(fasta_file, _NonClosingTextIOWrapper):
        fasta_file = click.open_file('-')
    main(fasta_file, char_num, outfile)


if __name__ == '__main__':
    run()
