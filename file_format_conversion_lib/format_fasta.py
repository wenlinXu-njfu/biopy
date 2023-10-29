#!/usr/bin/env python
"""
File: format_fasta.py
Description: Make each sequence to be displayed in a single line or in multiple lines and sort sequence by length or ID.
CreateDate: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from re import findall
from natsort import natsort_key
import click
from pybioinformatic import Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fasta_file: Union[str, TextIOWrapper],
         parse_seqids: bool = True,
         char_num: int = 60,
         sort_by_len: bool = True,
         sort_by_id: bool = True,
         out_file: TextIOWrapper = None):
    with Fasta(fasta_file) as fa:
        seq_obj_generator = fa.merge_sequence(parse_seqids) if char_num == 0 else \
            fa.split_sequence(parse_seqids, char_num)
        seq_obj_list = [seq_obj for seq_obj in seq_obj_generator]
        if sort_by_len and not sort_by_id:
            seq_obj_list.sort(key=lambda i: i.len)
        elif not sort_by_len and sort_by_id:
            seq_obj_list.sort(key=lambda seq_obj: natsort_key(seq_obj.id))
        elif sort_by_len and sort_by_id:
            seq_obj_list.sort(key=lambda i: (i.len, findall(r'[a-zA-Z]+', i.id)[0], int(findall(r'\d+', i.id)[0])))
        for seq_obj in seq_obj_list:
            click.echo(seq_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input FASTA file.')
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-n', '--char_num', 'char_num',
              metavar='<int>', type=int, default=60, show_default=True,
              help='Specify how many character show in per line, 0 presents one line show per sequence.')
@click.option('-L', '--sort_by_len', 'sort_by_len',
              is_flag=True, flag_value=True,
              help='Sort sequence by length. If both "-L, --sort_by_len" and "-I, --sort_by_id" options are specified, '
                   'sort by length first then id by default.')
@click.option('-I', '--sort_by_id', 'sort_by_id',
              is_flag=True, flag_value=True,
              help='Sort sequence by id. If both "-L, --sort_by_len" and "-I, --sort_by_id" options are specified, '
                   'sort by length first then id by default.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<fasta file>', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, char_num, sort_by_len, sort_by_id, outfile):
    """Make each sequence to be displayed in a single line or in multiple lines and sort sequence by length or ID."""
    main(fasta_file, parse_seqids, char_num, sort_by_len, sort_by_id, outfile)


if __name__ == '__main__':
    run()
