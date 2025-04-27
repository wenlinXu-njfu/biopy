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
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(fasta_file: Union[str, TextIOWrapper],
         parse_seqids: bool = True,
         char_num: int = 60,
         sort_by_len: bool = True,
         sort_by_id: bool = True,
         reverse_id: bool = False,
         reverse_len: bool = False,
         out_file: TextIOWrapper = None):
    with Fasta(fasta_file) as fa:
        seq_obj_generator = fa.merge_sequence(parse_seqids) if char_num == 0 else \
            fa.split_sequence(parse_seqids, char_num)
        seq_obj_list = [seq_obj for seq_obj in seq_obj_generator]
        if sort_by_len and not sort_by_id:
            seq_obj_list.sort(key=lambda i: i.len, reverse=reverse_len)
        elif not sort_by_len and sort_by_id:
            seq_obj_list.sort(key=lambda seq_obj: natsort_key(seq_obj.id), reverse=reverse_id)
        elif sort_by_len and sort_by_id:
            if reverse_len:
                seq_obj_list.sort(
                    key=lambda seq_obj: (-seq_obj.len, natsort_key(seq_obj.id), int(findall(r'\d+', seq_obj.id)[0])),
                    reverse=reverse_id
                )
            else:
                seq_obj_list.sort(
                    key=lambda seq_obj: (seq_obj.len, natsort_key(seq_obj.id), int(findall(r'\d+', seq_obj.id)[0])),
                    reverse=reverse_id
                )
        for seq_obj in seq_obj_list:
            click.echo(seq_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta_file', metavar='<fasta file|stdin>', nargs=1, required=True, type=click.File('r'))
@click.option('-p', '--parse_seqids', 'parse_seqids',
              is_flag=True, flag_value=True,
              help='Parse sequence id.')
@click.option('-n', '--char_num', 'char_num',
              metavar='<int>', type=int, default=60, show_default=True,
              help='Specify how many character show in per line, 0 presents one line show per sequence.')
@click.option('-l', '--sort_by_len', 'sort_by_len',
              is_flag=True, flag_value=True,
              help='Sort sequence by length. If both "-L, --sort_by_len" and "-I, --sort_by_id" options are specified, '
                   'sort by length first then id by default.')
@click.option('-i', '--sort_by_id', 'sort_by_id',
              is_flag=True, flag_value=True,
              help='Sort sequence by id. If both "-L, --sort_by_len" and "-I, --sort_by_id" options are specified, '
                   'sort by length first then id by default.')
@click.option('-r', '--reverse-id', 'reverse_id',
              is_flag=True, flag_value=True,
              help='Reverse the result of sort by id.')
@click.option('-R', '--reverse-len', 'reverse_len',
              is_flag=True, flag_value=True,
              help='Reverse the result of sort by length.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<fasta file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_file, parse_seqids, char_num, sort_by_len, sort_by_id, reverse_id, reverse_len, outfile):
    """Make each sequence to be displayed in a single line or in multiple lines and sort sequence by length or ID."""
    main(
        fasta_file=fasta_file,
        parse_seqids=parse_seqids,
        char_num=char_num,
        sort_by_len=sort_by_len,
        sort_by_id=sort_by_id,
        reverse_id=reverse_id,
        reverse_len=reverse_len,
        out_file=outfile
    )


if __name__ == '__main__':
    run()
