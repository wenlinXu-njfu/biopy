#!/usr/bin/env python
"""
File: get_seq_len.py
Description: Get each sequence length of FASTA file
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.fasta import Fasta
from Biolib.show_info import Displayer


def main(fa_file, parse_seqids: bool, out_file):
    content = ''
    for nucl_obj in Fasta(fa_file).parse(parse_seqids):
        if out_file:
            content += nucl_obj.get_seq_len_info() + '\n'
        else:
            print(nucl_obj.get_seq_len_info())
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_file', 'fasta_file', help='Input FASTA file.')
@click.option('-P', '--parse_seqids', 'parse_seqids', is_flag=True, flag_value=True, help='Parse sequence id.')
@click.option('-o', '--output_file', 'outfile',
              help='Output file(seq_id\\tseq_len\\n), if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_file, parse_seqids, outfile):
    """Get each sequence length of FASTA file"""
    main(fasta_file, parse_seqids, outfile)


if __name__ == '__main__':
    run()
