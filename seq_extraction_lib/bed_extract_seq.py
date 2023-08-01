#!/usr/bin/env python
"""
File: bed_extract_seq.py
Description: Extract sequences from BED file
Date: 2022/3/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.bed import Bed
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(bed_file, fa_file, non_id, up, down, both, extension, out_file: str = None):
    bed_file_obj = Bed(bed_file)
    nucl_obj_generator = bed_file_obj.bed_extract_seq(fa_file, non_id, up, down, both, extension)
    content = ''
    for nucl_obj in nucl_obj_generator:
        if out_file:
            content += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
        else:
            print(f">{nucl_obj.id}\n{nucl_obj.seq}")
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--bed_file', 'bed_file', help='Input BED file. (chr_num\\tstart\\tend\\tframe\\tname\\tstrand\\tetc)')
@click.option('-f', '--ref_fasta', 'ref_fasta_file', help='Input reference sequence FASTA file.')
@click.option('-non_id/-id', type=bool, default=True,
              help='''If "-id" specified, the sequence ID will be fourth column content in the BED file. 
              Otherwise, "Chr_num:start-end(strand)" by default.''')
@click.option('-up', '--upstream', 'upstream', type=int, default=0,
              help='Make sequence in bed file to extent upstream the specified length. {default: 0}')
@click.option('-down', '--downstream', 'downstream', type=int, default=0,
              help='Make sequence in bed file to extent downstream the specified length. {default: 0}')
@click.option('-b', '--both', 'booth', type=int, default=0,
              help='''Make sequence in bed file to extent both end the specified length {default: 0}. 
              If "up" "down" and "both" are specified, by default, only both is specified.''')
@click.option('-extension/-non_extension', type=bool, default=True,
              help='''If "up" "down" or "both" is specified, specify whether extension 
              (including sequence self in BED file). {default: extension}''')
@click.option('-o', '--output_file', 'output_file',
              help='[optional] Output FASTA file, if not specified, print results to terminal as stdout.')
def run(bed_file, ref_fasta_file, non_id, upstream, downstream, both, extension, output_file):
    """Extract sequences from BED file."""
    main(bed_file, ref_fasta_file, non_id, upstream, downstream, both, extension, output_file)


if __name__ == '__main__':
    run()
