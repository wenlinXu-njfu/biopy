#!/usr/bin/env python
"""
File: bed_extract_seq.py
Description: Extract sequences from BED file.
Date: 2022/3/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from _io import TextIOWrapper
import click
from Biolib.bed import Bed
from Biolib.show_info import Displayer


def main(bed_file: TextIOWrapper,
         fa_file: TextIOWrapper,
         use_id: bool,
         up: int,
         down: int,
         both: int,
         extension: bool,
         out_file: TextIOWrapper = None):
    content = []
    for nucl_obj in Bed(bed_file).bed_extract_seq(fa_file, use_id, up, down, both, extension):
        content.append(f">{nucl_obj.id}\n{nucl_obj.seq}\n") if out_file else print(f">{nucl_obj.id}\n{nucl_obj.seq}")
    if out_file:
        with out_file as o:
            o.write(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--bed_file', 'bed_file', type=click.File('r'),
              help='Input BED file.\n(chr_num\\tstart\\tend\\tframe\\tname\\tstrand\\tetc)')
@click.option('-r', '--ref_fasta', 'ref_fasta_file', type=click.File('r'), help='Input reference sequence FASTA file.')
@click.option('-u', '--upstream', 'upstream', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent upstream the specified length.')
@click.option('-d', '--downstream', 'downstream', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent downstream the specified length.')
@click.option('-b', '--both_end', 'both_end', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent both end the specified length. '
                   'If "-u --upstream", "-d --downstream", and "-b --both_end" are specified, '
                   'by default, only "--both_end" is valid.')
@click.option('-U', '--use_id', 'use_id', is_flag=True, flag_value=True,
              help='Use the fourth column content in the BED file as sequence ID. '
                   'Otherwise, "Chr_num:start-end(strand)" by default.')
@click.option('--extension/--non_extension', type=bool, default=True, show_default=True,
              help='If "-u --upstream", "-d --downstream" or "-b --both_end" is specified, '
                   'specify whether extension (including sequence self in BED file).')
@click.option('-o', '--output_file', 'output_file', type=click.File('w'),
              help='Output FASTA file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(bed_file, ref_fasta_file, use_id, upstream, downstream, both_end, extension, output_file):
    """Extract sequences from BED file."""
    main(bed_file, ref_fasta_file, use_id, upstream, downstream, both_end, extension, output_file)


if __name__ == '__main__':
    run()
