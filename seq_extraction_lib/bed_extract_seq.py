#!/usr/bin/env python
"""
File: bed_extract_seq.py
Description: Extract sequences from BED file.
CreateDate: 2022/3/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Bed, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(bed_file: Union[str, TextIOWrapper],
         fa_file: Union[str, TextIOWrapper],
         use_id: bool = True,
         up: int = 0,
         down: int = 0,
         both: int = 0,
         extension: bool = True,
         out_file: TextIOWrapper = None):
    with Bed(bed_file) as bed:
        for nucl_obj in bed.extract_seq(fa_file, use_id, up, down, both, extension):
            click.echo(nucl_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--bed_file', 'bed_file',
              metavar='<bed file|stdin>', type=click.File('r'), required=True,
              help=r'Input BED file (Chr_num\tStart\tEnd\tFrame\tName\tStrand\tEtc).')
@click.option('-r', '--ref_fasta', 'ref_fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input reference sequence FASTA file.')
@click.option('-up', '--upstream', 'upstream',
              metavar='<int>', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent upstream the specified length.')
@click.option('-down', '--downstream', 'downstream',
              metavar='<int>', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent downstream the specified length.')
@click.option('-both', '--both_end', 'both_end',
              metavar='<int>', type=int, default=0, show_default=True,
              help='Make sequence in bed file to extent both end the specified length. '
                   'If "-u --upstream", "-d --downstream", and "-b --both_end" are specified, '
                   'by default, only "--both_end" is valid.')
@click.option('-U', '--use_id', 'use_id',
              is_flag=True, flag_value=True,
              help='Use the fourth column content in the BED file as sequence ID. '
                   'Otherwise, "Chr_num:start-end(strand)" by default.')
@click.option('--extension/--non_extension',
              default=True, show_default=True,
              help='If "-u --upstream", "-d --downstream" or "-b --both_end" is specified, '
                   'specify whether extension (including sequence self in BED file).')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file|stdout>', type=click.File('w'),
              help='Output FASTA file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(bed_file, ref_fasta_file, use_id, upstream, downstream, both_end, extension, output_file):
    """Extract sequences from BED file."""
    main(bed_file, ref_fasta_file, use_id, upstream, downstream, both_end, extension, output_file)


if __name__ == '__main__':
    run()
