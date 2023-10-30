#!/usr/bin/env python
"""
File: gff_extract_seq.py
Description: Extract sequences from GFF file.
CreateDate: 2022/3/24
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Gff, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gff_file: Union[str, TextIOWrapper],
         fa_file: Union[str, TextIOWrapper],
         feature_type: str,
         id_file: TextIOWrapper = None,
         out_file: TextIOWrapper = None):
    with Gff(gff_file) as gff:
        if id_file:
            id_list = set(i.strip() for i in id_file if i.strip())
            for nucl_obj in gff.extract_seq(fa_file, feature_type, id_list):
                click.echo(nucl_obj, out_file)
        else:
            for nucl_obj in gff.extract_seq(fa_file, feature_type):
                click.echo(nucl_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gff_file', 'gff_file',
              metavar='<gff file>', type=click.File('r'), required=True,
              help='Input GFF file.')
@click.option('-r', '--ref_fasta', 'ref_fasta_file',
              metavar='<fasta file>', type=click.File('r'), required=True,
              help='Input reference sequence FASTA file.')
@click.option('-t', '--feature_type', 'feature_type',
              metavar='<str>', required=True,
              help='Specify feature type. (eg: gene, mRNA, etc)')
@click.option('-d', '--id_file', 'id_file',
              metavar='<id file>', type=click.File('r'),
              help='Provides an ID file (one id per line) that extracts sequences from the GFF file that '
                   '\033[1mmatch\033[0m the IDs in the ID file.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file>', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gff_file, ref_fasta_file, feature_type, id_file, output_file):
    """Extract sequences from GFF file."""
    main(gff_file, ref_fasta_file, feature_type, id_file, output_file)


if __name__ == '__main__':
    run()
