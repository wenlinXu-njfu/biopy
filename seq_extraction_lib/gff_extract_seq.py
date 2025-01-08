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
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(gff_file: Union[str, TextIOWrapper],
         fa_file: Union[str, TextIOWrapper],
         feature_type: str,
         id_field: str,
         id_file: TextIOWrapper = None,
         remain_attr: list = None,
         num_processing: int = 4,
         out_file: TextIOWrapper = None):
    id_set = set(i.strip() for i in id_file if i.strip()) if id_file else None
    with Gff(gff_file) as gff:
        generator = gff.extract_seq(
            fasta_file=fa_file,
            feature_type=feature_type,
            id_field=id_field,
            feature_id_set=id_set,
            remain_attr=remain_attr,
            num_processing=num_processing
        )
        for nucl_obj in generator:
            click.echo(nucl_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gff-file', 'gff_file',
              metavar='<gff file|stdin>', type=click.File('r'), required=True,
              help='Input GFF file.')
@click.option('-r', '--ref-fasta', 'ref_fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input reference sequence FASTA file.')
@click.option('-t', '--feature-type', 'feature_type',
              metavar='<str>', required=True,
              help='Specify feature type of sequence. (Column 3 of the gff file, eg: gene, mRNA, etc)')
@click.option('-I', '--id-field', 'id_field',
              metavar='<str>', default='ID', show_default=True,
              help='Specify the gff file column 9 field as the sequence id.')
@click.option('-i', '--id-file', 'id_file',
              metavar='<id file|stdin>', type=click.File('r'),
              help='Provides an ID file (one id per line) that extracts sequences from the GFF file that '
                   '\033[1mmatch\033[0m the IDs in the ID file.')
@click.option('-a', '--other-infor', 'other_infor', metavar='<str>',
              help='Information to retain in column 9 of the gff file, multiple messages are separated by commas. (eg. Parent,pacid)')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=4, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<fasta file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gff_file, ref_fasta_file, feature_type, id_file, id_field, other_infor, num_processing, output_file):
    """Extract sequences from GFF file."""
    other_infor = other_infor.split(',') if other_infor else None
    main(
        gff_file=gff_file,
        fa_file=ref_fasta_file,
        feature_type=feature_type,
        id_field=id_field,
        id_file=id_file,
        remain_attr=other_infor,
        num_processing=num_processing,
        out_file=output_file
    )


if __name__ == '__main__':
    run()
