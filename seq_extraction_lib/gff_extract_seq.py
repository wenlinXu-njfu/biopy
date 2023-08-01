#!/usr/bin/env python
"""
File: gff_extract_seq.py
Description: Extract sequences from GFF file
Date: 2022/3/24
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gff import Gff
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(gff_file, fa_file, feature_type, id_file, out_file):
    content = ''
    if id_file:
        id_list = [i.strip() for i in open(id_file).readlines() if i.strip()]
        for nucl_obj in Gff(gff_file).gff_extract_seq(fa_file, feature_type, id_list):
            if out_file:
                content += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
            else:
                print(nucl_obj)
    else:
        for nucl_obj in Gff(gff_file).gff_extract_seq(fa_file, feature_type):
            if out_file:
                content += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
            else:
                print(nucl_obj)
    if content:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-g', '--gff_file', 'gff_file', help='Input GFF file.')
@click.option('-f', '--ref_fasta', 'ref_fasta_file', help='Input reference sequence FASTA file.')
@click.option('-t', '--feature_type', 'feature_type', help='Specify feature type. (eg: gene, mRNA, etc)')
@click.option('-d', '--id_file', 'id_file',
              help='[optional] Provides an ID file (one id per line) that extracts sequences from the GFF file that '
                   '\033[1mmatch\033[0m the IDs in the ID file.')
@click.option('-o', '--output_file', 'output_file',
              help='[optional] Output file, if not specified, print results to terminal as stdout.')
def run(gff_file, ref_fasta_file, feature_type, id_file, output_file):
    """Extract sequences from GFF file."""
    main(gff_file, ref_fasta_file, feature_type, id_file, output_file)


if __name__ == '__main__':
    run()
