#!/usr/bin/env python
"""
File: gtf_extract_seq.py
Description: Extract cDNA sequence from GTF file
Date: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gtf import Gtf
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(gtf_file, fasta_file, out_file):
    if out_file:
        with open(out_file, 'a') as o:
            for cDNA_nucl_obj in Gtf(gtf_file).get_cDNA(fasta_file):
                o.write(f">{cDNA_nucl_obj.id}\n{cDNA_nucl_obj.seq}\n")
    else:
        for cDNA_nucl_obj in Gtf(gtf_file).get_cDNA(fasta_file):
            print(cDNA_nucl_obj)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-g', '--gtf_file', 'gtf_file', help='Input GTF file.')
@click.option('-f', '--ref_fasta', 'ref_fasta_file', help='Input reference sequence FASTA file.')
@click.option('-o', '--output_file', 'output_file',
              help='[optional] Output FASTA file, if not specified, print results to terminal as stdout.')
def run(gtf_file, ref_fasta_file, output_file):
    """Extract cDNA sequence from GTF file."""
    main(gtf_file, ref_fasta_file, output_file)


if __name__ == '__main__':
    run()
