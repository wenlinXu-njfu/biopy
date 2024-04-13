#!/usr/bin/env python
"""
File: gtf_extract_seq.py
Description: Extract cDNA sequence from GTF file
CreateDate: 2022/3/25
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import Gtf, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gtf_file: TextIOWrapper,
         fasta_file: TextIOWrapper,
         out_file: TextIOWrapper):
    with Gtf(gtf_file) as gtf:
        for cDNA_nucl_obj in gtf.get_cDNA(fasta_file):
            click.echo(cDNA_nucl_obj, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gtf_file', 'gtf_file',
              metavar='<gtf file|stdin>', type=click.File('r'), required=True,
              help='Input GTF file.')
@click.option('-r', '--ref_fasta', 'ref_fasta_file',
              metavar='<fasta file|stdin>', type=click.File('r'), required=True,
              help='Input reference sequence FASTA file.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file|stdout>', type=click.File('w'),
              help='Output FASTA file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gtf_file, ref_fasta_file, output_file):
    """Extract cDNA sequence from GTF file."""
    main(gtf_file, ref_fasta_file, output_file)


if __name__ == '__main__':
    run()
