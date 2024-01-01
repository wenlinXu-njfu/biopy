#!/usr/bin/env python
"""
File: fq2fa.py
Description: Convert FASTQ to FASTA
CreateDate: 2022/4/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from gzip import GzipFile
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(fq_file: str, out_file: TextIOWrapper):
    try:
        with open(fq_file) as f:
            while 1:
                read_id = f.readline().replace('@', '>')
                if not read_id:
                    break
                seq = f.readline()
                f.readline()
                f.readline()
                click.echo(read_id.strip(), out_file)
                click.echo(seq.strip(), out_file)
    except UnicodeDecodeError:
        with GzipFile(fq_file) as f:
            while 1:
                read_id = str(f.readline(), 'utf8').replace('@', '>')
                if not read_id:
                    break
                seq = str(f.readline(), 'utf8')
                f.readline()
                f.readline()
                click.echo(read_id.strip(), out_file)
                click.echo(seq.strip(), out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fastq_file', 'fastq_file',
              metavar='<fastq file>', required=True,
              help='Input FASTQ file(XXX.fq) or FASTQ compressed files(XXX.fq.gz).')
@click.option('-o', '--output_fasta', 'fasta_file',
              metavar='<fasta file|stdout>', type=click.File('r'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fastq_file, fasta_file):
    """Convert FASTQ to FASTA."""
    main(fastq_file, fasta_file)


if __name__ == '__main__':
    run()
