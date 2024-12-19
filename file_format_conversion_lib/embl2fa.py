#!/usr/bin/env python
"""
File: embl2fa.py
Description: Convert EMBL to FASTA.
CreateDate: 2024/12/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from Bio import SeqIO
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(embl_file: str):
    if not embl_file.endswith('embl'):
        fasta_file = f'{embl_file}.fasta'
    else:
        fasta_file = '.'.join(embl_file.split('.')[:-1]) + '.fasta'
    records = SeqIO.parse(embl_file, "embl")
    count = SeqIO.write(records, fasta_file, "fasta")
    print("Converted %i records" % count)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('embl_file', nargs=1, metavar='<embl file>', required=True)
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(embl_file):
    """Convert EMBL to FASTA."""
    main(embl_file)


if __name__ == '__main__':
    run()
