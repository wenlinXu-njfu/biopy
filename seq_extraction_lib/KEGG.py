#!/usr/bin/env python
"""
File: KEGG.py
Description: Extract sequences from KEGG url.
Date: 2022/6/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os.path import exists
import click
import requests
from tqdm import tqdm
from Biolib import Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(in_file, seq_type: click.Choice(['aaseq', 'ntseq']), out_file):
    set1 = set([line.strip() for line in open(in_file) if line.strip()])
    if exists(out_file):
        set2 = set(line.split(' ')[0].replace('>', '') for line in open(out_file) if line.startswith('>'))
        l = list(set1 - set2)
    else:
        l = list(set1)
    l.sort()
    with open(out_file, 'a') as o:
        with tqdm(total=len(l)) as pbar:
            for gene_id in l:
                url = f'https://rest.kegg.jp/get/{gene_id}/{seq_type}'
                response = requests.get(url)
                o.write(response.text)
                pbar.update()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--seq_id', 'id_file',
              metavar='<id file>', required=True,
              help='Input KEGG list file, one id per line. (eg. pop:112323434\\npop:112323435)')
@click.option('-t', '--seq_type', 'seq_type',
              metavar='<aaseq|ntseq>', type=click.Choice(['aaseq', 'ntseq']), default='aaseq', show_default=True,
              help='Specified sequence type.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file>',
              help='Output FASTA file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(id_file, seq_type, output_file):
    """Extract sequences from KEGG url."""
    main(id_file, seq_type, output_file)


if __name__ == '__main__':
    run()
