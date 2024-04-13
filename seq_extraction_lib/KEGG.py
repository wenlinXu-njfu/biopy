#!/usr/bin/env python
"""
File: KEGG.py
Description: Extract sequences from KEGG url.
CreateDate: 2022/6/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
import requests
from natsort import natsort_key
from tqdm import tqdm
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(in_file: TextIOWrapper,
         seq_type: click.Choice(['aaseq', 'ntseq']),
         out_file: TextIOWrapper):
    id_list = list({line.strip() for line in in_file if line.strip()})
    id_list.sort(key=natsort_key)
    with tqdm(total=len(id_list)) as pbar:
        for gene_id in id_list:
            url = f'https://rest.kegg.jp/get/{gene_id}/{seq_type}'
            response = requests.get(url)
            click.echo(response.text, out_file)
            pbar.update()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--seq_id', 'id_file',
              metavar='<id file|stdin>', type=click.File('r'), required=True,
              help=r'Input KEGG list file, one id per line. (eg. pop:112323434\npop:112323435)')
@click.option('-t', '--seq_type', 'seq_type',
              metavar='<aaseq|ntseq>', type=click.Choice(['aaseq', 'ntseq']), default='aaseq', show_default=True,
              help='Specified sequence type.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<fasta file>', type=click.File('w'), required=True,
              help='Output FASTA file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(id_file, seq_type, output_file):
    """Extract sequences from KEGG url."""
    main(id_file, seq_type, output_file)


if __name__ == '__main__':
    run()
