#!/usr/bin/env python
"""
File: KEGG.py
Description: Extract sequences from KEGG url.
CreateDate: 2022/6/19
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import time
import click
import requests
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def download(url):
    response = requests.get(url=url,
                            headers={
                                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"
                            })
    click.echo(response.text)
    time.sleep(1.5)


def main(in_file: TextIOWrapper,
         seq_type: click.Choice(['aaseq', 'ntseq']),
         num_processing: int = 10):
    id_list = list({line.strip() for line in in_file if line.strip()})
    params = ((f'https://rest.kegg.jp/get/{gene_id}/{seq_type}',) for gene_id in id_list)
    tkm = TaskManager(num_processing=num_processing, params=params)
    tkm.parallel_run_func(func=download)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--seq_id', 'id_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input KEGG list file, one id per line. (eg. pop:112323434\npop:112323435)')
@click.option('-t', '--seq_type', 'seq_type',
              metavar='<aaseq|ntseq>', type=click.Choice(['aaseq', 'ntseq']), default='aaseq', show_default=True,
              help='Specified sequence type.')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=10, show_default=True,
              help='Number of processing.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(id_file, seq_type, num_processing):
    """Extract sequences from KEGG url."""
    main(id_file, seq_type, num_processing)


if __name__ == '__main__':
    run()
