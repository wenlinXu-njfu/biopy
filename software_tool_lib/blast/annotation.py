#!/usr/bin/env python
"""
File: annotation.py
Description: Annotate the blast results for species and functions.
CreateDate: 2024/8/9
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from collections import defaultdict
from gzip import GzipFile
from json import loads
import click
from pybioinformatic import Blast, Displayer
parent_dir = '/'.join(__file__.split('/')[:-1])
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(blast_result: TextIOWrapper,
         nr_database: str,
         output_file: str,
         align_rate: float = 90.0,
         align_len: int = 100,
         e_value: float = 10e-5):
    if nr_database.endswith('gz'):
        with GzipFile(nr_database) as f:
            d = loads(f.read())
    else:
        with open(nr_database) as f:
            d = loads(f.read())

    nr_database = defaultdict(lambda _=None: {'species': 'NA', 'annotation': 'NA'})
    for k, v in d.items():
        nr_database[k] = v

    with Blast(blast_result) as blast:
        df = blast.to_dataframe()
        filter_df = df[(df.AlignRate >= align_rate) & (df.AlignLength >= align_len) & (df.Evalue <= e_value)]
        filter_df['Species'] = filter_df.apply(func=lambda row: nr_database[row['Sbject']]['species'], axis=1)
        filter_df['Annotation'] = filter_df.apply(func=lambda row: nr_database[row['Sbject']]['annotation'], axis=1)
        filter_df.to_csv(output_file, sep='\t', index=False, na_rep='NA')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--blast-result', 'blast_result',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Enter blast result file with format 6.')
@click.option('-d', '--nr-database', 'nr_database',
              metavar='<file>', default=f'{parent_dir}/accession_species_anno.json.gz', show_default=True,
              help='Enter nr database json file.')
@click.option('-r', '--align-rate', 'align_rate',
              metavar='<float>', type=float, default=90.0, show_default=True,
              help='Align ratio threshold.')
@click.option('-l', '--align-length', 'align_len',
              metavar='<int>', type=int, default=100, show_default=True,
              help='Align length threshold.')
@click.option('-e', '--e-value', 'e_value',
              metavar='<float>', type=float, default=10e-5, show_default=True,
              help='E value threshold.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', required=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(blast_result, nr_database, align_rate, align_len, e_value, output_file):
    """Annotate the blast results for species and functions."""
    main(
        blast_result=blast_result,
        nr_database=nr_database,
        output_file=output_file,
        align_rate=align_rate,
        align_len=align_len,
        e_value=e_value
    )


if __name__ == '__main__':
    run()
