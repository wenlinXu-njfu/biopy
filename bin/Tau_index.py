#!/usr/bin/env python
"""
File: Tau_index.py
Description: Calculate gene expression tissue-specificity based on Tau index.
CreateDate: 2022/10/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from pybioinformatic import read_in_gene_expression_as_dataframe, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(exp_file: str, out_file_prefix: str):
    gene_exp_matrix = read_in_gene_expression_as_dataframe(exp_file)
    if isinstance(gene_exp_matrix, str):
        click.echo(gene_exp_matrix, err=True)
        exit()
    df = gene_exp_matrix.copy()
    df['max'] = df.max(axis=1)
    for i in range(len(df.index.tolist())):
        df.iloc[i] = 1 - df.iloc[i] / df.iloc[i].max()
    df['max'] = df.sum(axis=1)
    df['Tau'] = df['max'] / (len(gene_exp_matrix.columns.tolist()) - 1)
    gene_exp_matrix['Tau'] = df['Tau']
    gene_exp_matrix.to_csv(f'./{out_file_prefix}.csv')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--expression_file', 'expression_file',
              metavar='<exp file>', required=True,
              help='Input gene expression profile file. (support format: txt, xls, xlsx, and csv)')
@click.option('-o', '--output_prefix', 'output_prefix',
              metavar='<str>', default='Tau_index', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(expression_file, output_prefix):
    """Calculate gene expression tissue-specificity based on Tau index."""
    main(expression_file, output_prefix)


if __name__ == '__main__':
    run()
