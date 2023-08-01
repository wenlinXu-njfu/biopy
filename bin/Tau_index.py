#!/usr/bin/env python
"""
File: Tau_index.py
Description: Calculate gene expression tissue-specificity based on Tau index.
Date: 2022/10/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.statistics import read_in_gene_expression_as_dataframe
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(exp_file, out_file_prefix):
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


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--expression_file', 'expression_file',
              help='Input gene expression profile file. (support format: txt, xls, xlsx, csv)')
@click.option('-o', '--output_file', 'outfile', help='Output file prefix.')
def run(expression_file, outfile):
    """Calculate gene expression tissue-specificity based on Tau index."""
    main(expression_file, outfile)


if __name__ == '__main__':
    run()
