#!/usr/bin/env python
"""
File: Tau_index.py
Description: Tissue specific expression analysis of circRNAs (Tau index).
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
@click.option('-i', '--circ_exp_file', 'circ_exp_file',
              help='Input circRNA expression profile file. (support format: txt, xls, xlsx, csv)')
@click.option('-o', '--output_prefix', 'out_prefix', default='circ_Tau',
              help='[optional] Output file prefix. {default=circ_Tau}')
def run(circ_exp_file, out_prefix):
    """Tissue specific expression analysis of circRNAs (Tau index)."""
    main(circ_exp_file, out_prefix)


if __name__ == '__main__':
    run()
