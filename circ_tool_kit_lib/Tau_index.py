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
from Biolib.show_info import Displayer


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


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--circ_exp_file', 'circ_exp_file',
              help='Input circRNA expression profile file. (support format: txt, xls, xlsx, csv)')
@click.option('-o', '--output_prefix', 'out_prefix', default='circ_Tau', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(circ_exp_file, out_prefix):
    """Tissue specific expression analysis of circRNAs (Tau index)."""
    main(circ_exp_file, out_prefix)


if __name__ == '__main__':
    run()
