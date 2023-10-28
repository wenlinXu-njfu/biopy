#!/usr/bin/env python
"""
File: get_circ_exp.py
Description: Standardize circRNAs expression with CPM.
CreateDate: 2023/3/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from pybioinformatic import read_in_gene_expression_as_dataframe, Displayer
displayer = Displayer(__file__.split('/')[-1], version='1.0.0')


def main(BSJ_matrix_file, out_file):
    df = read_in_gene_expression_as_dataframe(BSJ_matrix_file)
    if isinstance(df, str):
        click.echo(df, err=True)
        exit()
    read_sum = df.sum(axis=0)
    df = df.div(read_sum, axis=1) * 10 ** 6
    df.to_csv(f'./{out_file}.csv')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--bsj_matrix', 'bsj_matrix',
              metavar='<exp file>', required=True,
              help='Input BSJ matrix file. (support format: txt, xls, xlsx, csv)')
@click.option('-o', '--output_prefix', 'output_prefix',
              metavar='<str>', default='circ_CPM', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(bsj_matrix, output_prefix):
    """Standardize circRNAs expression with CPM."""
    main(bsj_matrix, output_prefix)


if __name__ == '__main__':
    run()
