#!/usr/bin/env python
"""
File: PCC.py
Description: Calculation of Pearson correlation coefficient from gene expression.
CreateDate: 2022/9/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from scipy.stats import pearsonr
import click
from pybioinformatic import read_in_gene_expression_as_dataframe, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(exp_matrix_file: str, output_prefix: str):
    with open(f'./{output_prefix}.fmt1.xls', 'w') as fmt1:
        df = read_in_gene_expression_as_dataframe(exp_matrix_file)
        if isinstance(df, str):
            click.echo(df, err=True)
            exit()
        else:
            df2 = df.T
            df2.columns.name = None
            df2 = df2.corr()
            i = 0
            j = 1
            click.echo('Gene1\tGene2\tPCC\tP_value', fmt1)
            while i < len(df.index.tolist()):
                while j < len(df.index.tolist()):
                    x = df.iloc[i]
                    y = df.iloc[j]
                    r, p = pearsonr(x, y)
                    click.echo(f'{x.name}\t{y.name}\t{r}\t{p}', fmt1)
                    df2.iloc[i, j] = ''
                    j += 1
                i += 1
                j = i + 1
            df2.to_csv(f'./{output_prefix}.fmt2.xls', sep='\t')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--exp_file', 'input_file',
              metavar='<exp file>', required=True,
              help='Input gene expression matrix file (including header). Supported formats: txt, xls, xlsx, csv')
@click.option('-o', '--output_prefix', 'output_prefix',
              metavar='<str>', default='PCC', show_default=True,
              help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, output_prefix):
    """Calculation of Pearson correlation coefficient from gene expression."""
    main(input_file, output_prefix)


if __name__ == '__main__':
    run()
