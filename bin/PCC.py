#!/usr/bin/env python
"""
File: PCC.py
Description: Calculation of Pearson correlation coefficient from gene expression.
CreateDate: 2022/9/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union
from tqdm import tqdm
from scipy.stats import pearsonr
import click
from pybioinformatic import TaskManager, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def sub_process(x_name, x_data, y_name, y_data):
    r, p = pearsonr(x_data, y_data)
    return f'{x_name}\t{y_name}\t{r}\t{p}'


def main(exp_matrix_file1: Union[str, TextIOWrapper],
         exp_matrix_file2: Union[str, TextIOWrapper],
         min_exp1: float = 0.5,
         min_exp2: float = 0.5,
         num_processing: int = 10,
         output_file: TextIOWrapper = None):

    exp_dict1 = {
        line.strip().split('\t')[0]: [float(i) for i in line.strip().split('\t')[1:]]
        for line in exp_matrix_file1
        if not line.startswith('Geneid') and line.strip() and
           all([float(i) >= min_exp1 for i in line.strip().split('\t')[1:]])
    }

    exp_dict2 = {
        line.strip().split('\t')[0]: [float(i) for i in line.strip().split('\t')[1:]]
        for line in exp_matrix_file2
        if not line.startswith('Geneid') and line.strip() and
           all([float(i) >= min_exp2 for i in line.strip().split('\t')[1:]])
    }

    with tqdm(total=len(exp_dict1) * len(exp_dict2)) as pbar:
        params = ((k1, v1, k2, v2) for k1, v1 in exp_dict1.items() for k2, v2 in exp_dict2.items())
        tkm = TaskManager(num_processing=num_processing, params=params)
        rets = tkm.parallel_run_func(func=sub_process, call_back_func=lambda _: pbar.update())
        for ret in rets:
            click.echo(ret.get(), output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--exp-file1', 'exp_file1',
              metavar='<file|stdin>', required=True, type=click.File('r'),
              help=r'Input gene expression matrix file. (header must start with "Geneid")')
@click.option('-I', '--exp-file2', 'exp_file2',
              metavar='<file|stdin>', required=True, type=click.File('r'),
              help='Input another gene expression matrix file. (header must start with "Geneid")')
@click.option('-e', '--min-exp1', 'min_exp1',
              metavar='<float>', type=float, default=0.5, show_default=True,
              help="Min expression of gene in exp file1, if some one's expression is less than specified value in all samples, then filter out.")
@click.option('-E', '--min-exp2', 'min_exp2',
              metavar='<float>', type=float, default=0.5, show_default=True,
              help="Min expression of gene in exp file2, if some one's expression is less than specified value in all samples, then filter out.")
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=10, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(exp_file1, exp_file2, min_exp1, min_exp2, num_processing, output_file):
    """Calculation of Pearson correlation coefficient from gene expression."""
    main(
        exp_matrix_file1=exp_file1,
        exp_matrix_file2=exp_file2,
        min_exp1=min_exp1,
        min_exp2=min_exp2,
        num_processing=num_processing,
        output_file=output_file
    )


if __name__ == '__main__':
    run()
