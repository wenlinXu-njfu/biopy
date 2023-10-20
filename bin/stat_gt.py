#!/usr/bin/env python
"""
File: stat_gt.py
Description: Calculate the miss rate, heterozygosity rate and MAF of SNP sites from GT files.
Date: 2023/10/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from re import sub
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1])


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-gt', '--gt-file', 'gt_file',
              type=click.File('r'), metavar='<GT file>', required=True,
              help='Input GT file.')
@click.option('-n', '--num-processing', 'num_processing',
              type=int, metavar='<int>', default=1, show_default=True, help='Number of Processing.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>',
              help='Output file path and name. Print to terminal by default.')
@Timer('Calculating MissRate, HetRate, and MAF.')
def run(gt_file, num_processing, output_file):
    gt = GenoType(gt_file)
    stat_df = gt.parallel_stat_MHM(num_processing)
    if output_file:
        stat_df.to_csv(output_file, sep='\t')
    else:
        text = stat_df.to_string(index=False)
        text = sub(r' +', '\t', text)
        text = sub(r'\n\t', '\n', text).strip()
        click.echo(text)


if __name__ == '__main__':
    run()
