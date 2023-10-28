#!/usr/bin/env python
"""
File: stat_gt.py
Description: Calculate the miss rate, heterozygosity rate and MAF of SNP sites from GT files.
Date: 2023/10/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from matplotlib.pyplot import rcParams, style, figure, subplots_adjust, savefig
from seaborn import histplot
import click
from pybioinformatic import GenoType, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gt_file: Union[str, TextIOWrapper],
         num_processing: int,
         output_path: str):
    # stat MHM
    gt = GenoType(gt_file)
    stat_df = gt.parallel_stat_MHM(num_processing)
    # plot figure
    style.use('ggplot')
    rcParams['font.family'] = 'Arial'
    rcParams['font.size'] = 8
    fig = figure(figsize=(12, 5), dpi=300)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    histplot(stat_df, x='MissRate(%)', ax=ax1, bins=20)
    histplot(stat_df, x='HetRate(%)', ax=ax2, bins=20)
    histplot(stat_df, x='MAF', ax=ax3, bins=20)
    ax1.set_xticks(range(0, 110, 10), range(0, 110, 10), rotation=45)
    ax2.set_xticks(range(0, 110, 10), range(0, 110, 10), rotation=45)
    subplots_adjust(wspace=0.3)
    # output file
    stat_df['MissRate(%)'] = stat_df.loc[:, ['MissRate(%)']].applymap(lambda value: '%.2f' % value)
    stat_df['HetRate(%)'] = stat_df.loc[:, ['HetRate(%)']].applymap(lambda value: '%.2f' % value)
    stat_df['MAF'] = stat_df.loc[:, ['MAF']].applymap(lambda value: '%.3f' % value)
    stat_df.to_csv(f'{output_path}/site_stat.xls', sep='\t', index=False)
    savefig(f'{output_path}/distribution.png', bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-gt', '--gt-file', 'gt_file',
              type=click.File('r'), metavar='<GT file>', required=True,
              help='Input GT file.')
@click.option('-n', '--num-processing', 'num_processing',
              type=int, metavar='<int>', default=1, show_default=True, help='Number of Processing.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default='./', show_default=True,
              help='Output file path.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
@Timer('Calculating MissRate, HetRate, and MAF.')
def run(gt_file, num_processing, output_path):
    main(gt_file, num_processing, output_path)


if __name__ == '__main__':
    run()
