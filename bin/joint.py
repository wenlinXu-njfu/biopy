#!/usr/bin/env python
"""
File: joint.py
Description: Joint two table files by different ways ('left', 'right', 'outer', 'inner', or 'cross').
CreateDate: 2024/10/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union, Literal
from pandas import read_table
from natsort import natsort_key
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(left_table: Union[str, TextIOWrapper],
         right_table: Union[str, TextIOWrapper],
         left_column: int,
         right_column: int,
         joint_method: Literal['left', 'right', 'outer', 'inner', 'cross'],
         output_file: TextIOWrapper):
    left_table = read_table(left_table, index_col=left_column - 1, dtype=str)
    right_table = read_table(right_table, index_col=right_column - 1, dtype=str)
    merge = left_table.join(other=right_table, how=joint_method, lsuffix='_left_table', rsuffix='_right_table')
    merge.sort_index(key=natsort_key, inplace=True)
    merge.drop_duplicates(inplace=True)
    merge.to_csv(output_file, sep='\t', na_rep='NA')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--left-table', 'left_table',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input left table file.')
@click.option('-I', '--right-table', 'right_table',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input right table file.')
@click.option('-lc', '--left-column', 'left_column',
              metavar='<int>', default=1, show_default=True,
              help='Left table column to join on the right table.')
@click.option('-rc', '--right-column', 'right_column',
              metavar='<int>', default=1, show_default=True,
              help='Right table column to join on the left table.')
@click.option('-j', '--joint-method', 'joint_method',
              type=click.Choice(['left', 'right', 'outer', 'inner', 'cross']),
              default='left', show_default=True,
              help='How to handle the operation of the two tables.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', type=click.File('w'), default='joint.xls', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(left_table, right_table, left_column, right_column, joint_method, output_file):
    """Joint two table files by different ways ('left', 'right', 'outer', 'inner', or 'cross')."""
    main(
        left_table=left_table,
        right_table=right_table,
        left_column=left_column,
        right_column=right_column,
        joint_method=joint_method,
        output_file=output_file
    )


if __name__ == '__main__':
    run()
