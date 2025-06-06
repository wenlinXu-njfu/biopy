#!/usr/bin/env python
"""
File: joint2.py
Description: Perform inner or outer joins on a group of table files based on the common field names of the table files.
CreateDate: 2024/10/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union, Tuple, Literal
from pandas import read_table,concat
from natsort import natsort_key
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(master_table: Union[str, TextIOWrapper],
         other_tables: Tuple[Union[str, TextIOWrapper]],
         output_file: TextIOWrapper,
         joint_method: Literal['left', 'outer', 'inner'] = 'left'):
    # read in raw table
    master_table = read_table(master_table, dtype=str)
    other_tables = [
        read_table(other_table, dtype=str)
        for other_table in other_tables
    ]

    # reset table index
    all_columns = [set(other_table.columns.tolist()) for other_table in other_tables]
    intersection = set(master_table.columns.tolist()).intersection(*all_columns)
    index = [i for i in master_table.columns if i in intersection]
    master_table.set_index(keys=index, drop=True, inplace=True)
    other_tables = [
        other_table.set_index(keys=index, drop=True)
        for other_table in other_tables
    ]

    # merge tables
    other_tables.insert(0, master_table)
    if joint_method == 'left':
        merged_table = concat(other_tables, axis=1, join='outer')
        merged_table = merged_table.loc[master_table.index.tolist()]
    else:
        merged_table = concat(other_tables, axis=1, join=joint_method)
    merged_table.reset_index(inplace=True)
    merged_table.sort_values(by=index, key=natsort_key, inplace=True)
    merged_table.to_csv(output_file, sep='\t', na_rep='NA', index=False)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('table_files', nargs=-1, metavar='<table files|stdin>', type=click.File('r'), required=True)
@click.option('-i', '--master-table', 'master_table_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input master table file.')
@click.option('-j', '--joint-method', 'joint_method',
              type=click.Choice(['left', 'outer', 'inner']),
              default='left', show_default=True,
              help='How to handle the operation of these tables.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', type=click.File('w'), default='joint.xls', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(master_table_file, table_files, joint_method, output_file):
    """Perform left inner or outer joins on a group of table files based on the common field names of the table files."""
    main(
        master_table=master_table_file,
        other_tables=table_files,
        output_file=output_file,
        joint_method=joint_method
    )


if __name__ == '__main__':
    run()
