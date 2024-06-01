#!/usr/bin/env python
"""
File: merge.py
Description: Merges gt files based on sites intersection or union.
CreateDate: 2024/6/1
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union
from typing_extensions import Literal
import click
from pybioinformatic import GenoType, dataframe_to_str, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(gt_file: Union[str, TextIOWrapper],
         another_gt: Union[str, TextIOWrapper],
         mode: Literal['inner', 'outer'] = 'inner',
         output_file: Union[TextIOWrapper, None] = None):
    with GenoType(gt_file) as gt1, GenoType(another_gt) as gt2:
        merge_gt = gt1.merge(other=gt2, how=mode)
        merge_gt.reset_index(inplace=True)
        merge_gt = dataframe_to_str(df=merge_gt, index=False)
        click.echo(merge_gt, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1', '--gt-file', 'gt_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Enter GT file.')
@click.option('-2', '--another-gt', 'another_gt',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Enter another GT file.')
@click.option('-m', '--mode', 'mode',
              type=click.Choice(['inner', 'outer']), default='inner', show_default=True,
              help='Merge mode.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gt_file, another_gt, mode, output_file):
    """Merges gt files based on sites intersection or union."""
    main(
        gt_file=gt_file,
        another_gt=another_gt,
        mode=mode,
        output_file=output_file
    )


if __name__ == '__main__':
    run()
