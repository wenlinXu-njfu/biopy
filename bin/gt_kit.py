#!/usr/bin/env python
"""
File: gt_kit.py
Description: GT file tools.
CreateDate: 2024/6/1
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from gt_lib import gs, merge, stat
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def gtkit():
    """GT file tools."""
    pass


gtkit.add_command(gs, 'gs')
gtkit.add_command(merge, 'merge')
gtkit.add_command(stat, 'stat')


if __name__ == '__main__':
    gtkit()
