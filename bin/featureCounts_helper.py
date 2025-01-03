#!/usr/bin/env python
"""
File: featureCounts_helper.py
Description: FeatureCounts software helper.
CreateDate: 2022/3/29
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.featureCounts.normalization import run as run1
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def featureCounts_helper():
    """FeatureCounts software helper."""
    pass


featureCounts_helper.add_command(run1, 'normalization')

if __name__ == '__main__':
    featureCounts_helper()
