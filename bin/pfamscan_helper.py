#!/usr/bin/env python
"""
File: pfamscan_helper.py
Description: PfamScan helper.
Date: 2022/4/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.PfamScan.batch_pfamscan import run as run1
from Biolib.show_info import Displayer


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def pfamscan_helper():
    """PfamScan helper."""
    pass


pfamscan_helper.add_command(run1, 'batch')

if __name__ == '__main__':
    pfamscan_helper()
