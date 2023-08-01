#!/usr/bin/env python
"""
File: pfamscan_helper.py
Description: PfamScan helper
Date: 2022/4/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.PfamScan.batch_pfamscan import run as run1
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def pfamscan_helper():
    """PfamScan helper."""
    pass


pfamscan_helper.add_command(run1, 'batch')

if __name__ == '__main__':
    pfamscan_helper()
