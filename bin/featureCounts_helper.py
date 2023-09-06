#!/usr/bin/env python
"""
File: featureCounts_helper.py
Description: Standardize gene expression with FPKM or TPM based on featureCounts results.
Date: 2022/3/29
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.featureCounts.get_FPKM import run as run1
from software_tool_lib.featureCounts.get_TPM import run as run2
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def featureCounts_helper():
    """Standardize gene expression with FPKM or TPM based on featureCounts results."""
    pass


featureCounts_helper.add_command(run1, 'FPKM')
featureCounts_helper.add_command(run2, 'TPM')

if __name__ == '__main__':
    featureCounts_helper()
