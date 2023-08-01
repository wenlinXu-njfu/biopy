#!/usr/bin/env python
"""
File: featureCounts_helper.py
Description: Standardize gene expression with FPKM or TPM based on featureCounts results
Date: 2022/3/29
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.featureCounts.get_FPKM import run as run1
from software_tool_lib.featureCounts.get_TPM import run as run2
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def featureCounts_helper():
    """Standardize gene expression with FPKM or TPM based on featureCounts results."""
    pass


featureCounts_helper.add_command(run1, 'FPKM')
featureCounts_helper.add_command(run2, 'TPM')

if __name__ == '__main__':
    featureCounts_helper()
