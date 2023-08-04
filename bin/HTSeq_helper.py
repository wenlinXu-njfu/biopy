#!/usr/bin/env python
"""
File: HTSeq_helper.py
Description: Standardize gene expression with FPKM or TPM based on HTSeq results.
Date: 2023/6/15
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.HTseq.get_FPKM import run as run1
from software_tool_lib.HTseq.get_TPM import run as run2
from Biolib.show_info import Displayer


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def HTSeq_helper():
    """Standardize gene expression with FPKM or TPM based on HTSeq results."""
    pass


HTSeq_helper.add_command(run1, 'FPKM')
HTSeq_helper.add_command(run2, 'TPM')

if __name__ == '__main__':
    HTSeq_helper()
