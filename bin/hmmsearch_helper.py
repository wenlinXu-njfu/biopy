#!/usr/bin/env python
"""
File: hmmsearch_helper.py
Description: Hmmsearch helper
Date: 2022/4/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.hmmsearch.extract_seq_id import run as run1
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def hmmsearch_helper():
    """Hmmsearch helper."""
    pass


hmmsearch_helper.add_command(run1, 'extract_id')

if __name__ == '__main__':
    hmmsearch_helper()
