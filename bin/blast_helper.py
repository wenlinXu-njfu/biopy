#!/usr/bin/env python
"""
File: blast_helper.py
Description: Blast software helper
Date: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.blast.blast2bed import run as run1
from software_tool_lib.blast.reciprocal_blast import run as run2
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def blast_helper():
    """Blast software helper."""
    pass


blast_helper.add_command(run1, 'blast2bed')
blast_helper.add_command(run2, 'reciprocal_blast')

if __name__ == '__main__':
    blast_helper()
