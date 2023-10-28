#!/usr/bin/env python
"""
File: blast_helper.py
Description: Blast software helper
CreateDate: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.blast.blast2bed import run as run1
from software_tool_lib.blast.reciprocal_blast import run as run2
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def blast_helper():
    """Blast software helper."""
    pass


blast_helper.add_command(run1, 'blast2bed')
blast_helper.add_command(run2, 'reciprocal_blast')

if __name__ == '__main__':
    blast_helper()
