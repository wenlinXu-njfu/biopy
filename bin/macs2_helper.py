#!/usr/bin/env python
"""
File: macs2_helper.py
Description: 
CreateDate: 2025/6/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from software_tool_lib.macs2.cent_identifier import run as cent_identifier
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def macs2_helper():
    """Macs2 software helper."""
    pass


macs2_helper.add_command(cent_identifier, 'cent_identifier')

if __name__ == '__main__':
    macs2_helper()
