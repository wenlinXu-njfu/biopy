#!/usr/bin/env python
"""
File: SSR_finder.py
Description: 
Date: 2023/9/9
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1])
def main(input_file, output_file):
    pass
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', 'input_file', type=click.File('r'), help='Input file.')
@click.option('-o', '--output_file', 'output_file', type=click.File('w'), help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, output_file):
    """Description."""
    main(input_file, output_file)
if __name__ == '__main__':
    run()
