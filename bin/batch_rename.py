#!/usr/bin/env python
"""
File: batch_rename.py
Description: Batch rename files.
Date: 2021-10-06
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import listdir, rename
from re import sub
import click
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(in_dir, old, new):
    files = listdir(in_dir)
    for file in files:
        if new:
            s = sub(old, new, file)
        else:
            s = sub(old, '', file)
        rename(f'{in_dir}/{file}', f'{in_dir}/{s}')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-d', '--input_dir_path', 'input_dir',
              metavar='<dir>', required=True,
              help='The directory where the file to be renamed resides.')
@click.option('-old', '--old_name', 'old',
              metavar='<str>',
              help='The string to be replaced, it supports for regular expressions')
@click.option('-new', '--new_name', 'new',
              metavar='<str>',
              help='Replacement string, it supports for regular expressions')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_dir, old, new):
    """Batch rename files."""
    main(input_dir, old, new)


if __name__ == '__main__':
    run()
