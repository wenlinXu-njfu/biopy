#!/usr/bin/env python
"""
File: execute_duplicate_command.py
Description: Execute commands in a file line by line.
Date: 2022/1/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from multiprocessing import Pool
from io import TextIOWrapper
import click
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def sub_processing(command: str):
    displayer.echo_and_execute_command(command)


def main(command_file: TextIOWrapper, processing_num: int):
    with Pool(processing_num) as pool:
        pool.map(sub_processing, (line.strip() for line in command_file))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-f', '--command_file', 'command_file', type=click.File('r'), required=True,
              help='Input file including command (one command per line).')
@click.option('-n', '--processing_num', 'processing_num', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(command_file, processing_num):
    """Execute commands in a file line by line."""
    main(command_file, processing_num)


if __name__ == '__main__':
    run()
