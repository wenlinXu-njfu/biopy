#!/usr/bin/env python
"""
File: exec_cmds.py
Description: Execute commands asynchronously.
CreateDate: 2022/1/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(command_file: TextIOWrapper, processing_num: int):
    tkm = TaskManager([line.strip() for line in command_file], processing_num)
    tkm.parallel_run_cmd()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-f', '--command_file', 'command_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input file including command (one command per line).')
@click.option('-n', '--processing_num', 'processing_num',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(command_file, processing_num):
    """Execute commands asynchronously."""
    main(command_file, processing_num)


if __name__ == '__main__':
    run()
