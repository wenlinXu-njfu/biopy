#!/usr/bin/env python
"""
File: execute_duplicate_command.py
Description: Execute commands in a file line by line.
Date: 2022/1/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import system
from datetime import datetime
import click
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(file):
    for line in file:
        if not line.strip():
            continue
        click.echo(f"\033[1m[{datetime.now().replace(microsecond=0)}]\033[0m {line.strip()}", err=True)
        system(command=line.strip())


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-f', '--command_file', 'f', type=click.File('r'), required=True,
              help='Input file including command (one command per line).')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(f):
    """Execute commands in a file line by line."""
    main(f)


if __name__ == '__main__':
    run()
