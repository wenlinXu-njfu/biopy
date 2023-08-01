#!/usr/bin/env python
"""
File: execute_duplicate_command.py
Description: Execute commands in a file line by line
Date: 2022/1/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import system
from datetime import datetime
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(file):
    for line in open(file):
        if not line.strip():
            continue
        click.echo(f"\033[1m[{datetime.now().replace(microsecond=0)}]\033[0m {line.strip()}", err=True)
        system(command=line.strip())


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--command_file', 'f', help='Input file including command (one command per line).')
def run(f):
    """Execute commands in a file line by line."""
    main(f)


if __name__ == '__main__':
    run()
