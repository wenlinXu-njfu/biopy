#!/usr/bin/env python
"""
File: biogrep.py
Description: This program works just like fishing, it helps you to get things that you wanted from a target file.
CreateDate: 2023/12/31
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from re import sub
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(bait_file: TextIOWrapper,
         fish_file: TextIOWrapper,
         bait_column: int,
         fish_column: int,
         match: bool,
         invert_match: bool,
         output_file: TextIOWrapper):
    baits = set(sub(r' {2,}', '\t', line).split('\t')[bait_column - 1].strip()
                for line in bait_file if line.strip())
    for line in fish_file:
        if line.strip():
            fish = sub(r' {2,}', '\t', line).split('\t')[fish_column - 1].strip()
            if match and not invert_match and (fish in baits):
                click.echo(line.strip(), output_file)
            elif match and invert_match and (fish not in baits):
                click.echo(line.strip(), output_file)
            elif not match and not invert_match:
                for bait in baits:
                    if (bait in fish) or (fish in bait):
                        click.echo(line.strip(), output_file)
            elif not match and invert_match:
                for bait in baits:
                    if (bait not in fish) and (fish not in bait):
                        click.echo(line.strip(), output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-b', '--bait-file', 'bait_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input bait file.')
@click.option('-f', '--fish-file', 'fish_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input fish file.')
@click.option('-bc', '--bait-column', 'bait_column',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Column of bait file that used to match the fish file.')
@click.option('-fc', '--fish-column', 'fish_column',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Column of fish file that used to match the bait file.')
@click.option('--match/--contain', default=True, show_default=True,
              help='Match mode.')
@click.option('-v', '--invert-match', 'invert_match',
              is_flag=True, flag_value=True,
              help='Select non-matching lines.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, print results to terminal as stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(bait_file, fish_file, bait_column, fish_column, match, invert_match, output_file):
    """This program works just like fishing, it helps you to get things that you wanted from a target file."""
    main(bait_file, fish_file, bait_column, fish_column, match, invert_match, output_file)


if __name__ == '__main__':
    run()
