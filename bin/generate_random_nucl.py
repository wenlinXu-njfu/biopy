#!/usr/bin/env python
"""
File: generate_random_nucl.py
Description: Generate random nucleotide sequence.
CreateDate: 2025/5/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Union, List
import click
from pybioinformatic import Nucleotide, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def check_params(
    names: Union[List[str], str, None],
    length: Union[List[int], int, None]
):
    if all([isinstance(names, list), isinstance(length, list)]) and len(names) != len(length):
        click.echo('\033[31mError: The name and the quantity of the length do not match.\033[0m', err=True)
        exit()


def main(names: Union[List[str], str, None],
         length: Union[List[int], str, None],
         bias: Union[List[float], str, None],
         output_file: TextIOWrapper = None):
    check_params(names, length)
    if isinstance(names, list) and isinstance(length, list):
        for i in range(len(names)):
            click.echo(Nucleotide.random_nucl(name=names[i], length=length[i], bias=bias), output_file)
    else:
        click.echo(Nucleotide.random_nucl(name=names, length=length, bias=bias), output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-n', '--names', 'names', metavar='<str>',
              help='Nucleotide sequence name. Multiple names should be separated by commas (eg. seq1,se2),'
                   'and the number of names should be consistent with the "-l, --length" option values.')
@click.option('-l', '--length', 'length', metavar='<str>',
              help='Nucleotide sequence length. Multiple length should be separated by commas,'
                   'and the number of length should be consistent with the "-n, --names" option values (eg. 100,200).')
@click.option('-b', '--bias', 'bias', metavar='<str>',
              help='Nucleotide sequence bias of AGCT. Multiple bias should be separated by commas (eg. 0.25,0.5,0.1,0.15).')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(names, length, bias, output_file):
    """Generate random nucleotide sequence."""
    if names is not None and ',' in names:
        names = names.split(',')
    if length is not None and ',' in length:
        length = [int(i) for i in length.split(',')]
    elif isinstance(length, str):
        length = int(length)
    if bias is not None and ',' in bias:
        bias = [float(i) for i in bias.split(',')]
    main(names, length, bias, output_file)


if __name__ == '__main__':
    run()
