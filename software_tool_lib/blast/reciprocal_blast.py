#!/usr/bin/env python
"""
File: reciprocal_blast.py
Description: By reciprocal blast, obtain sequence pair that best match each other.
CreateDate: 2022/4/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Blast, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(blast1: Union[str, TextIOWrapper],
         blast2: Union[str, TextIOWrapper],
         top: int,
         out_file: TextIOWrapper = None):
    with Blast(blast1) as blast_obj1:
        pair_dict1 = blast_obj1.get_pair_dict(top)  # {query1: {sbject1: [], sbject2: [], ...}, query2: {}, ...}
    with Blast(blast2) as blast_obj2:
        pair_dict2 = blast_obj2.get_pair_dict(top)  # {query1: {sbject1: [], sbject2: [], ...}, query2: {}, ...}
    for query, d1 in pair_dict1.items():
        for sbject, info in d1.items():
            try:
                queries = list(pair_dict2[sbject].keys())  # queries that sbject best align
            except KeyError:
                pass
            else:
                if query in queries:
                    click.echo(f'{query}\t{sbject}', out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--blast_result1', 'blast_result1',
              metavar='<blast file|stdin>', type=click.File('r'), required=True,
              help='Input blast result file.')
@click.option('-I', '--blast_result2', 'blast_result2',
              metavar='<blast file|stdin>', type=click.File('r'), required=True,
              help='Input another blast result file.')
@click.option('-t', '--top', 'top',
              metavar='<int>', type=int, default=3, show_default=True,
              help='Specify max alignment num of each sequence.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file|stdout>', type=click.File('w'),
              help='Output file, if not specified, print result to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(blast_result1, blast_result2, top, output_file):
    """By reciprocal blast, obtain sequence pair that best match each other."""
    main(blast_result1, blast_result2, top, output_file)


if __name__ == '__main__':
    run()
