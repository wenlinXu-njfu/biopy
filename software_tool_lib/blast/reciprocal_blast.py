#!/usr/bin/env python
"""
File: reciprocal_blast.py
Description: By reciprocal blast, obtain sequence pair that best match each other
Date: 2022/4/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.blast import Blast
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(blast1, blast2, top, out_file):
    pair_dict1 = Blast(blast1).get_pair_dict(top)  # {query1: {sbject1: [], sbject2: [], ...}, query2: {}, ...}
    pair_dict2 = Blast(blast2).get_pair_dict(top)  # {query1: {sbject1: [], sbject2: [], ...}, query2: {}, ...}
    content = ''
    for query, d1 in pair_dict1.items():
        for sbject, info in d1.items():
            try:
                queries = list(pair_dict2[sbject].keys())  # queries that sbject best align
            except KeyError:
                pass
            else:
                if query in queries:
                    if out_file:
                        content += f'{query}\t{sbject}\n'
                    else:
                        print(f'{query}\t{sbject}')
    if out_file and content:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--blast_result1', 'blast_result1', help='Input blast result file.')
@click.option('-I', '--blast_result2', 'blast_result2', help='Input another blast result file.')
@click.option('-t', '--top', 'top', type=int, default=3,
              help='[optional] Specify max alignment num of each sequence. {default: 3}')
@click.option('-o', '--output_file', 'output_file',
              help='[optional] Output file, if not specified, print result to terminal as stdout.')
def run(blast_result1, blast_result2, top, output_file):
    """By reciprocal blast, obtain sequence pair that best match each other."""
    main(blast_result1, blast_result2, top, output_file)


if __name__ == '__main__':
    run()
