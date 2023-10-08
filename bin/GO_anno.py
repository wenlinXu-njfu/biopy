#!/usr/bin/env python
"""
File: GO_anno.py
Description: Preprocess go-basic.obo file
Date: 2022/1/14
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(go_basic_obo_file, out_file):
    content = []
    d = {}
    is_term = None
    for line in open(go_basic_obo_file):
        if line.startswith('data-version'):
            content.append(f'##{line}##ID\tChild_id\tName\tNamespace\tDefinition')
        elif line.startswith('[Term]'):
            is_term = True
        elif line.startswith('id: '):
            if is_term:
                d[1] = f"{line.strip().replace('id: ', '')}"
        elif line.startswith('name: '):
            if is_term:
                d[3] = f"{line.strip().replace('name: ', '')}"
        elif line.startswith('namespace: '):
            if is_term:
                d[4] = f"{line.strip().replace('namespace: ', '')}"
        elif line.startswith('alt_id: '):
            if is_term:
                if 2 in d:
                    d[2] = d[2] + f",{line.strip().replace('alt_id: ', '')}"
                else:
                    d[2] = f"{line.strip().replace('alt_id: ', '')}"
        elif line.startswith('def: '):
            if is_term:
                d[5] = f"{line.strip().replace('def: ', '')}"
        elif not line.strip():
            if d:
                try:
                    content.append(f'{d[1]}\t{d[2]}\t{d[3]}\t{d[4]}\t{d[5]}')
                except KeyError:
                    content.append(f'{d[1]}\tNull\t{d[3]}\t{d[4]}\t{d[5]}')
            is_term = False
            d = {}
    content = '\n'.join(content)
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--obo_file', 'obo_file',
              metavar='<obo file>', required=True,
              help='Input GO annotation file (go-basic.obo).')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file>',
              help='Output file (ID\\tChild_id\\tName\\tNamespace\\tDefinition), '
                   'if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(obo_file, outfile):
    """Preprocess go-basic.obo file."""
    main(obo_file, outfile)


if __name__ == '__main__':
    run()
