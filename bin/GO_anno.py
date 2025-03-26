#!/usr/bin/env python
"""
File: GO_anno.py
Description: GO annotate.
CreateDate: 2022/1/14
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import StringIO, TextIOWrapper
from os import getcwd, makedirs
from re import findall
from pandas import read_table
import click
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def parse_obo(go_basic_obo_file: str):
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
    return content


def main(anno_file: TextIOWrapper,
         go_basic_obo_file: str = None,
         output_path: str = getcwd()):
    makedirs(output_path, exist_ok=True)
    if not go_basic_obo_file:
        cmd = f'wget -c -O {output_path}/go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo'
        tkm = TaskManager(num_processing=1)
        tkm.echo_and_exec_cmd(cmd=cmd, show_cmd=True)
    content = parse_obo(go_basic_obo_file) if go_basic_obo_file else parse_obo(f'{output_path}/go-basic.obo')
    df = read_table(StringIO(content), comment='#', header=None, names=['ID', 'Child_id', 'Name', 'Namespace', 'Definition'])
    df.set_index(keys='ID', drop=True, inplace=True)
    with open(f'{output_path}/GO_anno.xls', 'w') as o:
        for line in anno_file:
            query_id = line.strip().split('\t')[0]
            GO_ids = findall(r'GO:\d{7}', line)
            if GO_ids:
                for GO_id in GO_ids:
                    try:
                        definition = df.loc[GO_id, 'Name']
                    except KeyError:
                        pass
                    else:
                        click.echo(f'{query_id}\t{GO_id}\t{definition}', o)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('anno_file', metavar='<anno file | stdin>', nargs=1, required=True, type=click.File('r'))
@click.option('-i', '--obo-file', 'obo_file',
              metavar='<obo file>',
              help='Input GO annotation file (go-basic.obo). If not specified, download from url automatically.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help=r'Output path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(obo_file, anno_file, output_path):
    """GO annotate."""
    main(
        anno_file=anno_file,
        go_basic_obo_file=obo_file,
        output_path=output_path
    )


if __name__ == '__main__':
    run()
