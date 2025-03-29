#!/usr/bin/env python
"""
File: KEGG_anno.py
Description: Preprocess xxx00001.keg file.
CreateDate: 2022/6/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import requests
from json import loads
from re import sub
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.3.0')


def download_json(species_name: str):
    response = requests.get(url=f'https://www.kegg.jp/kegg-bin/download_htext?htext={species_name}00001&format=json&filedir=kegg/brite/{species_name}',
                            headers={
                                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"
                            })
    return response.text


def parse_json(content: str,
               out_file: TextIOWrapper):
    items = loads(content)['children']
    for level_1 in items:
        level_1_name = level_1['name'].split(' ')
        level_1_ko = f'ko{level_1_name[0]}'
        level_1_name = ' '.join(level_1_name[1:])
        level_1_children = level_1['children']
        for level_2 in level_1_children:
            level_2_name = level_2['name'].split(' ')
            level_2_ko = f'ko{level_2_name[0]}'
            level_2_name = ' '.join(level_2_name[1:])
            level_2_children = level_2['children']
            for level_3 in level_2_children:
                level_3_name = level_3['name'].split(' ')
                level_3_ko = f'ko{level_3_name[0]}'
                level_3_name = sub(r'\[.*]', '', ' '.join(level_3_name[1:]))
                try:
                    level_3_children = level_3['children']
                except KeyError:
                    pass
                else:
                    for i in level_3_children:
                        name = i['name'].split('\t')
                        K = name[1].split(' ')[0]
                        id = 'pop:' + name[0].split(' ')[0] + '\t' + K + '\t' + ' '.join(name[0].split(' ')[1:])
                        click.echo(f'{id}\t{level_1_ko}\t{level_1_name}\t{level_2_ko}\t{level_2_name}\t{level_3_ko}\t{level_3_name}', out_file)


def main(species_name: str,
         out_file: TextIOWrapper):
        content = download_json(species_name=species_name)
        parse_json(content=content, out_file=out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-n', '--species-name', 'species_name',
              metavar='<str>', required=True,
              help='Species name. (eg: pop)')
@click.option('-o', '--output-file', 'outfile',
              metavar='<file|stdout>', type=click.File('w'),
              help=r'Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(species_name, outfile):
    """Preprocess xxx00001.keg file"""
    main(species_name, outfile)


if __name__ == '__main__':
    run()
