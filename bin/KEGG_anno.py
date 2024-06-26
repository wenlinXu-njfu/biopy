#!/usr/bin/env python
"""
File: KEGG_anno.py
Description: Preprocess xxx00001.keg file.
CreateDate: 2022/6/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from re import findall
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(keg_file: TextIOWrapper,
         out_file: TextIOWrapper):
    ko = pathway = species = ''
    for line in keg_file:
        line = line.strip()
        if not line.startswith('A09160'):
            if line.startswith('C'):
                ko = findall(r'[a-z]{3}\d{5}', line)
                if ko:
                    ko = ko[0]
                    species = findall(r'[a-z]{3}', ko)[0]
                    ko_num = findall(r'\d{5}', ko)[0]
                    ko = f"ko{ko_num}"
                    pathway = ' '.join(line.split('    ')[1].split(' ')[1:])
                else:
                    ko = f"ko{line.split('    ')[1].split(' ')[0]}"
                    pathway = ' '.join(line.split('    ')[1].split(' ')[1:])
                if '[' in pathway:
                    _id = findall(r' \[.+]', pathway)[0]
                    pathway = pathway.replace(_id, '')
            elif line.startswith('D'):
                gene_id = f"{species}{line.split('      ')[1].split(' ')[0]}"
                K = line.split('      ')[1].split('\t')[1].split(' ')[0]
                EC = findall(r' \[EC:.+]', line)
                if EC:
                    EC = EC[0]
                    KO = line.split('      ')[1].split('\t')[1].split('; ')[1].replace(EC, '')
                else:
                    KO = line.split('      ')[1].split('\t')[1].split('; ')[1]
                click.echo(f'{gene_id}\t{K}\t{KO}\t{ko}\t{pathway}', out_file)
        else:
            break


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--htext_file', 'htext_file',
              metavar='<htext file|stdin>', required=True,
              help='Input kegg file. (format: htext)')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file>',
              help=r'Output file (eg. pop7465650\tK18835\tWRKY transcription factor 2\tko04626\tPlant-pathogen interaction), '
                   'stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(htext_file, outfile):
    """Preprocess xxx00001.keg file"""
    main(htext_file, outfile)


if __name__ == '__main__':
    run()
