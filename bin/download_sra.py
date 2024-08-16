#!/usr/bin/env python
"""
File: download_sra.py
Description: Download SRA data from SRA database and convert it to fq.
CreateDate: 2024/7/22
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from os import getcwd, makedirs
from os.path import abspath
import click
from pybioinformatic import TaskManager, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(sra_list: TextIOWrapper,
         num_processing: int,
         output_path: str):
    output_path = abspath(output_path)
    makedirs(f'{output_path}/fq', exist_ok=True)
    cmds = []
    for line in sra_list:
        sra_accession = line.strip()
        cmd = (f'prefetch {sra_accession} -O {output_path}; '
               f'fastq-dump --gzip --split-3 {output_path}/{sra_accession}/{sra_accession}.sra -O {output_path}/fq')
        cmds.append(cmd)
    tkm = TaskManager(commands=cmds, num_processing=num_processing)
    tkm.parallel_run_cmd()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-l', '--sra-list', 'sra_list',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Input sra accession id file. (one accession per line)')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=1, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(sra_list, num_processing, output_path):
    """Download SRA data from SRA database and convert it to fq."""
    main(sra_list, num_processing, output_path)


if __name__ == '__main__':
    run()
