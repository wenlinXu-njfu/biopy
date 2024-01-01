#!/usr/bin/env python
"""
File: batch_pfamscan.py
Description: Batch pfamscan with multiple FASTA files
CreateDate: 2022/4/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import listdir, mkdir
from os.path import exists
import click
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(input_dir: str,
         pfamscan_database: str,
         num_processing: int,
         out_dir: str):
    if not exists(out_dir):
        mkdir(out_dir)
    files = listdir(input_dir)
    cmds = []
    for file in files:
        cmds.append(f"pfam_scan.pl -fasta {input_dir}/{file} "
                    f"-dir {pfamscan_database} "
                    f"-outfile {out_dir}/{file}_pfamsacn_out.txt")
    tkm = TaskManager(commands=cmds, num_processing=num_processing)
    tkm.parallel_run_cmd()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_dir', 'fasta_dir',
              metavar='<dir>', required=True,
              help='Input FASTA file directory.')
@click.option('-d', '--database', 'pfam_database',
              metavar='<path>', required=True,
              help='PfamScan database path.')
@click.option('-o', '--output_dir', 'output_dir',
              metavar='<dir>', default='./', show_default=True,
              help='Output directory, if not exits, create automatically.')
@click.option('-n', '--num_processing', 'num_processing',
              )
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fasta_dir, pfam_database, num_processing, output_dir):
    """Batch pfamscan with multiple FASTA files."""
    main(fasta_dir, pfam_database, num_processing, output_dir)


if __name__ == '__main__':
    run()
