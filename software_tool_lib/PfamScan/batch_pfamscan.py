#!/usr/bin/env python
"""
File: batch_pfamscan.py
Description: Batch pfamscan with multiple FASTA files
Date: 2022/4/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import listdir, mkdir, system
from os.path import exists
import click
from tqdm import tqdm
from Biolib.show_info import Displayer


def main(in_dir, pfamscan_database, out_dir):
    if not exists(out_dir):
        mkdir(out_dir)
    files = listdir(in_dir)
    with tqdm(total=len(files)) as pbar:
        for file in files:
            system(command=f"pfam_scan.pl -fasta {in_dir}/{file} -dir {pfamscan_database} "
                              f"-outfile {out_dir}/{file}_pfamsacn_out.txt")
            pbar.update()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta_dir', 'fasta_dir', help='Input FASTA file directory.')
@click.option('-d', '--database', 'pfam_database', help='Input PfamScan database.')
@click.option('-o', '--output_dir', 'output_dir', default='./', show_default=True,
              help='Output directory, if not exits, create automatically.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(fasta_dir, pfam_database, output_dir):
    """Batch pfamscan with multiple FASTA files."""
    main(fasta_dir, pfam_database, output_dir)


if __name__ == '__main__':
    run()
