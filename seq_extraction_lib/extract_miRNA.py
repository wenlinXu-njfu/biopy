#!/usr/bin/env python
"""
File: extract_miRNA.py
Description: Extract miRNA sequence from miRNA GFF file
Date: 2022/3/24
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import mkdir
from os.path import exists
import click
from Biolib import Gff, Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(miRNA_gff_file, out_dir):
    if out_dir and not exists(out_dir):
        mkdir(out_dir)
    elif not out_dir:
        out_dir = './'
    gff_file_obj = Gff(miRNA_gff_file)
    primary = star = mature = ''
    for nucl_obj in gff_file_obj.miRNA_extraction():
        if 'star' in nucl_obj.id:
            star += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
        elif 'mature' in nucl_obj.id:
            mature += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
        else:
            primary += f">{nucl_obj.id}\n{nucl_obj.seq}\n"
    with open(f"{out_dir}/primary.fa", 'w') as o1:
        o1.write(primary)
    with open(f"{out_dir}/star.fa", 'w') as o2:
        o2.write(star)
    with open(f"{out_dir}/mature.fa", 'w') as o3:
        o3.write(mature)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gff_file', 'gff_file', required=True, help='Input miRNA GFF file.')
@click.option('-o', '--output_dir', 'output_dir', default='./', show_default=True,
              help='Output directory, if the output directory does not exist, it will be created automatically.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(gff_file, output_dir):
    """Extract miRNA sequence from miRNA GFF file"""
    main(gff_file, output_dir)


if __name__ == '__main__':
    run()
