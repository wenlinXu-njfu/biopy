#!/usr/bin/env python
"""
File: file_format_conversion.py
Description: File format conversion tool (version=3.0)
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from file_format_conversion_lib.gff_sort import run as run1
from file_format_conversion_lib.gff2gtf import run as run2
from file_format_conversion_lib.gff2bed import run as run3
from file_format_conversion_lib.gtf2bed import run as run4
from file_format_conversion_lib.fasta_conversion import run as run5
from file_format_conversion_lib.fq2fa import run as run6


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def file_format_conversion():
    """
    Program: File format conversion tool\n
    Version: 1.0.0\n
    Contact: WenlinXu \033[1m(wenlinxu.njfu@outlook.com)\033[0m
    """
    pass


file_format_conversion.add_command(run1, 'gff_sort')
file_format_conversion.add_command(run2, 'gff2gtf')
file_format_conversion.add_command(run3, 'gff2bed')
file_format_conversion.add_command(run4, 'gtf2bed')
file_format_conversion.add_command(run5, 'fasta')
file_format_conversion.add_command(run6, 'fq2fa')


if __name__ == '__main__':
    file_format_conversion()
