#!/usr/bin/env python
"""
File: seq_extraction_tool.py
Description: Sequence extraction tool (version=3.0)
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from seq_extraction_lib.extract_single_sequence import run as run1
from seq_extraction_lib.bed_extract_seq import run as run2
from seq_extraction_lib.gtf_extract_seq import run as run3
from seq_extraction_lib.gff_extract_seq import run as run4
from seq_extraction_lib.extract_miRNA import run as run5
from circ_tool_kit_lib.extract_circRNA import run as run6
from seq_extraction_lib.KEGG import run as run7
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def seq_extraction_tool():
    """
    Program: Sequence extraction tool\n
    Version: 1.0.0\n
    Contact: WenlinXu \033[1m(wenlinxu.njfu@outlook.com)\033[0m
    """
    pass


seq_extraction_tool.add_command(run1, 'single')
seq_extraction_tool.add_command(run2, 'bed')
seq_extraction_tool.add_command(run3, 'gtf')
seq_extraction_tool.add_command(run4, 'gff')
seq_extraction_tool.add_command(run5, 'miRNA')
seq_extraction_tool.add_command(run6, 'circRNA')
seq_extraction_tool.add_command(run7, 'KEGG')

if __name__ == '__main__':
    seq_extraction_tool()
