#!/usr/bin/env python
"""
File: seq_extraction_tool.py
Description: Sequence extraction tool.
Date: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from circ_tool_kit_lib import extract_circrna
from seq_extraction_lib import (extract_single_sequence, bed_extract_seq, gtf_extract_seq, gff_extract_seq,
                                extract_mirna, kegg, __version__)
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version=__version__)


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def seq_extraction_tool():
    """Sequence extraction tool."""
    pass


seq_extraction_tool.add_command(extract_single_sequence, 'single')
seq_extraction_tool.add_command(bed_extract_seq, 'bed')
seq_extraction_tool.add_command(gtf_extract_seq, 'gtf')
seq_extraction_tool.add_command(gff_extract_seq, 'gff')
seq_extraction_tool.add_command(extract_mirna, 'miRNA')
seq_extraction_tool.add_command(extract_circrna, 'circRNA')
seq_extraction_tool.add_command(kegg, 'KEGG')

if __name__ == '__main__':
    seq_extraction_tool()
