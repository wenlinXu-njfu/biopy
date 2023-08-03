#!/usr/bin/env python
"""
File: CircToolKit.py
Description: CircRNAs analysis tools
Date: 2022/4/3
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from circ_tool_kit_lib.extract_circRNA import run as run1
from circ_tool_kit_lib.get_circ_exp import run as run2
from circ_tool_kit_lib.Tau_index import run as run3
from circ_tool_kit_lib.reverse_complementary_analysis import run as run4
from circ_tool_kit_lib.repeat_seq_analysis import run as run5
from circ_tool_kit_lib.alternative_cyclization_analysis import run as run6
from circ_tool_kit_lib.circular_translation import run as run7
from circ_tool_kit_lib.ceRNA_identification import run as run8


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def CircToolKit():
    """
    Program: CircRNAs analysis tools\n
    Version: 1.0.0\n
    Contact: WenlinXu \033[1m(wenlinxu.njfu@outlook.com)\033[0m
    """
    pass


CircToolKit.add_command(run1, 'get_circ_seq')
CircToolKit.add_command(run2, 'get_circ_exp')
CircToolKit.add_command(run3, 'Tau_index')
CircToolKit.add_command(run4, 'rca')
CircToolKit.add_command(run5, 'rsa')
CircToolKit.add_command(run6, 'aca')
CircToolKit.add_command(run7, 'translation')
CircToolKit.add_command(run8, 'ceRNA_identification')

if __name__ == '__main__':
    CircToolKit()
