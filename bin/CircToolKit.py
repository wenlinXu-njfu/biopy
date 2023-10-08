#!/usr/bin/env python
"""
File: CircToolKit.py
Description: CircRNAs analysis tools.
Date: 2022/4/3
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from pybioinformatic import Displayer
from circ_tool_kit_lib import (extract_circrna, get_circ_exp, tau_index, reverse_complementary_analysis,
                               repeat_seq_analysis, alternative_cyclization_analysis, circular_translation,
                               cerna_identification, __version__)
displayer = Displayer(__file__.split('/')[-1], version=__version__)


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def CircToolKit():
    """CircRNAs analysis tools."""
    pass


CircToolKit.add_command(extract_circrna, 'get_circ_seq')
CircToolKit.add_command(get_circ_exp, 'get_circ_exp')
CircToolKit.add_command(tau_index, 'Tau_index')
CircToolKit.add_command(reverse_complementary_analysis, 'rca')
CircToolKit.add_command(repeat_seq_analysis, 'rsa')
CircToolKit.add_command(alternative_cyclization_analysis, 'aca')
CircToolKit.add_command(circular_translation, 'translation')
CircToolKit.add_command(cerna_identification, 'ceRNA_identification')

if __name__ == '__main__':
    CircToolKit()
