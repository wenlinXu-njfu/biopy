#!/usr/bin/env python
"""
File: file_format_conversion.py
Description: File format conversion tool.
CreateDate: 2022/3/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from file_format_conversion_lib import fq2fa, vcf2gt, embl2fa, __version__
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version=__version__)


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def file_format_conversion():
    """File format conversion tool."""
    pass


file_format_conversion.add_command(fq2fa, 'fq2fa')
file_format_conversion.add_command(vcf2gt, 'vcf2gt')
file_format_conversion.add_command(embl2fa, 'embl2fa')


if __name__ == '__main__':
    file_format_conversion()
