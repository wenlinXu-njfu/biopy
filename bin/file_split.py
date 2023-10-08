#!/usr/bin/env python
"""
File: file_split.py
Description: Divide a large file into several smaller files.
Date: 2022/1/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import mkdir
from os.path import exists
from typing import Union
from _io import TextIOWrapper
import click
from click._compat import _NonClosingTextIOWrapper
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(input_file: Union[TextIOWrapper, list],
         sub_file_line_num: int,
         output_dir: str,
         header: bool,
         header_line_num: int):
    if not exists(output_dir):
        mkdir(output_dir)
    line_count, subfile_num = 0, 1
    line_num = sum(1 for _ in input_file)
    header_content = content = ''
    # Read content from a file.
    if isinstance(input_file, TextIOWrapper):
        output_prefix = f"{input_file.name.split('/')[-1]}."
        input_file = open(input_file.name)
    # Read content from stdout.
    else:
        output_prefix = ''
    for line in input_file:
        line_count += 1
        if header:
            if line_count == line_num:
                content += line
                with open(f"{output_dir}/{output_prefix}part{subfile_num}", 'w') as o:
                    o.write(header_content + content)
            elif line_count <= header_line_num:
                header_content += line
            else:
                if (line_count - header_line_num) % sub_file_line_num != 0:
                    content += line
                else:
                    content += line
                    with open(f"{output_dir}/{output_prefix}part{subfile_num}", 'w') as o:
                        o.write(header_content + content)
                    content = ''
                    subfile_num += 1
        else:
            if line_count == line_num:
                content += line
                with open(f"{output_dir}/{output_prefix}part{subfile_num}", 'w') as o:
                    o.write(content)
            elif line_count % sub_file_line_num != 0:
                content += line
            else:
                content += line
                with open(f"{output_dir}/{output_prefix}part{subfile_num}", 'w') as o:
                    o.write(content)
                content = ''
                subfile_num += 1


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', 'input_file',
              metavar='<file>', type=click.File('r'), required=True,
              help='Input file.')
@click.option('-n', '--line_num', 'line_num',
              metavar='<int>', type=int,  required=True,
              help='The line num of each sub file, not including header.')
@click.option('-H', '--header', is_flag=True, flag_value=True,
              help='If specified header, each sub file will contain header content.')
@click.option('-N', '--header_num', 'header_num',
              metavar='<int>', type=int, default=1, show_default=True,
              help='If header specified, specify the number of header line.')
@click.option('-o', '--output_dir', 'output_dir',
              metavar='<dir>', required=True,
              help='Output directory, if not exists, it will be created automatically.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, line_num, header, header_num, output_dir):
    """Divide a large file into several smaller files."""
    if isinstance(input_file, _NonClosingTextIOWrapper):
        input_file = click.open_file('-').readlines()
    main(input_file, line_num, output_dir, header, header_num)


if __name__ == '__main__':
    run()
