#!/usr/bin/env python
"""
File: file_split.py
Description: Divide a large file into several smaller files
Date: 2022/1/23
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import mkdir
from os.path import exists
import click


def file_split(in_file: str, sub_file_line_num: int, out_dir: str, header: bool = True, header_line_num: int = 1):
    if not exists(out_dir):
        mkdir(out_dir)
    line_num = sum(1 for _ in open(in_file))
    line_count = 0
    subfile_num = 1
    header_content = content = ''
    out_prefix = in_file.split('/')[-1]
    for line in open(in_file):
        line_count += 1
        if header:
            if line_count == line_num:
                content += line
                with open(f"{out_dir}/{out_prefix}.part{subfile_num}", 'w') as o:
                    o.write(header_content + content)
            elif line_count <= header_line_num:
                header_content += line
            else:
                if (line_count - header_line_num) % sub_file_line_num != 0:
                    content += line
                else:
                    content += line
                    with open(f"{out_dir}/{out_prefix}.part{subfile_num}", 'w') as o:
                        o.write(header_content + content)
                    content = ''
                    subfile_num += 1
        else:
            if line_count == line_num:
                content += line
                with open(f"{out_dir}/{out_prefix}.part{subfile_num}", 'w') as o:
                    o.write(content)
            elif line_count % sub_file_line_num != 0:
                content += line
            else:
                content += line
                with open(f"{out_dir}/{out_prefix}.part{subfile_num}", 'w') as o:
                    o.write(content)
                content = ''
                subfile_num += 1


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', 'input_file', help='Input file.')
@click.option('-n', '--line_num', 'line_num', type=int, help='The line num of each sub file, not including header.')
@click.option('-H', '--header', is_flag=True, flag_value=True,
              help='If specified header, each sub file will contain header content.')
@click.option('-N', '--header_num', 'header_num', type=int, default=1, show_default=True,
              help='If header specified, specify the number of header line.')
@click.option('-o', '--output_dir', 'output_dir',
              help='Output directory, if not exists, it will be created automatically.')
def run(input_file, line_num, header, header_num, output_dir):
    """Divide a large file into several smaller files."""
    file_split(input_file, line_num, output_dir, header, header_num)


if __name__ == '__main__':
    run()
