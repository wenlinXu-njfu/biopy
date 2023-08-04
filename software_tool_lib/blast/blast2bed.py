#!/usr/bin/env python
"""
File: blast2bed.py
Description: Convert balst result to BED format
Date: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.blast import Blast
from Biolib.show_info import Displayer


def main(blast_file, out_file):
    blast_file_obj = Blast(blast_file)
    content = blast_file_obj.blast_to_BED()
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--blast_result', 'blast_result_file', help='Input blast result file with format 6.')
@click.option('-o', '--output_file', 'output_file',
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(blast_result_file, output_file):
    """Convert balst result to BED format."""
    main(blast_result_file, output_file)


if __name__ == '__main__':
    run()
