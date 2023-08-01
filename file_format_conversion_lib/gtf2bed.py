#!/usr/bin/env python
"""
File: gtf2bed.py
Description: Convert the file format from GTF to BED
Date: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gtf import Gtf
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(gtf_file, out_file):
    content = []
    for line in Gtf(gtf_file).gtf_to_bed():
        if out_file:
            content.append(line)
        else:
            print(line.strip())
    if out_file:
        with open(out_file, 'w') as o:
            o.write(''.join(content))


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--gtf_file', 'gtf_file', help='Input GTF file.')
@click.option('-o', '--output_bed_file', 'output_bed_file',
              help='[optional] Output BED file, if not specified, print results to terminal as stdout.')
def run(gtf_file, output_bed_file):
    """Convert the file format from GTF to BED."""
    main(gtf_file, output_bed_file)


if __name__ == '__main__':
    run()
