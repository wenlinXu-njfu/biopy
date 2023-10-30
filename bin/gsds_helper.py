#!/usr/bin/env python
"""
File: gsds_helper.py
Description: Convert the file format from GFF or GTF to GSDS.
CreateDate: 2022/3/7
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
import click
from pybioinformatic import Gff, Gtf, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(anno_file: Union[str, TextIOWrapper],
         file_format: click.Choice(['gff', 'gtf']),
         feature_type: click.Choice(['gene', 'transcript']),
         out_file: TextIOWrapper = None):
    if file_format == 'gff':
        with Gff(anno_file) as gff:
            for line in gff.to_gsds():
                click.echo(line, out_file)
    else:
        with Gtf(anno_file) as gtf:
            for line in gtf.to_gsds(feature_type):
                click.echo(line, out_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--anno_file', 'anno_file', metavar='<file>',
              type=click.File('r'), required=True, help='Input GFF or GTF file.')
@click.option('-f', '--format', 'file_format', metavar='<gff|gtf>',
              type=click.Choice(['gff', 'gtf']), default='gff', show_default=True,
              help='Specify input annotation file format.')
@click.option('-t', '--feature_type', 'feature_type', metavar='<gene|transcript>',
              type=click.Choice(['gene', 'transcript']), default='transcript', show_default=True,
              help='If input file is GTF, specify feature type.')
@click.option('-o', '--output_file', 'outfile', metavar='<file>', type=click.File('w'),
              help='Output file (ID\\tStart\\tEnd\\tFeature\\tFrame), '
                   'if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(anno_file, file_format, feature_type, outfile):
    """Convert the file format from GFF or GTF to GSDS."""
    if file_format == 'gff' and feature_type:
        click.echo('\033[33mWarning: input file is GFF, "-t --feature_type" option is invalid.\033[0m', err=True)
    main(anno_file, file_format, feature_type, outfile)


if __name__ == '__main__':
    run()
