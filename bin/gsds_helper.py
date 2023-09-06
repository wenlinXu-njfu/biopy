#!/usr/bin/env python
"""
File: gsds_helper.py
Description: Convert the file format from GFF or GTF to GSDS.
Date: 2022/3/7
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from _io import TextIOWrapper
import click
from Biolib import Gff, Gtf, Displayer, __version__
displayer = Displayer(__file__.split('/')[-1], version=__version__)


def main(in_file, file_format, feature_type, out_file: TextIOWrapper):
    if file_format == 'gff':
        gff_file_obj = Gff(in_file)
        content = gff_file_obj.gff_to_gsds()
        if out_file:
            with out_file as o:
                o.write(content)
        else:
            print(content)
    else:
        gtf_file_obj = Gtf(in_file)
        content = gtf_file_obj.gtf_to_gsds(feature_type)
        if out_file:
            with out_file as o:
                o.write(content)
        else:
            print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--anno_file', 'anno_file', type=click.File('r'), required=True, help='Input GFF or GTF file.')
@click.option('-f', '--format', 'file_format', type=click.Choice(['gff', 'gtf']), default='gff', show_default=True,
              help='Specify input annotation file format.')
@click.option('-t', '--feature_type', 'feature_type',
              type=click.Choice(['gene', 'transcript']), default='transcript', show_default=True,
              help='If input file is GTF, specify feature type.')
@click.option('-o', '--output_file', 'outfile', type=click.File('w'),
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
