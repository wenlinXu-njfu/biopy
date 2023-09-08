#!/usr/bin/env python
"""
File: get_feature_density.py
Description: Get feature density from GFF file.
Date: 2022/10/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
import click
from Biolib import Gff, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(chr_len_file: TextIOWrapper,
         gff_file: TextIOWrapper,
         feature: str,
         span: int,
         output_file: TextIOWrapper):
    chr_len_dict = {line.split('\t')[0]: int(line.strip().split('\t')[1]) for line in chr_len_file}
    for line in Gff(gff_file).get_feature_density(chr_len_dict, feature, span):
        click.echo(line, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-l', '--chr_len_file', 'chr_len_file', type=click.File('r'), required=True,
              help='Input chromosome length file. (Chr_num\\tLength\\n)')
@click.option('-g', '--gff_file', 'gff_file', type=click.File('r'), required=True, help='Input genome annotation GFF file.')
@click.option('-f', '--feature_type', 'feature_type', required=True, help='Specify feature type. (eg. exon)')
@click.option('-s', '--span', 'span', type=int, default=100000, show_default=True, help='Density statistical span.')
@click.option('-o', '--output_file', 'out', type=click.File('w'),
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(chr_len_file, gff_file, feature_type, span, out):
    """Get feature density from GFF file."""
    main(chr_len_file, gff_file, feature_type, span, out)


if __name__ == '__main__':
    run()
