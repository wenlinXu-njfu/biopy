#!/usr/bin/env python
"""
File: get_feature_density.py
Description: Get feature density from GFF file.
Date: 2022/10/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gff import Gff
from Biolib.show_info import Displayer


def main(chr_len_file, gff_file, feature, span: int, out_file):
    chr_len_dict = {line.split('\t')[0]: int(line.strip().split('\t')[1]) for line in open(chr_len_file)}
    content = Gff(gff_file).get_feature_density(chr_len_dict, feature, span)
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)
    else:
        print(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-l', '--chr_len_file', 'chr_len_file', help='Input chromosome length file. (Chr_num\\tLength\\n)')
@click.option('-a', '--anno_file', 'anno', help='Input GFF or GTF file.')
@click.option('-f', '--feature_type', 'feature_type', help='Specify feature type. (eg. exon)')
@click.option('-s', '--span', 'span', type=int, default=100000, show_default=True, help='Density statistical span.')
@click.option('-o', '--output_file', 'out',
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(chr_len_file, anno, feature_type, span, out):
    """Get feature density from GFF file."""
    main(chr_len_file, anno, feature_type, span, out)


if __name__ == '__main__':
    run()
