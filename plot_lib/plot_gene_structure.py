#!/usr/bin/env python
"""
File: plot_gene_structure.py
Description: Plot gene structure based on annotation file
Date: 2022/3/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from plot_lib.gene_structure.gff import plot_mRNA_structure
from plot_lib.gene_structure.gtf import plot_gene_structure
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(in_file: str, in_file_format: click.Choice(['gff', 'gtf']),
         utr_color: str, cds_color: str, exon_color: str, edge_color: str, out_path: str,
         figure_width: float, figure_height: float,
         out_format: click.Choice(['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']),
         utr_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None,
         cds_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None,
         exon_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None):
    if in_file_format == 'gff':
        plot_mRNA_structure(in_file, utr_color, cds_color, edge_color, figure_width, figure_height, out_path,
                            out_suffix=out_format, utr_hatch=utr_hatch, cds_hatch=cds_hatch)
    else:
        plot_gene_structure(in_file, exon_color, edge_color=edge_color,
                            figure_width=figure_width, figure_height=figure_height,
                            out_path=out_path, out_suffix=out_format, exon_hatch=exon_hatch)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', 'input_file', help='Input file.')
@click.option('-t', '--format', 't', type=click.Choice(['gff', 'gtf']), default='gff', show_default=True,
              help='Specify the format of input file.')
@click.option('-u', '--utr_color', 'utr_color', default='salmon', show_default=True,
              help='[optional] If input GFF file, specify color of utr, it supports color code.')
@click.option('-uh', '--utr_hatch', 'utr_hatch',
              type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']), default=None, show_default=True,
              help='[optional] If input GFF file, specify hatch of utr.')
@click.option('-c', '--cds_color', 'cds_color', default='skyblue', show_default=True,
              help='[optional] If input GFF file, specify color of utr, it supports color code.')
@click.option('-ch', '--cds_hatch', 'cds_hatch',
              type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']), default=None, show_default=True,
              help='[optional] If input GFF file, specify hatch of cds.')
@click.option('-e', '--exon_color', 'exon_color', default='salmon', show_default=True,
              help='[optional] If input GTF file, specify color of exon, it supports color code.')
@click.option('-eh', '--exon_hatch', 'exon_hatch',
              type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']), default=None, show_default=True,
              help='[optional] If input GTF file, specify hatch of exon.')
@click.option('-E', '--edge_color', 'edge_color', default=None, show_default=True, help='[optional] Set edge color.')
@click.option('-W', '--figure_width', 'figure_width', type=float, default=20.0, show_default=True,
              help='Output figure width.')
@click.option('-H', '--figure_height', 'figure_height', type=float, default=10.0, show_default=True,
              help='Output figure height.')
@click.option('-o', '--output_path', 'output_path', default='./', show_default=True,
              help='Output file path. {default: ./}')
@click.option('-O', '--output_format', 'output_format', default='pdf', show_default=True,
              type=click.Choice(['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']),
              help='Output file format.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, t, utr_color, utr_hatch, cds_color, cds_hatch, exon_color, exon_hatch, edge_color,
        figure_width, figure_height, output_path, output_format):
    """Plot gene structure based on annotation file."""
    if t == 'gff':
        if exon_hatch or exon_color != 'salmon':
            click.echo('\033[33mThere are conflicting options, ignore "exon_color" and "exon_hatch".\033[0m')
    elif t == 'gtf':
        if utr_color != 'salmon' or utr_hatch or cds_color != 'skyblue' or cds_hatch:
            click.echo('\033[33mThere are conflicting options, ignore "utr_color", "utr_hatch", "cds_color" and '
                       '"cds_hatch".\033[0m')
    main(input_file, t, utr_color, cds_color, exon_color, edge_color, output_path, figure_width, figure_height,
         output_format, utr_hatch, cds_hatch, exon_hatch)


if __name__ == '__main__':
    run()
