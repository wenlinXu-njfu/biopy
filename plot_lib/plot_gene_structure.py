#!/usr/bin/env python
"""
File: plot_gene_structure.py
Description: Plot gene structure based on annotation file
CreateDate: 2022/3/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from plot_lib.gene_structure.gff import plot_mRNA_structure
from plot_lib.gene_structure.gtf import plot_gene_structure
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(in_file: str,
         in_file_format: click.Choice(['gff', 'gtf']),
         utr_color: str,
         cds_color: str,
         exon_color: str,
         edge_color: str,
         figure_width: float,
         figure_height: float,
         output_file: str,
         utr_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None,
         cds_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None,
         exon_hatch: click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']) = None):
    if in_file_format == 'gff':
        plot_mRNA_structure(gff_file=in_file,
                            utr_color=utr_color,
                            cds_color=cds_color,
                            edge_color=edge_color,
                            figure_width=figure_width,
                            figure_height=figure_height,
                            out_file=output_file,
                            utr_hatch=utr_hatch,
                            cds_hatch=cds_hatch)
    else:
        plot_gene_structure(gtf_file=in_file,
                            exon_color=exon_color,
                            edge_color=edge_color,
                            figure_width=figure_width,
                            figure_height=figure_height,
                            out_file=output_file,
                            exon_hatch=exon_hatch)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input-file', 'input_file',
              metavar='<anno file>', required=True,
              help='Input GFF or GTF file.')
@click.option('-t', '--file-format', 'file_format',
              metavar='<gff|gtf>', type=click.Choice(['gff', 'gtf']), default='gff', show_default=True,
              help='Specify the format of input file.')
@click.option('-u', '--utr-color', 'utr_color',
              metavar='<str>', default='salmon', show_default=True,
              help='If input GFF file, specify color of utr, it supports color code.')
@click.option('-uh', '--utr-hatch', 'utr_hatch',
              metavar='<str>', type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']),
              help='If input GFF file, specify hatch of utr.')
@click.option('-c', '--cds-color', 'cds_color',
              metavar='<str>', default='skyblue', show_default=True,
              help='If input GFF file, specify color of utr, it supports color code.')
@click.option('-ch', '--cds-hatch', 'cds_hatch',
              metavar='<str>', type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']),
              help='If input GFF file, specify hatch of cds.')
@click.option('-e', '--exon-color', 'exon_color',
              metavar='<str>', default='salmon', show_default=True,
              help='[optional] If input GTF file, specify color of exon, it supports color code.')
@click.option('-eh', '--exon-hatch', 'exon_hatch',
              metavar='<str>', type=click.Choice(['/', '|', '\\', '+', '-', 'x', '*', 'o', 'O', '.']),
              help='If input GTF file, specify hatch of exon.')
@click.option('-ec', '--edge-color', 'edge_color',
              metavar='<str>', default='black', show_default=True,
              help='Set edge color.')
@click.option('-width', '--figure-width', 'figure_width',
              metavar='<float>', type=float, default=20.0, show_default=True,
              help='Output figure width.')
@click.option('-height', '--figure-height', 'figure_height',
              metavar='<float>', type=float, default=10.0, show_default=True,
              help='Output figure height.')
@click.option('-o', '--output-file', 'output_file',
              metavar='<file>', default='gene_structure.pdf', show_default=True,
              help='Output file. (support formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, and tiff)')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(input_file, file_format, utr_color, utr_hatch, cds_color, cds_hatch, exon_color, exon_hatch, edge_color, figure_width, figure_height, output_file):
    """Plot gene structure based on annotation file."""
    if file_format == 'gff':
        if exon_hatch or exon_color != 'salmon':
            click.echo('\033[33mThere are conflicting options, ignore "exon_color" and "exon_hatch".\033[0m')
    elif file_format == 'gtf':
        if utr_color != 'salmon' or utr_hatch or cds_color != 'skyblue' or cds_hatch:
            click.echo('\033[33mThere are conflicting options, ignore "utr_color", "utr_hatch", "cds_color" and '
                       '"cds_hatch".\033[0m')
    main(in_file=input_file,
         in_file_format=file_format,
         utr_color=utr_color,
         cds_color=cds_color,
         exon_color=exon_color,
         edge_color=edge_color,
         figure_width=figure_width,
         figure_height=figure_height,
         output_file=output_file,
         utr_hatch=utr_hatch,
         cds_hatch=cds_hatch,
         exon_hatch=exon_hatch)


if __name__ == '__main__':
    run()
