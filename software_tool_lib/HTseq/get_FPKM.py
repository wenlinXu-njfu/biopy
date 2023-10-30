#!/usr/bin/env python
"""
File: get_FPKM.py
Description: Standardize gene expression with FPKM based on HTSeq results.
CreateDate: 2023/6/15
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from pandas import read_table
import click
from pybioinformatic import get_FPKM, Gtf, Gff, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def main(header_file: TextIOWrapper,
         htseq_file: TextIOWrapper,
         anno_file: Union[str, TextIOWrapper],
         anno_format: click.Choice(['gff', 'gtf']),
         min_exp: float,
         output_file: Union[str, TextIOWrapper]):
    parse_anno_file = {'gff': Gff, 'gtf': Gtf}
    with parse_anno_file[anno_format](anno_file) as anno:
        length_dict = {}
        for line in anno.parse():
            if line[2] == 'exon':
                length = int(line[4]) - int(line[3]) + 1
                if isinstance(anno, Gff):
                    transcript_id = line[-1]['Parent']
                else:
                    transcript_id = line[-1]['transcript_id']
                if transcript_id in length_dict:
                    length_dict[transcript_id] += length
                else:
                    length_dict[transcript_id] = length
        columns = [line.strip() for line in header_file]
        df = read_table(htseq_file, index_col=0, names=columns).iloc[:-5]
        df.insert(0, 'length', 0)
        for transcript_id in df.index.tolist():
            df.loc[transcript_id, 'length'] = length_dict[transcript_id]
        FPKM = get_FPKM(df, min_exp)
        FPKM.to_csv(output_file, sep='\t')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--header_info', 'header_info_file',
              metavar='<header file>', type=click.File('r'), required=True,
              help='Input header file. (eg. Gene_id\\nSample1\\nSample2\\netc\\n)')
@click.option('-I', '--htseq_result', 'htseq_result_file',
              metavar='<htseq file>', type=click.File('r'), required=True,
              help='Input Htseq results file.')
@click.option('-a', '--anno_file', 'anno_file',
              metavar='<anno file>', type=click.File('r'), required=True,
              help='Input genome annotation GFF or GTF file, must contain exon information.')
@click.option('-f', '--format', 'anno_format',
              metavar='<gff|gtf>', type=click.Choice(['gff', 'gtf']), default='gff', show_default=True,
              help='Specify the format of anno file.')
@click.option('-m', '--min_exp', 'min_exp',
              metavar='<float>', type=float, default=0, show_default=True,
              help='Gene minimum expression threshold in all samples.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<file>', type=click.File('w'), default='htseq_FPKM.xls', show_default=True,
              help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(header_info_file, htseq_result_file, anno_file, anno_format, min_exp, output_file):
    """Standardize gene expression with FPKM based on HTSeq results."""
    main(header_info_file, htseq_result_file, anno_file, anno_format, min_exp, output_file)


if __name__ == '__main__':
    run()
