#!/usr/bin/env python
"""
File: get_TPM.py
Description: Standardize gene expression with TPM based on HTSeq results.
Date: 2023/6/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from Biolib.statistics import click, pd, get_TPM
from Biolib.gtf import Gtf
from Biolib.gff import Gff
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(header_file: str, htseq_file: str, anno_file: str, min_exp: float, out_prefix: str):
    parse_anno_file = {'gff': Gff, 'gff3': Gff, 'gtf': Gtf}
    file_obj = parse_anno_file[anno_file.split('.')[-1]](anno_file)
    length_dict = {}
    for line in file_obj.parse():
        if line[2] == 'exon':
            length = int(line[4]) - int(line[3]) + 1
            if isinstance(file_obj, Gff):
                transcript_id = line[-1]
            else:
                transcript_id = line[-2]
            if transcript_id in length_dict:
                length_dict[transcript_id] += length
            else:
                length_dict[transcript_id] = length
    columns = [line.strip() for line in open(header_file)]
    df = pd.read_table(htseq_file, index_col=0, names=columns).iloc[:-5]
    df.insert(0, 'length', 0)
    for transcript_id in df.index.tolist():
        df.loc[transcript_id, 'length'] = length_dict[transcript_id]
    TPM = get_TPM(df, min_exp)
    TPM.to_csv(f'./{out_prefix}.csv')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--header_infor', 'header_infor_file', help='Input header file. (eg. Gene_id\\nSample1\\nSample2\\netc\\n)')
@click.option('-I', '--htseq_result', 'htseq_result_file', help='Input Htseq results file.')
@click.option('-a', '--anno_file', 'anno_file', help='Input genome annotation GFF or GTF file, must contain exon information.')
@click.option('-m', '--min_exp', 'min_exp', type=float, default=0,
              help='[optional] Gene minimum expression threshold in all samples. {default: 0}')
@click.option('-o', '--output_prefix', 'output_prefix', default='htseq_TPM',
              help='[optional] Prefix of output file. {default: htseq_TPM}')
def run(header_infor_file, htseq_result_file, anno_file, min_exp, output_prefix):
    """Standardize gene expression with TPM based on HTSeq results."""
    main(header_infor_file, htseq_result_file, anno_file, min_exp, output_prefix)


if __name__ == '__main__':
    run()
