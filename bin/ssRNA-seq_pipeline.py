#!/usr/bin/env python
"""
File: ssRNA-seq_pipeline.py
Description: Strand specific RNA pair end sequencing analysis pipeline (including lncRNA and its target gene prediction).
Date: 2023/2/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from typing import Literal
from os import getcwd, makedirs, system
from os.path import abspath
from shutil import which
import click
from pybioinformatic import (
    parse_sample_info,
    MergeSamples,
    StrandSpecificRNASeqAnalyser,
    LncRNAPredictor,
    Displayer
)
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def main(sample_info: TextIOWrapper,
         genome: str,
         gff: str,
         feature_type: str,
         count_unit: str,
         module: Literal['ve', 'pl'],
         pfamscan_database: str,
         num_threads: int,
         num_processing: int,
         output_path: str):
    output_path = abspath(output_path)
    storer = MergeSamples(set(), set())
    d = parse_sample_info(sample_info=sample_info)

    # single sample analyse
    for sample_name, fq_list in d.items():
        makedirs(f'{output_path}/shell/normal', exist_ok=True)
        analyser = StrandSpecificRNASeqAnalyser(
            read1=fq_list[0],
            read2=fq_list[1],
            ref_gff=gff,
            ref_genome=genome,
            output_path=output_path,
            num_threads=num_threads,
            sample_name=sample_name
        )
        with open(f'{output_path}/shell/normal/{sample_name}.sh', 'w') as f:
            fastp = analyser.run_fastp(q_value=20, fastp_exe='fastp')
            hisat2 = analyser.run_hisat2(hisat2_exe='hisat2', samtools_exe='samtools', storer=storer)
            stringtie = analyser.run_stringtie(stringtie_exec='stringtie', storer=storer)
            cmds = [fastp, hisat2, stringtie]
            f.write('\n'.join(cmds))
    system(
        f'for i in `ls {output_path}/shell/normal`;'
        f'do echo sh {output_path}/shell/normal/$i;'
        f'done > {output_path}/shell/merge_normal.sh'
    )

    # merge results of each sample
    a = StrandSpecificRNASeqAnalyser(
        read1='',
        read2='',
        ref_genome=genome,
        output_path=output_path,
        num_threads=num_threads,
        ref_gff=gff
    )

    stringtie_merge = a.run_stringtie_merge(
        gtf_list=list(storer.stringtie_gtf),
        m=200,
        c=1.0,
        F=0.5,
        g=0,
        stringtie_exec='stringtie'
    )

    cuffcompare = a.run_cuffcompare(cuffcompare_exe='cuffcompare', gffread_exe='gffread')

    featureCounts = a.run_featureCounts(
        bam_list=list(storer.hisat2_bam),
        feature_type=feature_type,
        count_unit=count_unit,
        featureCounts_exec='featureCounts'
    )

    # lncRNA prediction
    makedirs(f'{output_path}/05.lncRNA_prediction', exist_ok=True)
    lp = LncRNAPredictor(
        nucl_fasta_file=f'{output_path}/03.assembly/novel_transcript.fa',
        output_path=f'{output_path}/05.lncRNA_prediction',
        pep_fasta_file=None,
        num_thread=num_threads,
        module=module
    )
    CNCI_cmd = lp.run_CNCI(CNCI_exec='CNCI.py')
    CPC2_cmd = lp.run_CPC2(CPC2_exec='CPC2.py')
    PLEK_cmd = lp.run_PLEK(PLEK_exec='PLEK')
    lncRNA_prediction_cmd = (
        f'{CNCI_cmd}\n{CPC2_cmd}\n{PLEK_cmd}\n'
        f'pfamscan_helper batch '
        f'-i {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input '
        f'-d {pfamscan_database} '
        f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out '
        f'-n {num_threads}')
    with open(f'{output_path}/shell/lncRNA_prediction.sh', 'w') as o:
        o.write(lncRNA_prediction_cmd)

    # lncRNA target gene prediction
    makedirs(f'{output_path}/06.lncRNA_target_prediction', exist_ok=True)
    target_prediction_script = f'''#!/usr/bin/env python
from pybioinformatic import LncRNATargetPredictor

ltp = LncRNATargetPredictor(
    lncRNA_gtf_file='{output_path}/06.lncRNA_target_prediction/lncRNA.gtf',
    mRNA_gtf_file='{output_path}/06.lncRNA_target_prediction/target.gtf',
    lncRNA_exp_file='{output_path}/06.lncRNA_target_prediction/lncRNA_exp.xls',
    mRNA_exp_file='{output_path}/06.lncRNA_target_prediction/target_exp.xls',
    output_path='{output_path}/06.lncRNA_target_prediction',
    lncRNA_min_exp=1,
    mRNA_min_exp=1,
    r=0.8,
    FDR=0.05,
    q_value=0.05,
    distance=100000,
    num_processing={num_threads},
)

if __name__ == '__main__':
    ltp.co_location()
    ltp.co_expression()
'''
    with open(f'{output_path}/shell/lncRNA_target_prediction.py', 'w') as o:
        o.write(target_prediction_script)

    # write all step commands
    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        exec_cmds = which('exec_cmds')
        featureCounts_helper = which('featureCounts_helper')
        seqkit = which('seqkit')
        ORF_finder = which('ORF_finder')
        file_split = which('file_split')
        lncRNA_results = lp.merge_results(
            CNCI_results=f'{output_path}/05.lncRNA_prediction/CNCI/CNCI.index',
            CPC2_results=f'{output_path}/05.lncRNA_prediction/CPC2/CPC2.txt',
            PLEK_results=f'{output_path}/05.lncRNA_prediction/PLEK/PLEK.xls',
            PfamScan_results=f'{output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out/all_results.xls',
            seqkit_exec=seqkit
        )
        cmd = (
            f'#!/bin/bash\n\n'
            f'{exec_cmds} -f {output_path}/shell/merge_normal.sh -n {num_processing}\n\n'
            f'{stringtie_merge}\n\n'
            f'{cuffcompare}\n\n'
            f'{featureCounts}\n\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/04.expression/featureCounts.xls '
            f'-o {output_path}/04.expression\n\n'
            f'{ORF_finder} -l 30 -P -log {output_path}/03.assembly/ORF.log '
            f'-o {output_path}/03.assembly {output_path}/03.assembly/novel_transcript.fa -n {num_threads}\n\n'
            f'{file_split} -n 200 '
            f'-i {output_path}/03.assembly/novel_transcript_pep.fa '
            f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input\n\n'
            f'{exec_cmds} -f {output_path}/shell/lncRNA_prediction.sh -n 4\n\n'
            f'{lncRNA_results}\n\n'
            f'grep ">" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -wf - {output_path}/03.assembly/All.gtf > {output_path}/06.lncRNA_target_prediction/lncRNA.gtf\n\n'
            f'grep ">" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -vwf - {output_path}/03.assembly/All.gtf > {output_path}/06.lncRNA_target_prediction/target.gtf\n\n'
            f'grep ">" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -wf - {output_path}/04.expression/FPKM.fc.xls | '
            f'cat <(head -n 1 {output_path}/04.expression/FPKM.fc.xls) - '
            f'> {output_path}/06.lncRNA_target_prediction/lncRNA_exp.xls\n\n'
            f'grep ">" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -vwf - {output_path}/04.expression/FPKM.fc.xls '
            f'> {output_path}/06.lncRNA_target_prediction/target_exp.xls\n\n'
            f'python {output_path}/shell/lncRNA_target_prediction.py'
        )
        o.write(cmd)
    system(f'chmod 755 {output_path}/shell/All_step.sh')
    click.echo(f'Commands created successfully, please run "bash {output_path}/shell/All_step.sh".', err=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-l', '--sample-info', 'sample_info',
              metavar='<file|stdin>', required=True, type=click.File('r'),
              help=r'Sample information file. (Sample_name\tFastq_path\tEtc)')
@click.option('-r', '--ref-genome', 'ref_genome',
              metavar='<file>', required=True,
              help='Reference genome fasta file.')
@click.option('-g', '--gff', 'gff',
              metavar='<file>', required=True,
              help='Reference genome gff annotation file.')
@click.option('-f', '--feature-type', 'feature_type',
              metavar='<str>', default='exon', show_default=True,
              help='Genome feature type, such as exon, CDS, etc. Multiple features are separated by commas.')
@click.option('-c', '--count-unit', 'count_unit',
              metavar='<str>', default='transcript_id', show_default=True,
              help='Raw reads count unit, such as transcript_id or gene_id.')
@click.option('-m', '--module', 'module',
              metavar='<str>', type=click.Choice(['ve', 'pl']), default='pl', show_default=True,
              help='')
@click.option('-d', '--pfamscan-database', 'pfamscan_database',
              metavar='<path>', required=True,
              help='PfamScan database path.')
@click.option('-t', '--num-threads', 'num_threads',
              metavar='<int>', type=int, default=10, show_default=True,
              help='The number of threads.')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=4, show_default=True,
              help='The number of processing.')
@click.option('-o', '--out-path', 'out_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(sample_info,
        ref_genome,
        gff,
        feature_type,
        count_unit,
        module,
        pfamscan_database,
        num_threads,
        num_processing,
        out_path):
    """Strand specific RNA pair end sequencing analysis pipeline (including lncRNA and its target gene prediction)."""
    main(
        sample_info=sample_info,
        genome=ref_genome,
        gff=gff,
        feature_type=feature_type,
        count_unit=count_unit,
        module=module,
        pfamscan_database=pfamscan_database,
        num_threads=num_threads,
        num_processing=num_processing,
        output_path=out_path
    )


if __name__ == '__main__':
    run()
