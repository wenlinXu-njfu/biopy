#!/usr/bin/env python
"""
File: ssRNA-seq_pipeline.py
Description: Strand specific (dUTP constructed library) pair end RNA-seq analysis pipeline (including lncRNA and its target gene prediction).
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
    LncRNAClassification,
    Displayer
)
displayer = Displayer(__file__.split('/')[-1], version='0.4.0')


def check_dependency():
    software_list = ['fastp', 'hisat2', 'gffread', 'stringtie', 'cuffcompare', 'featureCounts',
                     'CNCI.py', 'CPC2.py', 'PLEK', 'pfam_scan.pl', 'bedtools']
    click.echo('Dependency check.\n', err=True)
    for software in software_list:
        path = which(software)
        if path:
            click.echo(f'{software}: {path}', err=True)
        else:
            click.echo(f'{software}: command not found', err=True)
    if not all([which(i) for i in software_list]):
        exit()


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
    check_dependency()
    output_path = abspath(output_path)
    storer = MergeSamples(set(), set())
    d = parse_sample_info(sample_info=sample_info)
    script_template_path = '/'.join(__file__.split('/')[:-2]) + '/script_template'

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
        g=500,
        stringtie_exec='stringtie'
    )

    cuffcompare = a.run_cuffcompare(cuffcompare_exe='cuffcompare', gffread_exe='gffread')

    anno_file = 'gtf' if feature_type == 'exon' else 'gff'
    all_featureCounts = a.run_featureCounts(
        anno_file=anno_file,
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
    pfamscan_helper = which('pfamscan_helper')
    lncRNA_prediction_cmd = (
        f'{CNCI_cmd}\n{CPC2_cmd}\n{PLEK_cmd}\n'
        f'{pfamscan_helper} batch '
        f'-i {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input '
        f'-d {pfamscan_database} '
        f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out '
        f'-n {num_threads}')
    with open(f'{output_path}/shell/lncRNA_prediction.sh', 'w') as o:
        o.write(lncRNA_prediction_cmd)

    # lncRNA target gene prediction
    makedirs(f'{output_path}/06.lncRNA_target_prediction', exist_ok=True)
    with open(f'{script_template_path}/lncRNA_target_prediction') as f, \
        open(f'{output_path}/shell/lncRNA_target_prediction.py', 'w') as o:
        args = [output_path for _ in range(5)]
        args.append(num_threads * num_processing)
        args = tuple(args)
        target_prediction_script = f.read() % args
        o.write(target_prediction_script)

    # lncRNA classification
    lc = LncRNAClassification(
        mRNA_gff_file=gff,
        lncRNA_gtf_file=f'{output_path}/06.lncRNA_target_prediction/lncRNA.gtf',
        out_dir=f'{output_path}/07.lncRNA_classification'
    )
    classify_script = lc.classification()
    with open(f'{output_path}/shell/lncRNA_classification.sh', 'w') as o:
        o.write(classify_script)

    # write all step commands
    python = which('python')
    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        exec_cmds = which('exec_cmds')
        featureCounts = which('featureCounts')
        featureCounts_helper = which('featureCounts_helper')
        seqkit = which('seqkit')
        ORF_finder = which('ORF_finder')
        file_split = which('file_split')
        plot = which('plot')
        CNCI_results = r'''awk '{if($3 =="noncoding"){print $1}}' %s/05.lncRNA_prediction/CNCI/CNCI.index | sort -uV''' % output_path
        CPC2_results = r'''awk '{if($9 == "noncoding"){print $1}}' %s/05.lncRNA_prediction/CPC2/CPC2.txt | sort -uV''' % output_path
        PfamScan_results = r'''cut -f1 %s/05.lncRNA_prediction/PfamScan/pfamscan_out/all_results.xls | grep -vFw -f - <(grep '>' %s/03.assembly/novel_transcript.fa | sed 's/>//') | cut -d' ' -f 1 | sort -uV''' % (
        output_path, output_path)
        PLEK_results = r'''awk '{if($1 == "Non-coding"){print $3}}' %s/05.lncRNA_prediction/PLEK/PLEK.xls | sed 's/>//' | sort -uV''' % output_path
        plot_venn = f'''{plot} venn -g CNCI,CPC2,PfamScan,PLEK <({CNCI_results}) <({CPC2_results}) <({PfamScan_results}) <({PLEK_results}) -o {output_path}/05.lncRNA_prediction/venn.pdf'''
        lncRNA_results = lp.merge_results(
            CNCI_results=f'{output_path}/05.lncRNA_prediction/CNCI/CNCI.index',
            CPC2_results=f'{output_path}/05.lncRNA_prediction/CPC2/CPC2.txt',
            PLEK_results=f'{output_path}/05.lncRNA_prediction/PLEK/PLEK.xls',
            PfamScan_results=f'{output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out/all_results.xls',
            seqkit_exec=seqkit
        )
        create_lncRNA_gtf = (
            f'grep "^>" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -Fwf - {output_path}/03.assembly/All.gtf > {output_path}/06.lncRNA_target_prediction/lncRNA.gtf'
        )
        makedirs(f'{output_path}/06.lncRNA_target_prediction/lncRNA_exp', exist_ok=True)
        makedirs(f'{output_path}/06.lncRNA_target_prediction/target_exp', exist_ok=True)
        lncRNA_exp = (
            f'{featureCounts} -t exon -g transcript_id -fMO -p --countReadPairs -s 2 -T {num_threads} '
            f'-a {output_path}/06.lncRNA_target_prediction/lncRNA.gtf '
            f'-o {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.xls '
            f'{output_path}/02.mapping/*/*.ht2.sort.bam '
            f'2> {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.log\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.xls '
            f'-o {output_path}/06.lncRNA_target_prediction/lncRNA_exp'
        )
        create_target_gtf = (
            f'grep "^>" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -vFwf - {output_path}/03.assembly/All.gtf > {output_path}/06.lncRNA_target_prediction/target.gtf'
        )
        target_exp = (
            f'{featureCounts} -t exon -g transcript_id -fMO -p --countReadPairs -s 2 -T {num_threads} '
            f'-a {output_path}/06.lncRNA_target_prediction/target.gtf '
            f'-o {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.xls '
            f'{output_path}/02.mapping/*/*.ht2.sort.bam '
            f'2> {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.log\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.xls '
            f'-o {output_path}/06.lncRNA_target_prediction/target_exp'
        )
        cmd = (
            f'#!/bin/bash\n\n'
            f'{exec_cmds} -f {output_path}/shell/merge_normal.sh -n {num_processing}\n\n'
            f'{stringtie_merge}\n\n'
            f'{cuffcompare}\n\n'
            f'{all_featureCounts}\n\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/04.expression/featureCounts.xls '
            f'-o {output_path}/04.expression\n\n'
            f'{ORF_finder} -l 30 -F -log {output_path}/03.assembly/ORF.log '
            f'-o {output_path}/03.assembly {output_path}/03.assembly/novel_transcript.fa -n {num_threads}\n\n'
            f'{file_split} -n 200 '
            f'-i {output_path}/03.assembly/novel_transcript_pep.fa '
            f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input\n\n'
            f'{exec_cmds} -f {output_path}/shell/lncRNA_prediction.sh -n 4\n\n'
            f'{plot_venn}\n\n'
            f'{lncRNA_results}\n\n'
            f'{create_lncRNA_gtf}\n\n'
            f'{create_target_gtf}\n\n'
            f'{lncRNA_exp}\n\n'
            f'{target_exp}\n\n'
            f'{python} {output_path}/shell/lncRNA_target_prediction.py\n\n'
            f'bash {output_path}/shell/lncRNA_classification.sh\n'
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
              help='CNCI -m option.')
@click.option('-d', '--pfamscan-database', 'pfamscan_database',
              metavar='<path>', required=True,
              help='PfamScan database path.')
@click.option('-t', '--num-threads', 'num_threads',
              metavar='<int>', type=int, default=10, show_default=True,
              help='The number of threads for each sample.')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=4, show_default=True,
              help='The number of processing. It means how many samples are analyzed in parallel at a time.')
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
    """Strand specific (dUTP constructed library) pair end RNA-seq analysis pipeline (including lncRNA and its target gene prediction)."""
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
