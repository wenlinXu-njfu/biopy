#!/usr/bin/env python
"""
File: ssRNA-seq_pipeline.py
Description: ssRNA-seq pipeline.
Date: 2023/2/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from os import getcwd, makedirs, system
from shutil import which
import click
from pybioinformatic import parse_sample_info, MergeSamples, SpecificStrandRNASeqAnalyser
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(sample_info: TextIOWrapper,
         genome: str,
         gff: str,
         feature_type: str,
         count_unit: str,
         num_threads: int,
         num_processing: int,
         output_path: str):
    storer = MergeSamples(set(), set())
    d = parse_sample_info(sample_info=sample_info)

    # single sample analyse
    for sample_name, fq_list in d.items():
        makedirs(f'{output_path}/shell/normal', exist_ok=True)
        analyser = SpecificStrandRNASeqAnalyser(
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
    a = SpecificStrandRNASeqAnalyser(
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

    with open(f'{output_path}/shell/run_all_normal.sh', 'w') as o:
        exec_cmds = which('exec_cmds')
        featureCounts_helper = which('featureCounts_helper')
        cmd = (f'{exec_cmds} -f {output_path}/shell/merge_normal.sh -n {num_processing}\n'
               f'{stringtie_merge}\n'
               f'{cuffcompare}\n'
               f'{featureCounts}\n'
               f'{featureCounts_helper} normalization -i {output_path}/04.expression/featureCounts.xls -o {output_path}/04.expression')
        o.write(cmd)


@click.command(context_settings=CONTEXT_SETTINGS)
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
@click.option('-t', '--num-threads', 'num_threads',
              metavar='<int>', type=int, default=10, show_default=True,
              help='The number of threads.')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=4, show_default=True,
              help='The number of processing.')
@click.option('-o', '--out-path', 'out_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output path, if not exist, automatically created.')
def run(sample_info, ref_genome, gff, feature_type, count_unit, num_threads, num_processing, out_path):
    """Strand specific RNA pair end sequencing analysis pipeline."""
    main(
        sample_info=sample_info,
        genome=ref_genome,
        gff=gff,
        feature_type=feature_type,
        count_unit=count_unit,
        num_threads=num_threads,
        num_processing=num_processing,
        output_path=out_path
    )


if __name__ == '__main__':
    run()
