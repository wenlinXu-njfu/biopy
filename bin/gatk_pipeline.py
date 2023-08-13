#!/usr/bin/env python
"""
File: gatk_pipeline.py
Description: Variation analysis pipeline of GATK.
Date: 2023/8/13
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import system
from typing import List
import click
from Biolib.timer import Timer
from Biolib.show_info import Displayer


@Timer('Build genome index.')
def build_genome_index(genome_fasta_file: str):
    command = f'samtools faidx {genome_fasta_file}'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    command = f'gatk CreateSequenceDictionary -R {genome_fasta_file}'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)


@Timer('Merge bam files.')
def merge_bam_files(sorted_bam_files: list, output_file_prefix: str):
    bam_files = ' '.join(sorted_bam_files)
    command = f'samtools merge {output_file_prefix}.merge.bam {bam_files}'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{output_file_prefix}.merge.bam'


@Timer('Mark duplicates.')
def MarkDuplicates(sorted_bam_file: str):
    markdup_bam_file = sorted_bam_file.replace('bam', 'markdup.bam')
    command = f'gatk MarkDuplicates ' \
              f'-I {sorted_bam_file} ' \
              f'-M {sorted_bam_file}.metrics ' \
              f'-O {markdup_bam_file}'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    command = f'samtools index {markdup_bam_file}'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return markdup_bam_file


@Timer('Variation detection.')
def HaplottypeCaller(genome_fasta_file: str,
                     bam_file: str,
                     sample_name: str):
    command = f'gatk HaplotypeCaller ' \
              f'--emit-ref-confidence GVCF ' \
              f'-R {genome_fasta_file} ' \
              f'-I {bam_file} ' \
              f'--sample-name {sample_name} ' \
              f'-O {sample_name}.HC.gvcf.gz'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    command = f'gatk IndexFeatureFile -F {sample_name}.HC.gvcf.gz'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{sample_name}.HC.gvcf.gz'


@Timer('SNP calling.')
def SNP_calling(genome_fasta_file: str,
                vcf_file: str,
                sample_name: str):
    command = f'gatk SelectVariants ' \
              f'-R {genome_fasta_file} ' \
              f'--select-type-to-include SNP ' \
              f'-V {vcf_file} ' \
              f'-O {sample_name}.HC.snp.vcf.gz'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{sample_name}.HC.snp.vcf.gz'


@Timer('SNP filtering.')
def SNP_filter(snp_vcf_file: str, filter_expression: str, sample_name: str):
    command = f'gatk SelectVariantsFilteration ' \
              f'-V {snp_vcf_file} ' \
              f'-O {sample_name}.HC.snp.filter.vcf.gz ' \
              f'--filter-expression {filter_expression} ' \
              f'--filter-name PASS'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{sample_name}.HC.snp.filter.vcf.gz'


@Timer('INDEL calling.')
def INDEL_calling(genome_fasta_file: str,
                  vcf_file: str,
                  sample_name: str):
    command = f'gatk SelectVariants ' \
              f'-R {genome_fasta_file} ' \
              f'--select-type-to-include INDEL ' \
              f'-V {vcf_file} ' \
              f'-O {sample_name}.HC.indel.vcf.gz'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{sample_name}.HC.indel.vcf.gz'


@Timer('INDEL filtering.')
def INDEL_filter(indel_vcf_file: str, filter_expression: str, sample_name: str):
    command = f'gatk SelectVariantsFilteration ' \
              f'-V {indel_vcf_file} ' \
              f'-O {sample_name}.HC.indel.filter.vcf.gz ' \
              f'--filter-expression {filter_expression} ' \
              f'--filter-name PASS'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)
    return f'{sample_name}.HC.indel.filter.vcf.gz'


def merge_variant(snp_vcf_file: str,
                  indel_vcf_file: str,
                  sample_name: str):
    command = f'gatk mergeVcfs -I {snp_vcf_file} -I {indel_vcf_file} -O {sample_name}.HC.filter.vcf.gz'
    click.echo(f'\033[36m{command}\033[0m', err=True)
    system(command)


def main(genome_fasta_file: str,
         sample_list: List[str],
         snp_filter_expression: str,
         indel_filter_expression: str,
         output_prefix: str,
         sorted_bam_files: List[str]):
    """Variation analysis pipeline of GATK."""
    build_genome_index(genome_fasta_file)
    merge_bam_file = merge_bam_files(sorted_bam_files, output_prefix)
    markdup_bam_file = MarkDuplicates(merge_bam_file)
    for sample_name in sample_list:
        gvcf_file = HaplottypeCaller(genome_fasta_file, markdup_bam_file, sample_name)
        snp_vcf_file = SNP_calling(genome_fasta_file, gvcf_file, sample_name)
        snp_vcf_file = SNP_filter(snp_vcf_file, snp_filter_expression, sample_name)
        indel_vcf_file = INDEL_calling(genome_fasta_file, gvcf_file, sample_name)
        indel_vcf_file = INDEL_filter(indel_vcf_file, indel_filter_expression, sample_name)
        merge_variant(snp_vcf_file, indel_vcf_file, sample_name)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-r', '--ref_genome_file', 'genome_fasta_file', type=click.File('r'), required=True,
              help='Reference sequence file.')
@click.option('-samples', '--sample_list', 'sample_list', type=click.File('r'), required=True,
              help='Sample name list. (One sample name per line, eg. Control\\nThreat\\n)')
@click.option('-snp', '--snp_filter', 'snp_filter_expression', help='SNP filter expression.')
@click.option('-indel', '--indel_filter', 'indel_filter_expression', help='INDEL filter expression.')
@click.option('-o', '--output_prefix', 'output_prefix', required=True, help='Output file prefix.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
@click.argument('sorted_bam_files', nargs=-1, type=click.File('r'), required=True)
def run(genome_fasta_file, sample_list, snp_filter_expression, indel_filter_expression, output_prefix, sorted_bam_files):
    """Variation analysis pipeline of GATK."""
    from datetime import datetime
    start_time = datetime.now().replace(microsecond=0)
    if sample_list.name == '<stdin>':
        sample_list = [line.strip() for line in click.open_file('-').readlines()]
    else:
        sample_list = [line.strip() for line in sample_list]
    genome_fasta_file = genome_fasta_file.name
    sorted_bam_files = [bam_file.name for bam_file in sorted_bam_files]
    main(genome_fasta_file, sample_list, snp_filter_expression, indel_filter_expression, output_prefix, sorted_bam_files)
    end_time = datetime.now().replace(microsecond=0)
    click.echo(f'[{datetime.now().replace(microsecond=0)}] Total time spent {end_time - start_time}.', err=True)


if __name__ == '__main__':
    run()
