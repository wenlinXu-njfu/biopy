#!/usr/bin/env python
"""
File: gatk_pipeline.py
Description: Variation analysis pipeline of GATK.
CreateDate: 2023/8/13
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from os import makedirs, system, getcwd
from os.path import abspath
from shutil import which
from datetime import datetime
import click
from pybioinformatic import Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def check_dependency():
    software_list = ['fastp', 'bwa', 'samtools', 'gatk']
    click.echo('Dependency check.\n', err=True)
    for software in software_list:
        path = which(software)
        if path:
            click.echo(f'{software}: {path}', err=True)
        else:
            click.echo(f'{software}: command not found', err=True)
    if not all([which(i) for i in software_list]):
        exit()


def build_genome_index(genome_fasta_file: str):
    bwa = which('bwa')
    bwa_index = f'{bwa} index {genome_fasta_file}'
    samtools = which('samtools')
    samtools_index = f'{samtools} faidx {genome_fasta_file}'
    gatk = which('gatk')
    gatk_index = f'{gatk} CreateSequenceDictionary -R {genome_fasta_file}'
    cmds = [bwa_index, samtools_index, gatk_index]
    return '\n'.join(cmds)


def fastp(sample_name: str, fq1: str, fq2: str, num_threads: int, output_path: str):
    fastp = which('fastp')
    cmd = (f'{fastp} -w {num_threads} -i {fq1} -I {fq2} '
           f'-o {output_path}/{sample_name}_clean_1.fq.gz '
           f'-O {output_path}/{sample_name}_clean_2.fq.gz '
           f'-h {output_path}/{sample_name}.fastp.html '
           f'-j {output_path}/{sample_name}.fastp.json '
           f'-n 10 -q 20 -u 40 2> {output_path}/{sample_name}.fastp.log')
    return cmd


def align(sample_name: str,
          genome_fasta_file: str,
          fq1: str,
          fq2: str,
          num_threads: int,
          output_path: str):
    cmds = []
    bwa = which('bwa')
    samtools = which('samtools')
    gatk = which('gatk')
    # Align command
    cmd1 = (fr'{bwa} mem -t {num_threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" {genome_fasta_file} {fq1} {fq2} | '
            fr'{samtools} view -b -S -T {genome_fasta_file} - > {output_path}/{sample_name}.bam')
    cmds.append(cmd1)
    # map=30 filter
    awk = r'''awk '{if($1~/@/){print}else{if( $7 == "=" &&  $5 >= 30 ){print $0}}}' '''
    cmd2 = (fr'{samtools} view -h {output_path}/{sample_name}.bam | {awk} | '
            fr'{samtools} view -b -S -T {genome_fasta_file} - > {output_path}/{sample_name}.map30.bf.bam')
    cmds.append(cmd2)
    # Sort bam file
    cmd3 = f'{samtools} sort -@ {num_threads} -o {output_path}/{sample_name}.map30.sort.bam {output_path}/{sample_name}.map30.bf.bam'
    cmds.append(cmd3)
    cmd4 = f'{samtools} index {output_path}/{sample_name}.map30.sort.bam'
    cmds.append(cmd4)
    cmd5 = f'rm -rf {output_path}/{sample_name}.map30.bf.bam'
    cmds.append(cmd5)
    # MarkDuplicates reads
    cmd6 = (f'{gatk} MarkDuplicates '
            f'-I {output_path}/{sample_name}.map30.sort.bam '
            f'-M {output_path}/{sample_name}.map30.sort.bam.metrics '
            f'-O {output_path}/{sample_name}.map30.sort.markdup.bam')
    cmds.append(cmd6)
    cmd7 = f'{samtools} index {output_path}/{sample_name}.map30.sort.markdup.bam'
    cmds.append(cmd7)
    # Stat reads depth
    cmd8 = f'{samtools} depth {output_path}/{sample_name}.map30.sort.bam > {output_path}/{sample_name}.map30.depth'
    cmds.append(cmd8)
    return cmds


def HaplotypeCaller(genome_fasta_file: str,
                     bam_file: str,
                     sample_name: str,
                     output_path: str):
    cmds = []
    gatk = which('gatk')
    cmd1 = (f'{gatk} HaplotypeCaller '
             f'-I {bam_file} '
             f'-R {genome_fasta_file} -ERC GVCF '
             f'-O {output_path}/{sample_name}.map30.gvcf')
    cmds.append(cmd1)
    cmd2 = (f'{gatk} GenotypeGVCFs '
            f'-R {genome_fasta_file} '
            f'-V {output_path}/{sample_name}.map30.gvcf '
            f'-O {output_path}/{sample_name}.map30.vcf')
    cmds.append(cmd2)
    cmd3 = (f'{gatk} VariantFiltration '
            f'--filter-name  "HARD_TO_VALIDATE" '
            f'--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" '
            f'-R {genome_fasta_file} '
            f'-V {output_path}/{sample_name}.map30.vcf '
            f'-O {output_path}/{sample_name}.map30.filt.vcf')
    cmds.append(cmd3)
    return cmds


def main(fq_path: str,
         genome_fasta_file: str,
         build_index: bool,
         sample_list: TextIOWrapper,
         num_threads: int,
         num_processing: int,
         output_path: str,
         read_depth: int = 5):
    """Variation analysis pipeline of GATK."""
    check_dependency()

    fq_path = abspath(fq_path)
    genome_fasta_file = abspath(genome_fasta_file)
    output_path = abspath(output_path)

    makedirs(f'{output_path}/shell/normal', exist_ok=True)
    for line in sample_list:
        fq_prefix = line.strip().split('\t')[0]
        sample_name = line.strip().split('\t')[1]
        makedirs(f'{output_path}/01.QC/{sample_name}', exist_ok=True)
        makedirs(f'{output_path}/02.mapping/{sample_name}', exist_ok=True)
        makedirs(f'{output_path}/03.variant/{sample_name}', exist_ok=True)
        fq1 = f'{fq_path}/{fq_prefix}_1.fq.gz'
        fq2 = f'{fq_path}/{fq_prefix}_2.fq.gz'
        with open(f'{output_path}/shell/normal/{sample_name}.sh', 'w') as o:
            cmds = [
                fastp(sample_name=sample_name,
                      fq1=fq1, fq2=fq2,
                      num_threads=num_threads,
                      output_path=f'{output_path}/01.QC/{sample_name}')
            ]
            cmds.extend(
                align(sample_name=sample_name,
                      genome_fasta_file=genome_fasta_file,
                      fq1=fq1, fq2=fq2,
                      num_threads=num_threads,
                      output_path=f'{output_path}/02.mapping/{sample_name}')
            )
            cmds.extend(
                HaplotypeCaller(genome_fasta_file=genome_fasta_file,
                                bam_file=f'{output_path}/02.mapping/{sample_name}/{sample_name}.map30.sort.markdup.bam',
                                sample_name=sample_name,
                                output_path=f'{output_path}/03.variant/{sample_name}')
            )
            o.write('\n'.join(cmds))

    system(f'for i in `ls {output_path}/shell/normal`; do echo "sh {output_path}/shell/normal/$i"; done > {output_path}/shell/run_normal.sh')

    # get genotype
    file_format_conversion = which('file_format_conversion')
    depth_dir = f'{output_path}/02.mapping'
    vcf_dir = f'{output_path}/03.variant'
    output_file = f'{output_path}/03.variant/All.GT.xls'
    vcf2gt = f'{file_format_conversion} vcf2gt -d {depth_dir} -D {read_depth} -s "map30.depth" -o {output_file} {vcf_dir}/*/*.map30.filt.vcf'

    # merge genotype
    gatk = which('gatk')
    CombineGVCFs = f'{gatk} CombineGVCFs -R {genome_fasta_file} $(for i in `ls {output_path}/03.variant/*/*.gvcf`; do echo "-V $i" ;done) -O {output_path}/03.variant/cohort.gvcf'
    GenotypeGVCFs = f'{gatk} GenotypeGVCFs -R {genome_fasta_file} -V {output_path}/03.variant/cohort.gvcf -O {output_path}/03.variant/cohort.vcf'

    # write all step commands
    exec_cmds = which('exec_cmds')
    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        cmds = [
            f'{exec_cmds} -f {output_path}/shell/run_normal.sh -n {num_processing}',
            vcf2gt, CombineGVCFs, GenotypeGVCFs
        ]
        if build_index:
            with open(f'{output_path}/shell/build_index.sh', 'w') as o2:
                build_index_cmd = build_genome_index(genome_fasta_file)
                o2.write(build_index_cmd)
            cmds.insert(0, f'{exec_cmds} -f {output_path}/shell/build_index.sh -n 3')
        o.write('\n\n'.join(cmds))
    system(f'chmod 755 {output_path}/shell/All_step.sh')
    click.echo(f'Commands created successfully, please run "bash {output_path}/shell/All_step.sh".', err=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-fq', '--fq-path', 'fq_path',
              metavar='<dir>', required=True,
              help='Input fastq path.')
@click.option('-r', '--ref_genome_file', 'genome_fasta_file',
              metavar='<fasta file>', required=True,
              help='Reference sequence file.')
@click.option('--build-index', 'build_index',
              is_flag=True, flag_value=True,
              help='Build genome index.')
@click.option('-d', '--read-depth', 'read_depth',
              metavar='<int>', type=int, default=5, show_default=True,
              help='The number of reads supported variant.')
@click.option('-l', '--sample_list', 'sample_list',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input sample list file. (FastqPrefix\tSampleName)')
@click.option('-t', '--num-threads', 'num_threads',
              metavar='<int>', type=int, default=10, show_default=True,
              help='The number of threads for each sample.')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=5, show_default=True,
              help='Number of processing. It means how many samples are analyzed in parallel at a time.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<str>', default=getcwd(), show_default=True,
              help='Output path, if not exist, automatically created.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fq_path: str,
        genome_fasta_file: str,
        build_index: bool,
        read_depth: int,
        sample_list: TextIOWrapper,
        num_threads: int,
        num_processing: int,
        output_path: str):
    """Variation analysis pipeline of GATK."""
    main(
        fq_path=fq_path,
        genome_fasta_file=genome_fasta_file,
        build_index=build_index,
        sample_list=sample_list,
        num_threads=num_threads,
        num_processing=num_processing,
        output_path=output_path,
        read_depth=read_depth
    )


if __name__ == '__main__':
    run()
