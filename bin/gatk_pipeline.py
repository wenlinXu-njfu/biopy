#!/usr/bin/env python
"""
File: gatk_pipeline.py
Description: Variation analysis pipeline of GATK.
CreateDate: 2023/8/13
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import makedirs, system, getcwd
from os.path import abspath
from shutil import which
import click
from pybioinformatic import parse_sample_info, build_ref_index, GatkSNPCalling, Displayer
displayer = Displayer(__file__.split('/')[-1], version='1.0.1')


def check_dependency():
    software_list = ['fastp', 'bwa', 'samtools', 'gatk']
    click.echo('\033[36mDependency check.\033[0m', err=True)
    for software in software_list:
        path = which(software)
        if path:
            click.echo(f'{software}: {path}', err=True)
        else:
            click.echo(f'{software}: command not found', err=True)
    if not all([which(i) for i in software_list]):
        exit()


def main(genome_fasta_file: str,
         build_index: bool,
         sample_info: str,
         num_threads: int,
         num_processing: int,
         output_path: str):
    """Variation analysis pipeline of GATK."""
    check_dependency()
    output_path = abspath(output_path)
    sample_info_dict = parse_sample_info(sample_info)
    makedirs(f'{output_path}/shell/normal', exist_ok=True)

    # Single sample SNA calling
    for sample_name, fq_path in sample_info_dict.items():
        fq1, fq2 = fq_path
        with open(f'{output_path}/shell/normal/{sample_name}.sh', 'w') as o:
            sc = GatkSNPCalling(
                read1=fq1,
                read2=fq2,
                ref_genome=genome_fasta_file,
                output_path=output_path,
                num_threads=num_threads,
                sample_name=sample_name
            )
            pipeline = sc.pipeline()
            o.write(pipeline)
    system(
        f'for i in `ls {output_path}/shell/normal`; '
        f'do echo "sh {output_path}/shell/normal/$i"; '
        f'done > {output_path}/shell/run_normal.sh'
    )

    # Merge vcf
    gatk = which('gatk')
    CombineGVCFs = (
        f'{gatk} CombineGVCFs '
        f'-R {genome_fasta_file} '
        f'$(for i in `ls {output_path}/03.variant/*/*.gvcf`; do echo "-V $i" ;done) '
        f'-O {output_path}/03.variant/cohort.gvcf.gz'
    )
    GenotypeGVCFs = (
        f'{gatk} GenotypeGVCFs '
        f'-R {genome_fasta_file} '
        f'-V {output_path}/03.variant/cohort.gvcf.gz '
        f'-O {output_path}/03.variant/cohort.vcf.gz'
    )

    # Separate SNP and INDEL
    select_snp = (
        f'{gatk} SelectVariants '
        f'-R {genome_fasta_file} '
        f'--select-type-to-include SNP '
        f'-V {output_path}/03.variant/cohort.vcf.gz '
        f'-O {output_path}/03.variant/cohort.snp.vcf.gz'
    )
    select_indel = (
        f'{gatk} SelectVariants '
        f'-R {genome_fasta_file} '
        f'--select-type-to-include INDEL '
        f'-V {output_path}/03.variant/cohort.vcf.gz '
        f'-O {output_path}/03.variant/cohort.indel.vcf.gz'
    )

    # Vcf2GT
    file_format_conversion = which('file_format_conversion')
    vcf2gt = (
        f'{file_format_conversion} vcf2gt '
        f'-o {output_path}/03.variant/All.GT.xls '
        f'{output_path}/03.variant/cohort.vcf.gz'
    )

    # Stat genotype
    gt_kit = which('gt_kit')
    makedirs(f'{output_path}/04.stats', exist_ok=True)
    stat_gt = (
        f'{gt_kit} stat -gt {output_path}/03.variant/All.GT.xls -n {num_processing * num_threads} -o {output_path}/04.stats'
    )

    # Write all step commands
    exec_cmds = which('exec_cmds')
    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        cmds = [
            f'{exec_cmds} -f {output_path}/shell/run_normal.sh -n {num_processing}',
            CombineGVCFs,
            GenotypeGVCFs,
            select_snp,
            select_indel,
            vcf2gt,
            stat_gt
        ]
        if build_index:
            with open(f'{output_path}/shell/build_index.sh', 'w') as o2:
                build_index_cmd = build_ref_index(
                    fasta_file=genome_fasta_file,
                    bwa_exe='bwa',
                    gatk_exe='gatk',
                    samtools_exe='samtools'
                )
                o2.write(build_index_cmd)
            cmds.insert(0, f'{exec_cmds} -f {output_path}/shell/build_index.sh -n 3')
        o.write('\n\n'.join(cmds))
    system(f'chmod 755 {output_path}/shell/All_step.sh')
    click.echo(
        message='\033[32mCommands created successfully, please run "bash {output_path}/shell/All_step.sh".\033[0m',
        err=True
    )


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-r', '--ref_genome_file', 'genome_fasta_file',
              metavar='<fasta file>', required=True,
              help='Reference sequence file.')
@click.option('--build-index', 'build_index',
              is_flag=True, flag_value=True,
              help='Build genome index.')
@click.option('-l', '--sample-info', 'sample_info',
              metavar='<file>', required=True,
              help=r'Input sample information file. (SampleName\tRead1Path\tRead2Path)')
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
def run(genome_fasta_file: str,
        build_index: bool,
        sample_info: str,
        num_threads: int,
        num_processing: int,
        output_path: str):
    """Variation analysis pipeline of GATK."""
    main(
        genome_fasta_file=genome_fasta_file,
        build_index=build_index,
        sample_info=sample_info,
        num_threads=num_threads,
        num_processing=num_processing,
        output_path=output_path
    )


if __name__ == '__main__':
    run()
