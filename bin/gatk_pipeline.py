#!/usr/bin/env python
"""
File: gatk_pipeline.py
Description: Variation analysis pipeline of GATK.
CreateDate: 2023/8/13
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper, StringIO
from os import makedirs, system, listdir
from collections import defaultdict
from pandas import read_table
import click
from pybioinformatic import TaskManager, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.1.0')


def build_genome_index(genome_fasta_file: str):
    from pybioinformatic import TaskManager
    bwa_index = f'bwa index {genome_fasta_file}'
    samtools_index = f'samtools faidx {genome_fasta_file}'
    gatk_index = f'gatk CreateSequenceDictionary -R {genome_fasta_file}'
    cmds = [bwa_index, samtools_index, gatk_index]
    tkm = TaskManager(commands=cmds, num_processing=3)
    tkm.parallel_run_cmd()


def fastp(sample_name: str, fq1: str, fq2: str, output_path: str):
    cmd = (f'fastp -i {fq1} -I {fq2} '
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
          output_path: str):
    cmds = []
    # Align command
    cmd1 = fr'bwa mem -t 5 -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" {genome_fasta_file} {fq1} {fq2} \
    | samtools view -b -S -T {genome_fasta_file} - > {output_path}/{sample_name}.bam'
    cmds.append(cmd1)
    # map=30 filter
    awk = r'''awk '{if($1~/@/){print}else{if( $7 =="=" &&  $5>=30 ){print $0}}}' '''
    cmd2 = (fr'samtools view -h {output_path}/{sample_name}.bam | {awk} | '
            fr'samtools view -b -S -T {genome_fasta_file} - > {output_path}/{sample_name}.map30.bf.bam')
    cmds.append(cmd2)
    # Sort bam file
    cmd3 = f'samtools sort -@ 5 -o {output_path}/{sample_name}.map30.sort.bam {output_path}/{sample_name}.map30.bf.bam'
    cmds.append(cmd3)
    cmd4 = f'samtools index {output_path}/{sample_name}.map30.sort.bam'
    cmds.append(cmd4)
    cmd5 = f'rm -rf {output_path}/{sample_name}.map30.bf.bam'
    cmds.append(cmd5)
    # MarkDuplicates reads
    cmd6 = (f'gatk MarkDuplicates '
            f'-I {output_path}/{sample_name}.map30.sort.bam '
            f'-M {output_path}/{sample_name}.map30.sort.bam.metrics '
            f'-O {output_path}/{sample_name}.map30.sort.markdup.bam')
    cmds.append(cmd6)
    cmd7 = f'samtools index {output_path}/{sample_name}.map30.sort.markdup.bam'
    cmds.append(cmd7)
    # Stat reads depth
    cmd8 = f'samtools depth {output_path}/{sample_name}.map30.sort.bam > {output_path}/{sample_name}.map30.depth'
    cmds.append(cmd8)
    return cmds


def HaplotypeCaller(genome_fasta_file: str,
                     bam_file: str,
                     sample_name: str,
                     output_path: str):
    cmds = []
    cmd1 = (f'gatk HaplotypeCaller '
             f'-I {bam_file} '
             f'-R {genome_fasta_file} -ERC GVCF '
             f'-O {output_path}/{sample_name}.map30.gvcf')
    cmds.append(cmd1)
    cmd2 = (f'gatk GenotypeGVCFs -R {genome_fasta_file} '
             f'-V {output_path}/{sample_name}.map30.gvcf '
             f'-O {output_path}/{sample_name}.map30.vcf')
    cmds.append(cmd2)
    cmd3 = (f'gatk VariantFiltration '
             f'--filter-name  "HARD_TO_VALIDATE" '
             f'--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" '
             f'-R {genome_fasta_file} '
             f'-V {output_path}/{sample_name}.map30.vcf '
             f'-O {output_path}/{sample_name}.map30.filt.vcf')
    cmds.append(cmd3)
    return cmds


def get_gt(depth_dir: str,
           vcf_dir: str,
           depth: int = 5,
           tkm: TaskManager = None):
    depth_dict = defaultdict(list)  # {sample: [positions]}
    for sample in listdir(depth_dir):
        with open(f'{depth_dir}/{sample}/{sample}.map30.depth') as f:
            for line in f:
                split = line.strip().split('\t')
                if int(split[2]) >= depth:
                    depth_dict[sample].append('_'.join(split[:2]))

    if tkm is None:
        tkm = TaskManager(num_processing=1)
    cmd = f'file_format_conversion vcf2gt {vcf_dir}/*/*.map30.vcf'
    stdout = tkm.echo_and_exec_cmd(cmd)
    gt_df = read_table(StringIO(stdout), index_col=0)
    for sample in gt_df.columns[3:]:
        mask = (gt_df.index.isin(depth_dict[sample])) & (gt_df[sample].isna())
        gt_df.loc[mask, sample] = gt_df.loc[mask, 'Ref'] * 2
    return gt_df


def main(fq_path: str,
         genome_fasta_file: str,
         build_index: bool,
         sample_list: TextIOWrapper,
         num_processing: int,
         output_path: str):
    """Variation analysis pipeline of GATK."""
    if build_index:
        build_genome_index(genome_fasta_file)
    makedirs(f'{output_path}/shell', exist_ok=True)
    for line in sample_list:
        fq_prefix = line.strip().split('\t')[0]
        sample_name = line.strip().split('\t')[1]
        makedirs(f'{output_path}/QC/{sample_name}', exist_ok=True)
        makedirs(f'{output_path}/align/{sample_name}', exist_ok=True)
        makedirs(f'{output_path}/variant/{sample_name}', exist_ok=True)
        fq1 = f'{fq_path}/{fq_prefix}_1.fq.gz'
        fq2 = f'{fq_path}/{fq_prefix}_2.fq.gz'
        with open(f'{output_path}/shell/{sample_name}.sh', 'w') as o:
            cmds = [
                fastp(sample_name=sample_name,
                      fq1=fq1, fq2=fq2,
                      output_path=f'{output_path}/QC/{sample_name}')
            ]
            cmds.extend(
                align(sample_name=sample_name,
                      genome_fasta_file=genome_fasta_file,
                      fq1=fq1, fq2=fq2,
                      output_path=f'{output_path}/align/{sample_name}')
            )
            cmds.extend(
                HaplotypeCaller(genome_fasta_file=genome_fasta_file,
                                bam_file=f'{output_path}/align/{sample_name}/{sample_name}.map30.sort.markdup.bam',
                                sample_name=sample_name,
                                output_path=f'{output_path}/variant/{sample_name}')
            )
            o.write('\n'.join(cmds))
    system(f'for i in `ls {output_path}/shell`; do echo "sh {output_path}/shell/$i"; done > {output_path}/All.sh')
    system(f'exec_cmds -f {output_path}/All.sh -n {num_processing}')
    gt = get_gt(depth_dir=f'{output_path}/align/', vcf_dir=f'{output_path}/variant/')
    gt.to_csv(f'{output_path}/All.GT.xls', sep='\t', na_rep='NA')


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
@click.option('-samples', '--sample_list', 'sample_list',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input sample list file. (FastqPrefix\tSampleName)')
@click.option('-p', '--num-processing', 'num_processing',
              metavar='<int>', type=int, default=5, show_default=True,
              help='Number of processing.')
@click.option('-o', '--output-path', 'output_path',
              metavar='<str>', required=True,
              help='Output path.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(fq_path: str,
        genome_fasta_file: str,
        build_index: bool,
        sample_list: TextIOWrapper,
        num_processing: int,
        output_path: str):
    """Variation analysis pipeline of GATK."""
    from datetime import datetime
    start_time = datetime.now().replace(microsecond=0)
    main(fq_path, genome_fasta_file, build_index, sample_list, num_processing, output_path)
    end_time = datetime.now().replace(microsecond=0)
    click.echo(f'[{datetime.now().replace(microsecond=0)}] Total time spent {end_time - start_time}.', err=True)


if __name__ == '__main__':
    run()
