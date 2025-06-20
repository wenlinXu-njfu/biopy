#!/usr/bin/env python
"""
File: centromere_identification.py
Description: Centromere identification based on ChIP-seq.
CreateDate: 2025/6/1
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from os import makedirs, system
from os.path import abspath
from shutil import which
from yaml import safe_load
from schema import Schema, And, Use, SchemaError
import click
from pybioinformatic import check_cmds, parse_sample_info, build_ref_index, Macs2PeakCalling, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


def check_config(yaml_file: TextIOWrapper):
    click.echo('\033[36mCheck params.\033[0m', err=True)
    config_schema = Schema(
        {
            "input": {
                "sample_info": str
            },
            "output": {
                "dir": str
            },
            "global_params": {
                "build_index": bool,
                "num_threads": And(Use(int), lambda x: 0 < x, error="num_threads must be positive integer."),
                "num_processing": And(Use(int), lambda x: 0 < x, error="num_processing must be positive integer.")
            },
            "centromere_identification_params": {
                "min_len": And(Use(int), lambda x: 0 < x, error="min_len must be positive integer."),
                "max_len": And(Use(int), lambda x: 0 < x, error="max_len must be positive integer."),
                "z_threshold": float,
                "homogenize": bool
            }
        }
    )

    config = safe_load(yaml_file)
    try:
        config_schema.validate(config)
        click.echo("\033[32mConfiguration verification passed.\n\033[0m", err=True)
    except SchemaError as e:
        click.echo(f"\033[31mConfig error: {e}\033[0m", err=True)
        exit()
    else:
        return config


def main(config_file: Union[TextIOWrapper, str]):
    # check dependency and config
    check_cmds(cmds_list=['fastp', 'bowtie2', 'samtools', 'gatk', 'macs2'])
    config = check_config(config_file)

    # input
    sample_info = config['input']['sample_info']
    sample_info_dict = parse_sample_info(sample_info=sample_info)

    # output
    output_path = abspath(config['output']['dir'])

    # global params
    build_index = config['global_params']['build_index']
    num_threads = config['global_params']['num_threads']
    num_processing = config['global_params']['num_processing']

    # centromere identification params
    min_len = config['centromere_identification_params']['min_len']
    max_len = config['centromere_identification_params']['max_len']
    z_threshold = config['centromere_identification_params']['z_threshold']
    homogenize = config['centromere_identification_params']['homogenize']

    # commands path
    exec_cmds = which('exec_cmds')
    get_seq_len = which('get_seq_len')
    macs2_helper = which('macs2_helper')

    makedirs(f'{output_path}/shell/samples', exist_ok=True)
    makedirs(f'{output_path}/04.centromere', exist_ok=True)
    for sample_name, fq_list in sample_info_dict.items():
        with open(f'{output_path}/shell/samples/{sample_name}.sh', 'w') as o:
            ChIP_1, ChIP_2, ref_genome, Input_1, Input_2 = fq_list[:5]

            cmd1 = Macs2PeakCalling(
                ChIP_read1=ChIP_1,
                ChIP_read2=ChIP_2,
                Input_read1=Input_1,
                Input_read2=Input_2,
                ref_genome=ref_genome,
                output_path=output_path,
                num_threads=num_threads,
                sample_name=sample_name
            ).pipeline(macs2_callpeak_options={'--broad': ''})

            cmd2 = (
                f'{macs2_helper} cent_identifier '
                f'-i {ref_genome} '
                f'-p {output_path}/03.peaks/{sample_name}/{sample_name}_peaks.xls '
                f'-min {min_len} '
                f'-max {max_len} '
                f'-z {z_threshold} '
                f'-o {output_path}/04.centromere/{sample_name}'
            )
            if homogenize:
                cmd2 += ' --homogenize'
            pipeline = [cmd1, cmd2]

            if build_index:
                pipeline.insert(
                    __index=0,
                    __object=build_index(
                        fasta_file=ref_genome,
                        bowtie2_build='bowtie2-build',
                        large=True
                    )
                )

            o.write('\n'.join(pipeline))

    system(
        f'for i in `ls {output_path}/shell/samples`;'
        f'do echo sh {output_path}/shell/samples/$i;'
        f'done > {output_path}/shell/All_samples.sh'
    )

    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        o.write(f'{exec_cmds} -f {output_path}/shell/All_samples.sh -n {num_processing}')

    system(f'chmod 755 {output_path}/shell/All_step.sh')
    click.echo(
        message=f'\033[32mCommands created successfully, please run "bash {output_path}/shell/All_step.sh".\033[0m',
        err=True
    )


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('config', nargs=1, metavar='<yaml file|stdin>', type=click.File('r'), required=True)
def run(config):
    """Centromere identification based on ChIP-seq."""
    main(config)


if __name__ == '__main__':
    run()
