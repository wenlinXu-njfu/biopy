#!/usr/bin/env python
"""
File: multi_vcf2gt.py
Description: Merge genotypes for all samples from vcf files.
CreateDate: 2023/10/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Tuple
from io import TextIOWrapper, StringIO
from os.path import exists
from re import sub
from natsort import natsort_key
from pandas import read_table, concat
from tqdm import tqdm
import click
from pybioinformatic import VCF, Timer, Displayer
displayer = Displayer(__file__.split('/')[-1], version='0.2.0')


@Timer('Converting the vcf file to gt format.')
def main(vcf_files: Tuple[TextIOWrapper],
         output_file: TextIOWrapper,
         read_depth: int = 5,
         depth_files_dir: str = None,
         depth_file_suffix: str = 'depth'):
    if len(vcf_files) == 1 and depth_files_dir is None:
        with VCF(vcf_files[0]) as vcf:
            for line in vcf.to_genotype():
                click.echo(line, output_file)
    else:
        # merge genotype
        dfs = []
        for vcf_file in vcf_files:
            with VCF(vcf_file) as vcf:
                gt = '\n'.join([line for line in vcf.to_genotype()])
                df = read_table(StringIO(gt), index_col=[0, 1, 2, 3])
                dfs.append(df)
        merge = concat(dfs, axis=1)
        merge.reset_index(inplace=True)
        merge.sort_values(by=[merge.columns[1], merge.columns[2]],
                          key=natsort_key,
                          inplace=True)
        columns = merge.columns.tolist()
        head = columns[:4]
        samples = columns[4:]
        samples.sort(key=natsort_key)
        head.extend(samples)
        merge = merge.loc[:, head]

        # fill NA
        if depth_files_dir:
            miss_depth = []
            with tqdm(total=len(samples)) as pbar:
                for sample in samples:
                    depth_file1 = f'{depth_files_dir}/{sample}.{depth_file_suffix}'
                    depth_file2 = f'{depth_files_dir}/{sample}/{sample}.{depth_file_suffix}'
                    if exists(depth_file1):
                        depth = [
                            '_'.join(line.strip().split('\t')[:2])
                            for line in open(depth_file1)
                            if int(line.strip().split('\t')[2]) >= read_depth
                        ]
                        # SNP
                        merge.loc[merge[sample].isnull(), sample] = \
                            merge.loc[
                                (merge[sample].isnull()) &
                                (merge.loc[merge[sample].isnull(), 'ID'].isin(depth)) &
                                (merge['Ref'].str.len() == 1), 'Ref'
                            ] * 2
                        # non-SNP
                        merge.loc[
                            (merge[sample].isnull()) &
                            (merge.loc[merge[sample].isnull(), 'ID'].isin(depth)) &
                            (merge['Ref'].str.len() > 1), sample
                        ] = '-'
                    elif exists(depth_file2):
                        depth = [
                            '_'.join(line.strip().split('\t')[:2])
                            for line in open(depth_file2)
                            if int(line.strip().split('\t')[2]) >= read_depth
                        ]
                        # SNP
                        merge.loc[merge[sample].isnull(), sample] = \
                            merge.loc[
                                (merge[sample].isnull()) &
                                (merge.loc[merge[sample].isnull(), 'ID'].isin(depth)) &
                                (merge['Ref'].str.len() == 1), 'Ref'
                            ] * 2
                        # non-SNP
                        merge.loc[
                            (merge[sample].isnull()) &
                            (merge.loc[merge[sample].isnull(), 'ID'].isin(depth)) &
                            (merge['Ref'].str.len() > 1), sample
                        ] = '-'
                    else:
                        miss_depth.append(sample)
                    pbar.update()
            if miss_depth:
                stderr = '\n'.join([f'\033[33mWaring: Missing reads depth file of {sample}.' for sample in miss_depth])
                click.echo(stderr, err=True)

        # output results
        merge = merge.to_string(index=False, na_rep='NA').strip()
        merge = sub(r'\n +', '\n', merge)
        merge = sub(r' +', '\t', merge)
        click.echo(merge, output_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('vcf_files', nargs=-1, metavar='<vcf files|stdin>', type=click.File('r'), required=True)
@click.option('-d', '--depth-files-dir', 'depth_files_dir',
              metavar='<dir>', help='Reads depth file for each sample.')
@click.option('-D', '--read-depth', 'read_depth',
              metavar='<int>', default=5, show_default=True,
              help='Read depth threshold value.')
@click.option('-s', '--depth-file-suffix', 'depth_file_suffix',
              metavar='<str>', default='depth', show_default=True,
              help='Suffix of depth file for each sample.')
@click.option('-o', '--output_file', 'output_file',
              metavar='<gt file|stdout>', type=click.File('w'),
              help='Output file, stdout by default.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(vcf_files, depth_files_dir, depth_file_suffix, read_depth, output_file):
    """Merge genotypes for all samples from vcf files."""
    main(
        vcf_files=vcf_files,
        output_file=output_file,
        read_depth=read_depth,
        depth_files_dir=depth_files_dir,
        depth_file_suffix=depth_file_suffix
    )


if __name__ == '__main__':
    run()
