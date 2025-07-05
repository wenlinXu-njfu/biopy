#!/usr/bin/env python
"""
File: cent_identifier.py
Description: Centromere identification based on ChIP-seq.
CreateDate: 2024/11/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from functools import partial
from os import makedirs, getcwd
from tqdm import tqdm
import numpy as np
from pandas import read_table, DataFrame, Series, concat
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import click
from pybioinformatic import Fasta, TaskManager, Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.3.0')


def cluster_and_filter_intervals(
    df: DataFrame,
    max_gap: int = 1000000,
    z_threshold: float = 2.5
) -> DataFrame:
    """
    Use hierarchical clustering to group the intervals and remove the discrete intervals.
    """
    df['midpoint'] = (df['start'] + df['end']) / 2

    results = []
    for chrom in df['chr'].unique():
        chrom_df = df[df['chr'] == chrom].copy().reset_index(drop=True)

        if len(chrom_df) <= 1:
            chrom_df['group'] = 0
            results.append(chrom_df)
            continue

        midpoints = chrom_df['midpoint'].values.reshape(-1, 1)

        dist_matrix = pdist(midpoints, metric='cityblock')

        Z = linkage(dist_matrix, method='single')
        clusters = fcluster(Z, t=max_gap, criterion='distance')

        chrom_df['group'] = clusters

        filtered_dfs = []
        for group_id, group_df in chrom_df.groupby('group'):
            if len(group_df) == 1:
                filtered_dfs.append(group_df)
                continue

            points = group_df['midpoint'].values
            median = np.median(points)
            mad = np.median(np.abs(points - median))

            if mad == 0:
                mad = np.std(points) / 1.4826 if np.std(points) > 0 else 1
                z_scores = np.abs(points - median) / (np.std(points) + 1e-9)
            else:
                z_scores = np.abs(points - median) / (1.4826 * mad)

            is_not_outlier = z_scores <= z_threshold
            filtered_df = group_df[is_not_outlier].copy()

            if len(filtered_df) >= 1:
                filtered_dfs.append(filtered_df)

        if filtered_dfs:
            results.append(concat(filtered_dfs))

    result_df = concat(results).reset_index(drop=True)
    result_df = result_df.drop(columns=['midpoint'])

    return result_df


def apply_func(
    grouped_df: DataFrame,
    chr_seq_dict: dict,
    min_cent_len: int,
    max_cent_len: int,
    homogenize: bool,
    out_path: str,
    pbar
) -> Series:
    data_list = []
    chr_name = grouped_df['chr'].unique()[0]
    chr_seq = chr_seq_dict[chr_name]
    chr_len = len(chr_seq)
    centromere_id = f'Cent_{chr_name}'
    for _, sub_df in grouped_df.groupby('group'):
        if sub_df.iloc[-1, 2] - sub_df.iloc[0, 1] > min_cent_len:
            for i in range(1, 11):
                for j in [k / 10 for k in range(1, 10)]:
                    fold_enrichment = i + j
                    filtered_df = sub_df[sub_df['fold_enrichment'] >= fold_enrichment]
                    centromere_start = filtered_df['start'].min()
                    centromere_end = filtered_df['end'].max()
                    centromere_len = centromere_end - centromere_start + 1
                    q_p = (
                            max([chr_len - centromere_end, centromere_start - 1]) /
                            min([chr_len - centromere_end, centromere_start - 1])
                    )
                    if 1 < q_p <= 1.7:
                        karyotype = 'M'
                    elif 1.7 < q_p <= 3:
                        karyotype = 'SM'
                    elif 3 < q_p <= 7:
                        karyotype = 'ST'
                    else:
                        karyotype = 'T'
                    if min_cent_len <= centromere_len <= max_cent_len:
                        num_interval = len(filtered_df)
                        data = [centromere_id, chr_name, chr_len, centromere_start, centromere_end,
                                centromere_len, fold_enrichment, q_p, karyotype, num_interval]
                        data_list.append(data)

    columns = ['Centromere_id', 'Chromosome', 'Chromosome_length', 'Centromere_start', 'Centromere_end',
               'Centromere_length', 'fold_enrichment', 'q:p', 'Karyotype', 'num_interval']
    raw = DataFrame(data=data_list, columns=columns)
    Centromere_length_median = min_cent_len + (max_cent_len - min_cent_len) / 2
    raw['Diff'] = abs(raw['Centromere_length'] - Centromere_length_median)
    if raw.empty:
        best = Series(data=[centromere_id, chr_name, chr_len, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'], index=columns[:-1])
    else:
        best = raw.iloc[raw['Diff'].idxmin(), :-2] if homogenize else \
            raw.iloc[raw['num_interval'].idxmax(), :-2]

    # plot
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    for line in grouped_df.iterrows():
        width = line[1]['length']
        x = line[1]['start'] + width / 2
        height = line[1]['fold_enrichment']
        ax1.bar(x, height, ec='k')
        ax2.bar(x, height, ec='k')
    ax1.set_xticks(
        [i for i in range(10 ** 6, chr_len, 10 ** 6)],
        [i for i in range(1, int(chr_len / 1000000) + 1)],
        rotation=45
    )
    if not raw.empty:
        ax2.plot([0, chr_len], [best['fold_enrichment'], best['fold_enrichment']], 'r--', linewidth=0.8)
        ax2.vlines(
            x=[best['Centromere_start'], best['Centromere_end']],
            ymin=0,
            ymax=grouped_df['fold_enrichment'].max(),
            colors='orange',
            linestyles='dashed',
            linewidth=0.8
        )
        ax2.set_xticks(
            [i for i in range(10 ** 6, chr_len, 10 ** 5)],
            [str(i / 1000000) for i in range(10 ** 6, chr_len, 10 ** 5)],
            rotation=45
        )
        range_start = best['Centromere_start'] - 10 ** 6
        range_end = best['Centromere_end'] + 10 ** 6
        ax2.set_xlim([range_start, range_end])
    else:
        ax2.set_xticks(
            [i for i in range(10 ** 6, chr_len, 10 ** 6)],
            [i for i in range(1, int(chr_len / 1000000) + 1)],
            rotation=45
        )
    ax2.set_xlabel('length (Mb)')
    ax1.set_ylabel('fold_enrichment')
    ax2.set_ylabel('fold_enrichment')
    ax1.set_title(chr_name)
    plt.savefig(f'{out_path}/{chr_name}_peaks.pdf')
    pbar.update(1)

    return best


def main(
    chr_fasta_file: str,
    peaks_file: TextIOWrapper,
    min_cent_len: int,
    max_cent_len: int,
    z_threshold: float,
    homogenize: bool,
    out_path: str
):
    makedirs(out_path, exist_ok=True)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.size'] = 6
    with Fasta(chr_fasta_file) as fa:
        chr_seq_dict = fa.to_dict(parse_id=True)

    df = read_table(peaks_file, index_col=-1, comment='#')
    df = cluster_and_filter_intervals(df=df, max_gap=int((min_cent_len + max_cent_len) / 2), z_threshold=z_threshold)

    with tqdm(total=len(chr_seq_dict)) as pbar:
        ret = df.groupby(by='chr').apply(func=partial(
            apply_func,
            chr_seq_dict=chr_seq_dict,
            min_cent_len=min_cent_len,
            max_cent_len=max_cent_len,
            homogenize=homogenize,
            out_path=out_path,
            pbar=pbar
        ))
    ret.set_index(keys='Centromere_id', drop=True, inplace=True)
    ret.to_csv(f'{out_path}/centromere_pos.xls', sep='\t')

    tkm = TaskManager(num_processing=1)
    bed = r'''sed '1d' %s/centromere_pos.xls | awk -F'\t' '{print $2,$4-1,$5,$1}' OFS='\t' ''' % out_path
    cmd = f'bedtools getfasta -fi {chr_fasta_file} -bed <({bed}) -nameOnly > {out_path}/centromere_seq.fa'
    tkm.echo_and_exec_cmd(cmd, show_cmd=False)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--chr-seq', 'chr_seq_file',
              metavar='<file>', required=True,
              help=r'Input chromosome sequence fasta file.')
@click.option('-p', '--peaks-file', 'peaks_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Peaks file generated by macs2.')
@click.option('-min', '--min-cent-len', 'min_cent_len',
              metavar='<int>', type=int, default=400000, show_default=True,
              help='Minimal centromere length.')
@click.option('-max', '--max-cent-len', 'max_cent_len',
              metavar='<int>', type=int, default=650000, show_default=True,
              help='Maximum centromere length.')
@click.option('-z', '--z-threshold', 'z_threshold',
              metavar='<float>', type=float, default=10.0, show_default=True,
              help='Discrete interval Z-value threshold.')
@click.option('--homogenize', 'homogenize',
              is_flag=True, flag_value=True,
              help='Maintaining the consistency of the telomere length.')
@click.option('-o', '--output-path', 'out_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output path.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(chr_seq_file, peaks_file, min_cent_len, max_cent_len, z_threshold, homogenize, out_path):
    """Centromere identification based on ChIP-seq."""
    main(
        chr_fasta_file=chr_seq_file,
        peaks_file=peaks_file,
        min_cent_len=min_cent_len,
        max_cent_len=max_cent_len,
        z_threshold=z_threshold,
        homogenize=homogenize,
        out_path=out_path
    )


if __name__ == '__main__':
    run()
