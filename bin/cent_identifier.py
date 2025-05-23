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
from pandas import read_table, DataFrame, Series
import matplotlib.pyplot as plt
import click
from pybioinformatic import Displayer

displayer = Displayer(__file__.split('/')[-1], version='0.2.1')


def apply_func(grouped_df: DataFrame,
               chr_len_dict: dict,
               min_cent_len: int,
               max_cent_len: int,
               out_path: str):
    makedirs(out_path, exist_ok=True)
    data_list = []
    chr_name = grouped_df['chr'].unique()[0]
    chr_len = chr_len_dict[chr_name]
    centromere_id = f'Cent_{chr_name}'
    for i in range(1, 11):
        for j in [k / 10 for k in range(1, 10)]:
            fold_enrichment = i + j
            filtered_df = grouped_df[grouped_df['fold_enrichment'] >= fold_enrichment]
            centromere_start = filtered_df['start'].min()
            centromere_end = filtered_df['end'].max()
            centromere_len = centromere_end - centromere_start + 1
            q_p = max([chr_len - centromere_end, centromere_start -1]) / min([chr_len - centromere_end, centromere_start -1])
            if 1 < q_p <= 1.7:
                karyotype = 'M'
            elif 1.7 < q_p <= 3:
                karyotype = 'SM'
            elif 3 < q_p <= 7:
                karyotype = 'ST'
            else:
                karyotype = 'T'
            if min_cent_len <= centromere_len <= max_cent_len and \
                    filtered_df.iloc[1, 1] - filtered_df.iloc[0, 2] <= min_cent_len * 0.1 and \
                    filtered_df.iloc[-1, 1] - filtered_df.iloc[-2, 2] <= min_cent_len * 0.1:
                data = [centromere_id, chr_name, chr_len, centromere_start, centromere_end, centromere_len, fold_enrichment, q_p, karyotype]
                data_list.append(data)
    columns = ['Centromere_id', 'Chromosome', 'Chromosome_length', 'Centromere_start', 'Centromere_end', 'Centromere_length', 'fold_enrichment', 'q:p', 'Karyotype']
    raw = DataFrame(data=data_list, columns=columns)
    Centromere_length_median = min_cent_len + (max_cent_len - min_cent_len) / 2
    raw['Diff'] = abs(raw['Centromere_length'] - Centromere_length_median)
    if raw.empty:
        best = Series(data=[centromere_id, chr_name, chr_len, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'], index=columns)
    else:
        best = raw.iloc[raw['Diff'].idxmin(), :-1]

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

    return best


def main(chr_len_file: TextIOWrapper,
        peaks_file: TextIOWrapper,
        min_cent_len,
        max_cent_len,
        out_path: str):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6
    chr_len_dict = {
        line.strip().split('\t')[0]: int(line.strip().split('\t')[1])
        for line in chr_len_file
    }

    df = read_table(peaks_file, index_col=-1, comment='#')
    ret = df.groupby(by='chr').apply(func=partial(
        apply_func,
        chr_len_dict=chr_len_dict,
        min_cent_len=min_cent_len,
        max_cent_len=max_cent_len,
        out_path=out_path
    ))
    ret.set_index(keys='Centromere_id', drop=True, inplace=True)
    ret.to_csv(f'{out_path}/Centromere_pos.xls', sep='\t')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-l', '--chr-len', 'chr_len_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help=r'Input chromosome length file. (Chr_name\tChr_length\tEtc)')
@click.option('-p', '--peaks-file', 'peaks_file',
              metavar='<file|stdin>', type=click.File('r'), required=True,
              help='Peaks file generated by macs2.')
@click.option('-min', '--min-cent-len', 'min_cent_len',
              metavar='<int>', type=int, default=400000, show_default=True,
              help='Minimal centromere length.')
@click.option('-max', '--max-cent-len', 'max_cent_len',
              metavar='<int>', type=int, default=650000, show_default=True,
              help='Maximum centromere length.')
@click.option('-o', '--output-path', 'out_path',
              metavar='<path>', default=getcwd(), show_default=True,
              help='Output path.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(chr_len_file, peaks_file, min_cent_len, max_cent_len, out_path):
    """Centromere identification based on ChIP-seq."""
    main(
        chr_len_file=chr_len_file,
        peaks_file=peaks_file,
        min_cent_len=min_cent_len,
        max_cent_len=max_cent_len,
        out_path=out_path
    )


if __name__ == '__main__':
    run()
