#!/usr/bin/env python
"""
File: qPCR.py
Description: Calculate relative expression based on qPCR results.
Date: 2022/5/7
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from warnings import filterwarnings
from typing import Tuple
import click
import numpy as np
import pandas as pd
from scipy.stats import f_oneway
import matplotlib.pyplot as plt
from Biolib.show_info import Displayer


def main(in_file: str,
         ref_gene_name: str,
         control_sample_name: str,
         out_file: str,
         figsize: Tuple[float, float]):
    # step1: Read in qPCR results file.
    df = pd.read_excel(in_file, 'Results', header=39, usecols=['Sample Name', 'Target Name', 'CT'])
    df = df[df['Sample Name'].notnull()]
    all_sample_name = list(set(list(df['Sample Name'])))
    all_sample_name.sort()

    # step2: Isolate the reference gene.
    ref_gene_list = ref_gene_name.split(',')
    all_gene_name = list(set([i for i in df['Target Name']]))
    all_gene_name.sort()
    all_target_gene_name = [i for i in all_gene_name if i not in ref_gene_list]

    # step3: Calculate geometric mean of CT value of reference genes.
    d = {'Sample Name': [], 'Ct Mean': []}
    for sample_name in all_sample_name:
        ct_mean_product = 1
        for ref_gene in ref_gene_list:
            ref_gene_df = df[(df['Target Name'] == ref_gene) &
                             (df['Sample Name'] == sample_name) &
                             (df['CT'] != 'Undetermined')]
            if ref_gene_df.empty:
                err_msg = f'Error: Reference gene "{ref_gene}" lacks sample "{sample_name}".'
                click.echo(f'\033[31m{err_msg}\033[0m', err=True)
                exit()
            ct_mean = ref_gene_df['CT'].mean(0)
            ct_mean_product *= ct_mean
        ct_mean = ct_mean_product ** (1 / len(ref_gene_list))
        d['Sample Name'].append(sample_name)
        d['Ct Mean'].append(ct_mean)
    ref_gene_df = pd.DataFrame(d)

    # step4: Calculate 2^(-ΔΔCT) of target genes and significance between control sample and threat sample.
    target_gene_df_list = []
    for target_gene in all_gene_name:
        control_sample_ct_mean = None
        df_list = []
        for sample_name in all_sample_name:
            target_gene_df = df[(df['Target Name'] == target_gene) & (df['Sample Name'] == sample_name)]
            if target_gene_df.empty:
                if sample_name == control_sample_name:
                    err_msg = f'Error: Gene "{target_gene}" lacks control sample "{sample_name}".'
                    click.echo(f'\033[31m{err_msg}\033[0m', err=True)
                else:
                    warning_msg = f'Warning: Gene "{target_gene}" lacks sample "{sample_name}".'
                    click.echo(f'\033[33m{warning_msg}\033[0m', err=True)
            else:
                target_gene_df.loc[target_gene_df['CT'] == 'Undetermined', 'CT'] = np.inf
                ref_gene_ct_mean = float(ref_gene_df.loc[ref_gene_df['Sample Name'] == sample_name, 'Ct Mean'])
                target_gene_df['ΔCt'] = target_gene_df['CT'] - ref_gene_ct_mean
                target_gene_df['ΔCt Mean'] = target_gene_df['ΔCt'].mean()
                if sample_name == control_sample_name:
                    control_sample_ct_mean = float(list(set(target_gene_df['ΔCt Mean']))[0])
                df_list.append(target_gene_df)
        if df_list and control_sample_ct_mean is not None:
            target_gene_df = pd.concat(df_list)
            target_gene_df['ΔΔCt'] = target_gene_df['ΔCt'] - control_sample_ct_mean
            target_gene_df['2^(-ΔΔCt)'] = 2 ** (-target_gene_df['ΔΔCt'])
            target_gene_df['F'] = np.NAN
            target_gene_df['PR(>F)'] = np.NAN
            for sample_name in all_sample_name:
                CK = target_gene_df.loc[target_gene_df['Sample Name'] == control_sample_name, '2^(-ΔΔCt)']
                if sample_name == control_sample_name:
                    args = [CK, CK]
                    f, p = f_oneway(*args)
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, 'F'] = f
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, 'PR(>F)'] = p
                else:
                    threat = target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, '2^(-ΔΔCt)']
                    args = [CK, threat]
                    f, p = f_oneway(*args)
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, 'F'] = f
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, 'PR(>F)'] = p
                target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, '2^(-ΔΔCt) Mean'] = \
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, '2^(-ΔΔCt)'].mean()
                target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, 'Std'] = \
                    target_gene_df.loc[target_gene_df['Sample Name'] == sample_name, '2^(-ΔΔCt)'].std()
            target_gene_df_list.append(target_gene_df)
    target_gene_df = pd.concat(target_gene_df_list)

    # step5: Save relative expression results.
    target_gene_df.to_csv(f'{out_file}.csv', float_format="%.3f", index=False, encoding='gbk')

    # step6: Draw bar graph and save figure.
    target_gene_df = target_gene_df.loc[:,
                     ['Sample Name', 'Target Name', 'PR(>F)', '2^(-ΔΔCt) Mean', 'Std']].drop_duplicates()
    if len(all_target_gene_name) % 3 == 0:
        nrows = int(len(all_target_gene_name) / 3)
    else:
        nrows = int(len(all_target_gene_name) / 3) + 1
    fig, axes = plt.subplots(nrows, 3, figsize=figsize)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    i = j = 0
    for gene_name in all_target_gene_name:
        df = target_gene_df[target_gene_df['Target Name'] == gene_name]
        if nrows == 1:
            ax = axes[j]
        else:
            ax = axes[i, j]
        ax.bar('Sample Name', '2^(-ΔΔCt) Mean', yerr='Std', data=df, color='lightgrey', edgecolor='k',
               error_kw=dict(elinewidth=0.5, capsize=10))
        y = list(df['2^(-ΔΔCt) Mean'] + df['Std'] + 0.01)
        p = list(df['PR(>F)'])
        df2 = pd.DataFrame(dict(y=y, p=p))
        df2.loc[(df2['p'] > 0.05) | (df2['p'].isna()), 'text'] = ''
        df2.loc[(df2['p'] > 0.01) & (df2['p'] <= 0.05), 'text'] = '*'
        df2.loc[(df2['p'] > 0.001) & (df2['p'] <= 0.01), 'text'] = '**'
        df2.loc[df2['p'] <= 0.001, 'text'] = '***'
        k = 0
        while k < len(df.index.tolist()):
            temp_df = df2.iloc[[k]]
            text = ''.join(temp_df['text'].values)
            ax.text(k, float(temp_df['y']), text, ha='center', va='center')
            k += 1
        ax.set_title(gene_name, fontstyle='italic')
        ax.tick_params('x', color='w', rotation=45)
        if j == 2:
            j = 0
            i += 1
        else:
            j += 1
    if j == 2:
        if nrows == 1:
            plt.delaxes(axes[j])
        else:
            plt.delaxes(axes[i, j])
    elif j == 1:
        if nrows == 1:
            plt.delaxes(axes[j])
            plt.delaxes(axes[j + 1])
        else:
            plt.delaxes(axes[i, j])
            plt.delaxes(axes[i, j + 1])
    plt.subplots_adjust(wspace=0.5, hspace=0.4)
    plt.savefig(f'{out_file}.pdf', bbox_inches='tight')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', 'input_file', required=True, help='Input qPCR result excel file.')
@click.option('-r', '--ref_name', 'ref_name', required=True,
              help='Specify reference gene name. Multiple reference genes are separated by commas. (eg. 18s,ef1)')
@click.option('-c', '--control_name', 'control_name', required=True, help='Specify name of control sample.')
@click.option('-f', '--figure_size', 'figure_size', default='10x10', show_default=True, help='Figure size.')
@click.option('-o', '--output_file_prefix', 'out_prefix', help='Prefix of output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(input_file, ref_name, control_name, figure_size, out_prefix):
    """Calculate relative expression based on qPCR results."""
    filterwarnings("ignore")
    figure_size = tuple([float(i) for i in figure_size.split('x')])
    main(input_file, ref_name, control_name, out_prefix, figure_size)


if __name__ == '__main__':
    run()
