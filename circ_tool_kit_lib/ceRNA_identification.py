#!/usr/bin/env python
"""
File: ceRNA_identification.py
Description: Identify the expression correlation between ceRNA (mRNA-circRNA) with Pearson coefficient.
Date: 2023/4/11
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import mkdir
from os.path import exists
from warnings import filterwarnings
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import click
from Biolib.statistics import pd, read_in_gene_expression_as_dataframe
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def get_ceRNA_dict(ceRNA_file: str):
    ceRNA_dict = {}  # {miRNA: [circRNA_list, mRNA_list]}
    for line in open(ceRNA_file):
        split = line.strip().split('\t')
        miRNA, target, target_type = split[0], split[1], split[2]
        if miRNA not in ceRNA_dict:
            ceRNA_dict[miRNA] = [[], []]
        if target_type == 'circRNA':
            ceRNA_dict[miRNA][0].append(target)
        else:
            ceRNA_dict[miRNA][1].append(target)
    return ceRNA_dict


def plot_ceRNA_exp(miRNA_exp: pd.Series,
                   mRNA_exp: pd.Series,
                   circRNA_exp: pd.Series,
                   r, p, out_path):
    miRNA_plot, = plt.plot(range(1, len(miRNA_exp) + 1),
                           (miRNA_exp - miRNA_exp.min()) / (miRNA_exp.max() - miRNA_exp.min()),
                           marker='.', label='miRNA')
    mRNA_plot, = plt.plot(range(1, len(miRNA_exp) + 1),
                          (mRNA_exp - mRNA_exp.min()) / (mRNA_exp.max() - mRNA_exp.min()),
                          marker='.', label='mRNA')
    circRNA_plot, = plt.plot(range(1, len(miRNA_exp) + 1),
                             (circRNA_exp - circRNA_exp.min()) / (circRNA_exp.max() - circRNA_exp.min()),
                             marker='.', label='circRNA')
    plt.xticks(range(1, len(miRNA_exp) + 1), miRNA_exp.index.tolist())
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0], color='#3e3a39')
    plt.ylabel('Normalized Expression', color='#3e3a39')
    plt.tick_params('both', color='w')
    plt.title(f'{mRNA_exp.name}_{miRNA_exp.name}_{circRNA_exp.name}\nr={r}, p={p}', color='#3e3a39')
    plt.legend()
    plt.savefig(f"{out_path}ceRNA_exp_plot/"
                f"{mRNA_exp.name}_{miRNA_exp.name}_{str(circRNA_exp.name).replace(':', '_').replace('|', '_')}.pdf",
                bbox_inches='tight')
    plt.clf()
    return miRNA_plot, mRNA_plot, circRNA_plot


def main(mRNA_exp_file: str,
         circRNA_exp_file: str,
         ceRNA_file: str,
         miRNA_exp_file: str = None,
         out_path: str = './',
         ceRNA_r_cutoff: float = 0.9,
         ceRNA_p_cutoff: float = 0.05,
         miRNA_target_r_cutoff: float = -0.9,
         miRNA_target_P_cutoff: float = 0.05):
    if not exists(f'{out_path}ceRNA_exp_plot'):
        mkdir(f'{out_path}ceRNA_exp_plot')
    plt.style.use('ggplot')
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    filterwarnings("ignore")
    miRNA_exp_matrix = read_in_gene_expression_as_dataframe(miRNA_exp_file)
    mRNA_exp_matrix = read_in_gene_expression_as_dataframe(mRNA_exp_file)
    circRNA_exp_matrix = read_in_gene_expression_as_dataframe(circRNA_exp_file)
    content = '## miRNA\tcircRNA\tmRNA\tceRNA_r\tceRNA_p\tmiRNA_target_r\tmiRNA_target_p\n'
    ceRNA_dict = get_ceRNA_dict(ceRNA_file)
    for k, v in ceRNA_dict.items():
        miRNA_id = k
        circRNAs, mRNAs = v[0], v[1]
        miRNA_exp = miRNA_exp_matrix.loc[miRNA_id]
        for circRNA in circRNAs:
            circRNA_exp = circRNA_exp_matrix.loc[circRNA]
            for mRNA in mRNAs:
                mRNA_exp = mRNA_exp_matrix.loc[mRNA]
                r1, p1 = pearsonr(circRNA_exp, mRNA_exp)
                r2, p2 = pearsonr(miRNA_exp, mRNA_exp)
                r3, p3 = pearsonr(miRNA_exp, circRNA_exp)
                if r1 >= ceRNA_r_cutoff and p1 <= ceRNA_p_cutoff:
                    if r2 <= miRNA_target_r_cutoff and p2 <= miRNA_target_P_cutoff:
                        plot_ceRNA_exp(miRNA_exp, mRNA_exp, circRNA_exp, r1, p1, out_path)
                        content += f'{miRNA_id}\t{circRNA}\t{mRNA}\t{r1}\t{p1}\t{r2}\t{p2}\n'
                    elif r3 <= miRNA_target_r_cutoff and p3 <= miRNA_target_P_cutoff:
                        plot_ceRNA_exp(miRNA_exp, mRNA_exp, circRNA_exp, r1, p1, out_path)
                        content += f'{miRNA_id}\t{circRNA}\t{mRNA}\t{r1}\t{p1}\t{r3}\t{p3}\n'
    with open(f'{out_path}PCC_identified_ceRNA.txt', 'w') as o:
        o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-m', '--mRNA_exp_file', 'mrna', help='Input mRNA expression file. (Supported formats: txt, xls, xlsx and csv)')
@click.option('-c', '--circ_exp_file', 'circ', help='Input circRNA expression file. (Supported formats: txt, xls, xlsx and csv)')
@click.option('-C', '--ceRNA', 'cerna', help='Input ceRNA file as followed content:'
                                             '\nmiRNA156\\tChr01:100|1000\\tcircRNA\\n\nmiRNA156\\tPotri.001G001000.1\\tmRNA\\n')
@click.option('-M', '--miRNA_exp_file', 'mirna', help='[optional] Input miRNA expression file. '
                                                      '(Supported formats: txt, xls, xlsx and csv)')
@click.option('-o', '--output_path', 'out_path', default='./', help='[optional] Output path. {default=./}')
@click.option('-r1', type=float, default=0.9,
              help='[optional] Specify the cutoff of Pearson coefficient between ceRNA. {default=0.9}')
@click.option('-p1', type=float, default=0.05,
              help='[optional] Specify the cutoff of p value between ceRNA. {default=0.05}')
@click.option('-r2', type=float, default=-0.3,
              help='[optional] Specify the cutoff of Pearson coefficient between miRNA and target (mRNA/circRNA). '
                   '{default=-0.3}')
@click.option('-p2', type=float, default=1,
              help='[optional] Specify the cutoff of p value between miRNA and target (mRNA/circRNA). {default=1}')
def run(mrna, circ, cerna, mirna, out_path, r1, p1, r2, p2):
    """Identify the expression correlation between ceRNA (mRNA-circRNA) with Pearson coefficient."""
    main(mrna, circ, cerna, mirna, out_path, r1, p1, r2, p2)


if __name__ == '__main__':
    run()
