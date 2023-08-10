#!/usr/bin/env python
"""
File: plot_circos.py
Description: Plot circos graphic
Date: 2022/10/12
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import matplotlib.pyplot as plt
from plot_lib.circos.draw_gene_density import draw_chr, draw_gene_density, draw_outer_bar, draw_bezier_curve
from Biolib.show_info import Displayer


def main(chr_len_file, gene_density_file, stat_file, link_file, out_file):
    # Set matplotlib params.
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'

    # Draw chromosome bar.
    raw_chr_len_dict, chr_theta_dict, chr_width_dict, ax, bar_container = draw_chr(chr_len_file)

    # Draw gene density.
    ax, line2D = draw_gene_density(gene_density_file, chr_theta_dict, chr_width_dict, ax, color='#1f77b4')

    # Draw outer bar.
    ax, artist_list = draw_outer_bar(chr_theta_dict, chr_width_dict, ax, stat_file)

    # Draw bezier curve.
    used_color = []
    for line in open(link_file):
        split = line.strip().split('\t')
        chr1, start1, end1 = split[0], int(split[1]), int(split[2])
        chr2, start2, end2 = split[3], int(split[4]), int(split[5])
        color, description = split[-2], split[-1]
        angle1 = chr_theta_dict[chr1] - chr_width_dict[chr1] / 2
        angle2 = chr_theta_dict[chr2] - chr_width_dict[chr2] / 2
        x1 = (end1 + start1) / 2 / raw_chr_len_dict[chr1] * chr_width_dict[chr1] + angle1
        x2 = (end2 + start2) / 2 / raw_chr_len_dict[chr2] * chr_width_dict[chr2] + angle2
        if color not in used_color:
            used_color.append(color)
            _, line_2D = draw_bezier_curve(ax, x1, 7.6, x2, 7.6, 0.6, color, description)
            artist_list.append(line_2D)
        else:
            draw_bezier_curve(ax, x1, 7.6, x2, 7.6, 0.6, color)

    # Set legend.
    artist_list.insert(0, bar_container)
    artist_list.append(line2D)
    plt.legend(artist_list, [i.get_label() for i in artist_list], fontsize=6, loc=(0.9, 0.9), ncol=1)

    # Save graphic.
    plt.savefig(out_file, bbox_inches='tight', dpi=300)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-c', '--chr_len', 'chr_len_file', help='Input chromosome length file.\n(Chr_num\\tLength\\n)')
@click.option('-d', '--gene_density', 'gene_density_file',
              help='Input gene density file.\n(Chr_num\\tStart\\tEnd\\tCount\\n)')
@click.option('-s', '--feature_stat', 'feature_stat', help='Input statistical file.\n(Chr_num\\tType\\tCount\\tColor\\n)')
@click.option('-l', '--link_file', 'link_file',
              help='Input associated site file.\n(Chr_num\\tStart\\tEnd\\tChr_num\\tStart\\tEnd\\tColor\\tDescription\\n)')
@click.option('-o', '--output_file', 'outfile', default='circos.pdf', show_default=True, help='Output file.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=Displayer(__file__.split('/')[-1]).version_info)
def run(chr_len_file, gene_density_file, feature_stat, link_file, outfile):
    """Draw the circos graph."""
    main(chr_len_file, gene_density_file, feature_stat, link_file, outfile)


if __name__ == '__main__':
    run()
