#!/usr/bin/env python
"""
File: plot_figure.py
Description: Plot tools.
Date: 2022/2/24
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from plot_lib import plot_venn, plot_volcano, plot_heatmap, plot_cluster_heatmap, plot_gene_structure, plot_circos
from Biolib import Displayer
displayer = Displayer(__file__.split('/')[-1], version='1.1.0')


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def plot():
    """Plot tools."""
    pass


plot.add_command(plot_venn, name='venn')
plot.add_command(plot_volcano, name='volcano')
plot.add_command(plot_heatmap, name='heatmap')
plot.add_command(plot_cluster_heatmap, name='cluster_heatmap')
plot.add_command(plot_gene_structure, name='gene_structure')
plot.add_command(plot_circos, name='circos')


if __name__ == '__main__':
    plot()
