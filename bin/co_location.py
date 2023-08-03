#!/usr/bin/env python
"""
File: co_location.py
Description: Target genes were predicted according to the location relationship between lncRNA and genes.
Date: 2023/4/16
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gff import Gff


def judge_distance_location(lncRNA_start: int,
                            lncRNA_end: int,
                            lncRNA_strand: str,
                            gene_start: int,
                            gene_end: int,
                            gene_strand: str) -> tuple:
    """
    According to the position relationship between lncRNA and target gene, the type of lncRNA was determined.
    :param lncRNA_start: lncRNA start site on chromosome (type=int)
    :param lncRNA_end: lncRNA end site on chromosome (type=int)
    :param lncRNA_strand: lncRNA chain direction (type=str, value='+' or '-')
    :param gene_start: gene start site on chromosome (type=int)
    :param gene_end: gene end site on chromosome (type=int)
    :param gene_strand: gene chain direction (type=str, value='+' or '-')
    :return: (distance between lncRNA and target gene, location of lncRNA relative to target genes,
              whether lncRNA and target gene are in the same chain) (type=tuple(int,str,str))
    """
    if lncRNA_start < lncRNA_end < gene_start < gene_end:
        if lncRNA_strand == gene_strand == '+':
            return gene_start - lncRNA_end, 'upstream', 'sense'
        elif lncRNA_strand == gene_strand == '-':
            return gene_start - lncRNA_end, 'downstream', 'sense'
        elif lncRNA_strand != gene_strand and gene_strand == '+':
            return gene_start - lncRNA_end, 'upstream', 'antisense'
        else:
            return gene_start - lncRNA_end, 'downstream', 'antisense'
    elif lncRNA_start < gene_start < lncRNA_end < gene_end:
        if lncRNA_strand == gene_strand:
            return lncRNA_end - gene_start, 'overlap', 'sense'
        else:
            return lncRNA_end - gene_start, 'overlap', 'antisense'
    elif gene_start < lncRNA_start < lncRNA_end < gene_end:
        if lncRNA_strand == gene_strand:
            return lncRNA_end - lncRNA_start, 'overlap', 'sense'
        else:
            return lncRNA_end - lncRNA_start, 'overlap', 'antisense'
    elif gene_start < lncRNA_start < gene_end < lncRNA_end:
        if lncRNA_strand == gene_strand:
            return gene_end - lncRNA_start, 'overlap', 'sense'
        else:
            return gene_end - lncRNA_start, 'overlap', 'antisense'
    elif gene_start < gene_end < lncRNA_start < lncRNA_end:
        if lncRNA_strand == gene_strand == '+':
            return lncRNA_start - gene_end, 'downstream', 'sense'
        elif lncRNA_strand == gene_strand == '-':
            return lncRNA_start - gene_end, 'upstream', 'sense'
        elif lncRNA_strand != gene_strand and gene_strand == '+':
            return lncRNA_start - gene_end, 'downstream', 'antisense'
        else:
            return lncRNA_start - gene_end, 'upstream', 'antisense'
    else:
        if lncRNA_strand == gene_strand:
            return gene_end - gene_start, 'overlap', 'sense'
        else:
            return gene_end - gene_start, 'overlap', 'antisense'


def lncRNA_target_gene_prediction(gtf_file: str,
                                  gff_file: str,
                                  feature: str,
                                  target_range: int,
                                  out_file: str = None):
    """
    LncRNA target genes were predicted according to lncRNA annotation files (GTF) and gene annotation files (GFF).
    :param gtf_file: lncRNA annotation files(GTF)
    :param gff_file: gene annotation files(GFF)
    :param out_file: output file
    :param feature: feature type (type=str)
    :param target_range: the farthest distance between lncRNA and target gene (type=int, unit=1000)
    :return: None
    """
    gff_dict = Gff(gff_file).get_gff_dict(feature)
    content = ['Chr_num\tLncRNA_id\tTarget_id\tDistance\tLocation\tLncRNA_strand\n']
    if out_file is None:
        print(''.join(content).strip())
    for line in open(gtf_file):
        if not line.strip():
            continue
        if not line.startswith('#'):
            split = line.strip().split('\t')
            if split[2] == 'transcript':
                chr_num = split[0]
                start = int(split[3]) - target_range * 1000
                end = int(split[4]) + target_range * 1000
                strand = split[-3]
                attr_list = split[-1].split(' ')
                id_index = attr_list.index('transcript_id')
                transcript_id = attr_list[id_index + 1].replace(';', '').replace('''"''', '''''')
                if chr_num in gff_dict:
                    gene_list = gff_dict[chr_num]
                    for gene in gene_list:
                        dis = loc = None
                        if gene['start'] < start < gene['end'] < end:
                            dis, loc, strand = judge_distance_location(int(split[3]), int(split[4]), strand,
                                                                       gene['start'], gene['end'], gene['strand'])
                        elif start < gene['start'] < gene['end'] < end:
                            dis, loc, strand = judge_distance_location(int(split[3]), int(split[4]), strand,
                                                                       gene['start'], gene['end'], gene['strand'])
                        elif start < gene['start'] < end:
                            dis, loc, strand = judge_distance_location(int(split[3]), int(split[4]), strand,
                                                                       gene['start'], gene['end'], gene['strand'])
                        if dis is not None:
                            if out_file:
                                content.append(f"{chr_num}\t{transcript_id}\t{gene['id']}\t{dis}\t{loc}\t{strand}\n")
                            else:
                                print(f"{chr_num}\t{transcript_id}\t{gene['id']}\t{dis}\t{loc}\t{strand}")
    if out_file:
        with open(out_file, 'w') as o:
            o.write(''.join(content))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-g', '--gtf_file', 'gtf', help='LncRNA annotation file (format: GTF)')
@click.option('-a', '--gff_file', 'gff', help='Reference annotation file (format: GFF)')
@click.option('-f', '--feature_type', 'feature_type', default='mRNA', help='Feature type. {default: mRNA}')
@click.option('-d', '--distance', 'distance', type=int, default=100,
              help='Genes within a certain range of upstream and downstream of lncRNA were selected as target genes. '
                   '{default: 100 unit: kb}')
@click.option('-o', '--output_file', 'outfile',
              help='[optional] Output file, if not specified, print results to terminal as stdout.')
def run(gtf, gff, feature_type, distance, outfile):
    """
    Target genes were predicted according to the location relationship between lncRNA and genes.
    """
    lncRNA_target_gene_prediction(gtf, gff, feature_type, distance, outfile)


if __name__ == '__main__':
    run()
