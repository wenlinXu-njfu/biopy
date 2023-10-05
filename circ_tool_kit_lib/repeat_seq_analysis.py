#!/usr/bin/env python
"""
File: repeat_seq_analysis.py
Description: Repeat sequence analysis of circRNAs flanking sequences.
Date: 2023/2/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib import Gff, Bed, Fasta, Displayer
displayer = Displayer(__file__.split('/')[-1], version='1.0.0')


def main(circ_bed_file, repeat_seq_gff_file, genome_fasta_file, distance: int, out_file):
    repeat_loci: dict = Gff(repeat_seq_gff_file).to_dict()
    repeat_class: dict = {line.strip().split('\t')[8].replace('=', ';').split(';')[1]:
                              line.strip().split('\t')[8].replace('=', ';').split(';')[-1]
                          for line in open(repeat_seq_gff_file) if not line.startswith('#')}
    circ_loci: dict = Bed(circ_bed_file).get_bed_dict()
    content = '#Chr_num\tCirc_id\tCirc_start\tCirc_end\tCirc_strand\tRepeat_seq_id\tRepeat_seq_start\tRepeat_seq_end\t' \
              'Repeat_seq_strand\tRepeat_seq_class\tRepeat_seq\n'
    if not out_file:
        print(content)
    for chromosome in Fasta(genome_fasta_file).parse():
        try:
            circ_list = circ_loci[chromosome.id]
            repeats = repeat_loci[chromosome.id]
        except KeyError:
            pass
        else:
            for circ in circ_list:
                for repeat in repeats:
                    if repeat['strand'] == circ['strand']:
                        if 0 < circ['start'] - repeat['end'] <= distance or 0 < repeat['start'] - circ['end'] <= distance:
                            repeat_seq = chromosome.seq[repeat['start'] - 1:repeat['end']]
                            if out_file:
                                content += '\t'.join(
                                    [chromosome.id, circ['id'], str(circ['start']), str(circ['end']), circ['strand'],
                                     repeat['id'], str(repeat['start']), str(repeat['end']), repeat['strand'],
                                     repeat_class[repeat['id']], repeat_seq]
                                ) + '\n'
                            else:
                                print(chromosome.id, circ['id'], circ['start'], circ['end'], circ['strand'],
                                      repeat['id'],
                                      repeat['start'], repeat['end'], repeat['strand'], repeat_class[repeat['id']],
                                      repeat_seq, sep='\t')
    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-b', '--circ_bed_file', 'circ_bed_file',
              metavar='<bed file>', required=True,
              help='Input circRNA BED file\n(Chr_num\\tStart\\tEnd\\tID\\tFrame\\tStrand).')
@click.option('-g', '--gff_file', 'repeat_gff',
              metavar='<gff file>', required=True,
              help='Input repeat sequence GFF file generated from RepeatMasker software.')
@click.option('-f', '--genome_fasta', 'genome_fasta',
              metavar='<fasta file>', required=True,
              help='Input genome FASTA file.')
@click.option('-d', '--distance', 'distance',
              metavar='<int>', type=int, default=100, show_default=True,
              help='Specify the minimum distance between circRNA and repeat sequence.')
@click.option('-o', '--output_file', 'outfile',
              metavar='<file>',
              help='Output file, if not specified, print results to terminal as stdout.')
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(circ_bed_file, repeat_gff, genome_fasta, distance, outfile):
    """Repeat sequence analysis of circRNAs flanking sequences."""
    main(circ_bed_file, repeat_gff, genome_fasta, distance, outfile)


if __name__ == '__main__':
    run()
