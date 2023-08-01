#!/usr/bin/env python
"""
File: extract_circRNA.py
Description: Extract circRNA sequence.
Date: 2022/10/29
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
from Biolib.gtf import Gtf
from Biolib.bed import Bed
from Biolib.fasta import Fasta
from Biolib.sequence import Nucleotide
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(genome_fasta_file, genome_gtf_file, circRNA_bed_file, out_file):
    # Step 1: Read in exon, gene and circRNA locus information.
    gtf_file = Gtf(genome_gtf_file)
    exon_dict = gtf_file.get_non_redundant_exon()  # {chr_num: [{gene_id: str, start: int, end: int, strand: str}, ...], ...}
    gene_dict = {}  # {chr_num: [{gene_id: str, start: int, end: int, strand: str}, ...], ...}
    for line in gtf_file.parse():
        if line[2] == 'gene':
            chr_num, start, end, strand, gene_id = line[0], int(line[3]), int(line[4]), line[6], line[8]
            if chr_num in gene_dict:
                gene_dict[chr_num].append({'gene_id': gene_id, 'start': start, 'end': end, 'strand': strand})
            else:
                gene_dict[chr_num] = [{'gene_id': gene_id, 'start': start, 'end': end, 'strand': strand}]
    circ_dict = Bed(circRNA_bed_file).get_bed_dict()  # {chr_num: [{circ_id: str, start: int, end: int, strand: str}, ...], ...}

    # Step 2: Extract circRNA sequence.
    content = ''
    for nucl in Fasta(genome_fasta_file).parse():
        try:
            circs: list = circ_dict[nucl.id]
            exons: list = exon_dict[nucl.id]
            genes: list = gene_dict[nucl.id]
        except KeyError:
            pass
        else:
            for circ in circs:
                circ_id = circ['id']
                circ_seq = ''
                circ_type = None
                host_gene = ''
                start = end = None  # Mark whether the BSJ site is boundary of exon
                exon_number = 0

                for exon in exons:
                    if exon['strand'] == circ['strand']:
                        # single exon circRNA
                        if exon['start'] <= circ['start'] <= circ['end'] <= exon['end']:
                            circ_seq += nucl.seq[circ['start'] - 1:circ['end']]
                            exon_number += 1
                            host_gene = exon['id']
                            start = end = True
                            break
                        # multiple exon circRNA
                        elif exon['start'] <= circ['start'] < exon['end'] < circ['end']:
                            circ_seq = nucl.seq[exon['start'] - 1:exon['end']]
                            exon_number += 1
                            if exon['id'] not in host_gene and host_gene:
                                host_gene += f";{exon['id']}"
                            elif not host_gene:
                                host_gene = exon['id']
                            start = True
                        elif circ['start'] < exon['start'] < exon['end'] < circ['end']:
                            circ_seq += nucl.seq[exon['start'] - 1:exon['end']]
                            exon_number += 1
                            if exon['id'] not in host_gene and host_gene:
                                host_gene += f";{exon['id']}"
                            elif not host_gene:
                                host_gene = exon['id']
                        elif circ['start'] < exon['start'] < circ['end'] <= exon['end']:
                            circ_seq += nucl.seq[exon['start'] - 1:circ['end']]
                            exon_number += 1
                            if exon['id'] not in host_gene and host_gene:
                                host_gene += f";{exon['id']}"
                            elif not host_gene:
                                host_gene = exon['id']
                            end = True
                            break
                # exonic circRNA
                if start and end:
                    if len(host_gene.split(';')) == 2:
                        circ_type = 'read_through'
                    elif len(host_gene.split(';')) == 1:
                        circ_type = 'exonic'
                    circ_id += f' type={circ_type} exon_number={exon_number} host_gene={host_gene}'
                elif exon_number == 0:
                    for gene in genes:
                        # single intronic circRNA
                        if gene['strand'] == circ['strand']:
                            if gene['start'] < circ['start'] < circ['end'] < gene['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'intronic'
                                host_gene = gene['gene_id']
                                circ_id += f' type={circ_type} host_gene={host_gene}'
                                break
                        # antisense circRNA
                        else:
                            if circ['start'] < gene['start'] < circ['end'] < gene['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'antisense'
                                antisense_gene = gene['gene_id']
                                circ_id += f' type={circ_type} antisense_gene={antisense_gene}'
                                break
                            elif gene['start'] <= circ['start'] < circ['end'] <= gene['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'antisense'
                                antisense_gene = gene['gene_id']
                                circ_id += f' type={circ_type} antisense_gene={antisense_gene}'
                                break
                            elif gene['start'] < circ['start'] <= gene['end'] < circ['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'antisense'
                                antisense_gene = gene['gene_id']
                                circ_id += f' type={circ_type} antisense_gene={antisense_gene}'
                                break
                            elif circ['start'] < gene['start'] < gene['end'] < circ['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'antisense'
                                antisense_gene = gene['gene_id']
                                circ_id += f' type={circ_type} antisense_gene={antisense_gene}'
                                break
                    # intergenic circRNA
                    if not circ_type:
                        circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                        circ_type = 'intergenic'
                        circ_id += f' type={circ_type}'
                # single exonic-intronic circRNA
                elif start or end:
                    if exon_number == 1 and len(host_gene.split(';')) == 1:
                        for gene in genes:
                            if gene['strand'] == circ['strand'] and gene['start'] <= circ['start'] <= circ['end'] <= gene['end']:
                                circ_seq = nucl.seq[circ['start'] - 1:circ['end']]
                                circ_type = 'exonic_intronic'
                                circ_id += f' type={circ_type} host_gene={host_gene}'

                if circ_type:
                    circ_obj = Nucleotide(circ_id, circ_seq)
                    if circ['strand'] == '-':
                        circ_obj = -circ_obj
                        circ_obj.id = circ_obj.id.replace(' reverse_complementary_chain', '')
                    if not out_file:
                        print(f'>{circ_obj.id}\n{circ_obj.seq}')
                    else:
                        content += f'>{circ_obj.id}\n{circ_obj.seq}\n'

    if out_file:
        with open(out_file, 'w') as o:
            o.write(content)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--genome_fasta', 'genome_fasta', help='Input genome FASTA file.')
@click.option('-a', '--gtf_file', 'gtf_file', help='Input genome GTF file.')
@click.option('-b', '--circ_bed_file', 'circ_bed_file', help='Input circRNA locus BED file, the site of first nucleotide is one. '
                                                             '(Chr_num\\tStart\\tEnd\\tCircRNA_ID\\tFrame\\tStrand\\n)')
@click.option('-o', '--output_file', 'outfile',
              help='[optional] Output FASTA file, if not specified, print results to terminal as stdout.')
def run(genome_fasta, gtf_file, circ_bed_file, outfile):
    """Extract circRNA sequence."""
    main(genome_fasta, gtf_file, circ_bed_file, outfile)


if __name__ == '__main__':
    run()
