#!/usr/bin/env python
"""
File: reverse_complementary_analysis.py
Description: Reverse complementary analysis of circRNAs flanking sequences.
Date: 2022/4/2
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from os import system, mkdir
from datetime import datetime
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(bed_file, ref_seq_file, flanking_seq_len: int, out_file):
    click.echo(f"[{datetime.now().replace(microsecond=0)}] Extract flanking sequence.", err=True)
    system(command=f"seq_extraction_tool bed -i {bed_file} -f {ref_seq_file} -id -both {flanking_seq_len} "
                   f"-non_extension -o flanking.fa")
    click.echo(f"[{datetime.now().replace(microsecond=0)}] Analysis flanking sequence with blastn.", err=True)
    with open('flanking.fa') as f:
        while True:
            seq1_id = f.readline()
            seq1 = f.readline()
            seq2_id = f.readline()
            seq2 = f.readline()
            if not seq1_id:
                break
            else:
                mkdir('temp')
                with open('temp/1.fa', 'w') as o1:
                    o1.write(seq1_id + seq1)
                with open('temp/2.fa', 'w') as o2:
                    o2.write(seq2_id + seq2)
                system(command="makeblastdb -in temp/1.fa -dbtype nucl -out temp/1.fa -parse_seqids -logfile temp/1.fa.log")
                system(command=f"blastn -task blastn-short -strand minus -query temp/2.fa -db temp/1.fa -outfmt 6 >> {out_file}")
                system(command="rm -r temp/")
    click.echo(f"[{datetime.now().replace(microsecond=0)}] Done, please see results in {out_file}.", err=True)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--circ_bed_file', 'circ_bed_file',
              help='Input circRNA BED file (Chr_num\\tStart\\tEnd\\tID\\tFrame\\tStrand\\n).')
@click.option('-f', '--ref_fasta_file', 'ref_fasta_file', help='Input reference sequence FASTA file.')
@click.option('-l', '--flank_seq_len', 'flank_seq_len', type=int, default=100,
              help='[optional] Specify length of flanking sequence. {default=100}')
@click.option('-o', '--output_file', 'outfile', help='Output file.')
def run(circ_bed_file, ref_fasta_file, flank_seq_len, outfile):
    """Reverse complementary analysis of circRNAs flanking sequences."""
    main(circ_bed_file, ref_fasta_file, flank_seq_len, outfile)


if __name__ == '__main__':
    run()
