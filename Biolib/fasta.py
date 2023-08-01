"""
File: fasta.py
Description: Instantiate a FASTA file object
Date: 2021/11/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import re
from typing import Dict
from gzip import GzipFile
import click
from itertools import groupby
from Biolib.sequence import Nucleotide, Protein


class Fasta:
    def __init__(self, path: str):
        self.path = path

# Basic method==========================================================================================================
    def parse(self, parse_id: bool = True) -> Nucleotide:  # return Nucleotide generator
        """
        A FASTA file generator that returns one Nucleotide object at a time.
        """
        try:
            fa_generator = (ret[1] for ret in groupby(open(self.path), lambda line: line.startswith('>')))
            for g in fa_generator:
                seq_id = g.__next__().strip()
                seq = ''.join(line.strip() for line in fa_generator.__next__())
                if parse_id:
                    if '|' in seq_id:
                        seq_id = seq_id.split('|')[0]
                    else:
                        seq_id = seq_id.split(' ')[0]
                if 'M' not in seq and '*' not in seq:
                    yield Nucleotide(seq_id, seq)
                else:
                    yield Protein(seq_id, seq)
        except UnicodeDecodeError:
            l = []
            for line in GzipFile(self.path):
                l.append(str(line, 'utf8'))
            fa_generator = (ret[1] for ret in groupby(l, lambda line: line.startswith('>')))
            for g in fa_generator:
                seq_id = g.__next__().strip()
                seq = ''.join(line.strip() for line in fa_generator.__next__())
                if parse_id:
                    if '|' in seq_id:
                        seq_id = seq_id.split('|')[0]
                    else:
                        seq_id = seq_id.split(' ')[0]
                if 'M' not in seq and '*' not in seq:
                    yield Nucleotide(seq_id, seq)
                else:
                    yield Protein(seq_id, seq)

    def get_seq_dict(self, parse_id: bool = False) -> dict:
        """Get sequence dict from FASTA file"""
        seq_dict = {}
        for nucl_obj in self.parse(parse_id):
            if nucl_obj.id not in seq_dict:
                seq_dict[nucl_obj.id] = nucl_obj.seq
            else:
                click.echo(f'\033[31mError: FASTA file has repeat id {nucl_obj.id}.', err=True)
                exit()
        return seq_dict

# File format conversion method=========================================================================================
    def check_FASTA(self) -> bool:
        """Check whether a file is formal FASTA format."""
        try:
            with open(self.path) as f:
                f.readline()
                f.readline()
                line = f.readline()
                if line.startswith('>'):
                    return True
                else:
                    return False
        except UnicodeDecodeError:
            with GzipFile(self.path) as f:
                f.readline()
                f.readline()
                line = str(f.readline(), 'utf8')
                if line.startswith('>'):
                    return True
                else:
                    return False

    def merge_sequence(self) -> str:  # return str generator
        """Make each sequence to be displayed on a single line."""
        is_fa = self.check_FASTA()
        if not is_fa:
            fa_generator = self.parse(False)
            for nucl_obj in fa_generator:
                yield f">{nucl_obj.id}\n{nucl_obj.seq}\n"
        else:
            click.echo('\033[33mThe input FASTA file does not need to be formatted.\033[0m', err=True)
            exit()

    def split_sequence(self, char_num: int) -> str:  # return str generator
        """Make each sequence to be displayed in multiple lines."""
        fa_generator = self.parse(False)
        for nucl_obj in fa_generator:
            nucl_obj = nucl_obj.display_set(char_num)
            yield f">{nucl_obj.id}\n{nucl_obj.seq}"

# Other method==========================================================================================================
    def get_longest_seq(self, regular_exp: str = r'\w+.\w+', inplace_id: bool = False) -> Dict[str, str]:
        """Get the longest transcript of each gene."""
        all_seq_dict = {}  # {seq_id: seq}
        longest_seq_dict = {}
        id_map_dict = {}  # {'Potri.001G000100': 'Potri.001G000100.3', 'Potri.001G000200': 'Potri.001G000200.1', ...}
        if inplace_id:
            for seq_obj in self.parse():
                all_seq_dict[seq_obj.id] = seq_obj.seq
                gene_id = re.findall(regular_exp, seq_obj.id)[0]
                if gene_id not in id_map_dict:
                    id_map_dict[gene_id] = seq_obj.id
                else:
                    if len(seq_obj) >= len(all_seq_dict[id_map_dict[gene_id]]):
                        id_map_dict[gene_id] = seq_obj.id
            for locus, longest_seq_id in id_map_dict.items():
                longest_seq_dict[locus] = all_seq_dict[longest_seq_id]
        else:
            for seq_obj in self.parse(False):
                all_seq_dict[seq_obj.id] = seq_obj.seq
                gene_id = re.findall(regular_exp, seq_obj.id)[0]
                if gene_id not in id_map_dict:
                    id_map_dict[gene_id] = seq_obj.id
                else:
                    if seq_obj.len >= len(all_seq_dict[seq_obj.id]):
                        id_map_dict[gene_id] = seq_obj.id
            for locus, longest_seq_id in id_map_dict.items():
                longest_seq_dict[longest_seq_id] = all_seq_dict[longest_seq_id]
        return longest_seq_dict

    def filter_n(self, max_num=1) -> Nucleotide:  # return Nucleotide object generator
        for nucl_obj in self.parse():
            if nucl_obj.seq.count('n') < max_num or nucl_obj.seq.count('N') < max_num:
                yield nucl_obj
            else:
                click.echo(f'{nucl_obj.id} has been filtered out.', err=True)
