"""
File: fasta.py
Description: Instantiate a FASTA file object.
Date: 2021/11/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from re import findall
from _io import TextIOWrapper
from typing import Union
from gzip import GzipFile
from click import echo, open_file
from itertools import groupby
from Biolib.sequence import Nucleotide, Protein


class Fasta:
    def __init__(self, path: Union[str, TextIOWrapper]):
        if isinstance(path, str):
            self.path = path
        else:
            if path.name == '<stdin>':
                self.path = open_file('-').readlines()
            else:
                self.path = path.name

# Basic method==========================================================================================================
    def parse(self, parse_id: bool = True) -> Nucleotide:  # return Nucleotide generator
        """A FASTA file generator that returns one Nucleotide or Protein object at one time."""
        # Parse FASTA format from a file.
        if isinstance(self.path, str):
            # Parse uncompressed FASTA file (xx.fa).
            try:
                fa_generator = (ret[1] for ret in groupby(open(self.path), lambda line: line.startswith('>')))
                for g in fa_generator:
                    seq_id = g.__next__().strip()
                    seq = ''.join(line.strip() for line in fa_generator.__next__())
                    if parse_id:
                        if '\t' in seq_id:
                            seq_id = seq_id.split('\t')[0]
                        elif '|' in seq_id:
                            seq_id = seq_id.split('|')[0]
                        else:
                            seq_id = seq_id.split(' ')[0]
                    if 'M' not in seq and '*' not in seq:
                        yield Nucleotide(seq_id, seq)
                    else:
                        yield Protein(seq_id, seq)
            # Parse compressed FASTA file (xx.fa.gz).
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
        # Parse FASTA format from command line stdin.
        else:
            # Parse uncompressed FASTA file (xx.fa).
            try:
                fa_generator = (ret[1] for ret in groupby(self.path, lambda line: line.startswith('>')))
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
            # Parse compressed FASTA file (xx.fa.gz).
            except UnicodeDecodeError:
                l = []
                for line in GzipFile(self.path.name):
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
        """Get sequence dict from FASTA file."""
        seq_dict = {}
        for nucl_obj in self.parse(parse_id):
            if nucl_obj.id not in seq_dict:
                seq_dict[nucl_obj.id] = nucl_obj.seq
            else:
                echo(f'\033[31mError: FASTA file has repeat id {nucl_obj.id}.', err=True)
                exit()
        return seq_dict

# File format conversion method=========================================================================================
    def check_FASTA(self) -> bool:
        """Check whether a file is formal FASTA format."""
        if isinstance(self.path, str):
            try:
                with open(self.path) as f:
                    f.readline()
                    f.readline()
                    line = f.readline()
                    return True if line.startswith('>') else False
            except UnicodeDecodeError:
                with GzipFile(self.path) as f:
                    f.readline()
                    f.readline()
                    line = str(f.readline(), 'utf8')
                    return True if line.startswith('>') else False
        else:
            try:
                with self.path as f:
                    f.readline()
                    f.readline()
                    line = f.readline()
                    return True if line.startswith('>') else False
            except UnicodeDecodeError:
                with GzipFile(self.path.name) as f:
                    f.readline()
                    f.readline()
                    line = str(f.readline(), 'utf8')
                    return True if line.startswith('>') else False

    def merge_sequence(self) -> Union[Nucleotide, Protein]:
        """Make each sequence to be displayed on a single line."""
        is_fa = self.check_FASTA()
        if not is_fa:
            for seq_obj in self.parse(False):
                yield seq_obj
        else:
            echo('\033[33mThe input FASTA file does not need to be formatted.\033[0m', err=True)
            exit()

    def split_sequence(self, char_num: int) -> Union[Nucleotide, Protein]:
        """Make each sequence to be displayed in multiple lines."""
        for seq_obj in self.parse(False):
            seq_obj = seq_obj.display_set(char_num)
            yield seq_obj

# Other method==========================================================================================================
    def get_longest_seq(self, regular_exp: str = r'\w+.\w+', inplace_id: bool = False) -> Union[Nucleotide, Protein]:
        """Get the longest transcript of each gene locus."""
        all_seq_dict = {seq_obj.id: seq_obj.seq for seq_obj in self.parse(False)}  # {seq_id: seq}
        longest_seq_dict = {}  # {locus_id: seq}
        id_map_dict = {}  # {'Potri.001G000100': 'Potri.001G000100.3', 'Potri.001G000200': 'Potri.001G000200.1', ...}
        if inplace_id:
            for seq_obj in self.parse(False):
                gene_id = findall(regular_exp, seq_obj.id)[0]
                if gene_id not in id_map_dict:
                    id_map_dict[gene_id] = seq_obj.id
                else:
                    if len(seq_obj) >= len(all_seq_dict[id_map_dict[gene_id]]):
                        id_map_dict[gene_id] = seq_obj.id
            for locus, longest_seq_id in id_map_dict.items():
                longest_seq_dict[locus] = all_seq_dict[longest_seq_id]
        else:
            for seq_obj in self.parse(False):
                gene_id = findall(regular_exp, seq_obj.id)[0]
                if gene_id not in id_map_dict:
                    id_map_dict[gene_id] = seq_obj.id
                else:
                    if seq_obj.len >= len(all_seq_dict[id_map_dict[gene_id]]):
                        id_map_dict[gene_id] = seq_obj.id
            for locus, longest_seq_id in id_map_dict.items():
                longest_seq_dict[longest_seq_id] = all_seq_dict[longest_seq_id]
        for seq_id, seq in longest_seq_dict.items():
            yield Nucleotide(seq_id, seq) if 'M' not in seq and '*' not in seq else Protein(seq_id, seq)

    def filter_n(self, max_num=1) -> Nucleotide:
        for nucl_obj in self.parse():
            if nucl_obj.seq.count('n') < max_num or nucl_obj.seq.count('N') < max_num:
                yield nucl_obj
            else:
                echo(f'{nucl_obj.id} has been filtered out.', err=True)
