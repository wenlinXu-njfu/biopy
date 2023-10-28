"""
File: fasta.py
Description: Instantiate a FASTA file object.
Date: 2021/11/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from re import findall
from io import TextIOWrapper
from typing import Union
from os.path import abspath
from gzip import GzipFile
from click import echo, open_file
from itertools import groupby
from pybioinformatic.sequence import Nucleotide, Protein


class Fasta:
    def __init__(self, path: Union[str, TextIOWrapper]):
        if isinstance(path, str):
            self.name = abspath(path)
            if path.endswith('gz'):
                self.__open = GzipFile(path)
                self.seq_num = sum(1 for line in self.__open if str(line, 'utf8').startswith('>'))
                self.__open.seek(0)
            else:
                self.__open = open(path)
                self.seq_num = sum(1 for line in open(path) if line.startswith('>'))
        else:
            if path.name == '<stdin>':
                self.name = 'stdin'
                self.__open = open_file('-').readlines()
                self.seq_num = sum(1 for line in self.__open if line.startswith('>'))
            else:
                self.name = abspath(path.name)
                if path.name.endswith('gz'):
                    self.__open = GzipFile(path.name)
                    self.seq_num = sum(1 for line in self.__open if str(line, 'utf8').startswith('>'))
                    self.__open.seek(0)
                else:
                    self.__open = path
                    self.seq_num = sum(1 for line in self.__open if line.startswith('>'))
                    self.__open.seek(0)

# Basic method==========================================================================================================
    def parse(self, parse_id: bool = True) -> Nucleotide:  # return Nucleotide generator
        """A FASTA file generator that returns one Nucleotide or Protein object at one time."""
        if self.name.endswith('gz'):
            fa_generator = (ret[1] for ret in groupby(self.__open, lambda line: str(line, 'utf8').startswith('>')))
        else:
            fa_generator = (ret[1] for ret in groupby(self.__open, lambda line: line.startswith('>')))
        for g in fa_generator:
            if self.name.endswith('gz'):
                seq_id = str(g.__next__(), 'utf8').strip()
                seq = ''.join(str(line, 'utf8').strip() for line in fa_generator.__next__())
            else:
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

    def to_dict(self, parse_id: bool = False) -> dict:
        """Get sequence dict from FASTA file."""
        seq_dict = {}
        for nucl_obj in self.parse(parse_id):
            if nucl_obj.id not in seq_dict:
                seq_dict[nucl_obj.id] = nucl_obj.seq
            else:
                echo(f'\033[31mError: FASTA file has repeat id {nucl_obj.id}.', err=True)
                exit()
        return seq_dict

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.__open.close()
        except AttributeError:
            pass

# File format conversion method=========================================================================================
    def merge_sequence(self, parse_id: bool = False) -> Union[Nucleotide, Protein]:
        """Make each sequence to be displayed on a single line."""
        for seq_obj in self.parse(parse_id):
            yield seq_obj

    def split_sequence(self, parse_id: bool = False, char_num: int = 60) -> Union[Nucleotide, Protein]:
        """Make each sequence to be displayed in multiple lines."""
        for seq_obj in self.parse(parse_id):
            seq_obj = seq_obj.display_set(char_num)
            yield seq_obj

    def fa2tab(self, parse_id: bool = False):
        for seq_obj in self.parse(parse_id):
            yield f'{seq_obj.id}\t{seq_obj.seq}'

# Other method==========================================================================================================
    def get_longest_seq(self, regular_exp: str = r'\w+.\w+', inplace_id: bool = False) -> Union[Nucleotide, Protein]:
        """Get the longest transcript of each gene locus."""
        all_seq_dict = self.to_dict(False)  # {seq_id: seq}
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
