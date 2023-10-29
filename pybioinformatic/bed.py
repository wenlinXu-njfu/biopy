"""
File: bed.py
Description: Instantiate a BED file object.
CreateDate: 2021/12/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from os.path import abspath
from click import open_file
from pybioinformatic.sequence import Nucleotide
from pybioinformatic.fasta import Fasta


class Bed:
    def __init__(self, path: Union[str, TextIOWrapper]):
        if isinstance(path, str):
            self.name = abspath(path)
            self.__open = open(path)
            self.line_num = sum(1 for _ in self.__open)
            self.__open.seek(0)
        else:
            self.name = abspath(path.name)
            if path.name == '<stdin>':
                self.__open = open_file('-').readlines()
                self.line_num = sum(1 for _ in self.__open)
            else:
                self.__open = path
                self.line_num = sum(1 for _ in self.__open)
                self.__open.seek(0)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.__open.close()
        except AttributeError:
            pass

# Basic method==========================================================================================================
    def get_bed_dict(self) -> dict:
        bed_dict = {}  # {Chr1: [{start: int, end: int, id: str, frame: str, strand: str}, ...], ...}
        for line in self.__open:
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                item = {'start': int(split[1]), 'end': int(split[2]), 'id': split[3],
                        'frame': split[4], 'strand': split[5]}
                if split[0] in bed_dict:
                    bed_dict[split[0]].append(item)
                else:
                    bed_dict[split[0]] = [item]
        return bed_dict

# Sequence extraction method============================================================================================
    @staticmethod
    def __judge_range(start: int, end: int, strand: str, up: int = 0, down: int = 0, both: int = 0) -> tuple:
        if both:
            start -= both
            end += both
            if start < 0:
                start = 0
            return start, end
        else:
            if strand == '+':
                start -= up
                end += down
                if start < 0:
                    start = 0
                return start, end
            else:
                start -= down
                end += up
                return start, end

    def extract_seq(self,
                    fasta_file: Union[str, TextIOWrapper],
                    use_id: bool = True,
                    up: int = 0,
                    down: int = 0,
                    both: int = 0,
                    extension: bool = True) -> Nucleotide:
        """
        Extract the sequence in the BED file from the reference sequence file
        :param fasta_file: Reference sequence FASTA file
        :param use_id: Use fourth column content in the BED file as sequence ID. Otherwise,
                       "chr_num:start-end(strand)" by default.
        :param up: Make start site of sequence to extent to upstream specified length.
        :param down: Make end site of sequence to extent to downstream specified length.
        :param both: Make sequence to extend specified length both end.
        :param extension: Specify whether the extended sequence includes the sequence itself
        :return: Nucleotide_obj_generator
        """
        # bed_dict = {Chr_num: [{start: int, end: int, id: str, frame: str, strand: str}]}
        bed_dict = self.get_bed_dict()
        for nucl_obj in Fasta(fasta_file).parse():
            if nucl_obj.id in bed_dict:
                seqs: list = bed_dict[nucl_obj.id]  # [{start: int, end: int, id: str, frame: str, strand: str}]
                for d in seqs:
                    new_start, new_end = self.__judge_range(d['start'], d['end'], d['strand'], up, down, both)
                    if extension:
                        sub_seq_obj = nucl_obj[new_start:new_end]
                        if d['strand'] == '-':
                            sub_seq_obj = -sub_seq_obj
                        if use_id:
                            sub_seq_obj.id = d['id']
                        else:
                            sub_seq_obj.id = f"{nucl_obj.id}:{new_start}-{new_end}({d['strand']})"
                        yield sub_seq_obj
                    else:
                        if both:
                            if d['strand'] == '+':
                                up_seq_obj = nucl_obj[new_start:d['start']]
                                down_seq_obj = nucl_obj[d['end']:new_end]
                                if use_id:
                                    up_seq_obj.id = f"{d['id']} upstream_{both}_nt"
                                    down_seq_obj.id = f"{d['id']} downstream_{both}_nt"
                                else:
                                    up_seq_obj.id = f"{nucl_obj.id}:{new_start}-{d['start']}({d['strand']})"
                                    down_seq_obj.id = f"{nucl_obj.id}:{d['end']}-{new_end}({d['strand']})"
                            else:
                                up_seq_obj = nucl_obj[d['end']:new_end]
                                up_seq_obj = -up_seq_obj
                                down_seq_obj = nucl_obj[new_start:d['start']]
                                down_seq_obj = -down_seq_obj
                                if use_id:
                                    up_seq_obj.id = f"{d['id']} upstream_{both}_nt"
                                    down_seq_obj.id = f"{d['id']} downstream_{both}_nt"
                                else:
                                    up_seq_obj.id = f"{nucl_obj.id}:{d['end']}-{new_end}({d['strand']})"
                                    down_seq_obj.id = f"{nucl_obj.id}:{new_start}-{d['start']}({d['strand']})"
                            yield up_seq_obj
                            yield down_seq_obj
                        else:
                            if up:
                                if d['strand'] == '+':
                                    up_seq_obj = nucl_obj[new_start:d['start']]
                                    if use_id:
                                        up_seq_obj.id = f"{d['id']} upstream_{up}_nt"
                                    else:
                                        up_seq_obj.id = f"{nucl_obj.id}:{new_start}-{d['start']}({d['strand']})"
                                else:
                                    up_seq_obj = nucl_obj[d['end']:new_end]
                                    up_seq_obj = -up_seq_obj
                                    if use_id:
                                        up_seq_obj.id = f"{d['id']} upstream_{up}_nt"
                                    else:
                                        up_seq_obj.id = f"{nucl_obj.id}:{d['end']}-{new_end}({d['strand']})"
                                yield up_seq_obj
                            if down:
                                if d['strand'] == '+':
                                    down_seq_obj = nucl_obj[d['end']:new_end]
                                    if use_id:
                                        down_seq_obj.id = f"{d['id']} downstream_{down}_nt"
                                    else:
                                        down_seq_obj.id = f"{nucl_obj.id}:{d['end']}-{new_end}({d['strand']})"
                                else:
                                    down_seq_obj = nucl_obj[new_start:d['start']]
                                    down_seq_obj = -down_seq_obj
                                    if use_id:
                                        down_seq_obj.id = f"{d['id']} downstream_{down}_nt"
                                    else:
                                        down_seq_obj.id = f"{nucl_obj.id}:{new_start}-{d['start']}({d['strand']})"
                                yield down_seq_obj
