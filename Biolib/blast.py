"""
File: blast.py
Description: Instantiate a blast result file object
Date: 2021/12/4
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, IO, Tuple


class Blast:
    def __init__(self, path: Union[IO, str]):
        self.path = path
        self.line_num = sum(1 for _ in open(self.path))

# Basic method==========================================================================================================
    def parse_BLAST(self) -> Tuple[str]:
        """Parse information of each column of BLAST6 file line by line."""
        for line in open(self.path):
            split = line.strip().split('\t')
            querry_id, sbject_id, align_rate = split[0], split[1], split[2]
            align_len, mismatch, gap = split[3], split[4], split[5]
            querry_start, querry_end = split[6], split[7]
            sbject_start, sbject_end = split[8], split[9]
            e_value, score = split[10], split[11]
            yield querry_id, sbject_id, align_rate, align_len, querry_start, querry_end, sbject_start, sbject_end, e_value, score

    def get_pair_dict(self, top: int = 3):
        pair_dict = {}  # {query1: {sbject1: [], sbject2: [], ...}, query2: {}, ...}
        for line in self.parse_BLAST():
            if line[0] not in pair_dict:
                pair_dict[line[0]] = {line[1]: '\t'.join(line[2:])}
            else:
                if len(pair_dict[line[0]]) < top:
                    pair_dict[line[0]][line[1]] = '\t'.join(line[2:])
        return pair_dict

# File format conversion method=========================================================================================
    def blast_to_BED(self, query_is_chr: bool = False) -> str:
        """Transform the result of align the query sequence with the reference sequence into a BED file."""
        content = ''
        for line in self.parse_BLAST():
            if not query_is_chr:
                if int(line[6]) - int(line[7]) > 0:
                    content += f'{line[1]}\t{int(line[7]) - 1}\t{line[6]}\t{line[0]}\t.\t-\n'
                else:
                    content += f'{line[1]}\t{int(line[6]) - 1}\t{line[7]}\t{line[0]}\t.\t+\n'
            else:
                if int(line[4]) - int(line[5]) > 0:
                    content += f'{line[0]}\t{int(line[5]) - 1}\t{line[4]}\t{line[1]}\t.\t-\n'
                else:
                    content += f'{line[0]}\t{int(line[4]) - 1}\t{line[5]}\t{line[1]}\t.\t+\n'
        return content
