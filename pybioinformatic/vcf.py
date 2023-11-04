"""
File: vcf.py
Description: Instantiate a VCF file object.
CreateDate: 2023/10/28
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union
from io import TextIOWrapper
from os.path import abspath
from gzip import GzipFile


class VCF:
    def __init__(self, path: Union[str, TextIOWrapper]):
        self.name = abspath(path) if isinstance(path, str) else abspath(path.name)
        if isinstance(path, str):
            if path.endswith('gz'):
                self.__open = GzipFile(path)
            else:
                self.__open = open(path)
        else:
            if path.name.endswith('gz'):
                self.__open = GzipFile(path.name)
            else:
                self.__open = path

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.__open.close()
        except AttributeError:
            pass

# Format conversion method==============================================================================================
    def to_genotype(self):
        """Convert VCF to GenoType."""
        for line in self.__open:
            if not isinstance(line, str):
                line = str(line, 'utf8')
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_name, position, ID, ref, alt, = split[:5]
                ID = f'{chr_name}_{position}'
                samples = split[9:]
                samples_gt = []
                for sample in samples:
                    gt = sample.split(':')[0]
                    if gt == './.':
                        gt = 'NA'
                    elif gt == '0/0':
                        gt = ref * 2
                    elif gt == '0/1':
                        gt = ''.join(sorted(ref + alt))
                    elif gt == '1/1':
                        gt = alt * 2
                    samples_gt.append(gt)
                samples_gt = '\t'.join(samples_gt)
                yield f'{ID}\t{chr_name}\t{position}\t{ref}\t{samples_gt}'
            elif not line.startswith('##') and line.startswith('#'):
                split = line.strip().split('\t')
                header = 'ID\tChrom\tPosition\tRef\t'
                samples = '\t'.join(split[9:])
                yield header + samples
