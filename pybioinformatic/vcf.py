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
    def to_genotype(self) -> str:
        """Convert VCF to GenoType."""
        for line in self.__open:
            if not isinstance(line, str):
                line = str(line, 'utf8')
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_name, position, ID, ref, alts, = split[:5]
                alts = alts.split(',')
                alts.insert(0, ref)
                ID = f'{chr_name}_{position}'
                samples = split[9:]
                samples_gt = []
                for sample in samples:
                    gt = sample.split(':')[0]
                    if gt == './.':
                        gt = 'NA'
                    else:
                        alleles_index = gt.replace('|', '/').split('/')
                        # SNP
                        if len(alts[int(alleles_index[0])]) == len(alts[int(alleles_index[1])]) == len(ref) == 1:
                            gt = alts[int(alleles_index[0])] + alts[int(alleles_index[1])]
                        # INDEL
                        else:
                            if len(ref) == len(alts[int(alleles_index[0])]) == 1:
                                alt = alts[int(alleles_index[1])][1:]
                                gt = f'-/ins{alt}'
                            elif len(ref) == len(alts[int(alleles_index[0])]) > 1:
                                alt = alts[int(alleles_index[0])][1:]
                                gt = f'-/del{alt}'
                            elif ref == alts[int(alleles_index[0])] == alts[int(alleles_index[1])]:
                                gt = '-'
                            else:
                                alt1 = f'del{ref.replace(alts[int(alleles_index[0])], "", 1)}' \
                                    if len(ref) > len(alts[int(alleles_index[0])]) else \
                                    f'ins{alts[int(alleles_index[0])].replace(ref, "", 1)}'
                                alt2 = f'del{ref.replace(alts[int(alleles_index[1])], "", 1)}' \
                                    if len(ref) > len(alts[int(alleles_index[1])]) else \
                                    f'ins{alts[int(alleles_index[1])].replace(ref, "", 1)}'
                                gt = f'{alt1}/{alt2}' if alt1 != alt2 else alt1
                    samples_gt.append(gt)
                samples_gt = '\t'.join(samples_gt)
                yield f'{ID}\t{chr_name}\t{position}\t{ref}\t{samples_gt}'
            elif not line.startswith('##') and line.startswith('#'):
                split = line.strip().split('\t')
                header = 'ID\tChrom\tPosition\tRef\t'
                samples = '\t'.join(split[9:])
                yield header + samples
