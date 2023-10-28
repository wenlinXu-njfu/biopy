#!/usr/bin/env python
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
from natsort import natsort_key
from pandas import read_table, Series, DataFrame, concat
from pybioinformatic import TaskManager


def __record_to_genotype(row: Series):
    ref, alt, record = row[0], row[1], row[2]
    gt = record.split(':')[0]
    if gt == '0/0':
        allele = ref * 2
    elif gt == '0/1':
        allele = ref + alt
    elif gt == '1/1':
        allele = alt * 2
    else:
        allele = 'NA'
    return allele


def vcf_to_genotype(df: DataFrame):
    new_col = ['Chrom', 'Position', 'Ref', 'Alt']
    for sample in df.columns[5:]:
        new_col.append(sample)
        df[sample] = df.loc[:, ['REF', 'ALT', sample]].apply(__record_to_genotype, axis=1)
    df['ID'] = df['#CHROM'] + '_' + df['POS']
    df.set_index('ID', drop=True, inplace=True)
    df.columns = new_col
    df.drop(columns='Alt', inplace=True)
    return df


class VCF:
    def __init__(self, path: Union[str, TextIOWrapper]):
        self.name = abspath(path) if isinstance(path, str) else abspath(path.name.replace('<', '').replace('>', ''))
        self.__open = open(path) if isinstance(path, str) else path

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.__open.close()
        except AttributeError:
            pass

    @staticmethod
    def __use_col(column_name: str):
        except_col = {'QUAL', 'FILTER', 'INFO', 'FORMAT'}
        return True if column_name not in except_col else False

    def to_dataframe(self):
        if 'gz' in self.name:  # read xxx.vcf.gz file
            df = read_table(self.name, usecols=self.__use_col, dtype=str,
                            skiprows=sum(1 for line in GzipFile(self.name) if str(line, 'utf8').startswith('##')))
        else:
            skip_rows = sum(1 for line in self.__open if line.startswith('##'))
            self.__open.seek(0)
            df = read_table(self.__open, usecols=self.__use_col, dtype=str, skiprows=skip_rows)
        return df

    def to_dataframes(self, chunk_size: int = 10000):
        if 'gz' in self.name:  # read xxx.vcf.gz file
            dfs = read_table(self.name, usecols=self.__use_col, dtype=str, chunksize=chunk_size,
                             skiprows=sum(1 for line in GzipFile(self.name) if str(line, 'utf8').startswith('##')))
        else:
            skip_rows = sum(1 for line in self.__open if line.startswith('##'))
            self.__open.seek(0)
            dfs = read_table(self.__open, usecols=self.__use_col, dtype=str, skiprows=skip_rows, chunksize=chunk_size)
        return dfs

    def to_genotype(self, num_processing: int):
        params = [(df,) for df in self.to_dataframes()]
        tkm = TaskManager(processing_num=num_processing, params=params)
        results = tkm.parallel_run_func(vcf_to_genotype)
        results = [i.get() for i in results]
        gt = concat(results)
        gt.sort_index(key=natsort_key, inplace=True)
        return gt
