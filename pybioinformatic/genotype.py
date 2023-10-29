"""
File: genotype.py
Description: Instantiate a GT file object.
CreateDate: 2023/10/26
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, List
from io import TextIOWrapper
from os.path import abspath
from warnings import filterwarnings
from tqdm import tqdm
from natsort import natsort_key
from pandas import Series, DataFrame, read_table, read_excel, concat, cut
from click import echo
from pybioinformatic.task_manager import TaskManager
from pybioinformatic.biopandas import read_file_as_dataframe_from_stdin
filterwarnings("ignore")


# =====================================================================================================================#
# NOTICE:                                                                                                              #
# The "__check_hom" and "__get_all_allele" functions cannot be defined in the "GenoType" class.                        #
# Because the "GenoType.parallel_stat_MHM" method runs with multiprocessing,                                           #
# the subprocess cannot call the private method of the "GenoType" class.                                               #
# =====================================================================================================================#
def __check_hom(row: Series) -> Series:
    """Check which samples are homozygous genotypes at specific loci."""
    genotype_set = {'AT', 'TA',
                    'AG', 'GA',
                    'AC', 'CA',
                    'GC', 'CG',
                    'GT', 'TG',
                    'CT', 'TC'}
    data = [str(value) not in genotype_set and '/' not in str(value) and str(value) != ''
            for value in row]
    return Series(data, row.index)


def __get_all_allele(row: Series) -> str:
    """Get all allele at specified loci."""
    all_allele = ''
    for value in row:
        if '/' in str(value):
            left, right = str(value).split('/')[0], str(value).split('/')[1]
            if 'ins' in left:
                all_allele += 'I'
            elif 'del' in left:
                all_allele += 'D'
            else:
                all_allele += left[0]
            if 'ins' in right:
                all_allele += 'I'
            elif 'del' in right:
                all_allele += 'D'
            else:
                all_allele += right[0]
        elif 'ins' in str(value):
            all_allele += 'II'
        elif 'del' in str(value):
            all_allele += 'DD'
        else:
            all_allele += str(value)
    return all_allele


def stat_MHM(df: DataFrame) -> DataFrame:
    """
    Calculate MissRate, HetRate, and MAF (MHM) from specified DataFrame.
    The top 4 columns of DataFrame must be SNP ID, chromosome, position, and reference sequence info respectively,
    and the header of DataFrame is not None.
    """
    snp_ref = df.iloc[:, :4]
    # Calculate MissRate
    sample_num = len(df.columns.tolist()) - 4
    df['MissRate(%)'] = df.isnull().sum(axis=1) / sample_num * 100
    # Calculate HetRate
    df.fillna('', inplace=True)  # fill NA
    df['HetRate(%)'] = df.iloc[:, 4:-1].apply(lambda row: (1 - __check_hom(row).sum() / (row != '').sum()) * 100, axis=1)
    # Calculate MAF
    df['all_gt'] = df.iloc[:, 4:-2].apply(lambda row: __get_all_allele(row), axis=1)
    df['total'] = df['all_gt'].apply(len)
    df['A'] = df['all_gt'].str.count('A') / df['total']
    df['G'] = df['all_gt'].str.count('G') / df['total']
    df['C'] = df['all_gt'].str.count('C') / df['total']
    df['T'] = df['all_gt'].str.count('T') / df['total']
    df['D'] = df['all_gt'].str.count('D') / df['total']
    df['I'] = df['all_gt'].str.count('I') / df['total']
    df['MAF'] = df.loc[:, ['A', 'G', 'C', 'T', 'D', 'I']].apply(lambda row: sorted(row, reverse=True)[1], axis=1)
    # Output results
    df = df.loc[:, ['MissRate(%)', 'HetRate(%)', 'MAF']]
    merge = snp_ref.join(df)
    return merge


class GenoType:
    """
    Standard genotype file object.
    """

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
    def __allele_sort(allele: str) -> str:
        if '/' in allele:
            return '/'.join(sorted(allele.split('/')))
        elif 'ins' in allele or 'del' in allele:
            return allele
        else:
            return ''.join(sorted(allele))

    def allele_sort(self, df: DataFrame) -> DataFrame:
        snp_ref = df.iloc[:, :4]
        sorted_allele = df.iloc[:, 4:].applymap(self.__allele_sort, na_action='ignore')
        df = snp_ref.join(sorted_allele)
        return df

    def to_dataframe(self,
                     sheet: Union[str, int, List[Union[str, int]]] = 0,
                     index_col: int = None,
                     sort_allele: bool = True) -> DataFrame:
        try:  # read from text file
            if 'gz' in self.name:
                df = read_table(self.name, index_col=index_col)
            else:
                df = read_table(self.__open, index_col=index_col)
            if sort_allele:
                df = self.allele_sort(df)
            return df
        except UnicodeDecodeError:  # read from Excel file
            df = read_excel(self.name, sheet, index_col=index_col)
            if sort_allele:
                df = self.allele_sort(df)
            return df

    def to_dataframes(self,
                      sheet: Union[str, int, List[Union[str, int]]] = 0,
                      index_col: int = None,
                      chunk_size: int = 10000) -> DataFrame:
        if 'stdin' in self.name:  # read from stdin
            dfs = read_file_as_dataframe_from_stdin(chunk_size=chunk_size, index_col=index_col)
        elif self.name.endswith('gz'):  # read from compressed text (xx.gz) file
            dfs = read_table(self.name, chunksize=chunk_size, index_col=index_col)
        else:
            try:  # read from text file
                dfs = read_table(self.__open, chunksize=chunk_size, index_col=index_col)
            except UnicodeDecodeError:  # read from Excel file
                echo('\033[33mWarning: Excel file cannot be processed in chunks.\033[0m')
                return self.to_dataframe(sheet, index_col)
        return dfs

    def parallel_stat_MHM(self, processing_num: int):
        """Calculate the MissRate, HetRate and MAF (MHM) of SNP sites from GT files parallely."""
        dfs = self.to_dataframes()
        if not isinstance(dfs, DataFrame):  # Calculate with multiprocessing
            params = ((df,) for df in dfs)  # Index of DataFrame is 0, 1, 2, ...
            tkm = TaskManager(processing_num=processing_num, params=params)
            ret = tkm.parallel_run_func(stat_MHM)
            stat_dfs = [i.get() for i in ret]
            # Merge results of each multiprocessing
            stat_df = concat(stat_dfs)
        else:
            stat_df = stat_MHM(dfs)
        stat_df.sort_index(inplace=True, key=natsort_key)
        return stat_df

    def self_compare(self, other,
                     sheet1: Union[str, int, List[Union[str, int]]] = None,
                     sheet2: Union[str, int, List[Union[str, int]]] = None,
                     output_path: str = './') -> None:
        """Genotype consistency of different test batches in a single sample."""
        df1 = self.to_dataframe(sheet1, index_col=0)
        df2 = other.to_dataframe(sheet2, index_col=0)
        # Fill NA.
        df1.fillna('', inplace=True)
        df2.fillna('', inplace=True)
        # Sort by SNP ID (index of DataFrame).
        df1.sort_index(key=natsort_key, inplace=True)
        df2.sort_index(key=natsort_key, inplace=True)
        # Check whether the two GT files contain the same loci.
        if df1.index.tolist() != df2.index.tolist():
            echo('\033[31mError: The two GT file loci to be compared are inconsistent.\033[0m', err=True)
            exit()
        total_loci_num = len(df1)
        # Calculate genotype consistency of each sample under different test batches.
        bins = range(0, 110, 10)  # Set the consistency statistics interval.
        sample_count = df1.iloc[:, 3:].eq(df2.iloc[:, 3:]).sum(axis=0) / total_loci_num * 100
        bins_count = cut(sample_count, bins).value_counts(sort=False).to_string()
        sample_count = sample_count.to_string()
        # Write results to output file.
        echo(sample_count, open(f'{output_path}/Sample.consistency.xls', 'w'))
        echo(bins_count, open(f'{output_path}/Interval.stat.xls', 'w'))

    def compare(self, other,
                sheet1: Union[str, int, List[Union[str, int]]] = None,
                sheet2: Union[str, int, List[Union[str, int]]] = None,
                output_path: str = './') -> None:
        """Calculate genotype consistency."""
        # Step1: Read GT file as DataFrame.
        df1 = self.to_dataframe(sheet1)  # index = 0, 1, 2, ...
        df2 = other.to_dataframe(sheet2)  # index = 0, 1, 2, ...
        # Step2: Select the site intersection of two GT files.
        left_on = df1.columns[0]
        right_on = df2.columns[0]
        merge = df1.merge(df2, left_on=left_on,  # index = 0, 1, 2, ...
                          right_on=right_on)  # Avoid inconsistency between the two GT file ID fields
        merge.fillna('', inplace=True)  # fill NA
        # Step3: Calculate genotype consistency.
        loci_num = len(merge)
        df1_sample_num = len(df1.columns.tolist()) - 4
        left_sample_range = list(range(4, 4 + df1_sample_num))
        right_sample_range = list(range(len(df1.columns.tolist()) + 3, len(merge.columns.tolist())))
        with open(f'{output_path}/TestSample.consistency.xls', 'w') as o:
            echo('Database_sample\tTest_sample\tConsensus_number\tTotal_number\tRatio(%)', o)
            with tqdm(total=len(left_sample_range), unit='sample') as pbar:  # Show process bar
                for index1 in left_sample_range:
                    gt1 = merge.iloc[:, index1]
                    gt1.name = gt1.name.replace('_x', '')
                    pbar.set_description(f'Processing {gt1.name}')
                    for index2 in right_sample_range:
                        gt2 = merge.iloc[:, index2]
                        gt2.name = gt2.name.replace('_y', '')
                        consistency_count = gt1.eq(gt2).sum()
                        ratio = '%.2f' % (consistency_count / loci_num * 100)
                        echo(f'{gt1.name}\t{gt2.name}\t{consistency_count}\t{loci_num}\t{ratio}', o)
                    pbar.update(1)
        # Step4: Output GT file of test sample.
        right_sample_range.insert(0, 0)  # Only output site ID
        test_sample_gt_df = merge.iloc[:, right_sample_range]
        test_sample_gt_df.rename(columns=lambda i: str(i).replace('_y', '').replace('_x', ''), inplace=True)
        test_sample_gt_df.to_csv(f'{output_path}/TestSample.GT.xls', sep='\t', index=False)
