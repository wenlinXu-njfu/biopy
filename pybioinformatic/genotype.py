from typing import Union, List
from io import TextIOWrapper
from os.path import abspath
from warnings import filterwarnings
from natsort import natsort_key
from pandas import Series, DataFrame, read_table, read_excel, concat
from pybioinformatic.task_manager import TaskManager
from pybioinformatic.statistics import read_file_as_dataframe_from_stdin
filterwarnings("ignore")


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


def __allele_count(row: Series) -> str:
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
    snp_ref = df.copy().iloc[:, :4]
    # Calculate MissRate
    sample_num = len(df.columns.tolist()) - 4
    df['MissRate'] = df.isnull().sum(axis=1) / sample_num
    # Calculate HetRate
    df.fillna('', inplace=True)  # fill NA
    df['HetRate'] = df.iloc[:, 4:-1].apply(lambda row: 1 - __check_hom(row).sum() / (row != '').sum(), axis=1)
    # Calculate MAF
    df['all_gt'] = df.iloc[:, 4:-2].apply(lambda row: __allele_count(row), axis=1)
    df['total'] = df['all_gt'].apply(len)
    df['A'] = df['all_gt'].str.count('A') / df['total']
    df['G'] = df['all_gt'].str.count('G') / df['total']
    df['C'] = df['all_gt'].str.count('C') / df['total']
    df['T'] = df['all_gt'].str.count('T') / df['total']
    df['D'] = df['all_gt'].str.count('D') / df['total']
    df['I'] = df['all_gt'].str.count('I') / df['total']
    df['MAF'] = df.loc[:, ['A', 'G', 'C', 'T', 'D', 'I']].apply(lambda row: sorted(row, reverse=True)[1], axis=1)
    # Output results
    df = df.loc[:, ['MissRate', 'HetRate', 'MAF']]
    merge = snp_ref.join(df)
    return merge


class GenoType:
    """
    Standard genotype file object.
    """

    def __init__(self, path: Union[str, TextIOWrapper]):
        self.name = abspath(path) if isinstance(path, str) else abspath(path.name.replace('<', '').replace('>', ''))
        self.__open = path

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
                     sort_allele: bool = True) -> DataFrame:
        try:  # read from text file
            df = read_table(self.__open)
            if sort_allele:
                df = self.allele_sort(df)
            return df
        except UnicodeDecodeError:  # read from Excel file
            df = read_excel(self.__open, sheet)
            if sort_allele:
                df = self.allele_sort(df)
            return df

    def to_dataframes(self,
                      sheet: Union[str, int, List[Union[str, int]]] = None,
                      chunk_size: int = 10000) -> DataFrame:
        if 'stdin' in self.name:  # read from stdin
            dfs = read_file_as_dataframe_from_stdin(chunk_size=chunk_size, index_col=None)
        else:
            try:  # read from text file
                dfs = read_table(self.__open, chunksize=chunk_size)
            except UnicodeDecodeError:  # read from Excel file
                dfs = read_excel(self.__open, sheet)
        yield from dfs

    def parallel_stat_MHM(self, processing_num: int):
        """Calculate the MissRate, HetRate and MAF (MHM) of SNP sites from GT files parallely."""
        # Calculate with multiprocessing
        params = ((df,) for df in self.to_dataframes())
        tkm = TaskManager(processing_num=processing_num, params=params)
        ret = tkm.parallel_run_func(stat_MHM)
        stat_dfs = [i.get() for i in ret]
        # Merge results of each multiprocessing
        stat_df = concat(stat_dfs)
        stat_df.sort_index(inplace=True, key=natsort_key)
        return stat_df

    def compare(self, other,
                sheet1: Union[str, int, List[Union[str, int]]] = None,
                sheet2: Union[str, int, List[Union[str, int]]] = None,
                output_prefix: str = 'TestSample') -> str:
        """Calculate genotype consistency."""
        df1 = self.to_dataframe(sheet1)
        df2 = other.to_dataframe(sheet2)
        df1_sample_num = len(df1.columns.tolist()) - 4
        # 取出两个GT文件的位点交集
        left_on = df1.columns[0]
        right_on = df2.columns[0]
        merge = df1.merge(df2, left_on=left_on, right_on=right_on)
        merge.fillna('', inplace=True)
        loci_num = len(merge)
        left_sample_range = list(range(4, 4 + df1_sample_num))
        right_sample_range = list(range(len(df1.columns.tolist()) + 3, len(merge.columns.tolist())))
        # 计算基因型一致率
        yield 'Sample1\tSample2\tConsensus_number\tTotal_number\tRatio(%)'
        for index1 in left_sample_range:
            gt1 = merge.iloc[:, index1]
            gt1.name = gt1.name.replace('_x', '')
            for index2 in right_sample_range:
                gt2 = merge.iloc[:, index2]
                gt2.name = gt2.name.replace('_y', '')
                consistency_count = gt1.eq(gt2).sum()
                ratio = '%.2f' % (consistency_count / loci_num * 100)
                yield f'{gt1.name}\t{gt2.name}\t{consistency_count}\t{loci_num}\t{ratio}'
        # 输出测试样本GT文件
        right_sample_range.insert(0, 0)  # 只输出位点ID，不输出位点所在染色体及位置
        merge.iloc[:, right_sample_range].to_csv(f'{output_prefix}.GT.xls', sep='\t', index=False)
