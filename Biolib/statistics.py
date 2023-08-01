"""
File: statistics.py
Description: 
Date: 2022/1/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
import click
import pandas as pd
from typing import Union


def display_set(decimal: int = 2) -> None:
    from warnings import filterwarnings
    filterwarnings("ignore")
    pd.set_option('display.float_format', lambda x: f'%.{decimal}f' % x)
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)


def read_in_gene_expression_as_dataframe(gene_exp_file: str) -> Union[str, pd.DataFrame]:
    try:
        function_dict = {'txt': pd.read_table, 'xlsx': pd.read_excel, 'xls': pd.read_excel, 'csv': pd.read_csv}
        df = function_dict[gene_exp_file.split('.')[-1]](gene_exp_file, index_col=0)
    except KeyError:
        return '\033[31mError: Unrecognised file formats. Only txt, xls, xlsx, and csv formats are supported.\033[0m'
    else:
        return df


def merge_duplicate_indexes(data_frame: pd.DataFrame) -> pd.DataFrame:
    """Sums the column values of rows with the same index"""
    index_name = data_frame.index.name
    data_frame = data_frame.groupby(index_name).sum()
    return data_frame


def filter_by_min_value(data_frame: pd.DataFrame, min_value: float, start_column_num: int = 1,
                        end_column_num: int = None) -> pd.DataFrame:
    """Delete rows where all column values are less than the specified value"""
    if end_column_num is None:
        data_frame = data_frame[~(data_frame.iloc[:, start_column_num - 1:] <= min_value).all(axis=1)]
        return data_frame
    else:
        data_frame = data_frame[~(data_frame.iloc[:, start_column_num - 1:end_column_num - 1] <= min_value).all(axis=1)]
        return data_frame


def get_FPKM(read_count_DataFrame: pd.DataFrame,
             min_value: float = None) -> pd.DataFrame:
    """
    Standardize gene expression with FPKM.
    :param read_count_DataFrame: read count DataFrame. (type=pd.DataFrame)

            Geneid     Length   Sample1   Sample2   Sample3
            gene1      200         1         0         0
            gene2      1000       21        51        34
            ......     .......   .......   .......   .......
            gene1000   2100      2345      2137      1987

    :param min_value: Gene minimum expression (genes whose expression is less than the specified value in all samples
                      are filtered out).
    :return: FPKM matrix. (type=pd.DataFrame)
    """
    read_count_DataFrame = merge_duplicate_indexes(read_count_DataFrame)
    div_gene_length = read_count_DataFrame.div(read_count_DataFrame.iloc[:, 0], axis=0)
    div_gene_length.loc['sum'] = read_count_DataFrame.sum(axis=0)
    div_read_num = div_gene_length.div(div_gene_length.loc['sum'].T, axis=1) * 10 ** 9
    FPKM = div_read_num.iloc[0:-1, 1:]
    raw_gene_num = len(FPKM.index.tolist())
    if min_value:
        FPKM = filter_by_min_value(FPKM, min_value=min_value)
        new_gene_num = len(FPKM.index.tolist())
        click.echo(f"{raw_gene_num - new_gene_num} genes have been filtered out", err=True)
    return FPKM


def get_TPM(read_count_DataFrame: pd.DataFrame,
            min_value: float = None) -> pd.DataFrame:
    """
    Standardize gene expression with TPM.
    :param read_count_DataFrame: read count DataFrame. (type=pd.DataFrame)

            Geneid     Length   Sample1   Sample2   Sample3
            gene1      200         1         0         0
            gene2      1000       21        51        34
            ......     .......   .......   .......   .......
            gene1000   2100      2345      2137      1987

    :param min_value: Gene minimum expression (genes whose expression is less than the specified value in all samples
                      are filtered out).
    :return: TPM matrix. (type=pd.DataFrame)
    """
    read_count_DataFrame = merge_duplicate_indexes(read_count_DataFrame)
    read_count_DataFrame = read_count_DataFrame.div(read_count_DataFrame.iloc[:, 0], axis=0)
    read_count_DataFrame.loc['sum'] = read_count_DataFrame.sum(axis=0)
    read_count_DataFrame = read_count_DataFrame.div(read_count_DataFrame.iloc[-1], axis=1)
    read_count_DataFrame = read_count_DataFrame * 10 ** 6
    TPM = read_count_DataFrame.iloc[:-1, 1:]
    raw_gene_num = len(TPM.index.tolist())

    # FPKM = get_FPKM(read_count_DataFrame)
    # FPKM.loc['sum'] = FPKM.sum(axis=0)
    # TPM = FPKM.div(FPKM.iloc[-1], axis=1).iloc[:-1] * 10 ** 6
    # raw_gene_num = len(TPM.index.tolist())

    if min_value:
        TPM = filter_by_min_value(TPM, min_value=min_value)
        new_gene_num = len(TPM.index.tolist())
        click.echo(f"{raw_gene_num - new_gene_num} genes have been filtered out", err=True)
    return TPM
