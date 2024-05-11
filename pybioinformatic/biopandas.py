"""
File: biopandas.py
Description: Common functions to handle dataframe.
CreateDate: 2022/1/10
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from typing import Union, List, Iterable
from io import StringIO
from warnings import filterwarnings
from re import sub
from click import echo, open_file
from pandas import cut, set_option, read_table, read_excel, read_csv, Series, DataFrame, ExcelWriter


def display_set(decimal: int = 2) -> None:
    filterwarnings("ignore")
    set_option('display.float_format', lambda x: f'%.{decimal}f' % x)
    set_option('display.width', 1000)
    set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)


def read_file_as_dataframe_from_stdin(sep: str = '\t',
                                      line_terminator: str = '\n',
                                      index_col: Union[int, None] = 0,
                                      header: int = 0,
                                      skip_rows: int = 0,
                                      use_cols: List[Union[str, int]] = None,
                                      names: List[str] = None,
                                      chunk_size: int = None):
    df = read_table(StringIO(''.join(open_file('-').readlines())),
                    sep=sep, header=header, names=names, index_col=index_col, usecols=use_cols,
                    lineterminator=line_terminator, skiprows=skip_rows, chunksize=chunk_size)
    return df


def read_in_gene_expression_as_dataframe(gene_exp_file: str) -> Union[str, DataFrame]:
    function_dict = {'txt': read_table, 'xlsx': read_excel, 'xls': read_excel, 'csv': read_csv}
    try:
        df = function_dict[gene_exp_file.split('.')[-1]](gene_exp_file, index_col=0)
    except KeyError:
        return '\033[31mError: Unrecognised file formats. Only txt, xls, xlsx, and csv formats are supported.\033[0m'
    else:
        return df


def merge_duplicate_indexes(df: DataFrame) -> DataFrame:
    """Sums the column values of rows with the same index"""
    index_name = df.index.name
    data_frame = df.groupby(index_name).sum()
    return data_frame


def filter_by_min_value(df: DataFrame,
                        min_value: float,
                        start_column_num: int = 1,
                        end_column_num: int = None) -> DataFrame:
    """Delete rows where all column values are less than the specified value"""
    if end_column_num is None:
        data_frame = df[~(df.iloc[:, start_column_num - 1:] <= min_value).all(axis=1)]
        return data_frame
    else:
        data_frame = df[~(df.iloc[:, start_column_num - 1:end_column_num - 1] <= min_value).all(axis=1)]
        return data_frame


def get_FPKM(read_count: DataFrame,
             min_value: float = None) -> DataFrame:
    """
    Standardize gene expression with FPKM.
    :param read_count: read count DataFrame. (type=pd.DataFrame)

            Geneid     Length   Sample1   Sample2   Sample3
            gene1      200         1         0         0
            gene2      1000       21        51        34
            ......     .......   .......   .......   .......
            gene1000   2100      2345      2137      1987

    :param min_value: Gene minimum expression (genes whose expression is less than the specified value in all samples
                      are filtered out).
    :return: FPKM matrix. (type=pd.DataFrame)
    """
    read_count_DataFrame = merge_duplicate_indexes(read_count)
    div_gene_length = read_count_DataFrame.div(read_count_DataFrame.iloc[:, 0], axis=0)
    div_gene_length.loc['sum'] = read_count_DataFrame.sum(axis=0)
    div_read_num = div_gene_length.div(div_gene_length.loc['sum'].T, axis=1) * 10 ** 9
    FPKM = div_read_num.iloc[0:-1, 1:]
    raw_gene_num = len(FPKM.index.tolist())
    if min_value:
        FPKM = filter_by_min_value(FPKM, min_value=min_value)
        new_gene_num = len(FPKM.index.tolist())
        echo(f"{raw_gene_num - new_gene_num} genes have been filtered out", err=True)
    return FPKM


def get_TPM(read_count: DataFrame,
            min_value: float = None) -> DataFrame:
    """
    Standardize gene expression with TPM.
    :param read_count: read count DataFrame. (type=pd.DataFrame)

            Geneid     Length   Sample1   Sample2   Sample3
            gene1      200         1         0         0
            gene2      1000       21        51        34
            ......     .......   .......   .......   .......
            gene1000   2100      2345      2137      1987

    :param min_value: Gene minimum expression (genes whose expression is less than the specified value in all samples
                      are filtered out).
    :return: TPM matrix. (type=pd.DataFrame)
    """
    read_count_DataFrame = merge_duplicate_indexes(read_count)
    read_count_DataFrame = read_count_DataFrame.div(read_count_DataFrame.iloc[:, 0], axis=0)
    read_count_DataFrame.loc['sum'] = read_count_DataFrame.sum(axis=0)
    read_count_DataFrame = read_count_DataFrame.div(read_count_DataFrame.iloc[-1], axis=1)
    read_count_DataFrame = read_count_DataFrame * 10 ** 6
    TPM = read_count_DataFrame.iloc[:-1, 1:]
    raw_gene_num = len(TPM.index.tolist())

    if min_value:
        TPM = filter_by_min_value(TPM, min_value=min_value)
        new_gene_num = len(TPM.index.tolist())
        echo(f"{raw_gene_num - new_gene_num} genes have been filtered out", err=True)
    return TPM


def dfs_to_excel(df_list: list, output_file: str, sheets_name: list = None):
    with ExcelWriter(output_file, engine='xlsxwriter') as writer:
        if not sheets_name:
            sheets_name = [f'Sheet {i}' for i in range(1, len(df_list) + 1)]
        workbook = writer.book
        workbook.nan_inf_to_errors = True  # 添加 'nan_inf_to_errors' 选项
        default_cell_format = workbook.add_format({'bold': False, 'border': 0, 'align': 'left'})
        for df, sheet_name in zip(df_list, sheets_name):
            df.fillna('NA', inplace=True)
            df = df.reset_index()
            df.to_excel(writer, sheet_name=sheet_name, na_rep='NA', startrow=1, header=False, index=False)
            worksheet = writer.sheets[sheet_name]  # 获取工作表对象
            # 将标题写入工作表，并应用默认的单元格格式
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, default_cell_format)
            # 将数据写入工作表，并应用默认的单元格格式
            for row_num, row_data in enumerate(df.values):
                for col_num, col_data in enumerate(row_data):
                    worksheet.write(row_num + 1, col_num, col_data, default_cell_format)


def dataframe_to_str(df: DataFrame,
                     index: bool = True,
                     header: bool = True):
    if index:
        df = df.reset_index()
    string_df = sub(r'\n +', '\n', df.to_string(index=False, header=header).strip())
    string_df = sub(r' +', '\t', string_df)
    return string_df


def interval_stat(ser: Series, bins: Iterable, precision: int = 3, name: str = 'count'):
    ret: DataFrame = cut(x=ser, bins=bins, precision=precision, include_lowest=True).value_counts(sort=False).to_frame(name)
    ret.index.name = ret.index.name + '_interval'
    ret.rename(index={ret.index[0]: sub(r'.+,', '[0,', str(ret.index[0]))}, inplace=True)
    ret.loc['total', name] = ret[name].sum()
    return ret
