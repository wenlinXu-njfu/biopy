"""
File: gff.py
Description: Instantiate a GFF file object.
CreateDate: 2021/11/27
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper, StringIO
from typing import Union, List, Dict, Tuple, Generator
from os.path import abspath
from re import findall
from gzip import GzipFile
from click import echo
from pandas import DataFrame, read_table
from pybioinformatic.fasta import Fasta
from pybioinformatic.sequence import Nucleotide


class Gff:
    def __init__(self, path: Union[str, TextIOWrapper]):
        if isinstance(path, str):
            self.name = abspath(path)
            if path.endswith('gz'):
                self.__open = GzipFile(self.name, 'rb')
                self.line_num = sum(1 for line in self.__open if not str(line, 'utf8').startswith('#'))
                self.__open.seek(0)
                self.anno_line_num = sum(1 for line in self.__open if str(line, 'utf8').startswith('#'))
                self.__open.seek(0)
            else:
                self.__open = open(path)
                self.line_num = sum(1 for line in self.__open if not line.startswith('#'))
                self.__open.seek(0)
                self.anno_line_num = sum(1 for line in self.__open if line.startswith('#'))
                self.__open.seek(0)
        else:
            self.name = abspath(path.name)
            if path.name == '<stdin>':
                string_io = StringIO()
                for line in path:
                    string_io.write(line)
                string_io.seek(0)
                self.__open = string_io
                self.line_num = sum(1 for line in self.__open if not line.startswith('#'))
                self.__open.seek(0)
                self.anno_line_num = sum(1 for line in self.__open if line.startswith('#'))
                self.__open.seek(0)
            else:
                if self.name.endswith('gz'):
                    self.__open = GzipFile(self.name, 'rb')
                    self.line_num = sum(1 for line in self.__open if not str(line, 'utf8').startswith('#'))
                    self.__open.seek(0)
                    self.anno_line_num = sum(1 for line in self.__open if str(line, 'utf8').startswith('#'))
                    self.__open.seek(0)
                else:
                    self.__open = path
                    self.line_num = sum(1 for line in self.__open if not line.startswith('#'))
                    self.__open.seek(0)
                    self.anno_line_num = sum(1 for line in self.__open if line.startswith('#'))
                    self.__open.seek(0)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.__open.close()
        except AttributeError:
            pass

    def __seek_zero(self):
        try:
            self.__open.seek(0)
        except AttributeError:
            pass

# Basic method==========================================================================================================
    def parse(self) -> Generator[Tuple[str, str, str, str, str, str, str, str, Dict[str, str]], None, None]:
        """Parse information of each column of GFF file line by line."""
        for line in self.__open:
            line = str(line, 'utf8') if isinstance(line, bytes) else line
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_num, source, feature = split[0], split[1], split[2]
                start, end, score, strand, frame = split[3], split[4], split[5], split[6], split[7]
                attr_list = [attr for attr in split[8].split(';') if '=' in attr]
                attr_dict: Dict[str, str] = {attr.split('=')[0]: attr.split('=')[1] for attr in attr_list if attr}
                yield chr_num, source, feature, start, end, score, strand, frame, attr_dict
        self.__seek_zero()

    def to_dataframe(self) -> DataFrame:
        names = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute']
        if self.name.endswith('gz'):
            df = read_table(self.name, header=None, names=names, dtype=str, comment='#')
        else:
            df = read_table(self.__open, header=None, names=names, dtype=str, comment='#')
        df[['Start', 'End']] = df[['Start', 'End']].astype(int)
        self.__seek_zero()
        return df

    def __check_feature(self, feature: str) -> Tuple[bool, str]:
        """Check whether specified feature is included in GFF file."""
        features = set(line[2] for line in self.parse())
        if feature in features:
            return True, f'"{feature}" have found.'
        else:
            return False, f'"{feature}" not found.'

    def to_dict(self, feature_type: str = None) -> Dict[str, List[Dict[str, Union[str, int]]]]:
        """Save the feature information in the GFF file into the dictionary."""
        # gff_dict = {
        #             Chr_num: [{id: str, start: int, end: int, strand: str}, {}, ...],
        #             Chr_num: [{}, {}, ...], ...
        #             }
        is_in_gff, msg = self.__check_feature(feature_type)
        if not is_in_gff:
            echo(f'\033[31mError: {msg}\033[0m', err=True)
            exit()
        gff_dict = {}
        for line in self.parse():
            if line[2] == feature_type or feature_type is None:
                item = {'id': line[8]['ID'], 'start': int(line[3]), 'end': int(line[4]), 'strand': line[6]}
                if line[0] in gff_dict:
                    gff_dict[line[0]].append(item)
                else:
                    gff_dict[line[0]] = [item]
        return gff_dict

    def get_mRNA_dict(self) -> Dict[str, List[Dict[str, Union[str, int]]]]:
        """Get mRNA dict. The start and end based on mRNA length, not based on chromosome length."""
        mRNA_dict = {}  # {mRNA_id: [{feature_type: str, start: int, end: int, strand: str}, {}, ...], ...}
        mRNA_id = None
        mRNA_start = 0
        mRNA_len = 0
        for line in self.parse():
            if line[2] == 'mRNA':
                mRNA_id = line[8]['ID']
                mRNA_start = int(line[3])
                mRNA_len = int(line[4]) - int(line[3]) + 1
                if mRNA_id not in mRNA_dict:
                    mRNA_dict[mRNA_id] = []
                else:
                    echo(f'\033[31mError: The GFF file has repeat id {mRNA_id}.\033[0m')
                    exit()
            elif line[2] == 'CDS' or 'UTR' in line[2]:
                if line[8]['Parent'] == mRNA_id:
                    if line[6] == '-':
                        start = mRNA_len - (int(line[4]) - mRNA_start) - 1
                        end = mRNA_len - (int(line[3]) - mRNA_start) - 1
                        item = {'feature_type': line[2], 'start': start, 'end': end, 'strand': line[6]}
                    else:
                        item = {'feature_type': line[2], 'start': int(line[3]) - mRNA_start,
                                'end': int(line[4]) - mRNA_start, 'strand': line[6]}
                    mRNA_dict[mRNA_id].append(item)
                else:
                    echo('\033[31mError: GFF file is not sorted by mRNA ID.\033[0m', err=True)
                    exit()
        return mRNA_dict

    def summary(self) -> str:
        """A peek at the genome."""
        df = self.to_dataframe()
        df['Length'] = df['End'] - df['Start'] + 1
        content = [f'# Feature\tTotal\tMin_len\tMax_len\tMedian_len\tMean_len\n']
        features = df['Feature'].drop_duplicates().values
        for feature in features:
            total = len(df.loc[df['Feature'] == feature])
            min_len = df.loc[df['Feature'] == feature, 'Length'].min()
            max_len = df.loc[df['Feature'] == feature, 'Length'].max()
            median_len = '%.0f' % df.loc[df['Feature'] == feature, 'Length'].median()
            mean_len = '%.0f' % df.loc[df['Feature'] == feature, 'Length'].mean()
            content.append(f'{feature}\t{total}\t{min_len}\t{max_len}\t{median_len}\t{mean_len}\n')
        content = ''.join(content)
        return content

# GFF file sorted by id method==========================================================================================
    @staticmethod
    def __gff_sort(line: str) -> tuple:
        if not line.startswith('#') and line.strip():
            split = line.strip().split('\t')
            alpha = findall(r'[a-zA-Z]+', split[0])[0]
            number = int(findall(r'\d+', split[0])[0])
            if alpha == 'Chr' or alpha == 'chr':
                alpha = 0
            else:
                alpha = 1
            attr = split[8].replace('=', ';').split(';')
            if split[2] == 'gene':
                locus_id = findall(r'\w+.\w+', attr[attr.index('ID') + 1])[0]
                parent_id = '0'
            elif split[2] == 'mRNA':
                locus_id = attr[attr.index('ID') + 1]
                parent_id = '1'
            else:
                parent_id = attr[attr.index('Parent') + 1]
                locus_id = parent_id
            start, end = int(split[3]), int(split[4])
            return alpha, number, locus_id, parent_id, start, -end
        else:
            return 0, 0, '0', '0', 0, line

    def sort(self) -> Generator[str, None, None]:
        """Sort the GFF file by sequence ID."""
        if self.name.endswith('gz'):
            l = [str(line, 'utf8').strip() for line in self.__open
                 if not str(line, 'utf8').startswith('#') and line.strip()]
        else:
            l = [line.strip() for line in self.__open if not line.startswith('#') and line.strip()]
        l.sort(key=lambda line: self.__gff_sort(line))
        if not isinstance(self.__open, list):
            self.__open.seek(0)
        yield from l

# Sequence extraction method============================================================================================
    def extract_seq(self,
                    fasta_file: Union[str, TextIOWrapper],
                    feature_type: str = 'gene',
                    feature_id_set: set = None) -> Nucleotide:
        """Extract sequences of specified feature type from GFF file."""
        is_in_gff, msg = self.__check_feature(feature_type)
        if not is_in_gff:
            echo(f'\033[31mError: {msg}\033[0m', err=True)
            exit()
        gff_dict = self.to_dict(feature_type)
        with Fasta(fasta_file) as fa:
            for nucl_obj in fa.parse():
                try:
                    features = gff_dict[nucl_obj.id]  # features = [{feature1}, {feature2}, ...]
                except KeyError:
                    pass  # Some sequences (eg. scaffold, contig) may not have annotation
                else:
                    for feature in features:  # feature = {id: str, start: int, end: int, strand: str}
                        if feature_id_set and feature['id'] in feature_id_set:
                            sub_seq_obj = nucl_obj[feature['start'] - 1:feature['end']]
                            sub_seq_obj.id = feature['id']
                            yield sub_seq_obj
                        elif not feature_id_set:
                            sub_seq_obj = nucl_obj[feature['start'] - 1:feature['end']]
                            sub_seq_obj.id = feature['id']
                            yield sub_seq_obj

    def miRNA_extraction(self) -> Nucleotide:
        """Extract miRNA sequence from GFF file."""
        for line in self.parse():
            attr_dict = line[8]
            seq_id = attr_dict['ID']
            seq = attr_dict['seq']
            yield Nucleotide(seq_id, seq)

# File format conversion method=========================================================================================
    def to_gtf(self) -> Generator[str, None, None]:
        """Convert the file format from GFF to GTF."""
        last_line = None
        gene_id = transcript_id = None
        i = 0
        for line in self.parse():
            i += 1
            current_line = list(line)
            if current_line[2] == 'gene':
                if last_line:
                    yield '\t'.join(last_line[:8]) + f'''\tgene_id "{gene_id}"; transcript_id "{transcript_id}";'''
                    last_line = None
                    gene_id = current_line[8]['ID']
                else:
                    gene_id = current_line[8]['ID']
                yield '\t'.join(current_line[:8]) + f'''\tgene_id "{gene_id}";'''
            elif current_line[2] == 'mRNA' or current_line[2] == 'transcript':
                current_line[2] = 'transcript'
                if last_line:
                    if gene_id is not None:
                        yield '\t'.join(last_line[:8]) + f'''\tgene_id "{gene_id}"; transcript_id "{transcript_id}";'''
                    else:
                        yield '\t'.join(last_line[:8]) + f'''\ttranscript_id "{transcript_id}";'''
                    last_line = None
                    transcript_id = current_line[8]['ID']
                else:
                    transcript_id = current_line[8]['ID']
                if gene_id is not None:
                    yield '\t'.join(current_line[:8]) + f'''\tgene_id "{gene_id}"; transcript_id "{transcript_id}";'''
                else:
                    yield '\t'.join(current_line[:8]) + f'''\ttranscript_id "{transcript_id}";'''
            elif 'UTR' in current_line[2] or 'CDS' in current_line[2]:
                current_line[2] = 'exon'
                current_line[7] = '.'
                if last_line:
                    if last_line[8]['Parent'] == current_line[8]['Parent'] == transcript_id:
                        if int(last_line[4]) + 1 == int(current_line[3]):
                            current_line[3] = last_line[3]
                            last_line = current_line
                        else:
                            if gene_id is not None:
                                yield ('\t'.join(last_line[:8]) +
                                       f'''\tgene_id "{gene_id}"; transcript_id "{transcript_id}";''')
                            else:
                                yield '\t'.join(last_line[:8]) + f'''\ttranscript_id "{transcript_id}";'''
                            last_line = current_line
                else:
                    last_line = current_line
            if self.line_num == i:
                if gene_id is not None:
                    yield '\t'.join(last_line[:8]) + f'''\tgene_id "{gene_id}"; transcript_id "{transcript_id}";'''
                else:
                    yield '\t'.join(last_line[:8]) + f'''\ttranscript_id "{transcript_id}";'''

    def to_bed(self, feature_type: Union[str, list] = None) -> str:
        """Convert the file format from GFF to BED."""
        for line in self.parse():
            if feature_type:
                if line[2] == feature_type or line[2] in feature_type:
                    yield f"{line[0]}\t{int(line[3]) - 1}\t{line[4]}\t{line[8]['ID']}\t{line[7]}\t{line[6]}"
            elif not feature_type:
                yield f"{line[0]}\t{int(line[3]) - 1}\t{line[4]}\t{line[8]['ID']}\t{line[7]}\t{line[6]}"

    def to_gsds(self) -> str:
        """Convert the file format from GFF to GSDS."""
        transcript_id, transcript_start = None, 0
        for line in self.parse():
            line = list(line)
            if line[2] == 'mRNA':
                transcript_id = line[8]['ID']
                transcript_start = int(line[3])
            elif line[2] == 'CDS' or 'UTR' in line[2]:
                if line[8]['Parent'] == transcript_id:
                    line[3] = str(int(line[3]) - transcript_start)
                    line[4] = (int(line[4]) - transcript_start)
                    yield f"{transcript_id}\t{line[3]}\t{line[4]}\t{line[2]}\t{line[7]}"
                else:
                    echo(f'\033[33mWarning: The order of GFF file is wrong, '
                               f'this will cause some information to be lost.\033[0m', err=True)

# Feature density count=================================================================================================
    def get_feature_density(self,
                            chr_len_dict: Dict[str, int],
                            feature_type: str = 'gene',
                            span: int = 100000) -> Generator[str, None, None]:
        """Get feature density."""
        is_in_gff, msg = self.__check_feature(feature_type)
        if not is_in_gff:
            echo(f'\033[31mError: {msg}\033[0m', err=True)
            exit()
        if min(list(chr_len_dict.values())) / span < 1:
            echo('\033[33mError: Density statistical interval is too large.\033[0m', err=True)
            exit()
        gff_dataframe = self.to_dataframe()
        gff_dataframe['Site'] = (gff_dataframe['Start'] + gff_dataframe['End']) / 2
        for chr_num, length in chr_len_dict.items():
            df = gff_dataframe[(gff_dataframe['Chromosome'] == chr_num) & (gff_dataframe['Feature'] == feature_type)]
            sites = df['Site']
            for i in range(span, length, span):
                count = len(sites[(sites <= i) & (sites >= i - span + 1)])
                yield f'{chr_num}\t{i - span + 1}\t{i}\t{count}'
            if length % span != 0:
                count = len(sites[(sites <= length) & (sites >= length // span * span + 1)])
                yield f'{chr_num}\t{length // span * span + 1}\t{length}\t{count}'
