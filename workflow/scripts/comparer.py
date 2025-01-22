import ast
import os
from dataclasses import dataclass
import os.path
from pathlib import Path
from typing import List, Optional, Tuple, Union, Any
import psutil
import pandas as pd
from Bio.Seq import Seq
from csvsort import csvsort
import typer
import csv
from itertools import combinations

from pandas import Series, DataFrame

from count_fastq import read_csv, count_reads


class CSV:
    # used for filtering when extracting specific read types
    all_columns = ['flags', 'pos', 'chr', 'CIGAR',
                   'quality', 'edit_dist', 'MD', 'type', 'proper_pair']

    def __init__(self, path_to_csv: str, is_original: bool = True):
        self.path = path_to_csv
        self.is_single_ended = True
        self.is_real_data = False
        self.is_opened = False
        self.label = self.get_label(self.path)
        with open(self.path, 'r') as csv_file:
            reader = csv.reader(csv_file, delimiter=',')
            self.header = next(reader)
            if self.header[-1][-8:] == 'sequence':
                self.is_real_data = True
        if self.is_real_data:
            if is_original:
                self.sub = self.create_sub_csv()
            else:
                self.sub = self.path[:-7] + ".sub" + self.path[-4:]
            # buffer size (in number of reads) for I/O operations
            self.sub_BUFFER_SIZE = 50_000
            self.sub_buffered_reads = []

    @staticmethod
    def get_label(path_to_csv):
        with open(path_to_csv, 'r') as csv_file:
            reader = csv.reader(csv_file, delimiter=',')
            header = next(reader)
            return header[0][:-5]

    def _open(self):
        self.opened = open(self.path, 'r')
        self.is_opened = True
        return self.opened

    def _close(self):
        self.opened.close()
        self.is_opened = False

    def format_output_path(self, suffix: str):
        head, tail = os.path.split(self.path)
        tails = tail.split('.')
        tails.insert(-1, suffix)
        tail = '.'.join(tails)
        output_path = os.path.join(head, tail)
        return output_path

    def create_sub_csv(self, clone_header=True):
        new_csv_path = self.format_output_path('sub')
        with open(self.path, 'r') as csv_master:
            with open(new_csv_path, 'w') as csv_sub:
                if clone_header:
                    reader = csv.reader(csv_master, delimiter=',')
                    writer = csv.writer(csv_sub, delimiter=',')
                    header = next(reader)
                    writer.writerow(header[:-1])
        return new_csv_path

    def create_rv_csv(self, column: int = 9):
        BUFFER_SIZE = 100_000
        new_csv_path = self.format_output_path(suffix='rv')
        with open(self.path, 'r') as csv_master:
            with open(new_csv_path, 'w') as csv_sub:
                reader = csv.reader(csv_master, delimiter=',')
                writer = csv.writer(csv_sub, delimiter=',')
                header = next(reader)
                writer.writerow(header)
                buffered_rows = []
                for row in reader:
                    reversed_seq = Seq(row[column]).reverse_complement()
                    new_row = row
                    new_row[-1] = reversed_seq
                    buffered_rows.append(new_row)
                    if len(buffered_rows) == BUFFER_SIZE:
                        writer.writerows(buffered_rows)
                        buffered_rows = []
                if len(buffered_rows) != 0:
                    writer.writerows(buffered_rows)
                    buffered_rows = []
        csvsort(new_csv_path, [column], parallel=False, max_size=1000)
        return new_csv_path

    def _append(self):
        """Appends rows to the SUB(!) CSV file"""
        with open(self.sub, 'a') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerows(self.sub_buffered_reads)
        self.sub_buffered_reads = []

    def append_row(self, row):
        """Appends rows to the SUB(!) CSV file"""
        self.sub_buffered_reads.append(row)
        if len(self.sub_buffered_reads) == self.sub_BUFFER_SIZE:
            self._append()

    def finalize_append(self):
        """Appends rows to the SUB(!) CSV file"""
        if len(self.sub_buffered_reads) != 0:
            self._append()

    def get_reader(self):
        stream = self._open()
        reader = csv.reader(stream, delimiter=',')
        return reader


@dataclass
class ComparisonObject:
    """Input Data processing class

    Arguments
    ---------
    data : dict
        Dictionary with data for every sample. Can be accessed through
        input_data dictionary, e.g. input_data['original'],
        input_data['reversed'], input_data['shuffled']

    Attributes
    ----------
    label : str
        Label of the sample.
    path_to_csv : str
        Path to CSV table which contains columns with the corresponding label.
    is_comparable : bool
        True if all the required data is given, False if left blank"""

    data: dict

    def __post_init__(self):
        self.label = self.data['label']
        self.path_to_csv = self.data['path_to_csv']
        compare_condition = self.label.strip() and self.path_to_csv.strip()
        self.is_comparable = bool(compare_condition)


class Comparer:
    """
    Arguments
    ---------
    input_data : dict
        Dictionary that contains labels and path to the corresponding CSV
        for original, reversed and shuffled data.
    output_path : str (OPTIONAL)
        Name of the output CSV file, if left blank creates
        'compare_output.csv' in the same folder by default.

    Attributes
    ----------
    original_data : ComparisonObject
        Instance of ComparisonObject. Label, path to CSV and whether it should
        be compared can be accessed through attributes original_data.label,
        original_data.path_to_csv, original_data.is_comparable respectively.
    reversed_data : ComparisonObject
        The same as original_data but contains info about reversed sample
    shuffled_data : ComparisonObject
        The same as original_data but contains info about shuffled sample
    samples : list[ComparisonObject]
        List of all three ComparisonObject instances
    is_single_ended : bool
        Checks whether given data is single ended or not
    comparable : list[bool]
        List of booleans just to check how many samples to compare
    """

    def __init__(self, input_data: dict,
                 fastq_reads: int,
                 output_path: str = '',
                 reads_to_extract: list = [],
                 filtered_columns: list = []
                 ) -> None:
        self.input_data = input_data
        self.fastq_reads = fastq_reads
        self.output_path = self._parse_output_path(output_path)
        self.reads_to_extract = reads_to_extract
        self.filtered_columns = filtered_columns
        with open(self.output_path, 'w') as output:
            headers = ['FEATURE', 'READS', 'PERCENTAGE']
            writer = csv.DictWriter(output, fieldnames=headers)
            writer.writeheader()
        self._check_input()

    def write_row(self, feature: str, reads: int, total_reads: int):
        """
        Appends new row to the csv output file

        Arguments
        ---------
        feature : str
            Name of the feature, i.e. 'Unambiguous_reads_original_sample'
        reads : int
            Number of filtered reads that met the requirements
        total_reads : int
            Number of unfiltered reads to calculate the percentage
            (almost always number of unambiguous reads)
        """
        if isinstance(reads, str) and reads == "NA":
            percentage = "NA"
        elif total_reads == 0:
            percentage = 0
        else:
            percentage = round(100 * reads / total_reads, 3)

        with open(self.output_path, 'a') as f:
            header = ['FEATURE', 'READS', 'PERCENTAGE']
            row = {
                "FEATURE": feature,
                "READS": reads,
                "PERCENTAGE": percentage
            }
            dictwriter_object = csv.DictWriter(f, header)
            dictwriter_object.writerow(row)

    @staticmethod
    def _parse_output_path(output_path) -> str:
        """Generates output filepath in case it is not given or invalid"""
        DEFAULT_FILENAME = 'compare_output.csv'
        if not output_path:
            output_path = f'./{DEFAULT_FILENAME}'
        else:
            head, tail = os.path.split(output_path)
            if not tail:
                tail = DEFAULT_FILENAME
            if not tail.endswith('.csv'):
                tail += '.csv'
            output_path = os.path.join(head, tail)
        return output_path

    def _get_extract_path(self, read_type: str):
        return self.output_path[:-4] + f".{read_type}.csv"

    def save_to_csv(self, df: pd.DataFrame,
                    read_type: str
                    ):
        columns = "|".join(self.filtered_columns)
        csv_path = self._get_extract_path(read_type)
        df.filter(regex=columns).to_csv(csv_path)

    def _check_input(self) -> None:
        """Raises exception if input lacks required data"""
        self.original_data = ComparisonObject(self.input_data['original'])
        self.reversed_data = ComparisonObject(self.input_data['reversed'])
        self.shuffled_data = ComparisonObject(self.input_data['shuffled'])
        self.samples = [self.original_data,
                        self.reversed_data,
                        self.shuffled_data]
        self.comparable = [sample.is_comparable for sample in self.samples]
        if not self.original_data.is_comparable:
            print(self.original_data)
            raise Exception('Original sample data not provided')
        if sum(self.comparable) < 2:
            raise Exception('At least 2 samples including original required')

    def merge_dataframes(self) -> pd.DataFrame:
        """Merges 2 or 3 dataframes into one in case of multiple tables"""
        paths_to_csv = []
        for sample in self.samples:
            if sample.is_comparable:
                paths_to_csv.append(sample.path_to_csv)
        norm_paths = [os.path.normpath(path) for path in paths_to_csv]
        unique_paths = list(dict.fromkeys(norm_paths))
        # create list of dataframes
        print('Importing CSV(s)...')
        self.imported_dfs = [
            pd.read_csv(path, index_col=0)
            for path in unique_paths
        ]
        print('\nTables successfully imported')
        print(f'RAM usage: {self.get_current_memory_usage()}')
        # merging if there are multiple tables
        if len(self.imported_dfs) == 1:
            return self.imported_dfs[0]
        elif len(self.imported_dfs) == 2:
            return pd.merge(
                *self.imported_dfs,
                left_index=True,
                right_index=True
            )
        elif len(self.imported_dfs) == 3:
            return pd.merge(self.imported_dfs[0],
                            self.imported_dfs[1],
                            left_index=True,
                            right_index=True
                            ).merge(self.imported_dfs[2],
                                    left_index=True,
                                    right_index=True)

    @staticmethod
    def is_single_ended(dataframe: pd.DataFrame) -> bool:
        first_row_values = dataframe.iloc[0].values[-1]
        return not isinstance(first_row_values, str)

    def get_labels(self) -> list:
        labels = [self.original_data.label]
        if self.reversed_data.is_comparable:
            labels.append(self.reversed_data.label)
        if self.shuffled_data.is_comparable:
            labels.append(self.shuffled_data.label)
        return labels

    def count_total_reads(self, df: pd.DataFrame):
        """Counts total number of reads for the unfiltered raw dataframe,
        writes the first row to the csv"""
        # total_reads = df.shape[0]
        if sum(self.comparable) == 2:
            # total_reads_type1 = self.imported_dfs[0].shape[0]
            # total_reads_type2 = self.imported_dfs[1].shape[0]

            self.write_row(feature='FASTQ_reads',
                           reads=self.fastq_reads,
                           total_reads=self.fastq_reads)

    @staticmethod
    def get_current_memory_usage(in_gigabytes=True) -> float:
        """Checks the amount of RAM used by the script at the time"""
        process = psutil.Process(os.getpid())
        usage = process.memory_info().rss  # in bytes
        usage_GB = usage * 9.31 * 10 ** (-10)  # in Gigabytes
        return round(usage_GB, 3) if in_gigabytes else usage

    def count_unmapped(self, df: pd.DataFrame) -> None:
        if sum(self.comparable) != 2:
            raise ValueError("Comparison requires exactly two comparable datasets.")

        labels = self.get_labels()
        mapping_types = {
            'unmapped': ['U', 'U'],
            'forward_mapped': ['P', 'U'],
            'backward_mapped': ['U', 'P']
        }

        # Initialize a dictionary to store dataframes for each type and label
        dfs = {mtype: {label: None for label in labels} for mtype in mapping_types}

        # Iterate over each mapping type and label
        for mtype, mvalue in mapping_types.items():
            for label in labels:
                dfs[mtype][label] = df[df[f'{label}_type'].apply(lambda x: ast.literal_eval(x) == mvalue)]

                if '_g' in label:
                    # Log the count of reads per type and label
                    self.write_row(feature=f"{mtype}_type1",
                                   reads=dfs[mtype][label].shape[0],
                                   total_reads=self.fastq_reads)

                else:
                    self.write_row(feature=f"{mtype}_type2",
                                   reads=dfs[mtype][label].shape[0],
                                   total_reads=self.fastq_reads)

    # Define a function to check for unmapped to mapped transitions
    def count_inconsistent(self, df: pd.DataFrame):

        l1, l2 = self.get_labels()

        df_inconsistent_type1 = df.loc[
            df[f'{l1}_type'].apply(lambda x: ast.literal_eval(x) == ['P', 'P']) &
            df[f'{l2}_type'].apply(
                lambda x: (
                        ast.literal_eval(x) == ['U', 'P'] or
                        ast.literal_eval(x) == ['U', 'P'] or
                        ast.literal_eval(x) == ['P', 'U']
                )
            )
            ]

        df_inconsistent_type2 = df.loc[
            df[f'{l1}_type'].apply(lambda x: ast.literal_eval(x) == ['U', 'U']) &
            df[f'{l2}_type'].apply(
                lambda x: (
                        ast.literal_eval(x) == ['P', 'U'] or
                        ast.literal_eval(x) == ['U', 'P'] or
                        ast.literal_eval(x) == ['P', 'P']
                )
            )
            ]

        self.write_row(feature='inconsistent_type1',
                       reads=df_inconsistent_type1.shape[0],
                       total_reads=self.fastq_reads)

        self.write_row(feature='inconsistent_type2',
                       reads=df_inconsistent_type2.shape[0],
                       total_reads=self.fastq_reads)

        if len(df_inconsistent_type1) > 0: print('inconsistent_type1: ', df_inconsistent_type1.index[0])
        if len(df_inconsistent_type2) > 0: print('inconsistent_type2: ', df_inconsistent_type2.index[0])

        # self.write_row(feature='backward_to_forward',
        #                reads=df_backward_to_forward.shape[0],
        #                total_reads=self.fastq_reads)

    def extract_proper(self, df: pd.DataFrame) -> pd.DataFrame:

        l1, l2 = self.get_labels()

        df_proper_g = df[df[f'{l1}_proper_pair'] == 1]
        df_proper_rep = df[df[f'{l2}_proper_pair'] == 1]

        # Log the count of reads per type and label
        self.write_row(feature=f"proper_type1",
                       reads=df_proper_g.shape[0],
                       total_reads=self.fastq_reads)

        self.write_row(feature=f"proper_type2",
                       reads=df_proper_rep.shape[0],
                       total_reads=self.fastq_reads)

        df_common_proper = df_proper_g.loc[df_proper_g.index.isin(df_proper_rep.index)]

        self.write_row(feature=f"common_proper",
                       reads=df_common_proper.shape[0],
                       total_reads=self.fastq_reads)

        return df_common_proper

    def extract_common_mapped(self, df: pd.DataFrame) -> pd.DataFrame:

        l1, l2 = self.get_labels()

        df_mapped_g = df.loc[df[f'{l1}_type'].apply(lambda x: ast.literal_eval(x) == ['P', 'P'])]
        df_mapped_rep = df.loc[df[f'{l2}_type'].apply(lambda x: ast.literal_eval(x) == ['P', 'P'])]

        # Log the count of reads per type and label
        self.write_row(feature=f"mapped_type1",
                       reads=df_mapped_g.shape[0],
                       total_reads=self.fastq_reads)

        self.write_row(feature=f"mapped_type2",
                       reads=df_mapped_rep.shape[0],
                       total_reads=self.fastq_reads)

        df_common_mapped = df_mapped_g.loc[df_mapped_g.index.isin(df_mapped_rep.index)]

        self.write_row(feature=f"common_mapped",
                       reads=df_common_mapped.shape[0],
                       total_reads=self.fastq_reads)

        if len(df_common_mapped) > 0 : print('common_mapped: ', df_common_mapped.index[0])

        return df_common_mapped

    def count_nonidentical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Counts identical reads by position and edit distance,
        writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_nonidentical = df.loc[
                (df[f'{l1}_pos'] !=
                 df[f'{l2}_pos']) |
                (df[f'{l1}_edit_dist'] !=
                 df[f'{l2}_edit_dist'])
                ]
            # if len(df_nonidentical) > 0:
            #     print(df_nonidentical.index)

            self.write_row(feature='Non-identical',
                           reads=df_nonidentical.shape[0],
                           total_reads=df.shape[0])

            if len(df_nonidentical) > 0 : print('non_identical: ', df_nonidentical.index[0])

            return df_nonidentical

    def count_identical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Counts identical reads by position and edit distance,
        writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_identical = df.loc[
                (df[f'{l1}_pos'] ==
                 df[f'{l2}_pos']) &
                (df[f'{l1}_edit_dist'] ==
                 df[f'{l2}_edit_dist'])
                ]
            # self.write_row(feature='Non-identical',
            #                reads=df_identical.shape[0],
            #                total_reads=self.fastq_reads)
        return df_identical

    def count_CG_IL(self, df: pd.DataFrame, df_non: pd.DataFrame):
        """Counts reads with consistent global and inconsistent local
                alignment (common proper reads mapped to the same position
                with different edit distance), writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_CG_IL = df.loc[
                (df[f'{l1}_pos'] ==
                 df[f'{l2}_pos']) &
                (df[f'{l1}_edit_dist'] !=
                 df[f'{l2}_edit_dist'])
                ]
            self.write_row(feature='Consistent_global_inconsistent_local',
                           reads=df_CG_IL.shape[0],
                           total_reads=df_non.shape[0])
            if len(df_CG_IL) > 0: print('CG_IL', df_CG_IL.index[0])

        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()

            def count_CG_IL_between(label1: str, label2: str):
                df_CG_IL = df.loc[
                    (df[f'{label1}_pos'] ==
                     df[f'{label2}_pos']) &
                    (df[f'{label1}_edit_dist'] !=
                     df[f'{label2}_edit_dist'])
                    ]
                self.write_row(feature=f'CG_IL_{label1}_AND_{label2}',
                               reads=df_CG_IL.shape[0],
                               total_reads=df_non.shape[0])

            count_CG_IL_between(l1, l2)
            count_CG_IL_between(l1, l3)
            count_CG_IL_between(l2, l3)

        if any(x in self.reads_to_extract for x in ["all", "CG_IL"]):
            self.save_to_csv(df_CG_IL, "CG_IL")

    def count_IG(self, df: pd.DataFrame, df_non: pd.DataFrame):
        """Counts reads with inconsistent global alignment (common
        unambiguous reads mapped to the different position with
        different edit distance), writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_IG = df.loc[
                (df[f'{l1}_pos'] !=
                 df[f'{l2}_pos']) &
                (df[f'{l1}_edit_dist'] !=
                 df[f'{l2}_edit_dist'])
                ]
            self.write_row(feature='Inconsistent_global',
                           reads=df_IG.shape[0],
                           total_reads=df_non.shape[0])

            if len(df_IG) > 0: print('IG', df_IG.index[0])

        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()

            def count_IG_between(label1: str, label2: str):
                df_IG = df.loc[
                    (df[f'{label1}_pos'] !=
                     df[f'{label2}_pos']) &
                    (df[f'{label1}_edit_dist'] !=
                     df[f'{label2}_edit_dist'])
                    ]
                self.write_row(feature=f'IG_{label1}_AND_{label2}',
                               reads=df_IG.shape[0],
                               total_reads=df_non.shape[0])

            count_IG_between(l1, l2)
            count_IG_between(l1, l3)
            count_IG_between(l2, l3)

        if any(x in self.reads_to_extract for x in ["all", "IG"]):
            self.save_to_csv(df_IG, "IG")

    def count_MM(self, df: pd.DataFrame, df_non: pd.DataFrame):
        """Counts multi-mapped reads (common unambiguous reads mapped
        to different positions with the same edit distance),
        writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_MM = df.loc[
                (df[f'{l1}_pos'] !=
                 df[f'{l2}_pos']) &
                (df[f'{l1}_edit_dist'] ==
                 df[f'{l2}_edit_dist'])
                ]
            self.write_row(feature='multi_mapped',
                           reads=df_MM.shape[0],
                           total_reads=df_non.shape[0])

            if len(df_MM) > 0: print('MM: ', df_MM.index[0])


        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()

            def count_MM_between(label1: str, label2: str):
                df_MM = df.loc[
                    (df[f'{label1}_pos'] !=
                     df[f'{label2}_pos']) &
                    (df[f'{label1}_edit_dist'] ==
                     df[f'{label2}_edit_dist'])
                    ]
                self.write_row(feature='multi_mapped_' +
                                       f'{label1}_AND_{label2}',
                               reads=df_MM.shape[0],

                               total_reads=df.shape[0])

            count_MM_between(l1, l2)
            count_MM_between(l1, l3)
            count_MM_between(l2, l3)

        if any(x in self.reads_to_extract for x in ["all", "MM"]):
            self.save_to_csv(df_MM, "MM")

    def count_diff_qual(self, df: pd.DataFrame):
        """Counts identical reads by position and edit distance,
        writes output row to the csv"""
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()

            df_qual = df.loc[
                (df[f'{l1}_quality'] !=
                 df[f'{l2}_quality'])
            ]
            self.write_row(feature='quality_change',
                           reads=df_qual.shape[0],
                           total_reads=self.fastq_reads)

            if len(df_qual) > 0: print('quality_change', df_qual.index[0])

    def find_and_count_highest_quality_proper_mapped_reads(self, df: pd.DataFrame):
        # Retrieve labels dynamically
        l1, l2 = self.get_labels()

        # Check if quality columns contain valid data
        if df[f'{l1}_quality'].isnull().all() and df[f'{l2}_quality'].isnull().all():
            self.write_row(feature=f'high_quality_mapped_type1', reads="NA", total_reads=self.fastq_reads)
            self.write_row(feature=f'high_quality_mapped_type2', reads="NA", total_reads=self.fastq_reads)
            print("All quality values are None for both datasets. Skipping calculation.")
            return

        # Assuming quality columns store string representations of lists
        df[f'{l1}_parsed_quality'] = df[f'{l1}_quality'].apply(
            lambda x: ast.literal_eval(x) if pd.notnull(x) else []
        )
        df[f'{l2}_parsed_quality'] = df[f'{l2}_quality'].apply(
            lambda x: ast.literal_eval(x) if pd.notnull(x) else []
        )

        # Initialize highest_quality
        highest_quality = None

        # Find the highest quality for l1
        highest_quality_l1 = max(
            (max(row) for row in df[f'{l1}_parsed_quality'] if row),
            default=None
        )

        # Find the highest quality for l2
        highest_quality_l2 = max(
            (max(row) for row in df[f'{l2}_parsed_quality'] if row),
            default=None
        )

        # Determine the overall highest quality
        if highest_quality_l1 is not None or highest_quality_l2 is not None:
            highest_quality = max(
                filter(None, [highest_quality_l1, highest_quality_l2])
            )

        # Define a lambda to check if all items in the list equal the highest quality
        is_both_highest = lambda x: all(item == highest_quality for item in x)

        # Filter for rows where both quality scores match the highest quality and are marked as properly paired
        df_high_quality_proper_l1 = df[
            df[f'{l1}_parsed_quality'].apply(is_both_highest) & (df[f'{l1}_proper_pair'] == 1)]
        df_high_quality_proper_l2 = df[
            df[f'{l2}_parsed_quality'].apply(is_both_highest) & (df[f'{l2}_proper_pair'] == 1)]

        # Log the number of high-quality, properly mapped reads for each dataset
        self.write_row(feature=f'high_quality_mapped_type1',
                       reads=df_high_quality_proper_l1.shape[0],
                       total_reads=self.fastq_reads)
        self.write_row(feature=f'high_quality_mapped_type2',
                       reads=df_high_quality_proper_l2.shape[0],
                       total_reads=self.fastq_reads)

        if len(df_high_quality_proper_l1) > 0:
            print('high_quality_mapped_type1: ', df_high_quality_proper_l1.index[0])
        if len(df_high_quality_proper_l2) > 0:
            print('high_quality_mapped_type2: ', df_high_quality_proper_l2.index[0])

    def compare(self):
        df = self.merge_dataframes()
        self.is_single_ended = self.is_single_ended(df)
        # 2. Counting total number of reads in the raw dataframe
        self.count_total_reads(df)
        self.count_unmapped(df)
        # df_common_proper = self.extract_proper(df)
        df_common_mapped = self.extract_common_mapped(df)

        self.count_inconsistent(df)
        df_nonidentical = self.count_nonidentical(df_common_mapped)
        self.count_CG_IL(df_common_mapped, df_nonidentical)
        self.count_IG(df_common_mapped, df_nonidentical)
        self.count_MM(df_common_mapped, df_nonidentical)

        df_identical = self.count_identical(df_common_mapped)
        self.count_diff_qual(df_identical)
        self.find_and_count_highest_quality_proper_mapped_reads(df_common_mapped)


app = typer.Typer(add_completion=False)


@app.command()
def main(input_csvs: List[Path] = typer.Argument(
    ...,
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
    show_default=False,
    help="List of CSV files to compare; min: 2 files"
),
        fastq_file: Path = typer.Argument(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            show_default=False,
            help="Path to FASTQ File of the sample"
        ),
        output_path: Path = typer.Option(
            "./comparer_output.csv",
            "--output",
            "-o",
            help="Path to output directory or file."
        ),
        reads_to_extract: Optional[List[str]] = typer.Option(
            [],
            "--extract",
            "-x",
            help="Extract specific type of common reads"
        ),
        filtered_columns: Optional[List[str]] = typer.Option(
            CSV.all_columns,
            "--filter",
            "-f",
            help="Filter columns for extraction of specific read types"
        )

):
    """
    Subsampling and comparing tool for parsed CSV files.

    """
    read_types_notation = {
        "all": "Extract every read type listed below "
               "to an individual file.",
        "ID": "Identical. Common reads identical by "
              "position and edit distance.",
        "CG_IL": "Consistent Global and Inconsistent Local. "
                 "Common unambiguous reads mapped to the same "
                 "position with different edit distance.",
        "IT1": "Inconsistent Type 1. Common unambiguous reads "
               "mapped only with original data (Type 1).",
        "IT2": "Inconsistent Type 1. Common unambiguous reads "
               "only with replicated data (Type 2).",
        "IG": "Inconsistent Global. Common unambiguous reads "
              "mapped to the different position with "
              "different edit distance.",
        "MM": "Hidden 'Multi-Mapped'. Common unambiguous reads mapped "
              "to different positions with the same edit distance."
    }

    column_names_notation = {
        "flags": "SAM format flag for each read",
        "pos": "1-based leftmost mapping POSition",
        "chr": "Chromosome",
        "CIGAR": "CIGAR string",
        "quality": "ASCII of Phred-scaled base QUALity+33",
        "edit_dist": "Edit distance",
        "MD": "Mismatch - MD tag (if available)",
        "type": "mapping (P - primary mapped / S - secondary mapped / U - unmapped)",
        "proper_pair": "proper pair mapping (1- mapped in proper pair / 0- not mapped in proper pair)",
    }

    # extraction - input check
    for read_type in reads_to_extract:
        if read_type not in read_types_notation.keys():
            print("INCORRECT INPUT. No such read type:", read_type)
            print("Only the following types of reads",
                  "are available for extraction:",
                  ", ".join(read_types_notation.keys()), "\n")
            for _type, description in read_types_notation.items():
                print(f"{_type}: \t{description}")
            print()
            raise typer.Abort()

    # filtering - input check
    for column in filtered_columns:
        if column not in column_names_notation.keys():
            print("INCORRECT INPUT. No such column name:", column)
            print("Only the following column names",
                  "are available for filtering:",
                  ", ".join(column_names_notation.keys()), "\n")
            for _type, description in column_names_notation.items():
                print(f"{_type}: \t{description}")
            print()
            raise typer.Abort()

    fastq_reads = count_reads(fastq_file)

    csvs = [CSV(str(input_csv)) for input_csv in input_csvs]
    paths = input_csvs

    input_data = {
        'original': {
            'label': csvs[0].label,
            'path_to_csv': str(paths[0])
        },

        'reversed': {
            'label': csvs[1].label,
            'path_to_csv': str(paths[1])
        },

        'shuffled': {
            'label': '',
            'path_to_csv': ''
        },
    }

    if len(csvs) > 2:
        input_data['shuffled'] = {
            'label': csvs[2].label,
            'path_to_csv': str(paths[2])
        }

    comparer = Comparer(input_data,
                        output_path=str(output_path),
                        reads_to_extract=reads_to_extract,
                        filtered_columns=filtered_columns,
                        fastq_reads=fastq_reads)
    comparer.compare()


if __name__ == "__main__":
    app()
