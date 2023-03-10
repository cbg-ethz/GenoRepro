import os
from dataclasses import dataclass
import os.path
from pathlib import Path
from typing import List
import psutil
import pandas as pd
from Bio.Seq import Seq
from csvsort import csvsort
import typer
import csv
from itertools import combinations


class CSV:
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


class Subsampler:
    def __init__(self, *csv_files: CSV) -> None:
        self.csvs_g = csv_files
        self.csv_rv_paths = [
            csv_file.create_rv_csv() 
            for csv_file in csv_files
                ]
        self.csvs_rv = [
            CSV(path, is_original=False)
            for path in self.csv_rv_paths
            ]


    @staticmethod
    def find_common_reads(csvs: List[CSV], column: int = -1):
        def next_row(readers):
            return [next(reader) for reader in readers]
        
        readers = [i_csv.get_reader() for i_csv in csvs]
        # starting with header row
        current_rows = next_row(readers)
        # skipping header row to process 1st read
        current_rows = next_row(readers)

        while True:
            try:
                for i, reader in enumerate(readers):
                    while current_rows[i-1][column] > current_rows[i][column]:
                        current_rows[i] = next(reader)
                    if len(set(row[column] for row in current_rows)) == 1:
                        yield current_rows
                        current_rows = [next(reader) for reader in readers]
            except StopIteration:
                break
        
        # closing previously opened CSV files 
        for i_csv in csvs:
            i_csv._close()


    @staticmethod
    def get_combinations(iterable1, iterable2) -> list:
        """Generator that outputs combination of given iterables"""
        yield list(iterable1)
        length = len(iterable1)
        if length == 2:
            yield [iterable1[0], iterable2[1]]
            return
        length_range = range(length)
        iterable1_indeces = []
        current_length = len(iterable1)
        while current_length != 2:
            combs = combinations(range(current_length), current_length - 1)
            iterable1_indeces += list(combs)
            current_length -= 1
        iterable2_indeces = [length_range for _ in iterable1_indeces]
        iterable2_indeces_updated = []
        for i1, i2 in zip(iterable1_indeces, iterable2_indeces):
            difference = set(i2).difference(i1)
            new_i2 = tuple(i for i in i2 if i in difference)
            iterable2_indeces_updated.append(new_i2)
        for i1, i2 in zip(iterable1_indeces, iterable2_indeces_updated):
            yield [iterable1[i] for i in i1] + [iterable2[i] for i in i2]


    def run(self):
        for combinaiton in self.get_combinations(self.csvs_g, self.csvs_rv):
            for common_reads in self.find_common_reads(combinaiton):
                # reduce the names of the reads to one common format
                read_names = sorted(read[0] for read in common_reads)
                combined_read_name = "_".join(read_names)
                # iterate and append common read to the corresponding sub-CSV
                for read, csv_file in zip(common_reads, combinaiton):
                    read[0] = combined_read_name
                    csv_file.append_row(read[:-1]) # excluding sequence column
        # append rows left in the buffer
        for csv_file in tuple(self.csvs_g) + tuple(self.csvs_rv):
            csv_file.finalize_append()


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
    df_unambiguous : pd.DataFrame
        Pandas DataFrame of unambiguous reads
    """


    def __init__(self, input_data: dict, output_path: str = '') -> None:
        self.input_data = input_data
        self.output_path = self._parse_output_path(output_path)
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
        with open(self.output_path, 'a') as f:
            header = ['FEATURE', 'READS', 'PERCENTAGE']
            row = {
                "FEATURE": feature,
                "READS": reads,
                "PERCENTAGE": round(100 * reads / total_reads, 3)
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
        dataframes = [pd.read_csv(path, index_col=0) for path in unique_paths]
        print('\nTables successfully imported')
        print(f'RAM usage: {self.get_current_memory_usage()}')
        # merging if there are multiple tables
        if len(dataframes) == 1:
            return dataframes[0]
        elif len(dataframes) == 2:
            return pd.merge(
                *dataframes,
                left_index=True, 
                right_index=True
                )
        elif len(dataframes) == 3:
            return pd.merge(dataframes[0],
                            dataframes[1], 
                            left_index=True, 
                            right_index=True
                ).merge(dataframes[2], 
                        left_index=True, 
                        right_index=True)


    @staticmethod
    def is_single_ended(dataframe: pd.DataFrame) -> bool:
        first_row_values = dataframe.iloc[0].values[-1]
        return not isinstance(first_row_values, str)

    
    @staticmethod
    def get_current_memory_usage(in_gigabytes=True) -> float:
        """Checks the amount of RAM used by the script at the time"""
        process = psutil.Process(os.getpid())
        usage = process.memory_info().rss # in bytes
        usage_GB = usage * 9.31 * 10 ** (-10) # in Gigabytes
        return round(usage_GB, 3) if in_gigabytes else usage


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
        total_reads = df.shape[0]
        self.write_row(feature='Total_reads', 
                       reads=total_reads,
                       total_reads=total_reads)


    def extract_mapped(self, df: pd.DataFrame) -> pd.DataFrame:
        """Removes those reads which were mapped neither with original
        nor replicated data"""
        # excluded values for '_pos' column
        if self.is_single_ended:
            val = [0, 'None']
        else:
            val = ['[0, None]', '[None, 0]', '[0, 0]', '[None, None]']
        
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_mapped = df.loc[
                (~df[f'{l1}_pos'].isin(val)) |
                (~df[f'{l2}_pos'].isin(val))
                ]
            self.write_row(feature=f'Mapped_reads', 
                           reads=df_mapped.shape[0],
                           total_reads=df.shape[0])
        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()
            df_mapped = df.loc[
                (~df[f'{l1}_pos'].isin(val)) |
                (~df[f'{l2}_pos'].isin(val)) |
                (~df[f'{l3}_pos'].isin(val))
                ]
            self.write_row(feature=f'Mapped_reads', 
                           reads=df_mapped.shape[0],
                           total_reads=df.shape[0])
        return df_mapped


    def remove_inconsistent(self, df: pd.DataFrame) -> pd.DataFrame:
        """Counts reads belonging to Inconsistent type 1 and 2, writes 
        corresponding rows to the csv, removes unambiguous reads mapped only 
        with original data (type 1) and only with replicated data (type 2)"""
        # excluded values for '_pos' column
        if self.is_single_ended:
            val = [0, 'None']
        else:
            val = ['[0, None]', '[None, 0]', '[0, 0]', '[None, None]']
        
        def find_inconsistent(label1: str, label2: str, type_: int,
                              triple_wise: bool):
            """type_ can only take values 1 and 2 (inconsistent type)"""
            df_IT = df.loc[
                (~df[f'{label1}_pos'].isin(val)) &
                (df[f'{label2}_pos'].isin(val))
                ]
            feature = f'Inconsistent_type{type_}'
            if triple_wise:
                sorted_labels = sorted([label1, label2])
                feature += f'_{sorted_labels[0]}_AND_{sorted_labels[1]}'
            self.write_row(feature=feature, 
                           reads=df_IT.shape[0],
                           total_reads=df.shape[0])

        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()

            find_inconsistent(l1, l2, type_=1, triple_wise=False)
            find_inconsistent(l2, l1, type_=2, triple_wise=False)

            df_only_mapped = df.loc[
                (~df[f'{l1}_pos'].isin(val)) &
                (~df[f'{l2}_pos'].isin(val))
                ]
        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()

            find_inconsistent(l1, l2, type_=1, triple_wise=True)
            find_inconsistent(l2, l1, type_=2, triple_wise=True)

            find_inconsistent(l1, l3, type_=1, triple_wise=True)
            find_inconsistent(l3, l1, type_=2, triple_wise=True)

            df_only_mapped = df.loc[
                (~df[f'{l1}_pos'].isin(val)) &
                (~df[f'{l2}_pos'].isin(val)) &
                (~df[f'{l3}_pos'].isin(val))
                ]
        return df_only_mapped


    def count_unamiguous(self, df: pd.DataFrame):
        """Counts nubmer of unamiguous reads for every sample separetely,
        writes two or three output rows to the csv"""
        # allowed values for '_multi' flag
        if self.is_single_ended:
            val = [0]
        else:
            val = ['[0, None]', '[None, 0]', '[0, 0]']
        
        l1 = self.original_data.label
        v1 = df.loc[df[f'{l1}_multi'].isin(val)].shape[0]
        self.write_row(f'Unambiguous_{l1}', v1, df.shape[0])

        if self.reversed_data.is_comparable:
            l2 = self.reversed_data.label
            v2 = df.loc[df[f'{l2}_multi'].isin(val)].shape[0]
            self.write_row(f'Unambiguous_{l2}', v2, df.shape[0])

        if self.shuffled_data.is_comparable:
            l3 = self.shuffled_data.label
            v3 = df.loc[df[f'{l3}_multi'].isin(val)].shape[0]
            self.write_row(f'Unambiguous_{l3}', v3, df.shape[0])


    def extract_unambiguous(self, df: pd.DataFrame) -> pd.DataFrame:
        """Returns dataframe with only common unambiguous reads,
        writes output row to the csv"""
        # allowed values for '_multi' flag
        if self.is_single_ended:
            val = [0]
        else:
            val = ['[0, None]', '[None, 0]', '[0, 0]']
        
        if sum(self.comparable) == 2:
            l1, l2 = self.get_labels()
            df_unambiguous = df.loc[
                df[f'{l1}_multi'].isin(val) & 
                df[f'{l2}_multi'].isin(val)
                ]
        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()
            df_unambiguous = df.loc[
                df[f'{l1}_multi'].isin(val) & 
                df[f'{l2}_multi'].isin(val) &
                df[f'{l3}_multi'].isin(val)
                ]
        self.write_row(feature='Common_unambiguous', 
                       reads=df_unambiguous.shape[0],
                       total_reads=df.shape[0])
        self.df_unambiguous = df_unambiguous
        return df_unambiguous


    def count_identical(self, df: pd.DataFrame):
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
            self.write_row(feature='Identical', 
                           reads=df_identical.shape[0],
                           total_reads=self.df_unambiguous.shape[0])

        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()
            df_identical = df.loc[
                (df[f'{l1}_pos'] == 
                df[f'{l2}_pos']) &
                (df[f'{l2}_pos'] == 
                df[f'{l3}_pos']) &
                (df[f'{l1}_edit_dist'] == 
                df[f'{l2}_edit_dist']) &
                (df[f'{l2}_edit_dist'] == 
                df[f'{l3}_edit_dist'])
                ]
            self.write_row(feature='Identical_among_all_three', 
                           reads=df_identical.shape[0],
                           total_reads=self.df_unambiguous.shape[0])
            
            def count_identical_between(label1: str, label2: str):
                df_identical = df.loc[
                    (df[f'{label1}_pos'] == 
                    df[f'{label2}_pos']) &
                    (df[f'{label1}_edit_dist'] == 
                    df[f'{label2}_edit_dist'])
                    ]
                self.write_row(feature=f'Identical_{label1}_AND_{label2}',
                            reads=df_identical.shape[0],
                            total_reads=self.df_unambiguous.shape[0])
            
            count_identical_between(l1, l2)
            count_identical_between(l1, l3)
            count_identical_between(l2, l3)


    def count_CG_IL(self, df: pd.DataFrame):
        """Counts reads with consistent global and inconsistent local
        alignment (common unambiguous reads mapped to the same position
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
                           total_reads=self.df_unambiguous.shape[0])
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
                               total_reads=self.df_unambiguous.shape[0])

            count_CG_IL_between(l1, l2)
            count_CG_IL_between(l1, l3)
            count_CG_IL_between(l2, l3)


    def count_IG(self, df: pd.DataFrame):
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
                           total_reads=self.df_unambiguous.shape[0])
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
                               total_reads=self.df_unambiguous.shape[0])

            count_IG_between(l1, l2)
            count_IG_between(l1, l3)
            count_IG_between(l2, l3)


    def count_MM(self, df: pd.DataFrame):
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
            self.write_row(feature='Multi_mapped',
                           reads=df_MM.shape[0],
                           total_reads=self.df_unambiguous.shape[0])
        elif sum(self.comparable) == 3:
            l1, l2, l3 = self.get_labels()

            def count_MM_between(label1: str, label2: str):
                df_MM = df.loc[
                    (df[f'{label1}_pos'] != 
                    df[f'{label2}_pos']) &
                    (df[f'{label1}_edit_dist'] == 
                    df[f'{label2}_edit_dist'])
                    ]
                self.write_row(feature=f'Multi_mapped_{label1}_AND_{label2}',
                               reads=df_MM.shape[0],
                               total_reads=self.df_unambiguous.shape[0])

            count_MM_between(l1, l2)
            count_MM_between(l1, l3)
            count_MM_between(l2, l3)


    def compare(self):
        # 1. Merging dataframes from different csv if needed
        df = self.merge_dataframes()
        self.is_single_ended = self.is_single_ended(df)
        # 2. Counting total number of reads in the raw dataframe
        self.count_total_reads(df)
        # 3. Removing reads which were mapped neither with original
        #    nor replicated data
        df_mapped = self.extract_mapped(df)
        # 4. Counting unambiguous reads and creating dataframe with them
        self.count_unamiguous(df_mapped)
        df_unambiguous = self.extract_unambiguous(df_mapped)
        # 5. Counting inconsistent type 1 and 2
        df_without_IT1_IT2 = self.remove_inconsistent(df_unambiguous)
        # 6. Counting other features among filtered reads
        self.count_identical(df_without_IT1_IT2)
        self.count_CG_IL(df_without_IT1_IT2)
        self.count_IG(df_without_IT1_IT2)
        self.count_MM(df_without_IT1_IT2)
        
        print('\nComparision complete')
        print(f'RAM usage: {self.get_current_memory_usage()}')


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
         output_path: Path = typer.Option(
            "./comparer_output.csv",
            "--output",
            "-o",
            help="Path to output directory or file."
            )
            ):
    """
    Python script that parses Binary Alignment Map File,
    keeping only primary alignments.
    
    """

    csvs = [CSV(str(input_csv)) for input_csv in input_csvs]
    paths = input_csvs

    if csvs[0].is_real_data:
        subsampler = Subsampler(*csvs)
        subsampler.run()
        paths = [csv_file.sub for csv_file in csvs]

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
    
    comparer = Comparer(input_data, output_path=str(output_path))
    comparer.compare()

if __name__ == "__main__":
    app()