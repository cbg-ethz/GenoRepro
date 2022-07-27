from dataclasses import dataclass
import os.path
import psutil
import pandas as pd
import csv


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


    @staticmethod
    def crop_read_id(index: str):
        return index.split('.')[1]


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
        dataframes = [pd.read_csv(path) for path in unique_paths]
        for df in dataframes:
            df['Unnamed: 0'] = df['Unnamed: 0'].map(self.crop_read_id)
        print('\nTables successfully imported')
        print(f'RAM usage: {self.get_current_memory_usage()}')
        # merging if there are multiple tables
        if len(dataframes) == 1:
            return dataframes[0]
        elif len(dataframes) == 2:
            return pd.merge(*dataframes)
        elif len(dataframes) == 3:
            return pd.merge(*dataframes[:2]).merge(dataframes[2])


    @staticmethod
    def is_single_ended(dataframe: pd.DataFrame) -> bool:
        first_row_values = dataframe.loc[0].values[-1]
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
