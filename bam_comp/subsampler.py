import csv
import pandas as pd
from Bio.Seq import Seq
import psutil
from os import getpid


class CSV:
    def __init__(self, path_to_csv) -> None:
        self.path_to_csv = path_to_csv
        self.is_single_ended = True
        with open(self.path_to_csv, 'r') as csv_file:
            reader = csv.reader(csv_file, delimiter=',')
            next(reader) # skip header row
            if '[' in next(reader)[-1]:
                self.is_single_ended = False


    def create_sub_csv(self, clone_header=True):
        new_csv_path = self.path_to_csv[:-4] + '_sub.csv'
        with open(self.path_to_csv, 'r') as csv_master:
            with open(new_csv_path, 'w') as csv_sub:
                if clone_header:
                    reader = csv.reader(csv_master, delimiter=',')
                    writer = csv.writer(csv_sub, delimiter=',')
                    header = next(reader)
                    writer.writerow(header[:-1])
        return new_csv_path
    


class Subsampler:
    def __init__(self, csv_1: CSV, csv_2: CSV) -> None:
        self.csv_1 = csv_1
        self.csv_2 = csv_2
        self.is_single_ended = False
        if csv_1.is_single_ended and csv_2.is_single_ended:
            self.is_single_ended = True

    @staticmethod
    def get_current_memory_usage(in_gigabytes=True) -> float:
        """Checks the amount of RAM used by the script at the time"""
        process = psutil.Process(getpid())
        usage = process.memory_info().rss # in bytes
        usage_GB = usage * 9.31 * 10 ** (-10) # in Gigabytes
        return round(usage_GB, 3) if in_gigabytes else usage

    @staticmethod
    def import_dataframe(path_to_csv) -> pd.DataFrame:
        df = pd.DataFrame()
        for chunk in pd.read_csv(path_to_csv, chunksize=1000000):
            df = pd.concat([df, chunk], ignore_index=True)
        return df

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        return Seq(sequence).reverse_complement()

    @staticmethod
    def get_paired_end_sequence(string: str) -> str:
        return string.split("'")[1]
    
    @staticmethod
    def format_paired_end(sequence: str) -> str:
        return f"['{sequence}', None]"

    def run(self):
        # 1) assigning specific method to parse sequence from SE/PE data
        if self.is_single_ended:
            get_sequence = lambda string: string
            format_sequence = lambda sequence: sequence
        else:
            get_sequence = self.get_paired_end_sequence
            format_sequence = self.format_paired_end
        
        def get_reverse_complement(sequence: str) -> str:
            return str(format_sequence(self.reverse_complement(get_sequence(sequence))))

        
        df_csv_1 = self.import_dataframe(self.csv_1.path_to_csv)
        seq_column_name_1 = df_csv_1.columns.values.tolist()[-1]

        df_csv_2 = self.import_dataframe(self.csv_2.path_to_csv)
        seq_column_name_2 = df_csv_2.columns.values.tolist()[-1]
        
        df_csv_1_rv = df_csv_1.copy()
        df_csv_1_rv[seq_column_name_1] = df_csv_1_rv[seq_column_name_1].apply(get_reverse_complement)

        df_csv_1_g_in_2 = df_csv_1[df_csv_1.iloc[:, -1].isin(df_csv_2.iloc[:, -1])].sort_values(seq_column_name_1)
        df_csv_1_rv_in_2 = df_csv_1_rv[df_csv_1_rv.iloc[:, -1].isin(df_csv_2.iloc[:, -1])].sort_values(seq_column_name_1)
        df_csv_2_g_in_1_g = df_csv_2[df_csv_2.iloc[:, -1].isin(df_csv_1_g_in_2.iloc[:, -1])].sort_values(seq_column_name_2)
        df_csv_2_g_in_1_rv = df_csv_2[df_csv_2.iloc[:, -1].isin(df_csv_1_rv_in_2.iloc[:, -1])].sort_values(seq_column_name_2)
        
        sub_csv_1_path = self.csv_1.create_sub_csv()
        sub_csv_2_path = self.csv_2.create_sub_csv()

        df_csv_1_g_in_2.iloc[: , :-1].to_csv(sub_csv_1_path, mode='w', index=False)
        df_csv_1_rv_in_2.iloc[: , :-1].to_csv(sub_csv_1_path, mode='a', header=False, index=False)
        df_csv_2_g_in_1_g.iloc[: , :-1].to_csv(sub_csv_2_path, mode='w', index=False)
        df_csv_2_g_in_1_rv.iloc[: , :-1].to_csv(sub_csv_2_path, mode='a', header=False, index=False)