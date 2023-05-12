import csv
import os
from collections import namedtuple
from pathlib import Path
import chardet
import numpy as np
import typer
import re

from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp

""" The namedTuple represents a row of data in the input csv file """

FileRow = namedtuple("FileRow",
                     ["directory_path",
                      "file_names",
                      "pair_type",
                      "replicate_type",
                      "replicate_number"]
                     )

""" CSVFileParser parse the input csv_file which stores information for the fastq files"""


class CSVFileParser:
    def __init__(self, input_file):
        # Initialize an empty list to store the parsed rows
        self.rows = []
        self.csv_file = input_file

        # Detect the encoding of the input file using the chardet library
        with open(self.csv_file, 'rb') as f:
            result = chardet.detect(f.read())
            encoding = result['encoding']

        # Parse the CSV file using the csv module and the detected encoding
        with open(self.csv_file, encoding=encoding) as f:
            reader = csv.reader(f)
            for row in reader:
                # Extract the directory path from the first column
                directory_path = Path(row[0])

                # Extract the file names from the second column and split by comma for paired-end files
                file_names = [name.strip() for name in row[1].split(",")] if "," in row[1] else [row[1]]
                pair_type = self._check_fastq_files(directory_path, file_names)

                # Determine the paired-end type (single, paired-end 1, paired-end 2) based on the file names
                replicate_type = row[2].lower() if len(row) >= 3 and row[2] else "shuffle"

                if replicate_type not in ["shuffle", "reverse_complement", "both"]:
                    raise ValueError(f"Invalid replicate type {replicate_type} in CSV file row number {row}")

                # Extract the replicate number from the fourth column or set it to 1 if not provided
                replicate_number = int(row[3]) if len(row) >= 4 and row[3] else 1
                if len(file_names) > 2:
                    raise ValueError(f"the number of files cannot exceed 2")

                # Create a tuple of the parsed information and append it to the list of rows
                row_tuple = FileRow(directory_path,
                                    file_names,
                                    pair_type,
                                    replicate_type,
                                    replicate_number)

                self.rows.append(row_tuple)

    @staticmethod
    def _check_fastq_files(file_path, fastq_files):

        missing_files = [f for f in fastq_files if
                         not (file_path / f).is_file() or not f.endswith(("fastq", "fq"))]

        if missing_files:
            raise FileNotFoundError(
                f"FASTQ files not found in {file_path}: {', '.join(os.path.join(file_path, f) for f in missing_files)}")

        # Check the file names to determine whether they are single-end or paired-end reads
        if len(fastq_files) == 1:
            pair_type = 'single'

        else:
            pair_type = 'paired'

            # regex pattern to determine the file name pattern of the fastq files
            pattern = r'^.+(?:_r|R)[12]\.(fq|fastq)$'

            # Check if file names match pattern
            if all([re.match(pattern, fastq_files[0]), re.match(pattern, fastq_files[1])]):
                base1, ext1 = fastq_files[0].split('.')[0].lower().rsplit('_r', 1)
                base2, ext2 = fastq_files[1].split('.')[0].lower().rsplit('_r', 1)

                if base1 != base2:
                    raise ValueError(f"{fastq_files[0]} and {fastq_files[1]} do not match")
                else:
                    if ext1 != ext2:
                        pair_type = 'paired'
            else:
                raise ValueError(
                    f"Paired-end files {fastq_files[0]} and {fastq_files[1]} in {file_path} "
                    f"do not follow the expected naming convention '_R1' and '_R2'")
        return pair_type


"""FastqFileReplicator takes a CSVFileParser argument and replicate 
FASTQ files based on the information from this parser"""


class FastqFileReplicator:
    def __init__(self, parser: CSVFileParser, seed):
        # Initialize the class with a CSVFileParser instance and a seed for random number generation
        self.parser = parser
        self.rows = parser.rows
        self.records = []  # Initialize an empty list to hold FASTQ records
        self.file_path = []
        self.seed = seed
        self.generated_orders = []  # Initialize an empty list to hold the orders of shuffling

    def fastq_replicator(self, cores):
        # Create a multiprocessing pool with a number of cores
        pool = mp.Pool(cores)
        # Use the pool to process each row of the CSV file by calling the _process_row method
        pool.map(self._process_row, self.rows)
        # Close the pool to prevent further tasks from being submitted to it
        pool.close()
        # Wait for the worker processes to exit
        pool.join()

    def _process_row(self, row):

        # Process a single row of the CSV file by taking it as input
        for file_name in row.file_names:
            self.file_path = os.path.join(row.directory_path, file_name)

            # If the pair type is single, process the single-end file
            if row.pair_type == 'single':
                with open(self.file_path, 'r') as handle1:
                    self.records = list(SeqIO.parse(handle1, "fastq"))

                if row.replicate_type == 'shuffle':
                    np.random.seed(self.seed)

                    # If the replicate_number is greater than the maximum possible shuffle
                    # than assign the max possible shuffle to rep_number
                    rep_number = min(row.replicate_number, np.math.factorial(len(self.records)) - 1)

                    for i in range(rep_number):
                        shuffled_records = self._shuffle_records(self.records, self.generated_orders)
                        self._write_records(shuffled_records, row, i)

                elif row.replicate_type == "reverse_complement":
                    self.records = self._reverse_complement_records(self.records)
                    self._write_records(self.records, row)

                elif row.replicate_type == 'both':

                    np.random.seed(self.seed)

                    # If the replicate_number is greater than the maximum possible shuffle
                    # than assign the max possible shuffle to rep_number

                    rep_number = min(row.replicate_number, np.math.factorial(len(self.records)) - 1)
                    rc_records = self._reverse_complement_records(self.records)

                    for i in range(rep_number):
                        shuffled_records = self._shuffle_records(rc_records, self.generated_orders)
                        self._write_records(shuffled_records, row, i)

            elif row.pair_type == 'paired':
                pass

    @staticmethod
    def _shuffle_records(records, generated_orders):

        """Shuffle a list of records using the Fisher-Yates shuffle algorithm."""
        original_order = list(range(len(records)))
        new_order = original_order.copy()
        while new_order == original_order or new_order in generated_orders:
            np.random.shuffle(new_order)

        generated_orders.append(new_order)

        for i in range(len(records)):
            yield records[new_order[i]]

    @staticmethod
    def _reverse_complement_records(records):
        rev_record = []
        for rec in records:
            seq = str(rec.seq.reverse_complement())
            seq_qual = rec.letter_annotations["phred_quality"][::-1]
            rev_record.append(rec.__class__(id=rec.id,
                                            seq=Seq(seq),
                                            description=rec.description,
                                            letter_annotations={"phred_quality": seq_qual}))
        return rev_record

    @staticmethod
    def _get_output_file_path(row, idx=0):
        file_name = os.path.basename(row.file_names[0])
        output_file = os.path.join(row.directory_path, file_name)

        if row.pair_type == 'single':
            if row.replicate_type == 'shuffle':
                rep_type = 'sh'
                output_file = output_file.replace('.fastq', f'_{rep_type}{idx + 1}.fastq')

            elif row.replicate_type == 'reverse_complement':
                rep_type = 'rc'
                output_file = output_file.replace('.fastq', f'_{rep_type}.fastq')

            elif row.replicate_type == 'both':
                output_file = output_file.replace('.fastq', f'_{row.replicate_type}{idx + 1}.fastq')

        if row.pair_type == 'paired':
            pass

        return output_file

    def _write_records(self, recs, row, idx=0):

        output_file = self._get_output_file_path(row=row, idx=idx)
        with open(output_file, 'wt') as f:
            SeqIO.write(recs, f, "fastq")


# Define the command line interface using Typer
app = typer.Typer(add_completion=False)


@app.command()
def main(csv_file: str = typer.Option(
    ...,
    "--csv_file", "-i",
    exists=True,
    file_okay=True,
    dir_okay=False,
    show_default=False,
    help="Path to the tsv file containing FASTQ files information"
),
        output_folder: str = typer.Option(
            None,
            "--output", "-o",
            exists=True,
            help="Path to the output of replicated FASTQ file"
        ),
        seed: int = typer.Option(
            1,
            "--seed", "-s",
            help="Seed number for reproducibility"
        ),

        cores: int = typer.Option(
            1,
            "--cores", "-c",
            help="Number of cores to use"

        )
):
    parser = CSVFileParser(csv_file)
    replicator = FastqFileReplicator(parser, seed)
    replicator.fastq_replicator(cores=cores)


if __name__ == "__main__":
    app()
