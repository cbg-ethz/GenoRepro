import csv
import os
from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Union, List
import chardet
import numpy as np
import typer
import re
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import gzip

# The namedTuple represents a row of data in the input csv file
FileRow = namedtuple("FileRow",
                     ["directory_path",
                      "file_names",
                      "parsed_names",
                      "pair_type",
                      "replicate_type",
                      "replicate_number"],
                     )

# All the allowed spellings of replicate types in the input csv file
spellings = {
    "shuffle": ["shuffle", "sh", "s"],
    "reverse_complement": ["reverse_complement", "rv", "rc", "r"],
    "both": ["both", "b", "rs", "rcs", "rcsh", "shrc", "sr"]
}

# List of column names in the input header
header_row = [
    "directory_path",
    "read1",
    "read2",
    "replicate_type",
    "replicate_number",
]


@dataclass
class ParsedFastqFileName:
    """
    ParsedFastqFileName stores parsed by Regex groups of filename for each row
    to build output filepath
    """
    base: str
    read: List[str]
    ext: List[str]
    pair_type: Literal['single', 'paired']


class CSVFileParser:
    """
    CSVFileParser parse the input csv_file which stores 
    information for the fastq files
    """

    def __init__(self, input_file: Union[str, bytes, os.PathLike]):
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
            for idx, row in enumerate(reader):
                # Skip header row if it is present in the input CSV file
                if idx == 0 and \
                        row == \
                        list(map(lambda x: x.lower().strip(), header_row)):
                    continue

                # Extract the directory path from the first column
                directory_path = Path(row[0])

                # Extract the file names from the 2nd and 3rd columns
                file_names = [name.strip() for name in row[1:3] \
                              if name.strip() != '']

                # Throw an error if no files are provided
                if len(file_names) < 1:
                    raise ValueError(
                        f"No FASTQ files were provided "
                        f"in CSV file row number {idx + 1}"
                    )

                # Determine the paired-end type based on the file names
                parsed_names = self._check_fastq_files(
                    directory_path, file_names
                )
                pair_type = parsed_names.pair_type

                # Set shuffling by default if replicate type not provided
                replicate_type = row[3].lower() \
                    if len(row) >= 4 and row[3] else "shuffle"

                # Normalize type spelling to one of 3 keys in spellings dict
                replicate_type = [k for k, v in spellings.items() \
                                  if replicate_type.lower() in v][0]

                # Check for replicate type allowed spelling from the dict
                if replicate_type not in sum(spellings.values(), []):
                    raise ValueError(
                        f"Invalid replicate type '{replicate_type}' "
                        f"in CSV file row number {idx + 1}"
                    )

                # Extract the replicate number or set it to 1 if not provided
                replicate_number = int(row[4]) \
                    if len(row) >= 4 and row[4] else 1

                row_tuple = FileRow(directory_path,
                                    file_names,
                                    parsed_names,
                                    pair_type,
                                    replicate_type,
                                    replicate_number)

                self.rows.append(row_tuple)

    @staticmethod
    def _check_fastq_files(
            file_path: Path,
            fastq_files: List[str]) -> ParsedFastqFileName:

        missing_files = [f for f in fastq_files if
                         not (file_path / f).is_file()
                         or not f.endswith(("fastq", "fq", "gz"))]

        if missing_files:
            raise FileNotFoundError(
                f"FASTQ files: {', '.join(f for f in missing_files)} - "
                f"not found in the direcotry {file_path}"
            )

        # Check the file names to determine whether they are 
        # single-end or paired-end reads
        if len(fastq_files) == 1:
            # Regex pattern to determine the file name pattern of the file
            pattern = (
                # Sample basename (i.e. 'ERR001010')
                r'(?P<base>^.+)'
                # Extension of the fastq file (i.e. '.fq.gz')
                r'(?P<ext>\.(?:fq|fastq)(?:.gz)?)$'
            )
            base, ext = re.search(pattern, fastq_files[0]).groups()
            parsed_names = ParsedFastqFileName(
                base=base,
                read=[],
                ext=[ext],
                pair_type='single'
            )

        else:
            # Regex pattern to determine the filename pattern of the files
            pattern = (
                # Sample basename (i.e. 'ERR001010')
                r'(?P<base>^.+)'
                # Number of the read in the pair (i.e. '_2')
                r'(?P<read>_(?:r|R)?[12])'
                # Extension of the fastq file (i.e. '.fq.gz')
                r'(?P<ext>\.(?:fq|fastq)(?:.gz)?)$'
            )

            # Check if file names match pattern
            if all(re.match(pattern, f) for f in fastq_files):
                matches = list(map(lambda x: re.search(pattern, x), \
                                   fastq_files))
                base1, base2 = map(lambda x: x.group('base'), matches)
                read1, read2 = map(lambda x: x.group('read'), matches)
                ext1, ext2 = map(lambda x: x.group('ext'), matches)

                if base1 != base2:
                    raise ValueError(
                        f"{fastq_files[0]} and {fastq_files[1]} do not match"
                    )

                parsed_names = ParsedFastqFileName(
                    base=base1,
                    read=[read1, read2],
                    ext=[ext1, ext2],
                    pair_type='paired'
                )
            else:
                raise ValueError(
                    f"Paired-end files {fastq_files[0]} and {fastq_files[1]} "
                    f"in {file_path} do not follow the expected "
                    f"naming convention '_R1' and '_R2'"
                )
        return parsed_names


class FastqFileReplicator:
    """
    FastqFileReplicator takes a CSVFileParser argument and replicate 
    FASTQ files based on the information from this parser
    """

    def __init__(self, parser: CSVFileParser, seed: int, output_folder: Path):
        # Initialize the class with a CSVFileParser instance 
        # and a seed for random number generation
        self.parser = parser
        self.rows = parser.rows
        self.seed = seed
        self.output_folder = output_folder
        if output_folder and not os.path.isdir(output_folder):
            os.mkdir(output_folder)

    def fastq_replicator(self, cores: int):
        # Create a multiprocessing pool with a number of cores
        pool = mp.Pool(cores)
        # Use the pool to process each row of the CSV file 
        # by calling the _process_row method
        pool.map(self._process_row, self.rows)
        # Close the pool to prevent further tasks from being submitted to it
        pool.close()
        # Wait for the worker processes to exit
        pool.join()

    def _process_row(self, row):

        # Get list of records
        def _get_records(file_path) -> list:
            if file_path.endswith(".gz"):
                with gzip.open(file_path, 'rt') as handle:
                    return list(SeqIO.parse(handle, "fastq"))
            else:
                with open(file_path, 'r') as handle1:
                    return list(SeqIO.parse(handle1, "fastq"))

        # If the replicate_number is greater than the maximum possible shuffle
        # then assign the max possible shuffle to rep_number
        def _limit_number_of_replicates(row) -> int:
            number = min(
                row.replicate_number,
                np.math.factorial(len(records)) - 1
            )
            return number

        def shuffle(records: list, orders: list, keep_orders: bool = False):
            np.random.seed(self.seed)
            for i in range(_limit_number_of_replicates(row)):
                shuffled_records = self._shuffle_records(
                    records, orders, keep_orders, i
                )
                self._write_records(shuffled_records, row, i, int(keep_orders))

        def rev_comp(records: list, write_to_file: bool = True, read: int = 0):
            rc_records = self._reverse_complement_records(records)
            if write_to_file:
                self._write_records(rc_records, row, idx=0, read=read)
            return rc_records

        # Initialize an empty list to hold the orders of shuffling
        generated_orders = []
        # Initialize an empty list to hold FASTQ records
        records = []
        # Process a single row of the CSV file by taking it as input
        for idx, file_name in enumerate(row.file_names):
            file_path = os.path.join(row.directory_path, file_name)

            # If the pair type is single, process the single-end file
            if row.pair_type == 'single':
                records = _get_records(file_path)

                if row.replicate_type == 'shuffle':
                    shuffle(records, generated_orders)

                elif row.replicate_type == "reverse_complement":
                    rev_comp(records, write_to_file=True)

                elif row.replicate_type == 'both':
                    rc_records = rev_comp(records, write_to_file=False)
                    shuffle(rc_records, generated_orders)

            # Process two paired-end files, one per iteration
            elif row.pair_type == 'paired':
                records = _get_records(file_path)

                if row.replicate_type == 'shuffle':
                    # If it's the first reads in pair shuffle anyway
                    if idx == 0:
                        shuffle(records, generated_orders)
                    # If it's the second reads in pair, use already
                    # generated orders for the first reads in pair
                    elif idx == 1:
                        shuffle(records, generated_orders, keep_orders=True)

                elif row.replicate_type == "reverse_complement":
                    rev_comp(records, write_to_file=True, read=idx)

                elif row.replicate_type == 'both':
                    rc_records = rev_comp(records, write_to_file=False)
                    # If it's the first list of records in pair shuffle anyway
                    if idx == 0:
                        shuffle(rc_records, generated_orders)
                    # If it's the second list of records in pair, use already
                    # generated orders for the first list of records
                    elif idx == 1:
                        shuffle(rc_records, generated_orders, keep_orders=True)

    @staticmethod
    def _shuffle_records(
            records: list,
            generated_orders: list,
            keep_orders: bool = False,
            i: int = -1
    ):
        """Shuffle a list of records using 
        the Fisher-Yates shuffle algorithm."""
        if keep_orders:
            for k in range(len(records)):
                yield records[generated_orders[i][k]]
        else:
            original_order = list(range(len(records)))
            new_order = original_order.copy()
            while new_order == original_order or \
                    new_order in generated_orders:
                np.random.shuffle(new_order)

            generated_orders.append(new_order)

            for k in range(len(records)):
                yield records[new_order[k]]

    @staticmethod
    def _reverse_complement_records(records):
        rev_record = []
        for rec in records:
            seq = str(rec.seq.reverse_complement())
            seq_qual = rec.letter_annotations["phred_quality"][::-1]
            rev_record.append(
                rec.__class__(
                    id=rec.id,
                    seq=Seq(seq),
                    description=rec.description,
                    letter_annotations={"phred_quality": seq_qual}
                )
            )
        return rev_record

    def _get_output_file_path(self, row, read=0, idx=0):
        def build_name(rep_type: str) -> Path:
            filename = f"{row.parsed_names.base}_{rep_type}"
            if rep_type != "rc":
                filename += f"{idx + 1}"
            if row.parsed_names.pair_type == "paired":
                filename += f"{row.parsed_names.read[read]}"
            filename += f"{row.parsed_names.ext[read]}"
            if not self.output_folder:
                self.output_folder = row.directory_path
            return os.path.join(self.output_folder, filename)

        if row.replicate_type == 'shuffle':
            return build_name("sh")
        elif row.replicate_type == 'reverse_complement':
            return build_name("rc")
        elif row.replicate_type == 'both':
            return build_name(row.replicate_type)

    def _write_records(self, recs, row, idx=0, read=0):
        output_file = self._get_output_file_path(row=row, idx=idx, read=read)
        if output_file.endswith("gz"):
            with gzip.open(output_file, 'wt') as f:
                SeqIO.write(recs, f, "fastq")
        else:
            with open(output_file, 'wt') as f:
                SeqIO.write(recs, f, "fastq")


# Define the command line interface using Typer
app = typer.Typer(add_completion=False)


@app.command()
def main(

        csv_file: Path = typer.Option(
            ...,
            "--csv_file", "-i",
            exists=True,
            file_okay=True,
            dir_okay=False,
            show_default=False,
            help="Path to the CSV file containing FASTQ files information"
        ),

        output_folder: Path = typer.Option(
            None,
            "--output", "-o",
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=False,
            help="Path to the output folder of replicated FASTQ file"
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
    replicator = FastqFileReplicator(parser, seed, output_folder)
    replicator.fastq_replicator(cores=cores)


if __name__ == "__main__":
    app()
