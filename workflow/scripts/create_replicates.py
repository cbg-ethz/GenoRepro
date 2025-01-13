import base64
import hashlib
import math
import os
import time
from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, List
import numpy as np
import typer
import re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

# The namedTuple represents a row of data in the input csv file
InfoTuple = namedtuple("InfoTuple",
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

# Define the command line interface using Typer
app = typer.Typer(add_completion=False)

def generate_shuffles(records, base_seed, num_shuffles):
    """
    Generate a list of unique shuffle orders for a set of records.
    Dynamically adjusts the seed offset based on the number of reads.
    Excludes the original order.
    """
    orders = []
    original_order = list(range(len(records)))
    num_reads = len(records)
    scale = max(1000, int(1e6 / num_reads))  # Larger scale for smaller datasets

    for i in range(num_shuffles):
        rng = np.random.default_rng(base_seed + i * scale)
        while True:
            shuffle_order = rng.permutation(len(records))
            if not np.array_equal(shuffle_order, original_order):
                break
        orders.append(shuffle_order)
    return orders

class InputFileParser:
    def __init__(self, fastq_type: str,
                 rep_type: str,
                 rep_num: int,
                 fastq1: str,
                 fastq2: str = None):

        self.fastq_type = fastq_type
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.rep_type = rep_type
        self.rep_num = rep_num

        self.fastq1 = os.path.basename(fastq1)
        if fastq2 is not None:
            self.fastq2 = os.path.basename(fastq2)
        self.fastq_path = Path(os.path.dirname(fastq1))

        file_names, parsed_names = self._check_fastq_files(self.fastq_path,
                                                           self.fastq_type,
                                                           self.fastq1,
                                                           self.fastq2)

        self.rep_type = [k for k, v in spellings.items() if self.rep_type.lower() in v][0]

        self.file_tuple = InfoTuple(self.fastq_path,
                                    file_names,
                                    parsed_names,
                                    self.fastq_type,
                                    self.rep_type,
                                    self.rep_num)

        print('-----------------------------------')
        print(f'directory: {self.fastq_path}')
        print(f'file1: {self.fastq1}')
        print(f'file2: {self.fastq2}')
        print(f'fastq_type: {self.fastq_type}')
        print(f'rep_type: {self.rep_type}')
        print(f'rep_num: {self.rep_num}')
        print('-----------------------------------')


    @staticmethod
    def _check_fastq_files(
            file_path: Path,
            fastq_type: str,
            file1: str,
            file2: str = None) -> tuple[list[str], ParsedFastqFileName]:

        if file2 is not None:
            if not (file_path / file1).is_file() or not (file_path / file2).is_file() or not file1.endswith(
                    ("fastq", "fq", "gz")) or not file2.endswith(("fastq", "fq", "gz")):
                raise FileNotFoundError(f"One or both of the provided files do not exist "
                                        f"or have an invalid extension.")
            else:
                fastq_files = [file1, file2]
        else:
            if not (file_path / file1).is_file() or not file1.endswith(("fastq", "fq", "gz")):
                raise FileNotFoundError(f"The provided file does not exist or "
                                        f"has an invalid extension.")
            else:
                fastq_files = [file1]

        if fastq_type == 'single':
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
                matches = list(map(lambda x: re.search(pattern, x), fastq_files))
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
        return fastq_files, parsed_names

class FastqFileReplicator:

    def __init__(self, parser: InputFileParser, seed: int, out_folder: Path):
        self.parser = parser
        self.info_tuple = parser.file_tuple
        self.seed = seed
        self.output_folder = out_folder

        if out_folder and not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        self.output_filepaths = []

        # Create a subdirectory based on the seed number if it doesn't exist
        if out_folder:
            seed_subfolder = out_folder / f"seed_{self.seed}"
            print(seed_subfolder)
            if not os.path.isdir(seed_subfolder):
                os.mkdir(seed_subfolder)
            self.output_folder = seed_subfolder

        # print(out_folder)

        self.output_filepaths = []

    def fastq_replicator(self):
        self._process_tuple(self.info_tuple)

    def _process_tuple(self, info):

        def _get_records(file_path):
            if file_path.endswith(".gz"):
                with gzip.open(file_path, 'rt') as handle:
                    return list(SeqIO.parse(handle, "fastq"))
            else:
                with open(file_path, 'r') as handle:
                    return list(SeqIO.parse(handle, "fastq"))

        records = _get_records(os.path.join(info.directory_path, info.file_names[0]))
        num_shuffles = info.replicate_number

        shuffle_orders = generate_shuffles(records, self.seed, num_shuffles)

        for idx, order in enumerate(shuffle_orders):
            shuffled_records = [records[i] for i in order]
            self._write_records(shuffled_records, info, idx)

    def _write_records(self, recs, row, idx=0):
        output_file = self._get_output_file_path(row, idx)
        if output_file.endswith("gz"):
            with gzip.open(output_file, 'wt') as f:
                SeqIO.write(recs, f, "fastq")
        else:
            with open(output_file, 'wt') as f:
                SeqIO.write(recs, f, "fastq")
        self.output_filepaths.append(output_file)

    def _get_output_file_path(self, row, idx=0):
        filename = f"{row.parsed_names.base}_sh{idx + 1}{row.parsed_names.ext[0]}"
        if not self.output_folder:
            self.output_folder = row.directory_path
        return os.path.join(self.output_folder, filename)

@app.command()
def main(
        fastq_file1: str = typer.Option(..., "--fastq_file1", "-f1", help="First FASTQ file"),
        fastq_file2: str = typer.Option(None, "--fastq_file2", "-f2", help="Second FASTQ file (optional for paired data)"),
        output_folder: Path = typer.Option(None, "--output", "-o", help="Output folder"),
        replicate_type: str = typer.Option("shuffle", "--rep_type", "-r", help="Type of replicates"),
        replicate_number: int = typer.Option(1, "--rep_num", "-n", help="Number of replicates"),
        seed: int = typer.Option(1, "--seed", "-s", help="Seed for reproducibility"),
        pair_type: str = typer.Option("paired", "--pair_type", "-p", help="Pair type of data")):

    parser = InputFileParser(pair_type, replicate_type, replicate_number, fastq_file1, fastq_file2)
    replicator = FastqFileReplicator(parser, seed, output_folder)
    replicator.fastq_replicator()

if __name__ == "__main__":
    app()
