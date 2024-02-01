import csv
import math
import os
from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Union, List, Tuple
import chardet
import numpy as np
import typer
import re
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
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


class InputFileParser:
    def __init__(self, fastq_type: str,
                 rep_type: str,
                 rep_num: int,
                 fastq1: str,
                 fastq2: str = None,
                 ):

        # self.directory_path = Path(directory_path)
        self.fastq_type = fastq_type
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.rep_type = rep_type
        self.rep_num = rep_num
        # self.file_tuple = ()

        self.fastq1 = os.path.basename(fastq1)
        self.fastq2 = os.path.basename(fastq2)
        self.fastq_path = Path(os.path.dirname(fastq1))

        # print('-----------------------------------')
        # print(f'directory: {self.fastq_path}')
        # print(f'file1: {self.fastq1}')
        # print(f'file2: {self.fastq2}')
        # print('-----------------------------------')

        file_names, parsed_names = self._check_fastq_files(self.fastq_path,
                                                           self.fastq_type,
                                                           self.fastq1,
                                                           self.fastq2)

        # print(file_names, parsed_names)
        # Normalize type spelling to one of 3 keys in spellings dict
        self.rep_type = [k for k, v in spellings.items() \
                         if self.rep_type.lower() in v][0]

        self.file_tuple = InfoTuple(self.fastq_path,
                                    file_names,
                                    parsed_names,
                                    self.fastq_type,
                                    self.rep_type,
                                    self.rep_num
                                    )

    @staticmethod
    def _check_fastq_files(
            file_path: Path,
            fastq_type: str,
            file1: str,
            file2: str = None) -> tuple[list[str], ParsedFastqFileName]:

        # print(file1)
        # print(file_path)
        # print(f"{file_path} / {file1}")

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

    def __init__(self, parser: InputFileParser,
                 seed: int,
                 out_folder: Path
                 ):
        self.parser = parser
        self.info_tuple = parser.file_tuple
        self.seed = seed
        self.output_folder = out_folder
        self.filenames = []

        if out_folder and not os.path.isdir(out_folder):
            os.mkdir(out_folder)

    def fastq_replicator(self):
        self._process_tuple(self.info_tuple)

    def _process_tuple(self, info):

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
        def _limit_number_of_replicates(info) -> int:
            number = min(
                info.replicate_number,
                math.factorial(len(records)) - 1
            )
            return number

        def shuffle(records: list, orders: list, keep_orders: bool = False):
            np.random.seed(self.seed)
            for i in range(_limit_number_of_replicates(info)):
                shuffled_records = self._shuffle_records(
                    records, orders, keep_orders, i
                )
                self._write_records(shuffled_records, info, i, int(keep_orders))

        def rev_comp(records: list, write_to_file: bool = True, read: int = 0):
            rc_records = self._reverse_complement_records(records)
            if write_to_file:
                self._write_records(rc_records, info, idx=0, read=read)
            return rc_records

        # Initialize an empty list to hold the orders of shuffling
        generated_orders = []
        # Initialize an empty list to hold FASTQ records
        records = []

        for idx, file_name in enumerate(info.file_names):
            # print(idx, file_name)
            file_path = os.path.join(info.directory_path, file_name)

            # If the pair type is single, process the single-end file
            if info.pair_type == 'single':
                records = _get_records(file_path)

                if info.replicate_type == 'shuffle':
                    shuffle(records, generated_orders)

                elif info.replicate_type == "reverse_complement":
                    rev_comp(records, write_to_file=True)

                elif info.replicate_type == 'both':
                    rc_records = rev_comp(records, write_to_file=False)
                    shuffle(rc_records, generated_orders)

            # Process two paired-end files, one per iteration
            elif info.pair_type == 'paired':
                records = _get_records(file_path)

                if info.replicate_type == 'shuffle':
                    # If it's the first reads in pair shuffle anyway
                    if idx == 0:
                        shuffle(records, generated_orders)
                    # If it's the second reads in pair, use already
                    # generated orders for the first reads in pair
                    elif idx == 1:
                        shuffle(records, generated_orders, keep_orders=True)

                elif info.replicate_type == "reverse_complement":
                    rev_comp(records, write_to_file=True, read=idx)

                elif info.replicate_type == 'both':
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
            self.filenames.append(os.path.join(self.output_folder, filename))
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


@app.command()
def main(

        # input_folder: Path = typer.Option(
        #     ...,
        #     "--input_folder", "-i",
        #     file_okay=False,
        #     dir_okay=True,
        #     writable=True,
        #     readable=True,
        #     resolve_path=True,
        #     show_default=False,
        #     help="Path to the input FASTQ files"
        # ),

        fastq_file1: str = typer.Option(
            ...,
            "--fastq_file1", "-f1",
            help="Name of the first FASTQ file for paired data"
        ),

        fastq_file2: str = typer.Option(
            ...,
            "--fastq_file2", "-f2",
            help="Name of to the second FASTQ file for paired data"
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
            help="Path to the output folder of replicated FASTQ files"
        ),

        replicate_type: str = typer.Option(
            "shuffle",
            "--rep_type", "-r",
            help="The type of replicates"
            # callback=validate_replicate_type,
        ),

        replicate_number: int = typer.Option(
            1,
            "--rep_num", "-n",
            help="Number of replicates"
        ),

        seed: int = typer.Option(
            1,
            "--seed", "-s",
            help="Seed number for reproducibility"
        ),

        pair_type: str = typer.Option(
            "paired",
            "--pair_type", "-p",
            help="Pair type of the data"
        )

):
    if pair_type == "paired":
        if not (fastq_file1 and fastq_file2):
            typer.echo("Error: Both fastq_file1 and fastq_file2 are required for paired data.")
            raise typer.Abort()
    elif pair_type == "single":
        if not fastq_file1:
            typer.echo("Error: fastq_file1 is required for single data.")
            raise typer.Abort()
    else:
        typer.echo("Error: Invalid pair_type. Must be 'paired' or 'single'.")
        raise typer.Abort()

    parser = InputFileParser(pair_type,
                             replicate_type,
                             replicate_number,
                             fastq_file1,
                             fastq_file2 if pair_type == "paired" else None
                             )

    replicator = FastqFileReplicator(parser,
                                     seed,
                                     output_folder)
    replicator.fastq_replicator()


if __name__ == "__main__":
    app()
