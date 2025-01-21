import ast
import os
from pathlib import Path
import gzip
import numpy as np
import typer
from Bio.Seq import Seq
from Bio import SeqIO
import re
from dataclasses import dataclass
from typing import List, Optional, Literal


# ============================
# Constants and Configurations
# ============================
class Constants:
    SPELLINGS = {
        "shuffle": ["shuffle", "sh", "s"],
        "reverse_complement": ["reverse_complement", "rv", "rc", "r"],
        "both": ["both", "b", "rs", "rcs", "rcsh", "shrc", "sr"],
    }
    SINGLE_PATTERN = re.compile(r"(?P<base>^.+)(?P<ext>\.(?:fq|fastq)(?:.gz)?)$")
    PAIRED_PATTERN = re.compile(
        r"(?P<base>^.+)(?P<read>_(?:[rR]?[12]))(?P<ext>\.(?:fq|fastq)(?:.gz)?)$"
    )


# ============================
# Data Structures
# ============================
@dataclass
class ParsedFastqFileName:
    base: str
    read: List[str]
    ext: List[str]
    pair_type: Literal["single", "paired"]


# ============================
# Utility Classes
# ============================
class RecordTransformer:
    @staticmethod
    def shuffle(records: List[Seq], base_seed: int, num_shuffles: int) -> List[List[int]]:
        """
        Generate a list of unique shuffle orders for a set of records.
        """
        orders = []
        original_order = list(range(len(records)))
        scale = max(1000, int(1e6 / len(records)))

        for i in range(num_shuffles):
            rng = np.random.default_rng(base_seed + i * scale)
            while True:
                shuffle_order = rng.permutation(len(records))
                if not np.array_equal(shuffle_order, original_order):
                    break
            orders.append(shuffle_order.tolist())
        return orders

    @staticmethod
    def reverse_complement(records: List[Seq]) -> List[Seq]:
        """
        Create reverse complement of FASTQ records.
        """
        rc_records = []
        for rec in records:
            rc_seq = str(rec.seq.reverse_complement())
            rc_record = rec.__class__(
                id=rec.id,
                seq=Seq(rc_seq),
                description=rec.description,
                letter_annotations={"phred_quality": rec.letter_annotations["phred_quality"][::-1]}
            )
            rc_records.append(rc_record)
        return rc_records


# ============================
# Core Classes
# ============================
class InputFileParser:
    def __init__(self, pair_type: str, rep_type: str, rep_num: int, fastq1: str, fastq2: Optional[str] = None):
        self.pair_type = pair_type
        self.rep_type = rep_type  # Allow rep_type to be empty initially
        self.rep_num = rep_num
        self.fastq1 = Path(fastq1)
        self.fastq2 = Path(fastq2) if fastq2 else None

        # Parse and store file names immediately
        self.parsed_names = self._parse_file_names()

    def validate_rep_type(self):
        """Validate the replicate type."""
        self.rep_type = next(
            (k for k, v in Constants.SPELLINGS.items() if self.rep_type.lower() in v), None
        )
        if self.rep_type is None:
            raise ValueError(f"Invalid replicate type: {self.rep_type}")

    def _parse_file_names(self) -> ParsedFastqFileName:
        files = [self.fastq1, self.fastq2] if self.fastq2 else [self.fastq1]

        if len(files) == 1:
            match = Constants.SINGLE_PATTERN.search(str(files[0].name))
            if not match:
                raise ValueError(f"Invalid single-end file name: {files[0].name}")

            return ParsedFastqFileName(
                base=match.group("base" ),
                read=[],
                ext=[match.group("ext")],
                pair_type="single",
            )

        else:
            matches = [Constants.PAIRED_PATTERN.search(str(f.name)) for f in files]
            if not all(matches):
                raise ValueError(f"Invalid paired-end file names: {files}")

            base1, base2 = [m.group("base") for m in matches]
            if base1 != base2:
                raise ValueError(f"Paired-end file names do not match: {base1} vs {base2}")

            return ParsedFastqFileName(
                base=base1,
                read=[m.group("read") for m in matches],
                ext=[m.group("ext") for m in matches],
                pair_type="paired",
            )



class FastqFileProcessor:
    def __init__(self, parser: InputFileParser, seed: int, output_folder: Optional[Path]):
        self.parser = parser
        self.seed = seed
        self.output_folder = output_folder or parser.fastq1.parent

        if self.output_folder:
            seed_subfolder = self.output_folder / f"seed_{self.seed}"
            if not os.path.isdir(seed_subfolder):
                os.mkdir(seed_subfolder)
            self.output_folder = seed_subfolder

    def process(self, replicate_types: List[str]):
        for rtype in replicate_types:
            self.parser.rep_type = rtype
            if self.parser.parsed_names.pair_type == "single":
                self._process_single_end()
            else:
                self._process_paired_end()

    def _process_single_end(self):
        records = self._read_fastq(self.parser.fastq1)
        if self.parser.rep_type == "shuffle":
            self._shuffle_and_write(records)
        elif self.parser.rep_type == "reverse_complement":
            self._reverse_complement_and_write(records)
        elif self.parser.rep_type == "both":
            self._shuffle_reverse_complement_and_write(records)

    def _process_paired_end(self):
        records1 = self._read_fastq(self.parser.fastq1)
        records2 = self._read_fastq(self.parser.fastq2)
        if len(records1) != len(records2):
            raise ValueError("Paired-end files have unequal records.")

        if self.parser.rep_type == "shuffle":
            self._shuffle_and_write(records1, records2)
        elif self.parser.rep_type == "reverse_complement":
            self._reverse_complement_and_write(records1, records2)
        elif self.parser.rep_type == "both":
            self._shuffle_reverse_complement_and_write(records1, records2)

    def _read_fastq(self, file_path: Path) -> List[Seq]:
        with gzip.open(file_path, 'rt') if file_path.suffix == ".gz" else open(file_path, 'r') as handle:
            return list(SeqIO.parse(handle, "fastq"))

    def _shuffle_and_write(self, records1: List[Seq], records2: Optional[List[Seq]] = None):
        num_shuffles = min(
            self.parser.rep_num, np.math.factorial(len(records1)) - 1
        )
        orders = RecordTransformer.shuffle(records1, self.seed, num_shuffles)

        read_suffix1 = self.parser.parsed_names.read[0] if self.parser.parsed_names.read else ""

        for idx, order in enumerate(orders):
            shuffled1 = [records1[i] for i in order]
            self._write_fastq(shuffled1, suffix=f"_sh{idx + 1}", read_suffix=read_suffix1)

            if records2:
                read_suffix2 = self.parser.parsed_names.read[1]
                shuffled2 = [records2[i] for i in order]
                self._write_fastq(shuffled2, suffix=f"_sh{idx + 1}", read_suffix=read_suffix2)

    def _reverse_complement_and_write(self, records1: List[Seq], records2: Optional[List[Seq]] = None):
        rc1 = RecordTransformer.reverse_complement(records1)
        read_suffix1 = self.parser.parsed_names.read[0] if self.parser.parsed_names.read else ""
        self._write_fastq(rc1, suffix="_rc", read_suffix=read_suffix1)

        if records2:
            rc2 = RecordTransformer.reverse_complement(records2)
            read_suffix2 = self.parser.parsed_names.read[1]
            self._write_fastq(rc2, suffix="_rc", read_suffix=read_suffix2)

    def _shuffle_reverse_complement_and_write(self, records1: List[Seq], records2: Optional[List[Seq]] = None):
        num_shuffles = min(
            self.parser.rep_num, np.math.factorial(len(records1)) - 1
        )
        orders = RecordTransformer.shuffle(records1, self.seed, num_shuffles)

        read_suffix1 = self.parser.parsed_names.read[0] if self.parser.parsed_names.read else ""

        for idx, order in enumerate(orders):
            shuffled1 = [records1[i] for i in order]
            rc_shuffled1 = RecordTransformer.reverse_complement(shuffled1)
            self._write_fastq(rc_shuffled1, suffix=f"_both{idx + 1}", read_suffix=read_suffix1)

            if records2:
                read_suffix2 = self.parser.parsed_names.read[1]
                shuffled2 = [records2[i] for i in order]
                rc_shuffled2 = RecordTransformer.reverse_complement(shuffled2)
                self._write_fastq(rc_shuffled2, suffix=f"_both{idx + 1}", read_suffix=read_suffix2)

    def _write_fastq(self, records: List[Seq], suffix: str, read_suffix: str):
        """
        Writes the FASTQ records to an output file.

        Parameters:
            records (List[Seq]): The FASTQ records to write.
            suffix (str): The suffix for the output file (e.g., `_sh1`, `_rc`).
            read_suffix (str): The read suffix for paired-end files (e.g., `_r1`, `_1`).
        """

        output_file = self.output_folder / f"{self.parser.parsed_names.base}{suffix}{read_suffix}{self.parser.parsed_names.ext[0]}"
        with gzip.open(output_file, 'wt') if output_file.suffix == ".gz" else open(output_file, 'w') as handle:
            SeqIO.write(records, handle, "fastq")


# ============================
# Command-Line Interface
# ============================
app = typer.Typer()


@app.command()
def main(
        fastq_file1: str = typer.Option(..., "--fastq_file1", "-f1", help="First FASTQ file"),
        fastq_file2: Optional[str] = typer.Option(None, "--fastq_file2", "-f2",
                                                  help="Second FASTQ file (optional for paired data)"),
        output_folder: Optional[Path] = typer.Option(None, "--output_folder", "-o", help="Output folder"),
        replicate_types: Optional[List[str]] = typer.Option(None, "--rep_types", "-r",
                                                            help="List of replicate types (e.g., shuffle, rc, both)"),
        replicate_number: int = typer.Option(1, "--rep_num", "-n", help="Number of replicates"),
        seed: int = typer.Option(1, "--seed", "-s", help="Seed for reproducibility"),
        pair_type: str = typer.Option("paired", "--pair_type", "-p", help="Pair type of data"),
        all_replicates: bool = typer.Option(False, "--all", "-a", help="Process all replicate types"),
):
    # Default to all replicate types if `all_replicates` is True
    if all_replicates:
        replicate_types = ["shuffle", "reverse_complement", "both"]

    # Ensure `replicate_types` is not None
    if replicate_types is None:
        raise ValueError("You must provide replicate types with --rep_types or use --all.")

    # Validate replicate types
    replicate_types = [rtype.lower() for rtype in replicate_types]
    valid_types = Constants.SPELLINGS.keys()
    for rtype in replicate_types:
        if rtype not in valid_types:
            print('not valid: ', rtype)
            raise ValueError(f"Invalid replicate type: {rtype}. Allowed types: {list(valid_types)}")

    # Process each replicate type
    parser = InputFileParser(pair_type, "", replicate_number, fastq_file1, fastq_file2)
    processor = FastqFileProcessor(parser, seed, output_folder)
    processor.process(replicate_types)


if __name__ == "__main__":
    app()
