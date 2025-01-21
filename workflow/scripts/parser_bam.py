import csv
from pathlib import Path
import pysam
import typer
import re

class SameNameReadsArray(list):
    def __init__(self, *reads: pysam.AlignedSegment, include_secondary=False):
        self.name = reads[0].query_name if reads else None
        self.include_secondary = include_secondary
        super().__init__(reads)

    def process_reads(self, label: str) -> list:
        primary = [read for read in self if not read.is_secondary and not read.is_supplementary \
                   and not read.is_duplicate]
        primary_pairs = self.pair_reads(primary)

        rows = [self.create_row(pair, label) for pair in primary_pairs]

        if self.include_secondary:
            secondary = [read for read in self if read.is_secondary and not read.is_supplementary \
                         and not read.is_duplicate]
            secondary_pairs = self.pair_reads(secondary)
            rows.extend(self.create_row(pair, label) for pair in secondary_pairs if pair)

        return rows

    def pair_reads(self, reads):
        paired_reads = []
        first_in_pair = {}
        second_in_pair = {}

        for read in reads:
            if read.is_read1:
                first_in_pair[read.query_name] = read
            elif read.is_read2:
                second_in_pair[read.query_name] = read

        all_query_names = set(first_in_pair.keys()).union(second_in_pair.keys())
        for name in all_query_names:
            first_read = first_in_pair.get(name, None)
            second_read = second_in_pair.get(name, None)
            paired_reads.append((first_read, second_read))

        return paired_reads

    def create_row(self, reads, label):
        read1, read2 = reads
        row = {
            label + '_name': self.name,
            label + '_flags': [read.flag if read else None for read in reads],
            label + '_pos': [read.reference_start + 1 if read and not read.is_unmapped else None for read in reads],
            label + '_chr': [read.reference_name if read and not read.is_unmapped else None for read in reads],
            label + '_CIGAR': [read.cigarstring if read and not read.is_unmapped else None for read in reads],
            label + '_quality': [read.mapping_quality if read and not read.is_unmapped else None for read in reads],
            label + '_edit_dist': [self.calculate_edit_distance(read) if read and not read.is_unmapped else None for
                                   read in reads],
            label + '_MD': [self.get_tag_safe(read, 'MD') if read and not read.is_unmapped else None for read in reads],
            label + '_type': [self.determine_read_type(read) for read in reads],
            label + '_proper_pair': self.check_proper_pair(read1)  # Check if each read is in a proper pair
        }
        return row

    def check_proper_pair(self, read):
        """Check if the read is part of a proper pair."""
        if read and not read.is_unmapped:
            return 1 if read.is_proper_pair else 0
        return 0

    def get_tag_safe(self, read, tag):
        if read and not read.is_unmapped:
            try:
                return read.get_tag(tag)
            except KeyError:
                return None
        return None

    # def calculate_edit_distance(self, read):
    #     """Calculate the edit distance for a read. If the read is unmapped, return its length as the edit distance."""
    #     if read and read.is_unmapped:
    #         return read.query_length  # Return the total length of the read if it is unmapped.
    #     elif read and read.cigarstring:
    #         # Calculate the edit distance based on the CIGAR string.
    #         matches = sum(
    #             length for op, length in read.cigar if op == 0)  # Sum lengths where operation is 'M' (match/mismatch)
    #         return read.query_length - matches
    #     return None  # Return None if read is None or other unexpected cases.

    def calculate_edit_distance(self, read):
        """Calculate the edit distance for a read using the NM tag if available, else calculate manually."""
        if read and not read.is_unmapped:
            try:
                return read.get_tag('NM')
            except KeyError:
                # If NM tag is not present, fall back to manual calculation
                edit_distance = 0
                for op, length in read.cigar:
                    if op in {1, 2, 8}:  # Insertion (I), Deletion (D), and Mismatch (X) operations
                        edit_distance += length
                return edit_distance
        return None

    def determine_read_type(self, read):
        if not read:
            return 'U'
        elif read.is_secondary:
            return 'S'
        elif read.is_unmapped:
            return 'U'
        else:
            return 'P'

class CSV:
    def __init__(self, path: str, label: str):
        self.path = path
        self.label = label
        self.header = [f'{self.label}_{column}' for column in [
            'name', 'flags', 'pos', 'chr', 'CIGAR', 'quality', 'edit_dist', 'MD', 'type', 'proper_pair'
        ]]
        self.write_header()
        self.rows_buffer = []
        self.BUFFER_SIZE = 1_000_000

    def write_header(self):
        with open(self.path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(self.header)

    def append_row(self, row):
        self.rows_buffer.append(row)
        if len(self.rows_buffer) >= self.BUFFER_SIZE:
            self.flush()

    def flush(self):
        with open(self.path, 'a', newline='') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=self.header)
            writer.writerows(self.rows_buffer)
        self.rows_buffer = []

    def close(self):
        if self.rows_buffer:
            self.flush()

class Parser:
    def __init__(self, input_file: str, label: str, output_path: str, fastq_file: str = None, include_secondary: bool = False):
        self.input_file = input_file
        self.label = label
        self.output_path = output_path
        self.include_secondary = include_secondary
        self.output_csv = CSV(self.output_path, self.label)
        self.expected_reads = self.load_fastq_reads(fastq_file) if fastq_file else set()
        self.seen_reads = set()
        self.all_rows = []

    def load_fastq_reads(self, fastq_file):
        """Load read names from the FASTQ file."""
        read_names = set()
        with open(fastq_file, 'r') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:  # FASTQ headers are every 4 lines starting from 0
                    read_name = line.strip().split()[0][1:]  # Remove '@' and take the read name
                    read_names.add(read_name)
        return read_names

    def run(self):
        sorted_input_file = self.sort_bam_file(self.input_file)
        with pysam.AlignmentFile(sorted_input_file, 'rb') as reads:
            first_line = next(reads)
            reads_same_name = SameNameReadsArray(first_line, include_secondary=self.include_secondary)
            self.seen_reads.add(first_line.query_name)
            for read in reads:
                if read.query_name == reads_same_name.name:
                    reads_same_name.append(read)
                else:
                    rows = reads_same_name.process_reads(self.label)
                    self.all_rows.extend(rows)
                    reads_same_name = SameNameReadsArray(read, include_secondary=self.include_secondary)
                    self.seen_reads.add(read.query_name)
            rows = reads_same_name.process_reads(self.label)
            self.all_rows.extend(rows)

        self.add_missing_reads()
        self.write_sorted_csv()

    def add_missing_reads(self):
        """Add rows for reads that were expected but not seen in the BAM file."""
        missing_reads = self.expected_reads - self.seen_reads
        for read_name in missing_reads:
            row = {
                f'{self.label}_name': read_name,
                f'{self.label}_flags': [None, None],
                f'{self.label}_pos': [None, None],
                f'{self.label}_chr': [None, None],
                f'{self.label}_CIGAR': [None, None],
                f'{self.label}_quality': [None, None],
                f'{self.label}_edit_dist': [None, None],
                f'{self.label}_MD': [None, None],
                f'{self.label}_type': ['U', 'U'],
                f'{self.label}_proper_pair': 0
            }
            self.all_rows.append(row)

    def write_sorted_csv(self):
        """Sort the rows by the numeric part of the read name and write them to the CSV."""

        def extract_numeric_suffix(name):
            match = re.search(r'\.(\d+)$', name)
            return int(match.group(1)) if match else 0

        self.all_rows.sort(key=lambda x: extract_numeric_suffix(x[f'{self.label}_name']))
        for row in self.all_rows:
            self.output_csv.append_row(row)
        self.output_csv.close()

    def sort_bam_file(self, input_file):
        """Sort the BAM file by read name."""
        output_sorted = input_file.replace('.bam', '_sorted_by_name.bam')
        pysam.sort('-o', output_sorted, '-n', input_file)
        return output_sorted

app = typer.Typer()

@app.command()
def main(input_file: Path = typer.Argument(
    ...,
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
    show_default=False),
        output_path: Path = typer.Option(
            "./parser_output.csv",
            "--output",
            "-o",
            help="Specify the path to the output CSV file.",
            show_default=True),
        fastq_file: Path = typer.Option(
            None,
            "--fastq",
            "-f",
            help="Specify the path to the input FASTQ file to ensure all reads are included."),
        include_secondary: bool = typer.Option(
            False,
            "--include-secondary",
            "-s",
            help="Include secondary alignments in the output")):
    label = input_file.stem
    parser = Parser(input_file=str(input_file), label=label, output_path=str(output_path),
                    fastq_file=str(fastq_file) if fastq_file else None, include_secondary=include_secondary)
    parser.run()

if __name__ == "__main__":
    app()
