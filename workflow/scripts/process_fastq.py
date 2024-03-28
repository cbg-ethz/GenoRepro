import csv
import os
from itertools import zip_longest

import typer
from pathlib import Path

from Bio import SeqIO
from collections import Counter

from click.exceptions import Exit

app = typer.Typer(add_completion=False)

CACHE_FILE = "./reads.csv"


def write_csv(data, file_name):
    with open(file_name, 'w', newline='') as csvfile:  # 'w' mode to write anew
        fieldnames = ['sampleID', 'number_of_reads']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(data)


def process_reads(out_file, in_file1: Path, in_file2: Path = None, data_type: str = 'single'):
    print(out_file)
    generator1 = SeqIO.parse(in_file1, "fastq")
    if data_type == 'paired':
        generator2 = SeqIO.parse(in_file2, "fastq")
        read_iter = zip_longest(generator1, generator2, fillvalue=None)
    else:
        read_iter = ((rec1, None) for rec1 in generator1)

    read_counts1 = Counter()
    for rec1, rec2 in read_iter:
        if rec1 is None or (rec2 is None and data_type == 'paired'):
            raise ValueError("The input files have different number of records.")
        sample_name1 = rec1.id.split(".")[0]

        read_counts1[sample_name1] += 1

    if read_counts1:
        for sample_name, count in read_counts1.items():
            print(sample_name)

            write_csv({'sampleID': sample_name, 'number_of_reads': count}, out_file)


@app.command()
def main(

        fastq_file1: Path = typer.Option(
            ...,
            "--fastq_file1", "-f1",
            help="Name of the first FASTQ file for paired data"
        ),

        fastq_file2: Path = typer.Option(
            None,
            "--fastq_file2", "-f2",
            help="Name of to the second FASTQ file for paired data"
        ),

        pair_type: str = typer.Argument(
            ...,
            exists=True,
            file_okay=False,
            dir_okay=False,
            show_default=False,
            help="Paired or single"

        ),

        output_file: str = typer.Argument(
            ...,
            exists=True,
            file_okay=False,
            dir_okay=True,
            show_default=False,
            help="Output folder"

        ),
):
    if pair_type not in ["paired", "single"]:
        typer.echo("Error: Invalid pair_type. Must be 'paired' or 'single'.")
        raise typer.Abort()
    if pair_type == "paired" and (not fastq_file1.exists() or not fastq_file2.exists()):
        typer.echo(f"Error: One or both FASTQ files cannot be found.")
        raise typer.Abort()
    if pair_type == "single" and not fastq_file1.exists():
        typer.echo("Error: FASTQ file cannot be found.")
        raise typer.Abort()

    # read_csv(input_folder)
    process_reads(output_file, fastq_file1, fastq_file2 if pair_type == "paired" else None, pair_type)


if __name__ == "__main__":
    app()
