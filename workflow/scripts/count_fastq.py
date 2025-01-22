from Bio import SeqIO
from pathlib import Path
import typer
import csv
import os
from collections import Counter


CACHE_FILE = "./reads.csv"


# Function to read CSV file and return a set of existing sample IDs
def read_csv():
    existing_sample_ids = dict()
    if not os.path.isfile(CACHE_FILE):
        return existing_sample_ids  # File doesn't exist yet
    
    with open(CACHE_FILE, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            existing_sample_ids.update(
            {
                row['sampleID']: row['number_of_reads']
                }
    )
    
    return existing_sample_ids


def write_csv(data):
    with open(CACHE_FILE, 'a', newline='') as csvfile:
        fieldnames = ['sampleID', 'number_of_reads']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
        # Check if file is empty and write header if needed
        if os.path.getsize(CACHE_FILE) == 0:
            writer.writeheader()
        
        writer.writerow(data)


def count_reads(fastq_file: Path):
    # Get existing sample IDs from CSV file
    existing_sample_ids = read_csv()

    generator = SeqIO.parse(fastq_file, "fastq")
    read_counts = Counter()

    for rec in generator:
        sample_name = rec.id.split(".")[0]

        # Check if sample ID exists in the set of existing sample IDs
        if sample_name in existing_sample_ids:
            return int(existing_sample_ids[sample_name])

        read_counts[sample_name] += 1

    if read_counts[sample_name] != 0:
        write_csv(
            {
                'sampleID': sample_name, 
                'number_of_reads': read_counts[sample_name]
                }
        )
    return read_counts[sample_name]


if __name__ == "__main__":
    app = typer.Typer(add_completion=False)

    @app.command()
    def main(fastq_file: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        show_default=False,
        help="Path to FASTQ File to count reads within"
    ),
    ):
        count_reads(fastq_file=fastq_file)

    app()