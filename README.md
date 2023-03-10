# Binary Alignment Map Files Comparer

Set of 2 Python scripts to compare alignment data among replicates.

## Installation

First, create a virtual environemnt.
Then install dependencies via pip:

```bash
pip install -r requirements.txt
```

## Pipeline

![alt text](data/img/pipeline.svg "Title")

## I. BAM -> CSV

Parser `parser.py` is a command line Python application. \
To show available commands run:

```bash
python parser.py --help
```

It requires only path to the input BAM File as a sole argument. There is no need to sort and index BAM Files prior to parsing as well as to specify whether it contains paired end or single end reads. Simply run:

```bash
python parser.py <input_file.bam>
```

By default, it creates a CSV File named `parser_output.csv` in the current directory. The output can be specified by providing an absolute or relative path with `-o` option:

```bash
python parser.py -o <output_file.csv> <input_file.bam>
```

It is also recommended to provide a unique label corresponding to each replicate with `-l` option. For further comparison of replicates of different sizes flag `--real-data` or `-r` should be given.

### Output

Generated CSV File will contain only unambiguous (non-multimapped primary alignments) reads (pairs of reads) if flag `--keep-multimapped_reads` was not given.

All the supportive data (e.g. number of reads by type, output path, memory usage) is printed to the stdout.

### Time complexity

Time complexity is O(n) for the algorithm, where n is the number of reads in the BAM File.

### Space complexity

Space complexity is O(1) for the algorithm. Regardless of the replicate size, the script will use no more than 3 GB of memory at peak.

## II. Comparison

To compare previously obtained CSV files use `compare.py`, which is also a command line Python application. \

To show available commands run:

```bash
python comparer.py --help
```

The data will be automatically subsampled if tables of different size are given.

During subsampling `input_file.rv.csv` and `input_file.sub.csv` will be generated (corresponding to the sorted table of reversed complemented reads and table with only common reads) for each `input_file.csv`.

To run:
```bash
python comparer.py -o <features_counts.csv> <input_file1.csv> <input_file2.csv> ...
```

### Output

Generated CSV File contains list of counters of selected features.

### Time complexity

Time complexity is O(n) for the subsampling process and O(1) for comparison process.

### Space complexity

Space complexity is O(1) for the subsampling process and O(n) for comparison process.
