# Binary Alignment Map Files Comparer

Set of 3 Python scripts to compare alignment data among replicates.

## Installation

First, create a virtual environemnt.
Then install dependencies via pip:

```bash
pip install -r requirements.txt
```

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

## II. CSV Subsampling (only for real data samples)

1. Edit `run_subsampler.py` by specifying path to each CSV. Script can subsample CSVs only pairwise at the moment.

```python3
csv_1 = CSV('/home/user/csv_1.csv') # path to first CSV
csv_2 = CSV('/home/user/csv_2.csv') # path to second CSV
```

2. And run the script:

```bash
python run_subsampler.py
```

## III. Multiple CSVs -> Features counter

Finally, in order to compare CSVs between each other:

1. Edit `run_comparer.py`, specifying label and path to each CSV. Script can compare CSV both pairwise and triplewise. In the latter case just leave the third field blank.

2. And run the script:

```bash
python run_comparer.py
```
