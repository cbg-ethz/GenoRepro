# BAM Parser

## How to run BAM Parser

### 0. Requirements

Before running the BAM parser make sure that:

* BAM file is sorted and indexed, there should be 2 files in the working directory:
  * <...>.sorted.bam
  * <...>.sorted.bam.bai

In order to sort and index BAM file run

```bash
samtools sort -o <sorted.bam> <initial.bam>
```

Then

```bash
samtools index <sorted.bam>
```

* conda environment contains all the necessary libraries

To install required packages run:

```bash
conda install --file requirements.txt
```

### I. BAM -> CSV

Open `test.py` file and edit `input_data` array, specifying label, path to the BAM file and tool field in case BAM file was generated using bwa or ngm tool, otherwise it should be left empty.
Output path can be specified too, otherwise the output CSV will be written in the working directory.

#### Caching

In order to minimize memory usage while parsing large (>15 GB) files caching may be enabled.

```python3
enable_caching=True
```

Though there are some exceptions to bear in mind:

1. Caching is not supported in multiple tables mode, so make sure this feature is disabled at the bottom of the script.

```python3
test.save()
```

or

```python3
test.save(multiple_tables=0)
```

And run the script:

```bash
python test.py
```

2. If caching is enabled, parse only one sample per run.

3. For real data parsing an additional column with read sequence must be generated, so the corresponding argument should be passed as True:

```python3
is_real_data=True
```

### III. CSV Subsampling (only for real data samples)

1. Edit `run_subsampler.py` by specifying path to each CSV. Script can subsample CSVs only pairwise at the moment.

```python3
csv_1 = CSV('/home/user/csv_1.csv') # path to first CSV
csv_2 = CSV('/home/user/csv_2.csv') # path to second CSV
```

2. And run the script:

```bash
python run_subsampler.py
```


### IV. Multiple CSVs -> Features counter

Finally, in order to compare CSVs between each other:

1. Edit `run_comparer.py`, specifying label and path to each CSV. Script can compare CSV both pairwise and triplewise. In the latter case just leave the third field blank.

2. And run the script:

```bash
python run_comparer.py
```

## Links

* [Google Drive Folder](https://drive.google.com/drive/folders/1e54IloZcnRdownjMEaMoNSOCQtakF47z)
