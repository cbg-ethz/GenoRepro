# File Structure

Suggested project file structure draft

```bash
├── sample_data
│   ├── plots
│   │   ├── ...
│   ├── raw_data
│   │   ├── simulated_g.bam
│   │   ├── simulated_rv.bam
│   │   ├── simulated_rv_s.bam
│   │   ├── ...
│   ├── processed_data
│   │   ├── simulated_samples
│   │   │   ├── E0
│   │   │   │   ├── bowtie2.csv
│   │   │   │   ├── bwa.csv
│   │   │   │   ├── ...
│   │   │   ├── EH
│   │   │   │   ├── bowtie2.csv
│   │   │   │   ├── bwa.csv
│   │   │   │   ├── ...
│   │   │   ├── RL
│   │   │   │   ├── bowtie2.csv
│   │   │   │   ├── bwa.csv
│   │   │   │   ├── ...
│   │   ├── synthetic_replicates
│   │   │   ├── bowtie2
│   │   │   │   ├── ERR009308
│   │   │   │   │   ├── g_rv.csv
│   │   │   │   │   ├── g_s1.csv
│   │   │   │   │   ├── g_s2.csv
│   │   │   │   │   ├── g_s3.csv
│   │   │   │   ├── ERR009309
│   │   │   │   │   ├── g_rv.csv
│   │   │   │   │   ├── g_s1.csv
│   │   │   │   │   ├── g_s2.csv
│   │   │   │   │   ├── g_s3.csv
│   │   │   │   ├── ...
│   │   │   ├── bwa
│   │   │   │   ├── ERR009308
│   │   │   │   │   ├── g_rv.csv
│   │   │   │   │   ├── g_s1.csv
│   │   │   │   │   ├── g_s2.csv
│   │   │   │   │   ├── g_s3.csv
│   │   │   │   ├── ERR009309
│   │   │   │   │   ├── g_rv.csv
│   │   │   │   │   ├── g_s1.csv
│   │   │   │   │   ├── g_s2.csv
│   │   │   │   │   ├── g_s3.csv
│   │   │   │   ├── ...
│   │   │   ├── ...
│   │   ├── true_position
│   │   │   ├── bowtie2
│   │   │   │   ├── ...
│   │   │   ├── bwa
│   │   │   │   ├── ...
│   │   │   ├── ...
│   └── partials/template
```
