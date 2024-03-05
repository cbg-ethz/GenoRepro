#indexing
#!/bin/bash

# Assuming the script receives two arguments:
# 1. Genome index
# 2. output indices

genome_index=$1


# Your alignment command, redirecting stderr to the log file
bwa-mem2 mem index "$genome_index"
