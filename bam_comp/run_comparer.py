from comparer import Comparer


# REQUIREMENTS
# -------------------------------------
# Some packages may not be already installed. To install run:
#   conda install -c conda-forge dataclasses psutil
#   conda install -c anaconda pandas


# INPUT
# -------------------------------------
# Fields to be filled in the input_data:
#   label       | STRING | Should be the same label that was specified when parsing
#   path_to_csv | STRING | Path to the CSV file which contains columns with the corresponding label given above;
#                          Supports multiple tables as well as one large table with all labels
#
# Script generates another output CSV file. Its location and name can be specified below with output_path variable
#   output_path | STRING | (OPTIONAL) if left blank creates 'compare_output.csv' in the same folder by default


# BEHAVIOUR
# -------------------------------------
# Script can process any combination of samples:
#   Original vs Reversed, Original vs Shuffled, Original vs both Reversed and Shuffled.
#
# For pairwise comparison (Original vs Reversed, Original vs Shuffled):
#   Leave one of the fields in input_data blank.
#   For example, to compare Original sample with Shuffled leave input for Reversed sample blank.
# For triple-wise comparison (Original vs both Reversed and Shuffled):
#   All the fields in the input_data should be filled.


input_data = {
    'original': {
                    'label': 'bowtie_SE',
                    'path_to_csv': '/Users/mike/Documents/Projects/Reproducibility/real_data/10_5_no_duplicates_sub.csv'
                },

    'reversed': {
                    'label': 'bowtie_11',
                    'path_to_csv': '/Users/mike/Documents/Projects/Reproducibility/real_data/11_5_no_duplicates_sub.csv'
                },

    'shuffled': {
                    'label': '',
                    'path_to_csv': ''
                },
}

output_path = '' # e.g. '~/scratch/project/run_1/output.csv'


a = Comparer(input_data, output_path, is_real_data=True)
a.compare()

# OUTPUT
# -------------------------------------
# Script generates output CSV file.
# In case of pairwise comparison:
#   --------------------------------------------------------------------------------
#   FEATURE             |       NUMBER OF READS             |       PERCENTAGE     |
#   --------------------------------------------------------------------------------
#   Total_reads         | Initial number of reads           |     Always 100%      |
#                       | (unfiltered, raw)                 |                      |
#   --------------------------------------------------------------------------------
#   Mapped_reads        | Number of reads without those     |                      |
#                       | which were mapped neither with    |  from Total_reads    |
#                       | original nor with replicated data |                      |
#   --------------------------------------------------------------------------------
#   Unambiguous_{label} | Number of unambiguous reads       |  from Mapped_reads   |
#                       | by their label                    |                      |
#   --------------------------------------------------------------------------------
#   Common_unambiguous  | Number of common unambiguous reads|  from Mapped_reads   |
#   --------------------------------------------------------------------------------
#   Inconsistent_type1  | Unambiguous reads mapped only     |         from         |
#                       | with original data                |  Common_unambiguous  |
#   --------------------------------------------------------------------------------
#   Inconsistent_type2  | Unambiguous reads mapped only     |         from         |
#                       | with replicated data              |  Common_unambiguous  |
#   --------------------------------------------------------------------------------
#   Identical           | Identical reads in terms of       |         from         |
#                       | position and edit distance        |  Common_unambiguous  |
#   --------------------------------------------------------------------------------
#   Consistent_global_  | Common unambiguous reads mapped   |         from         |
#   _inconsistent_local | to the same position with         |  Common_unambiguous  |
#                       | different edit distance           |                      |
#   --------------------------------------------------------------------------------
#   Inconsistent_global | Common unambiguous reads mapped   |         from         |
#                       | to different positions            |  Common_unambiguous  |
#                       | with different edit distance      |                      |
#   --------------------------------------------------------------------------------
#   Multi_mapped        | Common unambiguous reads mapped   |         from         |
#                       | to different positions            |  Common_unambiguous  |
#                       | with the same edit distance       |                      |
#   --------------------------------------------------------------------------------
