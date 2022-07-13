from bam_comp import BamComp

# May need to install the pysam package for python: 
#     python3 -m pip install pysam

# Creating class object. 
#     input_data requires LIST of DICTIONARIES with folliwing STRING fields: label, path, tool. (can be empty and field later through append_data function)
#     output_path requires a STRING path to output directory or file
# -------
# Tool field in input_data can be empty ('tool':'') or contain tool tag (Should be field for ngm and bwa tools, because of specific tags). 
# It should work with SAM and BAM files. 
# Also should work with single end reads if "is_single_ended" set to 1 when creating class object (set to 0 by default).
# -------
test = BamComp(
    input_data=[
        {'label':'bowtie_1', 'path':'~/sample.sorted.bam', 'tool':''}
        ],
    output_path='',
    is_single_ended=0,
    enable_caching=True
    )

# As mention before new files can be added with append_data function. Function requires STRING label, path and tool parameters
# test.append_data(label='BOWTIE2_1', path='~/GW/bowtie2_res.bam', tool='')

# compare function will start parsing for each added file one by one
test.compare()

# After using compare function you can accsess created table by using *class_object_name*.comp_data
# it will return PANDAS DataFrame with following columns (ID and 8 columns for each label):
#   ID                  STRING          read index
#   flags               LIST (INT)      SAM format flag for each read
#   pos                 LIST (INT)      position of primary aligned reads in the pair
#   chr                 LIST (STR)      chromosome
#   CIGAR               LIST (STR)      cigar strings
#   edit_dist           LIST (INT)      edit distance
#   quality             LIST (INT)      quality of reads given by tool
#   MD                  LIST (STR)      mismatch - MD tag (if available)
#   multi               LIST (INT)      multimapping (0 - no multimapping / 1 - read multimapped)
df = test.comp_data

# Save result in output_path folder. If multiple_tables set to 1, multiple files will be created (one for each label). Parameters:
#   multiple_tables     INTEGER     0 (default) or 1. If set to 0 - one large table with all labels will be created.
#   output_path         STRING      optionall. Will replace output_path set when class object was created.
test.save()
