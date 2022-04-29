import pysam
import pandas as pd

class BamComp:
    def __init__(self, input_data=[], output_path='', is_single_ended=0):
        if len(output_path) > 0 and output_path[-4:] != '.csv' and '\\' != output_path[-1] != '/': output_path += '/'
        self.input_data = input_data #[{'label': 'Label1', 'path': '/path/input_1.bam', 'tool': 'tool_name'}, ...]
        self.output_path = output_path if output_path[-4:] == '.csv' else output_path + 'BamComp_res.csv' #'/path/output.csv'
        self.is_single_ended = is_single_ended
        self.comp_data = {}
        self.column_names = ['_flags', '_pos', '_chr', '_CIGAR', '_edit_dist', '_quality', '_MD', '_multi']

    # Append more files to existing class object
    def append_data(self, label, path, tool):
        self.input_data.append({'label':label, 'path':path, 'tool':tool})

    # Save result. If multiple_tables set to 1, multiple files will be created (one for each label)
    def save(self, multiple_tables = 0, output_path = ''):
        if len(output_path) > 0 and output_path[-4:] != '.csv' and '\\' != output_path[-1] != '/': output_path += '/'
        if len(output_path) > 0: self.output_path = output_path if output_path[-4:] == '.csv' else output_path + 'BamComp_res.csv'
        if not multiple_tables:
            self.comp_data.to_csv(path_or_buf=self.output_path)
            print('Saved as csv.', self.output_path)
        else:
            all_columns = self.comp_data.columns
            output_paths = []
            for file in self.input_data:
                output_paths.append(self.output_path[:-4] +  '_' + file['label'] + self.output_path[-4:])
                file_columns = [column for column in all_columns if column[:len(file['label'])] == file['label']]
                self.comp_data.to_csv(path_or_buf=output_paths[-1], columns=file_columns)
            print('Saved as csv. Path:\n', '\n'.join(output_paths), sep='')

    def compare(self):
        def get_edit_dist(cigar, reads_len):
            match = 0
            for el in cigar:
                if el[0] == 0: match+=el[1]
            return reads_len - match
        
        def find_tag(tags, target):
            for tag in tags:
                if tag[0] == target: return tag[-1]
            return None
        
        # Returns only primary aligned reads and checking for multimapping.
        def get_primary(reads, tool):
            primary_reads, mm = [None, None], [None, None] # for reads without a pair will be set to None
            for read in reads:
                if bin(read.flag)[-12:-11] == '1': continue # skipping supplementary aligned reads (0x800)
                place = 1 if bin(read.flag)[-8:-7] == '1' else 0 # position of read in pair (0x80)
                if tool[:3] == 'bwa' or tool == 'ngm':
                    primary_reads[place] = read
                    # In order to check if reads are multimapped looking for specific flags (BWA and NGM tools)
                    if (tool[:3] == 'bwa' and None != find_tag(read.tags, 'XA')) or (tool == 'ngm' and None != find_tag(read.tags, 'X0') > 0): 
                        mm[place] = 1
                    else: mm[place] = 0
                else:
                    if bin(read.flag)[-9:-8] == '1': mm[place] = 1 # Checking if read have secondary alignment (0x100)
                    else: primary_reads[place] = read
            for i, read in enumerate(primary_reads):
                if read != None and mm[i] == None: mm[i] = 0
            return primary_reads, mm

        for file in self.input_data:
            label, tool_name = file['label'], file['tool']
            if file['path'][-4:] == '.bam': qualifier = 'rb' # for .BAM format files
            elif file['path'][-4:] == '.sam': qualifier = 'r' # for .SAM format files

            with pysam.AlignmentFile(file['path'], qualifier) as f:
                next_line = next(f)
                for line in f: #reading file line by line 
                    reads = [line] #list of reads
                    if next_line.query_name != reads[0].query_name: reads[0], next_line = next_line, reads[0]
                    else: 
                        reads.append(next_line)
                        while 1:
                            try: next_line = next(f)
                            except StopIteration: break #in case end of file is reached
                            if next_line.query_name == reads[0].query_name : #looking for rows with same read ID
                                reads.append(next_line)
                            else: break
                    
                    index = reads[0].query_name 
                    reads, mm = get_primary(reads, tool_name)
                    if self.comp_data.get(index, None) == None: self.comp_data[index] = {}
                    if self.is_single_ended:
                        self.comp_data[index].update({
                            label + '_flags': reads[0].flag,
                            label + '_pos': reads[0].reference_start + 1, #indexing starts from 0
                            label + '_chr': reads[0].reference_name,
                            label + '_CIGAR': reads[0].cigarstring,
                            label + '_edit_dist': get_edit_dist(reads[0].cigar, len(reads[0].query_sequence)),
                            label + '_quality': reads[0].mapping_quality,
                            label + '_MD': find_tag(reads[0].tags, 'MD'),
                            label + '_multi': mm[0]
                        })
                    else:
                        parsed_read = {label + column_name : [None, None] for column_name in self.column_names}
                        for i, read in enumerate(reads):
                            if read == None: continue
                            parsed_read[label + '_flags'][i] = read.flag
                            parsed_read[label + '_pos'][i] = read.reference_start + 1
                            parsed_read[label + '_chr'][i] = read.reference_name
                            parsed_read[label + '_CIGAR'][i] = read.cigarstring
                            parsed_read[label + '_edit_dist'][i] = get_edit_dist(read.cigar, len(read.query_sequence))
                            parsed_read[label + '_quality'][i] = read.mapping_quality
                            parsed_read[label + '_MD'][i] = find_tag(read.tags, 'MD')
                            parsed_read[label + '_multi'][i] = mm[i]
                        self.comp_data[index].update(parsed_read)

        self.comp_data = pd.DataFrame.from_dict(self.comp_data).T
