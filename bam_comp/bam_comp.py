import pysam
import pandas as pd
import psutil
import glob
from os import path, makedirs, getcwd, getpid
import pickle
from shutil import rmtree


def cache(py_object, filename, *subdirectories):
    # creating subdirectories if needed
    abs_path = path.join(getcwd(), *subdirectories)
    makedirs(abs_path, exist_ok=True)
    # creating cache file and dumping objects into it
    filepath = path.join(abs_path, str(filename) + '.tmp')
    with open(filepath, 'wb+') as File:
        pickle.dump(py_object, File)


def decache(filename, *subdirectories):
    # importing stored objects from cache file
    abs_path = path.join(
        getcwd(), *subdirectories, str(filename) + '.tmp'
        )
    with open(abs_path, 'rb') as File:
        decached = pickle.load(File)
    return decached


def delete_temp_directory(folder_name: str='temp_cache'):
    """Deletes temp directory and all its contents"""
    if path.exists(folder_name):
        rmtree(folder_name)


class BamComp:
    def __init__(self, input_data=[], output_path='', is_single_ended=0, enable_caching=False):
        if len(output_path) > 0 and output_path[-4:] != '.csv' and '\\' != output_path[-1] != '/': output_path += '/'
        self.input_data = input_data #[{'label': 'Label1', 'path': '/path/input_1.bam', 'tool': 'tool_name'}, ...]
        self.output_path = output_path if output_path[-4:] == '.csv' else output_path + 'BamComp_res.csv' #'/path/output.csv'
        self.is_single_ended = is_single_ended
        self.comp_data = {}
        self.enable_caching = enable_caching
        self.column_names = ['_flags', '_pos', '_chr', '_CIGAR', '_edit_dist', '_quality', '_MD', '_multi']


    @staticmethod
    def get_current_memory_usage(in_gigabytes=True) -> float:
        """Checks the amount of RAM used by the script at the time"""
        process = psutil.Process(getpid())
        usage = process.memory_info().rss # in bytes
        usage_GB = usage * 9.31 * 10 ** (-10) # in Gigabytes
        return round(usage_GB, 3) if in_gigabytes else usage


    # Append more files to existing class object
    def append_data(self, label, path, tool):
        self.input_data.append({'label':label, 'path':path, 'tool':tool})


    def append_to_csv(self):
        pass


    # Save result. If multiple_tables set to 1, multiple files will be created (one for each label)
    def save(self, multiple_tables = 0, output_path = ''):
        if len(output_path) > 0 and output_path[-4:] != '.csv' and '\\' != output_path[-1] != '/': output_path += '/'
        if len(output_path) > 0: self.output_path = output_path if output_path[-4:] == '.csv' else output_path + 'BamComp_res.csv'
        if not multiple_tables:
            with open(self.output_path, 'w'):
                pass
            if self.enable_caching:
                for filename in glob.iglob(f"temp_cache_{self.input_data[0]['label']}/cached_reads_*"):
                    rows_dict = decache(filename.split('/')[-1].split('.')[-2], f"temp_cache_{self.input_data[0]['label']}")
                    pandas_rows = pd.DataFrame.from_dict(rows_dict).T
                    pandas_rows.to_csv(path_or_buf=self.output_path, mode='a')
                delete_temp_directory(f"temp_cache_{self.input_data[0]['label']}")
            else:
                self.comp_data.to_csv(path_or_buf=self.output_path)
            print('Saved as csv.', self.output_path)
        else:
            if not self.enable_caching:
                all_columns = self.comp_data.columns
                output_paths = []
                for file in self.input_data:
                    output_paths.append(self.output_path[:-4] +  '_' + file['label'] + self.output_path[-4:])
                    file_columns = [column for column in all_columns if column[:len(file['label'])] == file['label']]
                    self.comp_data.to_csv(path_or_buf=output_paths[-1], columns=file_columns)
                print('Saved as csv. Path:\n', '\n'.join(output_paths), sep='')
            elif self.enable_caching:
                raise Exception('Caching does not supported in multiple tables mode')


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
            reads_name_counter = {}
            for read in reads:
                reads_name_counter[read.query_name] = reads_name_counter.get(read.query_name, 0) + 1
                if bin(read.flag)[-12:-11] == '1': continue # skipping supplementary aligned reads (0x800)
                place = 1 if bin(read.flag)[-8:-7] == '1' else 0 # position of read in pair (0x80)
                if tool[:3] == 'bwa' or tool == 'ngm':
                    primary_reads[place] = read
                    # In order to check if reads are multimapped looking for specific flags (BWA and NGM tools)
                    if (tool[:3] == 'bwa' and None != find_tag(read.tags, 'XA')) or (tool == 'ngm' and None != find_tag(read.tags, 'X0') > 0):
                        mm[place] = 1
                    else: mm[place] = 0
                else:
                    if bin(read.flag)[-9:-8] == '1': mm[place] = 1 # Checking if read has secondary alignment (0x100)
                    else: primary_reads[place] = read
                
            for i, read in enumerate(primary_reads):
                if read != None and mm[i] == None: mm[i] = 0
                if read != None and reads_name_counter[read.query_name] > 2:
                    mm[i] = 1
            return primary_reads, mm
        read_counter = 0
        cached_files = 0

        for file in self.input_data:
            label, tool_name = file['label'], file['tool']
            if file['path'][-4:] == '.bam': qualifier = 'rb' # for .BAM format files
            elif file['path'][-4:] == '.sam': qualifier = 'r' # for .SAM format files

            with pysam.AlignmentFile(file['path'], qualifier) as f:
                next_line = next(f)
                unique_reads = {}
                for line in f: #reading file line by line
                    read_counter += 1
                    reads = [line] #list of reads
                    if next_line.query_name != reads[0].query_name: reads[0], next_line = next_line, reads[0]
                    else: 
                        reads.append(next_line)
                        while True:
                            try: next_line = next(f)
                            except StopIteration: break #in case end of file is reached
                            if next_line.query_name == reads[0].query_name : #looking for rows with same read ID
                                reads.append(next_line)
                            else: 
                                break
                    
                    for read in reads:
                        if read.query_name in unique_reads:
                            unique_reads[read.query_name].append(read)
                        else:
                            unique_reads[read.query_name] = [read]


                for read_ID, reads in unique_reads.items():
                    print(len(reads), read_ID, [read.reference_name for read in reads])
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
                    if self.enable_caching:
                        if read_counter % 1000000 == 0:
                            cached_files += 1
                            print(f'{read_counter} reads cached, {self.get_current_memory_usage()}G RAM used')
                            cache(self.comp_data, f'cached_reads_{str(cached_files).zfill(3)}', f"temp_cache_{self.input_data[0]['label']}")
                            self.comp_data = {}
                            

        if self.enable_caching:
            cached_files += 1
            cache(self.comp_data, f'cached_reads_{str(cached_files).zfill(3)}', f"temp_cache_{self.input_data[0]['label']}")
            self.comp_data = {}

        self.comp_data = pd.DataFrame.from_dict(self.comp_data).T
