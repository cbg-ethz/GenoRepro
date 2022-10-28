from os import path, remove, rename, getpid
from prettytable import PrettyTable
from csvsort import csvsort
from pathlib import Path
import resource
import psutil
import pysam
import typer
import time
import csv


class SameNameReadsArray(list):
    def __init__(self, *reads: pysam.AlignedSegment):
        self.name = reads[0].query_name if len(reads) > 0 else None
        super(SameNameReadsArray, self).__init__([*reads])
    

    def keep_only_primary_reads(self):
        filtered_reads = SameNameReadsArray()
        for read in self:
            if (
                not read.is_secondary and
                not read.is_supplementary and
                not read.is_duplicate
                ):
                filtered_reads.append(read)
        return filtered_reads
    

    def export(self, label: str,
                     is_paired_ended: bool,
                     is_multimapped: bool,
                     is_real_data: bool) -> dict:

        def get_edit_dist():
            """Calculates edit distance"""
            def calculate_edit_dist(read):
                cigar = read.cigar
                reads_len = read.query_length
                match = 0
                for el in cigar:
                    if el[0] == 0: match += el[1]
                return reads_len - match
            if is_paired_ended:
                return [calculate_edit_dist(read) for read in self]
            else:
                return calculate_edit_dist(self[0])
        
        def get_tag():
            def find_tag(read):
                tags, target = read.tags, 'MD'
                for tag in tags:
                    if tag[0] == target: return tag[-1]
                return None
            if is_paired_ended:
                return [find_tag(read) for read in self]
            else:
                return find_tag(self[0])

        
        def get_multimapping_status():
            if is_multimapped:
                return [1, 1] if is_paired_ended else 1
            else:
                return [0, 0] if is_paired_ended else 0


        def get(attribute: str):
            if is_paired_ended:
                result = [getattr(read, attribute) for read in self]
                if attribute == 'reference_start':
                    return [pos + 1 for pos in result]
                return result
            else:
                result = getattr(self[0], attribute)
                if attribute == 'reference_start':
                    result += 1
                return result


        new_row = {
            label + '_name': self[0].query_name,
            label + '_flags': get('flag'),
            label + '_pos': get('reference_start'),
            label + '_chr': get('reference_name'),
            label + '_CIGAR': get('cigarstring'),
            label + '_edit_dist': get_edit_dist(),
            label + '_quality': get('mapping_quality'),
            label + '_MD': get_tag(),
            label + '_multi': get_multimapping_status()
            }
        if is_real_data:
            new_row[label + '_sequence'] = self[0].query_sequence
        return new_row


class CSV:
    def __init__(self, path: str, label: str, is_real_data: bool):
        self.path = path
        self.label = label
        self.is_real_data = is_real_data
    
    def write_header(self):
        with open(self.path, 'w') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            columns = ['name', 'flags', 'pos', 'chr', 'CIGAR',
                       'edit_dist', 'quality', 'MD', 'multi']
            if self.is_real_data:
                columns.append('sequence')
            self.header = [f'{self.label}_{column}' for column in columns]
            writer.writerow(self.header)

    def append_rows(self, rows):
        with open(self.path, 'a') as csv_file:
            dictwriter_object = csv.DictWriter(csv_file, self.header)
            dictwriter_object.writerows(rows)

    def sort_csv(self):
        csvsort(self.path, [9], parallel=False)

    def drop_duplicates(self, column: int=9):
        """Takes SORTED csv file and deletes all duplicate rows inplace"""
        def format_output_path():
            SUFFIX = 'temp'
            head, tail = path.split(self.path)
            tails = tail.split('.')
            tails.insert(1, SUFFIX)
            tail = '.'.join(tails)
            output_path = path.join(head, tail)
            return output_path
        
        def delete_input_file():
            remove(self.path)
        
        def rename_temp_file():
            rename(format_output_path(), self.path)

        with open(self.path, 'r') as csv_master:
            with open(format_output_path(), 'w') as csv_filtered:
                reader = csv.reader(csv_master, delimiter=',')
                writer = csv.writer(csv_filtered, delimiter=',')
                
                header = next(reader)
                writer.writerow(header)

                prev_read = next(reader)
                prev_seq = prev_read[column]
                prev_prev_seq = None
                sequence_to_skip = None

                for read in reader:
                    if read[column] == prev_seq:
                        sequence_to_skip = prev_seq
                        continue
                    else:
                        if prev_seq == sequence_to_skip:
                            prev_prev_seq = prev_seq
                            prev_read = read
                            prev_seq = read[column]
                            continue
                        else:
                            writer.writerow(prev_read)
                            prev_prev_seq = prev_seq
                            prev_read = read
                            prev_seq = read[column]
                if prev_seq not in (sequence_to_skip, prev_prev_seq):
                    writer.writerow(prev_read)

        delete_input_file()
        rename_temp_file()


class Parser:

    def __init__(self, input_file: str,
                       label: str,
                       is_real_data: bool,
                       output_path: str,
                       discard_multimapped_reads: bool=True):
        self.input_file = input_file
        self.label = label
        self.is_real_data = is_real_data
        self.discard_multimapped_reads = discard_multimapped_reads
        self.is_paired_ended = self.is_paired_ended()
        self.max_non_miltimapped_reads = 2 if self.is_paired_ended else 1
        if self.discard_multimapped_reads:
            self.max_same_name_reads = self.max_non_miltimapped_reads
        else:
            self.max_same_name_reads = float('inf')
        self.output_path = self._parse_output_path(output_path)
        self.output_csv = CSV(self.output_path, label, is_real_data)
        self.output_csv.write_header()
        self.saved_reads_buffer = []
        self.BUFFER_SIZE = 1_000_000 # max number of reads in buffer
        self.all_reads_counter = 0
        self.unique_names_counter = 1
        self.multimapped_reads_counter = 0
        self.multimapped_reads_by_read_name_counter = 0
        self.unambiguous_reads_counter = 0
        self.primary_multimapped_counter = 0
    

    @staticmethod
    def _parse_output_path(output_path) -> str:
        """Generates output filepath in case it is not given or invalid"""
        DEFAULT_FILENAME = 'parser_output.csv'
        if not output_path:
            output_path = f'./{DEFAULT_FILENAME}'
        else:
            head, tail = path.split(output_path)
            if not tail:
                tail = DEFAULT_FILENAME
            if not tail.endswith('.csv'):
                tail += '.csv'
            output_path = path.join(head, tail)
        return output_path


    @staticmethod
    def get_current_memory_usage(in_gigabytes=True) -> float:
        """Checks the amount of RAM used by the script at the time"""
        process = psutil.Process(getpid())
        usage = process.memory_info().rss # in bytes
        usage_GB = usage * 9.31 * 10 ** (-10) # in Gigabytes
        return round(usage_GB, 3) if in_gigabytes else usage
    

    def is_paired_ended(self) -> bool:
        """Determines if data is paired ended by reading first 100 reads"""
        with pysam.AlignmentFile(self.input_file) as reads:
            paired_read_counter = 0
            for index, read in enumerate(reads):
                if read.is_paired:
                    paired_read_counter += 1
                if index == 100:
                    break
        return True if paired_read_counter > 0 else False


    def run(self):
        """Runs the script"""
        def filter_reads(reads: SameNameReadsArray):
            number_of_reads = len(reads)
            self.all_reads_counter += number_of_reads
            are_mm = number_of_reads > self.max_non_miltimapped_reads
            if are_mm:
                self.multimapped_reads_counter += number_of_reads
                self.multimapped_reads_by_read_name_counter += 1
            if number_of_reads > self.max_same_name_reads:
                return (SameNameReadsArray(), are_mm)
            return (reads.keep_only_primary_reads(), are_mm)
        
        def save_read(filtered_reads: SameNameReadsArray, are_mm: bool):
            new_row = filtered_reads.export(
                    label=self.label,
                    is_paired_ended=self.is_paired_ended,
                    is_multimapped=are_mm,
                    is_real_data=self.is_real_data
                    )
            if not are_mm:
                self.unambiguous_reads_counter += 1
            else:
                self.primary_multimapped_counter += 1
            self.saved_reads_buffer.append(new_row)
            if len(self.saved_reads_buffer) == self.BUFFER_SIZE:
                self.output_csv.append_rows(self.saved_reads_buffer)
                self.saved_reads_buffer = []
            
        with pysam.AlignmentFile(self.input_file, 'rb') as reads:
            first_line = next(reads)
            reads_same_name = SameNameReadsArray(first_line)
            for read in reads:
                if read.query_name == reads_same_name.name:
                    reads_same_name.append(read)
                else:
                    filtered_reads, are_mm = filter_reads(reads_same_name)
                    if len(filtered_reads) == 0:
                        reads_same_name = SameNameReadsArray(read)
                        self.unique_names_counter += 1
                        continue
                    save_read(filtered_reads, are_mm)
                    reads_same_name = SameNameReadsArray(read)
                    self.unique_names_counter += 1
            # process last read
            filtered_reads, are_mm = filter_reads(reads_same_name)
            if len(filtered_reads) != 0:
                save_read(filtered_reads, are_mm)
        # append remaining reads from buffer
        if len(self.saved_reads_buffer) != 0:
            self.output_csv.append_rows(self.saved_reads_buffer)
            self.saved_reads_buffer = []
            

class AlignmentFile:
    def __init__(self, path_to_file):
        self.input_file = path_to_file
        self.output_sorted_file = self._generate_output_path(path_to_file)


    def _generate_output_path(self, input_path: str) -> str:
        SUFFIX = 'sorted_by_name'
        head, tail = path.split(input_path)
        tails = tail.split('.')
        if len(tails) == 2:
            tails.insert(1, SUFFIX)
        elif len(tails) == 3:
            tails[1] = SUFFIX
        tail = '.'.join(tails)
        return path.join(head, tail)


    def sort_by_name(self) -> str:
        pysam.sort('-o', self.output_sorted_file, '-n', self.input_file)
        return self.output_sorted_file


class ResourceTracker:

    def __init__(self):
        self.times = [None, None, None]

    def start_tracking(self):
        for n, tracker in enumerate(self.times):
            if tracker == None:
                self.times[n] = time.time()
                break
    
    def end_tracking(self):
        for n, tracker in reversed(list(enumerate(self.times))):
            if tracker != None:
                time_format = "%M min %S sec"
                self.times[n] = time.strftime(
                    time_format, 
                    time.gmtime(time.time() - tracker)
                    ) 
                break
    
    def get_mem_usage(self) -> str:
        """Peak RAM usage, might be fair only for Linux systems"""
        memory_KB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        memory_MB = memory_KB / 1024
        memory_GB = memory_MB / 1024
        return f"{int(memory_MB):,} MB or {memory_GB:.2f} GB"


class MessageFormatter:

    @staticmethod
    def _write_heading(heading: str, body: str = "", new_line: bool = False):
        length = len(heading)
        equal_signs = ''.join(['=' for _ in range(length + 1)])
        print()
        if new_line:
            print(heading.upper())
            print(equal_signs)
            print(body)
        else:
            print(f"{heading.upper()}: {body}")
            print(equal_signs)

    @staticmethod
    def _wrap_between_asterisks(string: str, n_asterisks: int = 5):
        asterisks = "*" * n_asterisks
        print(f"\n{asterisks} {string} {asterisks}")
    
    def print_pre_sorting_message(self):
        self._wrap_between_asterisks('Sorting started...')

    def print_post_sorting_message(self):
        self._wrap_between_asterisks('Sorting complete')
    
    def print_pre_reading_message(self):
        self._wrap_between_asterisks('Reading input data started...')
    
    def print_post_reading_message(self, parser: Parser, input_file: Path):
        self._wrap_between_asterisks('Reading input data complete')
        self._write_heading('Input BAM File', input_file.absolute())
        self._write_heading('Output CSV File', 
                            Path(parser.output_path).absolute())
        self._write_heading(
            'Data Type',
            "Detected as PAIRED ENDED" if parser.is_paired_ended else
            "Detected as SINGLE ENDED"
            )
        self._write_heading('Label', parser.label)
        self._write_heading('Real Data', parser.is_real_data)
        self._write_heading(
            'Discard multimapped reads',
            parser.discard_multimapped_reads
            )
        self._wrap_between_asterisks('Parsing started...')
    
    def print_pre_real_data_processing_message(self):
        self._wrap_between_asterisks('Real data postprocessing started...')
    
    def print_post_real_data_processing_message(self):
        self._wrap_between_asterisks('Real data complete')

    def print_post_parsing_message(self, parser: Parser):
        self._wrap_between_asterisks('Parsing complete')
        is_PE = parser.is_paired_ended
        reads_per_name = 2 if is_PE else 1
        total = parser.all_reads_counter
        unique_names = parser.unique_names_counter
        mm = parser.multimapped_reads_counter
        mm_by_name = parser.multimapped_reads_by_read_name_counter

        unambiguous = parser.unambiguous_reads_counter
        unambiguous_singles = reads_per_name * unambiguous
        primary_mm = parser.primary_multimapped_counter
        primary_mm_singles = reads_per_name * primary_mm
        saved = unambiguous + primary_mm
        saved_singles = reads_per_name * saved

        discarded = total - saved * reads_per_name

        table = PrettyTable()
        table.field_names = ['Type of reads', 'Number of reads', 'Percentage']

        table.add_rows(
            [
                [
                    "Total reads in BAM File",
                    f"{total:,}",
                    ""
                ],
                [
                    "  - Multimapped total",
                    f"{mm:,} out of {total:,}",
                    f"{100*mm/total:.2f}%"
                ],
                [
                    "Unique read names",
                    f"{unique_names:,}",
                    ""
                ],
                [
                    "  - Multimapped read names",
                    f"{mm_by_name:,} out of {unique_names:,}",
                    f"{100*mm_by_name/unique_names:.2f}%"
                ],
                [
                    "Discarded reads",
                    f"{discarded:,} out of {total:,}",
                    f"{100*discarded/total:.2f}%"
                ],
                [
                    "Saved to CSV (reads)" if is_PE else "Saved to CSV",
                    f"{saved_singles:,} out of {total:,}",
                    f"{100*saved_singles/total:.2f}%"
                ],
                [
                    "  - Unambiguous",
                    f"{unambiguous_singles:,} out of {saved_singles:,}",
                    f"{100*unambiguous_singles/saved_singles:.2f}%"
                ],
                [
                    "  - Primary but multimapped",
                    f"{primary_mm_singles:,} out of {saved_singles:,}",
                    f"{100*primary_mm_singles/saved_singles:.2f}%"
                ]
            ]
        )
        if is_PE:
            table.add_rows(
                [
                    [
                        "Saved to CSV (pairs of reads)",
                        f"{saved:,}",
                        ""
                    ],
                    [
                        "  - Unambiguous",
                        f"{unambiguous:,} out of {saved:,}",
                        f"{100*unambiguous/saved:.2f}%"
                    ],
                    [
                        "  - Primary but multimapped",
                        f"{primary_mm:,} out of {saved:,}",
                        f"{100*primary_mm/saved:.2f}%"
                    ]
                ]
            )
        table.align['Type of reads'] = 'l'
        table.align['Number of reads'] = 'c'
        table.align['Percentage'] = 'r'
        print()
        print(table)

    def print_resourse_usage_message(self, tracker: ResourceTracker):
        self._write_heading(
            "Memory Usage", 
            f"Peak RAM consumption: {tracker.get_mem_usage()}",
            new_line=True
            )
        print("* Memory estimation might be fair only for Linux systems")
        self._write_heading("Elapsed Time")
        print(f"Sorting BAM execution time: {tracker.times[0]}")
        print(f"Parsing to CSV execution time: {tracker.times[1]}")
        if tracker.times[2] != None:
            print(
                "Postprocessing real data execution time:", 
                tracker.times[2]
                )


app = typer.Typer(add_completion=False)
@app.command()
def main(input_file: Path = typer.Argument(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            show_default=False
            ),
         output_path: Path = typer.Option(
            "./parser_output.csv",
            "--output",
            "-o",
            help="Path to output directory or file."
            ),
         label: str = typer.Option(
            "NoLabel",
            "--label", "-l",
            help=(
                "Prefix for each column name in the output CSV file. " + 
                "For further comparison it is " +
                "recommended to specify different labels for each sample."
                )
            ),
         is_real_data: bool = typer.Option(
            False,
            "--is-real-data", "-r", 
            show_default=False,
            is_flag=True,
            help = (
                "Prepare CSV for further subsampling. " +
                "Only for comparison samples of different size."
                )
            ),
         keep_multimapped_reads: bool = typer.Option(
            False, 
            "--keep-multimapped_reads", "-m", 
            show_default=True,
            is_flag=True,
            help=(
                "Do not discard primary multimapped reads. " +
                "If flag not given, only unambiguous reads will be " +
                "written to the CSV.")
            )
            ):
    """
    Python script that parses Binary Alignment Map File,
    keeping only primary alignments.
    
    """
    messager = MessageFormatter()
    tracker = ResourceTracker()

    # 1) Sort by read name
    messager.print_pre_sorting_message()
    tracker.start_tracking()

    alignment_file = AlignmentFile(str(input_file))
    sorted_alignment_file = alignment_file.sort_by_name()

    tracker.end_tracking()
    messager.print_post_sorting_message()

    # 2) Parse BAM file to CSV file read by read
    messager.print_pre_reading_message()
    tracker.start_tracking()

    pysam.set_verbosity(0) # set minimum verbosity to silence index warning
    parser = Parser(
        input_file = sorted_alignment_file,
        label = label,
        is_real_data = is_real_data,
        output_path = str(output_path),
        discard_multimapped_reads = not keep_multimapped_reads)
    
    messager.print_post_reading_message(parser, input_file)

    parser.run()
    generated_csv_path = parser.output_path

    messager.print_post_parsing_message(parser)
    tracker.end_tracking()

    # 3) Process real data (optional)
    if is_real_data:
        messager.print_pre_real_data_processing_message()
        tracker.start_tracking()

        generated_csv = CSV(generated_csv_path, label, is_real_data)
        generated_csv.sort_csv()
        generated_csv.drop_duplicates()

        tracker.end_tracking()
        tracker.get_peak_memory_usage()
        messager.print_post_real_data_processing_message()

    messager.print_resourse_usage_message(tracker)

if __name__ == "__main__":
    app()