import pysam
import tqdm
import csv
import functools
import collections
import pandas as pd

@functools.total_ordering
class bed_interval:
    """Class representing a bed interval"""
    # Note: using total_ordering has a performance penalty
    
    def __init__(self, ref, start, end, name):
        self.ref = ref
        self.start = start
        self.end = end
        self.name = name
        self.properties = collections.defaultdict(int)
    
    def _is_valid_operand(self, other):
        return isinstance(other, bed_interval)
        
    def __eq__(self, other):
        if not (self._is_valid_operand(other)):
            return NotImplemented
        return self.start == other.start
        
    def __lt__(self, other):
        if not (self._is_valid_operand(other)):
            return NotImplemented
        return self.start < other.start
    
    def __setitem__(self, k, v):
        self.properties[k] = v
    
    def __getitem__(self, k):
        return self.properties[k]
    
    def contains(self, pos):
        if pos >= self.start and pos <= self.end:
            return True
        else:
            return False


class bed_intervals:
    def __init__(self):
        self.intervals = []
    
    def load_bed(self, bedfilename):
        with open(bedfilename) as file:
            tsv_file = csv.reader(file, delimiter = '\t')
            for line in tsv_file:
                ref, start, end, name, *other = line
                start = int(start)
                end = int(end)
                interval = bed_interval(ref, start, end, name)
                self.intervals.append(interval)
        self.intervals.sort()
        
    def find_overlapping_interval(self, pos):
        """Find the index of the overlapping interval"""
        ## This is very inefficient, for a small bed file it doesn't matter
        for (i, cur_interval) in enumerate(self.intervals):
            if cur_interval.contains(pos):
                return i
        return None
    
    def get_overlapping_interval(self, pos):
        """Get the overlapping interval or none"""
        i = self.find_overlapping_interval(pos)
        if i is not None:
            return self.intervals[i]
        return None
    
def simple_count_sgRNA(bamfilename, bedfile, sgRNA_bam_tag_name = 'TS'):
    """
    Count sgRNA reads in the most simple way possible
    
    Consider first reads only to avoid dealing with read pairs
    In the defined windows get counts for TRS containing and non-trs containisn
    ignore orientations
    
    """
    
    trs_intervals = bed_intervals()
    trs_intervals.load_bed(bedfile)
    
    with pysam.AlignmentFile(bamfilename, "rb", require_index=True) as bamfile:
        total_input_reads = bamfile.mapped + bamfile.unmapped
        for read in tqdm.tqdm(bamfile, total=total_input_reads, smoothing=0):
            # Ignore reads without a sequence or that are unmapped
            if (
                read.seq is None or
                read.is_unmapped
            ):
                continue
            if (read.is_read1):
                read_interval= trs_intervals.get_overlapping_interval(read.reference_start)
                if read_interval is not None:
                    if read.has_tag(sgRNA_bam_tag_name):
                        read_interval['sg_count'] += 1
                    else:
                        read_interval['non_sg_count'] += 1
                else:
                    continue
            else:
                # Ignore non-read 1
                continue
    
    # Summarize to pandas df
    names = []
    refs = []
    starts = []
    ends = []
    sg_counts = []
    non_sg_counts = []
    
    for interval in trs_intervals.intervals:
        names.append(interval.name)
        refs.append(interval.ref)
        starts.append(interval.start)
        ends.append(interval.end)
        sg_counts.append(interval['sg_count'])
        non_sg_counts.append(interval['non_sg_count'])
                
    return pd.DataFrame({'name': names, 
                         'ref': refs,
                         'start': starts,
                         'end': ends,
                         'sg_count': sg_counts, 
                         'non_sg_counts': non_sg_counts})
    
