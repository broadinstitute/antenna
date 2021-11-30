#!/usr/bin/env python3

import pysam
import tqdm
import csv
import functools
import collections
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Masks for output annotation
TRS_3_PRIME_RC = 0x1 << 7
TRS_3_PRIME_C = 0x1 << 6
TRS_3_PRIME_R = 0x1 << 5
TRS_3_PRIME_O = 0x1 << 4

TRS_5_PRIME_RC = 0x1 << 3
TRS_5_PRIME_C = 0x1 << 2
TRS_5_PRIME_R = 0x1 << 1
TRS_5_PRIME_O = 0x1 << 0


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
            tsv_file = csv.reader(file, delimiter="\t")
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
        ## The list is already sorted for binary search
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

def load_sgRNA_scores(bamfilename, sgRNA_bam_score_tag_name="TO"):
    """Load the read scores in all possible orientations from a bam file"""
    
    df_data = {
        'read_name': [],
        'is_read1': [],
        'reverse': [],
        'p5_o_score': [],
        'p5_rc_score': [],
        'p5_r_score': [],
        'p5_c_score': [],
        'p3_o_score': [],
        'p3_rc_score': [],
        'p3_r_score': [],
        'p3_c_score': [],
    }
    
    with pysam.AlignmentFile(bamfilename, "rb", require_index=True) as bamfile:
        total_input_reads = bamfile.mapped + bamfile.unmapped
        for read in tqdm.tqdm(bamfile, total=total_input_reads, smoothing=0):
            if read.seq is None or read.is_unmapped:
                continue
            
            try:
                scores = read.get_tag(sgRNA_bam_score_tag_name)
                scores_array = scores.split(',')
                df_data['read_name'].append(read.query_name)
                df_data['is_read1'].append(read.is_read1) 
                df_data['reverse'].append(read.is_reverse)
                
                df_data['p5_o_score'].append(scores_array[0])
                df_data['p5_rc_score'].append(scores_array[1])
                df_data['p5_r_score'].append(scores_array[2])
                df_data['p5_c_score'].append(scores_array[3])
                
                df_data['p3_o_score'].append(scores_array[4])
                df_data['p3_rc_score'].append(scores_array[5])
                df_data['p3_r_score'].append(scores_array[6])
                df_data['p3_c_score'].append(scores_array[7])
            except KeyError:
                continue
            
        df =  pd.DataFrame(df_data)
        df = df.astype({
            'read_name': 'U', 'is_read1': 'bool', 
            'p5_o_score': 'int32', 'p5_rc_score': 'int32', 'p5_r_score': 'int32', 'p5_c_score': 'int32',
            'p3_o_score': 'int32', 'p3_rc_score': 'int32', 'p3_r_score': 'int32', 'p3_c_score': 'int32',
        })
        
        return df
            
    

def count_sgRNA(bamfilename, bedfile, sgRNA_bam_tag_name="TS"):
    """Count sgRNAs stratifying by orientation"""
    trs_intervals = bed_intervals()
    trs_intervals.load_bed(bedfile)

    def count_read(read, trs_intervals, sgRNA_bam_tag_name, double_count=True):
        """Count each read, note that we allow double counting here"""

        # TODO: examine shorter version of this function -- consider speed
        read_interval = trs_intervals.get_overlapping_interval(read.reference_start)
        if read_interval is not None:
            if read.has_tag(sgRNA_bam_tag_name):
                trs_tag_info = read.get_tag(sgRNA_bam_tag_name)
                read_interval["sg_count"] += 1
                if read.is_read1:
                    if trs_tag_info & TRS_5_PRIME_O:
                        read_interval["sg_count_read1_5p_o"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_R:
                        read_interval["sg_count_read1_5p_r"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_C:
                        read_interval["sg_count_read1_5p_c"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_RC:
                        read_interval["sg_count_read1_5p_rc"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_O:
                        read_interval["sg_count_read1_3p_o"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_R:
                        read_interval["sg_count_read1_3p_r"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_C:
                        read_interval["sg_count_read1_3p_c"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_RC:
                        read_interval["sg_count_read1_3p_rc"] += 1
                        if not double_count:
                            return
                if read.is_read2:
                    if trs_tag_info & TRS_5_PRIME_O:
                        read_interval["sg_count_read2_5p_o"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_R:
                        read_interval["sg_count_read2_5p_r"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_C:
                        read_interval["sg_count_read2_5p_c"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_5_PRIME_RC:
                        read_interval["sg_count_read2_5p_rc"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_O:
                        read_interval["sg_count_read2_3p_o"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_R:
                        read_interval["sg_count_read2_3p_r"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_C:
                        read_interval["sg_count_read2_3p_c"] += 1
                        if not double_count:
                            return
                    if trs_tag_info & TRS_3_PRIME_RC:
                        read_interval["sg_count_read2_3p_rc"] += 1
                        if not double_count:
                            return
            else:
                read_interval["non_sg_count"] += 1

    with pysam.AlignmentFile(bamfilename, "rb", require_index=True) as bamfile:
        total_input_reads = bamfile.mapped + bamfile.unmapped
        for read in tqdm.tqdm(bamfile, total=total_input_reads, smoothing=0):
            if read.seq is None or read.is_unmapped:
                continue
            count_read(read, trs_intervals, sgRNA_bam_tag_name)

    return trs_intervals


def summarize_trs_intervals(trs_intervals):
    # Summarize to pandas df
    names = []
    refs = []
    starts = []
    ends = []
    sg_counts = []
    sg_count_read1_5p_o_l = []
    sg_count_read1_5p_r_l = []
    sg_count_read1_5p_c_l = []
    sg_count_read1_5p_rc_l = []
    sg_count_read1_3p_o_l = []
    sg_count_read1_3p_r_l = []
    sg_count_read1_3p_c_l = []
    sg_count_read1_3p_rc_l = []
    sg_count_read2_5p_o_l = []
    sg_count_read2_5p_r_l = []
    sg_count_read2_5p_c_l = []
    sg_count_read2_5p_rc_l = []
    sg_count_read2_3p_o_l = []
    sg_count_read2_3p_r_l = []
    sg_count_read2_3p_c_l = []
    sg_count_read2_3p_rc_l = []
    non_sg_counts = []
    for interval in trs_intervals.intervals:
        names.append(interval.name)
        refs.append(interval.ref)
        starts.append(interval.start)
        ends.append(interval.end)
        sg_counts.append(interval["sg_count"])
        sg_count_read1_5p_o_l.append(interval["sg_count_read1_5p_o"])
        sg_count_read1_5p_r_l.append(interval["sg_count_read1_5p_r"])
        sg_count_read1_5p_c_l.append(interval["sg_count_read1_5p_c"])
        sg_count_read1_5p_rc_l.append(interval["sg_count_read1_5p_rc"])
        sg_count_read1_3p_o_l.append(interval["sg_count_read1_3p_o"])
        sg_count_read1_3p_r_l.append(interval["sg_count_read1_3p_r"])
        sg_count_read1_3p_c_l.append(interval["sg_count_read1_3p_c"])
        sg_count_read1_3p_rc_l.append(interval["sg_count_read1_3p_rc"])
        sg_count_read2_5p_o_l.append(interval["sg_count_read2_5p_o"])
        sg_count_read2_5p_r_l.append(interval["sg_count_read2_5p_r"])
        sg_count_read2_5p_c_l.append(interval["sg_count_read2_5p_c"])
        sg_count_read2_5p_rc_l.append(interval["sg_count_read2_5p_rc"])
        sg_count_read2_3p_o_l.append(interval["sg_count_read2_3p_o"])
        sg_count_read2_3p_r_l.append(interval["sg_count_read2_3p_r"])
        sg_count_read2_3p_c_l.append(interval["sg_count_read2_3p_c"])
        sg_count_read2_3p_rc_l.append(interval["sg_count_read2_3p_rc"])
        non_sg_counts.append(interval["non_sg_count"])

    return pd.DataFrame(
        {
            "name": names,
            "ref": refs,
            "start": starts,
            "end": ends,
            "sg_count": sg_counts,
            "sg_count_read1_5p_o": sg_count_read1_5p_o_l,
            "sg_count_read1_5p_r": sg_count_read1_5p_r_l,
            "sg_count_read1_5p_c": sg_count_read1_5p_c_l,
            "sg_count_read1_5p_rc": sg_count_read1_5p_rc_l,
            "sg_count_read1_3p_o": sg_count_read1_3p_o_l,
            "sg_count_read1_3p_r": sg_count_read1_3p_r_l,
            "sg_count_read1_3p_c": sg_count_read1_3p_c_l,
            "sg_count_read1_3p_rc": sg_count_read1_3p_rc_l,
            "sg_count_read2_5p_o": sg_count_read2_5p_o_l,
            "sg_count_read2_5p_r": sg_count_read2_5p_r_l,
            "sg_count_read2_5p_c": sg_count_read2_5p_c_l,
            "sg_count_read2_5p_rc": sg_count_read2_5p_rc_l,
            "sg_count_read2_3p_o": sg_count_read2_3p_o_l,
            "sg_count_read2_3p_r": sg_count_read2_3p_r_l,
            "sg_count_read2_3p_c": sg_count_read2_3p_c_l,
            "sg_count_read2_3p_rc": sg_count_read2_3p_rc_l,
            "non_sg_counts": non_sg_counts,
        }
    )


def simple_count_sgRNA(
    bamfilename, bedfile, sgRNA_bam_tag_name="TS", count_read_1=True, count_read_2=False
):
    """
    Count sgRNA reads in the most simple way possible
    
    Consider first reads only to avoid dealing with read pairs
    In the defined windows get counts for TRS containing and non-trs containisn
    ignore orientations
    
    """

    trs_intervals = bed_intervals()
    trs_intervals.load_bed(bedfile)

    def count_read(read, trs_intervals, sgRNA_bam_tag_name):
        read_interval = trs_intervals.get_overlapping_interval(read.reference_start)
        if read_interval is not None:
            if read.has_tag(sgRNA_bam_tag_name):
                read_interval["sg_count"] += 1
            else:
                read_interval["non_sg_count"] += 1

    with pysam.AlignmentFile(bamfilename, "rb", require_index=True) as bamfile:
        total_input_reads = bamfile.mapped + bamfile.unmapped
        for read in tqdm.tqdm(bamfile, total=total_input_reads, smoothing=0):
            # Ignore reads without a sequence or that are unmapped
            if read.seq is None or read.is_unmapped:
                continue
            if read.is_read1 and count_read_1:
                count_read(read, trs_intervals, sgRNA_bam_tag_name)
            if read.is_read2 and count_read_2:
                count_read(read, trs_intervals, sgRNA_bam_tag_name)

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
        sg_counts.append(interval["sg_count"])
        non_sg_counts.append(interval["non_sg_count"])

    return pd.DataFrame(
        {
            "name": names,
            "ref": refs,
            "start": starts,
            "end": ends,
            "sg_count": sg_counts,
            "non_sg_counts": non_sg_counts,
        }
    )


def main():
    argparser = argparse.ArgumentParser("antena_count_reads")
    argparser.add_argument("--bam", help="input bam file", required=True)
    argparser.add_argument("--bed", help="bed file with intervals", required=True)
    argparser.add_argument("--outcsv", help="output csv file", required=True)
    argparser.add_argument(
        "--progress", help="display progress bar", action="store_true", default=False
    )
    args = argparser.parse_args()

    # Main execution flow
    trs_intervals = count_sgRNA(args.bam, args.bed)
    intervals_counts = summarize_trs_intervals(trs_intervals)
    intervals_counts.to_csv(args.outcsv)
    
    
def plot_read_score_dist(scores_df):
    """Return a plot of TRS alignment score distributions by read orientation, read type, flanking side and TRS orientation"""
    
    scores_subset = {
        'r1': {
            'fd': scores_df[(scores_df['is_read1'] == True) & (scores_df['reverse'] == False)],
            'rv': scores_df[(scores_df['is_read1'] == False) & (scores_df['reverse'] == False)]
        },
        'r2': {
            'fd': scores_df[(scores_df['is_read1'] == True) & (scores_df['reverse'] == True)],
            'rv': scores_df[(scores_df['is_read1'] == False) & (scores_df['reverse'] == True)]
        }
    }
    
    abbr = {
        'fd': 'Forward',
        'rv': 'Reverse',
        'r1': 'R1',
        'r2': 'R2',
        'p5': '5\'',
        'p3': '3\'',
        'o': 'TRS',
        'rc': 'RC TRS',
        'r': 'R TRS',
        'c': 'C TRS',
    }
    
    nrows = 8
    ncols = 4

    fig, axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=(10,15), sharex=True, sharey=True, constrained_layout=True)
    
    i=0
    for read_type in ('fd','rv'):
        for read_order in ('r1','r2'):
            for flanking_side in ('p5','p3'):
                for orientation in ('o','rc','r','c'):
                    col_index = i//nrows
                    row_index = i%nrows
                    value = scores_subset[read_order][read_type][f'{flanking_side}_{orientation}_score']
                    axes[row_index,col_index].hist(value, log=True)
                    axes[row_index, col_index].set_title(f'{abbr[read_type]} {abbr[read_order]} {abbr[flanking_side]} {abbr[orientation]}')
                    axes[row_index, col_index].set_xlabel('score')
                    i+=1
            
    return fig
    


if __name__ == "__main__":
    main()
