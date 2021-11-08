import pysam
import tqdm
import csv
import functools
import collections
import pandas as pd

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


if __name__ == "__main__":
    main()
