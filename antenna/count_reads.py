#!/usr/bin/env python3

import pysam
import tqdm
import csv
import functools
import collections
import pandas as pd
import argparse
import matplotlib.pyplot as plt


@functools.total_ordering
class BedInterval:
    """Class representing a bed interval"""

    # Note: using total_ordering has a performance penalty

    def __init__(self, ref, start, end, name):
        self.ref = ref
        self.start = start
        self.end = end
        self.name = name
        self.properties = collections.defaultdict(int)

    def _is_valid_operand(self, other):
        return isinstance(other, BedInterval)

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


class BedIntervals:
    def __init__(self):
        self.intervals = []

    def load_bed(self, bedfilename):
        with open(bedfilename) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            for line in tsv_file:
                ref, start, end, name, *other = line
                start = int(start)
                end = int(end)
                interval = BedInterval(ref, start, end, name)
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
        "query_name": [],
        "reference_start": [],
        "reference_end": [],
        "reference_name": [],
        "is_read1": [],
        "reverse": [],
        "p5_o_score": [],
        "p5_rc_score": [],
        "p5_r_score": [],
        "p5_c_score": [],
        "p3_o_score": [],
        "p3_rc_score": [],
        "p3_r_score": [],
        "p3_c_score": [],
    }

    with pysam.AlignmentFile(bamfilename, "rb", require_index=True) as bamfile:
        total_input_reads = bamfile.mapped + bamfile.unmapped
        for read in tqdm.tqdm(bamfile, total=total_input_reads, smoothing=0):
            if read.seq is None or read.is_unmapped:
                continue

            try:
                scores = read.get_tag(sgRNA_bam_score_tag_name)
                scores_array = scores.split(",")

                df_data["query_name"].append(read.query_name)
                df_data["reference_start"].append(read.reference_start)
                df_data["reference_end"].append(read.reference_end)
                df_data["reference_name"].append(read.reference_name)

                df_data["is_read1"].append(read.is_read1)
                df_data["reverse"].append(read.is_reverse)

                df_data["p5_o_score"].append(scores_array[0])
                df_data["p5_rc_score"].append(scores_array[1])
                df_data["p5_r_score"].append(scores_array[2])
                df_data["p5_c_score"].append(scores_array[3])

                df_data["p3_o_score"].append(scores_array[4])
                df_data["p3_rc_score"].append(scores_array[5])
                df_data["p3_r_score"].append(scores_array[6])
                df_data["p3_c_score"].append(scores_array[7])
            except KeyError:
                continue

        df = pd.DataFrame(df_data)
        df = df.astype(
            {
                "query_name": "U",
                "reference_start": "int32",
                "reference_end": "int32",
                "reference_name": "U",
                "is_read1": "bool",
                "p5_o_score": "int32",
                "p5_rc_score": "int32",
                "p5_r_score": "int32",
                "p5_c_score": "int32",
                "p3_o_score": "int32",
                "p3_rc_score": "int32",
                "p3_r_score": "int32",
                "p3_c_score": "int32",
            }
        )

        return df


def summarize_trs_intervals(trs_intervals):
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


def merge_read_information(scores_df):
    """Merge read information for paired reads"""

    scores_df_r1 = scores_df[scores_df["is_read1"] == True]
    scores_df_r2 = scores_df[scores_df["is_read1"] == False]
    scores_df_r1 = scores_df_r1.set_index("query_name")
    scores_df_r2 = scores_df_r2.set_index("query_name")
    scores_read_pairs = scores_df_r1.join(scores_df_r2, lsuffix="_r1", rsuffix="_r2",)

    return scores_read_pairs


def count_sgRNA(merge_sg_read_info, bedfile, cutoff, sgRNA_bam_tag_name="TO"):
    """Count sgRNAs stratifying by orientation"""
    trs_intervals = BedIntervals()
    trs_intervals.load_bed(bedfile)

    columns_to_count = [
        "p5_o_score_r1",
        "p3_o_score_r1",
        "p5_o_score_r2",
        "p3_o_score_r2",
    ]

    for read_pair in merge_sg_read_info.itertuples():
        count_read_pair = False

        for cc in columns_to_count:
            if getattr(read_pair, cc) > cutoff:
                count_read_pair = True
                continue  # This should only jump out of inner

       ## check that I am looking at reads that point inwards and remove chimeric ones
    
       # consider finding insert sizes here and sd
       # 
       #     min(start_r1, end_r1, start_r2, end_r2)
       #     max(start_r1, end_r1, start_r2, end_r2)
       #     check if within range as valid pair
       #     check if whole thing overlaps the known interval
            
        # TODO: consider r1 orientation (that also dictates r2 orientation)
        read_interval = trs_intervals.get_overlapping_interval(
            read_pair.reference_start_r1
        )
        if read_interval:
            if count_read_pair:
                read_interval["sg_count"] += 1
            else:
                read_interval["non_sg_count"] += 1

    return trs_intervals


def plot_read_score_dist(scores_df, max_range = 80, n_bins = 50):
    """Return a plot of TRS alignment score distributions by read orientation, read type, flanking side and TRS orientation"""

    scores_subset = {
        "r1": {
            "fd": scores_df[
                (scores_df["is_read1"] == True) & (scores_df["reverse"] == False)
            ],
            "rv": scores_df[
                (scores_df["is_read1"] == False) & (scores_df["reverse"] == False)
            ],
        },
        "r2": {
            "fd": scores_df[
                (scores_df["is_read1"] == True) & (scores_df["reverse"] == True)
            ],
            "rv": scores_df[
                (scores_df["is_read1"] == False) & (scores_df["reverse"] == True)
            ],
        },
    }

    abbr = {
        "fd": "Forward",
        "rv": "Reverse",
        "r1": "R1",
        "r2": "R2",
        "p5": "5'",
        "p3": "3'",
        "o": "TRS",
        "rc": "RC TRS",
        "r": "R TRS",
        "c": "C TRS",
    }

    nrows = 8
    ncols = 4

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(10, 15),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    
    

    i = 0
    for read_type in ("fd", "rv"):
        for read_order in ("r1", "r2"):
            for flanking_side in ("p5", "p3"):
                for orientation in ("o", "rc", "r", "c"):
                    col_index = i // nrows
                    row_index = i % nrows
                    value = scores_subset[read_order][read_type][
                        f"{flanking_side}_{orientation}_score"
                    ]
                    axes[row_index, col_index].hist(value, log=True, range=(0,max_range), bins=n_bins)
                    axes[row_index, col_index].set_title(
                        f"{abbr[read_type]} {abbr[read_order]} {abbr[flanking_side]} {abbr[orientation]}"
                    )
                    axes[row_index, col_index].set_xlabel("score")
                    i += 1

    return fig, axes


def main():
    argparser = argparse.ArgumentParser("antena_count_reads")
    argparser.add_argument("--bam", help="input bam file", required=True)
    argparser.add_argument("--bed", help="bed file with intervals", required=True)
    argparser.add_argument("--outcsv", help="output csv file", required=True)
    argparser.add_argument("--cutoff", help="score cutoff", default=50)
    argparser.add_argument(
        "--progress", help="display progress bar", action="store_true", default=False
    )
    args = argparser.parse_args()

    # Main execution flow
    # TODO: Check bed file ok before starting the counting
    sgRNA_scores = load_sgRNA_scores(args.bam)
    merged_read_information = merge_read_information(sgRNA_scores)
    trs_intervals = count_sgRNA(merged_read_information, args.bed, args.cutoff)
    intervals_counts = summarize_trs_intervals(trs_intervals)
    intervals_counts.to_csv(args.outcsv)


if __name__ == "__main__":
    main()
