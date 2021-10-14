#!/usr/bin/env python3

import pysam
from Bio import pairwise2
from pybedtools import *
import numpy as np
import tqdm
import argparse
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import logging
from collections import defaultdict

# CONSTS
BAM_CSOFT_CLIP = 4


class ClassifiedRead:
    """Class representing a read that has been classified"""

    def __init__(
        self, sgRNA: bool, orf: str, read: pysam.AlignedRead, motif_orientation: int
    ):
        self.sgRNA = sgRNA
        self.orf = orf
        self.pos = read.pos
        self.motif_orientation = motif_orientation
        self.read_orientation = "reverse" if read.is_reverse else "forward"
        self.first_in_pair = "r1" if read.is_read1 else "r2"


def get_mapped_reads(bam):
    """Get total mapped reads from the bam file index"""
    mapped_reads = (
        int(pysam.flagstat(bam).split("\n")[4].split(" ")[0])
        - int(pysam.flagstat(bam).split("\n")[2].split(" ")[0])
        - int(pysam.flagstat(bam).split("\n")[1].split(" ")[0])
    )
    return mapped_reads


def get_coverage(start, end, inbamfile):
    """Get coverage in a specified region of a bam file"""
    coverage = []
    if start < 0:
        start = 1
    for pileupcolumn in inbamfile.pileup("MN908947.3", int(start), int(end)):
        coverage.append(pileupcolumn.n)
    return np.median(coverage)


def check_trs_alignment(bases_clipped, TRS_sequence):
    """Get the alignment of the TRS against the clipped region"""
    alignments = pairwise2.align.localms(
        bases_clipped, TRS_sequence, 2, -2, -20, -0.1, one_alignment_only=True
    )
    return alignments[0].score if len(alignments) > 0 else 0


def check_alignment_factory(TRS_sequence):
    TRS_sequence = Seq(TRS_sequence)
    TRS_sequence_rc = TRS_sequence.reverse_complement()
    TRS_sequence_c = TRS_sequence.complement()
    TRS_sequence_r = TRS_sequence.reverse_complement().complement()

    if __debug__:
        logging.debug(f"TRS {TRS_sequence}")
        logging.debug(f"TRS RC: {TRS_sequence_rc}")
        logging.debug(f"TRS C: {TRS_sequence_c}")
        logging.debug(f"TRS R: {TRS_sequence_r}")

    def check_alignment(bases_clipped, score_cutoff):
        subgenomic_read = False
        orientation = None

        if check_trs_alignment(bases_clipped, TRS_sequence) >= score_cutoff:
            subgenomic_read = True
            orientation = 1
        elif check_trs_alignment(bases_clipped, TRS_sequence_rc) >= score_cutoff:
            subgenomic_read = True
            orientation = 2
        elif check_trs_alignment(bases_clipped, TRS_sequence_r) >= score_cutoff:
            subgenomic_read = True
            orientation = 3
        elif check_trs_alignment(bases_clipped, TRS_sequence_c) >= score_cutoff:
            subgenomic_read = True
            orientation = 4
        return (subgenomic_read, orientation)

    return check_alignment


def count_reads(reads):
    """Convert reads into ORF counts"""
    if __debug__:
        logging.debug("Counting reads...")

    orfs = defaultdict(int)

    for id, pair in reads.items():
        if __debug__:
            logging.debug(f"Processing read query set {id} with {len(pair)} reads")
        if len(pair) == 2:
            # Process pair of reads
            left_read = min(pair, key=lambda x: x.pos)
            right_read = max(pair, key=lambda x: x.pos)
            if __debug__:
                logging.debug(
                    f"left_read orf: {left_read.orf} sgRNA: {left_read.sgRNA}"
                )
                logging.debug(
                    f"right_read orf: {right_read.orf} sgRNA: {right_read.sgRNA}"
                )
            if left_read.sgRNA:
                orfs[
                    (
                        left_read.orf,
                        left_read.read_orientation,
                        left_read.first_in_pair,
                        left_read.motif_orientation,
                    )
                ] += 1
        else:
            # Process a single read
            read = pair[0]
            if __debug__:
                logging.debug(f"read orf: {read.orf} sgRNA: {read.orf}")
            if read.sgRNA:
                orfs[
                    (
                        read.orf,
                        read.read_orientation,
                        read.first_in_pair,
                        read.motif_orientation,
                    )
                ] += 1
    return orfs


def run_antenna(
    orf_bed_filename,
    inbam_filename,
    TRS_sequence,
    score_cutoff=50,
    n_clipped_cutoff=6,
    n_clipped_overhang=3,
    process_3prime_clipped=True
):
    """The main entry point for the antenna algorithm"""

    check_alignment = check_alignment_factory(TRS_sequence)

    orf_bed_object = BedTool(orf_bed_filename)
    reads = defaultdict(list)

    with pysam.AlignmentFile(inbam_filename, "rb", require_index=True) as inbamfile:
        # Process Reads
        total_input_reads = inbamfile.mapped + inbamfile.unmapped
        for read in tqdm.tqdm(inbamfile, total=total_input_reads, smoothing=0):
            subgenomic_read = False
            orientation = 0

            if (
                read.seq == None
                or read.is_unmapped
                or read.is_supplementary
                or read.is_secondary
            ):
                if __debug__:
                    logging.debug("Skipping read")
                continue

            cigar = read.cigartuples

            if cigar[0][0] == BAM_CSOFT_CLIP:
                n_clipped = cigar[0][1]
                if n_clipped > n_clipped_cutoff:
                    bases_clipped = read.seq[0 : (n_clipped + n_clipped_overhang)]

                    subgenomic_read, orientation = check_alignment(
                        bases_clipped, score_cutoff
                    )

                    if subgenomic_read:
                        read_orf = None
                        for row in orf_bed_object:
                            if row.end >= read.reference_start >= row.start:
                                read_orf = row.name

                        if read_orf == None:
                            read_orf = "novel_" + str(read.reference_start)

                        reads[read.query_name].append(
                            ClassifiedRead(
                                sgRNA=subgenomic_read,
                                orf=read_orf,
                                read=read,
                                motif_orientation=orientation,
                            )
                        )

            if process_3prime_clipped:
                if (not subgenomic_read) and cigar[-1][0] == BAM_CSOFT_CLIP:
                    n_clipped = cigar[-1][1]
                    if n_clipped > n_clipped_cutoff:
                        read_length = read.template_length
                        bases_clipped = read.seq[
                            read_length - n_clipped - n_clipped_overhang : read_length
                        ]

                        subgenomic_read, orientation = check_alignment(
                            bases_clipped, score_cutoff
                        )

                        if subgenomic_read:
                            read_orf = None
                            for row in orf_bed_object:
                                if row.end >= read.reference_start >= row.start:
                                    read_orf = row.name

                                if read_orf == None:
                                    read_orf = "novel_" + str(read.reference_start)

                                reads[read.query_name].append(
                                    ClassifiedRead(
                                        sgRNA=subgenomic_read,
                                        orf=read_orf,
                                        read=read,
                                        motif_orientation=orientation,
                                    )
                                )

        orfs = count_reads(reads)

    #        # Get coverage of different orfs
    #        orf_coverage = defaultdict(int)
    #        for row in orf_bed_object:
    #            orf_coverage[row.name] = get_coverage(row.start, row.end, inbamfile)

    return orfs


def save_output(orfs, output_filename):
    """Save the sgRNA counts as csv"""

    with open(output_filename, "w") as output_file:
        output_file.write(
            ",".join(["orf", "read_orientation", "first_in_pair", "motif_orientation", "sgRNA counts", "\n"])
        )
        for orf, read_orientation, first_in_pair, motif_orientation in orfs:
            count = orfs[(orf, read_orientation, first_in_pair, motif_orientation)]
            output_file.write(
                ",".join(
                    [
                        orf,
                        read_orientation,
                        str(first_in_pair),
                        str(motif_orientation),
                        str(count),
                        "\n",
                    ]
                )
            )


def main():
    argparser = argparse.ArgumentParser("antenna: a tool for sgRNA read search")
    argparser.add_argument("--bam", help="input bam file", required=True)
    argparser.add_argument("--output-counts", help="output csv file", required=True)
    argparser.add_argument(
        "--orf-bed", help="bed file with orf starting positions", required=True
    )
    argparser.add_argument("--progress", help="display progress bar")
    argparser.add_argument(
        "--trs-sequence",
        help="TRS sequence to search",
        default="AACCAACTTTCGATCTCTTGTAGATCTGTTCTC",
    )

    args = argparser.parse_args()

    if __debug__:
        logging.basicConfig(level=logging.DEBUG)

    count_data = run_antenna(
        orf_bed_filename=args.orf_bed,
        inbam_filename=args.bam,
        TRS_sequence=args.trs_sequence,
    )
    save_output(count_data, output_filename=args.output_counts)


if __name__ == "__main__":
    main()
