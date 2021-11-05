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


def check_trs_alignment(bases_clipped, TRS_sequence):
    """Get the alignment of the TRS against the clipped region"""
    alignments = pairwise2.align.localms(
        bases_clipped, TRS_sequence, 2, -2, -20, -0.1, one_alignment_only=True
    )
    return alignments[0].score if len(alignments) > 0 else 0


def check_alignment_factory(TRS_sequence):
    """ Function factory for checking alignment of TRS sequence in different orientations"""
    TRS_sequence = Seq(TRS_sequence)
    TRS_sequence_rc = TRS_sequence.reverse_complement()
    TRS_sequence_c = TRS_sequence.complement()
    TRS_sequence_r = TRS_sequence.reverse_complement().complement()

    if __debug__:
        logging.debug(f"TRS {TRS_sequence}")
        logging.debug(f"TRS RC: {TRS_sequence_rc}")
        logging.debug(f"TRS C: {TRS_sequence_c}")
        logging.debug(f"TRS R: {TRS_sequence_r}")

    def check_alignment(bases_clipped, score_cutoff, check_all_orientations=False):
        subgenomic_read = False
        orientation = None

        if check_trs_alignment(bases_clipped, TRS_sequence) >= score_cutoff:
            subgenomic_read = True
            orientation = 1
        elif check_all_orientations:
            if check_trs_alignment(bases_clipped, TRS_sequence_rc) >= score_cutoff:
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
    inbam_filename,
    outbam_filename,
    TRS_sequence,
    score_cutoff=50,
    n_clipped_cutoff=6,
    n_clipped_overhang=3,
    process_3prime_clipped=True,
    check_all_orientations=True,
):
    """The main entry point for the antenna algorithm"""

    check_alignment = check_alignment_factory(TRS_sequence)

    with pysam.AlignmentFile(inbam_filename, "rb", require_index=True) as inbamfile:
        with pysam.AlignmentFile(outbam_filename, "wb", template=inbamfile) as outbamfile:
            total_input_reads = inbamfile.mapped + inbamfile.unmapped

            for read in tqdm.tqdm(inbamfile, total=total_input_reads, smoothing=0):

                if (
                    read.seq is None or
                    read.is_unmapped
                    ):
                        continue

                read_length = read.template_length
                cigar = read.cigartuples

                if cigar[0][0] == BAM_CSOFT_CLIP:
                    n_clipped = cigar[0][1]
                    if n_clipped > n_clipped_cutoff:
                        five_prime_bases_clipped = read.seq[
                            0 : (n_clipped + n_clipped_overhang)
                        ]
                        subgenomic_read, orientation = check_alignment(
                            five_prime_bases_clipped,
                            score_cutoff,
                            check_all_orientations=check_all_orientations,
                        )
                        if subgenomic_read:
                            read.tags += [('TS', orientation)]

                if cigar[-1][0] == BAM_CSOFT_CLIP:
                    n_clipped = cigar[-1][1]
                    if n_clipped > n_clipped_cutoff:
                        five_prime_bases_clipped = read.seq[
                            read_length - n_clipped - n_clipped_overhang : read_length
                        ]
                        subgenomic_read, orientation = check_alignment(
                            five_prime_bases_clipped,
                            score_cutoff,
                            check_all_orientations=check_all_orientations,
                        )
                        if subgenomic_read:
                            read.tags += [('TS', orientation)]

                outbamfile.write(read)


def main():
    argparser = argparse.ArgumentParser("antenna: a tool for sgRNA read search")
    argparser.add_argument("--bam", help="input bam file", required=True)
    argparser.add_argument("--outbam", help="output bam file", required=True)
    argparser.add_argument(
        "--progress", help="display progress bar", default=False, action="store_true"
    )
    argparser.add_argument(
        "--trs-sequence",
        help="TRS sequence to search",
        default="AACCAACTTTCGATCTCTTGTAGATCTGTTCTC",
    )
    argparser.add_argument(
        "--check-all-orientations",
        default=False,
        action="store_true",
        help="Check for TRS motif in all possible orientations",
    )
    argparser.add_argument(
        "--score-cutoff", help="Score cutoff", default=50, type=int,
    )

    args = argparser.parse_args()

    if __debug__:
        logging.basicConfig(level=logging.DEBUG)

    run_antenna(
        inbam_filename=args.bam,
        outbam_filename=args.outbam,
        TRS_sequence=args.trs_sequence,
        check_all_orientations=args.check_all_orientations,
        score_cutoff=args.score_cutoff,
    )


if __name__ == "__main__":
    main()
