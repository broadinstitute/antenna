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
    # TODO: Review these penatlies
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

    def check_alignment(bases_clipped):
        o_score = check_trs_alignment(bases_clipped, TRS_sequence)
        rc_score = check_trs_alignment(bases_clipped, TRS_sequence_rc)
        r_score = check_trs_alignment(bases_clipped, TRS_sequence_r)
        c_score = check_trs_alignment(bases_clipped, TRS_sequence_c)

        return (o_score, rc_score, r_score, c_score)

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
    n_clipped_cutoff=6,
    n_clipped_overhang=3,
    process_3prime_clipped=True,
):
    """The main entry point for the antenna algorithm"""

    check_alignment = check_alignment_factory(TRS_sequence)

    with pysam.AlignmentFile(inbam_filename, "rb", require_index=True) as inbamfile:
        with pysam.AlignmentFile(
            outbam_filename, "wb", template=inbamfile
        ) as outbamfile:
            total_input_reads = inbamfile.mapped + inbamfile.unmapped

            for read in tqdm.tqdm(inbamfile, total=total_input_reads, smoothing=0):
                trs_found = False
                orientation_flag = 0x0

                if read.seq is None or read.is_unmapped:
                    outbamfile.write(read)
                    continue

                read_length = read.template_length
                cigar = read.cigartuples
                
                p5_o_score = 0
                p5_rc_score = 0
                p5_r_score = 0
                p5_c_score = 0
                p3_o_score = 0
                p3_rc_score = 0
                p3_r_score = 0
                p3_c_score = 0

                if cigar[0][0] == BAM_CSOFT_CLIP:
                    n_clipped = cigar[0][1]
                    if n_clipped > n_clipped_cutoff:
                        five_prime_bases_clipped = read.seq[
                            0 : (n_clipped + n_clipped_overhang)
                        ]
                        (
                            p5_o_score,
                            p5_rc_score,
                            p5_r_score,
                            p5_c_score,
                        ) = check_alignment(five_prime_bases_clipped)

                if cigar[-1][0] == BAM_CSOFT_CLIP:
                    n_clipped = cigar[-1][1]
                    if n_clipped > n_clipped_cutoff:
                        three_prime_bases_clipped = read.seq[
                            read_length - n_clipped - n_clipped_overhang : read_length
                        ]
                        (
                            p3_o_score,
                            p3_rc_score,
                            p3_r_score,
                            p3_c_score,
                        ) = check_alignment(three_prime_bases_clipped)


                TO_string = f'{p5_o_score:.0f},{p5_rc_score:.0f},{p5_r_score:.0f},{p5_c_score:.0f},{p3_o_score:.0f},{p3_rc_score:.0f},{p3_r_score:.0f},{p3_c_score:.0f}'

                read.tags += [("TO", TO_string)]

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

    args = argparser.parse_args()

    if __debug__:
        logging.basicConfig(level=logging.DEBUG)

    run_antenna(
        inbam_filename=args.bam,
        outbam_filename=args.outbam,
        TRS_sequence=args.trs_sequence,
    )


if __name__ == "__main__":
    main()
