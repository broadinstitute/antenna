#!/usr/bin/env python3

import pysam
from Bio import pairwise2
from pybedtools import *
import numpy as np
import tqdm
import argparse

class ClassifiedRead():
    def __init__(self,sgRNA: bool,orf: str,read: pysam.AlignedRead):
        self.sgRNA = sgRNA
        self.orf = orf
        self.pos = read.pos
#        self.read = read.to_string()
        
def get_mapped_reads(bam):
    # use get_index_statistics
    mapped_reads = int(pysam.flagstat(bam).split("\n")[4].split(" ")[0])-int(pysam.flagstat(bam).split("\n")[2].split(" ")[0])-int(pysam.flagstat(bam).split("\n")[1].split(" ")[0])
    return mapped_reads


def get_coverage(start,end,inbamfile):
    coverage = []
    if start < 0:
        start=1
    for pileupcolumn in inbamfile.pileup("MN908947.3", int(start), int(end)):
        coverage.append(pileupcolumn.n)
    return np.median(coverage)


def run_antenna(orf_bed_filename, inbam_filename, TRS_sequence = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCTC'):
    orf_bed_object = BedTool(orf_bed_filename)
    reads = {}
    with pysam.AlignmentFile(inbam_filename, 'rb', require_index = True) as inbamfile:
        # Process Reads
        total_input_reads = inbamfile.mapped + inbamfile.unmapped
        for read in tqdm.tqdm(inbamfile, total = total_input_reads,smoothing=0):
            subgenomic_read = False
            if read.seq == None or read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            else:
                cigar = read.cigartuples
                if cigar[0][0] == 4:
                    # read is softclippled
                    n_clipped = cigar[0][1]
                    if n_clipped > 6:
                        perfect = 33 if n_clipped >= 33 else n_clipped * 2 - 2
                        bases_clipped = read.seq[0:n_clipped+3]
                        # Perform local alignment
                        align = pairwise2.align.localms(bases_clipped, TRS_sequence, 2, -2, -20, -.1, one_alignment_only = True)
                        align_score = align[0][2]
                        align_right_position = align[0][4]
                        if align_right_position >= len(TRS_sequence):
                            if perfect - align_score <= 0:
                                subgenomic_read = True
                        if read.query_name not in reads:
                            reads[read.query_name] = []
                        read_orf = None
                        for row in orf_bed_object:
                            if row.end >= read.reference_start >= row.start:
                                read_orf = row.name
                        if read_orf == None:
                            if subgenomic_read:
                                read_orf = "novel_" + str(read.reference_start)
                        reads[read.query_name].append(
                            ClassifiedRead(sgRNA = subgenomic_read, orf = read_orf, read=read)
                        )
        # Merge Pairs
        orfs = dict()
        for id, pair in reads.items():
            left_read = min(pair, key = lambda x:x.pos)
            if left_read.orf == None:
                continue
            if left_read.sgRNA == False:
                continue
            if left_read.orf not in orfs:
                orfs[left_read.orf] = [left_read]
            else:
                orfs[left_read.orf].append(left_read)
        # Get coverage of different orfs
        orf_coverage = {}
        for row in orf_bed_object:
            orf_coverage[row.name] = get_coverage(row.start, row.end, inbamfile)
    return {'orf_coverage': orf_coverage, 'orfs': orfs}

def save_output(count_data, output_filename):
    orfs = count_data['orfs']
    orf_coverage = count_data['orf_coverage']
    with open(output_filename,'w') as output_file:
        output_file.write(",".join(['orf','sgRNA counts','coverage','\n']))
        for orf in orfs:
            if "novel" not in orf:
                output_file.write(",".join( [orf, str(len(orfs[orf])), str(orf_coverage[orf]), '\n'] ))
        
def main():
    argparser = argparse.ArgumentParser('antenna: a tool for sgRNA read search')
    argparser.add_argument('--bam', help = 'input bam file', required = True)
    argparser.add_argument('--output-counts', help = "output csv file", required = True)
    argparser.add_argument('--orf-bed', help = "bed file with orf starting positions", required = True)
    argparser.add_argument('--progress', help = 'display progress bar')
    
    args = argparser.parse_args()
    count_data = run_antenna(orf_bed_filename = args.orf_bed, inbam_filename = args.bam)
    save_output(count_data, output_filename = args.output_counts)

    
if __name__ == "__main__":
    main()

