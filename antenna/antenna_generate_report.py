import antenna.antenna_count_reads
import argparse


def main():
    argparser = argparse.ArgumentParser("antenna_generate_report")
    argparser.add_argument("--bam", help="input bam file", required=True)
    argparser.add_argument("--output-image", help="output image", required=True)
    argparser.add_argument("--verbose", help="verbose output", default=True, action="store_true")
    args = argparser.parse_args()
    
    if args.verbose:
        print("Loading read sgRNA scores...")
    sgRNA_scores = antenna.antenna_count_reads.load_sgRNA_scores(args.bam)
        
    # Output image report
    antenna.antenna_count_reads.plot_read_score_dist(sgRNA_scores, filename=args.output_image)
    

if __name__ == "__main__":
    main()