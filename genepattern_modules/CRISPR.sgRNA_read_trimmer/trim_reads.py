import argparse
import sgrna_trimmer

#---------------------
#   Script
#---------------------


parser = argparse.ArgumentParser("CRISPR.dual_sgRNA_read_trimmer parser.")

parser.add_argument('reads')
parser.add_argument('prefix')
parser.add_argument('mismatches', type=int)
parser.add_argument('sgRNA_length', type=int)
parser.add_argument('max_num_reads', type=int)
parser.add_argument('basename')
parser.add_argument('-r', action='store_true', dest='reverse_reads')

args = parser.parse_args()

reads_file = args.reads
reverse_reads = args.reverse_reads
prefix_sequence = args.prefix
mismatches = args.mismatches
sgRNA_length = args.sgRNA_length
max_num_reads = args.max_num_reads
output_basename = args.basename

trimmerObj = sgrna_trimmer.SgRNAExtractor(reads_file, reverse_reads, prefix_sequence,
                                          mismatches, sgRNA_length, max_num_reads,
                                          output_basename)
trimmerObj.extract_sgRNAs()
trimmerObj.write_summary_stats()

