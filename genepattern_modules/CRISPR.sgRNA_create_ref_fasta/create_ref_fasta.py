import argparse
import csv
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#---------------------
#   Script
#---------------------

parser = argparse.ArgumentParser("CRISPR.sgRNA_create_ref_fasta parser.")

parser.add_argument("reference")
parser.add_argument("basename")

args = parser.parse_args()

reference_file = args.reference
basename = args.basename

sgRNA_dict = {}
sgRNA_length = None
sgRNA_seq_list = []

with open(reference_file, "rU") as input_handle:
    reader = csv.reader(input_handle)

    for row in reader:
        if row[0] in sgRNA_dict:
            raise Error("Duplicate sgRNA: (%s, %s)" % (row[0], row[1]))

        if sgRNA_length is None:
            sgRNA_length = len(row[0])
        else:
            assert (len(row[0]) == sgRNA_length), "sgRNA sequence length differs.  Record ({}, {})".format(row[0], row[1])
        sgRNA_seq_list.append(SeqRecord(Seq(row[0], DNAAlphabet()),id=row[1],description=""))

fasta_file = basename + ".fasta"
with open(fasta_file, "w") as output_handle:
    SeqIO.write(sgRNA_seq_list, output_handle, "fasta")

            
