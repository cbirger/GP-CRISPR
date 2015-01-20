import argparse
import csv
import pysam

#---------------------
#   Script
#---------------------

parser = argparse.ArgumentParser("CRISPR.dual_sgRNA_count parser.")

parser.add_argument("fwd_reads_sam_file")
parser.add_argument("rvs_reads_sam_file")
parser.add_argument("sample_name")
parser.add_argument("sgRNA_reference_file")
parser.add_argument("output_file_basename")
parser.add_argument("max_reads")

args = parser.parse_args()

fwd_reads_sam_file = args.fwd_reads_sam_file
rvs_reads_sam_file = args.rvs_reads_sam_file
sample_name = args.sample_name
sgRNA_reference_file = args.sgRNA_reference_file
output_file_basename = args.output_file_basename
max_reads = args.max_reads

# read sgRNA_reference file (a csv file containing reference sgRNA sequences and names)
UNMAPPED = "UNMAPPED"

sgRNA_name_list = []
with open(sgRNA_reference_file, "rU") as input_handle:
    reader = csv.reader(input_handle)
    for row in reader:
        sgRNA_name_list.append(row[1])
    sgRNA_name_list.append(UNMAPPED)
    
sgRNA_mapping_counts = {}
for sgRNA_name_fwd in sgRNA_name_list:
    sgRNA_mapping_counts[sgRNA_name_fwd] = {}
    for sgRNA_name_rvs in sgRNA_name_list:
        sgRNA_mapping_counts[sgRNA_name_fwd][sgRNA_name_rvs] = 0

# count aligned sgRNA reads

fwd_samfile = pysam.Samfile(fwd_reads_sam_file, "r")
rvs_samfile = pysam.Samfile(rvs_reads_sam_file, "r")

fwd_iterator = fwd_samfile.fetch()
rvs_iterator = rvs_samfile.fetch()

read_count = 0

fwd_reads_unexpected_flags = 0
rvs_reads_unexpected_flags = 0

for fwd_alignedRead in fwd_iterator:
    if read_count == max_reads:
        break
    read_count += 1
    rvs_alignedRead = rvs_iterator.next()

    # bowtie1 appeared to in some cases extend the qname with additonal field (e.g.,
    # "1:N:0:1") that is separated from the main qname with a single space. We don't
    # want to include this additional field in the following assertion.
    assert (fwd_alignedRead.qname.split()[0] == rvs_alignedRead.qname.split()[0]), "Unmated reads: read_count = %i, fwd_qname = %s, rvs_qname = %s" % (read_count, fwd_alignedRead.qname, rvs_alignedRead.qname)

    if fwd_alignedRead.flag == 0:
        fwd_name = fwd_samfile.getrname(fwd_alignedRead.rname)
    elif fwd_alignedRead.flag == 4:
        fwd_name = UNMAPPED
    else:
        print("unexpected fwd read flag: %s, %i, %s"% (fwd_alignedRead.qname, fwd_alignedRead.flag, fwd_samfile.getrname(fwd_alignedRead.rname)))
        fwd_name = UNMAPPED
        fwd_reads_unexpected_flags += 1
               
    if rvs_alignedRead.flag == 16:
        rvs_name = rvs_samfile.getrname(rvs_alignedRead.rname)
    elif rvs_alignedRead.flag == 4:
        rvs_name = UNMAPPED
    else:
        print("unexpected rvs flag: %s, %i, %s"% (rvs_alignedRead.qname, rvs_alignedRead.flag, rvs_samfile.getrname(rvs_alignedRead.rname)))
        fwd_name = UNMAPPED
        rvs_reads_unexpected_flags +=1
        
    sgRNA_mapping_counts[fwd_name][rvs_name] += 1 

fwd_samfile.close()
rvs_samfile.close()

print("fwd_reads_unexpected_flags = %i" % fwd_reads_unexpected_flags)
print("rvs_reads_unexpected_flags = %i" % rvs_reads_unexpected_flags)

# write csv file with quantitation results
handle = open(output_file_basename + ".csv", "wt") 
writer = csv.writer(handle)

writer.writerow(("NAME_1", "NAME_2", sample_name))

for first_sgRNA in sgRNA_name_list:
    for second_sgRNA in sgRNA_name_list:
        writer.writerow((first_sgRNA, second_sgRNA, sgRNA_mapping_counts[first_sgRNA][second_sgRNA]))
        
handle.close()

    
