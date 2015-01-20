import argparse
import csv
import pysam

#---------------------
#   Script
#---------------------

parser = argparse.ArgumentParser("CRISPR.single_sgRNA_count parser.")

parser.add_argument("sam_file")
parser.add_argument("sample_name")
parser.add_argument("sgRNA_reference_file")
parser.add_argument("output_file_basename")
parser.add_argument("max_reads")

args = parser.parse_args()

sam_file = args.sam_file
sample_name = args.sample_name
sgRNA_reference_file = args.sgRNA_reference_file
output_file_basename = args.output_file_basename
max_reads = args.max_reads

# read sgRNA_reference file (a csv file containing reference sgRNA sequences and names)
UNMAPPED = "UNMAPPED"
sgRNA_name_list = []
sgRNA_mapping_counts = {}
with open(sgRNA_reference_file, "rU") as input_handle:
    reader = csv.reader(input_handle)
    for row in reader:
        sgRNA_name_list.append(row[1])
    sgRNA_name_list.append(UNMAPPED)
    
sgRNA_mapping_counts = {}
for sgRNA_name in sgRNA_name_list:
    sgRNA_mapping_counts[sgRNA_name] = 0

# count aligned sgRNA reads

py_sam_file = pysam.Samfile(sam_file, "r")

sam_iterator = py_sam_file.fetch()

read_count = 0
unexpected_flags = 0

for alignedRead in sam_iterator:
    if read_count == max_reads:
        break
    read_count += 1

    if alignedRead.flag == 0:
        sgRNA_name = py_sam_file.getrname(alignedRead.rname)
    elif alignedRead.flag == 4:
        sgRNA_name = UNMAPPED
    else:
        print("unexpected read flag: %s, %i, %s"% (alignedRead.qname, alignedRead.flag, py_sam_file.getrname(alignedRead.rname)))
        sgRNA_name = UNMAPPED
        unexpected_flags += 1
               
    sgRNA_mapping_counts[sgRNA_name] += 1 

py_sam_file.close()

print("unexpected flag count = %i" % unexpected_flags)

# write csv file with count results
handle = open(output_file_basename + ".csv", "wt") 
writer = csv.writer(handle)

writer.writerow(("NAME", sample_name))

for sgRNA_name in sgRNA_name_list:
    writer.writerow((sgRNA_name, sgRNA_mapping_counts[sgRNA_name]))
        
handle.close()

    
