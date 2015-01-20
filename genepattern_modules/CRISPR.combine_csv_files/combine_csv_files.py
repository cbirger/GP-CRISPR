import argparse
import csv

parser = argparse.ArgumentParser("CRISPR.merge_csv_count_files")

parser.add_argument("csv_file_list")
parser.add_argument("output_file_basename")

args = parser.parse_args()

file_list_filename = args.csv_file_list
output_file_basename = args.output_file_basename

row_names = []
combined_counts = {}
combined_header_row = None
paired_sgRNA_screen = False

# read first file in list

with open(file_list_filename, "rU") as f:
    file_list = f.read().splitlines()
    
with open(file_list[0], "rU") as input_handle:
    reader = csv.reader(input_handle)
    combined_header_row = reader.next()
    if len(combined_header_row) == 3:
        paired_sgRNA_screen = True
        
    if not paired_sgRNA_screen:
        for row in reader:
            row_names.append(row[0:1])
            combined_counts[row[0]] = [row[1]]
    else:
        # is paried sgRNA screen
        for row in reader:
            row_names.append(row[0:2])
            combined_counts[row[0]+'#'+row[1]] = [row[2]]


# now read remaining files

for csv_file in file_list[1:]:
    with open(csv_file, "rU") as input_handle:
        reader = csv.reader(input_handle)

        header_row = reader.next()
        if not paired_sgRNA_screen:
            assert(header_row[0] == combined_header_row[0])
            combined_header_row.append(header_row[1])
            for row in reader:
                combined_counts[row[0]].append(row[1])
        else:
            # is a paired sgRNA screen
            assert(header_row[0] == combined_header_row[0] and 
                   header_row[1] == combined_header_row[1])
            combined_header_row.append(header_row[2])
            for row in reader:
                combined_counts[row[0]+'#'+row[1]].append(row[2])
                
# now write to single csv file

handle = open(output_file_basename + ".csv", "wt")
writer = csv.writer(handle)

writer.writerow(tuple(combined_header_row))

if not paired_sgRNA_screen:
    for name in row_names:
        # note that combined_counts contains lists of strings, so we are concatenating lists of strings
        row_data = name + combined_counts[name[0]]
        writer.writerow(tuple(row_data))
else:
    # is a paired sgRNA screen
    for name_pair in row_names:
        # note that combined_counts contains lists of strings, so we are concatenating lists of strings
        row_data = name_pair + combined_counts[name_pair[0]+"#"+name_pair[1]]
        writer.writerow(tuple(row_data))

handle.close()
