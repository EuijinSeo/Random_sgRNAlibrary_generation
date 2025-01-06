"""
2025-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., in revision, 2025

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: SAM file (format:text file(.txt)) (The FASTQ file consisting of sgRNA aligned to the reference genome undergoes a identical sequence compression process using FASTX-collapser, and then the SAM file obtained by aligning it to the reference genome again is used as input.)
gtf file (format:text file(.txt)) (The gtf file of the reference genome)

Output: sgRNA sequence/sgRNA reads/sgRNA-aligend gene name  (format:excel) (By providing mapped genes and reads numbers for each individual sgRNA, it is possible to create an input file used in the later 'Read mathicng' program.)
"""

# Import pandas library (Version = 2.0.3)
import pandas as pd

# Open the SAM file (replace 'Path_of_SAM_file.txt' with the actual path of your SAM file)
sam = open('Path_of_SAM_file.txt', 'r')

# Skip the first 3 lines of the SAM file (usually header information)
for _ in range(3):
    sam.readline()

# Function to process a line from the SAM file and extract relevant information
def sam_mining(x):
    xx = x.split('\t')  # Split the line by tab characters
    if xx[2] != '*':  # Check if the read is mapped
        # Extract the UMI, start position, end position, and sequence
        return [xx[0].split('-')[1], int(xx[3]), int(xx[3]) + int(xx[5][:-1]) - 1, xx[9]]
    else:
        # Return a placeholder if the read is unmapped
        return 'data empty'

# Function to process a line from the GTF file and extract relevant gene information
def gtf_mining(x):
    xx = x.split('\t')  # Split the line by tab characters
    # Extract feature type, start position, end position, and gene name
    y = [xx[2], int(xx[3]), int(xx[4]), xx[-1].split('"')[1]]
    return y

# Function to check if two intervals overlap
def in_or_not(x1, x2, y1, y2):
    if x2 < y1 or x1 > y2:  # No overlap condition
        return False
    return True  # Overlap exists

result = []  # List to store results

# Process each line in the SAM file
while True:
    sam_line = sam.readline()  # Read the next line from the SAM file
    if sam_line == '':  # End of file
        break
    else:
        sam_temp = sam_mining(sam_line)  # Extract information from the SAM line
        if sam_temp == 'data empty':  # Skip unmapped reads
            continue
        result_len = len(result)  # Track the current length of the result list
        # Open the GTF file (replace with the actual path of your GTF file)
        gtf = open('Path_of_gtf_file_of_reference_genome.txt', 'r')
        # Skip the first 3 lines of the GTF file (usually header information)
        for _ in range(3):
            gtf.readline()
        while True:
            gtf_line = gtf.readline()  # Read the next line from the GTF file
            if gtf_line == '###\n':  # End of relevant GTF section
                break
            else:
                gtf_temp = gtf_mining(gtf_line)  # Extract information from the GTF line
                if gtf_temp[0] == 'gene':  # Only consider "gene" entries
                    # Check if the SAM read overlaps with the gene region
                    if in_or_not(gtf_temp[1], gtf_temp[2], sam_temp[1], sam_temp[2]):
                        # Append the sequence, count, and gene name to the result
                        result.append([sam_temp[3], sam_temp[0], gtf_temp[-1]])
                        break
        if len(result) == result_len:  # If no overlap was found, add a placeholder
            result.append([sam_temp[3], sam_temp[0], '-'])

# Lists to store sequences, read counts, and gene names
seq_list = []
cnt_list = []
gene_list = []

# Extract data from the result for DataFrame creation
for i in result:
    seq_list.append(i[0])  # Sequence
    cnt_list.append(int(i[1]))  # Read counts
    gene_list.append(i[2])  # Gene name

# Define column names
name_list = ['Sequence', 'Read counts', 'Gene']

# Create a DataFrame from the extracted data
df = pd.DataFrame({name_list[0]: seq_list, name_list[1]: cnt_list, name_list[2]: gene_list})

# Save the DataFrame to an Excel file (specify the desired output path)
df.to_excel('reads_per_sgRNA_dCas9_220506_1m.xlsx')

# Close the GTF and SAM files
gtf.close()
sam.close()
