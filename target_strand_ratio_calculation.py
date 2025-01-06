"""
2025-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., in revision, 2025

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: SAM file (format:text file(.txt)) (The reference genome alignment result of your FASTQ obtained by NGS) GTF file (format:text file(.txt)) (The GTF file of reference genome)

Output: A, B, A/(A+B), C
A: The number of reads that target the non-template strand in the SAM file
B: The number of reads that target the template strand in the SAM file
C: The number of reads in the SAM file that are not matched in the reference gene
"""

# Function to extract read information from a SAM file
def sam_mining(x):
    xx = x.split('\t')  # Split the SAM line into fields using tab as the delimiter
    y = [xx[1], int(xx[3]), int(xx[3]) + int(xx[5][:-1]) - 1]  # Extract mapping direction, start, and end positions
    return y  # y[0]: mapping direction (0 = forward, 16 = reverse), y[1]: start position, y[2]: end position

# Function to extract information from a GTF file
def gtf_mining(x):
    xx = x.split('\t')  # Split the GTF line into fields using tab as the delimiter
    y = [xx[2], int(xx[3]), int(xx[4]), xx[-3]]  # Extract feature type, start and end positions, and genome direction
    return y  # y[0]: sequence type (e.g., gene, CDS), y[1]: start position, y[2]: end position, y[3]: genome direction (+ or -)

# Function to check if two DNA regions overlap
def in_or_not(x, y):
    # Conditions for overlap between regions x and y
    if x[1] >= y[1] and x[1] <= y[2]:  # x start within y region
        return True
    elif x[2] >= y[1] and x[2] <= y[2]:  # x end within y region
        return True
    elif y[1] >= x[1] and y[2] <= x[2]:  # y region completely within x region
        return True
    elif y[2] >= x[1] and y[2] <= x[2]:  # y end within x region
        return True
    else:  # No overlap
        return False

# Open the SAM file
sam = open('Path_of_alignment_1mismatch_aligned_only_sam.txt', 'r')  # Replace with your SAM file path

# Skip the first three lines of the SAM file (header information)
for _ in range(3):
    sam.readline()

# Initialize counters
A = 0  # Number of reads targeting the non-template strand
B = 0  # Number of reads targeting the template strand
C = 0  # Number of reads not mapped to any reference gene

# Process each read in the SAM file
while True:
    sam_temp = sam.readline()  # Read a line from the SAM file

    # Exit the loop if the end of the SAM file is reached
    if sam_temp == '':
        break

    sam_seq = sam_mining(sam_temp)  # Extract relevant information from the SAM line

    # Open the GTF file
    gtf = open('Path_of_gtf_file.txt', 'r')  # Replace with your GTF file path

    # Skip the first three lines of the GTF file (header information)
    for _ in range(3):
        gtf.readline()

    # Process each entry in the GTF file
    while True:
        get_temp = gtf.readline()  # Read a line from the GTF file

        # Exit the loop if the end of the GTF file is reached
        if get_temp == '###\n':
            break

        gtf_seq = gtf_mining(get_temp)  # Extract relevant information from the GTF line

        # Skip non-gene entries
        if gtf_seq[0] != 'gene':
            continue

        # Check if the SAM read overlaps with the gene region
        if in_or_not(gtf_seq, sam_seq):
            # Classify the read based on genome direction and read direction
            if gtf_seq[3] == '+' and sam_seq[0] == '16':  # Forward genome, reverse read
                A += 1
            elif gtf_seq[3] == '+' and sam_seq[0] == '0':  # Forward genome, forward read
                B += 1
            elif gtf_seq[3] == '-' and sam_seq[0] == '0':  # Reverse genome, forward read
                A += 1
            elif gtf_seq[3] == '-' and sam_seq[0] == '16':  # Reverse genome, reverse read
                B += 1
        else:
            # If no overlap is found, count it as unmatched
            C += 1

# Print the results
print(f'The number of reads that target the non-template strand in the SAM file (A): {A}')
print(f'The number of reads that target the template strand in the SAM file (B): {B}')
print(f'A/(A+B): {A/(A+B):.2f}')  # Proportion of reads targeting the non-template strand
print(f'The number of reads in the SAM file that are not matched in the reference gene: {C}')

# Close the SAM and GTF files
sam.close()
gtf.close()
