"""
2025-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., in revision, 2025

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: SAM file (format:text file(.txt)) (The reference genome alignment result of your FASTQ obtained by NGS) GCF file (format:text file(.txt)) (The GCF file of reference genome)

Output: A, B, A/(A+B), C
A: The number of reads that target the non-template strand in the SAM file
B: The number of reads that target the template strand in the SAM file
C: The number of reads in the SAM file
"""

# Obtain the read information from SAM file 
def sam_mining(x):
    xx = x.split('\t')
    y = [xx[1],int(xx[3]),int(xx[3])+int(xx[5][:-1])-1]
    return y # y[0]:, y[1]: start position, y[2]: end position

# Obtain the read information from GCF file
def gcf_mining(x):
    xx = x.split('\t')
    y = [xx[2],int(xx[3]),int(xx[4]),xx[-3]]
    return y # y[0]:, y[1]: start position, y[2]:end position

# Check if there is an overlapping region between the coordinates of two DNA sequences (SAM, GCF)
def in_or_not(x,y):
    if x[1] >= y[1] and x[1] <=y[2]:
        return True
    elif x[2] >= y[1] and x[2] <=y[2]:
        return True
    elif y[1] >= x[1] and y[2] <=x[2]:
        return True
    elif y[2] >= x[1] and y[2] <=x[2]:
        return True
    else:
        return False

tot  = 0

# Open and read the SAM file
sam = open('Path_of_alignment_1mismatch_aligned_only_sam.txt', 'r') # Replace 'Path_of_alignment_1mismatch_aligned_only_sam.txt' with your SAM file path

# Remove the first three lines of unnecessary information from a SAM file before running the program 
for _ in range(3):
    sam.readline()
    
#
A = 0 # (The number of reads that target the non-template strand in the SAM file)
B = 0 # (The number of reads that target the template strand in the SAM file)
C = 0 # (The number of reads in the SAM file)

# Perform the following operations on every line of the SAM file
while True:
    sam_temp = sam.readline()

    # Terminate the while loop when the SAM file has been read completely
    if sam_temp == '':
        break
        
    sam_seq = sam_mining(sam_temp)

    # Open and read the GCF file
    gcf = open('Path_of_GCF_file.txt', 'r') # Replace 'Path_of_GCF_file.txt' with your GCF file path
    
    # Remove the first three lines of unnecessary information from a GCF file before running the program  
    for _ in range(3):
        gcf.readline()
        
    # Perform the following operations on every line of the GCF file
    while True:
        get_temp = gcf.readline()

        # Terminate the while loop when the GCF file has been read completely
        if get_temp == '###\n':
            break
        else:
            gcf_seq = gcf_mining(get_temp)
            if gcf_seq[0] != 'gene':
                continue
            else:
                if in_or_not(gcf_seq,sam_seq):
                    if gcf_seq[3] == '+' and sam_seq[0] == '16':
                        A += 1
                    elif gcf_seq[3] == '+' and sam_seq[0] == '0':
                        B += 1
                    elif gcf_seq[3] == '-' and sam_seq[0] == '0':
                        A += 1
                    else:
                        B += 1
                else:
                    C+=1 

print(A,B,A/(A+B),C)

# Close the opened SAM file and GCF file after reading them
sam.close()
gcf.close()
