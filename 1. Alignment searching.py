"""

2024-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., Nat. Chem. Biol., 2024

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: SAM file (format:text file(.txt)) (The reference genome alignment result of your FASTQ obtained by NGS) GCF file (format:text file(.txt)) (The GCF file of reference genome)

Output: A, B, A/(A+B), C
A: The number of reads that target the non-template strand in the SAM file
B: The number of reads that target the template strand in the SAM file
C: The number of reads in the SAM file
"""

def sam_mining(x):
    xx = x.split('\t')
    y = [xx[1],int(xx[3]),int(xx[3])+int(xx[5][:-1])-1]
    return y

def gtf_mining(x):
    xx = x.split('\t')
    y = [xx[2],int(xx[3]),int(xx[4]),xx[-3]]
    return y

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

sam = open('Path_of_alignment_1mismatch_aligned_only_sam.txt', 'r') # Replace 'Path_of_alignment_1mismatch_aligned_only_sam.txt' with your SAM file path

for _ in range(3):
    sam.readline()
for _ in sam:
    tot += 1

sam = open('Path_of_alignment_1mismatch_aligned_only_sam.txt', 'r') # Replace 'Path_of_alignment_1mismatch_aligned_only_sam.txt' with your SAM file path

for _ in range(3):
    sam.readline()
    
cnt = 0

A = 0 # (The number of reads that target the non-template strand in the SAM file)
B = 0 # (The number of reads that target the template strand in the SAM file)
C = 0 # (The number of reads in the SAM file)

while True:
    cnt += 1
    sam_temp = sam.readline()
    if sam_temp == '':
        break
    sam_seq = sam_mining(sam_temp)
    gtf = open('Path_of_GCF_file.txt', 'r') # Replace 'Path_of_GCF_file.txt' with your GCF file path
    for _ in range(3):
        gtf.readline()
    while True:
        get_temp = gtf.readline()
        if get_temp == '###\n':
            break
        else:
            gtf_seq = gtf_mining(get_temp)
            if gtf_seq[0] != 'gene':
                continue
            else:
                if in_or_not(gtf_seq,sam_seq):
                    if gtf_seq[3] == '+' and sam_seq[0] == '16':
                        A += 1
                    elif gtf_seq[3] == '+' and sam_seq[0] == '0':
                        B += 1
                    elif gtf_seq[3] == '-' and sam_seq[0] == '0':
                        A += 1
                    else:
                        B += 1
                else:
                    C+=1 

print(A,B,A/(A+B),C)

sam.close()
gtf.close()
