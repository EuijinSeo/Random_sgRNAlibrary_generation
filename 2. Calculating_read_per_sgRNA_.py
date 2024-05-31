"""

2024-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., Nat. Chem. Biol., 2024

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: SAM file (format:text file(.txt)) (The FASTQ file consisting of sgRNA aligned to the reference genome undergoes a identical sequence compression process using FASTX-collapser, and then the SAM file obtained by aligning it to the reference genome again is used as input.)
GCF file (format:text file(.txt)) (The GCF file of the reference genome)

Output: sgRNA sequence/sgRNA reads/sgRNA-aligend gene name  (format:excel) (By providing mapped genes and reads numbers for each individual sgRNA, it is possible to create an input file used in the later 'Read mathicng' program.)

"""

import pandas as pd

sam = open('Path_of_SAM_file.txt','r') # Replace 'Path_of_SAM_file.txt' with your SAM file path

for _ in range(3):
    sam.readline()

def sam_mining(x):
  xx = x.split('\t')
  if xx[2] != '*':
      return [xx[0].split('-')[1],int(xx[3]),int(xx[3])+int(xx[5][:-1])-1,xx[9]]
  else:
      return 'data empty'

def gca_mining(x):
  xx = x.split('\t')
  y = [xx[2],int(xx[3]),int(xx[4]),xx[-1].split('"')[1]]
  return y

def in_or_not(x1,x2,y1,y2):
  if x2 < y1 or x1 > y2:
    return False
  return True

result = []

while True:
    sam_line = sam.readline()
    if sam_line=='':
        break
    else:
        sam_temp = sam_mining(sam_line)
        if sam_temp == 'data empty':
            continue
        result_len = len(result)
        gca = open('Path_of_GCF_file_of_reference_genome.txt','r') # Replace 'Path_of_GCF_file_of_reference_genome.txt' with your GCF file path
        for _ in range(3):
            gca.readline()
        while True:
            gca_line = gca.readline()
            if gca_line =='###\n':
                break
            else:
                gca_temp = gca_mining(gca_line)
                if gca_temp[0] == 'gene':
                    if in_or_not(gca_temp[1],gca_temp[2],sam_temp[1],sam_temp[2]):
                        result.append([sam_temp[3],sam_temp[0],gca_temp[-1]])
                        break
        if len(result) == result_len:
            result.append([sam_temp[3],sam_temp[0],'-'])
            

seq_list = []
cnt_list = []
gene_list = []

for i in result:
    seq_list.append(i[0])
    cnt_list.append(int(i[1]))
    gene_list.append(i[2])


name_list = ['Sequence','Read counts','gene']

df = pd.DataFrame({name_list[0]:seq_list, name_list[1]:cnt_list, name_list[2]:gene_list})
df.to_excel('reads_per_sgRNA_dCas9_220506_1m.xlsx') # Put your desired output file directory

gca.close()
sam.close()
