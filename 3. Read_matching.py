"""

2024-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., Nat. Chem. Biol., 2024

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: Modified version of the output of 'Calculating_read_per_sgRNA' (format:excel). 'Noramlized read counts' should be calculated using the below formular and inserted as a new column between 'Read counts' and 'Gene'
Normalized read count: (Read counts/Total read counts)*10^6+1
Output: sgRNA sequence/Normalized sgRNA reads in control condition/Normalized sgRNA reads in selection condition/sgRNA-aligend gene name  (format:excel)

This program matches the normalized read counts of each sgRNA in control, selection condition, enabling the calculation of the fold change of the abundance of each sgRNA.


"""

import pandas as pd

sel_data = pd.read_excel('Path_to_selected_sgRNA_file.xlsx') # Replace 'Path_to_selected_sgRNA_file.xlsx' with your selected sgRNA file path
ctrl_data = pd.read_excel('Path_to_control_sgRNA_file.xlsx') # Replace 'Path_to_control_sgRNA_file.xlsx' with your control sgRNA file path

seqsel = sel_data['Sequence'].tolist()
seqctrl = ctrl_data['Sequence'].tolist()
nsel = sel_data['Normalized read counts'].tolist()
nctrl = ctrl_data['Normalized read counts'].tolist()
gensel = sel_data['gene'].tolist()
genctrl = ctrl_data['gene'].tolist()

seqlist = []
nsellist = []
nctrllist = []
genlist = []

for i in range(len(seqsel)):
    original_seq = len(seqctrl)
    for j in range(len(seqctrl)):
        if seqsel[i] == seqctrl[j]:
            seqlist.append(seqsel[i])
            nsellist.append(nsel[i])
            nctrllist.append(nctrl[j])
            genlist.append(gensel[i])
            seqctrl.remove(seqctrl[j])
            nctrl.remove(nctrl[j])
            genctrl.remove(genctrl[j])
            break

    if len(seqctrl) == original_seq:
        seqlist.append(seqsel[i])
        nsellist.append(nsel[i])
        nctrllist.append(1)
        genlist.append(gensel[i])

seqsel = sel_data['Sequence'].tolist()
seqctrl = ctrl_data['Sequence'].tolist()
nsel = sel_data['Normalized read counts'].tolist()
nctrl = ctrl_data['Normalized read counts'].tolist()
gensel = sel_data['gene'].tolist()
genctrl = ctrl_data['gene'].tolist()

for i in range(len(seqctrl)):
    original_seq = len(seqsel)
    tf = 0
    for j in range(len(seqsel)):
        if seqctrl[i] == seqsel[j]:
            tf = 1
            break
    if tf == 0:
        seqlist.append(seqctrl[i])
        nsellist.append(1)
        nctrllist.append(nctrl[i])
        genlist.append(genctrl[i])

result = pd.DataFrame({'Sequence':seqlist, 'Normalized reads in control condition':nctrllist, 'Normalized reads in selection condition':nsellist,'Gene':genlist})

result.to_excel('output.xlsx') # Put your desired output file directory
