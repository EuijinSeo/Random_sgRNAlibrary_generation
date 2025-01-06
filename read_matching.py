"""

2025-x-xx

Author : Euijin Seo (ufoooo1391@postech.ac.kr)

Paper : Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application, Lee et al., in revision, 2025

Copyright by POSTECH Synthetic Biology Laboratory(SBL)

Inputs: Modified version of the output of 'Calculating_read_per_sgRNA' (format:excel). 'Noramlized read counts' should be calculated using the below formular and inserted as a new column between 'Read counts' and 'Gene'
Normalized read count: (Read counts/Total read counts)*10^6+1
Output: sgRNA sequence/Normalized sgRNA reads in control condition/Normalized sgRNA reads in selection condition/sgRNA-aligend gene name  (format:excel)

This program matches the normalized read counts of each sgRNA in control, selection condition, enabling the calculation of the fold change of the abundance of each sgRNA.


"""

# Import pandas library (Version = 2.0.3)
import pandas as pd

# Load the selected sgRNA data from an Excel file
sel_data = pd.read_excel('Path_to_selected_sgRNA_file.xlsx')  # Replace with the path to your selected sgRNA file

# Load the control sgRNA data from an Excel file
ctrl_data = pd.read_excel('Path_to_control_sgRNA_file.xlsx')  # Replace with the path to your control sgRNA file

# Extract relevant columns from the selected sgRNA data
seqsel = sel_data['Sequence'].tolist()  # List of sgRNA sequences in the selection condition
nsel = sel_data['Normalized read counts'].tolist()  # List of normalized read counts in the selection condition
gensel = sel_data['gene'].tolist()  # List of genes associated with sgRNAs in the selection condition

# Extract relevant columns from the control sgRNA data
seqctrl = ctrl_data['Sequence'].tolist()  # List of sgRNA sequences in the control condition
nctrl = ctrl_data['Normalized read counts'].tolist()  # List of normalized read counts in the control condition
genctrl = ctrl_data['gene'].tolist()  # List of genes associated with sgRNAs in the control condition

# Initialize lists to store results
seqlist = []    # Combined list of sequences
nsellist = []   # Combined list of normalized reads in the selection condition
nctrllist = []  # Combined list of normalized reads in the control condition
genlist = []    # Combined list of gene names

# Compare sgRNA sequences between selected and control conditions
for i in range(len(seqsel)):
    original_seq = len(seqctrl)  # Track the original number of control sequences
    for j in range(len(seqctrl)):
        # If the sequence from the selection condition matches a sequence in the control condition
        if seqsel[i] == seqctrl[j]:
            # Append matching data to the result lists
            seqlist.append(seqsel[i])
            nsellist.append(nsel[i])
            nctrllist.append(nctrl[j])
            genlist.append(gensel[i])
            
            # Remove the matched sequence and corresponding data from the control lists
            seqctrl.remove(seqctrl[j])
            nctrl.remove(nctrl[j])
            genctrl.remove(genctrl[j])
            break

    # If no match was found, append the selection sequence with a default control value
    if len(seqctrl) == original_seq:
        seqlist.append(seqsel[i])
        nsellist.append(nsel[i])
        nctrllist.append(1)  # Default value for unmatched sequences in the control condition
        genlist.append(gensel[i])

# Reset the lists for another pass to account for unmatched control sequences
seqsel = sel_data['Sequence'].tolist()
seqctrl = ctrl_data['Sequence'].tolist()
nsel = sel_data['Normalized read counts'].tolist()
nctrl = ctrl_data['Normalized read counts'].tolist()
gensel = sel_data['gene'].tolist()
genctrl = ctrl_data['gene'].tolist()

# Check for unmatched control sequences
for i in range(len(seqctrl)):
    original_seq = len(seqsel)  # Track the original number of selection sequences
    tf = 0  # Flag to track whether a match is found
    for j in range(len(seqsel)):
        if seqctrl[i] == seqsel[j]:
            tf = 1  # Match found
            break
    if tf == 0:  # If no match is found
        seqlist.append(seqctrl[i])
        nsellist.append(1)  # Default value for unmatched sequences in the selection condition
        nctrllist.append(nctrl[i])
        genlist.append(genctrl[i])

# Create a DataFrame to store the combined results
result = pd.DataFrame({
    'Sequence': seqlist,
    'Normalized reads in control condition': nctrllist,
    'Normalized reads in selection condition': nsellist,
    'Gene': genlist
})

# Save the results to an Excel file
result.to_excel('output.xlsx')  # Replace with your desired output file path
