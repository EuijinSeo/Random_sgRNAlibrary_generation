# Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application
In our study, 'Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application', we conducted high-throughput phenotypic screenings to find essential genes and auxotrophic genes of Escherichia coli by CRISPRi using enzymatically generated sgRNA. Here, we present the code used in this study for analysis of NGS data; matching the reads from control and selective conditions, calculating the reads per sgRNAs and characterization of the sgRNA library. 

![Figure 1_2](https://github.com/EuijinSeo/Random_sgRNAlibrary_generation/assets/97028331/e9f3c8a1-931f-4321-b6a0-31a2f9bd6f50)

---
## Table of contents
##### Target Strand Ratio Calculation (target_strand_ratio_caculation.py)
  - This custom script calculates the target strand ratio in a sgRNA library by determining whether each sgRNA targets the non-template strand or the template strand. The script utilizes sequencing results (SAM file) and gene annotation data (GTF file) as input.
##### Calculating Reads per sgRNA (calculating_read_per_sgRNA.py)
  - This custom script calculates read counts and mapped genes for each sgRNA in a CRISPR library.It takes a SAM file, generated by aligning a compressed FASTA file of sgRNAs (produced using FASTX-collapser) to the reference genome, as input and uses the genome's GTF annotation to identify mapped genes. The output includes sgRNA sequences, read counts, and associated gene names, providing a concise summary for downstream analysis.
##### Read matching (read_matching.py)
  - This custom script matches identical sgRNAs between selective and control conditions, providing their normalized read counts and mapped genes. It takes as input Excel files containing sgRNAs with normalized read counts for each condition (derived from the output of "Calculating Reads per sgRNA" with added normalization). The output includes normalized read counts for each sgRNA in both conditions along with their mapped genes.

---
## Prerequisites
- Python = 3.6.0
- Pandas = 2.0.3

---
## Citation
Lee, J., Jeon, H.H., Seo, E., Park, S., Choe, D., Cho, B.-K. and Lee, J.W.\* “Generation of an Ultra-dense Genome-scale Single-guide RNA Library for CRISPRi Application”, (submitted).

---
## E-mail
- 1st author : Jiseon Lee (ljs72@postech.ac.kr)\n
- Corresponding author : Jeong Wook Lee (jeongwook@postech.ac.kr)
- Program author : Euijin Seo (ufoooo1391@postech.ac.kr)
