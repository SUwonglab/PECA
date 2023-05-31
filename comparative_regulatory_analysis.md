# Comparative gene regulatory analysis via PECA
If you have three or more samples in each group, you can leverage PECA for performing comparative regulatory analysis. PECA is a powerful tool designed for analyzing bulk or pseudobulk RNA-seq and ATAC-seq data. It generates condition-specific transcriptional regulatory networks and identifies condition-specific driver regulators. Additionally, PECA allows for controlling covariates such as race, age, sex, and treatment to obtain more accurate results
## Simple run
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt
```
Here 
## Step by step
