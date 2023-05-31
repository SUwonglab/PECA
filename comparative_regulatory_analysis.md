# Comparative gene regulatory analysis via PECA
If you have three or more samples in each group, you can leverage PECA for performing comparative regulatory analysis. PECA is a powerful tool designed for analyzing bulk or pseudo-bulk RNA-seq and ATAC-seq data. It generates condition-specific transcriptional regulatory networks and identifies condition-specific driver regulators. Additionally, PECA allows for controlling covariates such as race, age, sex, and treatment to obtain more accurate results
## Simple run
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt
```
Here Group1.txt Group2.txt All_sample_name.txt are three one-column-text files containing sample names, which have to be consistent with the sample name given under the ./Input folder. For the input detail please see README.md.
## Additional Features and Customization
### Covariate control
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt --Design_Matrix covariates_knock.txt
```
Write all covariates into a design matrix (i.e. covariates_knock.txt). Here the design matrix (all elements have to be **numbers**) is a **tab-delimited** text file as this example:
```
Name  Is_treatment  Is_male  Age
Sample1  1  1  66
Sample2  0  1  54
Sample3  1  0  76
Sample4  0  0  62
Sample5  1  1  59
Sample6  0  0  68
```
### Without bam file (useful for sc-multiome)
For some studies, like single cell multiome, generating a bam file is not convenient. In this case, we can use predefined peaks and openness to run the model. For each cluster, we can use the aggregated RNA-seq and ATAC-seq over cells from the same individual to generate pseudo-bulk gene expression and chromatin accessibility for each individual and run PECA. High-resolution peaks (i.e. separate peak-calling by MACS2 for each cluster) are helpful for regulatory network analysis compared to the default broad peaks from the cellRanger.  Peak calling for each cluster after clustering is recommended (for detail, please see the scREG package: https://github.com/Durenlab/RegNMF).
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 \
                                  --peak_file /full/path/to/peak.bed \
                                  --openness_file /full/path/to/openness.txt
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt --Design_Matrix covariates_knock.txt
```
Here the peak is a three-column bed file, and the openness is a tab-delimited text file where the first three columns are (chr start end). The 4th column is the openness of the 1st sample, the 5th column is the openness of the 2nd sample,.....
### Predefined motif-scan
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 \
                                  --peak_file /full/path/to/peak.bed \
                                  --openness_file /full/path/to/openness.txt
                                  --motif_file /full/path/to/MotifTarget.txt
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt --Design_Matrix covariates_knock.txt
```
Here the MotifTarget.txt is a three-column tab-delimited text file. We can obtain this file from Homer by following the script.
```
findMotifsGenome.pl region.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ./Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt
```
