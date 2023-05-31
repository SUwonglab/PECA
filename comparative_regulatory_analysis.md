# Comparative gene regulatory analysis via PECA
If you have three or more samples in each group, you can leverage PECA for performing comparative regulatory analysis. PECA is a powerful tool designed for analyzing bulk or pseudobulk RNA-seq and ATAC-seq data. It generates condition-specific transcriptional regulatory networks and identifies condition-specific driver regulators. Additionally, PECA allows for controlling covariates such as race, age, sex, and treatment to obtain more accurate results
## Simple run
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt
```
Here Group1.txt Group2.txt All_sample_name.txt are three one-column-text-files contain sample names, which have to be consistent with the sample name given under the ./Input folder. For the input detail please see README.md.
## Additional Features and Customization
### Covariate control
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt --Design_Matrix covariates_knock.txt
```
Write all covariates into a design matrix (i.e. covariates_knock.txt). Here the design matrix (all element have to be **numbers**) is a **tab delimited** text file as this example:
```
Name  Is_treatment  Is_male  Age
Sample1  1  1  66
Sample2  0  1  54
Sample3  1  0  76
Sample4  0  0  62
Sample5  1  1  59
Sample6  0  0  68
```
### Without bam file
For some study, like single cell multiome, generating bam file is to convinient. In this case, we can use predefined peak and openness to run the model.
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 \
                                  --peak_file /full/path/to/peak.bed \
                                  --openness_file /full/path/to/openness.txt
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt --Design_Matrix covariates_knock.txt
```
Here the peak is a three-columns bed file, and the openness is a tab delimited text file where the first three columns are (chr start end). The 4th column is the openness of the 1-st sample, 5-th column is the openness of 2-nd sample,.....
