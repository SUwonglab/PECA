# Comparative gene regulatory analysis via PECA
If you have three or more samples in each group, you can leverage PECA for performing comparative regulatory analysis. PECA is a powerful tool designed for analyzing bulk or pseudo-bulk RNA-seq and ATAC-seq data. It generates condition-specific transcriptional regulatory networks and identifies condition-specific driver regulators. Additionally, PECA allows for controlling covariates such as race, age, sex, and treatment to obtain more accurate results
## Simple run
```
sh PECA_multi.sh /full/path/to/All_sample_name.txt hg19 
sh PECA_compare_withoutENCODE.sh Group1.txt Group2.txt All_sample_name.txt
```
Here Group1.txt Group2.txt All_sample_name.txt are three one-column-text files containing sample names, which have to be consistent with the sample name given under the ./Input folder. For the input detail please see README.md. The results will be in folder ./Results_multiSample/CompareGroup_Group1_Group2/
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
### Without bam file (could be used for sc-multiome)
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
Motif-scan is the most time and space consuming step in PECA so that user can scan motifs first and give it to PECA as input to avoid doing motif-scan agian and agian for multiple samples or multiple runs. Here the MotifTarget.txt is a three-column tab-delimited text file. We can obtain this file from Homer by following the script.
```
findMotifsGenome.pl region.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ./Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt
```
## Figures and Tables
### Filtering networks
The condition-specific network may contain many TF-TG pairs so it may difficult to be visualized. Users can filter networks (i.e. Group1_specific_network.txt ) based on fold-change in regulation score difference (column 9), Activity of TF-TG regulation (ranging from 0-1, column 10), and/or TF-TG correlation (column 3). 
```
cat ${Group1}_specific_network.txt |awk 'NR>1'|awk 'BEGIN{OFS="\t"}{if($9>1.1) print $0}'|awk 'BEGIN{OFS="\t"}{if($10>0.5) print $0}'|awk 'BEGIN{OFS="\t"}{if($3>0.2) print $0}' > ${Group1}_specific_network_filter.txt
```
The same filtering could be done for Group1_specific_TFnetwork.txt and Group1_specific_module.txt.
### TF-TG network by Cytoscape
For the filtered TF-TG network, TF-TF network, or TF-TG module, users can plot condition-specific regulatory networks in Cytoscape by using the first three columns of information from the filter network files.
Users also can take a union of networks from two conditions, and visualize them together in the same network in Cytoscape. The Node_label.txt file can be used for node attributes and give colors for nodes based on the expression patterns. 
### Condition-specific driver TFs
The Group1_specific_TF_outDegree.txt and Group2_specific_TF_outDegree.txt reflect the importance of TFs in two conditions. Users may visualize the outdegree (# of TGs in the filtered network) in each condition.
### differential TF-TG regulation score
The Diff_TF_TG_pairs_scores_normalized.txt (covariates controlled) or Diff_TF_TG_pairs_scores_raw.txt could be used to visualize the differential regulatory patterns by clustering/heatmap for chosen TF-TG pairs (i.e. target genes of a specific TF, multiple TFs for a specific target gene, top ranking differential regulations).

