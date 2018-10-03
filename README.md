# PECA

## Introduction:

PECA is a software for inferring context specific gene regulatory network from paired gene expression and chromatin accessibility data.
Please cite: 

Duren, Zhana, et al. "Modeling gene regulation from paired expression and chromatin accessibility data." Proceedings of the National Academy of Sciences 114.25 (2017): E4914-E4923.

## Quickly start:

wget https://github.com/durenzn/PECA/archive/master.zip

unzip master.zip

cd PECA-master/

bash install.sh

bash run_PECA.sh

## Install:

bash install.sh

## Run PECA:

1, edit the two input in run_PECA.sh files, sampleName and genome (line 10 and 11 in run_PECA.sh).

2, Put the input files in folder named ./Input. Three files: ${SampleName}.txt, ${SampleName}.bam, ${SampleName}.bam.bai.

${SampleName}.txt is gene expression file containing two columns (tab delimited), gene Symbol and FPKM (or TPM). 

${SampleName}.bam is chromatin accessibility data, DNase-seq or ATAC-seq. 

${SampleName}.bam.bai is the index file of bam file. 

Please see the example of RAd4 in the ./Input directory.

3, bash run_PECA.sh

Example: run_PECA.sh

The results will be ./Results/${SampleName}/ .
${SampleName}_network.txt is the tissue specific network.

TFTG_regulationScore.txt is regulation strength for the all TF to TG. Row name and column name are TFName.txt and TGName.txt

## Run PECA_net_dif:
If you have two samples and want to compare the two samples at network level, please do it by following steps:

1, Run PECA on two samples one by one by runing PECA

2, Edit the Sample1 and Sample2 in run_PECA_compare_dif.sh (line 3 and 4 in run_PECA_compare_dif.sh)

3, bash run_PECA_compare_dif.sh

The results will be ./Results/Compare_${Sample1}_${Sample2}. Containing six files:  

specific network of two samples: ${Sample1}_specific_network.txt and ${Sample2}_specific_network.txt

common network of two samples: ${Sample1}_${Sample2}_common_network.txt 

specific module of two networks:  ${Sample1}_specific_module.txt and ${Sample2}_specific_module.txt

common module of two samples: ${Sample1}_${Sample2}_common_module.txt 

## Run PECA_net_dif_multiple:
If you have two conditions (multiple samples in each conditions) and want to compare the two conditions at network level, please do it by following steps:

1, Run PECA on all the samples from two conditions one by one by "bash run_PECA.sh"

2, Write the sample names of Group1 into one file named $Group1, and the same in $Group2 (eg. create one file named "Control" and put the sample names of first condition to this file, create one file named "Case" and put the sample names of second condition to this file)

3, Edit the Group1 and Group2 in run_PECA_compare_dif_multiple.sh (line 3 and 4 in run_PECA_compare_dif_multiple.sh, eg. Group1=Control; Group2=Case )

4, bash run_PECA_compare_dif_multiple.sh

The results will be ./Results/CompareGroup_${Group1}_${Group2}. Containing six files:  

specific network of two conditions: ${Group1}_specific_network.txt and ${Group2}_specific_network.txt

common network of two conditions: ${Group1}_${Group2}_common_network.txt 

specific module of two conditions:  ${Group1}_specific_module.txt and ${Group2}_specific_module.txt

common module of two conditions: ${Group1}_${Group2}_common_module.txt

## Requirements:

* macs2

* bedtools

* homer

* samtools

* Matlab



