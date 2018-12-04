# PECA

## Introduction:

PECA is a software for inferring context specific gene regulatory network from paired gene expression and chromatin accessibility data.
Please cite: 

Duren, Zhana, et al. "Modeling gene regulation from paired expression and chromatin accessibility data." Proceedings of the National Academy of Sciences 114.25 (2017): E4914-E4923.

## Quickly start:

wget https://github.com/SUwonglab/PECA/archive/master.zip

unzip master.zip

cd PECA-master/

bash install.sh

bash PECA.sh ${sampleName} ${genome}

## Install:

bash install.sh

## Run PECA:

Run PECA by following two steps:

### Step 1: Input 
Put the input files in folder named ./Input. Three files: ${SampleName}.txt, ${SampleName}.bam, ${SampleName}.bam.bai.

${SampleName}.txt is gene expression file containing two columns (tab delimited), gene Symbol and FPKM (or TPM). 

${SampleName}.bam is chromatin accessibility data, DNase-seq or ATAC-seq. 

${SampleName}.bam.bai is the index file of bam file. 

Note that all the three files should have same before-dot-file-name ${SampleName},only difference is after dot ".txt", ".bam" or ".bam.bai". Please see the example of RAd4 in the ./Input directory.

### Step 2: Run 
sh PECA.sh ${sampleName} ${genome}

Example: sh PECA.sh RAd4 mm9

The results will be ./Results/${SampleName}/ .
${SampleName}_network.txt is the tissue specific network.

TFTG_score.txt is regulation strength for the all TF to TG. Each row represent one TF and each column represents one target gene. Higher value represents higher possibility of regulation.

CRB_pval.txt is the Chromatin regulators' (CR) binding site matrix, each column represent one CR, each row represent one region, the values are p-values.

## Run PECA_net_dif:
If you have two samples and want to compare the two samples at network level, please do it by following steps:

1, Prepare two networks: Run PECA on two samples one by one by "sh PECA.sh ${sampleName} ${genome}"

2, Run:  sh PECA_compare_dif.sh ${Sample1} ${Sample2} ${Organism}

Example: sh PECA_compare_dif.sh K562 GM12878 human ; sh PECA_compare_dif.sh mESC RAd4 mouse

The results will be ./Results/Compare_${Sample1}_${Sample2}. Containing six files:  

specific network of two samples: ${Sample1}_specific_network.txt and ${Sample2}_specific_network.txt

common network of two samples: ${Sample1}_${Sample2}_common_network.txt 

specific module of two networks:  ${Sample1}_specific_module.txt and ${Sample2}_specific_module.txt

common module of two samples: ${Sample1}_${Sample2}_common_module.txt 

Files PooledNetwork.txt or PooledModuole.txt can be used to visualize the network by cytoscype, and the node lable is given in file Node_lable.txt. "1" and "-1" in PooledNetwork.txt or PooledModuole.txt represent "Activation" and "Repression" respectively. "1" and "2" in Node_lable.txt represent the gene is Sample1 specific or Sample2 specific.

## Run PECA_net_dif_multiple:
If you have two conditions (multiple samples in each conditions) and want to compare the two conditions at network level, please do it by following steps:

1, Prepare networks: Run PECA on all the samples from two conditions one by one by "sh PECA.sh ${sampleName} ${genome}"

2, Construct lables: Write the sample names of Group1 and Group2 into text files named $Group1 and $Group2, respectively. (eg. create one text file named "Control" and put the sample names of one condition to this file, create other text file named "Case" and put the names of the other condition to this file. Note that the sample name files contain one sample name per line )

3, Run: sh PECA_compare_dif_multiple.sh $Group1 $Group2 ${Organism}
Example： sh PECA_compare_dif_multiple.sh Control Case human
 
The results will be ./Results/CompareGroup_${Group1}_${Group2}. Containing six files:  

specific network of two conditions: ${Group1}_specific_network.txt and ${Group2}_specific_network.txt

common network of two conditions: ${Group1}_${Group2}_common_network.txt 

specific module of two conditions:  ${Group1}_specific_module.txt and ${Group2}_specific_module.txt

common module of two conditions: ${Group1}_${Group2}_common_module.txt

Files PooledNetwork.txt or PooledModuole.txt can be used to visualize the network by cytoscype, and the node lable is given in file Node_lable.txt. "1" and "-1" in PooledNetwork.txt or PooledModuole.txt represent "Activation" and "Repression" respectively. "1" and "2" in Node_lable.txt represent the gene is Group1 specific or Group2 specific.

## Requirements:

* Matlab

* macs2

* homer

* samtools

* bedtools



