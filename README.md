# PECA

Introduction:
PECA is a software for inferring context specific gene regulatory network from paired gene expression and chromatin accessibility data.
Please cite: 

Duren, Zhana, et al. "Modeling gene regulation from paired expression and chromatin accessibility data." Proceedings of the National Academy of Sciences 114.25 (2017): E4914-E4923.

## Install:

bash install.sh

## Run PECA:

1, edit the two input in run_PECA.sh files, sampleName and genome.

2, Put the input files in folder named ./Input. Three files: ${SampleName}.txt, ${SSampleName}.bam, ${SSampleName}.bam.bai

3, bash run_PECA.sh

Example: run_PECA.sh

The results will be ./Results/${SampleName}/ .
${SampleName}_network.txt is the tissue specific network
TFTG_regulationScore.txt is regulation strength for the all TF to TG. Row name and column name are TFName.txt and TGName.txt

## Quickly start:

wget https://github.com/durenzn/PECA/archive/master.zip

unzip master.zip

cd PECA-master/

bash install.sh

bash run_PECA.sh

## Requirements:

* macs2

* bedtools

* homer

* samtools

* Matlab



