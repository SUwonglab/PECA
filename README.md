# PECA

Introduction
PECA is a software for inferring context specific gene regulatory network from paired gene expression and chromatin accessibility data.

Install
bash install.sh

Run PECA
1,edit the three input in run_PECA.sh files, sampleName, genome and single-end/pair-end.
2,Put the input files in folder named Input. Three files: ${SampleName}.txt, ${SSampleName}.bam, ${SSampleName}.bam.bai
3,bash run_PECA.sh
Example: run_PECA.sh

The results will be ./Results/${SampleName}/
${SampleName}_network.txt is the tissue specific network
TFTG_regulationScore.txt is regulation strength for the all TF to TG. Row name and column name are TFName.txt and TGName.txt


