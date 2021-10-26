#/bin/bash

#download the prior data 
cd ./Prior/
wget -O Opn_median_mm9.bed https://github.com/SUwonglab/PECA/raw/master/Prior/Opn_median_mm9.bed
wget -O Opn_median_mm10.bed https://github.com/SUwonglab/PECA/raw/master/Prior/Opn_median_mm10.bed
wget -O Opn_median_hg19.bed https://github.com/SUwonglab/PECA/raw/master/Prior/Opn_median_hg19.bed
wget -O Opn_median_hg38.bed https://github.com/SUwonglab/PECA/raw/master/Prior/Opn_median_hg38.bed
wget -O RE_gene_corr_mm9.bed https://github.com/SUwonglab/PECA/raw/master/Prior/RE_gene_corr_mm9.bed
wget -O RE_gene_corr_mm10.bed https://github.com/SUwonglab/PECA/raw/master/Prior/RE_gene_corr_mm10.bed
wget -O RE_gene_corr_hg19.bed https://github.com/SUwonglab/PECA/raw/master/Prior/RE_gene_corr_hg19.bed
wget -O RE_gene_corr_hg38.bed https://github.com/SUwonglab/PECA/raw/master/Prior/RE_gene_corr_hg38.bed
wget -O TFTG_corr_mouse.mat https://github.com/SUwonglab/PECA/raw/master/Prior/TFTG_corr_mouse.mat
wget -O TFTG_corr_human.mat https://github.com/SUwonglab/PECA/raw/master/Prior/TFTG_corr_human.mat
cd ../
mkdir ./Results/
mkdir ./Homer/
