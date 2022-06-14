#/bin/bash

# PECA2 v3.0.1 updated May 27th 2019
# step 1: call peak from bam file
# step 2: motif binding
# step 3: calculate opn
# step 4: corr+dist
# step 5: score(Exp,binding,Opn,weight)

input=$1
genome=$2

# Examples
#input = RAd4     your input file RAd4.txt (gene expression data), RAd4.bam (chromatin accessibility data), and RAd4.bam.bai (bam index file) are under folder./Input
#genome = mm9  Your refference genome is mm9, software is only available for mm9, mm10, hg19 and hg38 currently 

mkdir ./Results/${input}
cd ./Results/${input}/
echo step 1: call peak from bam file....
if [ `echo $genome|grep hg|wc -l` -gt 0 ]
then
	species=hs
	speciesFull=human
else
	species=mm
	speciesFull=mouse
fi
macs2 callpeak -t ../../Input/${input}.bam -f BAM -n ${input} -g ${species} --nomodel --shift -100 --extsize 200
cat ${input}_peaks.narrowPeak|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn > region.txt

echo step 2: motif binding....
findMotifsGenome.pl region.txt ${genome} ./. -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../../Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt

echo step 3: calculate opn...
cat region.txt|cut -f 1-3|sort -k1,1 -k2,2n > region.bed
cat region.bed|awk 'BEGIN{OFS="\t"}{print $1,$2-500000,$3+500000}'| awk '{$2=$2<0?1:$2}1'|tr ' ' '\t'>background.bed
cat region.bed >a1
samtools bedcov region.bed  ../../Input/${input}.bam >read.bed
samtools bedcov ../../Data/genome_100k_${genome}.bed ../../Input/${input}.bam > back_read
bedtools intersect -a background.bed -b back_read -wa -wb -sorted > background
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$6-$5}' background|sed 's/\t/\_/1'|sed 's/\t/\_/1'|awk '{a[$1]+=$2;b[$1]+=$3}END{for(i in a){print i,a[i],b[i];}}'|tr ' ' '\t'|tr '_' '\t'|sort -k1,1 -k2,2n > background_fc
paste -d '\t' read.bed  background_fc|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*$9/(($3-$2)*($8+5))}'|sed 's/\t/\_/1'|sed 's/\t/\_/1' > openness.bed 
rm a1
rm back*
rm read.bed
cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b ../../Prior/Opn_median_${genome}.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
mkdir Enrichment
cat openness2.bed|awk 'BEGIN{OFS="\t"}{print $1,($2+0.5)/($3+0.5)}'|sort -k2nr|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{if ($3-$2 < 2000) print $0}'|head -10000 > ./Enrichment/region.bed
sed "s/species/${speciesFull}/g" ../../scr/mf_collect.m > ./Enrichment/mf_collect.m 
cd ./Enrichment/
findMotifsGenome.pl region.bed ${genome} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir ../../../Homer/
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "mf_collect; exit"
cd ../


echo step 4: Prior....
bedtools intersect -a region.bed -b ../../Prior/RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
bedtools intersect -a region.bed -b ../../Prior/Enhancer_RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>>peak_gene_100k_corr.bed
cat ../../Prior/TFTG_corr_${speciesFull}.mat > TFTG_corr.mat

echo step 5: Network....
cp ../../scr/mfbs.m ./.
sed "s/toreplace/${input}/g" ../../scr/PECA_network_${genome}.m > PECA_network.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_network; exit"
echo region > CRbinding_region
cat region.txt|cut -f 4 > CRbinding_region1
cat CRbinding_region CRbinding_region1> CRbinding_region2
paste -d '\t' CRbinding_region2 CR_binding_pval.txt > CRB_pval.txt
rm CRbinding_region*
rm CR_binding_pval.txt
cat TGName.txt |tr '\n' '\t'> TG_head
echo >> TG_head
cat TG_head TFTG_regulationScore.txt > TG_score
cat TG_head TFTG_regulationScore_norm.txt > TG_score_norm
echo TFName > TG_1
cat TG_1 TFName.txt >TG_2
paste -d '\t' TG_2 TG_score > TFTG_score.txt
paste -d '\t' TG_2 TG_score_norm > TFTG_score_norm.txt
rm TG_*
echo ${input} PECA done
