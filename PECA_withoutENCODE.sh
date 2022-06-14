#/bin/bash

# Run PECA2 based_on_your_own_data without using ENCODE public data infromation: updated June 13th 2022
# usage PECA_withoutENCODE.sh sampleName_file genome
# sampleNmae_file: a txt file: one sample name in each line, consistent with the sample names in Input folder. 
# step 1: call peak from bam file
# step 2: motif binding
# step 3: calculate opn
# step 4: corr+dist
# step 5: score(Exp,binding,Opn,weight)



input_file=$1
genome=$2

numCore=`nproc`
echo $numCore core will be used_in motif analysis

resultFolder=`basename ${input_file}|awk 'BEGIN{FS="."}{print $1}'`
mkdir ./Results/${resultFolder}
cd ./Results/${resultFolder}/
echo step 1: call peak from bam file....

if [ `echo $genome|grep hg|wc -l` -gt 0 ]
then
	species=hs
	speciesFull=human
else
	species=mm
	speciesFull=mouse
fi

> peak_raw.bed
while read line
do
  macs2 callpeak -t ../../Input/${line}.bam -f BAM -n ${line} -g ${species} --nomodel --shift -100 --extsize 200
  cat ${line}_peaks.narrowPeak|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn > ${line}_peak.txt
  cat ${line}_peak.txt >> peak_raw.bed
done < ${input_file} 
sort -k1,1 -k2,2n peak_raw.bed|mergeBed > peak.bed
cat peak.bed|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn > region.txt

echo step 2: motif binding....
findMotifsGenome.pl region.txt ${genome} ./. -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../../Homer/ -p $numCore > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt

echo step 3: calculate opn...
cat region.txt|cut -f 1-3|sort -k1,1 -k2,2n > region.bed
cat region.bed|awk 'BEGIN{OFS="\t"}{print $1,$2-500000,$3+500000}'| awk '{$2=$2<0?1:$2}1'|tr ' ' '\t'>background.bed

cat region.bed >openness_allSamples.bed
while read input
do
cat region.bed >a1
samtools bedcov region.bed  ../../Input/${input}.bam >read.bed
samtools bedcov ../../Data/genome_100k_hg19.bed ../../Input/${input}.bam > back_read
bedtools intersect -a background.bed -b back_read -wa -wb -sorted > background
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$6-$5}' background|sed 's/\t/\_/1'|sed 's/\t/\_/1'|awk '{a[$1]+=$2;b[$1]+=$3}END{for(i in a){print i,a[i],b[i];}}'|tr ' ' '\t'|tr '_' '\t'|sort -k1,1 -k2,2n > background_fc
paste -d '\t' read.bed  background_fc|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*$9/(($3-$2)*($8+5))}'|sed 's/\t/\_/1'|sed 's/\t/\_/1' > openness_${input}.bed
mkdir ../${input}/
cat openness_${input}.bed > ../${input}/openness.bed
cp region.bed ../${input}/
cp region.txt ../${input}/
ln -s MotifTarget.txt ../${input}/
rm a1
rm read.bed
cat openness_${input}.bed|cut -f 2 > a2
paste -d '\t' openness_allSamples.bed  a2 > a1
mv a1 openness_allSamples.bed
done < ${input_file} 
awk 'BEGIN{OFS="\t"}{s=0; for (i=4; i<=NF;i++) s += $i; print $1,$2,$3,s/(NF-3)}' openness_allSamples.bed > openness_mean.bed

echo step 4: Prior....
cat ../../Data/promoter_${genome}.bed| awk 'BEGIN{OFS="\t"}{print $1,$2-200000,$3+200000,$4}'| awk '{$2=$2<0?1:$2}1'|tr ' ' '\t' > promoter_200k.bed
bedtools intersect -a region.txt -b  promoter_200k.bed -wa -wb|cut -f 4,8 > peak_gene.txt
bedtools intersect -a region.txt -b  promoter_200k.bed -wa -wb|awk '{print $3-$7+100000}'|awk '{$1=$1<0?-$1:$1}1' > dis
cat ${input_file} > SampleNameFile
sed "s/speciesFull/${speciesFull}/g" ../../scr/prior_calculation.m > prior_calculation.m
matlab -nodisplay -nosplash -nodesktop -r "prior_calculation; exit"
paste -d '\t' peak_gene.txt dis > a
paste -d '\t' a corr > RE_TG.bed
while read input
do
cat RE_TG.bed > ../${input}/peak_gene_100k_corr.bed
cp TFTG_corr.mat ../${input}/
cp openness_mean.bed ../${input}/
done < ${input_file} 
cd ../

echo step 5: Network....
while read input
do
cd ./${input}/
cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b openness_mean.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
mkdir Enrichment
cat openness2.bed|awk 'BEGIN{OFS="\t"}{print $1,($2+0.5)/($3+0.5)}'|sort -k2nr|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{if ($3-$2 < 2000) print $0}'|head -10000 > ./Enrichment/region.bed
sed "s/species/${speciesFull}/g" ../../scr/mf_collect.m > ./Enrichment/mf_collect.m 
cd ./Enrichment/
findMotifsGenome.pl region.bed ${genome} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir ../../../Homer/ -p $numCore
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "mf_collect; exit"
cd ../
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
echo ${input} PECA completed
cd ../
done < ${input_file} 
