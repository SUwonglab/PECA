#!/bin/bash

# Function to display script usage
display_usage() {
    echo "Usage: $0 <input_name> <genome> [--peak_file <value>] [--motif_file <value>] [--openness_file <value>]"
    echo "  <input_name>   : full path to a text file contains the name of all samples, no spatial chractors in the file name, /home/user/abc.txt"
    echo "  <genome>   : hg19,hg38,mm9,mm10"
    echo "  --peak_file (optional)    : full path to a bed file"
    echo "  --motif_file (optional)   : full path to a motif file, which is generated from homer findMotifsGenome.pl, cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt"
    echo "  --openness_file (optional): full path to text file, the first three columns are chr start end the same as peak file, from the 4th column put the openness file as the same order as input file"
    echo "  --each_sample_net (optional): False (default), or True"
}
# Check if the required inputs are provided
if [[ -z $1 || -z $2 ]]; then
    display_usage
    exit 1
fi

# Assign the required inputs to variables
input_file=$1
genome=$2
# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$3"
    case $key in
    --peak_file)
        peak_file="$4"
        shift
        shift
        ;;
    --motif_file)
        motif_file="$4"
        shift
        shift
        ;;
    --openness_file)
        openness_file="$4"
        shift
        shift
        ;;
    --each_sample_net)
        each_sample_net="$4"
        shift
        shift
        ;;
    *)
        break
        ;;
    esac
done
# Display the inputs
echo "input_file: $input_file"
echo "genome: $genome"
if [[ -v peak_file ]]; then
  echo "peak_file(optional): $peak_file"
else
  echo "no peak provided, we will do peak calling by MACS2"
fi
if [[ -v motif_file ]]; then
  echo "motif_file(optional): $motif_file"
else
  echo "no motif_scan_file provided, we will do motif scanning by Homer"
fi
if [[ -v openness_file ]]; then
  echo "openness_file(optional): $openness_file"
else
    echo "no openness_file provided, we will canculated openness from bam files"
fi
###folder
resultFolder=`basename ${input_file}|awk 'BEGIN{FS="."}{print $1}'`
mkdir ./Results_multiSample
mkdir ./Results_multiSample/${resultFolder}
cd ./Results_multiSample/${resultFolder}/
if [ `echo $genome|grep hg|wc -l` -gt 0 ]
then
	species=hs
	speciesFull=human
else
  species=mm
  speciesFull=mouse
fi
####step1: peak
if [[ -v peak_file ]]; then
  cat $peak_file|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn |sort -k1,1 -k2,2n > region.txt
else
  echo step 1: call peak from bam file....
  > peak_raw.bed
  while read line
  do
    macs2 callpeak -t ../../Input/${line}.bam -f BAM -n ${line} -g ${species} --nomodel --shift -100 --extsize 200
    cat ${line}_peaks.narrowPeak|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn > ${line}_peak.txt
    cat ${line}_peak.txt >> peak_raw.bed
  done < ${input_file} 
  sort -k1,1 -k2,2n peak_raw.bed|mergeBed > peak.bed
  cat peak.bed|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'|grep -v rand |grep -v chrUn |sort -k1,1 -k2,2n > region.txt
fi
####step2: motif 
if [[ -v motif_file ]]; then
  ln -s $motif_file MotifTarget.txt
else
  findMotifsGenome.pl region.txt ${genome} ./. -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../../Homer/ -p $numCore > MotifTarget.bed
  cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
  rm MotifTarget.bed
  rm motifFindingParameters.txt
fi
#####step3: openness
if [[ -v openness_file ]]; then
  cat $openness_file > openness_allSamples.bed
  let i=3;
  while read input
  do
    mkdir ../${input}/
    let i=i+1;
    cat openness_allSamples.bed|awk -v i="$i" 'BEGIN{OFS="\t"}{print $1,$2,$3,$i}'|sed 's/\t/_/1'|sed 's/\t/_/1' > openness_${input}.bed
    cat openness_${input}.bed > ../${input}/openness.bed
    cp region.txt ../${input}/
  done < ${input_file} 
else
  cat region.txt|cut -f 1-3 > region.bed
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
    rm a1
    rm read.bed
    cat openness_${input}.bed|cut -f 2 > a2
    paste -d '\t' openness_allSamples.bed  a2 > a1
    mv a1 openness_allSamples.bed
  done < ${input_file} 
fi
awk 'BEGIN{OFS="\t"}{s=0; for (i=4; i<=NF;i++) s += $i; print $1,$2,$3,s/(NF-3)}' openness_allSamples.bed > openness_mean.bed
####Step 4: prior and activity analysis
echo step 4: Prior....
cat ../../Data/promoter_${genome}.bed| awk 'BEGIN{OFS="\t"}{print $1,$2-200000,$3+200000,$4}'| awk '{$2=$2<0?1:$2}1'|tr ' ' '\t' > promoter_200k.bed
bedtools intersect -a region.txt -b  promoter_200k.bed -wa -wb|cut -f 4,8 > peak_gene.txt
bedtools intersect -a region.txt -b  promoter_200k.bed -wa -wb|awk '{print $3-$7+100000}'|awk '{$1=$1<0?-$1:$1}1' > dis
cat ${input_file} > SampleNameFile
cp ../../scr/mfbs.m ./.
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
cd ../../
####Step 5: network for each sample
if [[(-v each_sample_net) && ($each_sample_net == "True" )]]; then
echo step 5: Network....
cd ./Results_multiSample/
while read input
do
cd ./${input}/
cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b openness_mean.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
cp ../../scr/mfbs.m ./.
sed "s/toreplace/${input}/g" ../../scr/PECA_network_${genome}.m > PECA_network.m
matlab -nodisplay -nosplash -nodesktop -r "addpath('../$resultFolder/');PECA_network; exit"
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
fi