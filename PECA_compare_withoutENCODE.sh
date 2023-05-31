#!/bin/bash

# Function to display script usage
display_usage() {
    echo "Usage: $0 <Group1_name_file> <Group2_name_file> [--Design_Matrix <value>]"
    echo "  <Group1_name_file>   : A text file contains the name of samples from group1"
    echo "  <Group2_name_file>   : A text file contains the name of samples from group2"
    echo "  <input_name>   : the file given in PECA_multi.sh as input, a text file contain all sample names"
    echo "  --Design_Matrix (optional)    : Design matrix to provide all covairates"
}
# Check if the required inputs are provided
if [[ -z $1 || -z $2 ]]; then
    display_usage
    exit 1
fi

Group1=$1
Group2=$2
AllSampleName=`basename ${3}|awk 'BEGIN{FS="."}{print $1}'`

if [[ $4 == "--Design_Matrix" ]]; then
    if [[ -n $5 && ! $5 =~ ^-[^-] ]]; then
        DesignMat=$5
    else
        echo "Missing or invalid value for --optional_input"
        display_usage
        exit 1
    fi
fi

####pbuild folder and repare scripts

mkdir ./Results_multiSample/CompareGroup_noENCODE_${Group1}_${Group2}
cd ./Results_multiSample/CompareGroup_noENCODE_${Group1}_${Group2}
cp ../../scr/FindSubNetwork_bi1.m ./.
cp ../../scr/FindSubNetwork_bi2.m ./.
cp ../../scr/dif_test.m ./.
cp ../../${Group1} ./.
cp ../../${Group2} ./.
cp ../../${DesignMat} ./.
#########run comparison
if [[ -v DesignMat ]]; then
  echo "using ${DesignMat} as design matrix to control covariates"
  sed "s/Group1/${Group1}/g" ../../scr/compareNet.m|sed "s/Group2/${Group2}/g"|sed "s/AllSample/${AllSampleName}/g"|sed "s/DesignMat/${DesignMat}/g" > compareNet.m
else
  echo "no covariates to control"
  sed "s/Group1/${Group1}/g" ../../scr/compareNet.m|sed "s/Group2/${Group2}/g"|sed "s/AllSample/${AllSampleName}/g"  > compareNet.m
fi

matlab -nodisplay -nosplash -nodesktop -r "compareNet; exit"
cat *_specific_network.txt|grep -v log10P_Regulation|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledNetwork.txt
cat *_specific_module.txt|grep -v log10P_Regulation|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledModule.txt
cat ${Group1}_specific_network.txt |awk 'NR>1'|awk 'BEGIN{OFS="\t"}{if($7>1.1) print $0}'|awk 'BEGIN{OFS="\t"}{if($8>0.5) print $0}'|awk 'BEGIN{OFS="\t"}{if($3>0.2) print $0}' > ${Group1}_specific_network_filter.txt
cat ${Group1}_specific_network_filter.txt|cut -f 1|sort|uniq -c|awk 'BEGIN{OFS="\t"}{print $2,$1}'|sort -k2nr > ${Group1}_specific_TF_outDegree.txt
cat ${Group2}_specific_network.txt |awk 'NR>1'|awk 'BEGIN{OFS="\t"}{if($7>1.1) print $0}'|awk 'BEGIN{OFS="\t"}{if($8>0.5) print $0}'|awk 'BEGIN{OFS="\t"}{if($3>0.2) print $0}' > ${Group2}_specific_network_filter.txt
cat ${Group2}_specific_network_filter.txt|cut -f 1|sort|uniq -c|awk 'BEGIN{OFS="\t"}{print $2,$1}'|sort -k2nr > ${Group2}_specific_TF_outDegree.txt
paste -d '\t' Diff_TF_TG_pairs.txt Diff_net_raw.txt > Diff_TF_TG_pairs_scores.txt
paste -d '\t' Diff_TF_TG_pairs.txt Diff_net_normalized.txt > Diff_TF_TG_pairs_scores_normalized.txt
echo GroupComparison ${Group1} with ${Group2} with compareNet done
