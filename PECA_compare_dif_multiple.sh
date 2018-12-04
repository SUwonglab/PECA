#/bin/bash

Group1=$1
Group2=$2
organism=$3
# Examples
# Group1=Control
# Group2=Patients
# organism=human

mkdir ./Results/CompareGroup_${Group1}_${Group2}
cd ./Results/CompareGroup_${Group1}_${Group2}
cp ../../scr/FindSubNetwork_bi1.m ./.
cp ../../scr/FindSubNetwork_bi2.m ./.
cp ../../scr/dif_test.m ./.
cp ../../${Group1} ./.
cp ../../${Group2} ./.
sed "s/Sample1/${Group1}/g" ../../scr/PECA_dif_net_multiple.m|sed "s/Sample2/${Group2}/g"|sed "s/organism/${organism}/g" > PECA_dif_net_multiple.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_dif_net_multiple; exit"
cat *_specific_network.txt|grep -v log10P_Regulation|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledNetwork.txt
cat *_specific_module.txt|grep -v log10P_Regulation|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledModule.txt
echo GroupComparison ${Group1} with ${Group2} with PECA_dif_net_multiple done
