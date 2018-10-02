#/bin/bash

Group1=stem1
Group2=stem2
# Examples
# Group1=Control
# Group2=Patients

mkdir ./Results/CompareGroup_${Group1}_${Group2}
cd ./Results/CompareGroup_${Group1}_${Group2}
cp ../../scr/FindSubNetwork_bi1.m ./.
cp ../../scr/dif_test.m ./.
cp ../../${Group1} ./.
cp ../../${Group2} ./.
sed "s/Sample1/${Group1}/g" ../../scr/PECA_dif_net_multiple.m|sed "s/Sample2/${Group2}/g" > PECA_dif_net_multiple.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_dif_net_multiple; exit"
echo GroupComparison ${Group1} with ${Group2} with PECA_dif_net_multiple done
