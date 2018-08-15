#/bin/bash

Sample1=mESC
Sample2=RAd2
# Examples
#Sample1 = mESC     
#Sample1 = RAd2  

mkdir ./Results/Compare_${Sample1}_${Sample2}
cd ./Results/Compare_${Sample1}_${Sample2}
cp ../../scr/FindSubNetwork_bi1.m ./.
sed "s/Sample1/${Sample1}/g" ../../scr/PECA_dif_net.m|sed "s/Sample2/${Sample2}/g" > PECA_dif_net.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_dif_net; exit"
echo Compare ${Sample1} with ${Sample2} with PECA_dif_net done