#/bin/bash

Sample1=$1   
Sample2=$2 
organism=$3
# Examples
#Sample1 = mESC     
#Sample1 = RAd2  
#organism=human

mkdir ./Results/Compare_${Sample1}_${Sample2}
cd ./Results/Compare_${Sample1}_${Sample2}
cp ../../scr/FindSubNetwork_bi1.m ./.
sed "s/Sample1/${Sample1}/g" ../../scr/PECA_dif_net.m|sed "s/Sample2/${Sample2}/g"|sed "s/organism/${organism}/g" > PECA_dif_net.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_dif_net; exit"
cat *_specific_network.txt|grep -v Activity|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledNetwork.txt
cat *_specific_module.txt|grep -v Activity|cut -f 1-3|awk '{$3=$3>0?1:-1}1'|tr ' ' '\t'|sort|uniq > PooledModule.txt
echo Compare ${Sample1} with ${Sample2} with PECA_dif_net done
