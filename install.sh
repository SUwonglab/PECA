#/bin/bash

#download the prior data 
wget http://web.stanford.edu/~zduren/PECA/Thresholding-Based%20SVD_files/Prior.tar.gz
tar -zxvf Prior.tar.gz
rm Prior.tar.gz
mkdir ./Results/
mkdir ./Homer/
