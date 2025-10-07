
## Enviroment setup and install
Install MACS2 if not available.
allocate 128GB resource interactve job on quartz (replace your slurm-account-name)
```
srun -p debug -A slurm-account-name --time=01:00:00 --mem=128GB --pty bash
```
create conda enviroment
```
ml conda
conda create -n PECA
```
actiavte conda enviroment
```
conda activate PECA
```
install macs2 using bioconda
```
conda install -c bioconda macs2
```
check all requirements
```
ml matlab
ml homer
ml samtools
ml bedtools
macs2 --help
```
Install PECA
```
wget https://github.com/SUwonglab/PECA/archive/master.zip
unzip master.zip
cd PECA-master/
bash install.sh
```
test PECA
```
sh PECA.sh RAd4 mm10
```
## run PECA for multiple samples

```
srun -p debug -A slurm-account-name --time=01:00:00 --mem=128GB --pty bash
```
Input
textfile

