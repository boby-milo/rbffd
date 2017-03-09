#!/bin/bash -l
 
#SBATCH -A snic2015-6-159
#SBATCH -p core -n 1
#SBATCH -t 72:00:00
#SBATCH -J millovanovic
#SBATCH -C usage_mail

module load matlab
cd ~/phs/rbffd/slask/T0.1
matlab -nojvm < BSeuCall2DbasketReference.m > Out.txt

