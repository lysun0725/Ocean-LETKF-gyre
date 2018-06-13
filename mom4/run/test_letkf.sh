#!/bin/sh --login
#SBATCH -n 80         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

YYYYMMDDHH=1981011600

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4

ENS_NUM=40
days=1
NUM_DRIFTERS=50
OBS_ENS=41

sh letkf_prof.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${OBS_ENS} ${NUM_DRIFTERS}  
