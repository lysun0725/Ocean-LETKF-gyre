#!/bin/sh --login
#SBATCH -n 80      
#SBATCH -t 30:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J letkf_dr_mom4p1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

IYYYYMMDDHH=1981011600
MYYYYMMDDHH=1981011600
EYYYYMMDDHH=1981050100
YYYYMMDDHH=$IYYYYMMDDHH

echo "Start running the experiment..."

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4/OUTPUT/EXP8
root_run=${root}/run

days=1

while [ "${YYYYMMDDHH}" -lt "${EYYYYMMDDHH}" ]
do
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00 

  echo "${YYYYMMDDHH}"

  awk 'NR%2==1' ${root}/${TYYYY}${TMM}${TDD}${THH}/letkf/gsdr_me.txt > ${root}/${TYYYY}${TMM}${TDD}${THH}/letkf/gsdr_me_1.txt
  

  date=/bin/date
  sinc=$days
  sinc_units=days
  NYYYY=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%Y`
  NMM=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%m`
  NDD=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%d`
  NHH=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%H`
  NNN=$TNN
  NSS=$TSS

  YYYYMMDDHH=${NYYYY}${NMM}${NDD}${NHH}

done

echo "NATURAL END."

exit 0
























