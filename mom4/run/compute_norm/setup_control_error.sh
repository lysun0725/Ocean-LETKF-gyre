#!/bin/sh --login
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

set -e

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run

OBS_DIR=${root}/OBS
EXP_DIR=${root}/CONTROL

OBS_ENS=041
days=1
COUNT=1

PYYYYMMDDHH=1981010100
EYYYYMMDDHH=1981021400
YYYYMMDDHH=$PYYYYMMDDHH



echo "Starting to setup files for computing errors..."

if [ -e ${EXP_DIR}/norm ]; then
  rm -rf ${EXP_DIR}/norm
fi

mkdir -p ${EXP_DIR}/norm
cd ${EXP_DIR}/norm

echo "Collecting data from pre run..."

while [ "${YYYYMMDDHH}" -le "${EYYYYMMDDHH}" ]
do
  
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00

  COUNT=`printf %.4d "${COUNT}"`
  echo "COUNT = ${COUNT}"

  ln ${EXP_DIR}/${TYYYY}${TMM}${TDD}${THH}/CON/RESTART/drifters_inp.nc drif_me_${COUNT}.nc

  ln ${OBS_DIR}/${TYYYY}${TMM}${TDD}${THH}/${OBS_ENS}/DRIFTERS/drifters_inpc.txt obs_dr_${COUNT}.txt
  ln ${root}/OUTPUT/OUTPUT.3/${EYYYYMMDDHH}/letkf/grd.drifters_inp.txt grd.drifters_inp_${COUNT}.txt

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


  COUNT=`expr ${COUNT} + 1`
done




