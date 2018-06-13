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
EXP_DIR=${root}/OUTPUT/
EXP_DIR2=${root}/OUTPUT/
OBS_ENS=041
days=1
COUNT=1

PYYYYMMDDHH=1981010100
IYYYYMMDDHH=1981010200
EYYYYMMDDHH=1981010300
YYYYMMDDHH=$PYYYYMMDDHH



echo "Starting to setup files for computing errors..."

if [ -e ${EXP_DIR2}/norm ]; then
  rm -rf ${EXP_DIR2}/norm
fi

mkdir -p ${EXP_DIR2}/norm
cd ${EXP_DIR2}/norm

echo "Collecting data from pre run..."


while [ "${YYYYMMDDHH}" -lt "${IYYYYMMDDHH}" ]
do
  
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00

  COUNT=`printf %.4d "${COUNT}"`
  echo "COUNT = ${COUNT}"

  ln ${EXP_DIR}/${TYYYY}${TMM}${TDD}${THH}/model/DATA/MEAN_DR_${TMM}${TDD}.nc drif_me_${COUNT}.nc

  ln ${OBS_DIR}/${TYYYY}${TMM}${TDD}${THH}/${OBS_ENS}/DRIFTERS/drifters_inpc.txt obs_dr_${COUNT}.txt

  ln ${EXP_DIR2}/${IYYYYMMDDHH}/letkf/grd.drifters_inp.txt grd.drifters_inp_${COUNT}.txt

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

echo "Collecting data from Data Assimilation window...."

days=1
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
  
  if [ -f ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/andr_me.txt ]; then
    ln ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/andr_me.txt andr_me_${COUNT}.txt
  else
    echo "${TYYYY}${TMM}${TDD}${THH} dose not have andr_me.txt..."
  fi


  if [ -f ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/gsdr_me.txt ]; then
    ln ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/gsdr_me.txt gsdr_me_${COUNT}.txt
  else
    echo "${TYYYY}${TMM}${TDD}${THH} dose not have gsdr_me.txt..."
  fi

  if [ -f ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/andr_sp.txt ]; then  
    ln ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/andr_sp.txt andr_sp_${COUNT}.txt
  else
    echo "${TYYYY}${TMM}${TDD}${THH} dose not have andr_sp.txt..."
  fi 

  ln ${EXP_DIR2}/${TYYYY}${TMM}${TDD}${THH}/letkf/grd.drifters_inp.txt grd.drifters_inp_${COUNT}.txt

  ln ${OBS_DIR}/${TYYYY}${TMM}${TDD}${THH}/${OBS_ENS}/DRIFTERS/drifters_inpc.txt obs_dr_${COUNT}.txt

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

# copy the executive file and grads file in this directory.




