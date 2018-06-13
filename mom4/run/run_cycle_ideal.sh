#!/bin/sh --login
#SBATCH -n 20      
#SBATCH -t 20:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J letkf_dr_mom4p1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

IYYYYMMDDHH=1981011700
MYYYYMMDDHH=1981011700
EYYYYMMDDHH=1981030100
YYYYMMDDHH=$IYYYYMMDDHH

echo "Start running the experiment..."


EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run

ENS_NUM=40
days=10
isfirst=0
NUM_DRIFTERS=50
OBS_ENS=41


if [ "${isfirst}" -eq "1" ]; then
  sh pre_run.sh ${root} ${ENS_NUM} # generate obsop.001 obsop.DRIFTERS.001 letkf.DRIFTERS.005

  sh ${root_run}/perturbation/setup_ens_drift_perturbation.sh ${root} ${IYYYYMMDDHH} ${ENS_NUM} ${OBS_ENS} # Here we just let all the drifters start at the same positions and run for 10 days
fi

#NOTE: for this experiment the real experiment starts at 1981011000

while [ "${YYYYMMDDHH}" -lt "${EYYYYMMDDHH}" ]
do
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00 
  # sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst> <OBS_ENS>
  sh model_mom4p1.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${isfirst} ${OBS_ENS} ${MYYYYMMDDHH}
  
  wait

  if [ "${YYYYMMDDHH}" -ge "${MYYYYMMDDHH}" ]; then     
     # letkf.sh <root> <TIME> <days> <ENS_NUM> <OBS_ENS>
     echo "STARTS DOING DATA ASSIMILATION"
     sh letkf.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${OBS_ENS} ${NUM_DRIFTERS}
  else
     echo "NO DATA ASSIMILATION. PRE RUN."
  fi
  isfirst=0

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
























