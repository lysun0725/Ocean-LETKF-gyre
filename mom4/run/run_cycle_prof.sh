#!/bin/sh --login
#SBATCH -n 80      
#SBATCH -t 30:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J letkf_dr_mom4p1_ts
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

IYYYYMMDDHH=1981041500
MYYYYMMDDHH=1981011600
EYYYYMMDDHH=1981050100
YYYYMMDDHH=$IYYYYMMDDHH

echo "Start running the experiment..."


EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run

ENS_NUM=40
days=1
isfirst=0 # NOTE!!! Don't forget to check before run the experiment!!!
NUM_DRIFTERS=50
OBS_ENS=41
ISLOT=1
ddays=1


#if [ "${isfirst}" -eq "1" ]; then
#  sh pre_run.sh ${root} ${ENS_NUM} # generate obsop.001 obsop.DRIFTERS.001 letkf.DRIFTERS.040 

#  sh ${root_run}/perturbation/setup_ens_drift_perturbation.sh ${root} ${IYYYYMMDDHH} ${ENS_NUM} ${OBS_ENS} # generate ensemble members on drifters
#fi

while [ "${YYYYMMDDHH}" -lt "${EYYYYMMDDHH}" ]
do
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00 
  # sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst> <OBS_ENS> <MYYYYMMDDHH>
  sh model_mom4p1_prof.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${isfirst} ${OBS_ENS} ${MYYYYMMDDHH} ${ISLOT} ${ddays}
  
  wait

  if [ "${YYYYMMDDHH}" -ge "${MYYYYMMDDHH}" ]; then
    echo "The SLOT is ${ISLOT}."

    if [ "${ISLOT}" -eq "${ddays}" ]; then     
     # letkf.sh <root> <TIME> <days> <ENS_NUM> <OBS_ENS>
     echo "STARTS DOING DATA ASSIMILATION"
     sh letkf_prof.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${OBS_ENS} ${NUM_DRIFTERS}     
     # ISLOT=0
    else
     echo "SKIP DATA ASSIMILATION"
    fi
  else
     echo "NO DATA ASSIMILATION. PRE RUN."
  fi

  #ISLOT=`expr ${ISLOT} + 1`
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
























