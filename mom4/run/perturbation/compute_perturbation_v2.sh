#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

set -e

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
INPUT_INI=${root}/OUTPUT

IYYYYMMDDHH=$1
EYYYYMMDDHH=$1

days=30
ENS_NUM=40
OBS_ENS=41
OBS_ENS=`printf %.3d ${OBS_ENS}`
variable1=velocity
variable2=temp_salt

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00

while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do

echo "==================================================================="
echo "Running MOM4p1_drifters for NATURE RUN."
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/DATA ]; then
    rm -rf ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/DATA
  fi

  mkdir ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/DATA
  OUTDIR=${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/DATA

  # Compute ensemble mean 
  ncea ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/0??/RESTART/drifters_inp.nc $OUTDIR/MEAN_${IMM}${IDD}.nc

  ncdiff $OUTDIR/MEAN_${IMM}${IDD}.nc ${root}/OBS/${IYYYY}${IMM}${IDD}${IHH}/${OBS_ENS}/RESTART/drifters_inp.nc ${OUTDIR}/DIFF_${IMM}${IDD}.nc

  MEM3=1
  while [ "${MEM3}" -le "${ENS_NUM}" ]
  do 
    MEM3=`printf %.3d ${MEM3}`
    echo "Hello from member ${MEM3}."

    # Compute ensemble members as perturbations from zero mean
    ncbo --op_typ=subtract ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}/RESTART/drifters_inp.nc $OUTDIR/MEAN_${IMM}${IDD}.nc $OUTDIR/PERT_${MEM3}_${IMM}${IDD}.nc
    echo "Finish computing ensemble perturbation."

    # Square perturbations
    ncbo --op_typ=mlt $OUTDIR/PERT_${MEM3}_${IMM}${IDD}.nc $OUTDIR/PERT_${MEM3}_${IMM}${IDD}.nc $OUTDIR/PERT2_${MEM3}_${IMM}${IDD}.nc  
    echo "Finish computing ensemble perturbation square."

    MEM3=`expr ${MEM3} + 1`
  done # END the loop of ENS_NUM

  ncea $OUTDIR/PERT2_???_${IMM}${IDD}.nc $OUTDIR/VAR_${IMM}${IDD}.nc

  ncap2 -s 'sd_p=sqrt(positions*40/39)' $OUTDIR/VAR_${IMM}${IDD}.nc $OUTDIR/SD_${IMM}${IDD}.nc
  echo 'Finish computing sd'

  # Update the time
  date=/bin/date
  sinc=$days
  sinc_units='days'
  TYYYY=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%Y`
  TMM=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%m`
  TDD=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%d`
  THH=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%H`
  TNN=$INN
  TSS=$ISS

  IYYYY=${TYYYY}
  IMM=${TMM}
  IDD=${TDD}
  IHH=${THH}
  INN=${TNN}
  ISS=${TSS}

done
