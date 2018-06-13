#!/bin/sh
#SBATCH -n 80         
#SBATCH -t 30:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# This code is for generating initial conditions for different ensembles.

set -e

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
ENS_NUM=41
EXP_DIR=${root}/INPUT_INIT_SODA
INPUT_INI=${EXP_DIR}/INPUT_INIT
INPUT_TAU=${EXP_DIR}/INIT_TAU

IYYYYMMDDHH=1981010100
EYYYYMMDDHH=1981012900
days=30
isfirst=1

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00
MEM3=2

while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do
echo "==================================================================="
echo "Running MOM4p1_drifters"
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  
  while [ "${MEM3}" -le "$ENS_NUM" ]
  #while [ "${MEM3}" -le 1 ]
  do
    if [ "${MEM3}" -lt 100 ]; then
      MEM3=0${MEM3}
    fi
    if [ "${MEM3}" -lt 10 ]; then 
      MEM3=0${MEM3}
    fi
    
    if [ -e ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3} ]; then
      rm -rf ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}
    fi

    mkdir -p ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}

    echo "I am member ${MEM3}"
    workdir=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}
    cd ${workdir}

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART

    if [ "$isfirst" -eq "1" ]; then
      cp ${INPUT_INI}/input.nml ${workdir}/INPUT/
      cp ${INPUT_INI}/ocean_*.res.nc ${workdir}/INPUT/
         
      cp ${INPUT_INI}/*_table ${workdir}/INPUT/
      cp ${INPUT_TAU}/tau.${MEM3}.nc ${workdir}/INPUT/tau.nc    

      sed -i "106s/.*/      use_this_module=.false./" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    else
      # Update the date for the previous analysis cycle
      date=/bin/date
      pinc=$days
      pinc_units='days ago'
      PYYYY=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%Y`
      PMM=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%m`
      PDD=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%d`
      PHH=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%H`
      PNN=$INN
      PSS=$ISS

      # LINK all the analytic files in letkf of previous step or from the restart file from previous step. ****FINISH Later
      # First link the files in PREVIOUS INPUT folder:
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/input.nml ${workdir}/INPUT/
      #ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/drifters_inp.nc ${workdir}/INPUT/
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/ocean_*.res.nc ${workdir}/INPUT/         
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/*_table ${workdir}/INPUT/
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/tau.nc ${workdir}/INPUT/

    fi # END IF on isfirst


    # LUYU: Following is the code to update the initial data, find the position of this in input.nml and then update.
    # LUYU: For model gyre1
    ln ${INPUT_INI}/grid_spec.nc ${workdir}/INPUT
    ln ${INPUT_INI}/gotmturb.inp ${workdir}/INPUT 

    cp ${root_run}/mom4p1_solo_run.csh ${workdir}

    # sbatch mom4p1_solo_run.csh <full path to working directory> $MEM3
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}

    cd ${workdir}  
   
    MEM3=`expr $MEM3 + 1`
    
  done # END the loop of ENS_NUM

isfirst=0

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
done # END the loop of dates

echo 'Naturual Finish.'
exit 0
