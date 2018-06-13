#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 12:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst>
set -e

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
INPUT_INI=${root}/INPUT_INIT
INPUT_TAU=${root}/INPUT_INIT/INIT_TAU3

IYYYYMMDDHH=1981010100
EYYYYMMDDHH=1981063000

days=30
ENS_NUM=41
isfirst=1

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00

while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do
echo "==================================================================="
echo "Running MOM4p1_drifters for PRE RUN."
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}
  MEM3=1

  while [ "${MEM3}" -le "${ENS_NUM}" ]
  do

    MEM3=`printf %.3d ${MEM3}`
    echo "I am member ${MEM3}"
    workdir=${INPUT_INI}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}
    mkdir -p ${workdir}
    cd ${workdir}

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART

    if [ "$isfirst" -eq "1" ]; then
      cp ${INPUT_INI}/INPUT_INIT/input.nml ${workdir}/INPUT/
      cp ${INPUT_INI}/INPUT_INIT/ocean_*.res.nc ${workdir}/INPUT/
         
      cp ${INPUT_INI}/INPUT_INIT/*_table ${workdir}/INPUT/     

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
      ln ${INPUT_INI}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/input.nml ${workdir}/INPUT/
      ln ${INPUT_INI}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/ocean_*.res.nc ${workdir}/INPUT/         
      ln ${INPUT_INI}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/*_table ${workdir}/INPUT/

    fi # END IF on isfirst

    #LUYU: For model gyre1
    sed -i "246s/.*/      months = 0/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "247s/.*/      days = $days/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "248s/.*/      date_init = ${IYYYY},${IMM},${IDD},${IHH},${INN},${ISS}/" ${workdir}/INPUT/input.nml

    # LUYU: For model gyre1
    ln ${INPUT_INI}/INPUT_INIT/grid_spec.nc ${workdir}/INPUT
    ln ${INPUT_INI}/INPUT_INIT/gotmturb.inp ${workdir}/INPUT 
    cp ${INPUT_TAU}/tau.${MEM3}.nc       ${workdir}/INPUT/tau.nc

    cp ${root_run}/mom4p1_solo_run.csh ${workdir}
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}

    MEM3=`expr ${MEM3} + 1`

  done 

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

done
