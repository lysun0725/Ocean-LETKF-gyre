#!/bin/sh
#SBATCH -n 80         
#SBATCH -t 20:30:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_control
set -e

EXP_NAME=gyre1-test3
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
OBS_DIR=${root}/OBS
CON_DIR=${root}/CONTROL2
EXP_DIR=${root}/OUTPUT
INPUT_INI=${root}/INPUT_INIT/1981053100

IYYYYMMDDHH=1981010100
EYYYYMMDDHH=1981043000

days=1
ENS_NUM=1 # This is a nature run, so only one ensemble.
isfirst=1
NUM_DRIFTERS=50

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00
PEM3=041
MEM3=CON # ONLY have one ensemble, so set up MEM3 with 0('s)


mkdir -p ${CON_DIR}

while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do
echo "==================================================================="
echo "Running MOM4p1_drifters for NATURE RUN."
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${CON_DIR}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${CON_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${CON_DIR}/${IYYYY}${IMM}${IDD}${IHH}

    echo "I am member ${MEM3}"
    workdir=${CON_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}
    mkdir -p ${workdir}
    cd ${workdir}

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART
    mkdir -p ${workdir}/DRIFTERS

    if [ "$isfirst" -eq "1" ]; then
      cp ${INPUT_INI}/${PEM3}/RESTART/input.nml ${workdir}/INPUT/        
      cp ${INPUT_INI}/${PEM3}/RESTART/*_table ${workdir}/INPUT/
       
      # Now we are going to setup the control initial condition by taking the mean of the ensembles
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_age.res.nc ${workdir}/INPUT/ocean_age.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_barotropic.res.nc ${workdir}/INPUT/ocean_barotropic.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_bih_friction.res.nc ${workdir}/INPUT/ocean_bih_friction.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_con_temp.res.nc ${workdir}/INPUT/ocean_con_temp.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_density.res.nc ${workdir}/INPUT/ocean_density.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_passive.res.nc ${workdir}/INPUT/ocean_passive.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_psom_moments.res.nc ${workdir}/INPUT/ocean_psom_moments.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_thickness.res.nc ${workdir}/INPUT/ocean_thickness.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_velocity.res.nc ${workdir}/INPUT/ocean_velocity.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_sbc.res.nc ${workdir}/INPUT/ocean_sbc.res.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/ocean_velocity_advection.res.nc ${workdir}/INPUT/ocean_velocity_advection.res.nc
      #ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/tau.nc ${workdir}/INPUT/tau.nc
      ncea -O ${EXP_DIR}/${IYYYYMMDDHH}/model/0??/INPUT/drifters_inp.nc ${workdir}/INPUT/drifters_inp.nc

      echo 'Finish computing the initial condition of CONTROL RUN'
      cp ${EXP_DIR}/${IYYYYMMDDHH}/model/001/INPUT/tau.nc ${workdir}/INPUT/tau.nc ${workdir}/INPUT/

      sed -i "106s/.*/      use_this_module=.true./" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
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
      ln ${CON_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/input.nml ${workdir}/INPUT/
      ln ${CON_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/drifters_inp.nc ${workdir}/INPUT/
      ln ${CON_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/ocean_*.res.nc ${workdir}/INPUT/         
      ln ${CON_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/*_table ${workdir}/INPUT/

      ln ${CON_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/INPUT/tau.nc ${workdir}/INPUT/
    fi # END IF on isfirst

    #LUYU: For model gyre1
    sed -i "247s/.*/      days = $days/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "248s/.*/      date_init = ${IYYYY},${IMM},${IDD},${IHH},${INN},${ISS}/" ${workdir}/INPUT/input.nml

    # LUYU: For model gyre1
    ln ${INPUT_INI}/${PEM3}/INPUT/grid_spec.nc ${workdir}/INPUT
    ln ${INPUT_INI}/${PEM3}/INPUT/gotmturb.inp ${workdir}/INPUT  
    
    cp ${root_run}/mom4p1_solo_run.csh ${workdir}
 
    # sbatch mom4p1_solo_run.csh <full path to working directory> $MEM3
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}
 
    cp ${root_run}/drifters_code/* ${workdir}/DRIFTERS/

    cd ${workdir}/DRIFTERS/
    ln ${workdir}/INPUT/drifters_inp.nc INPUT_drifters_inp.nc
    sh make_rd.sh # LUYU: If you are going to change the number of processor running the model, you need to modify the corresponding number in make_rd.sh
    sh make_cd.sh
    
    if [ -e ${workdir}/RESTART/drifters_inp.nc ]; then
      rm -rf ${workdir}/RESTART/drifters_inp.nc
      cp drifters_inp.nc ${workdir}/RESTART/
    fi

    cd ${workdir}  

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

done # END the loop of time

echo 'Naturual Finish.'
exit 0
