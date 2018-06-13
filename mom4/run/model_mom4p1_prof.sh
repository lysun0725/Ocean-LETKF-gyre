#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst> <OBS_ENS> <MYYYYMMDDHH>
set -e

root=$1
root_run=${root}/run
EXP_DIR=${root}/OUTPUT/
INPUT_INI=${root}/INPUT_INIT
INPUT_INI_ENS2=${root}/INPUT_INIT/1981053100
INPUT_INI_ENS=${root}/OUTPUT/

YYYYMMDDHH=$2
days=$3
ENS_NUM=$4
isfirst=$5
OBS_ENS=$6
MYYYYMMDDHH=$7
ISLOT=$8
ddays=$9

if [ "${OBS_ENS}" -lt 100 ]; then
  OBS_ENS=0${OBS_ENS}
fi
if [ "${OBS_ENS}" -lt 10 ]; then 
  MEM3=0${OBS_ENS}
fi


IYYYY=${YYYYMMDDHH:0:4}
IMM=${YYYYMMDDHH:4:2}
IDD=${YYYYMMDDHH:6:2}
IHH=${YYYYMMDDHH:8:2}
INN=00
ISS=00
MEM3=1


echo "==================================================================="
echo "Running MOM4p1_drifters"
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  
  while [ "${MEM3}" -le "$ENS_NUM" ]
  do
    if [ "${MEM3}" -lt 100 ]; then
      MEM3=0${MEM3}
    fi
    if [ "${MEM3}" -lt 10 ]; then 
      MEM3=0${MEM3}
    fi

    echo "I am member ${MEM3}"
    workdir=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}
    mkdir -p ${workdir}
    cd ${workdir}
    echo "The workdir is ${workdir}...."

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART

    if [ "$isfirst" -eq "1" ]; then

      IMEM3=${MEM3}
       
      if [ "${IMEM3}" -ge "${OBS_ENS}" ]; then
        IMEM3=`expr ${IMEM3} + 1`        
        IMEM3=`printf %.3d ${IMEM3}`
      fi

      echo "Obtain initial data from ${IMEM3}"

      cp ${INPUT_INI_ENS2}/${IMEM3}/RESTART/input.nml ${workdir}/INPUT/
      cp ${INPUT_INI_ENS2}/${IMEM3}/RESTART/ocean_*.res.nc ${workdir}/INPUT/
      cp ${INPUT_INI_ENS2}/${IMEM3}/RESTART/*_table ${workdir}/INPUT/

      # drifters_inp need to be copied from OBS/INPUT/
      echo "The OBS_ENS is ${OBS_ENS}"

      # turn off the drifters module
      sed -i "106s/.*/      use_this_module=.false./" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
         
    else
      # Update the date for the previous analysis cycle
      date=/bin/date
      if [ "${YYYYMMDDHH}" -eq "${MYYYYMMDDHH}" ]; then
        pinc=1
      else 
        pinc=${days}
      fi
      pinc_units='days ago'
      PYYYY=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%Y`
      PMM=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%m`
      PDD=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%d`
      PHH=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%H`
      PNN=$INN
      PSS=$ISS
      echo "Obtainning initial data from ${PYYYY}${PMM}${PDD}${PHH}..."

      # LINK all the analytic files in letkf of previous step or from the restart file from previous step. ****FINISH Later
      # First link the files in PREVIOUS INPUT folder:
      cp ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/input.nml ${workdir}/INPUT/input.nml
      cp ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/*.res.nc ${workdir}/INPUT/         
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/*_table ${workdir}/INPUT/

      if [ "${ISLOT}" -eq "1" ]; then

        # Next, check whether the previous step has LETKF:
        echo "Checking whether the Standard-LETKF Analysis from previous step exists. If yes, link as initial conditions......"
        workdir_analysis=${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/letkf
        if [ -e ${workdir_analysis} ]; then

          ##### TEMPERATURE and SALINITY file:
          if [ -f ${workdir_analysis}/anal${MEM3}.ocean_temp_salt.res.nc ]; then
            echo "The anal${MEM3}.ocean_temp_salt.res.nc exists... and start linking...."
            rm -rf ${workdir}/INPUT/ocean_temp_salt.res.nc
            cp -f ${workdir_analysis}/anal${MEM3}.ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc
          else
            echo "ANALYSIS FILE DOES NOT EXIST: anal${MEM3}.ocean_temp_salt.res.nc..."
            exit 1
          fi

          ##### VELCOCITY file:
          if [ -f ${workdir_analysis}/anal${MEM3}.ocean_velocity.res.nc ]; then
            echo "The anal${MEM3}.ocean_velocity.res.nc exists... and start linking...."
            rm -rf ${workdir}/INPUT/ocean_velocity.res.nc
            cp -f ${workdir_analysis}/anal${MEM3}.ocean_velocity.res.nc ${workdir}/INPUT/ocean_velocity.res.nc
          else
            echo "ANALYSIS FILE DOES NOT EXIST: anal${MEM3}.ocean_velocity.res.nc..."
            exit 1
          fi

          ##### SBC file:
          if [ -f ${workdir_analysis}/anal${MEM3}.ocean_sbc.res.nc ]; then
            echo "The anal${MEM3}.ocean_sbc.res.nc exists... and start linking...."
            rm -rf ${workdir}/INPUT/ocean_sbc.res.nc
            cp -f ${workdir_analysis}/anal${MEM3}.ocean_sbc.res.nc ${workdir}/INPUT/ocean_sbc.res.nc
          else
            echo "ANALYSIS FILE DOES NOT EXIST: anal${MEM3}.ocean_sbc.res.nc..."
          fi
       
        fi
      fi
    fi # END IF on isfirst

    #LUYU: For model gyre1
    sed -i "106s/.*/      use_this_module=.false./" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "247s/.*/            days = $days/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "248s/.*/            date_init = ${IYYYY},${IMM},${IDD},${IHH},${INN},${ISS}/" ${workdir}/INPUT/input.nml
    sed -i '111s/^#"ocean_model"/"ocean_model"/' ${workdir}/INPUT/diag_table # Don't forget to add "-i" to save


    # LUYU: For model gyre1
    ln ${INPUT_INI_ENS2}/${MEM3}/INPUT/grid_spec.nc ${workdir}/INPUT
    ln ${INPUT_INI_ENS2}/${MEM3}/INPUT/gotmturb.inp ${workdir}/INPUT 
    ln ${INPUT_INI_ENS2}/${MEM3}/INPUT/tau.nc       ${workdir}/INPUT 
        
    cp ${root_run}/mom4p1_solo_run.csh ${workdir}

    # sbatch mom4p1_solo_run.csh <full path to working directory> $MEM3
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}
    
    cd ${workdir} 
   
    MEM3=`expr $MEM3 + 1`
    
  done # END the loop of ENS_NUM

echo 'Naturual Finish.'
exit 0
