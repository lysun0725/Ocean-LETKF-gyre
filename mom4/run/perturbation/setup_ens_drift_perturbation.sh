# LUYU: this code is to setup the ensembles of drifters.
# setup_ini_drift_perturbation.sh <root> <TIME> <ENS_NUM> <OBS_NUM>

set -e

root=$1
root_run=${root}/run
OBS_DIR=${root}/OBS
EXP_DIR=${root}/OUTPUT/OUTPUT
INPUT_INI=${root}/INPUT_INIT

YYYYMMDDHH=$2
ENS_NUM=$3
OBS_ENS=`printf %.3d $4`

MEM3=1

mkdir -p ${INPUT_INI}/INIT_DRIFTERS

cp ${OBS_DIR}/${YYYYMMDDHH}/${OBS_ENS}/INPUT/drifters_inp.nc ${INPUT_INI}/INIT_DRIFTERS

while [ "${MEM3}" -le ${ENS_NUM} ]
do

  MEM3=`printf %.3d ${MEM3}`
  cp -f ${INPUT_INI}/INIT_DRIFTERS/drifters_inp.nc ${INPUT_INI}/INIT_DRIFTERS/drifters_inp.${MEM3}.nc

  MEM3=`expr ${MEM3} + 1`

done

cp -f ${root_run}/perturbation/setup_ini_drift_perturbation.x ${INPUT_INI}/INIT_DRIFTERS/
cd ${INPUT_INI}/INIT_DRIFTERS/
#./setup_ini_drift_perturbation.x
