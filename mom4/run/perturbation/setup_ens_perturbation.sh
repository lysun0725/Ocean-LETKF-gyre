set -e

EXP_NAME=gyre1
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
INPUT_INI=${root}/INPUT_INIT

ENS_NUM=41
MEM3=1

mkdir -p ${INPUT_INI}/INIT_TAU

cp ${INPUT_INI}/INPUT_INIT/tau.nc ${INPUT_INI}/INIT_TAU/

while [ "${MEM3}" -le "${ENS_NUM}" ]
do

  MEM3=`printf %.3d ${MEM3}`
  cp ${INPUT_INI}/INIT_TAU/tau.nc ${INPUT_INI}/INIT_TAU/tau.${MEM3}.nc

  MEM3=`expr ${MEM3} + 1`

done

cp ${root_run}/perturbation/setup_ini_perturbation.x ${INPUT_INI}/INIT_TAU
cd ${INPUT_INI}/INIT_TAU/
./setup_ini_perturbation.x

