#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

set -e

#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 06:30:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst>
set -e

EXP_NAME=gyre1
root=/lustre/lysun/models/Ocean-LETKF-${EXP_NAME}/mom4
root_run=${root}/run
EXP_DIR=${root}/OBS
INPUT_INI=${root}/INPUT_INIT

IYYYYMMDDHH=1981010100
EYYYYMMDDHH=1981123100

days=1

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00
MEM3=041 # ONLY have one ensemble, so set up MEM3 with 0('s)



while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do
echo "==================================================================="
echo "Running MOM4p1_drifters for NATURE RUN."
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

cd ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}/DRIFTERS
cp ${root_run}/drifters_code/read_drifters.f90 ./
cp ${root_run}/drifters_code/make_rd.sh

sh make_rd.sh

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
