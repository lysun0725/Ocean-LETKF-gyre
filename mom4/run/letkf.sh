#!/bin/sh --login
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# letkf.sh <root> <TIME> <days> <ENS_NUM> <OBS_ENS> <NUM_DRIFTERS>

set -e

root=$1
root_run=${root}/run
#OBS_DIR=${root}/OBS
OBS_DIR=/lustre/lysun/models/Ocean-LETKF-gyre1-test3-025/mom4/OBS/
EXP_DIR=${root}/OUTPUT/
INPUT_INIT=${root}/INPUT_INIT

YYYYMMDDHH=$2
days=$3
ENS_NUM=$4
OBS_ENS=`printf %.3d $5`
NUM_DRIFTERS=$6

MEM3=1
USE_ALTIMETRY=0
oday=1  # observation time period

IYYYY=${YYYYMMDDHH:0:4}
IMM=${YYYYMMDDHH:4:2}
IDD=${YYYYMMDDHH:6:2}
IHH=${YYYYMMDDHH:8:2}
INN=00
ISS=00

ISLOT=01

date=/bin/date
sinc=${days} #` expr $day - $oday ` # "1" depends on the running period of nature run
sinc_units=days
NYYYY=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%Y`
NMM=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%m`
NDD=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%d`
NHH=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%H`

OBSDIR=${OBS_DIR}/${NYYYY}${NMM}${NDD}${NHH}/${OBS_ENS}/INPUT
echo "Read observation data from ${OBSDIR}..."
  #-----------------------------------------------------------------------------
  # (OPTIONAL),for simulation data, the following steps are necessary:
  # Converting observation netcdf to .dat file
  #-----------------------------------------------------------------------------
echo "Start coverting netcdf file to binary file...."
# Experiment 1: DRIFTERS DATA -> obsin.dat and obsin_drifters.dat
echo "Open the obs directory: ${OBSDIR}"
cd ${OBSDIR}/

if [ -e obsin_drifters.dat ]; then
  echo "obsin_drifters does exist"
else
  day2=`printf %.2d ${days}`
  cp -f ${root}/obs/NCEP_PROF/d2l_synth_${day2}.x ./
  ./d2l_synth_${day2}.x #Generating obsin_drifters.dat
  echo "obsin_drifters does not exit and we generate one by d2l_synth."
fi
	
echo 'Finished generating obsin_drifters.dat.'

while [ "${MEM3}" -le "${ENS_NUM}" ]
do 
  
  if [ "${MEM3}" -lt 100 ]; then
    MEM3="0${MEM3}"
  fi

  if [ "${MEM3}" -lt 10 ]; then 
    MEM3="0${MEM3}"
  fi

  echo "LETKF preparation step"
  echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"
  echo "I am member ${MEM3}"
  workdir_fcst=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}/RESTART
  workdir_inpt=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}/INPUT
  workdir_fcst_dr=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}/DRIFTERS
  workdir2=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/letkf
  OBSOPexe=obsop.001
  OBSOP_DRexe=obsop.DRIFTERS.001
  echo '*******letkf'
  mkdir -p ${workdir2}

  workdir=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/letkf_prep/${MEM3}/${ISLOT}
  echo '*******letkf_prep'
  mkdir -p ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/letkf_prep/
  echo '*******workdir'
  mkdir -p ${workdir}

  cd ${workdir} # pwd={workdir}

  #-----------------------------------------------------------------------------
  # FIRST, link the background model data
  #-----------------------------------------------------------------------------

  ln -f $workdir_fcst/ocean_temp_salt.res.nc     gs$ISLOT$MEM3.ocean_temp_salt.res.nc
  ln -f $workdir_fcst/ocean_velocity.res.nc      gs$ISLOT$MEM3.ocean_velocity.res.nc
  ln -f $workdir_fcst/ocean_sbc.res.nc           gs$ISLOT$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f $workdir_fcst/ocean_barotropic.res.nc  gs$ISLOT$MEM3.ocean_barotropic.res.nc
  fi
  
  cd ${workdir}
  ln -f ${workdir_fcst_dr}/drifters_inp.txt          gs$ISLOT$MEM3.drifters_inp.txt


  #STEVE: add 'fill value' to netcdf files for identification of missing values
#  cp $ncatted .
#  ncatted -O -a _FillValue,temp,o,f,-1.e+34 gs$MEM3.ocean_temp_salt.res.nc
#  ncatted -O -a _FillValue,salt,o,f,-1.e+34 gs$MEM3.ocean_temp_salt.res.nc
#  ncatted -O -a _FillValue,u,o,f,-1.e+34 gs$MEM3.ocean_velocity.res.nc
#  ncatted -O -a _FillValue,v,o,f,-1.e+34 gs$MEM3.ocean_velocity.res.nc
#  ncatted -O -a _FillValue,sea_lev,o,f,-1.e+34 gs$MEM3.ocean_sbc.res.nc
#  if [ "$USE_ALTIMETRY" -eq "1" ]; then
#    ncatted -O -a _FillValue,eta_t,o,f,-1.e+34 gs$MEM3.ocean_barotropic.res.nc
#  fi

  cp ${workdir}/gs$ISLOT$MEM3.ocean_temp_salt.res.nc ${workdir2}/anal$MEM3.ocean_temp_salt.res.nc
  cp ${workdir}/gs$ISLOT$MEM3.ocean_velocity.res.nc  ${workdir2}/anal$MEM3.ocean_velocity.res.nc
  cp ${workdir}/gs$ISLOT$MEM3.ocean_sbc.res.nc       ${workdir2}/anal$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    cp ${workdir}/gs$ISLOT$MEM3.ocean_barotropic.res.nc  ${workdir2}/anal$MEM3.ocean_barotropic.res.nc
  fi
  cp ${workdir}/gs$ISLOT${MEM3}.drifters_inp.txt      ${workdir2}/anal${MEM3}.drifters_inp.txt

  # Create grid.drifters_inp.txt for LETKF operator "set_common_drifters" (forcast)....
  cd ${workdir_fcst_dr}
  head -1 drifters_inp.txt > grd.drifters_inp.txt  
  cd ${workdir}

  # LUYU: This file is just for letkf_drifters_tools to set_common_drifters  
  if [ "${MEM3}" -eq 1 ]; then
    cp ${workdir_fcst_dr}/grd.drifters_inp.txt ${workdir2}/ 
    cp ${workdir_inpt}/grid_spec.nc ${workdir2}/
  fi
  #-----------------------------------------------------------------------------
  # SECOND, link the observations
  #-----------------------------------------------------------------------------
  cp -f $INPUT_INIT/INPUT_INIT/grid_spec.nc ./
  if [ ! -f grid_spec.nc ]; then
    echo "### NOTICE: grid_spec.nc does not exist. ###"
  fi

  cp  ${root}/obs/$OBSOP_DRexe .

  ln -f ${OBSDIR}/obsin_drifters.dat obsin_drifters.dat
  ln -f ${workdir}/gs$ISLOT$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs$ISLOT$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
  ln -f ${workdir}/gs$ISLOT$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f ${workdir}/gs$ISLOT$MEM3.ocean_barotropic.res.nc       gues.ocean_barotropic.res.nc
    echo "Linking: $altimetry_climatology_file to here..."
  fi
  ln -f ${workdir}/gs$ISLOT$MEM3.drifters_inp.txt         gues.drifters_inp.txt

  ####################################################################################################################################
  # Running LETKF Obs Operator executable to generate the observation innovations for each member at each timestep:
  ####################################################################################################################################
# $OBSOPexe -obsin $OBSDIR/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
  #STEVE: (perhaps make parallel and call with aprun)
  #$OBSOPexe
  $OBSOP_DRexe 

#  if [ -f "obsout.dat" ]; then
#    ln -f obsout.dat ${workdir2}/obs${ISLOT}${MEM3}.dat
#  else
#    echo "output obs2 formatted file not created by $OBSOPexe."
#    pwd
#    ls
#    echo "Exiting..."
#    exit 2
#  fi

  if [ -f "obsout_drifters.dat" ]; then
    ln -f obsout_drifters.dat ${workdir2}/obs${ISLOT}${MEM3}_drifters.dat
  else
    echo "output obs2_drifters formatted file not created by $OBSOP_DRexe."
    pwd
    ls
    echo "Exiting..."
    exit 2
  fi
  ln -f ${workdir}/gs${ISLOT}$MEM3.ocean_temp_salt.res.nc ${workdir2}/gs${ISLOT}$MEM3.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT}$MEM3.ocean_velocity.res.nc  ${workdir2}/gs${ISLOT}$MEM3.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT}$MEM3.ocean_sbc.res.nc       ${workdir2}/gs${ISLOT}$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc       ${workdir2}/gs${ISLOT}$MEM3.ocean_barotropic.res.nc
  fi
  ln -f ${workdir}/gs${ISLOT}$MEM3.drifters_inp.txt ${workdir2}/gs${ISLOT}$MEM3.drifters_inp.txt

  #STEVE: the hard-link limit is running out (65000 max), so best to delete unnecessary links
  rm -f grid_spec.nc
  rm -f ncatted
  
  MEM3=`expr $MEM3 + 1`

done #END loop of ENS_NUM

MEM3=$ENS_NUM
if [ "${MEM3}" -lt 100 ]; then
  MEM3="0${MEM3}"
fi

if [ "${MEM3}" -lt 10 ]; then 
  MEM3="0${MEM3}"
fi

echo "Running LETKF..."

LETKFexe=letkf.DRIFTERS.$MEM3
cp ${root}/letkf_drifters/$LETKFexe ${workdir2}
cd ${workdir2}

mpirun -n 20 $LETKFexe

###### LUYU: add perturbation to each analysis
MEM3=1
while [ "${MEM3}" -le "${ENS_NUM}" ]
do 
  MEM3=`printf %.3d $MEM3`
  cp anal${MEM3}.drifters_inp.nc analp${MEM3}.drifters_inp.nc

  MEM3=`expr ${MEM3} + 1`
done

#cp ${root_run}/perturbation/perturb_on_anal.x ./
#./perturb_on_anal.x

exit 0
