#!/bin/sh --login
#SBATCH -n 80      
#SBATCH -t 30:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J letkf_dr_mom4p1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

set -e

IYYYYMMDDHH=1981011500
EYYYYMMDDHH=1981011600
YYYYMMDDHH=$IYYYYMMDDHH

work_root=/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/OBS2
orig_root=/lustre/lysun/models/Ocean-LETKF-gyre1-test3-025/mom4/OBS
INPUT_dir=/lustre/lysun/models/Ocean-LETKF-gyre1-test3/mom4/INPUT_INIT/INPUT_INIT
days=1


while [ "${YYYYMMDDHH}" -lt "${EYYYYMMDDHH}" ]
do

  cd ${work_root}

  mkdir ${YYYYMMDDHH}
  mkdir ${YYYYMMDDHH}/041
  mkdir ${YYYYMMDDHH}/041/RESTART

  workdir=${work_root}/${YYYYMMDDHH}/041/RESTART
  origdir=${orig_root}/${YYYYMMDDHH}/041/RESTART

  cd ${workdir}

# create the input namelist for 'ocean_temp_salt'
  cat >input.nml <<!
    &regrid_nml
       src_data = '${origdir}/ocean_temp_salt.res.nc'
       src_grid = '${orig_root}/${YYYYMMDDHH}/041/INPUT/grid_spec.nc'
       dst_grid = '${INPUT_dir}/grid_spec.nc'
       dst_data = 'ocean_temp_salt.res.nc'
       num_flds = 2
       fld_name = 'temp', 'salt'
       fld_pos = 'T', 'T'
       vector_fld = .false.
       use_source_vertical_grid = .true.
       apply_mask = .true.
       debug      = .false./
!

  cp ../../../../run/regrid.x ./
  ./regrid.x
  rm -f input.nml

# create the input namelist for 'ocean_temp_salt'
  cat >input.nml <<!
    &regrid_nml
       src_data = '${origdir}/ocean_velocity.res.nc'
       src_grid = '${orig_root}/${YYYYMMDDHH}/041/INPUT/grid_spec.nc'
       dst_grid = '${INPUT_dir}/grid_spec.nc'
       dst_data = 'ocean_velocity.res.nc'
       num_flds = 2
       fld_name = 'u', 'v'
       fld_pos = 'T', 'T'
       vector_fld = .false.
       use_source_vertical_grid = .true.
       apply_mask = .true.
       debug      = .false./
!

  ./regrid.x
 
  cp ${origdir}/drifters_inp.nc ./ 

  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00 

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

exit 0
