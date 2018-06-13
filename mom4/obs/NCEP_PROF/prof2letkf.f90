PROGRAM prof2letkf

 !USE common_mom4 ! Use check subroutine
 USE netcdf
 USE common
 USE params_obs

 IMPLICIT NONE

!SUBROUTINE read_grid(infile,v3d,v2d)
  INTEGER :: i,j,k,m
  INTEGER :: ncid,istat,varid,dimid
  
  !INTEGER,PARAMETER :: nlon=192
  !INTEGER,PARAMETER :: nlat=189
  !INTEGER,PARAMETER :: nlev=31
  INTEGER :: nlon, nlat, nlev, dodebug=1
  REAL(r_size), ALLOCATABLE :: lon(:)
  REAL(r_size), ALLOCATABLE :: lat(:)
  REAL(r_size), ALLOCATABLE :: lev(:)
  REAL(r_size), ALLOCATABLE :: salt(:,:,:)
  REAL(r_size), ALLOCATABLE :: temp(:,:,:)
  REAL(r_size), ALLOCATABLE :: u(:,:,:)
  REAL(r_size), ALLOCATABLE :: v(:,:,:)
  REAL(r_size), ALLOCATABLE :: sst(:,:)
  REAL(r_size), ALLOCATABLE :: sss(:,:)

  REAL(r_size) :: time, err(1)
  REAL(r_sngl) :: wk(6)
  INTEGER, DIMENSION(5) :: ilat=(/ 26, 29, 32, 35, 38/) 

  CHARACTER (len=*), PARAMETER :: tsfile="ocean_temp_salt.res.nc"
  CHARACTER (len=*), PARAMETER :: grfile="grid_spec.nc"
!  CHARACTER (len=*), PARAMETER :: uvfile="ocean_velocity.res.nc"
!  CHARACTER (len=*), PARAMETER :: sffile="ocean_sbc.res.nc"
  CHARACTER (len=*), PARAMETER :: outfile="obsin.dat"

  call check(NF90_OPEN(grfile,NF90_NOWRITE,ncid))
  print *, 'finish open the file: grid_spec.nc.'
  ! Read Dimensions
  call check(NF90_INQ_DIMID(ncid,'grid_x_T',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlon))
  call check(NF90_INQ_DIMID(ncid,'grid_y_T',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlat))
  call check(NF90_INQ_DIMID(ncid,'zt',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlev))
    
  ALLOCATE(lon(nlon))
  ALLOCATE(lat(nlat))
  ALLOCATE(lev(nlev))
  ALLOCATE(salt(nlon,nlat,nlev))
  ALLOCATE(temp(nlon,nlat,nlev))
  ALLOCATE(u(nlon,nlat,nlev))
  ALLOCATE(v(nlon,nlat,nlev))
  ALLOCATE(sst(nlon,nlat))
  ALLOCATE(sss(nlon,nlat))  

  ! Read Data
  call check(NF90_INQ_VARID(ncid,'grid_x_T',varid))
  call check(NF90_GET_VAR(ncid,varid,lon))
  print *, 'finish setting up grid_x_T.'
  call check(NF90_INQ_VARID(ncid,'grid_y_T',varid))
  call check(NF90_GET_VAR(ncid,varid,lat))
  print *, 'finish setting up grid_y_T.'
  call check(NF90_INQ_VARID(ncid,'zt',varid))
  call check(NF90_GET_VAR(ncid,varid,lev))
  print *, 'finish setting up zt.'
  call check(NF90_CLOSE(ncid))
  
  call check(NF90_OPEN(tsfile,NF90_NOWRITE,ncid))
  print *, 'finish open the file: ocean_temp_salt.res.nc.'
  call check(NF90_INQ_VARID(ncid,'Time',varid))
  call check(NF90_GET_VAR(ncid,varid,time))
  print *, 'finish setting up Time.'
  call check(NF90_INQ_VARID(ncid,'salt',varid))
  call check(NF90_GET_VAR(ncid,varid,salt))
  print *, 'finish setting up salt.'
  call check(NF90_INQ_VARID(ncid,'temp',varid))
  call check(NF90_GET_VAR(ncid,varid,temp))
  print *, 'finish setting up temp.' 
  call check(NF90_CLOSE(ncid))

!  call check(NF90_OPEN(uvfile,NF90_NOWRITE,ncid))
!  print *, 'finish open the file: ocean_velocity.res.nc.'
!  call check(NF90_INQ_VARID(ncid,'Time',varid))
!  call check(NF90_GET_VAR(ncid,varid,time))
!  print *, 'finish setting up Time.'
!  call check(NF90_INQ_VARID(ncid,'u',varid))
!  call check(NF90_GET_VAR(ncid,varid,u))
!  print *, 'finish setting up u.'
!  call check(NF90_INQ_VARID(ncid,'v',varid))
!  call check(NF90_GET_VAR(ncid,varid,v))
!  print *, 'finish setting up v.' 
!  call check(NF90_CLOSE(ncid)) 

!  call check(NF90_OPEN(sffile,NF90_NOWRITE,ncid))
!  print *, 'finish open the file: ocean_sbc.res.nc.'
!  call check(NF90_INQ_VARID(ncid,'Time',varid))
!  call check(NF90_GET_VAR(ncid,varid,time))
!  print *, 'finish setting up Time.'
!  call check(NF90_INQ_VARID(ncid,'t_surf',varid))
!  call check(NF90_GET_VAR(ncid,varid,sst))
!  print *, 'finish setting up t_surf.'
!  call check(NF90_INQ_VARID(ncid,'s_surf',varid))
!  call check(NF90_GET_VAR(ncid,varid,sss))
!  print *, 'finish setting up s_surf.' 
!  call check(NF90_CLOSE(ncid)) 
 

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
       DO m=1,5
        if ( i .eq. 8 ) then ! Here I set up lon=2.5
          !if (dodebug) print *, "lon=", lon(i)
          if ( j .eq. ilat(m)) then
            if (dodebug) print *, "lat=", lat(j)          
            wk(1)=id_t_obs
            wk(2)=lon(i)
            wk(3)=lat(j)
            wk(4)=lev(k)
            wk(5)=temp(i,j,k)
            wk(6)=0.5
            WRITE(91) wk
            wk(1)=id_s_obs
            wk(5)=salt(i,j,k)
            wk(6)=0.05 
            WRITE(91) wk
!        wk(1)=id_u_obs
!        wk(5)=u(i,j,k)
!        WRITE(91) wk
!        wk(1)=id_v_obs
!        wk(5)=v(i,j,k)
!        WRITE(91) wk
!        if (k .eq. 1) then
!          wk(1)=id_sst_obs
!          wk(5)=sst(i,j)
!          WRITE(91) wk
!          wk(1)=id_sss_obs
!          wk(5)=sss(i,j)
!          WRITE(91) wk          
!        endif
        !print *, id_s_obs,lon(i),lat(j),lev(k),salt(i,j,k),err
          end if
        end if 
       END DO  
      END DO
    END DO
  END DO

  CLOSE(91)
!  RETURN
!END SUBROUTINE read_grid
CONTAINS

SUBROUTINE check(status)
!===============================================================================
! Check the error status of the netcdf command
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
END SUBROUTINE check
  
END PROGRAM prof2letkf
