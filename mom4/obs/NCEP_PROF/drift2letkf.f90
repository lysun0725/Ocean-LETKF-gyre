PROGRAM drift2letkf

 !USE common_mom4 ! Use check subroutine
 USE netcdf
 USE common
 USE params_obs

 IMPLICIT NONE

!SUBROUTINE read_grid(infile,v3d,v2d)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid,dimid
  
  INTEGER :: nd,np,nf !nf is number of fields

  REAL(r_size), ALLOCATABLE :: pos(:,:)
  REAL(r_size), ALLOCATABLE :: ids(:), noise(:)

  REAL(r_size) :: err(1)
  REAL(r_sngl) :: wk(7)

  !CHARACTER (len=*), PARAMETER :: tsfile="ocean_temp_salt.res.nc"
  !CHARACTER (len=*), PARAMETER :: uvfile="ocean_velocity.res.nc"
  CHARACTER (len=*), PARAMETER :: drffile="drifters_inp.nc"
  CHARACTER (len=*), PARAMETER :: outfile="obsin_drifters.dat"

  call check(NF90_OPEN(drffile,NF90_NOWRITE,ncid))
  print *, 'finish open the file: drifters_inp.nc.'
  ! Read Dim
  call check(NF90_INQ_DIMID(ncid,'nd',varid)) !Number of dimension
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nd))
  call check(NF90_INQ_DIMID(ncid,'np',varid)) !Number of drifters
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=np))

  ALLOCATE(pos(nd,np))
  ALLOCATE(ids(np))
  ALLOCATE(noise(nd))

  ! Read Data

  call check(NF90_INQ_VARID(ncid,'positions',varid))
  call check(NF90_GET_VAR(ncid,varid,pos))
  print *, 'finish setting up positions.'
  call check(NF90_INQ_VARID(ncid,'ids',varid))
  call check(NF90_GET_VAR(ncid,varid,ids))
  print *, 'finish setting up ids.'
  call check(NF90_CLOSE(ncid))

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  DO i=1,np
     CALL com_randn(nd,noise)
     wk(1)=id_x_obs
     wk(2)=pos(1,i) + 0*noise(1)
     wk(3)=pos(2,i) + 0*noise(2)
     wk(4)=pos(3,i)
     wk(5)=ids(i)
     wk(6)=0.01
     wk(7)=24 !Modified later
     WRITE(91) wk
     print *, wk
     !call com_randn(1,err)
     wk(1)=id_y_obs 
     !wk(6)=ABS(err(1))
     WRITE(91) wk
     print *, wk
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
  
END PROGRAM drift2letkf
