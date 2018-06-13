PROGRAM setup_ini_perturbation

 USE netcdf
 USE common
 USE params_obs
 USE params_letkf, ONLY: nbv

 IMPLICIT NONE

  INTEGER :: i,j,k,m,dist,jlat
  INTEGER :: ncid,istat,varid,dimid
  
  INTEGER :: nnlon, nnlat, nnlev, dodebug=1
  REAL(r_size), ALLOCATABLE :: tlon(:)
  REAL(r_size), ALLOCATABLE :: tlat(:)
  REAL(r_size), ALLOCATABLE :: tlev(:)
  REAL(r_size), ALLOCATABLE :: salt(:,:,:)
  REAL(r_size), ALLOCATABLE :: temp(:,:,:)
  REAL(r_size), ALLOCATABLE :: u(:,:,:)
  REAL(r_size), ALLOCATABLE :: v(:,:,:)
  REAL(r_size), ALLOCATABLE :: taux(:,:)
  REAL(r_size), ALLOCATABLE :: tauy(:,:)
  REAL(r_size), ALLOCATABLE :: taux_new(:,:)
  REAL(r_size), ALLOCATABLE :: tauy_new(:,:)
  REAL(r_size), ALLOCATABLE :: ones(:)

  REAL(r_size) :: time, err_x(1), err_y(1)

  CHARACTER (len=*), PARAMETER :: taufile="tau.nc"
  CHARACTER (10) :: outfile="tau.EEE.nc"

  call check(NF90_OPEN(taufile,NF90_NOWRITE,ncid))
  print *, 'Finish open the file: tau.nc.'
  ! Read Dimensions
  call check(NF90_INQ_DIMID(ncid,'grid_x_C',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nnlon))
  call check(NF90_INQ_DIMID(ncid,'grid_y_C',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nnlat))
    
  ALLOCATE(tlon(nnlon))
  ALLOCATE(tlat(nnlat))
  ALLOCATE(taux(nnlon,nnlat))
  ALLOCATE(tauy(nnlon,nnlat))
  ALLOCATE(taux_new(nnlon,nnlat))
  ALLOCATE(tauy_new(nnlon,nnlat))
  ALLOCATE(ones(nnlon))
   
  ! Read Data
  call check(NF90_INQ_VARID(ncid,'taux',varid))
  call check(NF90_GET_VAR(ncid,varid,taux))
  print *, 'Finish setting up taux.'
  call check(NF90_INQ_VARID(ncid,'tauy',varid))
  call check(NF90_GET_VAR(ncid,varid,tauy))
  print *, 'Finish setting up tauy.'
  call check(NF90_CLOSE(ncid))

  DO i=1,nbv

    CALL com_randn(1,err_y)
    CALL com_randn(1,err_x)

    DO j=1,nnlat
      !jlat = MODULO(j+dist,nnlat) + 1
      jlat = j
      taux_new(:,jlat)=taux(:,j) + err_x * 0.1d0 * ones    
      DO k=1,nnlon
        tauy_new(k,j) = tauy(k,j) 
        
      END DO 
    END DO

    WRITE(outfile(5:7), '(I3.3)') i

    CALL check(NF90_OPEN(outfile,NF90_WRITE,ncid))
    print *, 'Just opened file ', outfile

    CALL check(NF90_INQ_VARID(ncid,'taux',varid))
    CALL check(NF90_PUT_VAR(ncid,varid,taux_new))
    CALL check(NF90_INQ_VARID(ncid,'tauy',varid))
    CALL check(NF90_PUT_VAR(ncid,varid,tauy_new))

    CALL check( NF90_CLOSE(ncid) )    
    
  END DO

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
  
END PROGRAM setup_ini_perturbation

