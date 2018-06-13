PROGRAM perturb_on_anal

 USE netcdf
 USE common
 USE params_letkf, ONLY: nbv

 IMPLICIT NONE

  CHARACTER (24) :: drifile="analpEEE.drifters_inp.nc"

  INTEGER :: i,j,k
  INTEGER :: nd,np,ncid,varid
   
  REAL(r_size), ALLOCATABLE :: pos(:,:),pos_new(:,:)
  REAL(r_size), ALLOCATABLE :: noise(:)

  DO i=1,nbv

    WRITE(drifile(6:8),'(I3.3)') i
    call check(NF90_OPEN(drifile,NF90_NOWRITE,ncid))    
    print *, 'Finish open the file: drifters_inp.nc.'
    
    ! Read Dimensions
    call check(NF90_INQ_DIMID(ncid,'nd',varid))
    call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nd))
    call check(NF90_INQ_DIMID(ncid,'np',varid))
    call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=np))

    ALLOCATE(pos(nd,np))
    ALLOCATE(pos_new(nd,np))
    ALLOCATE(noise(np))

    ! Read Data
    call check(NF90_INQ_VARID(ncid,'positions',varid))
    call check(NF90_GET_VAR(ncid,varid,pos))
    print *, 'Finish setting up positions.'

    call check(NF90_CLOSE(ncid))

    CALL com_randn(np,noise)
    pos_new(1,:)=pos(1,:) + 0.01*noise

    CALL com_randn(np,noise)
    pos_new(2,:)=pos(2,:) + 0.01*noise 

    pos_new(3,:)=pos(3,:)

    call check(NF90_OPEN(drifile,NF90_WRITE,ncid))    

    CALL check(NF90_INQ_VARID(ncid,'positions',varid))
    CALL check(NF90_PUT_VAR(ncid,varid,pos_new)) 

    call check(NF90_CLOSE(ncid))

    DEALLOCATE(pos,pos_new,noise)

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

END
