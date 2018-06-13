MODULE letkf_drifters

USE netcdf
USE common
USE common_mpi
USE common_letkf
USE letkf_obs, ONLY : nobs
USE params_letkf
USE params_model
USE params_obs
USE letkf_drifters_tools
UsE letkf_drifters_local

CONTAINS

SUBROUTINE das_drifters(gues4d,anal4d)
  ! das_letkf went through all model grid points. Now we're going to go through
  ! each drifter_id, as if we had appended them to the model state vector.
  
  USE common_letkf, ONLY: letkf_core
  USE params_model

  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues4d(nid1,num_times,nbv,nv4d) ! background ensemble
  REAL(r_size),INTENT(OUT) :: anal4d(nid1,num_times,nbv,nv4d) ! analysis ensemble
  REAL(r_size),ALLOCATABLE :: mean4d(:,:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv4d)
  LOGICAL :: ex
  INTEGER :: id,it,n,m,i,j,k,nobsl,ierr
  !STEVE: for debugging
  LOGICAL :: dodebug = .true.
  INTEGER :: nn
  !STEVE: added for drifters:
  INTEGER :: nobstotal,itim

  WRITE(6,'(A)') 'Hello from das_drifters'
  nobstotal = nobs_dr
  WRITE(6,'(A,I8)') 'Target observation numbers : nobs_dr=',nobs_dr!,', NTVS=',ntvs

  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    anal4d = gues4d
    RETURN
  ELSE                   !(OCEAN)
    anal4d = 0.0d0       !(OCEAN)
  END IF

  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean4d(nid1,num_times,nv4d))
  CALL ensmean_drifters(nbv,nid1,gues4d,mean4d)
  !STEVE: nid1 is each one of the grid points

  DO n=1,nv4d
    DO m=1,nbv
      DO k=1,num_times
        DO i=1,nid1
          gues4d(i,k,m,n) = gues4d(i,k,m,n) - mean4d(i,k,n)
        END DO
      END DO
    END DO
  END DO

  !
  ! multiplicative inflation
  !
  !STEVE: (using simple constant inflation for now.)

  !
  ! MAIN ASSIMILATION LOOP
  !
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE(hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  WRITE(6,*) "... done."

  DO it=1,num_times !STEVE: go through every possible coordinate of the grid in list form...
    if (dodebug) WRITE(6,*) "it = ", it ! LUYU: note "it" is time index NOT Time
    
    DO i=1,nid1
      id = i + myrank * nid1 
      if (dodebug) WRITE(6,*) "id= ", id ! LUYU:: not "id" is id index NOT ID   

    ! For each coordinate, x,y,and z:
    DO n=1,nv4d
      ! Find the observations around this point.
      
      if (dodebug) WRITE(6,*) "Start calling obs_local..."
      CALL obs_local(id,it,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)
      ! LUYU:  IN: id, it var_local(n,:),nobstotal
      ! LUYU: OUT: hdxf, rdiag, rloc, dep, nobsl

      parm = cov_infl_mul !STEVE: keeping it simple
      if (dodebug) WRITE(6,*) "Finished calling obs_local..."
      if (dodebug) WRITE(6,*) "nobstotal=",nobstotal
      if (dodebug) WRITE(6,*) "nobsl=",nobsl 
      if (dodebug) WRITE(6,*) "hdxf=",hdxf
      if (dodebug) WRITE(6,*) "rdiag=",rdiag
      if (dodebug) WRITE(6,*) "rloc=",rloc
      if (dodebug) WRITE(6,*) "dep=",dep     

      CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))
      if (dodebug) WRITE(6,*) "Finished calling letkf_core"
      if (dodebug) WRITE(6,*) "trans(1,1,n)=", trans(1,1,n)
      ! The ":" in place of "itim" implies that all times are affected equally by
      ! the observed drifter location.
      DO m=1,nbv
        anal4d(i,it,m,n) = mean4d(i,it,n)
        DO k=1,nbv
          anal4d(i,it,m,n) = anal4d(i,it,m,n) + gues4d(i,it,k,n) * trans(k,m,n)
        END DO
      END DO

    ENDDO

    ENDDO

 ENDDO

  DEALLOCATE(hdxf,rdiag,rloc,dep)

! If there are observations that weren't on the model grid, add them to the
! output list.

END SUBROUTINE das_drifters

END MODULE letkf_drifters
