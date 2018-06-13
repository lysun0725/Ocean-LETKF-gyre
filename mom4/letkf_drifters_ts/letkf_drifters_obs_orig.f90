MODULE letkf_drifters_obs
!========================================================================

! obs_local_drifters
!========================================================================

  USE common
  USE common_mpi
  USE common_mom4
  USE common_obs_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE params_letkf
  USE params_obs
  USE vars_obs

  IMPLICIT NONE
  PUBLIC

  INTEGER :: cnt_obs_x, cnt_obs_y, cnt_obs_z
  INTEGER, DIMENSION(nv4d), SAVE :: cnt_obs = 0
  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0

CONTAINS

SUBROUTINE set_letkf_obs_drifters

  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm_dr(:)
  REAL(r_size),ALLOCATABLE :: tmplon_dr(:)
  REAL(r_size),ALLOCATABLE :: tmplat_dr(:)
  REAL(r_size),ALLOCATABLE :: tmplev_dr(:)
  REAL(r_size),ALLOCATABLE :: tmpid(:)
  REAL(r_size),ALLOCATABLE :: tmperr_dr(:)
  REAL(r_size),ALLOCATABLE :: tmpdep_dr(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf_dr(:,:)
  REAL(r_size),ALLOCATABLE :: tmptime(:)
  INTEGER,ALLOCATABLE :: tmpqc0_dr(:,:)
  INTEGER,ALLOCATABLE :: tmpqc_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2id(:)
  REAL(r_size),ALLOCATABLE :: tmp2err_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep_dr(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf_dr(:,:)
  REAL(r_size),ALLOCATABLE :: tmp2time(:)
  INTEGER,ALLOCATABLE :: tmp2qc_dr(:)

  INTEGER :: nobslots_dr(nslots)
  INTEGER :: islot,im,nn,l,ierr,n,i
  CHARACTER(21) :: obsfile_dr='obsTTNNN_drifters.dat'

  !STEVE: for obs qc:
  REAL(r_size) :: hdx2,mstd
  INTEGER :: gross_cnt,gross_2x_cnt

  do islot=1,nslots
    im = myrank+1
    WRITE(obsfile_dr(4:8),'(I2.2,I3.3)') islot,im
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading an obs2-formatted file ',obsfile_dr
    CALL get_nobs(obsfile_dr,obs2nrec,nobslots_dr(islot))
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nobslots_dr,nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  nobs_dr = SUM(nobslots_dr)
  WRITE(6,'(I10,A)') nobs_dr,' TOTAL OBSERVATIONS INPUT'

  if (nobs_dr == 0) then
    WRITE(6,'(A)') 'No observation assimilated'
    RETURN
  endif



  !-----------------------------------------------------------------------------
  ! INITIALIZE GLOBAL VARIABLES
  !-----------------------------------------------------------------------------

  ALLOCATE( tmpelm_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmplon_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmplat_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmplev_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmpid(nobs_dr) )           !(DRIFTERS)
  ALLOCATE( tmperr_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmpdep_dr(nobs_dr) )       !(DRIFTERS)
  ALLOCATE( tmphdxf_dr(nobs_dr,nbv) )  !(DRIFTERS)
  ALLOCATE( tmpqc0_dr(nobs_dr,nbv) )   !(DRIFTERS)
  ALLOCATE( tmpqc_dr(nobs_dr))         !(DRIFTERS)
  ALLOCATE( tmptime(nobs_dr) )      !(DRIFTERS)
  tmpqc0_dr = 0
  tmphdxf_dr = 0.0d0
  tmperr_dr = 0.0d0

  !-----------------------------------------------------------------------
  ! LOOP of timeslots
  !-----------------------------------------------------------------------
  nn = 0
  timeslots0: do islot=1,nslots
    if (nobslots_dr(islot) == 0) CYCLE
    l=0
    do 
      im = myrank+1+nprocs * l
      if (im > nbv) EXIT
      WRITE(obsfile_dr(4:8),'(I2.2,I3.3)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile_dr
      CALL read_obs2(obsfile_dr,nobslots_dr(islot),&
        & tmpelm_dr(nn+1:nn+nobslots_dr(islot)),tmplon_dr(nn+1:nn+nobslots_dr(islot)),&
        & tmplat_dr(nn+1:nn+nobslots_dr(islot)),tmplev_dr(nn+1:nn+nobslots_dr(islot)),&
        & tmpid(nn+1:nn+nobslots_dr(islot)),tmperr_dr(nn+1:nn+nobslots_dr(islot)),&
        & tmphdxf_dr(nn+1:nn+nobslots_dr(islot),im),tmpqc0_dr(nn+1:nn+nobslots_dr(islot),im),&
        & tmptime(nn+1:nn+nobslots_dr(islot)) )  
      l = l+1
    enddo
    nn = nn + nobslots_dr(islot)
  enddo timeslots0

  WRITE(6,*) "Commencing collecting obs on all procs..."
  !STEVE: broadcast the 1d arrays from root onto all procs
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_BCAST's..."
  CALL MPI_BCAST( tmpelm_dr, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  !STEVE: just to be safe, calling MPI_BARRIER after each MPI_BCAST
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplon_dr, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplat_dr, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplev_dr, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmpid, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmperr_dr, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !STEVE: compile the tmphdxf array on all procs
  ALLOCATE(wk2d(nobs_dr,nbv))
  wk2d = tmphdxf_dr
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_ALLREDUCE..."
  CALL MPI_ALLREDUCE(wk2d,tmphdxf_dr,nobs_dr*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)   
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(wk2d)

  !STEVE: compile the tmpqc0 array on all procs
  ALLOCATE(iwk2d(nobs_dr,nbv))
  iwk2d = tmpqc0_dr
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_ALLREDUCE..."
  CALL MPI_ALLREDUCE(iwk2d,tmpqc0_dr,nobs_dr*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(iwk2d)
  WRITE(6,*) "Finished collecting obs on all procs."

  WRITE(6,*) "STEVE: DEBUGGING..."
  WRITE(6,'(I10,A,I3.3)') nobs_dr,' OBSERVATIONS, MYRANK = ',myrank
  WRITE(6,*) "tmphdxf_dr(1,:) = ", tmphdxf_dr(1,:)
  WRITE(6,*) "tmphdxf_dr(2,:) = ", tmphdxf_dr(2,:)
  WRITE(6,*) "tmphdxf_dr(3,:) = ", tmphdxf_dr(3,:)
  WRITE(6,*) "..."
  WRITE(6,*) "tmphdxf_dr(5,:) = ", tmphdxf_dr(5,:)
  WRITE(6,*)
  n=1
  WRITE(6,*) "For n=1,"
  WRITE(6,*) "MINVAL(tmpqc0_dr(n,:)) = ",MINVAL(tmpqc0_dr(n,:))
  WRITE(6,*) "tmpelm_dr(n) = ", tmpelm_dr(n) 
  WRITE(6,*) "tmplon_dr(n) = ", tmplon_dr(n) 
  WRITE(6,*) "tmplat_dr(n) = ", tmplat_dr(n) 
  WRITE(6,*) "tmplev_dr(n) = ", tmplev_dr(n) 
  WRITE(6,*) "tmpid(n) = ", tmpid(n) 
  WRITE(6,*) "tmperr_dr(n) = ", tmperr_dr(n) 
  WRITE(6,*) "tmptime(n) = ", tmptime(n)
  WRITE(6,*)
  n=nobs_dr
  WRITE(6,*) "For n=nobs_dr=",nobs_dr
  WRITE(6,*) "MINVAL(tmpqc0_dr(n,:)) = ",MINVAL(tmpqc0_dr(n,:))
  WRITE(6,*) "tmpelm_dr(n) = ", tmpelm_dr(n) 
  WRITE(6,*) "tmplon_dr(n) = ", tmplon_dr(n) 
  WRITE(6,*) "tmplat_dr(n) = ", tmplat_dr(n) 
  WRITE(6,*) "tmplev_dr(n) = ", tmplev_dr(n) 
  WRITE(6,*) "tmpid(n) = ", tmpid(n) 
  WRITE(6,*) "tmperr_dr(n) = ", tmperr_dr(n) 
  WRITE(6,*) "tmptime(n) = ", tmptime(n)
  WRITE(6,*) "STEVE: END DEBUGGING."
  WRITE(6,*)



  !STEVE: After processing ensemble members, apply some actions based on
  !forecast mean, to all observations
  cnt_obs_x = 0
  cnt_obs_y = 0
  cnt_obs_z = 0
  gross_cnt = 0
  gross_2x_cnt = 0
  
  WRITE(6,*) "Processing tmphdxf for n=1 to n=nobs_dr=",nobs_dr 
  WRITE(6,*) "and filtering bad observations..."

 !STEVE: this is the original version
  tmpqc0_dr=1
  do n=1,nobs_dr
    tmpqc_dr(n) = MINVAL(tmpqc0_dr(n,:))
    if (tmpqc_dr(n) /= 1) CYCLE
    tmpdep_dr(n) = tmphdxf_dr(n,1) !note: tmpdep is just used as a dummy variable to compute the mean over the next few lines
    do i=2,nbv
      tmpdep_dr(n) = tmpdep_dr(n) + tmphdxf_dr(n,i)
    enddo
    tmpdep_dr(n) = tmpdep_dr(n) / REAL(nbv,r_size)
    do i=1,nbv
      tmphdxf_dr(n,i) = tmphdxf_dr(n,i) - tmpdep_dr(n) ! Hdxf (perturbations from mean)
    enddo
    ! Now, tmpdep is defined appropriately as the obs departure from mean background
    ! LUYU: for drifters, we directly take lon, lat and lev as y instead of odat. So we use select function to achieve this.
    SELECT CASE(INT(tmpelm_dr(n)))
    CASE(id_x_obs)
       tmpdep_dr(n) = tmplon_dr(n) - tmpdep_dr(n) ! y-Hx
    CASE(id_y_obs)
       tmpdep_dr(n) = tmplat_dr(n) - tmpdep_dr(n)
    CASE(id_z_obs)
       tmpdep_dr(n) = tmplev_dr(n) - tmpdep_dr(n)
    END SELECT
 
    if (ABS(tmpdep_dr(n)) > gross_error*tmperr_dr(n)) then !gross error
      tmpqc_dr(n) = 0
      gross_cnt = gross_cnt + 1
    endif

    !STEVE: as a check, count the number of each type of observation
    !(DRIFTERS)
    if (INT(tmpelm_dr(n)) .eq. id_x_obs) cnt_obs_x = cnt_obs_x + 1
    if (INT(tmpelm_dr(n)) .eq. id_y_obs) cnt_obs_y = cnt_obs_y + 1
    if (INT(tmpelm_dr(n)) .eq. id_z_obs) cnt_obs_z = cnt_obs_z + 1

  enddo
  WRITE(6,*) "tmpdep_dr"
  WRITE(6,*) tmpdep_dr
  WRITE(6,*) "tmperr_dr"
  WRITE(6,*) tmperr_dr
  DEALLOCATE(tmpqc0_dr)

  WRITE(6,*) "###luyu_debug"
  WRITE(6,*) tmpqc_dr

  WRITE(6,'(I10,A)') SUM(tmpqc_dr),' OBSERVATIONS TO BE ASSIMILATED'
  !STEVE:
  WRITE(6,*) "cnt_obs_x = ", cnt_obs_x
  WRITE(6,*) "cnt_obs_y = ", cnt_obs_y
  WRITE(6,*) "cnt_obs_z = ", cnt_obs_z
  WRITE(6,*) "gross_cnt = ", gross_cnt
  WRITE(6,*) "gross_2x_cnt = ", gross_2x_cnt

  cnt_obs(iv4d_x) = cnt_obs_x
  cnt_obs(iv4d_y) = cnt_obs_y
  cnt_obs(iv4d_z) = cnt_obs_z

  !CALL monit_dep(nobs_dr,tmpelm_dr,tmpdep_dr,tmpqc_dr)

  ALLOCATE( tmp2elm_dr(nobs_dr) )
  ALLOCATE( tmp2lon_dr(nobs_dr) )
  ALLOCATE( tmp2lat_dr(nobs_dr) )
  ALLOCATE( tmp2lev_dr(nobs_dr) )
  ALLOCATE( tmp2id(nobs_dr) )
  ALLOCATE( tmp2err_dr(nobs_dr) )
  ALLOCATE( tmp2dep_dr(nobs_dr) )
  ALLOCATE( tmp2hdxf_dr(nobs_dr,nbv) )
  ALLOCATE( tmp2qc_dr(nobs_dr) )
  ALLOCATE( tmp2time(nobs_dr) ) 

  nn = 0
  !STEVE: first, remove all of the Quality-Controlled data
  do n=1,nobs_dr
    if (tmpqc_dr(n) /= 1) CYCLE
    nn = nn+1
    tmp2elm_dr(nn) = tmpelm_dr(n)
    tmp2lon_dr(nn) = tmplon_dr(n)
    tmp2lat_dr(nn) = tmplat_dr(n)
    tmp2lev_dr(nn) = tmplev_dr(n)
    tmp2id(nn) = tmpid(n)
    tmp2err_dr(nn) = tmperr_dr(n)
    tmp2dep_dr(nn) = tmpdep_dr(n)
    tmp2hdxf_dr(nn,:) = tmphdxf_dr(n,:)
    tmp2qc_dr(nn) = tmpqc_dr(n)
    tmp2time(nn) = tmptime(n)
  enddo

  nobs_dr = nn
  WRITE(6,'(I10,A,I3.3)') nobs_dr,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank

!
! 
!
  DEALLOCATE( obsid )
  DEALLOCATE( obstime )
  ALLOCATE( obselm_dr(nobs_dr) )
  ALLOCATE( obslon_dr(nobs_dr) )
  ALLOCATE( obslat_dr(nobs_dr) )
  ALLOCATE( obslev_dr(nobs_dr) )
  ALLOCATE( obsid(nobs_dr) )
  ALLOCATE( obserr_dr(nobs_dr) )
  ALLOCATE( obsdep_dr(nobs_dr) )
  ALLOCATE( obshdxf_dr(nobs_dr,nbv) )
  ALLOCATE( obstime(nobs_dr) )

! REARRANGEMETN: LATER....

  obselm_dr = tmp2elm_dr(1:nobs_dr)
  obslon_dr = tmp2lon_dr(1:nobs_dr)
  obslat_dr = tmp2lat_dr(1:nobs_dr)
  obslev_dr = tmp2lev_dr(1:nobs_dr)
     obsid  = tmp2id(1:nobs_dr)
  obserr_dr = tmp2err_dr(1:nobs_dr)
  obsdep_dr = tmp2dep_dr(1:nobs_dr)
  do i=1,nbv
    obshdxf_dr(:,i) = tmp2hdxf_dr(1:nobs_dr,i)
  enddo
  obstime  = tmp2time(1:nobs_dr)

END SUBROUTINE set_letkf_obs_drifters

END MODULE letkf_drifters_obs
