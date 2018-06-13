MODULE letkf_drifters_obs
!========================================================================

! obs_local_drifters
!========================================================================

  USE common
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE params_letkf
  USE params_obs
  USE params_drift_extra
  USE letkf_drifters_tools, ONLY: get_nobs_drifters
  USE vars_obs

  IMPLICIT NONE
  PUBLIC
  
  LOGICAL :: debug_hdxf_0 = .true.
  INTEGER :: cnt_obs_x, cnt_obs_y, cnt_obs_z, cnt_obs_t, cnt_obs_s
  INTEGER, DIMENSION(nv4d), SAVE :: cnt_obs = 0
  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_tout=0, cnt_sout=0, cnt_triout=0
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
  REAL(r_size),ALLOCATABLE :: tmpval(:)
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
  REAL(r_size),ALLOCATABLE :: tmp2val(:)

  INTEGER :: nobslots_dr(nslots)
  INTEGER :: islot,im,nn,l,ierr,n,i,j
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(21) :: obsfile_dr='obsTTNNN_drifters.dat'

  !STEVE: for obs qc:
  REAL(r_size) :: hdx2,mstd
  INTEGER :: gross_cnt,gross_2x_cnt
  LOGICAL :: obs_drif_exists

  WRITE(6,'(A)') 'Hello from set_letkf_obs_drifters'
 
  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dlat_zero = dist_zero / pi / re * 180.0d0

  
  if (.not. ALLOCATED(dlon_zero)) ALLOCATE(dlon_zero(nij1))
  do i=1,nij1
    dlon_zero(i) = dlat_zero / COS(pi*lat1(i)/180.0d0)
  enddo


  do islot=1,nslots
    im = myrank+1
    WRITE(obsfile_dr(4:8),'(I2.2,I3.3)') islot,im
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading an obs2-formatted file ',obsfile_dr

    INQUIRE( FILE=obsfile_dr, EXIST=obs_drif_exists )
    if (obs_drif_exists) then
      CALL get_nobs_drifters(obsfile_dr,10,nobslots_dr(islot))
    else
      WRITE(6,*) 'The drifters observation file does not exist'
      nobs_dr = 0
      DO_DRIFTERS = .false.
      RETURN
    end if
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nobslots_dr,nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  nobs_dr = SUM(nobslots_dr)
  WRITE(6,'(I10,A)') nobs_dr,' TOTAL DRIF OBSERVATIONS INPUT'

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
  ALLOCATE( tmptime(nobs_dr) )         !(DRIFTERS)
  ALLOCATE( tmpval(nobs_dr) )          !(DRIFTERS)
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
        & tmptime(nn+1:nn+nobslots_dr(islot)), tmpval(nn+1:nn+nobslots_dr(islot)) )  
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
  CALL MPI_BCAST( tmpval, nobs_dr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
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
  WRITE(6,*) "tmphdxf_dr(4,:) = ", tmphdxf_dr(4,:)
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
  WRITE(6,*) "tmpval(n) = ", tmpval(n)
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
  WRITE(6,*) "tmpval(n) = ", tmpval(n)
  WRITE(6,*) "STEVE: END DEBUGGING."
  WRITE(6,*)



  !STEVE: After processing ensemble members, apply some actions based on
  !forecast mean, to all observations
  cnt_obs_x = 0
  cnt_obs_y = 0
  cnt_obs_z = 0
  cnt_obs_t = 0
  cnt_obs_s = 0
  gross_cnt = 0
  gross_2x_cnt = 0
  
  WRITE(6,*) "Processing tmphdxf for n=1 to n=nobs_dr=",nobs_dr 
  WRITE(6,*) "and filtering bad observations..."

 !STEVE: this is the original version
  do n=1,nobs_dr
    tmpqc_dr(n) = MINVAL(tmpqc0_dr(n,:))
    if (tmpqc_dr(n) /= 1) CYCLE
    tmpdep_dr(n) = tmphdxf_dr(n,1) !note: tmpdep is just used as a dummy variable to compute the ensemble means
    do i=2,nbv
      tmpdep_dr(n) = tmpdep_dr(n) + tmphdxf_dr(n,i)
    enddo
    tmpdep_dr(n) = tmpdep_dr(n) / REAL(nbv,r_size)
    do i=1,nbv
      tmphdxf_dr(n,i) = tmphdxf_dr(n,i) - tmpdep_dr(n) ! Hdxf (perturbations from mean)
    enddo
    ! Now, tmpdep is defined appropriately as the obs departure from mean background

    tmpdep_dr(n) = tmpval(n) - tmpdep_dr(n)
 
    if (ABS(tmpdep_dr(n)) > gross_error*tmperr_dr(n)) then !gross error
      tmpqc_dr(n) = 0
      gross_cnt = gross_cnt + 1
    endif

    !STEVE: as a check, count the number of each type of observation
    !(DRIFTERS)
    if (INT(tmpelm_dr(n)) .eq. id_x_obs) cnt_obs_x = cnt_obs_x + 1
    if (INT(tmpelm_dr(n)) .eq. id_y_obs) cnt_obs_y = cnt_obs_y + 1
    if (INT(tmpelm_dr(n)) .eq. id_z_obs) cnt_obs_z = cnt_obs_z + 1
    if (nv4d > 3) then
      if (INT(tmpelm_dr(n)) .eq. id_t_dobs) cnt_obs_t = cnt_obs_t + 1
      if (INT(tmpelm_dr(n)) .eq. id_s_dobs) cnt_obs_s = cnt_obs_s + 1
    endif 

  enddo

  DEALLOCATE(tmpqc0_dr)

  WRITE(6,'(I10,A)') SUM(tmpqc_dr),' OBSERVATIONS TO BE ASSIMILATED'
  !STEVE:
  WRITE(6,*) "cnt_obs_x = ", cnt_obs_x
  WRITE(6,*) "cnt_obs_y = ", cnt_obs_y
  WRITE(6,*) "cnt_obs_z = ", cnt_obs_z
  if (nv4d > 3) then
    WRITE(6,*) "cnt_obs_t = ", cnt_obs_t
    WRITE(6,*) "cnt_obs_s = ", cnt_obs_s
  endif
  WRITE(6,*) "gross_cnt = ", gross_cnt
  WRITE(6,*) "gross_2x_cnt = ", gross_2x_cnt

  !cnt_obs(iv4d_x) = cnt_obs_x
  !cnt_obs(iv4d_y) = cnt_obs_y
  !cnt_obs(iv4d_z) = cnt_obs_z
  !if (nv4d > 3) then
  ! cnt_obs(iv4d_t) = cnt_obs_t
  !  cnt_obs(iv4d_s) = cnt_obs_s
  !end if

  !CALL monit_dep(nobs_dr,tmpelm_dr,tmpdep_dr,tmpqc_dr)

  nn = 0
  !STEVE: first, remove all of the Quality-Controlled data
  do n=1,nobs_dr
    if (tmpqc_dr(n) /= 1) CYCLE
    nn = nn+1
    tmpelm_dr(nn) = tmpelm_dr(n)
    tmplon_dr(nn) = tmplon_dr(n)
    tmplat_dr(nn) = tmplat_dr(n)
    tmplev_dr(nn) = tmplev_dr(n)
    tmpid(nn) = tmpid(n)
    tmperr_dr(nn) = tmperr_dr(n)
    tmpdep_dr(nn) = tmpdep_dr(n)
    tmphdxf_dr(nn,:) = tmphdxf_dr(n,:)
    tmpqc_dr(nn) = tmpqc_dr(n)
    tmptime(nn) = tmptime(n)
    tmpval(nn) = tmpval(n)
  enddo

  nobs_dr = nn
  WRITE(6,'(I10,A,I3.3)') nobs_dr,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank
!
! SORT
! 
  ALLOCATE( tmp2elm_dr(nobs_dr) )
  ALLOCATE( tmp2lon_dr(nobs_dr) )
  ALLOCATE( tmp2lat_dr(nobs_dr) )
  ALLOCATE( tmp2lev_dr(nobs_dr) )
  ALLOCATE( tmp2id(nobs_dr) )
  ALLOCATE( tmp2err_dr(nobs_dr) )
  ALLOCATE( tmp2dep_dr(nobs_dr) )
  ALLOCATE( tmp2hdxf_dr(nobs_dr,nbv) )
  ALLOCATE( tmp2time(nobs_dr) ) 
  ALLOCATE( tmp2val(nobs_dr) )

  ALLOCATE( obselm_dr(nobs_dr) )
  ALLOCATE( obslon_dr(nobs_dr) )
  ALLOCATE( obslat_dr(nobs_dr) )
  ALLOCATE( obslev_dr(nobs_dr) )
  ALLOCATE( obsid(nobs_dr) )
  ALLOCATE( obserr_dr(nobs_dr) )
  ALLOCATE( obsdep_dr(nobs_dr) )
  ALLOCATE( obshdxf_dr(nobs_dr,nbv) )
  ALLOCATE( obstime(nobs_dr) )
  ALLOCATE( obsval(nobs_dr) )
  ALLOCATE( nobsgrd_dr(nlon,nlat))

  nobsgrd_dr = 0
  nj = 0
  ! Count the number of observations within each latitude to ALLOCATABLE
  if (.true.) WRITE(6,*) "Count the number of observations within each latitude to ALLOCATABLE."
  do j=1,nlat-1
    do n=1,nobs_dr
      if (tmplat_dr(n) < lat(j) .OR. lat(j+1) <= tmplat_dr(n)) CYCLE
      nj(j) = nj(j) + 1 ! LUYU: record the number of observation within [lat(j), lat(j+1)).
    enddo
  enddo
  ! Record cumulative sum of observations up to this latitude
  ! Creates the basis for an indexing of observations from lat to lat
  if (.true.) WRITE(6,*) "Record cumulative sum of observations up to this latitude."
  do j=1,nlat-1
    njs(j) = SUM(nj(0:j-1)) ! njs record the number of observations within [0,lat(j))
  enddo

! Rearrange observation by latitude
  if (.true.) WRITE(6,*) "Rearrange observation by latitude."

  do j=1,nlat-1
    nn = 0
    do n=1,nobs_dr
      if (tmplat_dr(n) < lat(j) .OR. lat(j+1) <= tmplat_dr(n)) CYCLE
!     if (tmplon(n) >= lon(nlon)-EPSILON(1.0d0)) CYCLE   !STEVE: I added this to align with the same condition in the code above
                                                        !       Otherwise, sometimes nn /= nj(j)
      nn = nn + 1
      tmp2elm_dr(njs(j)+nn) = tmpelm_dr(n)
      tmp2lon_dr(njs(j)+nn) = tmplon_dr(n)
      tmp2lat_dr(njs(j)+nn) = tmplat_dr(n)
      tmp2lev_dr(njs(j)+nn) = tmplev_dr(n)
      tmp2id(njs(j)+nn) = tmpid(n) 
      tmp2err_dr(njs(j)+nn) = tmperr_dr(n)
      tmp2dep_dr(njs(j)+nn) = tmpdep_dr(n)
      tmp2hdxf_dr(njs(j)+nn,:) = tmphdxf_dr(n,:)    
      tmp2time(njs(j)+nn) = tmptime(n)
      tmp2val(njs(j)+nn) = tmpval(n) 
    enddo
  enddo
 
  if (.true.) WRITE(6,*) "tmp2lat_dr=", tmp2lat_dr
  if (.true.) WRITE(6,*) "tmp2lon_dr=", tmp2lon_dr
  
  ! For each latitude, identify the number of obs per longitude.
  ! Then, rearrange observations by longitude within each latitude step
  if (.true.) WRITE(6,*) "For each latitude, identify the number of obs per longitude."
  do j=1,nlat-1
    if (nj(j) == 0) then
      nobsgrd_dr(:,j) = njs(j)
      CYCLE
    endif
    nn = 0
    do i=1,nlon
      do n=njs(j)+1,njs(j)+nj(j)         

        ! Find the correct longitude bin for this observation...
        if (i < nlon) then
          if (tmp2lon_dr(n) < lon(i) .OR. lon(i+1) <= tmp2lon_dr(n)) CYCLE
        else
! STEVE: this is causing nn /= nj(j), the error thrown below.
!        We need these points that are skipped, otherwise there are
!        blank entries in the obselm etc. arrays, and this will
!        lead to problems during the main letkf algorithm.
!        Another solution may be to cut out all the empty entries
!        by changing the obsxxx indicies.
!
          if (tmp2lon_dr(n) < lon(nlon)) CYCLE

          !STEVE: debugging
          if (.false.) then
            WRITE(6,*) "n, nn, njs(j), nj(j) = ", n, nn, njs(j), nj(j)
            WRITE(6,*) "KEEPING, i == nlon == ", i, nlon
            WRITE(6,*) "tmp2lon_dr(n) = ", tmp2lon_dr(n)
            WRITE(6,*) "lon(nlon) = ", lon(nlon)
            !WRITE(6,*) "either tmp2lon(n) >= lon(nlon) .OR. 360.0d0 > tmp2lon(n)"
            WRITE(6,*) "tmp2lon(n) >= lon(nlon)"
            WRITE(6,*) "========================================================"
          ENDIF
        endif
        nn = nn + 1

        obselm_dr(njs(j)+nn) = tmp2elm_dr(n)
        obslon_dr(njs(j)+nn) = tmp2lon_dr(n)
        obslat_dr(njs(j)+nn) = tmp2lat_dr(n)
        obslev_dr(njs(j)+nn) = tmp2lev_dr(n)
        obsid(njs(j)+nn) = tmp2id(n)
        obserr_dr(njs(j)+nn) = tmp2err_dr(n)
        obsdep_dr(njs(j)+nn) = tmp2dep_dr(n)
        obshdxf_dr(njs(j)+nn,:) = tmp2hdxf_dr(n,:)
        obstime(njs(j)+nn) = tmp2time(n)
        obsval(njs(j)+nn) = tmp2val(n)
      enddo ! end do i
      
      ! This now contains the accumulated count of obs up to this lat, up to this lon
      nobsgrd_dr(i,j) = njs(j) + nn
    enddo ! end do j

    if (nn /= nj(j)) then
      WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(6,'(F6.2,A,F6.2)') lat(j),'<= LAT <',lat(j+1)
      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2lat_dr(njs(j)+1:njs(j)+nj(j))),'<= OBSLAT <',MAXVAL(tmp2lat_dr(njs(j)+1:njs(j)+nj(j)))
      WRITE(6,*) "j = ", j
      WRITE(6,*) "njs(j) = ", njs(j)
      WRITE(6,*) "nj(j) = ", nj(j)
      WRITE(6,*) "obs"
      !STEVE: this is bad, something is wrong
      WRITE(6,*) "STEVE: this error will cause matrix eigenvalue < 0 error."
      STOP 3
    endif

  enddo

  if (.true.) WRITE(6,*) "DEALLOCATE..."

  DEALLOCATE( tmpelm_dr )       !(DRIFTERS)
  DEALLOCATE( tmplon_dr )       !(DRIFTERS)
  DEALLOCATE( tmplat_dr )       !(DRIFTERS)
  DEALLOCATE( tmplev_dr )       !(DRIFTERS)
  DEALLOCATE( tmpid )           !(DRIFTERS)
  DEALLOCATE( tmperr_dr )       !(DRIFTERS)
  DEALLOCATE( tmpdep_dr )       !(DRIFTERS)
  DEALLOCATE( tmphdxf_dr )      !(DRIFTERS)
  !DEALLOCATE( tmpqc0_dr )       !(DRIFTERS)
  DEALLOCATE( tmpqc_dr )        !(DRIFTERS)
  DEALLOCATE( tmptime )         !(DRIFTERS)
  DEALLOCATE( tmpval )

  DEALLOCATE( tmp2elm_dr )
  DEALLOCATE( tmp2lon_dr )
  DEALLOCATE( tmp2lat_dr )
  DEALLOCATE( tmp2lev_dr )
  DEALLOCATE( tmp2id )
  DEALLOCATE( tmp2err_dr )
  DEALLOCATE( tmp2dep_dr )
  DEALLOCATE( tmp2hdxf_dr )
  DEALLOCATE( tmp2time )
  DEALLOCATE( tmp2val )

  if (.true.) WRITE(6,*) "Finish deallocating..."


END SUBROUTINE set_letkf_obs_drifters


SUBROUTINE read_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,oval)
!===============================================================================
! Read in observations with appended H(xb) for each ob
!===============================================================================
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: ohx(nn)
  REAL(r_size),INTENT(OUT) :: obhr(nn)
  REAL(r_size),INTENT(OUT) :: oval(nn)
  INTEGER,INTENT(OUT) :: oqc(nn)
  REAL(r_sngl) :: wk(10)
  INTEGER :: n,iunit
  
  !ALLOCATE(wk(obs2nrec))
  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  do n=1,nn
    READ(iunit) wk
    print *,wk
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    ohx(n)  = REAL(wk(7),r_size)
    oqc(n)  = NINT(wk(8))
    obhr(n) = REAL(wk(9),r_size)
    oval(n) = REAL(wk(10),r_size)
    !print *, elem(n),rlon(n),rlat(n),rlev(n),odat(n),oerr(n),ohx(n),oqc(n),obhr(n),oval(n)
  enddo

  CLOSE(iunit)

END SUBROUTINE read_obs2

END MODULE letkf_drifters_obs
