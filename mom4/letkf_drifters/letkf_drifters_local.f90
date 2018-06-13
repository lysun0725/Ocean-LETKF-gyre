MODULE letkf_drifters_local
  !STEVE: the letkf_local had a few details that were specific to grid-based obs,
  !       so I copied it here and made the requried changes for the drifter Lagrangian obs.
  !
  USE common
  USE common_mpi
  USE common_mom4
  USE common_obs_mom4
  USE common_mpi_mom4
  USE common_letkf
  !USE letkf_drifters_obs, ONLY: dist_zero, dist_zerov, dlon_zero, dlat_zero !contains debug_hdxf_0, and nobsgrd
  USE params_letkf, ONLY: nbv, DO_NO_VERT_LOC, localization_method
  USE params_obs
  USE vars_obs
  USE letkf_drifters_obs
  USE letkf_drifters_tools

! Designed by Dr. Stephen G. Penny
! University of Maryland, College Park
! This module is an offshoot of letkf_tools, which previously contained
! all localization routines. The purpose of isolating these routines in
! an independent module is to allow for further development of localization
! approaches, in particular those utilizing computational geometry tools
! such as:
! kd-tree, kNN search and range search
! Graph representation and A* search

  !STEVE: Testing "Vertical Tube" localization:
  !       i.e. the localization is not applied vertically
  ! This provides the benefit that 
  ! (1) the analysis only has to be computed once
  ! per horizontal gridpoint, thus providing a nlevX (40X) speedup
  ! (2) the altimetry, SST, SSH, and bottom pressure (GRACE) can now be applied
  ! as direct constraints on the water column.
  !
  ! There is precedence for this as in the paper "Reconstructing the Ocean's
  ! Interior from Surface Data" Wang et al. (2013)
  !

  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d+nv4d,nid_dobs) = 1.0d0
  !NOTE: the obs are the rows and the model variables the columns
  !LUYU: if we are going to add profile observation, we need to change 2nd argument to nid_dobs+nid_obs

  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d+nv4d)

CONTAINS

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
!SUBROUTINE obs_local(ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl)
SUBROUTINE obs_local(ij,ilev,var_local,hdxf,rdiag,rloc,dep,nobsl,nobsl_dr,nobstotal)
  ! LUYU: obs_local: for each grid point in the MODEL SPACE, we select local observations of profile and drifters. we define the number of profiles as nn, while drifters as nn_dr.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nobstotal !LUYU: "ilev" = index of time 
  REAL(r_size),INTENT(IN) :: var_local(nid_dobs)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl, nobsl_dr
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:),nobs_use_dr(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,nn_dr,tvnn,iobs
  !STEVE: for (OCEAN):
  LOGICAL :: blocked_by_land     !Use these three lines to identify gulf vs. pacific points in localization
  REAL(r_size) :: xlat,xlon,lxpa
  REAL(r_size) :: olat,olon,olev,oelm,lopa
  REAL(r_size) :: f1lon,f1lat,f2lon,f2lat
  REAL(r_size) :: dist1,dist2,a,b,u,v,lu,lv
  REAL(r_size) :: ecc,ecc2,fcc,fcc2,theta,theta_rotate,sigma_min
  REAL(r_size) :: dist_zero_a,dist_zero_b,dist_min
  REAL(r_size) :: sigma_a,sigma_b
  REAL(r_size) :: minr,maxr,xdis,fcc0,fcc1,maxdN,maxdS,xrad
  REAL(r_size), PARAMETER :: cmpersec2kmperday=0.864d0, days = 5.0d0
  !STEVE: for initializing splits:
  LOGICAL :: splitlon = .false. ! initialize
  LOGICAL :: splitlonL = .false., splitlonR = .false. ! initialize
  LOGICAL :: splitlat = .false. ! initialize
  LOGICAL :: fulllon = .false. ! initialize
  !STEVE:
  LOGICAL :: dodebug = .false.
  !debug TEST:
  INTEGER :: nn_old,imin_old,imax_old,jmin_old,jmax_old
  REAL(r_size) :: minlon_old,maxlon_old,minlat_old,maxlat_old

!
! INITIALIZE
!
  IF( DO_OCEAN ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF

  IF( DO_DRIFTERS ) THEN
    ALLOCATE(nobs_use_dr(nobs_dr))
  END IF
  nn = 0
  nn_dr = 0

! Set model coordinates
  xlat = lat1(ij)
  xlon = lon1(ij)
  if (dodebug) WRITE(6,*) "---------------------------------- ij = ", ij
  if (dodebug) WRITE(6,*) "---------------------------------- ilev = ", ilev
  
!
! data search (NOTE: this only wraps in the longitude direction)
! Use bounding box first to reduce obs considered.
!
  minlon = xlon - dlon_zero(ij)
  maxlon = xlon + dlon_zero(ij)
  minlat = xlat - dlat_zero !STEVE: this should be changed to be (ij) dependent
  maxlat = xlat + dlat_zero

  if (dodebug) then 
    WRITE(6,*) "---minlon = ", minlon
    WRITE(6,*) "---maxlon = ", maxlon
    WRITE(6,*) "---minlat = ", minlat
    WRITE(6,*) "---maxlat = ", maxlat
    
    WRITE(6,*) "+++minlon = ", minlon_old
    WRITE(6,*) "+++maxlon = ", maxlon_old
    WRITE(6,*) "+++minlat = ", minlat_old
    WRITE(6,*) "+++maxlat = ", maxlat_old

    WRITE(6,*) "nn_old   = ", nn_old
    WRITE(6,*) "imin_old = ", imin_old
    WRITE(6,*) "imax_old = ", imax_old
    WRITE(6,*) "jmin_old = ", jmin_old
    WRITE(6,*) "jmax_old = ", jmax_old
  endif

  !STEVE: rewrite
  fulllon = .false.
  splitlonL = .false.
  splitlonR = .false.
  splitlat = .false.
  if (maxlon - minlon > lonf - lon0) then 
    fulllon = .true.
    if (dodebug) then 
      WRITE(6,*) "FullLon."
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "xlat = ", xlat
    endif
    minlon = lon0
    maxlon = lonf
  else
    if (minlon < lon0) then
      ! Local region overlaps on left
      if (dodebug) WRITE(6,*) "split minlon = ", minlon
      splitlonL = .true.
      minlon = REAL(lonf - abs(lon0 - minlon) + wrapgap,r_size)
    endif
    if (maxlon > lonf) then
      ! Local region overlaps on right
      if (dodebug) WRITE(6,*) "split maxlon = ", maxlon
      if (splitlonL) WRITE(6,*) "WARNING: lon already split. ij = ", ij
      splitlonR = .true.
      maxlon = REAL(lon0 + abs(maxlon - lonf) - wrapgap,r_size)
    endif 
    if (minlat < lat0) then
      ! Local region overlaps on bottom
      splitlat= .true.
      minlat = lat0
    endif
    if (maxlat > latf) then
      ! Local region overlaps on top
      if (splitlat) WRITE(6,*) "WARNING: lat already split. ij = ", ij
      splitlat = .true.
      maxlat = latf
    endif 
  endif

  if (splitlat) then
    ! Too close to pole, need to work out the details of this
    ! For now, cut off the local region at the pole 
!   WRITE(6,*) "letkf_local.f90::obs_local: Grid point is close to pole."
  endif
! else
  DO jmin=1,nlat-2
    IF(minlat < lat(jmin+1)) EXIT
  END DO
  DO jmax=1,nlat-2
    IF(maxlat < lat(jmax+1)) EXIT
  END DO
  nn = 1
  nn_dr = 1
! endif

  if ( (splitlonL .or. splitlonR) .and. .not. fulllon) then
    if (dodebug) WRITE(6,*) "splitlonL = ", splitlonL
    if (dodebug) WRITE(6,*) "splitlonR = ", splitlonR

    ! Do left
    imin = 1
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    if (dodebug) then
      WRITE(6,*) "splitlon left: nn = ", nn, "nn_dr = ", nn_dr
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif

    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)    
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)

    ! Do right
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    imax = nlon
    if (dodebug) then
      WRITE(6,*) "splitlon right: nn = ", nn, "nn_dr = ", nn_dr
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif

    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)

  elseif ( .not. fulllon ) then
    ! The local range is entirely contained within the global domain
    ! Find the min and max latitudes
    DO jmin=1,nlat-2
      IF(minlat < lat(jmin+1)) EXIT
    END DO
    DO jmax=jmin,nlat-2
      IF(maxlat < lat(jmax+1)) EXIT
    END DO
    ! Find the min and max longitudes
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    DO imax=imin,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO

    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)
  else
    ! The domain is the entire circle around the globe (zonally)
    ! Find the min and max latitudes
    DO jmin=1,nlat-2
      IF(minlat < lat(jmin+1)) EXIT
    END DO
    DO jmax=jmin,nlat-2
      IF(maxlat < lat(jmax+1)) EXIT
    END DO

    if (dodebug) then
      WRITE(6,*) "Doing FullLon..."
    endif
    imin=1
    imax=nlon

    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)
  endif
  nn = nn - 1
  nn_dr = nn_dr-1
  !STEVE: end rewrite

  ! If no observations remain, then we're done.
  ! LUYU: REMINDER: nobsl means # of profile obs around (ij,ilev); nobsl_dr means # of drifters obs around (ij,ilev).
  IF(nn < 1) THEN
    nobsl = 0  
    IF(nn_dr < 1) THEN
      nobsl_dr = 0
      RETURN
    END IF
  END IF

!
! CONVENTIONAL
!
!STEVE: most of the localization section has been completely edited for (OCEAN)
!LUYU: (**TBD**) This is a conventional routine... Need to consider later whether to construct a new subroutine for the localization!!
 
 nobsl = 0
 IF ( DO_OCEAN ) THEN ! profile observation exists
  DO n=1,nn
    !STEVE: (future) use custom localization with CGAL/BOOST algorithms
    !
    ! Observational localization Distance Cutoff
    !
    olat = obslat(nobs_use(n))
    olon = obslon(nobs_use(n))
    olev = obslev(nobs_use(n))
    oelm = obselm(nobs_use(n))

    ! variable localization
    !
    SELECT CASE(NINT(oelm))
    CASE(id_u_obs)
        iobs=1
    CASE(id_v_obs)
        iobs=2
    CASE(id_t_obs)
        iobs=3
    CASE(id_s_obs)   !(OCEAN)
        iobs=4
    CASE(id_ssh_obs) !(OCEAN)
        iobs=5
    CASE(id_sst_obs) !(OCEAN)
        iobs=6
    CASE(id_sss_obs) !(OCEAN)
        iobs=7
    CASE(id_eta_obs) !(OCEAN)
        iobs=8
    CASE DEFAULT
        WRITE(6,*) "letkf_local.f90 :: there is no variable localization for obs-type :: ", oelm
        WRITE(6,*) "letkf_local.f90 :: exitING..."
        STOP(95)
    END SELECT
    IF(var_local(iobs) < TINY(var_local)) CYCLE 
    !STEVE: skip obs that are set to "0" impact on this model variable.

    !--------------------------------------------------------------------
    ! vertical localization
    !--------------------------------------------------------------------
    !STEVE: make vertical localization depth dependent
    !STEVE: Could alternatively treat surface obs as if they occur at depth 0, regardless of model's top level
    IF(NINT(oelm) == id_ssh_obs) THEN
        !dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))        !(OCEAN)
        !STEVE: no need for log scale, but, don't have ssh obs right now... so I'll test this later
        !dlev = ABS( olev - lev(ilev) )
        dlev = 0.0d0 !STEVE: allow all levels to be influenced by ssh
    ELSE IF(NINT(oelm) == id_sst_obs) THEN
        dlev = ABS( olev - lev(ilev) )
        !WRITE(6,*) "obs_local:: SST ob, dlev = ",dlev
    ELSE IF(NINT(oelm) == id_sss_obs) THEN
        dlev = ABS( olev - lev(ilev) )
    ELSE IF(NINT(oelm) == id_u_obs .or. &
            NINT(oelm) == id_v_obs .or. &
            NINT(oelm) == id_t_obs .or. &
            NINT(oelm) == id_s_obs   ) THEN
        dlev = ABS( olev - lev(ilev) )
    ELSE
        dlev = 0.0d0
    END IF
    IF(dlev > dist_zerov) CYCLE

    !
    ! horizontal localization
    !
    horizontal_localization : if (localization_method .eq. 1) then
        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
        sigma_b = sigma_a
        dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    elseif (localization_method .eq. 2) then
        !STEVE: THIS LOCALIZATION METHOD WAS AN ERROR, however the results were
        !good, and better than a number of different horizontal localization
        !schemes that I tried. So I am keeping it here for further study.
        !
        ! The input sigma_obs should be 1000km, and simga_obs0 should be 200km

        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 !100.0 * 1000.0 !make the minimum distance 100 kilometers at 90ºN or 90ºS
                                  !WARNING: this should be bigger than the minimum grid cell width
!       dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad
        sigma_b = sigma_a
        dist_zero_a = sigma_a * (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_b = sigma_a
!       dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    else !STEVE: use the original (default) localization
        CALL com_distll_1(olon,olat,xlon,xlat,dist)
        IF(dist > dist_zero ) CYCLE
    endif horizontal_localization

    ! STEVE: ADD CHECK FOR OBSERVATIONS THAT ARE OCCLUDED BY LAND!
    ! STOLE atlpac from SODA, in the future, implement a general method
    ! for occluding points that can't be reached by a SPM in the local range
    blocked_by_land = .false.
    CALL atlpac(xlat,xlon,lxpa)
    CALL atlpac(olat,olon,lopa)
    IF (  lxpa .eq. 1 .and. lopa .eq. 2 &        !in caribbean, pacific and atlantic
       .or. lxpa .eq. 2 .and. lopa .eq. 1 &        !in caribbean, atlantic and pacific
         ) blocked_by_land = .true.
    IF (blocked_by_land) CYCLE                   !STEVE: the two points are in atl and pac

!--------- End of Observation Culling ---------!

    nobsl = nobsl + 1
    hdxf(nobsl,:) = obshdxf_dr(nobs_use(n),:)
    dep(nobsl) = obsdep_dr(nobs_use(n))
    tmperr=obserr_dr(nobs_use(n))
    rdiag(nobsl) = tmperr * tmperr
    if (ALLOCATED(obs_useidx)) then
      obs_useidx(nobsl) = nobs_use(n)
    endif
    !
    ! Observational localization (weighting)
    !
    ! Note: var_local scales the localization weighting based on the parameter
    ! set in letkf_tools.f90. A row corresponding to the model parameter is
    ! input to this subroutine, and the column indicates the proportion of that
    ! type of observation to use.
    observation_localization : if (localization_method .eq. 1 ) then
                                                       
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
        ELSE
          rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF

    elseif (localization_method .eq. 2) then
      
        !STEVE: trying the original approach that seemed to work before:
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_obs)**2) * var_local(iobs)
        ELSE
          rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF
                                                       
    elseif (localization_method .eq. 0) then
        !STEVE: this is R^2 and the R localization gaussian function
        rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
    else
        print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
        STOP(3)
    endif observation_localization
  END DO
 END IF ! test nobs > 0


 !###############################################################################
 !LUYU: start searching DRIFTERS OBS 
 nobsl_dr = 0
 IF (DO_DRIFTERS) THEN 
  DO n=1,nn_dr
    olat = obslat_dr(nobs_use_dr(n))
    olon = obslon_dr(nobs_use_dr(n))
    olev = obslev_dr(nobs_use_dr(n))
    oelm = obselm_dr(nobs_use_dr(n))
    
    SELECT CASE(NINT(oelm))    
    CASE(id_x_obs)
        iobs=1
    CASE(id_y_obs)
        iobs=2
    CASE(id_z_obs)
        iobs=3
    CASE(id_t_obs)   !(OCEAN)
        iobs=4
    CASE(id_s_obs) !(OCEAN)
        iobs=5
    END SELECT
    IF(var_local(iobs) < TINY(var_local)) CYCLE 

    ! vertical localization
    !
    !STEVE: make vertical localization depth dependent
    !STEVE: Could alternatively treat surface obs as if they occur at depth 0, regardless of model's top level
    IF(NINT(oelm) == id_t_obs .or. &
            NINT(oelm) == id_s_obs   ) THEN
        dlev = ABS( olev - lev(ilev) )
    ELSE
        dlev = 0.0d0
    END IF
    IF(dlev > dist_zerov) CYCLE

    !
    ! horizontal localization
    !
    horizontal_localization_dr : if (localization_method .eq. 1) then
        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
        sigma_b = sigma_a
        dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.true.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    elseif (localization_method .eq. 2) then
        !LUYU: (***TBD***) Need to modify this localization on drifters later!!

        dist_min = sigma_obs0 !100.0 * 1000.0 !make the minimum distance 100 kilometers at 90ºN or 90ºS
                                  !WARNING: this should be bigger than the minimum grid cell width
!       dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad
        sigma_b = sigma_a
        dist_zero_a = sigma_a * (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_b = sigma_a
!       dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    else !STEVE: use the original (default) localization
        CALL com_distll_1(olon,olat,xlon,xlat,dist)
        IF(dist > dist_zero ) CYCLE
    endif horizontal_localization_dr

    ! STEVE: ADD CHECK FOR OBSERVATIONS THAT ARE OCCLUDED BY LAND!
    ! STOLE atlpac from SODA, in the future, implement a general method
    ! for occluding points that can't be reached by a SPM in the local range
    blocked_by_land = .false.
    CALL atlpac(xlat,xlon,lxpa)
    CALL atlpac(olat,olon,lopa)
    IF (  lxpa .eq. 1 .and. lopa .eq. 2 &        !in caribbean, pacific and atlantic
       .or. lxpa .eq. 2 .and. lopa .eq. 1 &        !in caribbean, atlantic and pacific
         ) blocked_by_land = .true.
    IF (blocked_by_land) CYCLE                   !STEVE: the two points are in atl and pac

!--------- End of Observation Culling ---------!

    nobsl_dr = nobsl_dr + 1
    hdxf(nobsl+nobsl_dr,:) = obshdxf_dr(nobs_use_dr(n),:)
    dep(nobsl+nobsl_dr)    = obsdep_dr(nobs_use_dr(n))
    tmperr=obserr_dr(nobs_use_dr(n))
    rdiag(nobsl+nobsl_dr) = tmperr * tmperr
    if (ALLOCATED(obs_useidx)) then
      obs_useidx(nobsl+nobsl_dr) = nobs_use_dr(n)
    endif
    !
    ! Observational localization (weighting)
    !
    ! Note: var_local scales the localization weighting based on the parameter
    ! set in letkf_tools.f90. A row corresponding to the model parameter is
    ! input to this subroutine, and the column indicates the proportion of that
    ! type of observation to use.
    observation_localization_dr : if (localization_method .eq. 1 ) then
!       rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a  )**2 + (dlev/sigma_obsv)**2)) &
!                                                    & * var_local(iobs)
!STEVE: testing different localization functions:
!STEVE: doubing standard deviation distance
!       rloc(nobsl) =EXP( -0.5d0 * ( ( dist/(2.0d0*sigma_a) )**2 + ( dlev/(2.0d0*sigma_obsv) )**2 ) ) &
!                                                    & * var_local(iobs)
!STEVE: removing all localization other than culling as applied above based on
!absolute distance:
!       rloc(nobsl) = 1.0d0 * var_local(iobs)
                                                       
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
        ELSE
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF

    elseif (localization_method .eq. 2) then
      
        !STEVE: trying the original approach that seemed to work before:
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * (dist/sigma_obs)**2) * var_local(iobs)
        ELSE
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF
                                                       
    elseif (localization_method .eq. 0) then
        !STEVE: this is R^2 and the R localization gaussian function
        rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
    else
        print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
        STOP(3)
    endif observation_localization_dr
  END DO
 END IF ! test nobs_dr > 0 

  !STEVE: this should never happen, if it does something went wrong
  IF( nobsl_dr + nobsl> nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL_DR + NOBSL =',nobsl_dr + nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN_DR,NN,TVNN=', ij, nn_dr,nn, tvnn
    STOP 99
  END IF
 
  IF( nobs_dr > 0 ) THEN
    DEALLOCATE(nobs_use_dr)
  END IF

  RETURN
END SUBROUTINE obs_local

SUBROUTINE obs_local_drifters(id,it,var_local,v4d,hdxf,rdiag,rloc,dep,nobsl,nobsl_dr,nobstotal)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: id,it,nobstotal !,nvar
  REAL(r_size),INTENT(IN) :: var_local(nid_dobs)
  REAL(r_size),INTENT(IN) :: v4d(nv4d)  ! record ensemble mean positions of drifter "id" at time index "it"
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl, nobsl_dr
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:),nobs_use_dr(:),nobs_use_tmp(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: i,n,nn,nn_dr,tvnn,iobs,m,mm,isfirst
  !STEVE: for (OCEAN):
  LOGICAL :: blocked_by_land     !Use these three lines to identify gulf vs. pacific points in localization
  REAL(r_size) :: xlat,xlon,lxpa,ylon,ylat,xlev,ylev
  REAL(r_size) :: olat,olon,olev,oelm,lopa
  REAL(r_size) :: f1lon,f1lat,f2lon,f2lat
  REAL(r_size) :: dist1,dist2,a,b,u,v,lu,lv,dist_d
  REAL(r_size) :: ecc,ecc2,fcc,fcc2,theta,theta_rotate,sigma_min
  REAL(r_size) :: dist_zero_a,dist_zero_b,dist_min,dlon_zero_dr
  REAL(r_size) :: sigma_a,sigma_b
  REAL(r_size) :: minr,maxr,xdis,fcc0,fcc1,maxdN,maxdS,xrad
  REAL(r_size), PARAMETER :: cmpersec2kmperday=0.864d0, days = 5.0d0
  !STEVE: for initializing splits:
  LOGICAL :: splitlon = .false. ! initialize
  LOGICAL :: splitlonL = .false., splitlonR = .false. ! initialize
  LOGICAL :: splitlat = .false. ! initialize
  LOGICAL :: fulllon = .false. ! initialize
  !STEVE:
  LOGICAL :: dodebug = .true.
  !debug TEST:
  INTEGER :: nn_old,imin_old,imax_old,jmin_old,jmax_old,itime,ind
  REAL(r_size) :: minlon_old,maxlon_old,minlat_old,maxlat_old
  

  IF( DO_OCEAN ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF

  IF( DO_DRIFTERS ) THEN
    ALLOCATE(nobs_use_dr(nobs_dr))
  END IF
  nn = 0
  nn_dr = 0

! Set model coordinates
  xlon = v4d(1)
  xlat = v4d(2)
  xlev = v4d(3)
  if (dodebug) WRITE(6,*) "---------------------------------- id = ", id
  if (dodebug) WRITE(6,*) "---------------------------------- it = ", it
  if (dodebug) WRITE(6,*) "xlon = ", xlon
  if (dodebug) WRITE(6,*) "xlat = ", xlat
  if (dodebug) WRITE(6,*) "xlev = ", xlev
!
! data search (NOTE: this only wraps in the longitude direction)
! Use bounding box first to reduce obs considered.
!
! LUYU: Redefine dlon_zero_dr for drifters
  dlon_zero_dr = dlat_zero / COS(pi*xlat/180.0d0)

  minlon = xlon - dlon_zero_dr
  maxlon = xlon + dlon_zero_dr
  minlat = xlat - dlat_zero !STEVE: this should be changed to be (ij) dependent
  maxlat = xlat + dlat_zero

  if (dodebug) then 
    WRITE(6,*) "---minlon = ", minlon
    WRITE(6,*) "---maxlon = ", maxlon
    WRITE(6,*) "---minlat = ", minlat
    WRITE(6,*) "---maxlat = ", maxlat
    
    WRITE(6,*) "---lonf = ", lonf
    WRITE(6,*) "---lon0 = ", lon0
    WRITE(6,*) "---latf = ", latf
    WRITE(6,*) "---lon0 = ", lat0
  endif

  !STEVE: rewrite
  fulllon = .false.
  splitlonL = .false.
  splitlonR = .false.
  splitlat = .false.
  if (maxlon - minlon > lonf - lon0) then 
    fulllon = .true.
    if (dodebug) then 
      WRITE(6,*) "FullLon."
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "xlat = ", xlat
    endif
    minlon = lon0
    maxlon = lonf
  else
    if (minlon < lon0) then
      ! Local region overlaps on left
      if (dodebug) WRITE(6,*) "split minlon = ", minlon
      splitlonL = .true.
      minlon = REAL(lonf - abs(lon0 - minlon) + wrapgap,r_size)
    endif
    if (maxlon > lonf) then
      ! Local region overlaps on right
      if (dodebug) WRITE(6,*) "split maxlon = ", maxlon
      if (splitlonL) WRITE(6,*) "WARNING: lon already split. id = ", id
      splitlonR = .true.
      maxlon = REAL(lon0 + abs(maxlon - lonf) - wrapgap,r_size)
    endif 
    if (minlat < lat0) then
      ! Local region overlaps on bottom
      splitlat= .true.
      minlat = lat0
    endif
    if (maxlat > latf) then
      ! Local region overlaps on top
      if (splitlat) WRITE(6,*) "WARNING: lat already split. id = ", id
      splitlat = .true.
      maxlat = latf
    endif 
  endif

  if (splitlat) then
    ! Too close to pole, need to work out the details of this
    ! For now, cut off the local region at the pole 
!   WRITE(6,*) "letkf_local.f90::obs_local: Grid point is close to pole."
  endif
! else
  DO jmin=1,nlat-2
    IF(minlat < lat(jmin+1)) EXIT
  END DO
  DO jmax=1,nlat-2
    IF(maxlat < lat(jmax+1)) EXIT
  END DO
  nn_dr = 1
! endif

  if ( (splitlonL .or. splitlonR) .and. .not. fulllon) then
    if (dodebug) WRITE(6,*) "splitlonL = ", splitlonL
    if (dodebug) WRITE(6,*) "splitlonR = ", splitlonR

    ! Do left
    imin = 1
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    if (dodebug) then
      WRITE(6,*) "splitlon left: nn_dr = ", nn_dr
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif
    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)

    ! Do right
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    imax = nlon
    if (dodebug) then
      WRITE(6,*) "splitlon right: nn_dr = ", nn_dr
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif
    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)

  elseif ( .not. fulllon ) then
    ! The local range is entirely contained within the global domain
    ! Find the min and max latitudes
    DO jmin=1,nlat-2
      IF(minlat < lat(jmin+1)) EXIT
    END DO
    DO jmax=jmin,nlat-2
      IF(maxlat < lat(jmax+1)) EXIT
    END DO
    ! Find the min and max longitudes
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    DO imax=imin,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)
  else
    ! The domain is the entire circle around the globe (zonally)
    ! Find the min and max latitudes
    DO jmin=1,nlat-2
      IF(minlat < lat(jmin+1)) EXIT
    END DO
    DO jmax=jmin,nlat-2
      IF(maxlat < lat(jmax+1)) EXIT
    END DO

    if (dodebug) then
      WRITE(6,*) "Doing FullLon..."
    endif
    imin=1
    imax=nlon
    IF( DO_OCEAN ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( DO_DRIFTERS ) CALL obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)
  endif
  nn = nn -1 
  nn_dr = nn_dr-1
  !STEVE: end rewrite
  
  ! If no observations remain, then we're done.
  IF(nn < 1) THEN
    nobsl = 0
    IF(nn_dr < 1) THEN
      nobsl_dr = 0
      RETURN
    END IF
  END IF

  if (dodebug) WRITE(6,*) "Start doing localization..."
!
! CONVENTIONAL
!
 nobsl = 0
 IF ( DO_OCEAN ) THEN ! profile observation exists
  DO n=1,nn
    !STEVE: (future) use custom localization with CGAL/BOOST algorithms
    !
    ! Observational localization Distance Cutoff
    !
    olat = obslat(nobs_use(n))
    olon = obslon(nobs_use(n))
    olev = obslev(nobs_use(n))
    oelm = obselm(nobs_use(n))

    ! variable localization
    !
    SELECT CASE(NINT(oelm))
    CASE(id_u_obs)
        iobs=1
    CASE(id_v_obs)
        iobs=2
    CASE(id_t_obs)
        iobs=3
    CASE(id_s_obs)   !(OCEAN)
        iobs=4
    CASE(id_ssh_obs) !(OCEAN)
        iobs=5
    CASE(id_sst_obs) !(OCEAN)
        iobs=6
    CASE(id_sss_obs) !(OCEAN)
        iobs=7
    CASE(id_eta_obs) !(OCEAN)
        iobs=8
    CASE DEFAULT
        WRITE(6,*) "letkf_local.f90 :: there is no variable localization for obs-type :: ", oelm
        WRITE(6,*) "letkf_local.f90 :: exitING..."
        STOP(95)
    END SELECT
    IF(var_local(iobs) < TINY(var_local)) CYCLE 
    !STEVE: skip obs that are set to "0" impact on this model variable.

    !--------------------------------------------------------------------
    ! vertical localization
    !--------------------------------------------------------------------
    !STEVE: make vertical localization depth dependent
    !STEVE: Could alternatively treat surface obs as if they occur at depth 0, regardless of model's top level
    IF(NINT(oelm) == id_ssh_obs) THEN
        !dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))        !(OCEAN)
        !STEVE: no need for log scale, but, don't have ssh obs right now... so I'll test this later
        !dlev = ABS( olev - lev(ilev) )
        dlev = 0.0d0 !STEVE: allow all levels to be influenced by ssh
    ELSE IF(NINT(oelm) == id_sst_obs) THEN
        dlev = ABS( olev - xlev )
        !WRITE(6,*) "obs_local:: SST ob, dlev = ",dlev
    ELSE IF(NINT(oelm) == id_sss_obs) THEN
        dlev = ABS( olev - xlev )
    ELSE IF(NINT(oelm) == id_u_obs .or. &
            NINT(oelm) == id_v_obs .or. &
            NINT(oelm) == id_t_obs .or. &
            NINT(oelm) == id_s_obs   ) THEN
        dlev = ABS( olev - xlev )
    ELSE
        dlev = 0.0d0
    END IF

    
    IF(dlev > dist_zerov) CYCLE

    !
    ! horizontal localization
    !
    horizontal_localization : if (localization_method .eq. 1) then
        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
        sigma_b = sigma_a
        dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    elseif (localization_method .eq. 2) then
        !STEVE: THIS LOCALIZATION METHOD WAS AN ERROR, however the results were
        !good, and better than a number of different horizontal localization
        !schemes that I tried. So I am keeping it here for further study.
        !
        ! The input sigma_obs should be 1000km, and simga_obs0 should be 200km

        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 !100.0 * 1000.0 !make the minimum distance 100 kilometers at 90ºN or 90ºS
                                  !WARNING: this should be bigger than the minimum grid cell width
!       dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad
        sigma_b = sigma_a
        dist_zero_a = sigma_a * (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_b = sigma_a
!       dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    else !STEVE: use the original (default) localization
        CALL com_distll_1(olon,olat,xlon,xlat,dist)
        IF(dist > dist_zero ) CYCLE
    endif horizontal_localization

    ! STEVE: ADD CHECK FOR OBSERVATIONS THAT ARE OCCLUDED BY LAND!
    ! STOLE atlpac from SODA, in the future, implement a general method
    ! for occluding points that can't be reached by a SPM in the local range
    blocked_by_land = .false.
    CALL atlpac(xlat,xlon,lxpa)
    CALL atlpac(olat,olon,lopa)
    IF (  lxpa .eq. 1 .and. lopa .eq. 2 &        !in caribbean, pacific and atlantic
       .or. lxpa .eq. 2 .and. lopa .eq. 1 &        !in caribbean, atlantic and pacific
         ) blocked_by_land = .true.
    IF (blocked_by_land) CYCLE                   !STEVE: the two points are in atl and pac

!--------- End of Observation Culling ---------!

    nobsl = nobsl + 1
    hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
    dep(nobsl) = obsdep(nobs_use(n))
    tmperr=obserr(nobs_use(n))
    rdiag(nobsl) = tmperr * tmperr
    if (ALLOCATED(obs_useidx)) then
      obs_useidx(nobsl) = nobs_use(n)
    endif
    !
    ! Observational localization (weighting)
    !
    ! Note: var_local scales the localization weighting based on the parameter
    ! set in letkf_tools.f90. A row corresponding to the model parameter is
    ! input to this subroutine, and the column indicates the proportion of that
    ! type of observation to use.
    observation_localization : if (localization_method .eq. 1 ) then
                                                       
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
        ELSE
          rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF

    elseif (localization_method .eq. 2) then
      
        !STEVE: trying the original approach that seemed to work before:
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_obs)**2) * var_local(iobs)
        ELSE
          rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF
                                                       
    elseif (localization_method .eq. 0) then
        !STEVE: this is R^2 and the R localization gaussian function
        rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
    else
        print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
        STOP(3)
    endif observation_localization
  END DO
 END IF ! test DO_OCEAN

  !########################################################################
  !LUYU: select drifters observation
  !STEVE: most of the localization section has been completely edited for (OCEAN)
  nobsl_dr = 0
 IF ( DO_DRIFTERS ) THEN
  DO n=1,nn_dr
    !STEVE: (future) use custom localization with CGAL/BOOST algorithms
    !
    ! Observational localization Distance Cutoff
    !
    olat = obslat_dr(nobs_use_dr(n))
    olon = obslon_dr(nobs_use_dr(n))
    olev = obslev_dr(nobs_use_dr(n))
    oelm = obselm_dr(nobs_use_dr(n))

    !STEVE: debug
   !if ( i1(ij) .eq. 456 .and. j1(ij) .eq. 319 .and. ilev .eq. 5) then
!   if ( (((i1(ij)) .eq. 456 .AND. (j1(ij)) .eq. 319) .OR.  (ij .eq. 478)) .AND. ilev .eq. 5) then
!     ! Skipping salinity obs for this grid point to see if it fixes weird analysis result
!     WRITE(6,*) "letkf_local.f90:: i=456,j=319,k=5 :: oelm, NINT(oelm) = ", oelm, NINT(oelm)
!     if ( NINT(oelm) .eq. id_s_obs ) then
!       WRITE(6,*) "letkf_local.f90:: removing salinity observation..."
!       CYCLE
!     elseif ( NINT(oelm) .eq. id_sst_obs ) then
!       WRITE(6,*) "letkf_local.f90:: removing SST observation..."
!       CYCLE
!     endif
!   endif

    !
    ! variable localization
    !
    SELECT CASE(NINT(oelm))
    CASE(id_x_obs)
        iobs=1
    CASE(id_y_obs)
        iobs=2
    CASE(id_z_obs)
        iobs=3
    CASE(id_t_obs)   !(OCEAN)
        iobs=4
    CASE(id_s_obs) !(OCEAN)
        iobs=5
    END SELECT
    IF(var_local(iobs) < TINY(var_local)) CYCLE 
    !STEVE: skip obs that are set to "0" impact on this model variable.

    !
    ! vertical localization
    !
    if (dodebug) WRITE(6,*) "vertical localization..."
    !STEVE: make vertical localization depth dependent
    !STEVE: Could alternatively treat surface obs as if they occur at depth 0, regardless of model's top level
    IF(NINT(oelm) == id_ssh_obs) THEN
        !dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))        !(OCEAN)
        !STEVE: no need for log scale, but, don't have ssh obs right now... so I'll test this later
        !dlev = ABS( olev - lev(ilev) )
        dlev = 0.0d0 !STEVE: allow all levels to be influenced by ssh
    ELSE IF(NINT(oelm) == id_sst_obs) THEN
        dlev = ABS( olev - xlev )
        !WRITE(6,*) "obs_local:: SST ob, dlev = ",dlev
    ELSE IF(NINT(oelm) == id_sss_obs) THEN
        dlev = ABS( olev - xlev )
    ELSE IF(NINT(oelm) == id_x_obs .or. &
            NINT(oelm) == id_y_obs .or. &
            NINT(oelm) == id_z_obs .or. &
            NINT(oelm) == id_t_obs .or. &
            NINT(oelm) == id_s_obs   ) THEN
        dlev = ABS( olev - xlev )
    ELSE
        dlev = 0.0d0
    END IF
    IF(dlev > dist_zerov) CYCLE

    !
    ! horizontal localization
    !
    if (dodebug) WRITE(6,*) "horizontal localization..."
    horizontal_localization_dr : if (localization_method .eq. 1) then
        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
        sigma_b = sigma_a
        dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    elseif (localization_method .eq. 2) then
        !STEVE: THIS LOCALIZATION METHOD WAS AN ERROR, however the results were
        !good, and better than a number of different horizontal localization
        !schemes that I tried. So I am keeping it here for further study.
        !
        ! The input sigma_obs should be 1000km, and simga_obs0 should be 200km

        !STEVE: make horizontal localization latitude dependent
        ! STEVE: make sigma_obs a linear/lookup function of latitude
        dist_min = sigma_obs0 !100.0 * 1000.0 !make the minimum distance 100 kilometers at 90ºN or 90ºS
                                  !WARNING: this should be bigger than the minimum grid cell width
!       dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
        minr = dist_min
        maxr = dist_zero
        maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
        maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

        ! Shrink radius far from equator
        ! Shrink foci far from equator
        xdis = abs(xlat)*(pi/180.0d0)*re
        if (xlat >= 0) then
          xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
        else
          xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
        endif

        sigma_a = xrad
        sigma_b = sigma_a
        dist_zero_a = sigma_a * (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
!       sigma_b = sigma_a
!       dist_zero_a = xrad

        CALL com_distll_1(olon,olat,xlon,xlat,dist)

        if (.false.) then !(modulo(n,10000) .eq. 0) then
          !STEVE: DEBUG: print out the important values for a selection of obs
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
          print *, "DEBUG letkf_local.f90..."
          print *, "n = ", n
          print *, "dist_min = ", dist_min
          print *, "sigma_obs0 = ", sigma_obs0
          print *, "minr = ", minr
          print *, "maxr = ", maxr
          print *, "maxdN = ", maxdN
          print *, "maxdS = ", maxdS
          print *, "xlat = ", xlat
          print *, "xdis = ", xdis
          print *, "xrad = ", xrad
          print *, "sigma_a = ", sigma_a
          print *, "sigma_b = ", sigma_b
          print *, "sigma_obs = ", sigma_obs
          print *, "dist_zero_a = ", dist_zero_a
          print *, "dist = ", dist
!         print *, "dist > dist_zero_a = ", dist > dist_zero_a
          print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        endif

        IF(dist > dist_zero_a) CYCLE  !Points are outside of radius

    else !STEVE: use the original (default) localization
        CALL com_distll_1(olon,olat,xlon,xlat,dist)
        IF(dist > dist_zero ) CYCLE
    endif horizontal_localization_dr

    ! STEVE: ADD CHECK FOR OBSERVATIONS THAT ARE OCCLUDED BY LAND!
    ! STOLE atlpac from SODA, in the future, implement a general method
    ! for occluding points that can't be reached by a SPM in the local range
    blocked_by_land = .false.
    CALL atlpac(xlat,xlon,lxpa)
    CALL atlpac(olat,olon,lopa)
    IF (  lxpa .eq. 1 .and. lopa .eq. 2 &        !in caribbean, pacific and atlantic
       .or. lxpa .eq. 2 .and. lopa .eq. 1 &        !in caribbean, atlantic and pacific
         ) blocked_by_land = .true.
    IF (blocked_by_land) CYCLE                   !STEVE: the two points are in atl and pac

!--------- End of Observation Culling ---------!
    
    nobsl_dr = nobsl_dr + 1
    hdxf(nobsl+nobsl_dr,:) = obshdxf_dr(nobs_use_dr(n),:)
    dep(nobsl+nobsl_dr)    = obsdep_dr(nobs_use_dr(n))
    tmperr        = obserr_dr(nobs_use_dr(n))
    rdiag(nobsl+nobsl_dr) = tmperr * tmperr
    if (ALLOCATED(obs_useidx)) then
      obs_useidx(nobsl+nobsl_dr) = nobs_use_dr(n)
      if (dodebug) WRITE(6,*) "obs_useidx(nobsl_dr) = ", obs_useidx(nobsl+nobsl_dr)
    endif
    !
    ! Observational localization (weighting)
    !
    ! Note: var_local scales the localization weighting based on the parameter
    ! set in letkf_tools.f90. A row corresponding to the model parameter is
    ! input to this subroutine, and the column indicates the proportion of that
    ! type of observation to use.
    if (dodebug) WRITE(6,*) "observation localization..."
    observation_localization_dr : if (localization_method .eq. 1 ) then
!       rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a  )**2 + (dlev/sigma_obsv)**2)) &
!                                                    & * var_local(iobs)
!STEVE: testing different localization functions:
!STEVE: doubing standard deviation distance
!       rloc(nobsl) =EXP( -0.5d0 * ( ( dist/(2.0d0*sigma_a) )**2 + ( dlev/(2.0d0*sigma_obsv) )**2 ) ) &
!                                                    & * var_local(iobs)
!STEVE: removing all localization other than culling as applied above based on
!absolute distance:
!       rloc(nobsl) = 1.0d0 * var_local(iobs)
                                                       
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
        ELSE
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF

    elseif (localization_method .eq. 2) then
      
        !STEVE: trying the original approach that seemed to work before:
        IF (DO_NO_VERT_LOC) THEN
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * (dist/sigma_obs)**2) * var_local(iobs)
        ELSE
          rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
        ENDIF
                                                       
    elseif (localization_method .eq. 0) then
        !STEVE: this is R^2 and the R localization gaussian function
        rloc(nobsl+nobsl_dr) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
    else
        print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
        STOP(3)
    endif observation_localization_dr
  END DO
 END IF ! test DO_DRIFTERS

  !STEVE: this should never happen, if it does something went wrong
  IF( nobsl_dr + nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL_DR + NOBSL =',nobsl_dr+nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'ID,IT,NN_DR,NN,TVNN=', id, it, nn_dr, nn, tvnn
    STOP 99
  END IF
 
  IF( nobs_dr > 0 ) THEN
    DEALLOCATE(nobs_use_dr)
  END IF

  RETURN
  
END SUBROUTINE obs_local_drifters

PURE SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!===============================================================================
! Identify the observations within the local region
!===============================================================================
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  ! Cycle through each latitude
  do j=jmin,jmax
    ! Find the number of accumulated obs at the bottom of the range
    if(imin > 1) then
      ib = nobsgrd(imin-1,j)+1
    else
      ! Wrap around to the previous latitude at the last longitude
      if(j > 1) then
        ib = nobsgrd(nlon,j-1)+1
      else
        ib = 1
      endif
    endif
    ! Find the number of accumulated obs at the top of the range
    ie = nobsgrd(imax,j)

    ! Subtract to get the number of obs in this region
    n = ie - ib + 1

    if(n == 0) CYCLE !there are no obs here

    do ip=ib,ie
      if(nn > nobs) then
!       WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!       stop 1  !STEVE: (added)
      endif
      ! Index for observation used
      nobs_use(nn) = ip
      ! Count up the total obs used so far
      nn = nn + 1
    enddo
  enddo

  RETURN
END SUBROUTINE obs_local_sub

PURE SUBROUTINE obs_local_sub_dr(imin,imax,jmin,jmax,nn_dr,nobs_use_dr)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn_dr, nobs_use_dr(nobs_dr)
  INTEGER :: j,n,ib,ie,ip

  ! Cycle through each latitude
  DO j=jmin,jmax
    ! Find the number of accumulated obs at the bottom of the range
    IF(imin > 1) THEN
      ib = nobsgrd_dr(imin-1,j)+1
    ELSE
      ! Wrap around to the previous latitude at the last longitude
      IF(j > 1) THEN
        ib = nobsgrd_dr(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ! Find the number of accumulated obs at the top of the range
    ie = nobsgrd_dr(imax,j)

    ! Subtract to get the number of obs in this region
    n = ie - ib + 1

    IF(n == 0) CYCLE !there are no obs here

    DO ip=ib,ie
      IF(nn_dr > nobs_dr) THEN
!        WRITE(6,*) 'FATALERROR, NN_DR > NOBS + NOBS_DR: ', nn_dr,'>', nobs_dr+
!        stop 1  !STEVE: (added)
      END IF
      ! Index for observation used
      nobs_use_dr(nn_dr) = ip
      ! Count up the total obs used so far
      nn_dr = nn_dr + 1
    END DO
  END DO

  RETURN
END SUBROUTINE obs_local_sub_dr

! Set up kd-tree with observation data
SUBROUTINE kdtree

END SUBROUTINE kdtree

! Scan graph of model grid with search algorithm to avoid land in localization
SUBROUTINE cullBlocked
! Inputs:
!        model grid (grid or graph form)
!        land/sea map (or kmt data)
!        list of observations in range
!
! Outputs:
!        list of observations not blocked by land 
!


END SUBROUTINE cullBlocked

! Link to C++ Boost library for fast A* algorithm
SUBROUTINE Astar

END SUBROUTINE Astar

! Graph representation of model grid
SUBROUTINE grid2graph
! Just create a linked list that contains the info needed to access information on grid
! Use grid_graph from Boost: http://www.boost.org/doc/libs/1_46_1/libs/graph/doc/grid_graph.html
!

END SUBROUTINE grid2graph

!(OCEAN) STEVE: add checks for atlantic/pacific basin boundary
subroutine atlpac (xlat, xlon_in, lxap)
REAL(r_size), INTENT(IN) :: xlat, xlon_in
REAL(r_size), INTENT(OUT) :: lxap
REAL(r_size) :: xlon

! Ensure the comparisons are on the 0-360.0 longitude range
xlon = modulo(xlon_in,360.0)

!c STEVE: Stolen from SODA: use until a general method for managing land-blocked ocean basins...
!c=================================================================
!c X. Cao 12/9/99
!c
!c   to make a mark to the location of a point in Caribbean area
!c (xlat.gt.-2..and.xlat.lt.32..and.xlon.gt.245..and.xlon.lt.295.)
!c to indicate if it is in Atlantic ocean (1) or in Pacific ocean (2)
!c or out of Caribbean area (0)
!c=================================================================
!c
  lxap=0
!c
!c -- Atlantic ? Pacific?
!c
  if(xlat.gt.-2..and.xlat.le.8.5) then
    if(xlon.lt.285.) then
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.8.5.and.xlat.le.15.5) then
    if(xlon.lt.276.) then
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.15.5.and.xlat.le.19.5) then
    if(xlon.lt.270.) then
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.19.5.and.xlat.le.32.0) then
    if(xlon.lt.258.) then
      lxap=2
    else
      lxap=1
    endif
  endif
  return
end subroutine atlpac

END MODULE letkf_drifters_local
