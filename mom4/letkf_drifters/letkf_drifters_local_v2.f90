MODULE letkf_drifters_local
  !STEVE: the letkf_local had a few details that were specific to grid-based obs,
  !       so I copied it here and made the requried changes for the drifter Lagrangian obs.
  !
  USE common
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE letkf_obs, ONLY: dist_zero, dist_zerov !contains debug_hdxf_0, and nobsgrd
  USE params_letkf, ONLY: nbv, DO_NO_VERT_LOC, localization_method
  USE params_obs
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

  REAL(r_size),PARAMETER :: var_local(nv4d,nid_dobs) = 1.0d0
   !NOTE: the obs are the rows and the model variables the columns

  INTEGER,SAVE :: var_local_n2n(nv4d)

CONTAINS

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
!SUBROUTINE obs_local(ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl)
SUBROUTINE obs_local(ij,it,var_local,hdxf,rdiag,rloc,dep,nobsl,nobstotal)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,it,nobstotal !,nvar
  REAL(r_size),INTENT(IN) :: var_local(nid_dobs)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:),nobs_use_tmp(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: i,n,nn,tvnn,iobs,m,mm,isfirst
  !STEVE: for (OCEAN):
  LOGICAL :: blocked_by_land     !Use these three lines to identify gulf vs. pacific points in localization
  REAL(r_size) :: xlat,xlon,lxpa,ylon,ylat,xlev,ylev
  REAL(r_size) :: olat,olon,olev,oelm,lopa
  REAL(r_size) :: f1lon,f1lat,f2lon,f2lat
  REAL(r_size) :: dist1,dist2,a,b,u,v,lu,lv,dist_d
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
  LOGICAL :: dodebug = .true.
  !debug TEST:
  INTEGER :: nn_old,imin_old,imax_old,jmin_old,jmax_old,id,itime,ind
  REAL(r_size) :: minlon_old,maxlon_old,minlat_old,maxlat_old

!
! INITIALIZE

  IF( nobs_dr > 0 ) THEN
    ALLOCATE(nobs_use(nobs_dr))
    ALLOCATE(nobs_use_tmp(nobs_dr)) ! This array is designed to collect all the indices with the same time
  ELSE
    WRITE(6,*) "WARNING: there is no observation"
  END IF
  nn = 0

  id=NINT(drifter_ids(ij))              ! DRIFTERS id (DEFINE LATER)
  itime=NINT(drifter_times(it))         ! DRIFTERS time (DEFINE LATER)

  ! Use time and id to locate the same drifters in observation space
  ! and also the indices under the same time

  m = 0
  ind = 0
  DO i=1,nobs_dr

    if (dodebug) WRITE(6,*) "obsid(i)=",obsid(i)
    if (dodebug) WRITE(6,*) "obstime(i)=",obstime(i)
    IF ( NINT(obstime(i)) .eq. itime ) THEN
      m = m + 1
      nobs_use(m) = i 
      IF ( NINT(obsid(i) ) .eq. id ) THEN       
        ind = i
      END IF
    END IF
  ENDDO

  if (dodebug) WRITE(6,*) "nobstl_tmp=", m, "ind=", ind

  IF (m<1 .or. ind .eq. 0) THEN
    nobsl=0
    WRITE(6,*) "NOTE: There is no observation for id=", id, " and time = ", itime
    RETURN
  END IF

  xlon=obslon_dr(ind)
  xlat=obslat_dr(ind)
  
  IF (dodebug) WRITE(6,*) "Selecting verified observation..."
  IF (dodebug) WRITE(6,*) "dist_zero=",dist_zero,"dist_zerov=",dist_zerov
!  mm = 0
  ! Find all the indices in the ball with radius equal to sigma_zero  
!  DO i = 1,m     ! m, the number of observations at the time it.
!    ylon = obslon_dr(nobs_use_tmp(i))
!    ylat = obslat_dr(nobs_use_tmp(i))
!    IF (dodebug) WRITE(6,*) "ylon=",ylon,"ylat=",ylat
!    dist = SQRT((xlon-ylon)**2+(xlat-ylat)**2)
!    ylev = obslev_dr(nobs_use_tmp(i))
!    dlev = ABS(xlev - ylev)
!    IF (dodebug) WRITE(6,*) "dist_d=",dist_d
!    IF ( dist_d <= dist_zero ) THEN
!      mm = mm + 1
!      nobs_use(mm) = nobs_use_tmp(i)
!      if (dodebug) WRITE(6,*) "index=", nobs_use(mm)
!    END IF
!  ENDDO

!
! CONVENTIONAL
!

!  IF (mm <1) THEN
!    nobsl = 0
!    RETURN
!  END IF
  !STEVE: most of the localization section has been completely edited for (OCEAN)
  nobsl = 0
  DO n=1,m
    !STEVE: (future) use custom localization with CGAL/BOOST algorithms
    !
    ! Observational localization Distance Cutoff
    !
    if (dodebug) WRITE(6,*) "n=",n
    olat = obslat_dr(nobs_use(n))
    olon = obslon_dr(nobs_use(n))
    olev = obslev_dr(nobs_use(n))
    oelm = obselm_dr(nobs_use(n))

    !
    ! variable localization
    !
    SELECT CASE(INT(oelm))
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
     
    dlev=0
    if(dlev > dist_zerov) CYCLE

    !---------------------------------------------------------------------------
    ! horizontal localization
    !---------------------------------------------------------------------------
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

      if (dodebug) WRITE(6,*) "dist=",dist,"dist_zero_a=",dist_zero_a

      if (dist > 2 * dist_zero_a) CYCLE
   else

      CALL com_distll_1(olon,olat,xlon,xlat,dist)
      if (dist > 2 * dist_zero) CYCLE

    endif horizontal_localization


!--------- End of Observation Culling ---------!

    nobsl = nobsl + 1
    hdxf(nobsl,:) = obshdxf_dr(nobs_use(n),:)
    dep(nobsl)    = obsdep_dr(nobs_use(n))
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
        rloc(nobsl) =EXP(-0.5d0 * ((dist_d/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                     & * var_local(iobs)
    else
        print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
        STOP(3)
    endif observation_localization
  END DO

  !STEVE: this should never happen, if it does something went wrong
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN,TVNN=', ij, nn, tvnn
    STOP 99
  END IF
 
  IF( nobs_dr > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF

  RETURN
END SUBROUTINE obs_local

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  ! Cycle through each latitude
  DO j=jmin,jmax
    ! Find the number of accumulated obs at the bottom of the range
    IF(imin > 1) THEN
      ib = nobsgrd(imin-1,j)+1
    ELSE
      ! Wrap around to the previous latitude at the last longitude
      IF(j > 1) THEN
        ib = nobsgrd(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ! Find the number of accumulated obs at the top of the range
    ie = nobsgrd(imax,j)

    ! Subtract to get the number of obs in this region
    n = ie - ib + 1

    IF(n == 0) CYCLE !there are no obs here

    DO ip=ib,ie
      IF(nn > nobs_dr) THEN
        WRITE(6,*) 'FATALERROR, NN > NOBS', nn, nobs_dr
        stop 1  !STEVE: (added)
      END IF
      ! Index for observation used
      nobs_use(nn) = ip
      ! Count up the total obs used so far
      nn = nn + 1
    END DO
  END DO

  RETURN
END SUBROUTINE obs_local_sub

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
