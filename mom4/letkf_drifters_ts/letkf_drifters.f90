MODULE letkf_drifters

USE netcdf
USE common
USE common_mpi
USE common_mom4
USE common_mpi_mom4
USE common_letkf
USE letkf_drifters_obs, ONLY : debug_hdxf_0
USE params_letkf
USE params_model
USE params_obs
USE letkf_drifters_tools
UsE letkf_drifters_local

CONTAINS

SUBROUTINE das_drifters(gues4d,gues3d,gues2d,anal4d,anal3d,anal2d)
  ! das_letkf went through all model grid points. Now we're going to go through
  ! each drifter_id, as if we had appended them to the model state vector.
  
  USE common_letkf, ONLY: letkf_core
  USE params_model

  IMPLICIT NONE
  CHARACTER(12) :: inflinfile='infl_mul.grd'
  CHARACTER(12) :: infloutfile='infl_out.grd'
  REAL(r_size),INTENT(INOUT) :: gues4d(nid1,num_times,nbv,nv4d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d)
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)
  REAL(r_size),INTENT(OUT) :: anal4d(nid1,num_times,nbv,nv4d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean4d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans_dr(nbv,nbv,nv3d+nv2d+nv4d)
  LOGICAL :: ex
  INTEGER :: ij,ilev
  INTEGER :: id,it,n,m,i,j,k,nobsl,nobsl_dr,ierr
  !STEVE: for debugging
  LOGICAL :: debug_sfckmt = .false.
  LOGICAL :: dodebug = .true.
  INTEGER :: nn
  REAL(r_size) :: maxdep_val
  INTEGER :: maxdep_nn
  REAL(r_size) :: mindep_val
  INTEGER :: mindep_nn
  !STEVE: for DO_NO_VERT_LOC
  INTEGER :: klev
  !STEVE: added for drifters:
  INTEGER :: nobstotal,itim

  WRITE(6,'(A)') 'Hello from das_drifters'
  nobstotal = nobs + nobs_dr !LUYU: total number of observation from drifters and profiles
  WRITE(6,'(A,I8)') 'Target observation numbers : nobs_dr=',nobs_dr, 'nobs(profile)=', nobs !,', NTVS=',ntvs

  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    anal4d = gues4d      !(DRIFTERS)
    anal3d = gues3d      !(OCEAN)
    anal2d = gues2d      !(OCEAN)
    RETURN
  ELSE                   
    anal4d = 0.0d0       !(DRIFTERS)
    anal3d = 0.0d0       !(OCEAN)
    anal2d = 0.0d0       !(OCEAN)
  END IF

  !-----------------------------------------------------------------------------
  ! Variable localization
  !-----------------------------------------------------------------------------
  var_local_n2n(1) = 1
  do n=2,nv3d+nv2d+nv4d
    do i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    enddo
  enddo
  WRITE(6,*) "var_local_n2n = ", var_local_n2n

  !-------------------------------------------------------------------------
  ! FCST PERTURBATIONS.
  !-------------------------------------------------------------------------
  ! LUYU: From now on, gues*d = X^f, the perturbation matrix
  !

  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)

  do n=1,nv3d
    do m=1,nbv
      do k=1,nlev
        do i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        enddo
      enddo
    enddo
  enddo
  do n=1,nv2d
    do m=1,nbv
      do i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      enddo
    enddo
  enddo

  !
  ! multiplicative inflation
  !
  !STEVE: (using simple constant inflation for now.)
  !***************************************************************************
  !LUYU: (**TBD**, Just COPY following codes for now)
  if (cov_infl_mul > 0.0d0) then ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = cov_infl_mul
    work2d = cov_infl_mul
  endif
  if (cov_infl_mul <= 0.0d0) then ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    INQUIRE(FILE=inflinfile,EXIST=ex)
    if (ex) then
      if (myrank == 0) then
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflinfile
        CALL read_bingrd4(inflinfile,work3dg,work2dg)
      endif
      CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
    else
      WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflinfile
      work3d = -1.0d0 * cov_infl_mul
      work2d = -1.0d0 * cov_infl_mul
    endif
  endif

  !-----------------------------------------------------------------------------
  ! Reset inflation, if desired
  !-----------------------------------------------------------------------------
  if ( DO_INFL_RESET ) then
    work3d = 1.0d0
    !work2d = 1.0d0 !STEVE: don't reset the 2D field
    !STEVE: will this mess up SST, SSH and SSS? maybe should reset them also.
    !       not sure if it matters since they are not prognostic in the model.

    ! Eventually, I'd like to back up the adaptive inflation from the depths to
    ! the corresponding surface level, following the water mass generation, but
    ! this would be a lot more complicated.

    ! If using the Hybrid-LETKF, multiplicative inflation is not necessary
  endif
  !***************************************************************************

  !
  ! MAIN ASSIMILATION LOOP
  !
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE(hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  WRITE(6,*) "... done."

  ! LUYU: ***First, consider local observations around ij.
  ! LUYU: ***In this part, we directly copy codes from the part in letkf_tools.f90: das_letkf....(At least, most of it....)
  !##########################################################################
  do ilev=1,nlev
    WRITE(6,'(A,I3)') 'ilev = ',ilev

    if (DO_NO_VERT_LOC .and. ilev > 1) CYCLE
          
    do ij=1,nij1 !STEVE: go through every possible coordinate of the grid in list form...
      if (dodebug) WRITE(6,*) "ij = ", ij

      !STEVE: debug
!     if (.false. .AND. ilev == 1 .AND. NINT(i1(ij)) == 131 .AND. NINT(j1(ij)) == 76) then
!       do m=1,nbv
!         WRITE(6,*) "letkf_tools.f90"
!         WRITE(6,*) "kmt1(ij) = ", kmt1(ij)
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "i1(ij) = ", i1(ij)
!         WRITE(6,*) "j1(ij) = ", j1(ij)
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!         WRITE(6,*) "should be > 20 in this region"
!       enddo
!     endif

      !(OCEAN) STEVE: it's on land, so just assign undef values and CYCLE
      !STEVE: NEED to define kmt1 as ocean depth
      if (kmt1(ij) < ilev) then
        anal3d(ij,ilev,:,:) = 0.0
        work3d(ij,ilev,:) = 0.0
        if (ilev == 1) then
          anal2d(ij,:,:) = 0.0
          work2d(ij,:) = 0.0 !STEVE: added
        endif
        !STEVE: debug
        !WRITE(6,*) "CYCLE: kmt1(ij) = ", kmt1(ij)
        CYCLE
      endif

      !WRITE(6,*) "ASSIM: kmt1(ij) = ", kmt1(ij)
      !(OCEAN) STEVE:end
       
!     if (ilev == 1) then 
!     do m=1,nbv
!       if (NINT(i1(ij)) == 131 .AND. NINT(j1(ij)) == 76) then
!         WRITE(6,*) "letkf_tools.f90"
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "i1(ij) = ", i1(ij)
!         WRITE(6,*) "j1(ij) = ", j1(ij)
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!         WRITE(6,*) "should be > 20 in this region"
!       endif
!       if (debug_sfckmt .AND. ABS(gues2d(ij,m,iv2d_sst)-gues3d(ij,1,m,iv3d_t)) > TINY(1.0)) then
!         WRITE(6,*) "letkf_tools.f90:: SST does not equal SFC T" 
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!!        stop 4
!       endif
!     enddo
!     endif

      !-------------------------------------------------------------------------
      ! Loop through all prognostic variables (e.g. temp, salt, u, v, etc.)
      !-------------------------------------------------------------------------
      do n=1,nv3d
        if (var_local_n2n(n) < n) then
          trans_dr(:,:,n) = trans_dr(:,:,var_local_n2n(n))
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        else
          CALL obs_local(ij,ilev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobsl_dr,nobstotal)
          WRITE(6,*) "Finish calling obs_local..."

          parm = work3d(ij,ilev,n)
          WRITE(6,*) "ilev=",ilev,"n=", n, "work3d=", work3d(ij,ilev,n)

          debug_local_obs : if ( .false. ) then !NINT(i1(ij)) .eq. 456 .and. NINT(j1(ij)) .eq. 319 .and. ilev .eq. 5) then
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "nobsl = ", nobsl, "nobsl_dr = ", nobsl_dr
            WRITE(6,*) "i1(ij), j1(ij), ilev = ", i1(ij), j1(ij), ilev
            WRITE(6,*) "ij,n,var_local_n2n(n) = ", ij,n,var_local_n2n(n)
            maxdep_val = 0
            maxdep_nn = 0
            mindep_val = 0
            mindep_nn = 0
            do nn = 1 ,nobsl_dr
             !if (obselm(obs_useidx(nn)) .eq. id_t_obs .OR. obselm(obs_useidx(nn)) .eq. id_s_obs .and. ALLOCATED(obs_useidx)) then
               
                ! LUYU: NOTE: Here we only assume the obersavation data are all Lagrangian Data(Drifters) . If also including profile, we need to modify later...(**TBD**)
                WRITE(6,*) "---------- (letkf_tools.f90)"
                WRITE(6,*) "nn = ", nn
                if (NINT(obselm_dr(obs_useidx(nn))) .eq. id_x_obs) WRITE(6,*) "DRIFTERS_X_OB"
                if (NINT(obselm_dr(obs_useidx(nn))) .eq. id_y_obs) WRITE(6,*) "DRIFTERS_Y_OB"
                if (NINT(obselm_dr(obs_useidx(nn))) .eq. id_z_obs) WRITE(6,*) "DRIFTERS_Z_OB"
                WRITE(6,*) "obslon_dr/obslat_dr/obslev_dr(obs_useidx(nn)) = ", obslon_dr(obs_useidx(nn)), obslat_dr(obs_useidx(nn)), obslev_dr(obs_useidx(nn))
                WRITE(6,*) "obs_useidx(nn) = ", obs_useidx(nn)
                WRITE(6,*) "obselm_dr(obs_useidx(nn)) = ", obselm_dr(obs_useidx(nn))
                if (obsdat(obs_useidx(nn)) < 0) WRITE(6,*) "NEGOB!"
                WRITE(6,*) "obsid(obs_useidx(nn)) = ", obsid(obs_useidx(nn))
                WRITE(6,*) "dep(nn)    = ", dep(nn)
                WRITE(6,*) "obserr_dr(obs_useidx(nn)) = ", obserr_dr(obs_useidx(nn))
                WRITE(6,*) "rdiag(nn)    = ", rdiag(nn)
                WRITE(6,*) "rloc(nn)    = ", rloc(nn)
             
!               WRITE(6,*) "obsdep(obs_useidx(nn))    = ", obsdep(obs_useidx(nn))
                if (dep(nn) > maxdep_val) then
                  maxdep_val = dep(nn)
                  maxdep_nn = nn
                endif
                if (dep(nn) < mindep_val) then
                  mindep_val = dep(nn)
                  mindep_nn = nn
                endif
                WRITE(6,*) "hdxf(nn,:) = "
                WRITE(6,*) hdxf(nn,1:nbv)
!               WRITE(6,*) "obshdxf(obs_useidx(nn),:) = "
!               WRITE(6,*) obshdxf(obs_useidx(nn),1:nbv)
             !endif
            enddo
            WRITE(6,*) "maxdep_val = ", maxdep_val
            WRITE(6,*) "maxdep_nn  = ", maxdep_nn
            WRITE(6,*) "mindep_val = ", mindep_val
            WRITE(6,*) "mindep_nn  = ", mindep_nn
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
          endif debug_local_obs

          !STEVE: some critical checks to make sure letkf_core equations are valid:
          !STEVE: check to make sure input inflation parameter is valid
!         if ( isnan(parm) ) then
!           WRITE(6,*) "parm = work3d(ij,ilev,n) = ", parm
!           WRITE(6,*) "ij = ", ij
!           WRITE(6,*) "ilev = ", ilev
!           WRITE(6,*) "n = ", n
!         endif
          !STEVE: this shouldn't really happen, but fix it if it does...
          if ( parm == 0 ) then
            WRITE(6,*) "parm = work3d(ij,ilev,n) = ", parm
            WRITE(6,*) "ij = ", ij
            WRITE(6,*) "n = ", n
            WRITE(6,*) "letkf_tools.f90:: pre-letkf_core, parm changed to ABS(cov_infl_mul)"
            parm = ABS(cov_infl_mul)
          endif
          !STEVE: check rdiag for > 0
          if (MINVAL(rdiag(1:nobsl)) .le. 0) then
            WRITE(6,*) "letkf.f90:: after obs_local, before letkf_core, for 3D, MINVAL(rdiag(1:nobsl)) ≤ 0"
            WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
          endif

          !STEVE: debug
          if ( debug_hdxf_0 .AND. MINVAL(hdxf(1:nobsl,1:nbv)) == 0 ) then
            WRITE(6,*) "letkf_tools.f90:: (3D) ij = ", ij
            WRITE(6,*) "letkf_tools.f90:: inputs to letkf_core:"
            WRITE(6,*) "nobstotal = ", nobstotal
            WRITE(6,*) "nobsl = ", nobsl
            WRITE(6,*) "hdxf(1:nobsl,1:nbv) = ", hdxf(1:nobsl,:)
            WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
            WRITE(6,*) "rloc(1:nobsl) = ", rloc(1:nobsl)
            WRITE(6,*) "dep(1:nobsl) = ", dep(1:nobsl) 
            WRITE(6,*) "parm = ", parm
            WRITE(6,*) "trans_dr(:,:,n) = ", trans_dr(:,:,n)
          endif
          !STEVE: end

          !LUYU: debug
          if (dodebug) then
            WRITE(6,*) "nobstotal = ", nobstotal , "nobsl = ", nobsl, "nobsl_dr = ", nobsl_dr
            WRITE(6,*) "hdxf(1:nobsl_dr,1:nbv) = ", hdxf(1:nobsl_dr,:)
            WRITE(6,*) "rdiag(1:nobsl_dr) = ", rdiag(1:nobsl_dr)
            WRITE(6,*) "rloc(1:nobsl_dr) = ", rloc(1:nobsl_dr)
            WRITE(6,*) "dep(1:nobsl_dr) = ", dep(1:nobsl_dr) 
            WRITE(6,*) "parm = ", parm
          endif
          !LUYU: enddebug

          !-------------------------------------------------------------------------
          ! Call LETKF MAIN subroutine
          !-------------------------------------------------------------------------
          WRITE(6,*) "Start calling letkf_core...."
          CALL letkf_core(nobstotal,nobsl+ nobsl_dr,hdxf,rdiag,rloc,dep,parm,trans_dr(:,:,n))   !STEVE: need to change for RIP
          WRITE(6,*) "Finish calling letkf_core...."
          WRITE(6,*) "trans_dr(1,1,n)=", trans_dr(1,1,n)
          ! (if doing adaptive inflation)
          work3d(ij,ilev,n) = parm
        endif

        !STEVE: Use the trans_dr matrix computed in letkf_core to form the analysis
        do m=1,nbv
          anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)

          !STEVE: reset analysis to mean for all levels
          if(DO_NO_VERT_LOC .and. ilev .eq. 1) then
            do klev=2,nlev
              anal3d(ij,klev,m,n) = mean3d(ij,klev,n)
            enddo
          endif

          do k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) + gues3d(ij,ilev,k,n) * trans_dr(k,m,n)

            !STEVE: debug - check for bad values
!           if ( anal3d(ij,ilev,m,n) < -10 ) then
!             WRITE(6,*) "Problem in letkf_das after letkf_core. k = ", k
!             WRITE(6,*) "ij, ilev, m, n = ", ij,ilev,m,n
!             WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d(ij,ilev,m,n)
!             WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d(ij,ilev,m,n)
!             STOP 6 
!           endif

            if (DO_NO_VERT_LOC .and. ilev .eq. 1) then
              !STEVE: match up ij to ij at other vertical levels
              do klev=2,nlev
                if (kmt1(ij) < klev) then
                  anal3d(ij,klev,:,:) = 0.0
                  work3d(ij,klev,:) = 0.0
                else
                  anal3d(ij,klev,m,n) = anal3d(ij,klev,m,n) + gues3d(ij,klev,k,n) * trans_dr(k,m,n)
                  work3d(ij,klev,:) = work3d(ij,ilev,:)
                endif
              enddo
            endif
          enddo !k,nbv
          !STEVE: debug
!         if ( i1(ij) .eq. 456 .and. j1(ij) .eq. 319 .and. ilev .eq. 5) then
!           WRITE(6,*) "------------------------------------------------------------"
!           WRITE(6,*) "ij,ilev,m,n = ", ij,ilev,m,n
!           WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d(ij,ilev,m,n) + mean3d(ij,ilev,n)
!           WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d(ij,ilev,m,n)
!           WRITE(6,*) "A-B                 = ", anal3d(ij,ilev,m,n) - gues3d(ij,ilev,m,n) - mean3d(ij,ilev,n)
!           WRITE(6,*) "------------------------------------------------------------"
!         endif
        enddo !m, nbv

      enddo ! n=1,nv3d

      !-------------------------------------------------------------------------
      ! Go through the 2d variables
      !-------------------------------------------------------------------------
      if (ilev == 1) then !update 2d variable at ilev=1
        do n=1,nv2d
          if (var_local_n2n(nv3d+n) <= nv3d) then
            trans_dr(:,:,nv3d+n) = trans_dr(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n))
          elseif (var_local_n2n(nv3d+n) < nv3d+n) then
            trans_dr(:,:,nv3d+n) = trans_dr(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)
          else
            CALL obs_local(ij,ilev,var_local(nv3d+n,:),hdxf,rdiag,rloc,dep,nobsl,nobsl_dr,nobstotal)
            parm = work2d(ij,n)
            !STEVE: check rdiag for > 0
            if (MINVAL(rdiag(1:nobsl)) .le. 0) then
              WRITE(6,*) "letkf.f90:: after obs_local, before letkf_core, for 2D, MINVAL(rdiag(1:nobsl)) ≤ 0"
              WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
            endif

            CALL letkf_core(nobstotal,nobsl+nobsl_dr,hdxf,rdiag,rloc,dep,parm,trans_dr(:,:,nv3d+n)) !STEVE: change for RIP

            !Debugging:
!           print *, "Debugging SFC 2D adaptive inflation:"
!           print *, "pre letkf_core 2D, ilev=1: ij, n, work2d(ij,n) (parm out) = ", ij, n, work2d(ij,n)
            work2d(ij,n) = parm
!           print *, "post letkf_core 2D, ilev=1: ij, n, work2d(ij,n) (parm out) = ", ij, n, work2d(ij,n)
          endif

          !STEVE: process 2D SFC variables here:
          do m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            do k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans_dr(k,m,nv3d+n)
            enddo
          enddo
        enddo

      endif !(ilev == 1)

    enddo !ij
  enddo !ilev

  !-------------------------------------------------------------------------
  ! Write out the adaptive inflation
  !-------------------------------------------------------------------------
  adaptive_inflation : if (cov_infl_mul < 0.0d0) then
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    if (myrank == 0) then
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',infloutfile
      CALL write_bingrd4(infloutfile,work3dg,work2dg)
      !STEVE: check
!     do n=1,nv3d
!       do k=1,nlev
!         do j=1,nlat
!           do i=1,nlon
!             if (isnan4(work3dg(i,j,k,n))) then
!               WRITE(6,*) "writing work3dg(i,j,k,n) = ", work3dg(i,j,k,n)
!               WRITE(6,*) "i,j,k,n = ", i,j,k,n
!               stop 2
!             endif
!           enddo
!         enddo
!       enddo
!     enddo
      !STEVE: end
 
    endif
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  endif adaptive_inflation

  !-------------------------------------------------------------------------
  ! Compute and apply the additive inflation
  !-------------------------------------------------------------------------
  additive_inflation : if (sp_infl_add > 0.0d0) then
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    do n=1,nv3d
      do m=1,nbv
        do k=1,nlev
          do i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          enddo
        enddo
      enddo
    enddo
    do n=1,nv2d
      do m=1,nbv
        do i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        enddo
      enddo
    enddo

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    do n=1,nv3d
      do m=1,nbv
        do ilev=1,nlev
          do ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          enddo
        enddo
      enddo
    enddo
    do n=1,nv2d
      do m=1,nbv
        do ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        enddo
      enddo
    enddo
  endif additive_inflation

  DEALLOCATE(hdxf,rdiag,rloc,dep)

  !########################################################################## 
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE(hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  WRITE(6,*) "... done."

  ALLOCATE(mean4d(nid1,num_times,nv4d))
  CALL ensmean_drifters(nbv,nid1,gues4d,mean4d)

  DO n=1,nv4d
    DO m=1,nbv
      DO k=1,num_times
        DO i=1,nid1
          gues4d(i,k,m,n) = gues4d(i,k,m,n) - mean4d(i,k,n)
        END DO
      END DO
    END DO
  END DO

  ! LUYU: ***Second, find local observations around mean positions of (id,it)
    it = num_times ! LUYU: We just consider the last time index FOR NOW (**TBD**)
    if (dodebug) WRITE(6,*) "it = ", it ! LUYU: note "it" is time index NOT Time
    
    DO i=1,nid1
      id = i + myrank * nid1 
      if (dodebug) WRITE(6,*) "id= ", id ! LUYU:: not "id" is id index NOT ID   

    ! For each coordinate, x,y,and z:
    DO n=1,nv4d
      
      ! Find the observations around this point.
      
      if (dodebug) WRITE(6,*) "Start calling obs_local..."
      if (dodebug) WRITE(6,*) "mean4d(i,it,:) = ", mean4d(i,it,:)
      if (dodebug) WRITE(6,*) "nobstotal = ", nobstotal
      CALL obs_local_drifters(id,it,var_local(n,:),mean4d(i,it,:),hdxf,rdiag,rloc,dep,nobsl,nobsl_dr,nobstotal)
      ! LUYU:  IN: id, it var_local(n,:),nobstotal
      ! LUYU: OUT: hdxf, rdiag, rloc, dep, nobsl

      parm = 1.0d0 !STEVE: keeping it simple
      if (dodebug) WRITE(6,*) "Finished calling obs_local_drifters..."
      if (dodebug) WRITE(6,*) "nobstotal=",nobstotal
      if (dodebug) WRITE(6,*) "nobsl_dr=",nobsl_dr 
      if (dodebug) WRITE(6,*) "hdxf=",hdxf
      if (dodebug) WRITE(6,*) "rdiag=",rdiag
      if (dodebug) WRITE(6,*) "rloc=",rloc
      if (dodebug) WRITE(6,*) "dep=",dep     

      CALL letkf_core(nobstotal,nobsl+nobsl_dr,hdxf,rdiag,rloc,dep,parm,trans_dr(:,:,nv3d+nv2d+n))
      if (dodebug) WRITE(6,*) "Finished calling letkf_core"
      if (dodebug) WRITE(6,*) "trans_dr(1,1,n)=", trans_dr(1,1,nv3d+nv2d+n)
      ! The ":" in place of "itim" implies that all times are affected equally by
      ! the observed drifter location.
      DO m=1,nbv
        anal4d(i,:,m,n) = mean4d(i,:,n)
        DO k=1,nbv
          anal4d(i,it,m,n) = anal4d(i,it,m,n) + gues4d(i,it,k,n) * trans_dr(k,m,nv3d+nv2d+n)
        END DO
      END DO

    ENDDO

   ENDDO

  DEALLOCATE(hdxf,rdiag,rloc,dep)

! If there are observations that weren't on the model grid, add them to the
! output list.

END SUBROUTINE das_drifters

END MODULE letkf_drifters
