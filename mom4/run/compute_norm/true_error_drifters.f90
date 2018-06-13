PROGRAM true_error_drifters

  USE netcdf
  USE common
  USE params_model !nv4d

INTEGER,SAVE :: num_obdr, num_andr
INTEGER,SAVE :: num_max=100

! INPUT: obs_drifters_EEEE.txt and andr_me_EEEE.txt
  CHARACTER (15) :: mefile ="drif_me_EEEE.nc"  
  CHARACTER (16) :: mefile2="andr_me_EEEE.txt"
  CHARACTER (25) :: grdfile="grd.drifters_inp_EEEE.txt"
  CHARACTER (15) :: obfile="obs_dr_EEEE.txt"
  CHARACTER (len=*), PARAMETER :: outfile="true_error_drifters.txt"

! PARAMETER
  INTEGER :: pre_days=2
  INTEGER :: days=3 ! 365
  INTEGER :: nv4d_id
  REAL :: sigma=1
  INTEGER :: ncid,varid

! DATA STORAGE VAR
  INTEGER :: num_times,num_dr
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: ob_var, an_var, pos
  REAL(r_size),ALLOCATABLE,DIMENSION(:) :: ids 
  ! index   x   y   z  id 
  !     1   *   *   *  02
  !  ....
! DEBUG, VARIABLES
  INTEGER :: dodebug=0
  INTEGER :: i,j,k,n,fid=32
  REAL(r_size) :: error, norm, dist
  LOGICAL :: ex

  print *, "Hello from true_error_drifters..."

  nv4d_id = nv4d + 1 ! add another dimension for ids

  OPEN(fid,FILE=outfile,ACCESS='sequential')

  DO i=1,days

    norm =0

    ! STEP 1: Obtain data from observation
    WRITE(obfile(8:11),'(I4.4)') i
    WRITE(grdfile(18:21),'(I4.4)') i
    CALL read_dimension(grdfile)
    ALLOCATE(ob_var(num_obdr,nv4d_id))
    CALL read_obs(obfile,ob_var)

   ! if (dodebug) print *, "ob_var=", ob_var
   ! print *, "Finish STEP 1."


    ! STEP 2: Obain data from analysis
    if ( i .lt. pre_days )  then   
      WRITE(mefile(9:12),'(I4.4)') i
      call check(NF90_OPEN(mefile,NF90_NOWRITE,ncid))
      !print *, 'finish open the file: drifters_inp.nc.'
      ! Read Dim
      call check(NF90_INQ_DIMID(ncid,'np',varid)) !Number of drifters
      call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=num_andr))
      ALLOCATE(pos(nv4d,num_andr))
      ALLOCATE(ids(num_andr))
      ALLOCATE(an_var(num_andr,nv4d_id))

      ! Read Data
      call check(NF90_INQ_VARID(ncid,'positions',varid))
      call check(NF90_GET_VAR(ncid,varid,pos))
      call check(NF90_INQ_VARID(ncid,'ids',varid))
      call check(NF90_GET_VAR(ncid,varid,ids))
      !print *, 'finish setting up ids.'
      call check(NF90_CLOSE(ncid))  

      DO k=1,nv4d
        an_var(:,k)=pos(k,:)
      END DO
      an_var(:,k)=ids

      DEALLOCATE(pos,ids)
          
    else 
      WRITE(mefile2(9:12),'(I4.4)') i
      if ( i .eq. 1 ) then
        CALL read_andr_dim(mefile2)
      end if
      ALLOCATE(an_var(num_andr,nv4d_id))

      INQUIRE(FILE=mefile2,EXIST=ex)
      if (ex) then
        CALL read_andr(mefile2,an_var)
        if (dodebug) print *, "an_var=", an_var
      else
        print *, mefile2, "does not exist..."
        WRITE(fid,'(F16.6)') norm

        DEALLOCATE(an_var)
        DEALLOCATE(ob_var)
        cycle
      end if
      !print *, "Finish STEP 2." 
    end if  
 
    error=0
    ! STEP 3: Compute the error
    num_dr=num_obdr
    DO j=1,num_obdr
      ! LUYU: Need to modify later about matching id.

      ! Checking the whether the observation is inside the range
      if (ob_var(j,1) .LT. 0.166667 .AND. ob_var(j,1) .GT. 9.83333 ) then
        print *, "NO.",i,"drifter is out of the boundary."
        num_dr=num_dr-1
        cycle
      else
        if (ob_var(j,2) .LT. 15.1667 .AND. ob_var(j,2) .GT. 34.8333 ) then 
          print *, "NO.",i,"drifter is out of the boundary."
          num_dr=num_dr-1
          cycle
        else
          call com_distll_1(an_var(j,1),an_var(j,2),ob_var(j,1),ob_var(j,2),dist)
          error = error + dist
        end if
      end if
    END DO

    if (num_dr .ge. 0) then
      norm = SQRT(error / REAL(num_andr))/sigma
    else
      print *,"ERROR, NONE of the observation is reasonable" 
      norm = 0
      WRITE(fid,'(I16.6)') norm
      DEALLOCATE(an_var)
      DEALLOCATE(ob_var)     
      cycle
    end if

    !print *, error
    print *, norm

    WRITE(fid,'(I16.6)') norm

    DEALLOCATE(an_var)
    DEALLOCATE(ob_var)

  ENDDO

  CLOSE(fid)
  
 CONTAINS

SUBROUTINE read_andr_dim(infile)
  IMPLICIT NONE
  CHARACTER(16),INTENT(IN) :: infile
  INTEGER :: fid = 33
  INTEGER :: i,j,k,n,ios
  INTEGER :: dodebug=0, counter=1
  CHARACTER(slen) :: dummy

  OPEN(fid,FILE=infile,ACCESS='sequential')

  DO i=1,num_max
    READ(fid,*,IOSTAT=ios) dummy
    IF (ios /= 0) EXIT
    print *, i  
    num_andr=counter
    counter=counter+1
  END DO

  CLOSE(fid)

  RETURN

END SUBROUTINE read_andr_dim

SUBROUTINE read_andr(infile,var)
  IMPLICIT NONE
  CHARACTER(16),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: var(num_andr,nv4d_id) 
  INTEGER :: fid = 32
  INTEGER :: i,j,k,n,ios,id
  REAL(r_size) :: mlon, mlat, mlev
  INTEGER :: dodebug=0, counter=1
  CHARACTER(slen) :: dummy

  !print *, "Hello from read analysis data. Start reading ", infile

  OPEN(fid,FILE=infile,ACCESS='sequential')

  DO i=1,num_andr
    READ(fid, '(I16,3F16.6)') id,mlon,mlat,mlev
    var(i,4) = id
    var(i,1) = mlon
    var(i,2) = mlat
    var(i,3) = mlev
    if (dodebug) print *, id, mlon, mlat, mlev
  END DO

  CLOSE(fid)

END SUBROUTINE read_andr


SUBROUTINE read_dimension(infile)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  INTEGER :: num_times
  CHARACTER(16) :: dummy_char
  CHARACTER(8*16) :: header_line
  INTEGER :: fid = 33
  CHARACTER(slen) :: drfile
 
  !print *, 'Hello from read_dimension.'
  drfile = infile ! (DRIFTERS)

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_obdr, dummy_char, num_times
  !print *, 'num_drifters=', num_obdr, 'num_times=', num_times
  close(fid)

  RETURN
END SUBROUTINE read_dimension



SUBROUTINE read_obs(infile,var)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: var(num_obdr,nv4d_id)
  REAL :: dlon, dlat, ddepth, dtemp, dsalt, dtime
  INTEGER :: ditime, dids
  CHARACTER(8*16) :: header_line
  INTEGER :: di
  INTEGER :: fid = 33, dodebug = 0
  CHARACTER(24) :: drfile

  drfile = infile
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !print *, 'Hello from read_drifters.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,*) header_line
  read(fid,*) header_line

  ! Read all positions (and possibly temp and salt)  of each drifter:
  DO di=1,num_obdr
      read(fid,'(I16,6F16.6,I16)') dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      IF (dodebug .eq. 1) THEN
        print *, dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      END IF
      var(di,4) = dids
      var(di,1) = dlon
      var(di,2) = dlat
      var(di,3) = ddepth
      IF (nv4d .ge. 5) THEN
        ! If we have temperature and salinity observations at each position,
        ! we can assimilate this data too
        var(di,4) = dtemp
        var(di,5) = dsalt
      END IF
  END DO
  close(fid)

  RETURN

END SUBROUTINE read_obs

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

END PROGRAM true_error_drifters
