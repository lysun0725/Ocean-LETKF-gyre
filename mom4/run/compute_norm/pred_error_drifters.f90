PROGRAM pred_error_drifters

  USE common
  USE params_model !nv4d

INTEGER,SAVE :: num_andr
INTEGER,SAVE :: num_max=100

! INPUT: obs_drifters_EEEE.txt and andr_me_EEEE.txt
  CHARACTER (16) :: spfile="andr_sp_EEEE.txt"
  CHARACTER (len=*), PARAMETER :: outfile="pred_error_drifters.txt"

! PARAMETER
  INTEGER :: days=31 ! 365
  INTEGER :: nv4d_id
  INTEGER :: sigma=1


! DATA STORAGE VAR
  INTEGER :: num_times
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: an_var 
  ! index   x   y   z  id 
  !     1   *   *   *  02
  !  ....
! DEBUG, VARIABLES
  INTEGER :: dodebug=1
  INTEGER :: i,j,k,n,fid=32
  REAL(r_size) :: error, norm

  print *, "Hello from true_error_drifters..."

  nv4d_id = nv4d + 1 ! add another dimension for ids

  OPEN(fid,FILE=outfile)

  DO i=1,days

    WRITE(spfile(9:12),'(I4.4)') i

    ! STEP 1: Obain data from analysis
    if (i .eq. 1) then
      CALL read_andr_dim(spfile)
    end if
    ALLOCATE(an_var(num_andr,nv4d_id))
    CALL read_andr(spfile,an_var)
    
    error=0
    ! STEP 3: Compute the error
    DO j=1,num_andr
        ! LUYU: Need to modify later about matching id.
        error = error + (an_var(j,1))**2 + (an_var(j,2))**2
    END DO

    norm = SQRT(error / num_andr)

    WRITE(fid,'(I16.6)') norm
    print *, norm
    DEALLOCATE(an_var)

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
  INTEGER :: fid = 33
  INTEGER :: i,j,k,n,ios,id
  REAL(r_size) :: mlon, mlat, mlev
  INTEGER :: dodebug=0, counter=1
  CHARACTER(slen) :: dummy

  OPEN(fid,FILE=infile,ACCESS='sequential')

  DO i=1,num_andr
    READ(fid,'(I16,3F16.6)') id,mlon,mlat,mlev
    var(i,4) = id
    var(i,1) = mlon
    var(i,2) = mlat
    var(i,3) = mlev
    if (dodebug) print *, var(i,4), var(i,1:3)
  END DO

  CLOSE(fid)

END SUBROUTINE read_andr


END PROGRAM pred_error_drifters
