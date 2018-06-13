PROGRAM drif2prof2letkf

USE netcdf
USE common
USE params_obs

  REAL(r_size), ALLOCATABLE, DIMENSION(:,:), SAVE :: v4d 
  REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_ids 
  INTEGER, SAVE :: num_drifters 
  INTEGER, SAVE :: num_times
  INTEGER, SAVE :: nv4d=5

  CHARACTER (len=*), PARAMETER :: infile = "drifters_inpc.txt"
  CHARACTER (len=*), PARAMETER :: outfile="obsin.dat"
  REAL, DIMENSION(:,:), ALLOCATABLE :: positions
  INTEGER, DIMENSION(:), ALLOCATABLE :: ids
  INTEGER :: nd, np
  REAL(r_sngl) :: wk(6)

  call read_dimension(infile,num_drifters,num_times)
  nd = nv4d !LUYU: most likely this is fixed. And it means we only read the position of each drifter.
  np = num_drifters

  ALLOCATE(v4d(np,nd))

  call read_drifters(infile,v4d)

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')  
  DO i = 1,np
    wk(1) = id_t_obs
    wk(2) = v4d(i,1)   ! lon
    wk(3) = v4d(i,2)   ! lat
    wk(4) = v4d(i,3)   ! lev
    wk(5) = v4d(i,4)   ! temp
    wk(6) = 0.5        ! err
    print *, wk
    WRITE(91) wk
    wk(1) = id_s_obs
    wk(5) = v4d(i,5)   ! salt
    wk(6) = 0.05       ! err 
    WRITE(91) wk
    print *, wk
  END DO
  CLOSE(91)
CONTAINS

subroutine check(status)
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

SUBROUTINE read_dimension(infile,num_drifters,num_times)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  INTEGER, INTENT(OUT) :: num_drifters,num_times
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: fid = 33
  CHARACTER(slen) :: drfile
 
  print *, 'Hello from read_dimension.'
  drfile = infile ! (DRIFTERS)

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  print *, 'num_drifters=', num_drifters, 'num_times=', num_times
  close(fid)

  ALLOCATE(drifter_ids(num_drifters))

  RETURN
END SUBROUTINE read_dimension

SUBROUTINE read_drifters(infile,v4d_all)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile 
  REAL(r_size),INTENT(OUT) :: v4d_all(num_drifters,nv4d)
  REAL :: dlon, dlat, ddepth, dtemp, dsalt, dtime
  INTEGER :: ditime, dids
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: di
  INTEGER :: fid = 33, dodebug = 0
  CHARACTER(24) :: drfile

  drfile = infile
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *, 'Hello from read_drifters.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  read(fid,*) header_line

  ! Read all positions (and possibly temp and salt)  of each drifter:
  DO di=1,num_drifters
      read(fid,'(I12,6F12.4,I12)') dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      IF (dodebug .eq. 1) THEN
        print *, dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      END IF
      drifter_ids(di) = dids
      v4d_all(di,1) = dlon
      v4d_all(di,2) = dlat
      v4d_all(di,3) = ddepth
      IF (nv4d .ge. 5) THEN
        ! If we have temperature and salinity observations at each position,
        ! we can assimilate this data too
        v4d_all(di,4) = dtemp
        v4d_all(di,5) = dsalt
      END IF
  END DO
  close(fid)

END SUBROUTINE read_drifters

END PROGRAM drif2prof2letkf
