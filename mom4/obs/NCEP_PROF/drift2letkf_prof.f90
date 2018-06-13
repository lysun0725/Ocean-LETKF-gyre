PROGRAM drift2letkf_ts

 USE common
 USE params_obs

 IMPLICIT NONE

  INTEGER :: i,j,k 
  INTEGER,SAVE :: nids, ntimes 
  REAL(r_size), ALLOCATABLE ::var(:,:)
  REAL(r_sngl) :: wk(6)
  CHARACTER (len=*), PARAMETER :: drffile = "drifters_inpc.txt"
  CHARACTER (len=*), PARAMETER :: outfile = "obsin.dat"

  CALL read_dimension(drffile,nids,ntimes)
  ALLOCATE(var(nids,nid_dobs))
  CALL read_drifters(drffile,var)

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  DO i=1,nids
    wk(1) = id_t_obs ! elem
    wk(2) = var(i,1) ! rlon
    wk(3) = var(i,2) ! rlat
    wk(4) = var(i,3) ! rlev
    wk(5) = var(i,4)
    wk(6) = 0.01
    WRITE(91) wk
    PRINT *, wk


    wk(1) = id_s_obs
    wk(6) = 1.e-5
    WRITE(91) wk
    PRINT *, wk 
    
  END DO
  CLOSE(91)

CONTAINS

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

  !ALLOCATE(drifter_ids(num_drifters))

  RETURN
END SUBROUTINE read_dimension

SUBROUTINE read_drifters(infile,v4d_all)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile 
  REAL(r_size),INTENT(OUT) :: v4d_all(nids,nid_dobs)
  REAL :: dlon, dlat, ddepth, dtemp, dsalt, dtime
  INTEGER :: ditime, dids
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: di
  INTEGER :: dummy_num
  INTEGER :: fid = 33, dodebug = 0
  CHARACTER(24) :: drfile

  drfile = infile
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *, 'Hello from read_drifters.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, dummy_num, dummy_char, dummy_num
  read(fid,*) header_line

  ! Read all positions (and possibly temp and salt)  of each drifter:
  DO di=1,nids
      read(fid,'(I16,4F16.6,E16.4,F16.6,I16)') dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      IF (dodebug .eq. 1) THEN
        print *, dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      END IF
      !drifter_ids(di) = dids
      v4d_all(di,1) = dlon
      v4d_all(di,2) = dlat
      v4d_all(di,3) = ddepth
      IF (nid_dobs .ge. 5) THEN
        ! If we have temperature and salinity observations at each position,
        ! we can assimilate this data too
        v4d_all(di,4) = dtemp
        v4d_all(di,5) = dsalt
      END IF
  END DO
  close(fid)

END SUBROUTINE read_drifters


END PROGRAM drift2letkf_ts
