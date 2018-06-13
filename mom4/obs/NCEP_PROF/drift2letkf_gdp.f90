PROGRAM drift2letkf_gdp

 !USE common_mom4 ! Use check subroutine
 USE netcdf
 USE common
 USE params_obs

 IMPLICIT NONE


  INTEGER :: i,j,k
  
  INTEGER,SAVE :: nd
  INTEGER,SAVE :: nf=14 !nf is number of fields


  REAL(r_size), ALLOCATABLE :: var(:,:)

  REAL(r_sngl) :: wk(7)

  CHARACTER (len=*), PARAMETER :: drffile="drifters_inp.txt"
  CHARACTER (len=*), PARAMETER :: outfile="obsin_drifters.dat"

  call read_gdp_dim(drffile)
  ALLOCATE(var(nd,nf))
  print *, "Finish counting the dimension."

  call read_gdp(drffile,var)
  print *, "Finish reading drifters_inp.txt."

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  DO i=1,nd
     wk(1)=id_x_obs
     wk(2)=var(i,4) ! dlon
     wk(3)=var(i,5) ! dlat
     wk(4)=15       ! dlev
     wk(5)=var(i,1) ! id
     wk(6)=var(i,8)/(3.1824/sqrt(3.00)) ! dlon_err
     wk(7)=var(i,3) ! hour
     WRITE(91) wk
     !print *, wk

     wk(1)=id_y_obs 
     wk(6)=var(i,9)/(3.1824/sqrt(3.00)) ! dlat_err
     WRITE(91) wk
     !print *, wk
  END DO

  CLOSE(91)

CONTAINS

SUBROUTINE read_gdp_dim(infile)
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(IN) :: infile
  INTEGER :: fid = 33
  INTEGER :: i,j,k,n,ios=0
  INTEGER :: dodebug=0, counter=1
  CHARACTER(slen) :: dummy

  OPEN(fid,FILE=infile,ACCESS='sequential')

  DO WHILE (ios==0)
    READ(fid,*,IOSTAT=ios) dummy
    !print *, counter  
    nd=counter
    counter=counter+1
  END DO
  nd = nd - 1

  CLOSE(fid)

  RETURN

END SUBROUTINE read_gdp_dim

SUBROUTINE read_gdp(infile,var)
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: var(nd,nf) 
  INTEGER :: fid = 32
  INTEGER :: i,j,k,n,id,jday,hour,drogue
  REAL(r_size) :: dlon,dlat,du,dv,dlon_err,dlat_err
  REAL(r_size) :: du_err,dv_err,gap,rmsgap
  INTEGER :: dodebug=0, counter=1
  CHARACTER(slen) :: dummy

  !print *, "Hello from read analysis data. Start reading ", infile

  OPEN(fid,FILE=infile,ACCESS='sequential')

  DO i=1,nd
    READ(fid,*) id,jday,hour,dlon,dlat,du,dv,dlon_err,dlat_err,du_err,dv_err,gap,rmsgap,drogue
    var(i,1) = id
    var(i,2) = jday
    var(i,3) = hour
    var(i,4) = dlon
    var(i,5) = dlat
    var(i,6) = du
    var(i,7) = dv
    var(i,8) = dlon_err
    var(i,9) = dlat_err
    var(i,10)= du_err
    var(i,11)= dv_err
    var(i,12)= gap
    var(i,13)= rmsgap
    var(i,14)= drogue
  END DO

  CLOSE(fid)

  RETURN

END SUBROUTINE read_gdp

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
  
END PROGRAM drift2letkf_gdp
