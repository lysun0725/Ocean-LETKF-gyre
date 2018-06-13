PROGRAM setup_ini_drifters

  USE netcdf
  USE common


  REAL(r_size), ALLOCATABLE, DIMENSION(:,:), SAVE :: v4d 
  INTEGER, SAVE :: num_drifters=50 !set up number of drifters to generate 
  INTEGER, SAVE :: num_times
  INTEGER, SAVE :: nv4d=3
 
  CHARACTER(32) :: drfile = "drifters_inp.nc"
  INTEGER :: nd_dimid, np_dimid, dimids(2)
  INTEGER :: pos_varid, ids_varid, i
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: positions
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: rlon, rlat
  INTEGER, DIMENSION(:), ALLOCATABLE :: ids
  ! For netcdf:
  INTEGER :: ncid, nd, np, dodebug=0

  if (dodebug) print *, "Hello from setup_ini_drifters..."

  nd = nv4d !LUYU: most likely this is fixed. And it means we only read the position of each drifter.
  np = num_drifters

  if (dodebug) print *, "Finish setting up the dimension..."

  ALLOCATE(rlon(np))
  ALLOCATE(rlat(np))
  ALLOCATE(v4d(np,nd))
  ALLOCATE(positions(nd,np))
  ALLOCATE(ids(np))

  call com_rand(np,rlon)
  if (dodebug) print *, "rlon=", rlon
  rlon = 2.5 + rlon*2.5
  v4d(:,1)=rlon
  if (dodebug) print *, "v4d(:,1)=", v4d(:,1)

  call com_rand(np,rlat)
  rlat = 20 + rlat*7
  v4d(:,2)=rlat

  v4d(:,3)=15

  if (dodebug) print *, "Finish setting up the positions..."

  !WRITE(6,*) "#### BEFORE####"
  !WRITE(6,*) v4d_all(1:np,num_times,1:nd)
  positions = RESHAPE(TRANSPOSE(v4d),(/nd,np/)) 
  !WRITE(6,*) "#### AFTER####"
  !WRITE(6,*) positions
  
  DO i = 1, num_drifters
    ids(i)=i
  END DO

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  call check( nf90_create(drfile, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "nd", nd, nd_dimid) )
  call check( nf90_def_dim(ncid, "np", np, np_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ nd_dimid, np_dimid /)

  ! Define the variable.
  call check( nf90_def_var(ncid, "positions", NF90_DOUBLE, dimids, pos_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, pos_varid, "names", "lon lat depth") )
  call check( nf90_put_att(ncid, pos_varid, "units", "deg_E deg_N meters") )

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).
  call check( nf90_def_var(ncid, "ids", NF90_INT, np_dimid, ids_varid) )

  ! Add global attributes with NF90_GLOBAL
  call check( nf90_put_att(ncid, NF90_GLOBAL, "velocity_names", "u v w") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_names", "lon lat depth temp salt") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_units", "deg_E deg_N meters Celsius PSU") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "time_units", "seconds") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "title", "LETKF analyzed positions for drifters, for input into MOM4p1") )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the data to the file.
  call check( nf90_put_var(ncid, pos_varid, positions) )
  call check( nf90_put_var(ncid, ids_varid, ids) )

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

CONTAINS

subroutine check(status)
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

END PROGRAM setup_ini_drifters
