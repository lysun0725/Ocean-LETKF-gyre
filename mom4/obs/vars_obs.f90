MODULE vars_obs

USE common, ONLY: r_size, r_sngl

IMPLICIT NONE

PUBLIC

REAL(r_size),SAVE :: dist_zero
REAL(r_size),SAVE :: dist_zerov
REAL(r_size),ALLOCATABLE,SAVE :: dlon_zero(:)
REAL(r_size),SAVE :: dlat_zero
REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshour(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)

INTEGER, ALLOCATABLE, SAVE :: obs_useidx(:) !STEVE: general version of nobs_use()
INTEGER, ALLOCATABLE, SAVE :: nobsgrd(:,:) !(nlon,nlat)

!LUYU: for (DRIFTERS)
REAL(r_size),ALLOCATABLE,SAVE :: obselm_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslon_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslat_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslev_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsid(:)
REAL(r_size),ALLOCATABLE,SAVE :: obserr_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsdep_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshour_dr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshdxf_dr(:,:)
REAL(r_size),ALLOCATABLE,SAVE :: obstime(:)
REAL(r_sngl),ALLOCATABLE,SAVE :: obsval(:)
REAL(r_size),ALLOCATABLE,SAVE :: nobsgrd_dr(:,:)




END MODULE vars_obs
