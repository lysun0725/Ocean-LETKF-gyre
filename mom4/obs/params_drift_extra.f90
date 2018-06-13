MODULE params_drift_extra

USE common, ONLY: r_size
USE params_model, ONLY: nv4d

IMPLICIT NONE

PUBLIC

! EXTRA obs id for temperature and salinity
INTEGER, PARAMETER :: id_t_dobs=4444 !(OCEAN) (DRIFTERS)
INTEGER, PARAMETER :: id_s_dobs=5555 !(OCENA) (DRIFTERS)

! EXTRA model id for temperature and salinity on drifters
INTEGER,PARAMETER :: iv4d_t=4                !(OCEAN) (DRIFTERS)
INTEGER,PARAMETER :: iv4d_s=5                !(OCEAN) (DRIFTERS)

END MODULE params_drift_extra
