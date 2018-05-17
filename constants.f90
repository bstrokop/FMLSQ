MODULE constants
!
!     Purpose:
!     =======
!     Contains various assorted standard constants
!     
USE select_kinds
IMPLICIT NONE
! Potential 2-gaussian_atoms:
CHARACTER(LEN=4),   PARAMETER :: two_gaussian_atoms = 'HCNOS'
! Need this to interpret correctly scattering factor for common atoms:
CHARACTER(LEN=6),   PARAMETER :: common_atoms = 'HCNOPS'
INTEGER,            PARAMETER :: tsym_factor = 24
REAL(KIND=wp),      PARAMETER :: eps         = 0.1_wp ** 6 
REAL(KIND=wp),      PARAMETER :: pi          = 3.14159265358979_wp
REAL(KIND=wp),      PARAMETER :: twopi       = 2.0_wp * pi
REAL(KIND=wp),      PARAMETER :: small       = 0.1_wp ** 6 
REAL(KIND=wp),      PARAMETER :: pi4         = 4.0_wp * pi
REAL(KIND=wp),      PARAMETER :: pi_squared  = pi ** 2
REAL(KIND=wp),      PARAMETER :: pi_squared_8  = 8.0_wp * pi_squared
REAL(KIND=wp),      PARAMETER :: aniso_pdb_factor = 10.0_wp ** 4
REAL(KIND=wp),      PARAMETER :: default_genden_accuracy = 0.1_wp ** 7 
! Max length:
INTEGER,            PARAMETER :: wlen = 256    ! word                 
INTEGER,            PARAMETER :: llen = 400    ! line
INTEGER,            PARAMETER :: flen = 512    ! file name
! Matrix
INTEGER,            PARAMETER :: max_block_size = 11
! Xray:
LOGICAL                       :: xray          ! want to use xray data?
LOGICAL                       :: anom          ! with anomalous scattering?
! Debug:
INTEGER                       :: debug         ! for debugging
!  Save non-PARAMETER global variables:
REAL(KIND=wp)                 :: lambda_shift
!REAL(KIND=wp)                 :: MIN_ADP_RATIO = 0.15_wp
REAL(KIND=wp)                 :: MIN_ADP_RATIO = 0.05_wp
REAL(KIND=wp)                 :: max_abs_xyz_shift = 0.5_wp
REAL(KIND=wp), PARAMETER      :: max_abs_biso_shift = 20.0_wp
!REAL(KIND=wp)                 :: max_abs_U_shift = max_abs_biso_shift / pi_squared_8
SAVE
END MODULE constants
