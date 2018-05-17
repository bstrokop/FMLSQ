MODULE fft_util
!
!   Purpose:
!   =======
!   Small utility routines basically for FFT calculations
!
!   Date:         Programmer:        Description of changes:
!   ====          ==========         ======================
!   Oct 2003      B.Strokopytov      Original code
!   Oct 2005      B.Strokopytov      fftfac routine added
!   Nov 2005      B.Strokopytov      phase and phase_of made ELEMENTAL
!
USE select_kinds
USE vectors 
IMPLICIT NONE
CONTAINS
    FUNCTION hemi_l(hkl)
        LOGICAL                      :: hemi_l
        TYPE(vector_int), INTENT(IN) :: hkl
!       Local vars:
        INTEGER, DIMENSION(3)        :: v
!
        v = hkl
        hemi_l = ( v(3) > 0 .OR. ( v(3) == 0 .AND. ( v(1) > 0 .OR. ( v(1) == 0 .AND. v(2) > 0 ) ) ) )
    END FUNCTION hemi_l

    FUNCTION hemi_h(hkl)
        LOGICAL                      :: hemi_h
        TYPE(vector_int), INTENT(IN) :: hkl
!       Local vars:
        INTEGER, DIMENSION(3)        :: v

        v = hkl
        hemi_h = ( v(1) > 0 .OR. ( v(1) == 0 .AND. ( v(2) > 0 .OR. ( v(2) == 0 .AND. v(3) > 0 ) ) ) )
    END FUNCTION hemi_h

    ELEMENTAL FUNCTION phase ( phi )
        COMPLEX(KIND=wp)          :: phase
        REAL(KIND=wp), INTENT(IN) :: phi
        phase = CMPLX(COS(phi), SIN(phi), KIND=wp)
    END FUNCTION phase

    ELEMENTAL FUNCTION phase_of ( f )
        REAL(KIND=wp)                :: phase_of
        COMPLEX(KIND=wp), INTENT(IN) :: f
        phase_of = ATAN2(AIMAG(f), REAL(f, KIND=wp)) 
    END FUNCTION phase_of

    FUNCTION sphere_volume ( radius )
!
!       Date:          Programmer:        History of changes:
!       ====           ==========         ==================
!       Oct 2004       Strokopytov B.     Original code
!
        REAL(KIND=wp) :: sphere_volume
        REAL(KIND=wp) :: radius
        sphere_volume = (4.0_wp * ACOS(-1.0_wp) / 3.0_wp) * radius ** 3 
    END FUNCTION sphere_volume

    FUNCTION fftfac(n0, m0)
!
!       Purpose:
!       =======
!       Find the composite number greater n than n0 which
!       will optimise the speed of the fft and has m0 as factor
!       (naive calculation: TIME = (p1.p2.p3...pn).(p1+p2+p3+..pn) )
!
!       Lifted from K.Cowtan `dm' package. Adapted for F95 by B.Strokopytov   Oct 2003
!
        INTEGER               :: fftfac
        INTEGER, INTENT(IN)   :: n0
        INTEGER, INTENT(IN)   :: m0
!       Don't want to use prime factor 13:
        INTEGER, PARAMETER    :: nprime = 5
        INTEGER, DIMENSION(6) :: iprime = (/2,3,5,7,11,13/)
!
        INTEGER               :: i
        INTEGER               :: j
        INTEGER               :: n
        INTEGER               :: n1
        INTEGER               :: nmin
        INTEGER               :: time
        INTEGER               :: timemin
        INTEGER               :: sum
        INTEGER               :: prod
!
        n1      = ( ( n0 + m0 - 1 ) / m0 ) * m0
        nmin    = 0
        timemin = 1000000

!       Try a number:
        DO n = n1, 2*n1, m0
!
!           Decompose n and add its factors to sum and prod
            i    = n
            sum  = 0
            prod = 1

            DO j = 1, nprime
                DO WHILE ( MOD ( i, iprime(j) ) == 0 )
                    i    = i / iprime(j)
                    sum  = sum  + iprime(j)
                    prod = prod * iprime(j)
                ENDDO
            ENDDO
            time = sum * prod
            IF ( i == 1 .AND. time < timemin ) THEN
                 nmin = n
                 timemin = time
            ENDIF
        ENDDO
!
        fftfac = nmin
!
    END FUNCTION fftfac

END MODULE fft_util
