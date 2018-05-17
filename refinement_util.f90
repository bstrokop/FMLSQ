MODULE refinement_util
USE constants
USE fail
USE util
IMPLICIT NONE
CONTAINS
    SUBROUTINE which_params(mode, refine_xyz, refine_Biso, refine_occ, refine_Us)
        CHARACTER(LEN=*), INTENT(IN)              :: mode
        LOGICAL,          INTENT(INOUT)           :: refine_xyz
        LOGICAL,          INTENT(INOUT)           :: refine_Biso
        LOGICAL,          INTENT(INOUT)           :: refine_occ
        LOGICAL,          INTENT(INOUT)           :: refine_Us
!       Local variables:
        CHARACTER(LEN=32)                         :: srname

        srname = 'which_params'

        refine_xyz  = INDEX ( mode, 'X' ) > 0 .OR. INDEX ( mode, 'Y' ) > 0 .OR. INDEX ( mode, 'Z' ) > 0
        refine_Biso = INDEX ( mode, 'B' ) > 0
        refine_occ  = INDEX ( mode, 'Q' ) > 0 .OR. INDEX ( mode, 'O' ) > 0
        refine_Us   = INDEX ( mode, 'U' ) > 0

        IF ( .NOT. refine_xyz .AND. .NOT. refine_biso .AND. .NOT. refine_occ .AND. .NOT. refine_Us) THEN
            WRITE(*,*) ' mode= ', mode
            CALL die('Unreasonable mode. Nothing to refine...', srname)
        ENDIF

    END SUBROUTINE which_params

    SUBROUTINE cubic_line_search ( p, lambda1, lambda2, f1, f2, g1, g2, lambda_opt )
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: p
        REAL(KIND=wp),               INTENT(IN)    :: lambda1
        REAL(KIND=wp),               INTENT(IN)    :: lambda2
        REAL(KIND=wp),               INTENT(IN)    :: f1
        REAL(KIND=wp),               INTENT(IN)    :: f2
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: g1
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: g2
        REAL(KIND=wp),               INTENT(INOUT) :: lambda_opt
!       Local variables:
        REAL(KIND=wp)                              :: g1_slope
        REAL(KIND=wp)                              :: g2_slope
        REAL(KIND=wp)                              :: lambda_diff
!       Cubic polynomial:
        REAL(KIND=wp)                              :: a
        REAL(KIND=wp)                              :: b 
        REAL(KIND=wp)                              :: c
        REAL(KIND=wp)                              :: d
        CHARACTER(LEN=32),                    SAVE :: srname ='cubic_line_search'

        g1_slope = DOT_PRODUCT ( g1,  p)
        g2_slope = DOT_PRODUCT ( g2,  p)
        
        lambda_diff = lambda2 - lambda1

!       Cubic polynomial:
        a = ( -2.0_wp * (f2 - f1) + (         g1_slope + g2_slope) * lambda_diff ) / lambda_diff ** 3
        b = (  3.0_wp * (f2 - f1) - (2.0_wp * g1_slope + g2_slope) * lambda_diff ) / lambda_diff ** 2
        c = g1_slope
        d = f1

!       Obtain minimum:
        IF ( a /= 0.0_wp .AND. b ** 2 - 3.0_wp * a * c >  0.0_wp ) THEN
            lambda_opt = lambda1 + (-b + SQRT ( b ** 2 - 3.0_wp * a * c )) / ( 3.0_wp * a )
            CALL messag('Cubic polynomial has been used for calculation of optimal shifts...', srname)
        ELSE
            lambda_opt = lambda1 - c / ( 2.0_wp * b )
            CALL messag('Quadratic polynomial has been used for calculation of optimal shifts...', srname)
        ENDIF

    END SUBROUTINE cubic_line_search

    FUNCTION estimate_integration_radius ( b, xaccuracy )
!
!       Purpose:
!       =======
!       Estimates integration radius for Tronrud method. We basically use it
!       to calculate normal matrix diagonal to improve condition number of normal
!       matrix.
!
!       Date:                 Programmer:           History of changes:
!       ====                  ==========            ==================
!       November 14, 2007     B.Strokopytov         Original code
!       April 14, 2007        B.Strokopytov         Important addition to code has been made
!                                                   see fig.1
!                                                   Otherwise e.g. b=98.696 will suddenly
!                                                   give us abs(f) < 1.0E-5 for r=2.5 A
!
!       f(r)
!       |**
!       |  *    sign change 
!       |    *         monotonicity region             r
!       -------*--------------------------*----------->
!                *             *    *   
!                  *     *  *
!                      *  
!                      fig.1
!
!
        REAL(KIND=wp)                        :: estimate_integration_radius
        REAL(KIND=wp), INTENT(IN)            :: b
        REAL(KIND=wp), INTENT(IN), OPTIONAL  :: xaccuracy
!       Local variables:
        INTEGER, PARAMETER                   :: MAXRAD = 50
        REAL(KIND=wp)                        :: accuracy
        REAL(KIND=wp)                        :: f 
        REAL(KIND=wp)                        :: f1 
        REAL(KIND=wp)                        :: r
        INTEGER                              :: biso_start
!       Counters: 
        INTEGER                              :: i
        CHARACTER(LEN=32)                    :: srname = 'estimate_integration_radius'

        IF ( PRESENT ( xaccuracy ) ) THEN
            accuracy = xaccuracy
        ELSE
            accuracy = 0.1_wp ** 5 
        ENDIF
 
        IF ( b < 0.1_wp ** 6) THEN
            WRITE(*,*) ' b=', b
            CALL die('Programming error. Non-positive value of B on input.', srname)
        ENDIF

        biso_start = 2 * NINT ( SQRT ( 5.0_wp * b / 8.0_wp ) / pi  ) + 1

!       Due to this algorithm minimal integration radius is (0.5 * biso_start) A:
        DO i = biso_start, 2 * MAXRAD

            r = i / 2.0_wp

            f = -pi * r ** 3 * ( 5.0_wp - 8.0_wp * pi ** 2 * r ** 2 / b ) &
              * EXP ( -4.0_wp * pi ** 2 * r ** 2 / b )

!           To avoid problems when f changes sign we MUST check:
            r = r + 0.5_wp

            f1 = -pi * r ** 3 * ( 5.0_wp - 8.0_wp * pi ** 2 * r ** 2 / b ) &
              * EXP ( -4.0_wp * pi ** 2 * r ** 2 / b ) 

!           If function value is small and monotonically decreasing only then we may exit:
            IF ( ABS ( f ) < accuracy .AND. ABS ( f1 ) < ABS ( f ) ) EXIT
        ENDDO

        IF ( r > 50 ) THEN
            WRITE(*,*) ' b=', b, ' r=', r, ' f=', f
            CALL die('Abnormally high B-value on input...', srname)
        ELSE
            estimate_integration_radius = r - 0.5_wp
        ENDIF

    END FUNCTION estimate_integration_radius

    FUNCTION integration_radius ( a_gauss, b_gauss, number_of_electrons_left_out )
!
!       Purpose:
!       =======
!       Estimates integration radius for Agarwal method. It's important
!       to have common routine for electron densitry calculations and
!       convolution.
!      
!
!       Date:                 Programmer:           History of changes:
!       ====                  ==========            ==================
!       November 14, 2007     B.Strokopytov         Original code
!       March 31, 2009        B.Strokopytov         Precise calculation over all gaussians.
!
        REAL(KIND=wp)                                      :: integration_radius
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: a_gauss
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: b_gauss                     
!       Number of electrons left out:
        REAL(KIND=wp),               INTENT(IN), OPTIONAL  :: number_of_electrons_left_out
        REAL(KIND=wp)                                      :: C1 
        REAL(KIND=wp)                                      :: b
!       Local variables:
        INTEGER, PARAMETER                                 :: MAXRAD = 50
        REAL(KIND=wp)                                      :: accuracy
        REAL(KIND=wp)                                      :: f
        REAL(KIND=wp)                                      :: r
!       Counters: 
        INTEGER                                            :: i
        INTEGER                                            :: l
        CHARACTER(LEN=32)                                  :: srname = 'integration_radius'

        IF ( PRESENT ( number_of_electrons_left_out ) ) THEN
            accuracy =  number_of_electrons_left_out
        ELSE
            accuracy = default_genden_accuracy
        ENDIF




!       This will limit max radius to 50:
        DO i = 1, 2 * MAXRAD
            r = i / 2.0_wp
            f = 0.0_wp
            DO l = 1, SIZE ( a_gauss )
                C1 = a_gauss(l)
                b  = b_gauss(l)

                IF ( b < 0.1_wp ** 5 ) THEN
                    WRITE(*,*) ' b=', b
                    CALL die('Programming error. Non-positive value of B on input.', srname)
                ENDIF

                f = f + ( C1 *  b / (8.0_wp * pi ** 2) &
                      * ( SQRT ( pi * b ) * ( 1.0_wp - ERF ( 2.0_wp * pi * r / SQRT ( b ) ) ) &
                      + 4.0_wp * pi * r * EXP ( -4 * pi ** 2 * r ** 2 / b ) ) )

            ENDDO

!           Suitable radius has been found:         
            IF ( ABS ( f ) < accuracy ) EXIT
           
        ENDDO

        IF ( r > 49.9_wp ) THEN
            WRITE(*,"(' b=', F8.2, ' r=', F8.1, ' f=', ES9.2)") b, r, f
            CALL warn('Abnormally high B-value on input...', srname)
        ELSE
            integration_radius = r
        ENDIF

    END FUNCTION integration_radius

    SUBROUTINE convert_gg ( g, refine, gg, indxyz )
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: g
        LOGICAL,       DIMENSION(:), INTENT(IN)    :: refine
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: gg
        INTEGER,                     INTENT(IN)    :: indxyz
!       Local variables:
        INTEGER                                    :: l
        INTEGER                                    :: k

        l = 0
        DO k = 1, 3
            IF ( refine(k) ) THEN
                l = l + 1
                g(indxyz+l) = g(indxyz+l) + gg(k)
            ENDIF
        ENDDO

    END SUBROUTINE convert_gg

    SUBROUTINE get_num_threads(nthreads)
        INTEGER, INTENT(OUT):: nthreads
        INTEGER, EXTERNAL   :: omp_get_num_threads
!$OMP PARALLEL
!$OMP MASTER
      PRINT *,'----------------------------------------------'
      nthreads = OMP_GET_NUM_THREADS()
!$    PRINT *,'Number of Threads = ', nthreads
!$
!$OMP END MASTER
!$OMP END PARALLEL
    END SUBROUTINE get_num_threads

END MODULE refinement_util
