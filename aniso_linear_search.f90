MODULE aniso_linear_search
USE aniso_derivatives
USE mtz_io
USE refinement_util
IMPLICIT NONE
CONTAINS
    SUBROUTINE aniso_line_search(xold, np, fold, g, p, x, f, check,  mtz_1, nout)
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: xold      ! old params
        INTEGER,                     INTENT(IN)    :: np
        REAL(KIND=wp),               INTENT(IN)    :: fold
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: g
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: p         ! search direction
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: x         ! result
        REAL(KIND=wp),               INTENT(INOUT) :: f         ! optimal function value
        LOGICAL,                     INTENT(OUT)   :: check
        INTEGER,                     INTENT(IN)    :: nout
!       X-ray part:
        TYPE(mtz),                   INTENT(INOUT) :: mtz_1
!       Local variables:
        REAL(KIND=wp), DIMENSION(np,np)            :: AMAT
        INTEGER                                    :: n
        REAL(KIND=wp)                              :: stpmax
!        REAL(KIND=wp), PARAMETER                   :: ALF  = 1.0D-4
        REAL(KIND=wp), PARAMETER                   :: ALF  = 1.0D-2
        REAL(KIND=wp), PARAMETER                   :: TOLX = 1.0D-7
        REAL(KIND=wp)                              :: a
        REAL(KIND=wp)                              :: alam
        REAL(KIND=wp)                              :: alam2
        REAL(KIND=wp)                              :: alamin
        REAL(KIND=wp)                              :: b
        REAL(KIND=wp)                              :: disc
        REAL(KIND=wp)                              :: f2
        REAL(KIND=wp)                              :: rhs1
        REAL(KIND=wp)                              :: rhs2
        REAL(KIND=wp)                              :: slope
        REAL(KIND=wp)                              :: sum1
        REAL(KIND=wp)                              :: temp
        REAL(KIND=wp)                              :: test
        REAL(KIND=wp)                              :: tmplam
!       Counters:
        INTEGER                                    :: i
        CHARACTER(LEN=32),                   SAVE  :: srname = 'aniso_line_search'

!       Checks first:
        n = SIZE ( g )
        IF ( n /= SIZE ( p ) ) THEN
            WRITE(*,*) n, SIZE ( p )
            CALL die('Programming error. Inconsistent sizes of G and P arrays.', srname)
        ENDIF

!       Initialize check:
        check = .FALSE.
        sum1  = SQRT ( DOT_PRODUCT ( p, p ) )
        stpmax = 100.0_wp * MAX ( SQRT ( sum1 ), REAL(n, KIND=wp) )

        IF ( sum1 > stpmax ) THEN
            p = p * stpmax / sum1
        ENDIF

        g = -g
        slope = DOT_PRODUCT ( g, p )

        IF ( slope > 0.0_wp ) THEN
            WRITE(*,*) ' slope=', slope
            CALL die('Roundoff problem.', srname)
        ENDIF

!       Get params from current pdb:
!        CALL get_xold(xold, pdb_2)

        test = 0.0_wp
        DO i = 1, n
            temp = ABS ( p(i) ) / MAX ( ABS ( xold(i) ), 1.0_wp )
            IF ( temp > test ) test = temp
        ENDDO

        alamin = TOLX / test
        alam   = 1.0_wp

        DO WHILE (.TRUE.)
            x = xold + alam * p
            CALL bulk_aniso_deriv(x, np, mtz_1, f, .TRUE., AMAT, g)
            IF ( nout > 0 ) THEN
                WRITE(*,"(' ANISO_LINE_SEARCH> ', 'f=', ES11.4, ' fold=', ES11.4, ' alam=', F9.5)") f, fold, alam
            ENDIF
            IF ( alam < alamin ) THEN
                IF ( nout > 0 ) THEN
                    WRITE(*,*) alam, alamin
                    CALL warn('Zero shift has been applied. Refinement has converged...', srname)
                ENDIF
                x = xold
                check = .TRUE.
                RETURN

            ELSE IF ( f < fold + ALF * alam * slope ) THEN
                IF ( nout > 0 ) THEN
                    CALL messag('Sufficient function decrease. Line search stopped.', srname)
                ENDIF
                RETURN
            ELSE

                IF ( alam == 1.0_wp ) THEN

                    tmplam = -slope / ( 2.0_wp * ( f - fold - slope ) )

                ELSE

                    rhs1 = f  - fold - alam  * slope
                    rhs2 = f2 - fold - alam2 * slope
                    a    = (rhs1 / alam ** 2 - rhs2 / alam2 ** 2) / (alam - alam2)
                    b    = (-alam2 * rhs1 / alam ** 2 + alam * rhs2 / alam2 ** 2) / (alam - alam2)

                    IF ( a == 0.0_wp ) THEN

                        tmplam = -slope / (2.0 * b)

                    ELSE
                        disc = b * b - 3.0_wp * a * slope
 
                        IF ( disc < 0.0_wp ) THEN

                            tmplam = 0.5_wp * alam

                        ELSE IF ( b <= 0.0_wp ) THEN

                            tmplam = ( -b + SQRT ( disc ) ) / (3.0_wp * a)

                        ELSE

                            tmplam = -slope / ( b + SQRT ( disc ) )

                        ENDIF

                    ENDIF

                   IF ( tmplam > 0.5_wp * alam ) tmplam = 0.5_wp * alam

                ENDIF

            ENDIF

            alam2 = alam
            f2    = f
            alam  = MAX ( tmplam, 0.1_wp * alam )

        ENDDO
    END SUBROUTINE aniso_line_search
        
END MODULE aniso_linear_search
