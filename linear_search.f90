MODULE linear_search
USE basic_pdb_manip
USE first_derivatives
USE refinement_util
IMPLICIT NONE
CONTAINS
    SUBROUTINE book_line_search(xold, fold, g, p, x, f, check,  mtz_1, map_1, binary_map, logical_map, pdb_2, mode)
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: xold      ! old params
        REAL(KIND=wp),               INTENT(IN)    :: fold
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: g
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: p         ! search direction
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: x         ! result
        REAL(KIND=wp), INTENT(INOUT)               :: f         ! optimal function value
        LOGICAL,       INTENT(OUT)                 :: check
!       X-ray part:
        TYPE(mtz),                   INTENT(INOUT) :: mtz_1
        TYPE(map),                   INTENT(INOUT) :: map_1
!       Bulk solvent:
        TYPE(map),                   INTENT(INOUT) :: binary_map
        TYPE(map_logical),           INTENT(INOUT) :: logical_map
        TYPE(pdb),                   INTENT(INOUT) :: pdb_2
        CHARACTER(LEN=*)                           :: mode
!       Local variables:
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
        CHARACTER(LEN=32),                   SAVE  :: srname = 'book_line_search'

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
        CALL get_xold(xold, pdb_2)

        test = 0.0_wp
        DO i = 1, n
            temp = ABS ( p(i) ) / MAX ( ABS ( xold(i) ), 1.0_wp )
            IF ( temp > test ) test = temp
        ENDDO

        alamin = TOLX / test
        alam   = 1.0_wp

        DO WHILE (.TRUE.)
            x = xold + alam * p
            CALL apply_shifts (p, pdb_2, alam, final=.TRUE.)

!           f = func(x,n), just function call (calling slow old genden):

!           Important change. No images will computed if only_fun=.TRUE.:
            CALL new_lsq_first_deriv_wrt_all ( g, mtz_1, map_1, binary_map, logical_map, pdb_2, &
                                               mode, fun=f, only_fun=.TRUE.)

!           Restore params in PDB:
            CALL set_xold(xold, pdb_2)
!            CALL corr_coef_func_deriv_wrt_6p(x, dccdx, mtz_1, map_1, rotated_pdb, polar_axis, f=f)
            WRITE(*,"(' BOOK_LINE_SEARCH> ', 'f=', ES11.4, ' fold=', ES11.4, ' alam=', F9.5)") f, fold, alam

            IF ( alam < alamin ) THEN
                WRITE(*,*) alam, alamin
                CALL warn('Zero shift has been applied. Refinement has converged...', srname)
                x = xold
                CALL set_xold(xold, pdb_2)
                check = .TRUE.
                RETURN

            ELSE IF ( f < fold + ALF * alam * slope ) THEN

!                WRITE(*,"(' BOOK_LINE_SEARCH> ', 'f=', ES11.4, ' fold=', ES11.4, ' alam=', F9.5, ' ALF * alam * slope=', ES9.2 )") &
!                                                  f, fold, alam, ALF * alam * slope
!               Sufficient function decrease:
                CALL apply_shifts(p, pdb_2, alam, final=.TRUE.)
                CALL messag('Sufficient function decrease. Line search stopped.', srname)
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
    END SUBROUTINE book_line_search
        
END MODULE linear_search
