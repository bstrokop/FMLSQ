MODULE map_operations
USE map_fft
IMPLICIT NONE
INTERFACE MAXVAL
    MODULE PROCEDURE map_maxval
END INTERFACE
INTERFACE MINVAL
    MODULE PROCEDURE map_minval
END INTERFACE
CONTAINS
!   Special map operations:
    SUBROUTINE binary_map_multiplication ( map_result, map_1, map_2 )
!
!       Purpose:
!       =======
!       Calculate binary product of 2 maps.
!       resulting map = T3 map .x. packing map. 
!       If grid point in both maps has positive value
!       then value from T3 map is used to form a resulting map 
!     
!       *** Important: T3 map MUST correspond to map_1
!
!       Date:            Programmer:          Description of changes:
!       ====             ==========           ======================
!       Nov 2003         Strokopytov B.       Original code
!       Nov 2005             -"-              Cosmetic changes, comments
!
        TYPE(map), INTENT(INOUT) :: map_result
        TYPE(map), INTENT(IN)    :: map_1
        TYPE(map), INTENT(INOUT) :: map_2
!       Test:
        LOGICAL                  :: non_zero
!       Counters:
        INTEGER                  :: i
        INTEGER                  :: j
        INTEGER                  :: k

        non_zero = .FALSE.
        DO k = 0, map_1%nw - 1
            DO j = 0, map_1%nv - 1
                DO i = 0, map_1%nu - 1
                    IF ( map_1%array(i,j,k) > 0.0_wp .AND. map_2%array(i,j,k) > 0.0_wp ) THEN
                        map_result%array(i,j,k) = map_1%array(i,j,k)
                        non_zero = .TRUE.
                    ELSE
                        map_result%array(i,j,k) = 0.0_wp
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        IF (.NOT. non_zero ) CALL die('Not a single grid point in the map is positive.', 'multiply_maps')
    END SUBROUTINE binary_map_multiplication

    FUNCTION map_sigma ( map_1 )
!
!       Purpose:
!       =======
!       Calculates map sigma. Suitable for maps calculated in P1 only.
!
!       Date:        Programmer:       Description of changes:
!       ====         ==========        ======================
!       Nov 2005     Strokopytov B.    Original code
!
!
        REAL(KIND=wp)         :: map_sigma
        TYPE(map), INTENT(IN) :: map_1
!       Local variables:
        REAL(KIND=wp)         :: map_scale
        REAL(KIND=wp)         :: sumyc
!       Sums:
        REAL(KIND=wp)         :: sumc2

        sumc2 = SUM ( map_1%array(0:map_1%nu-1, 0:map_1%nv-1, 0:map_1%nw-1)**2 )
        sumyc = SUM ( map_1%array(0:map_1%nu-1, 0:map_1%nv-1, 0:map_1%nw-1) )

        map_scale = REAL(map_1%nu, KIND=wp) * REAL(map_1%nv, KIND=wp) * REAL(map_1%nw, KIND=wp)
        map_sigma = SQRT( ( sumc2 - sumyc ** 2 / map_scale ) / map_scale )
    END FUNCTION map_sigma

    FUNCTION map_positivity ( map_1, cutoff )
!
!       Purpose:
!       =======
!       Calculates number of grid points above cutoff value.
!       Suitable for maps calculated in P1 only
!
!       Date:        Programmer:       Description of changes:
!       ====         ==========        ======================
!       Nov 2005     Strokopytov B.    Original code
!
        INTEGER                   :: map_positivity
        TYPE(map),     INTENT(IN) :: map_1
        REAL(KIND=wp), INTENT(IN) :: cutoff

        map_positivity = COUNT ( map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) > cutoff ) 
    END FUNCTION map_positivity

    FUNCTION map_maxval ( map_1 )
!
!       Purpose:
!       =======
!       Calculates number of grid points above cutoff value.
!       Suitable for maps calculated in P1 only
!
!       Date:        Programmer:       Description of changes:
!       ====         ==========        ======================
!       Aug 2007     Strokopytov B.    Original code
!
        REAL(KIND=wp)             :: map_maxval
        TYPE(map),     INTENT(IN) :: map_1

        map_maxval =  MAXVAL ( map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) )
    END FUNCTION map_maxval

    FUNCTION map_minval ( map_1 )
!
!       Purpose:
!       =======
!       Calculates number of grid points above cutoff value.
!       Suitable for maps calculated in P1 only
!
!       Date:        Programmer:       Description of changes:
!       ====         ==========        ======================
!       Nov 2007     Strokopytov B.    Original code
!
        REAL(KIND=wp)             :: map_minval
        TYPE(map),     INTENT(IN) :: map_1

        map_minval =  MINVAL ( map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) )
    END FUNCTION map_minval

    SUBROUTINE negate_map ( map_1, c )
        TYPE(map),     INTENT(INOUT) :: map_1
        REAL(KIND=wp), INTENT(IN)    :: c

        map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) = &
        map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) - ABS(c)

    END SUBROUTINE negate_map

    FUNCTION count_map ( map_1 )
        INTEGER(KIND=eb)          :: count_map
        TYPE(map),     INTENT(IN) :: map_1

!        count_map = COUNT ( map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) == c, KIND=eb ) 
        count_map = COUNT ( map_1%array(0:map_1%nu-1,0:map_1%nv-1,0:map_1%nw-1) > 0.0_wp, KIND=eb ) 

    END FUNCTION count_map

    FUNCTION how_much_solvent ( map_1 )
!
!       Purpose:
!       =======
!       This allows to estimate volume occupied by the model.
!       Some overestimation is possible. Better routine is needed.
        REAL(KIND=wp)             :: how_much_solvent
        TYPE(map),     INTENT(IN) :: map_1

        how_much_solvent = 1.0_wp - map_1%sp_group%number_of_symops *  count_map(map_1) &
                         / REAL ( map_1%nu * map_1%nv * map_1%nw )

    END FUNCTION how_much_solvent

    FUNCTION convert_symop_to_ort_system_map ( sym_1, map_2 )
        TYPE(symop)             :: convert_symop_to_ort_system_map
        TYPE(symop), INTENT(IN) :: sym_1
        TYPE(map),   INTENT(IN) :: map_2

        convert_symop_to_ort_system_map = ( map_2%ORT * sym_1 ) * map_2%DEORT
    END FUNCTION convert_symop_to_ort_system_map

    SUBROUTINE get_sym_refl ( hCj, hdj, fCj, hkl, sp_group, j, map_2 )
        TYPE(vector_int),  INTENT(OUT) :: hCj
        REAL(KIND=wp),     INTENT(OUT) :: hdj
        COMPLEX(KIND=wp),  INTENT(OUT) :: fCj
        TYPE(vector_int),  INTENT(IN)  :: hkl
        TYPE(space_group), INTENT(IN)  :: sp_group
        INTEGER,           INTENT(IN)  :: j
        TYPE(map),         INTENT(IN)  :: map_2
!       Local variables:
        INTEGER, DIMENSION(3)          :: uvw

!       Some speed up:
        IF ( j == 1 ) THEN
            hCj = hkl
            hdj = 0.0_wp
        ELSE
            hCj = sp_group%SYM_HKL(j) * hkl
            hdj = hkl .DOT. .SYMV. sp_group%SYM(j)
        ENDIF

!       Note that FFTW produces complex conjugate:
        IF ( hemi_h ( hCj ) ) THEN
            uvw  = hCj .HTUMOD. map_2%map_size
            fCj  = CMPLX ( map_2%array ( uvw(1)    , uvw(2), uvw(3) ), &
                          -map_2%array ( uvw(1) + 1, uvw(2), uvw(3) ), KIND=wp )
        ELSE
            uvw = -hCj .HTUMOD. map_2%map_size

!           This is complex conjugate -- no plus sign below for imaginary part:
            fCj =  CMPLX ( map_2%array ( uvw(1)    , uvw(2), uvw(3) ), &
                           map_2%array ( uvw(1) + 1, uvw(2), uvw(3) ), KIND=wp )
        ENDIF

    END SUBROUTINE get_sym_refl

END MODULE map_operations
