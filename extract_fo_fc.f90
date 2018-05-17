MODULE extract_fo_fc 
!USE class_vector
USE constants
USE fail
USE map_fft
USE map_operations
USE fft_util
USE util
USE vectors
IMPLICIT NONE
INTERFACE abs_fobs_in_true_space_group
    MODULE PROCEDURE abs_fobs_in_true_space_group_wp
END INTERFACE

INTERFACE fcalc_in_true_space_group
    MODULE PROCEDURE fcalc_in_true_space_group_wp
END INTERFACE

INTERFACE fcalc_in_P1_space_group
    MODULE PROCEDURE fcalc_in_P1_space_group_wp
END INTERFACE

INTERFACE OPERATOR ( .MFOBS. )
    MODULE PROCEDURE get_mfobs
END INTERFACE
CONTAINS
    SUBROUTINE abs_fobs_in_true_space_group_sp ( fobs, mtz_1 )
!
!       Purpose:
!       =======
!       Converts double precision mtz_1%fo array to single precision array fobs
!
        REAL(KIND=sp), DIMENSION(:), ALLOCATABLE            :: fobs
        TYPE(mtz),                               INTENT(IN) :: mtz_1

        CALL messag(' ', 'abs_fobs_in_true_space_group_sp')
        IF ( .NOT. ALLOCATED (fobs) ) THEN
            CALL messag('Programming error: fobs array has not been initialized.', 'abs_fobs_in_true_space_group_sp')
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fo ) ) THEN
            CALL messag('Programming error: mtz_1%fo array has not been initialized.', 'abs_fobs_in_true_space_group_sp')
        ENDIF

        CALL messag('Converting fobs to single precision.', 'abs_fobs_in_true_space_group_sp')

        fobs = mtz_1%fo(1:mtz_1%current_number_of_reflections)
        CALL messag( ' ', 'abs_fobs_in_true_space_group_sp')

    END SUBROUTINE abs_fobs_in_true_space_group_sp 

    SUBROUTINE abs_fobs_in_true_space_group_wp ( fobs, mtz_1 )
!
!       Purpose:
!       =======
!       Converts double precision mtz_1%fo into allocatable double precision array fobs
!
!
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: fobs
        TYPE(mtz),                   INTENT(IN)  :: mtz_1

!       Checkz:
        IF ( .NOT. ALLOCATED (fobs) ) THEN
            CALL die('Programming error: fobs array has not been initialized.', 'abs_fobs_in_true_space_group_wp')
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fo ) ) THEN
            CALL messag('Programming error: mtz_1%fo array has not been initialized.', 'abs_fobs_in_true_space_group_wp')
        ELSE IF ( SIZE(fobs) /= SIZE(mtz_1%fo) ) THEN
            CALL die('Programming error: mtz_1%fo array has not been initialized.', 'abs_fobs_in_true_space_group_wp')
        ENDIF

        fobs = mtz_1%fo(1:mtz_1%current_number_of_reflections)
    END SUBROUTINE abs_fobs_in_true_space_group_wp

    FUNCTION get_mfobs ( mode, mtz_1 ) RESULT ( fobs )
!
!       Purpose:
!       =======
!       Converts  mtz_1%fo, mtz_1%phio and mtz_1%fomo arrays into allocatable complex fobs array.
!
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE            :: fobs
        CHARACTER(LEN=*),                           INTENT(IN) :: mode                                       
        TYPE(mtz) ,                                 INTENT(IN) :: mtz_1
!       Local variables:
        CHARACTER(LEN=LEN(mode))                               :: string

        string = ADJUSTL ( mode )
        IF ( string(1:2) == 'P1' ) THEN
            IF ( ALLOCATED ( mtz_1%fobs) ) THEN
                CALL allocate_array ( fobs, SIZE ( mtz_1%fobs ) )
                fobs = mtz_1%fobs
            ELSE
                CALL die ( 'Programming error. MTZ_1%FOBS array has not been allocated...', 'get_mfobs' )
            ENDIF 
        ELSE
            CALL allocate_array ( fobs, SIZE ( mtz_1%fo ) )
            fobs = mtz_1%fomo * CMPLX ( mtz_1%fo * COS ( mtz_1%phio ), mtz_1%fo * SIN ( mtz_1%phio ), KIND=wp )
        ENDIF
   END FUNCTION get_mfobs

    SUBROUTINE fcalc_in_true_space_group_sp ( fcalc, mtz_1, FC_map, txyz )
!
!
!       Purpose:
!       =======
!       Calculates COMPLEX SINGLE PRECISION fcalc ARRAY in true space group.
!       Supply mtz_1, and FC_map containing structure factors in P1 space group
!
!       WARNING: txyz must be in fractions of u.c.
!
        COMPLEX(KIND=sp), DIMENSION(:), INTENT(INOUT) :: fcalc
        TYPE(mtz),        INTENT(IN)                  :: mtz_1
        TYPE(map),        INTENT(IN)                  :: FC_map
        TYPE(vector),     INTENT(IN), OPTIONAL        :: txyz
!       Local variables:
        TYPE(vector_int)                              :: hkl
        TYPE(vector_int)                              :: hCj
        REAL(KIND=wp)                                 :: hdj
        REAL(KIND=wp)                                 :: fftscale
        COMPLEX(KIND=wp)                              :: fcalCj
        COMPLEX(KIND=wp)                              :: fc
!       Counters:
        INTEGER                                       :: i
        INTEGER                                       :: j

!       Checks first:
        IF ( SIZE ( fcalc ) /= mtz_1%current_number_of_reflections ) THEN
            CALL warn('Too small size of fcalc array has been detected.', 'fcalc_in_true_space_group') 
        ENDIF

!       Conversion factor:        
        fftscale = mtz_1%volume/(FC_map%nu * FC_map%nv * FC_map%nw )

!       Loop over ASU in reciprocal space:
        DO i = 1, SIZE ( fcalc )
            hkl = mtz_1%hkl(i)
            fc  = (0.0_wp, 0.0_wp)
            DO j = 1, mtz_1%sp_group%number_of_symops
                CALL get_sym_refl ( hCj, hdj, fcalCj, hkl, mtz_1%sp_group, j, FC_map )
                IF ( PRESENT ( txyz ) ) THEN
                    fc = fc + fcalCj * phase ( twopi * ( hdj + (hCj .DOT. txyz) ) )
                ELSE
                    fc = fc + fcalCj * phase( twopi * hdj )
                ENDIF
            ENDDO

!           Apply scale/B-factor correction to fcalc:
            fcalc(i) = ( fftscale * EXP ( mtz_1%b_scale * mtz_1%s(i) /4.0_wp ) ) * fc
        ENDDO

    END SUBROUTINE fcalc_in_true_space_group_sp

    SUBROUTINE fcalc_in_true_space_group_wp ( fcalc, mtz_1, FC_map, txyz )
!
!       Purpose:
!       =======
!       Calculates COMPLEX DOUBLE PRECISION fcalc array in true space group.
!       Supply mtz_1 and FC_map containing structure factors in P1 space group
!       BUG CORRECTED mtz_1%curren_number_Of_reflections has been replaced with SIZE (fcalc).
!       Some checks added. APR 2010 BVS
!
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT)        :: fcalc
        TYPE(mtz),                      INTENT(IN)           :: mtz_1
        TYPE(map),                      INTENT(IN)           :: FC_map
        TYPE(vector),                   INTENT(IN), OPTIONAL :: txyz
!       Local variables:
        TYPE(vector_int)                                     :: hkl
        TYPE(vector_int)                                     :: hCj
        REAL(KIND=wp)                                        :: hdj
        REAL(KIND=wp)                                        :: fftscale
        COMPLEX(KIND=wp)                                     :: fcalCj
        COMPLEX(KIND=wp)                                     :: fc
!       Counters:
        INTEGER                                              :: i
        INTEGER                                              :: j
        CHARACTER(LEN=32),                              SAVE :: srname='fcalc_in_true_space_group'     
!       Checks first:
        IF ( SIZE ( fcalc ) /= mtz_1%current_number_of_reflections ) THEN
            WRITE(*,*) SIZE(fcalc), mtz_1%current_number_of_reflections
            CALL warn('Too small size of fcalc array has been detected.', srname)
        ENDIF

        IF ( SIZE ( fcalc ) /=  SIZE ( mtz_1%hkl ) ) THEN
            WRITE(*,*) SIZE(fcalc),SIZE ( mtz_1%hkl ) 
            CALL die('Inconsistent array sizes.', srname)
        ENDIF

!       Conversion factor:
        fftscale = mtz_1%volume/(FC_map%nu * FC_map%nv * FC_map%nw )

!       Loop over ASU in reciprocal space:
        DO i = 1, SIZE ( fcalc )

            hkl = mtz_1%hkl(i)
            fc  = (0.0_wp, 0.0_wp)

            DO j = 1, mtz_1%sp_group%number_of_symops

                CALL get_sym_refl ( hCj, hdj, fcalCj, hkl, mtz_1%sp_group, j, FC_map )

                IF ( PRESENT ( txyz ) ) THEN
                    fc = fc + fcalCj * phase ( twopi * ( hdj + (hCj .DOT. txyz) ) )
                ELSE    
                    fc = fc + fcalCj * phase( twopi * hdj )
                ENDIF

            ENDDO

!           Apply scale/B-factor correction to fcalc:
            fcalc(i) = ( fftscale * EXP ( mtz_1%b_scale * mtz_1%s(i) / 4.0_wp) ) * fc
        ENDDO

    END SUBROUTINE fcalc_in_true_space_group_wp

    SUBROUTINE fcalc_in_P1_space_group_wp ( fcalc, mtz_1, FC_map, txyz, xnsym )
!
!       Purpose:
!       =======
!       Calculates fcalc array in hemisphere but for true space group.
!       Supply mtz_1, and FC_map containing structure factors in P1 space group
!
        COMPLEX(KIND=wp),      DIMENSION(:)          :: fcalc
        TYPE(mtz),             INTENT(INOUT)         :: mtz_1      ! BUG CORRECTED APR 2005
        TYPE(map),             INTENT(IN)            :: FC_map
        TYPE(vector),          INTENT(IN), OPTIONAL  :: txyz
        INTEGER,               INTENT(IN), OPTIONAL  :: xnsym 
!       Local variables:
        TYPE(vector_int)                             :: hkl
        TYPE(vector_int)                             :: hCj
        REAL(KIND=wp)                                :: hdj
        REAL(KIND=wp)                                :: fftscale
        COMPLEX(KIND=wp)                             :: fcalCj
!        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE  :: fc
        INTEGER                                      :: nsym
!       Counters:
        INTEGER                                      :: i
        INTEGER                                      :: j

!       Checks first:
        IF ( SIZE ( fcalc ) /= SIZE ( mtz_1%hkl_in_P1 ) ) THEN
            WRITE(*,*) SIZE ( fcalc ), SIZE ( mtz_1%hkl_in_P1 )
            CALL die ( 'Too small size of FCALC array has been detected...',&
                       'fcalc_in_P1_sp_group_wp' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%s_in_P1 ) ) THEN
            IF ( SIZE ( mtz_1%s_in_P1 ) /= SIZE ( fcalc ) ) THEN
                WRITE(*,*) SIZE ( mtz_1%s_in_P1 ), SIZE ( fcalc )
                CALL die('Programing error. Sizes of S_SQUARED and FCALC do not match.',&
                         'fcalc_in_P1_sp_group_wp')
            ENDIF
        ENDIF

!       FFT conversion factor:
        fftscale = mtz_1%volume / (FC_map%nu * FC_map%nv * FC_map%nw)

!       Allocate and initialize fc array:
!!        CALL allocate_array ( fc, SIZE ( mtz_1%hkl_in_P1 ) )
        fcalc = ( 0.0_wp, 0.0_wp )

!       
        IF ( PRESENT ( xnsym ) ) THEN
            nsym = xnsym
        ELSE
            nsym =  mtz_1%sp_group%number_of_symops
        ENDIF

!       Loop over P1 ASU in reciprocal space:
        DO i = 1, SIZE ( mtz_1%hkl_in_P1 )

            hkl = mtz_1%hkl_in_P1(i)

            DO j = 1, nsym

                CALL get_sym_refl ( hCj, hdj, fcalCj, hkl, mtz_1%sp_group, j, FC_map )

                IF ( PRESENT ( txyz ) ) THEN
                    fcalc(i) = fcalc(i) + fcalCj * phase ( twopi * ( hdj + (hCj .DOT. txyz) ) )
                ELSE
                    fcalc(i) = fcalc(i) + fcalCj * phase( twopi * hdj )
                ENDIF

            ENDDO

        ENDDO

!       Apply scale/B-factor correction to fcalc:
        IF ( ALLOCATED ( mtz_1%s_in_P1 ) ) THEN
            fcalc = ( fftscale * EXP ( mtz_1%b_scale * mtz_1%s_in_P1 / 4.0_wp ) ) * fcalc
        ELSE
            fcalc = fftscale * fcalc
        ENDIF

!       Free fc array:
!        CALL deallocate_array ( fc )

    END SUBROUTINE fcalc_in_P1_space_group_wp

END MODULE extract_fo_fc
