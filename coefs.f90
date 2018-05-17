MODULE coefs
USE aniso_symmetry_manip
USE basic_symmetry_operations
!USE bulk_solvent_manip
USE constants
USE fail
USE map_fft
USE map_operations
USE fft_util
USE symmetry_manip
USE util
USE vectors
IMPLICIT NONE
CONTAINS
    SUBROUTINE expand_fobs_to_hemisphere ( mtz_1 )
!
!       Purpose:
!       =======
!       Expands Fobs to P1 hemisphere
!       and multiplies Fobs with figure of merit.
!
        TYPE(mtz) ,INTENT(INOUT) :: mtz_1
!       Local variables:
        TYPE (vector_int)        :: hklm
        REAL(KIND=wp)            :: dphim
        COMPLEX(KIND=wp)         :: fobs
        COMPLEX(KIND=wp)         :: fobsm
        INTEGER, DIMENSION(3)    :: h
        INTEGER, DIMENSION(3)    :: mtzh
!       Counters:
        INTEGER                  :: i
        INTEGER                  :: m
        INTEGER                  :: total_number_in_P1

        total_number_in_P1 = 0
        DO i = 1, mtz_1%current_number_of_reflections
            fobs  = CMPLX ( mtz_1%fo(i) * COS(mtz_1%phio(i)), &
                            mtz_1%fo(i) * SIN(mtz_1%phio(i)), KIND=wp)

            DO m = 1, mtz_1%sp_group%number_of_primitive_symops
                total_number_in_P1 = total_number_in_P1 + 1
                hklm               = mtz_1%sp_group%SYM_HKL(m) * mtz_1%hkl(i)
                IF ( m > 1 .AND. hklm == mtz_1%hkl(i) ) THEN
                    mtzh = mtz_1%hkl(i)
                    h = hklm
                    WRITE(*,'('' EXPAND_FOBS_TO_HEMISPHERE: Duplicate detected'')')
                    WRITE(*,'('' EXPAND_FOBS_TO_HEMISPHERE: '','' hkl = '', 3I4, '' hklm = '',3I4)')&
                    mtzh, h
                ENDIF

!               This is correct (at least for transposed symmetry matrices): 
                dphim = twopi * ( hklm .DOT. .SYMV. mtz_1%sp_group%SYM(m) )
                fobsm = mtz_1%fomo(i) * fobs * phase ( -dphim )

                IF ( hemi_h ( hklm ) ) THEN
                    mtz_1%fobs(total_number_in_P1) = fobsm
                ELSE
                    hklm = -hklm
                    mtz_1%fobs(total_number_in_P1) = CONJG( fobsm )
                ENDIF

                IF ( MOD(i, 1000) == 1 ) THEN
                    h = hklm
                    WRITE(*, '('' EXPAND_FOBS_TO_HEMISPHERE> '',''hkl index'', 3I5, &
                              &'' Fom*Fobs = '', 2F10.2, &
                              &'' Fom = '',F10.6)') &
                                  h, mtz_1%fobs(total_number_in_P1), mtz_1%fomo(i)
                ENDIF

                IF ( m > 1 ) THEN
                    IF ( hklm == mtz_1%hkl(i) ) CYCLE   ! This duplicate
                    IF (-hklm == mtz_1%hkl(i) ) CYCLE   ! This is centric refl. producing duplicate
                ENDIF

            ENDDO
        ENDDO

        CALL messag('Expansion finished.','expand_fobs_to_hemisphere')
        WRITE(*,'('' EXPAND_FOBS_TO_HEMISPHERE> '',''Total number of F(obs) in P1 hemisphere ='',I8)')&
        total_number_in_P1
        CALL messag( ' ', 'expand_fobs_to_hemisphere')
    END SUBROUTINE expand_fobs_to_hemisphere

    SUBROUTINE phased_translation_fun_coefs ( mtz_1, map_1, phase_sign )
!
!       Purpose:
!       =======
!       Calculates phased translation function coefficients in P1 space group
!       and assign them to appopriate map_1 indices
!
!       Date:            Programmer:                 Description of changes:
!       ====             ==========                  ======================
!       Nov 2003         B.Strokopytov               Original code
!       Nov 2005         B.Strokopytov               2 long-stading bugs corrected:
!                                                    a) friedel logical variable introduced
!                                                    b) CONJG ( phased_translation_function ) if not
!                                                       in hemisphere
!
!       Nov 2005         B.Strokopytov               Almost complete rewrite. Using GET_SYM_REFL now
!                                                    Last part of the routine removed. Have to use
!                                                    simple_map_coefs_in_P1 from now on.
!       Note:
!       ====
!       The above bugs have not affected this routine since list of hkl was already in hemisphere
!       However, for general hkl list this would work incorrectly...
!
        TYPE (mtz),        INTENT(INOUT) :: mtz_1
        TYPE (map),        INTENT(INOUT) :: map_1
        INTEGER,           INTENT(IN)    :: phase_sign
!       Local variables: 
        TYPE (vector_int)                :: hklm
        REAL(KIND=wp)                    :: hdj
        REAL(KIND=wp)                    :: fftscale
        REAL(KIND=wp)                    :: cc_scale_factor
!       Counters:
        INTEGER                          :: i

        fftscale = mtz_1%volume / (map_1%nu * map_1%nv * map_1%nw )

        IF ( .NOT. ALLOCATED ( mtz_1%fc ) ) THEN
            CALL allocate_array ( mtz_1%fc, SIZE ( mtz_1%hkl_in_P1 ) )
        ENDIF

!       Get structure factor in P1 space group from a raw map:
        DO i = 1, SIZE ( mtz_1%hkl_in_P1 )
            CALL get_sym_refl ( hklm, hdj, mtz_1%fc(i), mtz_1%hkl_in_P1(i), mtz_1%sp_group, 1, map_1 )
        ENDDO

!       apply B-factor correction to fc (in P1) since we are getting FCs from raw map:
        mtz_1%fc = fftscale *  mtz_1%fc  *  EXP ( mtz_1%b_scale * 0.25_wp * mtz_1%s_in_P1 )        

!       factor of 0.5 ( can divide by 2 also ) to account for the whole sphere:
        cc_scale_factor = ( 0.5_wp * mtz_1%volume )  /  SQRT ( SUM ( REAL ( mtz_1%fobs * CONJG ( mtz_1%fobs ), KIND=wp ) ) &
                                                     * SUM ( REAL ( mtz_1%fc * CONJG ( mtz_1%fc ), KIND=wp ) ) )

!       Final scaling:
        IF ( phase_sign > 0 ) THEN
            mtz_1%fc_in_P1 = cc_scale_factor * mtz_1%fobs * CONJG ( mtz_1%fc )
        ELSE
            mtz_1%fc_in_P1 = cc_scale_factor * CONJG ( mtz_1%fobs ) * CONJG ( mtz_1%fc )
        ENDIF

!       Done. Ready to use SIMPLE_MAP_COEFS_IN_P1

    END SUBROUTINE phased_translation_fun_coefs

    FUNCTION simple_ptf_coefs ( mtz_1, ips ) RESULT ( ptf_coefs )
!
!       Purpose:
!       =======
!       Calculates phased translation function coefficients in P1 in most simple form:
!
!
!       Date:            Programmer:                  Description_of_changes:
!       ====             ==========                   ======================
!       Nov 2005         B.Strokopytov                Original code
!
!         
!       Note:
!       ====
!       This routine implicitly assumes that both fobs and fcalc are in the same hemisphere
!       i.e. MTZ_1%FOBS and MTZ_1%FC was calculated for the same list of reflections which
!       is normally the case. Should be used in conjuction with subroutine 
!       REDUCE_FOBS_TO_HEMISPHERE which moves MTZ_1%HKL_IN_P1 and MTZ_1%FOBS.
!
!       Currenly this simple approach is NOT suitable for calculation in 'true' space group.
!       However, if initial list of hkls is in the unique hemisphere this routine may and
!       shoud be used.
!
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE            :: ptf_coefs
        TYPE(mtz),                                  INTENT(IN) :: mtz_1
        INTEGER,                                    INTENT(IN) :: ips
!       Local variables:
        REAL(KIND=wp)                                          :: cc_scale

        CALL allocate_array ( ptf_coefs, SIZE ( mtz_1%hkl_in_P1 ) )

!       Form phased translation function coefs:
        IF ( ips > 0 ) THEN
            ptf_coefs = mtz_1%fobs * CONJG ( mtz_1%fc )
        ELSE
            ptf_coefs = CONJG ( mtz_1%fobs ) * CONJG ( mtz_1%fc )
        ENDIF

!       Scale the result:
        cc_scale =  0.5_wp * mtz_1%volume &
                 /  SQRT ( SUM ( REAL ( mtz_1%fobs * CONJG ( mtz_1%fobs ), KIND=wp ) ) &
                         * SUM ( REAL ( mtz_1%fc * CONJG ( mtz_1%fc ), KIND=wp ) ) )

        ptf_coefs = cc_scale * ptf_coefs
        
    END FUNCTION simple_ptf_coefs

    SUBROUTINE simple_map_coefs_in_P1 ( mtz_1, map_1, iprint )
!
!       Purpose:
!       =======
!
!       Fill map_1 with coefficients from mtz_1%fcalc in P1 space group:
!
!       Date:             Programmer:            History of changes:
!       ====              ==========             ==================
!       Oct 2004          B.Strokopytov          Original code
!       Oct 2007          B.Strokopytov          Cosmetic changes
!
        TYPE(mtz),             INTENT(IN)           :: mtz_1
        TYPE(map),             INTENT(INOUT)        :: map_1
        INTEGER,               INTENT(IN), OPTIONAL :: iprint
!
!       Local variables:
        TYPE(vector_int)                            :: hklm
        REAL(KIND=wp)                               :: fftscale
        INTEGER, DIMENSION(3)                       :: uvw
!
!       Counters:
        INTEGER                                     :: i
        CHARACTER(LEN=32),                     SAVE :: srname = 'simple_map_coefs_in_P1'

!       Checkz:
        IF (.NOT. ALLOCATED ( mtz_1%fc_in_P1 ) ) THEN
 
           CALL die('Programming error. mtz_1%fcalc array has not been allocated properly.',&
                    srname)

        ELSE IF (.NOT. ALLOCATED ( mtz_1%hkl_in_P1 ) ) THEN

            CALL die('Programming error. mtz_1%hkl_in_P1 array has not been allocated properly.',&
                     srname)

        ENDIF

        fftscale = 1.0_wp / map_1%volume

!       Initialize map:
        map_1%array = 0.0_wp

        DO i = 1, SIZE ( mtz_1%hkl_in_P1 )

            hklm = mtz_1%hkl_in_P1(i)
            uvw  = hklm

            IF ( PRESENT ( iprint ) ) THEN
                WRITE(*,"(' SIMPLE_MAP_COEFS_IN_P1> ', 'hkl ', 3I4, 2F10.2, I12)")&
                uvw, mtz_1%fc_in_P1(i), i
            ENDIF

            IF ( hemi_h ( hklm ) ) THEN

                uvw =  hklm .HTUMOD. map_1%map_size
                map_1%array(uvw(1),   uvw(2), uvw(3)) =  fftscale * REAL  ( mtz_1%fc_in_P1(i) )
                map_1%array(uvw(1)+1, uvw(2), uvw(3)) = -fftscale * AIMAG ( mtz_1%fc_in_P1(i) ) 

            ELSE

                uvw = -hklm .HTUMOD. map_1%map_size
                map_1%array(uvw(1),   uvw(2), uvw(3)) = fftscale * REAL  ( mtz_1%fc_in_P1(i) )
                map_1%array(uvw(1)+1, uvw(2), uvw(3)) = fftscale * AIMAG ( mtz_1%fc_in_P1(i) )

            ENDIF

        ENDDO

!       Ready to perform hermitian-to-real FFT...

    END SUBROUTINE simple_map_coefs_in_P1

    SUBROUTINE matrix_magic_coefs(mtz_1, map_1, mode)
!
!       Purpose:
!       =======
!       Calculates coefficients for MAX_BLOCK_SIZE maps to be used during convolution stage.
!
!       Date:             Programmer:                 History of changes:
!       ====              ==========                  ==================
!       Nov 2009          Boris Strokopytov           Original code
!       Mar 2009          Boris Strokopytov           ADPs added
!
        TYPE(mtz),                                    INTENT(INOUT) :: mtz_1
        TYPE(map),         DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: map_1
        CHARACTER(LEN=*),                             INTENT(IN)    :: mode
!       Local variables:
        INTEGER,           DIMENSION(3)                             :: uvw
        TYPE(vector_int)                                            :: hklm
        COMPLEX(KIND=wp),  DIMENSION(MAX_BLOCK_SIZE)                :: coefs
        REAL(KIND=wp),     DIMENSION(6)                             :: a
        REAL(KIND=wp),     DIMENSION(6)                             :: b
        REAL(KIND=wp)                                               :: hdm
        INTEGER                                                     :: nr
        INTEGER                                                     :: nsym
        REAL(KIND=wp)                                               :: fftscale
!       Refinement:
        LOGICAL                                                     :: refine_xyz
        LOGICAL                                                     :: refine_biso
        LOGICAL                                                     :: refine_occ
        LOGICAL                                                     :: refine_Us
        INTEGER                                                     :: my_first_allocated_map
!       Special matrices:
        TYPE(matrix),      DIMENSION(:), ALLOCATABLE                :: ORT_SYM_DEORT
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE                :: RU
!       Counters:
        INTEGER                                                     :: i
        INTEGER                                                     :: m  
        INTEGER                                                     :: imap      

!       Figure out what to refine:
        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

!       First allocated map:
        my_first_allocated_map = first_allocated_map(map_1)

        fftscale = 1.0_wp / map_1(my_first_allocated_map)%volume

        nsym = map_1(my_first_allocated_map)%sp_group%number_of_symops

!       Allocate special matrices:
        CALL allocate_array ( ORT_SYM_DEORT, nsym )
        CALL allocate_array ( RU, nsym )
!                                                                  T
!       Have to use TRANSPOSE because of identity: d[OAD]x/dx=[OAD] :
        DO m = 1, nsym
            ORT_SYM_DEORT(m) = map_1(my_first_allocated_map)%ORT                    &
                             * .SYMA. map_1(my_first_allocated_map)%sp_group%SYM(m) &
                             * map_1(my_first_allocated_map)%DEORT
            RU(m) = TRANSPOSE ( .SIXBYSIX. ORT_SYM_DEORT(m) )
!            WRITE(*,*) ' Number of zeroes in RU matrix:', COUNT(RU(m)%a == 0.0_wp)
            ORT_SYM_DEORT(m) = TRANSPOSE ( ORT_SYM_DEORT(m) )
        ENDDO

        IF ( refine_xyz ) THEN

!THis can be handled by the program automatically:
!$OMP PARALLEL DO !                    NUM_THREADS(3)
            DO imap = 1, 3
                map_1(imap)%array = 0.0_wp
            ENDDO
        ENDIF

        IF ( refine_biso ) THEN
            map_1(4)%array = 0.0_wp
        ENDIF

        IF ( refine_occ ) THEN
            map_1(5)%array = 0.0_wp
        ENDIF

        IF ( refine_Us ) THEN
!$OMP PARALLEL DO                      !NUM_THREADS(6)
            DO imap = 6, MAX_BLOCK_SIZE
                map_1(imap)%array = 0.0_wp
            ENDDO
        ENDIF

        nr = SIZE ( mtz_1%hkl_in_P1 )

!       These loops have been interchanged for parallelization (to avoid 'OMP CRITICAL' instructions):
        symmetry:DO m = 1, nsym

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (M, MTZ_1, MAP_1, MY_FIRST_ALLOCATED_MAP, REFINE_XYZ, REFINE_BISO, REFINE_OCC, &
!$OMP                                      REFINE_US, NR, ORT_SYM_DEORT, FFTSCALE, RU )

            DO i = 1, nr

!OMP CRITICAL
                hklm = mtz_1%sp_group%SYM_HKL(m) * mtz_1%hkl_in_P1(i)
!OMP END CRITICAL
                uvw  = hklm

!               Phase shift:
!OMP CRITICAL
                hdm  = -mtz_1%hkl_in_P1(i) .DOT. .SYMV. mtz_1%sp_group%SYM(m)
!OMP END CRITICAL

                IF ( refine_xyz ) THEN

                    coefs(1:3) = mtz_1%fcq_sym(i,1:3) * phase ( twopi * hdm )

!                   Super strange operation, actually summing up maps in a very special way:
                    a(1:3) = ORT_SYM_DEORT(m) * REAL  ( coefs(1:3), KIND=wp)
                    b(1:3) = ORT_SYM_DEORT(m) * AIMAG ( coefs(1:3) )

                    coefs(1:3) = CMPLX ( a(1:3), b(1:3), KIND=wp )
                ENDIF

                IF ( refine_biso ) THEN

                    coefs(4) = mtz_1%fcq_sym(i,4) * phase ( twopi * hdm )

                ENDIF


                IF ( refine_occ ) THEN

                    coefs(5) = mtz_1%fcq_sym(i,5) * phase ( twopi * hdm )                    

                ENDIF

                IF ( refine_Us ) THEN
                    coefs(6:MAX_BLOCK_SIZE) = mtz_1%fcq_sym(i,6:MAX_BLOCK_SIZE) * phase ( twopi * hdm )

!                    This is dangerous for parallel computing (need other function for overloading '*':
                    a(1:6) = RU(m) * REAL  ( coefs(6:MAX_BLOCK_SIZE), KIND=wp)
                    b(1:6) = RU(m) * AIMAG ( coefs(6:MAX_BLOCK_SIZE) ) 
!                   Multiply by TRANSPOSE according to Garib (see TRANSPOSE(RU) above):
!                    a(1:6) = MATMUL ( RU(m)%a, REAL ( coefs(6:MAX_BLOCK_SIZE), KIND=wp ) )
!                    b(1:6) = MATMUL ( RU(m)%a, AIMAG( coefS(6:MAX_BLOCK_SIZE) ) )
                    coefs(6:MAX_BLOCK_SIZE) = CMPLX ( a(1:6), b(1:6), KIND=wp )

                ENDIF
   
                hemisphere:IF ( hemi_h ( hklm ) ) THEN

                    uvw =  hklm .HTUMOD. map_1(my_first_allocated_map)%map_size

                    IF ( refine_xyz ) THEN
                        DO imap = 1,3
                            map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) =  map_1(imap)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                        + fftscale * REAL  ( coefs(imap) )
                            map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                        - fftscale * AIMAG ( coefs(imap) )
                        ENDDO
                    ENDIF

                    IF ( refine_biso ) THEN

                        map_1(4)%array(uvw(1),   uvw(2), uvw(3)) =  map_1(4)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                 +  fftscale * REAL  ( coefs(4) )
                        map_1(4)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(4)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                 - fftscale * AIMAG ( coefs(4) )

                    ENDIF

                    IF ( refine_occ ) THEN
                        map_1(5)%array(uvw(1),   uvw(2), uvw(3)) =  map_1(5)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                 +  fftscale * REAL  ( coefs(5) )
                        map_1(5)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(5)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                 - fftscale * AIMAG ( coefs(5) )
                    ENDIF

                    IF ( refine_Us ) THEN
                        DO imap = 6,MAX_BLOCK_SIZE
                            map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) =  map_1(imap)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                        + fftscale * REAL  ( coefs(imap) )
                            map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                        - fftscale * AIMAG ( coefs(imap) )
                        ENDDO
                    ENDIF

                ELSE

                    uvw = -hklm .HTUMOD. map_1(my_first_allocated_map)%map_size

                    IF ( refine_xyz ) THEN
                        DO imap = 1, 3
                            map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) = map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) &
                                                                        + fftscale * REAL  ( coefs(imap) )
                            map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                        + fftscale * AIMAG ( coefs(imap) )
                        ENDDO
                    ENDIF

                    IF ( refine_biso ) THEN
                        map_1(4)%array(uvw(1),   uvw(2), uvw(3)) = map_1(4)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                 + fftscale * REAL  ( coefs(4) )
                        map_1(4)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(4)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                 + fftscale * AIMAG ( coefs(4) )
                    ENDIF
                
                    IF ( refine_occ ) THEN
                        map_1(5)%array(uvw(1),   uvw(2), uvw(3)) =  map_1(5)%array(uvw(1),  uvw(2), uvw(3)) &
                                                                 +  fftscale * REAL  ( coefs(5) )
                        map_1(5)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(5)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                 + fftscale * AIMAG ( coefs(5) )
                    ENDIF

                    IF ( refine_Us ) THEN

                        DO imap = 6, MAX_BLOCK_SIZE
                            map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) = map_1(imap)%array(uvw(1),   uvw(2), uvw(3)) &
                                                                        + fftscale * REAL  ( coefs(imap) )
                            map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) = map_1(imap)%array(uvw(1)+1, uvw(2), uvw(3)) &
                                                                        + fftscale * AIMAG ( coefs(imap) )
                        ENDDO

                    ENDIF

                ENDIF hemisphere

            ENDDO
        ENDDO symmetry

!       Free memory:
        CALL deallocate_array(ORT_SYM_DEORT)
        CALL deallocate_array(RU)

    END SUBROUTINE matrix_magic_coefs

    SUBROUTINE Ax_matrix_vector_coefs(fcq, mtz_1, map_1, temp_complex, imap, mode, H2)
! 
!
!       Purpose:
!       =======
!       Calculates reciprocal space coefficients for general matrix-vector product in P1 space group.
!       Transformation to TRUE space group is done via MATRIX_MAGIC_COEFS.
!       Can be used for mixed second derivatives of X, Y, Z, B and Q
!
!       Date:              Programmer:           History of changes:
!       ====               ==========            ==================
!       Oct 2007           B.Strokopytov         Original code.
!       Nov 2007           B.Strokopytov         Added cross terms for XYZ contrary to Agarwal.
!                                                They maybe ~10% of diagonal terms.
!       Dec 2007           B.Strokopytov         Some loss of precision detected.
!                                                Added KIND=wp in CMPLX function for OCC refinement,
!                                                bug corrected.
!
!       
!       Note:
!       ====
!       a)Anisotropy code is on its way.
!
!       Theoretical note:
!       ================
!       Symmetries of these maps is still unknown (but not lower than P1).
!       Need to use P1 map alongside with this subroutine.
!    
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT)         :: fcq 
        TYPE(mtz),                      INTENT(INOUT)         :: mtz_1
!       Just one map is sufficient:
        TYPE(map),                      INTENT(IN)            :: map_1
!       Temporary array:
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT)         :: temp_complex
        INTEGER,                        INTENT(IN)            :: imap
        CHARACTER(LEN=*),               INTENT(IN)            :: mode
        LOGICAL,                        INTENT(IN)            :: H2
!       Local vars:
        REAL(KIND=wp)                                         :: fft_scale
        INTEGER                                               :: nrefl
!       Refinement:
        LOGICAL                                               :: refine_xyz
        LOGICAL                                               :: refine_biso
        LOGICAL                                               :: refine_occ
        LOGICAL                                               :: refine_Us
        CHARACTER(LEN=32),                               SAVE :: srname = 'matrix_vector_coefs_in_P1'


        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

        nrefl = SIZE ( mtz_1%fc_in_P1 )

!       Checkz:
        IF ( SIZE ( mtz_1%fo_in_P1 ) /= nrefl ) THEN

            WRITE(*,*) nrefl, SIZE ( mtz_1%fo_in_P1 )
            CALL messag('Programming error. Inconsistent sizes for fcalc and fo_in_P1',&
                        srname)
        ELSE IF ( SIZE ( fcq ) /= nrefl ) THEN

            WRITE(*,*) nrefl, SIZE ( fcq )
            CALL messag('Programming error. Inconsistent sizes for FC_IN_P1 and FCQ',&
                        srname)
        ELSE IF ( SIZE ( mtz_1%wphases ) /= nrefl ) THEN

            WRITE(*,*) nrefl, SIZE ( mtz_1%wphases )
            CALL messag('Programming error. Inconsistent sizes for FC_IN_P1 and PHASES',&
                        srname)

        ENDIF

!       Need this to scale fcalc properly otherwise the result will depend on
!       map grid:

        fft_scale = map_1%volume / ( map_1%nu * map_1%nv * map_1%nw )

        IF ( fft_scale <= 0.0_wp ) THEN
            WRITE(*,*) ' imap=', imap, ' fft_scale=', fft_scale, ' map_1%volume=', map_1%volume
            CALL die('Programming error. FFT scale factor is non-positive.', &
                     srname)
        ENDIF

!       Initialise FCQ:
        fcq = CMPLX ( 0.0_wp, 0.0_wp, KIND=wp )

!       Famous Agarwal(1978) - Tronrud(1999) - Strokopytov(2008) terms:
!       ---------------------------------------------------------------
        SELECT CASE ( imap )
       
        CASE (1)
        
!       4 pi ** 2 h ** 2 term:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ho ** 2, 0.0_wp, KIND=wp )

!       Multiply by fcq_in_P1 #1:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)                   ! H1 term
        IF ( H2 ) THEN 
            fcq = fcq  - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &  ! H2 term
                       * mtz_1%wphases                                                ! H2 term
        ENDIF

!       4 pi ** 2 h * k term:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ho * mtz_1%Ko, 0.0_wp, KIND=wp )

!       Multiply by fcq_in_P1 #2:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq  - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &  ! H2 term
                       * mtz_1%wphases                                                ! H2 term
        ENDIF


!       4 pi ** 2 h * l term:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ho * mtz_1%Lo, 0.0_wp, KIND=wp )

!       Multiply by fcq_in_P1 #3:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)                   ! H1 term
        IF ( H2 ) THEN
            fcq = fcq  - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &  ! H2 term
                       * mtz_1%wphases                                                ! H2 term
        ENDIF



        IF ( refine_biso ) THEN

!           pi * i * h * s ** 2 / 2 term:
            temp_complex = CMPLX ( 0.0_wp, pi * mtz_1%Ho * mtz_1%s_in_P1 / 2.0_wp, KIND=wp )

!           fcq_in_P1 #4
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)                  ! H1 term

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,4) )  & ! H2 term
                          * mtz_1%wphases                                                ! H2 term
            ENDIF

        ENDIF

        IF ( refine_occ ) THEN

!           2 * pi * i * h / O(j):
            temp_complex = CMPLX ( 0.0_wp, -twopi * mtz_1%Ho, KIND=wp ) 
        
!           fcq_in_P1 #5:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)
            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,5) ) &
                          * mtz_1%wphases
            ENDIF
        ENDIF

        aniso_row1:IF ( refine_Us ) THEN
!           4 pi ** 3 * i * h ** 3 (6):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ho ** 3, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF

!           4 pi ** 3 * i * h * k ** 2  (7):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko ** 2, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           4 pi ** 3 * i * h * l ** 2 (8):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Lo ** 2, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl, 8)
           
            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl, 8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 pi ** 3 * i * h ** 2 * k (9):
            temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Ko, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl, 9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl, 9) ) &
                          * mtz_1%wphases
            ENDIF

!           8 pi ** 3 * i * h ** 2 * l (10):
            temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Lo, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl, 10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl, 10) ) &
                          * mtz_1%wphases
            ENDIF

!           8 pi ** 3 * i * h * k * l (MAX_BLOCK_SIZE):
            temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl, MAX_BLOCK_SIZE)
                     

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl, MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF
 
        ENDIF aniso_row1
!       --------------------------------------------------------------------
        CASE (2)

!       Can't get here unless we refine Y coordinates:

!       4 pi ** 2 * H * K :
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ho * mtz_1%Ko, 0.0_wp, KIND=wp )

!       fcq_in_P1 #1:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF


!       4 pi ** 2 * k ** 2:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )

!       fcq_in_P1 #2:
        
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF

!       4 pi ** 2 * k * l:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )

!       fcq_in_P1 #3:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF

       
        IF ( refine_biso ) THEN

!           pi * i * k * s ** 2 / 2 term:
            temp_complex = CMPLX ( 0.0_wp, pi * mtz_1%Ko * mtz_1%s_in_P1 / 2.0_wp, KIND=wp ) 

!           fcq_in_P1 #4:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4) 

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,4) ) &
                          * mtz_1%wphases      
            ENDIF

        ENDIF

        IF ( refine_occ ) THEN

!           2 * pi * i * k / O(j):
            temp_complex = CMPLX ( 0.0_wp, -twopi * mtz_1%Ko, KIND=wp )

!           fcq_in_P1 #5:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,5) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF

        aniso_row2:IF ( refine_Us ) THEN

!           4 * pi ** 3 * i * h ** 2 * k (6):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Ko, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 3 * i * k ** 3 (7):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ko ** 3, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 3 * i * k * l ** 2 (8):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ko * mtz_1%Lo ** 2, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 3 * i * h * k ** 2 (9):
            temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko ** 2, KIND=wp )
            fcq = fcq + temp_complex *  mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!          8 * pi ** 3 * i * h * k * l (10):
           temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
           fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

           IF ( H2 ) THEN
               fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,10 ) ) &
                         * mtz_1%wphases
           ENDIF

!          8 * pi ** 3 * i * k ** 2 * l (MAX_BLOCK_SIZE):
           temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ko ** 2 * mtz_1%Lo, KIND=wp )
           fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

           IF ( H2 ) THEN
               fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                         * mtz_1%wphases
            ENDIF

        ENDIF aniso_row2

!       ---------------------------------------------------------------------
        CASE (3)

!       4 pi ** 2 * h * l :
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ho * mtz_1%Lo, 0.0_wp, KIND=wp )

!       fcq_in_P1 #1:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF

!       4 pi ** 2 * k * l:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )

!       fcq_in_P1 #2:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF

!       4 pi ** 2 * l ** 2:
        temp_complex = CMPLX ( 4.0_wp * pi_squared * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )

!       fcq_in_P1 #3:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)                   ! H1 term

        IF ( H2 ) THEN
            fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) )  &  ! H2 term
                      * mtz_1%wphases                                                 ! H2 term
        ENDIF

        IF ( refine_biso ) THEN

!           pi * i * l * s ** 2 / 2 term:
            temp_complex = CMPLX ( 0.0_wp, pi * mtz_1%Lo * mtz_1%s_in_P1 / 2.0_wp, KIND=wp )

!           fcq_in_P1 #4:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,4) ) &
                          * mtz_1%wphases     
            ENDIF
        ENDIF

        IF ( refine_occ ) THEN

!           2 * pi * i * l / O(j):
            temp_complex = CMPLX ( 0.0_wp, -twopi * mtz_1%Lo, KIND=wp )

!           fcq_in_P1 #5:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)           
            IF ( H2 ) THEN
                fcq = fcq  + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,5) ) &
                           * mtz_1%wphases
            ENDIF

        ENDIF

        aniso_row3:IF ( refine_Us ) THEN
            
!           4 * pi ** 3 * h ** 2 * l (6):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Lo, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 3 * k ** 2 * l (7):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Ko ** 2 * mtz_1%Lo, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF           

!           4 * pi ** 3 * l ** 3 (8):
            temp_complex = CMPLX ( 0.0_wp, 4.0_wp * pi ** 3 * mtz_1%Lo ** 3, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 3 * h * k * l (9):
            temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!          8 * pi ** 3 * h * l ** 2 (10):
           temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Lo ** 2, KIND=wp )
           fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

           IF ( H2 ) THEN
               fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,10) ) &
                         * mtz_1%wphases
           ENDIF

!          8 * pi ** 3 * k * l ** 2 (MAX_BLOCK_SIZE):
           temp_complex = CMPLX ( 0.0_wp, 8.0_wp * pi ** 3 * mtz_1%Ko * mtz_1%Lo ** 2, KIND=wp )
           fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)
           
           IF ( H2 ) THEN
               fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                         * mtz_1%wphases
           ENDIF

       ENDIF aniso_row3
        

!       --------------------------------------------------------------------
        CASE (4)

!       We may or may not refine XYZ during Biso refinement:
        IF ( refine_xyz ) THEN 

!           -pi * i * h * s ** 2 / 2:
            temp_complex = CMPLX ( 0.0_wp, -pi * mtz_1%Ho * mtz_1%s_in_P1 / 2.0_wp, KIND=wp )        

!           fcq_in_P1 #1:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                          * mtz_1%wphases
            ENDIF

!           -pi * i * k * s ** 2 / 2: 
            temp_complex = CMPLX ( 0.0_wp, -pi * mtz_1%Ko * mtz_1%s_in_P1 / 2.0_wp, KIND=wp )
        
!           fcq_in_P1 #2:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                          * mtz_1%wphases
            ENDIF

!           -pi * i * l * s ** 2 / 2: 
            temp_complex = CMPLX ( 0.0_wp, -pi * mtz_1%Lo * mtz_1%s_in_P1 / 2.0_wp, KIND=wp )

!           fcq_in_P1 #3:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF

!       s ** 4 / 16 (BB):
        temp_complex = CMPLX ( mtz_1%s_in_P1 ** 2 / 16.0_wp, 0.0_wp, KIND=wp )
 
!       fcq_in_P1 #4:
        fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
        IF ( H2 ) THEN
            fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,4) ) &
                      * mtz_1%wphases
        ENDIF

        IF ( refine_occ ) THEN

!           -s ** 2 / ( 4 * O(j) ):  
            temp_complex = CMPLX ( -mtz_1%s_in_P1 / 4.0_wp, 0.0_wp, KIND=wp )

!           fcq_in_P1 #5:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,5) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF

        aniso_row4:IF ( refine_Us ) THEN

!           pi ** 2  * h ** 2 * s ** 2 / 2 (6):
            temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Ho ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
!           BUG CORRECTED MAR 2009 BVS (was "7" here):
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF

!           pi ** 2 * k ** 2 * S ** 2 / 2 (7):
            temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Ko ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           pi ** 2 * l ** 2 * s ** 2 / 2 (8):
            temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Lo ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           pi ** 2 * h * k * s ** 2 (9):
            temp_complex = CMPLX ( pi ** 2 * mtz_1%Ho * mtz_1%Ko * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           pi ** 2 * h * l * s ** 2 (10):
            temp_complex = CMPLX ( pi ** 2 * mtz_1%Ho * mtz_1%Lo * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF            

!           pi ** 2 * k * l * s ** 2 (MAX_BLOCK_SIZE):
            temp_complex = CMPLX ( pi ** 2 * mtz_1%Ko * mtz_1%Lo * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF aniso_row4
!       --------------------------------------------------------------------
        CASE (5)

        IF ( refine_xyz ) THEN

!           2 * pi * i * h / O(i), OX term:
            temp_complex = CMPLX ( 0.0_wp, twopi * mtz_1%Ho, KIND=wp )

!           fcq_in_P1 #1:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                          * mtz_1%wphases
            ENDIF

!           2 * pi * i * k / O(i), OY term:
            temp_complex = CMPLX ( 0.0_wp, twopi * mtz_1%Ko, KIND=wp )

!           fcq_in_P1 #2:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                          * mtz_1%wphases
            ENDIF

!           2 * pi * i * l / O(i), OZ term:
            temp_complex = CMPLX ( 0.0_wp, twopi * mtz_1%Lo, KIND=wp )

!           fcq_in_P1 #3:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)

            IF ( H2 ) THEN
                fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                          * mtz_1%wphases
            ENDIF
        ENDIF

        IF ( refine_biso ) THEN

!           -s ** 2 / ( 4 * O(j) ), QB term:
            temp_complex = CMPLX ( -mtz_1%s_in_P1 / 4.0_wp, 0.0_wp, KIND=wp )

!           fcq_in_P1 #4:
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,4) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF 

!       1 / ( O(i) * O(j) ):
        fcq = fcq + mtz_1%fcq_in_P1(1:nrefl,5)

        IF ( H2 ) THEN
            fcq = fcq + CONJG ( mtz_1%fcq_in_P1(1:nrefl,5) ) &
                      * mtz_1%wphases
        ENDIF

        aniso_row5:IF ( refine_Us ) THEN

!           -2 * pi ** 2 * h ** 2 (6):
            temp_complex = CMPLX( -2.0_wp * pi ** 2 * mtz_1%Ho ** 2, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF

!           -2 * pi ** 2 * k ** 2 (7):
            temp_complex = CMPLX( -2.0_wp * pi ** 2 * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           -2 * pi ** 2 * l ** 2 (8):
            temp_complex = CMPLX( -2.0_wp * pi ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)
        
            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           -4 * pi ** 2 * h * k (9):
            temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ho * mtz_1%Ko, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF
            
!           -4 * pi ** 2 * h * l (10):
            temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ho * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)
        
            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           -4 * pi ** 2 * k * l (MAX_BLOCK_SIZE):
            temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq =fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

        ENDIF aniso_row5

!  
        CASE (6)


!           Row 6:
            IF ( refine_xyz ) THEN
            
!               -4 * pi ** 3 * h ** 3 (1):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ho ** 3, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                          * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * h ** 2 * k (2):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Ko, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * h ** 2 * l (3):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * h ** 2 * s ** 2 / 2 (4):
                temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Ho ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -2 * pi ** 2 * h ** 2 (5):
                temp_complex = CMPLX ( -2.0_wp * pi ** 2 * mtz_1%Ho ** 2, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           4 * pi ** 4 * h ** 4 (6):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ho ** 4, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           4 * pi ** 4 * h ** 2 * k ** 2 (7):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 4 * h ** 2 * l ** 2 (8):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h ** 3 * k (9):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 3 * mtz_1%Ko, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h ** 3 * l (10):
            temp_complex = CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 3 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h ** 2 * k * l (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

            CASE (7)
!           Row 7:

            IF ( refine_xyz ) THEN
            
!               -4 * pi ** 3 * h * k ** 2 (1):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                          * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * k ** 3  (2):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ko ** 3, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * k ** 2 * l (3):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ko ** 2 * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * k ** 2 * s ** 2 / 2 (4):
                temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Ko ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -2 * pi ** 2 * k ** 2 (5):
                temp_complex = CMPLX ( -2.0_wp * pi ** 2 * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           4 * pi ** 4 * h ** 2 * k ** 2 (6):
            temp_complex =  CMPLX( 4.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           4 * pi ** 4 * k ** 4 (7):
            temp_complex =  CMPLX( 4.0_wp * pi ** 4 * mtz_1%Ko ** 4, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 4 * k ** 2 * l ** 2 (8):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ko ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * k ** 3 (9):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * k ** 2 * l (10):
            temp_complex = CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 2 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * k ** 3 * l (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ko ** 3 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

            CASE (8)
!           Row 8 (7 without Biso):
            
            IF ( refine_xyz ) THEN
            
!               -4 * pi ** 3 * h * l ** 2 (1):
                temp_complex = CMPLX( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Lo ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                          * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * k * l ** 2 (2):  BUG CORRECTED (additional array mult eliminated)
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Ko * mtz_1%Lo ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -4 * pi ** 3 * l ** 3 (3):
                temp_complex = CMPLX ( 0.0_wp, -4.0_wp * pi ** 3 * mtz_1%Lo ** 3, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * l ** 2 * s ** 2 / 2 (4):
                temp_complex = CMPLX ( 0.5_wp * pi ** 2 * mtz_1%Lo ** 2 * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -2 * pi ** 2 * l ** 2 (5):
                temp_complex = CMPLX ( -2.0_wp * pi ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           4 * pi ** 4 * h ** 2 * l ** 2 (6):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           4 * pi ** 4 * k ** 2 * l ** 2 (7):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Ko ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           4 * pi ** 4 * l ** 4 (8):
            temp_complex =  CMPLX ( 4.0_wp * pi ** 4 * mtz_1%Lo ** 4, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * k * l ** 2  (9):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * l ** 3 (10):
            temp_complex = CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Lo ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * k * l ** 3 (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ko * mtz_1%Lo ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

            CASE(9)
!           Row 9:

            IF ( refine_xyz ) THEN
            
!               -8 * pi ** 3 * h ** 2 * k (1):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Ko, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * h * k ** 2 (2):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * h * k * l (3):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * h * k * s ** 2  (4):
                temp_complex = CMPLX ( pi ** 2 * mtz_1%Ho * mtz_1%Ko * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -4 * pi ** 2 * h * k (5):
                temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ho * mtz_1%Ko, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           8 * pi ** 4 * h ** 3 * k (6):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 3 * mtz_1%Ko, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           8 * pi ** 4 * h * k ** 3 (7):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * k * l ** 2 (8):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h ** 2 * k ** 2 (9):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h ** 2 * k * l (10):
            temp_complex = CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h * k ** 2 * l (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 2 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

!           Row 10:
            CASE(10)

            IF ( refine_xyz ) THEN
            
!               -8 * pi ** 3 * h ** 2 * l (1):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho ** 2 * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * h * k * l (2):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * h * l ** 2 (3):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Lo ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * h * l * s ** 2  (4):
                temp_complex = CMPLX ( pi ** 2 * mtz_1%Ho * mtz_1%Lo * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -4 * pi ** 2 * h * l (5):
                temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ho * mtz_1%Lo, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           8 * pi ** 4 * h ** 3 * l (6):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 3 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           8 * pi ** 4 * h * k ** 2 * l (7):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 2 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * h * l ** 3 (8):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Lo ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h ** 2 * k * l (9):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h ** 2 * l ** 2 (10):
            temp_complex = CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h * k * l ** 2 (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF

            CASE (MAX_BLOCK_SIZE)

!           Row MAX_BLOCK_SIZE:

            IF ( refine_xyz ) THEN
            
!               -8 * pi ** 3 * h * k * l (1):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,1)

                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,1) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * k ** 2 * l  (2):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ko ** 2 * mtz_1%Lo, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,2)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,2) ) &
                              * mtz_1%wphases
                ENDIF

!               -8 * pi ** 3 * k * l ** 2 (3):
                temp_complex = CMPLX ( 0.0_wp, -8.0_wp * pi ** 3 * mtz_1%Ko * mtz_1%Lo ** 2, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,3)
        
                IF ( H2 ) THEN
                    fcq = fcq - temp_complex * CONJG ( mtz_1%fcq_in_P1(1:nrefl,3) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_biso ) THEN

!               pi ** 2 * k * l * s ** 2  (4):
                temp_complex = CMPLX ( pi ** 2 * mtz_1%Ko * mtz_1%Lo * mtz_1%s_in_P1, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,4)
             
                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,4) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

            IF ( refine_occ ) THEN

!               -4 * pi ** 2 * k * l (5):
                temp_complex = CMPLX ( -4.0_wp * pi ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
                fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,5)

                IF ( H2 ) THEN
                    fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,5) ) &
                              * mtz_1%wphases
                ENDIF

            ENDIF

!           8 * pi ** 4 * h ** 2 * k *l  (6):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ho ** 2 * mtz_1%Ko * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,6)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,6) ) &
                          * mtz_1%wphases
            ENDIF
 
!           8 * pi ** 4 * k ** 3 * l (7):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ko ** 3 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,7)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,7) ) &
                          * mtz_1%wphases
            ENDIF

!           8 * pi ** 4 * k * l ** 3 (8):
            temp_complex =  CMPLX ( 8.0_wp * pi ** 4 * mtz_1%Ko * mtz_1%Lo ** 3, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,8)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,8) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h * k ** 2 * l  (9):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko ** 2 * mtz_1%Lo, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,9)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,9) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * h * k * l ** 2 (10):
            temp_complex = CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ho * mtz_1%Ko * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,10)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,10) ) &
                          * mtz_1%wphases
            ENDIF

!           16 * pi ** 4 * k ** 2 * l ** 2 (MAX_BLOCK_SIZE):
            temp_complex =  CMPLX ( 16.0_wp * pi ** 4 * mtz_1%Ko ** 2 * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            fcq = fcq + temp_complex * mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE)

            IF ( H2 ) THEN
                fcq = fcq + temp_complex * CONJG (  mtz_1%fcq_in_P1(1:nrefl,MAX_BLOCK_SIZE) ) &
                          * mtz_1%wphases
            ENDIF


        CASE DEFAULT

            WRITE(*,"(' MATRIX_VECTOR_COEFS_IN_P1> ','imap= ',I4)")  imap
            CALL die('Program is not ready for additional parameter refinement...', srname)

        END SELECT

!       Use weights:
!        fcq = mtz_1%weight_in_P1 * fcq
        fcq = mtz_1%weight_in_P1 * mtz_1%sc_in_P1 * fcq

!       Map based on these coefs will be used for convolution:
        fcq = (fft_scale * map_1%volume) * fcq * EXP ( mtz_1%b_scale * 0.25_wp * mtz_1%s_in_P1 ) 

!       Done.
    END SUBROUTINE Ax_matrix_vector_coefs

    SUBROUTINE tronrud_matrix_coefs_in_P1(fcq, mtz_1, map_1, wphases, H2, ai, bi, aj, bj)
! 
!
!       Purpose:
!       =======
!       Calculates reciprocal space coefficients for calculation of matrix diagonal elements.
!       Can be used for calculations of partial (mixed) second derivatives of X, Y, Z, B and Q
!
!       In general need to call it twice...
!       H2 = .FALSE. calculates H1 terms coefs
!       H2 = .TRUE.  calculates H2 terms coefs
!
!       Date:              Programmer:           History of changes:
!       ====               ==========            ==================
!       Oct 2007           B.Strokopytov         Original code based on Tronrud(1999) paper.
!       Apr 2010           B.Strokopytov         Scales array mtz_1%sc_in_P1 has been added.
!       Apr 2010           B.Strokopytov         Scales must be squared. Bug corrected. 
!       Note:
!       ====
!       According to results of D.Tronrud: 
!       a) H1 map has symmetry of corresponding Patterson group.
!       b) H2 map has symmetry of `squared space group'
!
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT)         :: fcq
        TYPE(mtz),                      INTENT(INOUT)         :: mtz_1
!       Just one map is sufficient:
        TYPE(map),                      INTENT(IN)            :: map_1
!       Temporary array:
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)            :: wphases
        LOGICAL,                        INTENT(IN)            :: H2
        REAL(KIND=wp), DIMENSION(:),    INTENT(IN), OPTIONAL  :: ai
        REAL(KIND=wp), DIMENSION(:),    INTENT(IN), OPTIONAL  :: bi
        REAL(KIND=wp), DIMENSION(:),    INTENT(IN), OPTIONAL  :: aj
        REAL(KIND=wp), DIMENSION(:),    INTENT(IN), OPTIONAL  :: bj
        
!       Local vars:
        REAL(KIND=wp)                                         :: fft_scale
!       Counters:
        INTEGER                                               :: i
        CHARACTER(LEN=32),                               SAVE :: srname = 'tronrud_matrix_coefs_in_P1'

        fft_scale = map_1%volume / ( map_1%nu * map_1%nv * map_1%nw )

!       Calculate maps based on weighting scheme chosen (see Tronrud, 1999):
        IF ( .NOT. H2 ) THEN
            fcq = CMPLX ( mtz_1%weight_in_P1 * mtz_1%sc_in_P1 ** 2, 0.0_wp, KIND=wp )
        ELSE
            fcq = -CMPLX ( mtz_1%weight_in_P1 * mtz_1%sc_in_P1 ** 2, 0.0_wp, KIND=wp ) * wphases
        ENDIF

!       Weights have been applied now usual scaling:
!       BUG CORRECTED mtz_1%b_scale replace with map_1%b_scale otherwise inconsistent results may arise:
        IF ( mtz_1%b_scale /= map_1%b_scale ) THEN
            WRITE(*,*) ' mtz_1%b_scale=', mtz_1%b_scale
            WRITE(*,*) ' map_1%b_scale=', map_1%b_scale
            CALL die('Programming error. Those values above should be equal. Inconsistent results may arise...',&
                     srname)
        ENDIF

!       Should we use 0.25 from the very beginning?: 
        fcq = (fft_scale * map_1%volume) * fcq * EXP ( (0.5_wp * mtz_1%b_scale) * mtz_1%s_in_P1 )

        IF ( PRESENT ( ai ) .AND. PRESENT ( aj ) ) THEN

!           Make smoothie:
            DO i = 1, SIZE ( fcq )
                fcq(i) = fcq(i) * SUM ( ai * EXP ( -0.25_wp * bi * mtz_1%s_in_P1(i) ) )
                fcq(i) = fcq(i) * SUM ( aj * EXP ( -0.25_wp * bj * mtz_1%s_in_P1(i) ) )
            ENDDO

!           THen this will not be necessary anymore:
            fcq = fcq * EXP ( (-0.25_wp * mtz_1%b_scale) * mtz_1%s_in_P1 )

        ENDIF

    END SUBROUTINE tronrud_matrix_coefs_in_P1

    SUBROUTINE matrix_diagonal_coefs_in_P1(mtz_1, map_1, wphases, H2)
! 
!
!       Purpose:
!       =======
!       Calculates reciprocal space coefficients for calculation of main diagonal.
!
!       Date:              Programmer:           History of changes:
!       ====               ==========            ==================
!       Oct 2007           B.Strokopytov         Original code based on Tronrud(1999) paper.
!
!       
        TYPE(mtz),                                   INTENT(INOUT) :: mtz_1
!       Just one map is sufficient:
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: map_1
!       Temporary array:
        COMPLEX(KIND=wp), DIMENSION(:),              INTENT(IN)    :: wphases
        LOGICAL,                                     INTENT(IN)    :: H2
!       Local vars:
        REAL(KIND=wp)                                              :: fft_scale
        INTEGER                                                    :: nrefl
!       Counters:
        INTEGER                                                    :: i
        CHARACTER(LEN=32),                               SAVE      :: srname = 'matrix_diagonal_coefs_in_P1'

        i = first_allocated_map ( map_1 )
        fft_scale = map_1(i)%volume / ( map_1(i)%nu * map_1(i)%nv * map_1(i)%nw )

        nrefl = SIZE ( mtz_1%fc_in_P1 )

        IF ( ALLOCATED ( map_1(1) ) ) THEN
            IF ( .NOT. H2 ) THEN
                mtz_1%fcq_in_P1(1:nrefl,1) =  CMPLX ( mtz_1%weight_in_P1 * 4.0_wp * pi_squared * mtz_1%Ho ** 2, 0.0_wp, KIND=wp )
            ELSE
                mtz_1%fcq_in_P1(1:nrefl,1) = -CMPLX ( mtz_1%weight_in_P1 * 4.0_wp * pi_squared * mtz_1%Ho ** 2, 0.0_wp, KIND=wp ) &
                                           *  wphases
            ENDIF
        ENDIF

        IF ( ALLOCATED ( map_1(2) ) ) THEN
            IF ( .NOT. H2 ) THEN
                mtz_1%fcq_in_P1(1:nrefl,2) =  CMPLX ( mtz_1%weight_in_P1 * 4.0_wp * pi_squared * mtz_1%Ko ** 2, 0.0_wp, KIND=wp )
            ELSE
                mtz_1%fcq_in_P1(1:nrefl,2) = -CMPLX ( mtz_1%weight_in_P1 * 4.0_wp * pi_squared * mtz_1%Ko ** 2, 0.0_wp, KIND=wp ) &
                                           *  wphases
            ENDIF
        ENDIF

        IF ( ALLOCATED ( map_1(3) ) ) THEN
            IF ( .NOT. H2 ) THEN
                mtz_1%fcq_in_P1(1:nrefl,3) = CMPLX ( mtz_1%weight_in_P1 * 4 * pi_squared * mtz_1%Lo ** 2, 0.0_wp, KIND=wp )
            ELSE
                mtz_1%fcq_in_P1(1:nrefl,3) = -CMPLX ( mtz_1%weight_in_P1 * 4.0_wp * pi_squared * mtz_1%Lo ** 2, 0.0_wp, KIND=wp ) &
                                           *  wphases
            ENDIF
        ENDIF

        IF ( ALLOCATED ( map_1(4) ) ) THEN
            IF ( .NOT. H2 ) THEN
                mtz_1%fcq_in_P1(1:nrefl,4) = CMPLX ( mtz_1%weight_in_P1 * mtz_1%s_in_P1 ** 2 / 16.0_wp, 0.0_wp, KIND=wp )
            ELSE
                mtz_1%fcq_in_P1(1:nrefl,4) = CMPLX ( mtz_1%weight_in_P1 * mtz_1%s_in_P1 ** 2 / 16.0_wp, 0.0_wp, KIND=wp ) &
                                           * wphases
            ENDIF
        ENDIF

        IF ( ALLOCATED ( map_1(5) ) ) THEN

            IF ( .NOT. H2 ) THEN
                mtz_1%fcq_in_P1(1:nrefl,5) = CMPLX ( mtz_1%weight_in_P1, 0.0_wp, KIND=wp )
            ELSE
                mtz_1%fcq_in_P1(1:nrefl,5) = CMPLX ( mtz_1%weight_in_P1, 0.0_wp, KIND=wp ) &
                                           * wphases
            ENDIF

        ENDIF

        DO i = 1, SIZE ( map_1 )
            IF ( mtz_1%b_scale /= map_1(i)%b_scale ) THEN
                WRITE(*,*) ' mtz_1%b_scale=', mtz_1%b_scale
                WRITE(*,*) ' map_1%b_scale=', map_1(i)%b_scale
                CALL die('Programming error. Those values above should be equal. Inconsistent results may arise...',&
                          srname)
            ENDIF

            IF ( ALLOCATED ( map_1(i) ) ) THEN
                mtz_1%fcq_in_P1(1:nrefl,i) = (fft_scale * map_1(i)%volume) * mtz_1%fcq_in_P1(1:nrefl,i) &
                                           * EXP ( (0.5_wp * mtz_1%b_scale) * mtz_1%s_in_P1 )
            ENDIF  
        ENDDO

   END SUBROUTINE matrix_diagonal_coefs_in_P1

   SUBROUTINE Bx_matrix_vector_coefs(fcq, mtz_1, map_1, wphases, imap)
!
!
!       Purpose:
!       =======
!       Coefficients for subsequent calculation of second derivatives during convolution step.
!       This might be slower than the original version but no problems with symmetry at all.
!
!       Note:
!       ====
!       All calculations MUST be done in P1 sp. group at this stage.
!       This is relatively cheap in reciprocal space.
!       Algorithm does not allow for transformation of these coefficients into coefs for TRUE sp.group.
!
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT) :: fcq
        TYPE(mtz),                      INTENT(INOUT) :: mtz_1
        TYPE(map),                      INTENT(IN)    :: map_1
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)    :: wphases
        INTEGER,                        INTENT(IN)    :: imap
!       Local vars:
        REAL(KIND=wp)                                 :: fft_scale
        INTEGER                                       :: nr
        CHARACTER(LEN=32),                       SAVE :: srname='Bx_matrix_vector_coefs'

        nr = SIZE ( mtz_1%fc_in_P1 )

!       Checkz:
        IF ( SIZE ( mtz_1%fo_in_P1 ) /= nr ) THEN
            CALL messag('Programming error. Inconsistent sizes for FC_IN_P1 and FO_IN_P1',&
                        srname)
        ENDIF

        fft_scale = map_1%volume / ( map_1%nu * map_1%nv * map_1%nw )

!       Note change of sign for last two cases:
        SELECT CASE ( imap )
            CASE (1)
                fcq(1:nr) = mtz_1%fcq_in_P1(1:nr,1) - CONJG ( mtz_1%fcq_in_P1(1:nr,1) ) * wphases
            CASE (2)
                fcq(1:nr) = mtz_1%fcq_in_P1(1:nr,2) - CONJG ( mtz_1%fcq_in_P1(1:nr,2) ) * wphases
            CASE (3)
                fcq(1:nr) = mtz_1%fcq_in_P1(1:nr,3) - CONJG ( mtz_1%fcq_in_P1(1:nr,3) ) * wphases
            CASE (4)
                fcq(1:nr) = mtz_1%fcq_in_P1(1:nr,4) + CONJG ( mtz_1%fcq_in_P1(1:nr,4) ) * wphases
            CASE (5)
                fcq(1:nr) = mtz_1%fcq_in_P1(1:nr,5) + CONJG ( mtz_1%fcq_in_P1(1:nr,5) ) * wphases
            CASE DEFAULT
                WRITE(*,"(' BX_MATRIX_VECTOR_COEFS> ','imap= ', I4)")  imap
                CALL die('Program is not ready for U''s refinement...', srname)
        END SELECT

!       Use weights:
        fcq = CMPLX ( mtz_1%weight_in_P1, 0.0_wp, KIND=wp) * fcq

!       Map based on these coefs will be used for convolution:
        fcq = (fft_scale * map_1%volume) * fcq * EXP ( mtz_1%b_scale * 0.25_wp * mtz_1%s_in_P1 )

!       Done.
   END SUBROUTINE Bx_matrix_vector_coefs

END MODULE coefs
