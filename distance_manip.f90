MODULE distance_manip
USE mtz_io
USE basic_symmetry_operations
IMPLICIT NONE
CONTAINS
    FUNCTION min_dist ( mtz_1, xyz, xyz_target, isym_start )

!       Purpose:
!       =======
!       Calculates minimum eucledian distance between xyz and xyz_target vectors
!       xyz and xyz_target MUST be in fractional coordinates
!
!       Date:        Programmer:          Description of changes:
!       ====         ==========           ======================
!       Nov 2003     B.Strokopytov        Original code
!       Nov 2005     B.Strokopytov        Bug in P1 space group removed
!
        REAL(KIND=wp)                      :: min_dist
        TYPE(mtz),    INTENT(IN)           :: mtz_1
        TYPE(vector), INTENT(IN)           :: xyz
        TYPE(vector), INTENT(IN)           :: xyz_target
        INTEGER,      INTENT(IN), OPTIONAL :: isym_start
!       Local variables:
        INTEGER, PARAMETER                 :: trmax = 2
        REAL(KIND=wp)                      :: r2
        REAL(KIND=wp)                      :: r2_min
        TYPE(vector)                       :: xyz_sym
        TYPE(vector)                       :: xyz_diff
        TYPE (vector_int)                  :: translation
        INTEGER                            :: nsym_start
!       Counters:
        INTEGER                            :: iu
        INTEGER                            :: iv
        INTEGER                            :: iw
        INTEGER                            :: m
        
        IF ( PRESENT( isym_start ) ) THEN
            nsym_start = isym_start
        ELSE
            nsym_start = 1
        ENDIF

!       P1 space group problem:
        IF ( nsym_start > mtz_1%sp_group%number_of_symops ) THEN 
            min_dist = 0.5_wp * MIN ( mtz_1%cell(1), mtz_1%cell(2), mtz_1%cell(3) ) 
            RETURN            
        ENDIF

        r2_min = 1.0E20_wp
        symmetry:DO m = nsym_start, mtz_1%sp_group%number_of_symops
            xyz_sym  = mtz_1%sp_group%SYM(m) * xyz 
            xyz_diff = xyz_sym - xyz_target           
            along_z:DO iw = -trmax, trmax
                along_y:DO iv = -trmax, trmax
                    along_x:DO iu = -trmax, trmax
                       translation = (/iu, iv, iw/)
                       r2 = mtz_1%REAL_TENSOR * ( xyz_diff + translation )
                       IF  ( r2 < r2_min ) THEN
                           r2_min = r2  
                       ENDIF
                    ENDDO along_x
                ENDDO along_y
            ENDDO along_z
        ENDDO symmetry
        min_dist = r2_min
    END FUNCTION min_dist

    FUNCTION min_dist_map ( sp_group, REAL_TENSOR, xyz, xyz_target, r_crit, isym_start )
!
!       Purpose:
!       =======
!       Calculates minimum eucledian distance between xyz and xyz_target vectors
!       using cheshirean operators...
!       xyz and xyz_target are assumed to be in fractional coordinates
!
!       Date:          Programmer:           Description of changes:
!       ====           ==========            ======================
!       Dec 2003       B. Strokopytov        Original code
!       Oct 2005       B. Strokopytov        map_1 removed from list of arguments
!                                            replaced with sp_group
!       Nov 2005       B. Strokopytov        Cosmestic changes 
!
        REAL(KIND=wp)                        :: min_dist_map
!        TYPE(map_double), INTENT(IN)         :: map_1
        TYPE(space_group)                    :: sp_group
        TYPE(tensor)                         :: REAL_TENSOR
        TYPE(vector),     INTENT(IN)         :: xyz
        TYPE(vector),     INTENT(IN)         :: xyz_target
        REAL(KIND=wp),    INTENT(IN)         :: r_crit
        INTEGER, INTENT(IN), OPTIONAL        :: isym_start
!       Local variables:
        INTEGER, PARAMETER                   :: trmax = 1
        REAL(KIND=wp)                        :: r2
        REAL(KIND=wp)                        :: r2_min
        REAL(KIND=wp)                        :: r2_crit
        TYPE (vector)                        :: xyz_sym
        TYPE (vector)                        :: xyz_diff
        TYPE (vector_int)                    :: translation
!       Counters:
        INTEGER                              :: iu
        INTEGER                              :: iv
        INTEGER                              :: iw
        INTEGER                              :: m

        INTEGER                              :: nsym_start

        IF ( PRESENT( isym_start ) ) THEN
            nsym_start = isym_start
        ELSE
            nsym_start = 1
        ENDIF
        r2_min = 10.0_wp ** 10

!       P1 space group problem:
        IF ( nsym_start > sp_group%number_of_symops ) THEN
!           This should be something else. Probably smallest cell dim divided by 2.
            min_dist_map = r2_min
            RETURN
        ENDIF

!       Checks:
        IF ( ALLOCATED ( sp_group%sym_cheshire ) ) THEN
            IF ( SIZE ( sp_group%sym_cheshire ) <  sp_group%number_of_symops ) THEN
                CALL die ( 'Number of cheshireian operators is smaller then number of symmetry operators.', &
                           'min_dist_map')
            ENDIF
        ELSE
            CALL die ( 'Programming error. Cheshirean operators has not been allocated.', &
                       'min_dist_map' )
        ENDIF

        r2_min  = 1.0E10_wp
        r2_crit = r_crit ** 2
        symmetry:DO m = nsym_start, SIZE(sp_group%sym_cheshire)
            IF ( .SYMA. sp_group%SYM_INT_cheshire(m) /= 1 ) CYCLE symmetry
            xyz_sym  = sp_group%SYM_cheshire(m) * xyz
            xyz_diff = xyz_sym - xyz_target
            along_z:DO iw = -trmax, trmax
                along_y:DO iv = -trmax, trmax
                    along_x:DO iu = -trmax, trmax
                       translation = (/iu, iv, iw/)
                       r2          = REAL_TENSOR * ( xyz_diff + translation )
                       IF  ( r2 < r2_min ) THEN
                           r2_min = r2

!                          No need to continue wasting CPU time. Goal is achieved.
                           IF ( r2_min < r2_crit ) EXIT symmetry
                       ENDIF
                    ENDDO along_x
                ENDDO along_y
            ENDDO along_z
        ENDDO symmetry
        min_dist_map = r2_min
    END FUNCTION min_dist_map

    FUNCTION min_dist_symop ( mtz_1, xyz, xyz_target, r_crit, isym_start )
!
!       Purpose:
!       =======
!       Calculats symmetry operator (including additional cell translation)
!       corresponding to  minimum eucledian distance between xyz and xyz_target vectors
!       xyz and xyz_target are assumed to be in fractional coordinates
!
!       This operator should be applied to molecule corresponding to xyz vector
!       on input of this routine
!       When this is done you can calculate rmsd between two set of coordinates
!       using fast_rmsd function (module pdb_manip).
!
        TYPE(symop)                               :: min_dist_symop
        TYPE(mtz), INTENT(IN)                     :: mtz_1
        TYPE(vector), INTENT(IN)                  :: xyz
        TYPE(vector), INTENT(IN)                  :: xyz_target
        REAL(KIND=wp), INTENT(IN)                 :: r_crit
        INTEGER, INTENT(IN), OPTIONAL             :: isym_start
!       Local vars:
        INTEGER, PARAMETER                        :: trmax = 2
        REAL(KIND=wp)                             :: r2
        REAL(KIND=wp)                             :: r2_min
        TYPE (vector)                             :: xyz_sym
        TYPE (vector)                             :: xyz_diff
        TYPE (vector_int)                         :: translation
        TYPE (symop)                              :: best_symop
        INTEGER                                   :: nsym_start
!
!       Counters
        INTEGER                                   :: iu
        INTEGER                                   :: iv
        INTEGER                                   :: iw
        INTEGER                                   :: m
!
        IF ( PRESENT( isym_start ) ) THEN
            nsym_start = isym_start
        ELSE
            nsym_start = 1
        ENDIF

!       P1 space group problem:
        IF ( nsym_start > mtz_1%sp_group%number_of_symops ) THEN
!           This should be something else. Probably smallest cell dim divided by 2:
            min_dist_symop = -10.0_wp
            RETURN
        ENDIF

        r2_min     = 1.0E10_wp
        best_symop = mtz_1%sp_group%SYM(nsym_start)
        symmetry:DO m = nsym_start, mtz_1%sp_group%number_of_symops
            xyz_sym  = mtz_1%sp_group%SYM(m) * xyz
            xyz_diff = xyz_sym - xyz_target
            along_z:DO iw = -trmax, trmax
                along_y:DO iv = -trmax, trmax
                    along_x:DO iu = -trmax, trmax
                       translation = (/iu, iv, iw/)
                       r2          = mtz_1%REAL_TENSOR * ( xyz_diff + translation )
                       IF  ( r2 < r2_min ) THEN
                           r2_min     = r2
                           best_symop = mtz_1%sp_group%SYM(m) + translation
                       ENDIF
                    ENDDO along_x
                ENDDO along_y
            ENDDO along_z
        ENDDO symmetry

        IF ( r2_min < r_crit * r_crit ) THEN
            min_dist_symop = best_symop
        ELSE

!           set error flag which means that the routine failed to find
!           a symmetry operator which brings to molecular centroids close enough
!           to each other:
            min_dist_symop = -10.0_wp
        ENDIF

    END FUNCTION min_dist_symop

    SUBROUTINE short_distance_vectors(xyz_close, xyz, xyz_target,  mtz_1, r2_cut_off, isym_start)
        TYPE(vector), DIMENSION(:), ALLOCATABLE          :: xyz_close
        TYPE(vector),               INTENT(IN)           :: xyz
        TYPE(vector),               INTENT(IN)           :: xyz_target
        TYPE(mtz),                  INTENT(IN)           :: mtz_1
        REAL(KIND=wp),              INTENT(IN)           :: r2_cut_off
        INTEGER,                    INTENT(IN), OPTIONAL :: isym_start
!       Local
        INTEGER, PARAMETER                               :: trmax = 2
        INTEGER                                          :: nsym_start
        TYPE(vector)                                     :: xyz_sym
        TYPE(vector)                                     :: xyz_diff
        TYPE(vector_int)                                 :: translation
        REAL(KIND=wp)                                    :: r2
!       Counters
        INTEGER                                          :: iu
        INTEGER                                          :: iv
        INTEGER                                          :: iw
        INTEGER                                          :: m
        INTEGER                                          :: run
        INTEGER                                          :: short_contact

        IF ( PRESENT(isym_start) ) THEN
            nsym_start = isym_start
        ELSE
            nsym_start = 1
        ENDIF

        DO run = 0, 1
            short_contact = 0
            symmetry:DO m = nsym_start, mtz_1%sp_group%number_of_symops
                xyz_sym  = mtz_1%sp_group%SYM(m) * xyz
                xyz_diff = xyz_sym - xyz_target
                along_z:DO iw = -trmax, trmax
                    along_y:DO iv = -trmax, trmax
                        along_x:DO iu = -trmax, trmax
                            translation = (/iu, iv, iw/)
                            r2 = mtz_1%REAL_TENSOR * ( xyz_diff + translation )
                            IF ( r2 <= r2_cut_off ) THEN
                                short_contact = short_contact + 1
                                IF ( run == 1 ) xyz_close(short_contact) = xyz_diff + translation
                            ENDIF
                        ENDDO along_x
                    ENDDO along_y
                ENDDO along_z
            ENDDO symmetry
            IF ( run == 0 ) CALL allocate_array(xyz_close, short_contact)
        ENDDO
    END SUBROUTINE short_distance_vectors

END MODULE distance_manip
