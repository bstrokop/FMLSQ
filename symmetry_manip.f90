MODULE symmetry_manip
!
!  Purpose:
!  =======
!  For reading of symop library + utilities
!  Intended for `space_group' type manipulations
!
!  Date:          Programmer:           History of changes:
!  ====           ==========            ==================
!  Oct 2003       B.Strokopytov         Original code
!  Sep 2005       B.Strokopytov         Module iso_varying_string added
!  Sep 2005       B.Strokopytov         `find_symop_in_library' modified
!  Oct 2005       B. Strokopytov        Moved symop allocation/deallocation routines to this module
!  Nov 2005          -"-                Added ALLOCATED and SIZE interfaces
!
USE basic_symmetry_operations
USE constants
USE fail
USE fft_util
USE parser_library
USE select_kinds
USE string_manip
USE util
USE vectors
IMPLICIT NONE

TYPE :: space_group
    INTEGER           :: space_group_number                        ! e.g., 20
    INTEGER           :: number_of_symops                          ! total number of symmetry operations, e.g., 4
    INTEGER           :: number_of_primitive_symops                ! number of primitive cell operators (nsymp)
    INTEGER           :: laue_group_number                         ! nlaue
!
!   character length chosen to simplify interaction with ccp4 package:
    CHARACTER(len=10) :: space_group_name                              ! true space group name
    CHARACTER(len=10) :: point_group_name                              ! corresponding point space group name, e.g., 'PG2'
    CHARACTER(len=32) :: syngony                                       ! e.g., 'MONOCLINIC' 
    CHARACTER(len=32) :: cif_space_group_name                          ! CIF space group?
    CHARACTER(len= 1) :: lattice_type                                  ! e.g., 'C'
    CHARACTER(len=10) :: laue_group_name                               ! e.g., '2/m'

    TYPE(symop_int),  DIMENSION(:), ALLOCATABLE :: sym_int           ! symmetry operators in integer form
    TYPE(symop),      DIMENSION(:), ALLOCATABLE :: sym               ! symmetry operators in wp form
    TYPE(matrix_int), DIMENSION(:), ALLOCATABLE :: sym_hkl           ! inverted symmetry operators for hkl operations
    INTEGER,          DIMENSION(:), ALLOCATABLE :: non_polar_axes    ! e.g., [1 3] in monoclinic sp. groups (Y-axis unique)
    LOGICAL                                     :: chiral            ! checks whether enantiomorphic sp. group exists...
    INTEGER                                     :: primitive_cheshire_symops_number
    TYPE(symop_int),  DIMENSION(:), ALLOCATABLE :: sym_int_cheshire  ! array of origins in symop_int form
    TYPE(symop),      DIMENSION(:), ALLOCATABLE :: sym_cheshire      ! mix of sym and sym_int_cheshire ( cheshirean sp. group)
    LOGICAL,          DIMENSION(3)              :: lpaxis
END TYPE

PRIVATE
PUBLIC :: ALLOCATED
PUBLIC :: allocate_array
PUBLIC :: ASSIGNMENT(=)
PUBLIC :: centric_and_epsilon_tests
PUBLIC :: centric_phase
PUBLIC :: deallocate_array
PUBLIC :: deallocate_space_group
PUBLIC :: equivalent_origins
PUBLIC :: find_symop_in_library
PUBLIC :: get_chiral_space_group 
PUBLIC :: get_polar_axes
PUBLIC :: gridsize
PUBLIC :: in_ASU
PUBLIC :: is_chiral
PUBLIC :: operator(*)
PUBLIC :: put_in_ASU
PUBLIC :: print_space_group_contents
PUBLIC :: SIZE
PUBLIC :: space_group

! Declare interface operators:
INTERFACE ALLOCATED
    MODULE PROCEDURE allocated_sp_group
END INTERFACE

INTERFACE ASSIGNMENT ( = )
    MODULE PROCEDURE copy_space_group
END INTERFACE

INTERFACE get_polar_axes
    MODULE PROCEDURE get_polar_axes_from_sp_group
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_1d_symop_int_array
    MODULE PROCEDURE allocate_1d_symop_array
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_1d_symop_int_array
    MODULE PROCEDURE deallocate_1d_symop_array
END INTERFACE

INTERFACE SIZE
    MODULE PROCEDURE size_sp_group
END INTERFACE

! Local module data:
INTEGER, DIMENSION(11), SAVE :: chiral_1 = (/76, 91, 92, 144, 151, 152, 169, 171, 178, 180, 212 /)
INTEGER, DIMENSION(11), SAVE :: chiral_2 = (/78, 95, 96, 145, 153, 154, 170, 172, 179, 181, 213 /)

CONTAINS
    FUNCTION allocated_sp_group ( sp_group )
        LOGICAL                       :: allocated_sp_group
        TYPE(space_group), INTENT(IN) :: sp_group
       
        allocated_sp_group =      ALLOCATED( sp_group%sym_int )  &
                             .OR. ALLOCATED ( sp_group%sym )     &
                             .OR. ALLOCATED ( sp_group%sym_hkl )

    END FUNCTION allocated_sp_group

    FUNCTION size_sp_group ( sp_group )
        INTEGER                       :: size_sp_group
        TYPE(space_group), INTENT(IN) :: sp_group

        IF ( ALLOCATED ( sp_group%sym_int ) ) THEN
            size_sp_group = SIZE ( sp_group%sym_int )
        ELSE IF ( ALLOCATED ( sp_group%sym ) ) THEN
            size_sp_group = SIZE ( sp_group%sym )
        ELSE IF ( ALLOCATED ( sp_group%sym_hkl ) ) THEN
            size_sp_group = SIZE ( sp_group%sym_hkl )
        ENDIF

    END FUNCTION size_sp_group

    FUNCTION is_chiral ( sp_group )
!      
!       Purpose:
!       =======
!       Determines whether current space group is chiral
!       Only standard space groups included
!
!       Date:           Programmer:         History of changes:
!       ====            ==========          ==================
!       Nov 2003        B. Strokopytov      Original code
!       Oct 2005        B. Strokopytov      Modified. `ANY' intrinsic added 
!
        LOGICAL                          :: is_chiral
        TYPE ( space_group ), INTENT(IN) :: sp_group
 
        is_chiral = ANY ( chiral_1 == sp_group%space_group_number ) .OR. &
                    ANY ( chiral_2 == sp_group%space_group_number )
    END FUNCTION is_chiral

    FUNCTION get_chiral_space_group ( sp_group )
!
!       Purpose:
!       =======
!       Reads chiral counterpart of current space group
!       Only standard space groups included
!
!       Date:           Programmer:         History of changes:
!       ====            ==========          ==================
!       Nov 2003        B. Strokopytov      Original code
!       Oct 2005        B. Strokopytov      chiral=.TRUE. flag added 
!
!
!
!
        TYPE ( space_group )             :: get_chiral_space_group
        TYPE ( space_group ), INTENT(IN) :: sp_group
!       Local variables:
        CHARACTER(len=240)               :: file_symop
!       Counters:
        INTEGER                          :: i

        IF ( is_chiral ( sp_group ) ) THEN
            DO i = 1, SIZE(chiral_1)
                IF ( chiral_1(i) == sp_group%space_group_number ) THEN
                    get_chiral_space_group%space_group_number = chiral_2(i)
                    EXIT
                ELSE IF ( chiral_2(i) == sp_group%space_group_number ) THEN
                    get_chiral_space_group%space_group_number = chiral_1(i)
                    EXIT
                ENDIF
            ENDDO

!           Read database of space groups:
            CALL ugtenv( 'SYMOP', file_symop )
            CALL find_symop_in_library ( get_chiral_space_group, file_symop )

!           flag this space group as chiral:
            get_chiral_space_group%chiral = .TRUE.
        ELSE
            get_chiral_space_group = sp_group
        ENDIF
        CALL print_space_group_contents ( get_chiral_space_group )
    END FUNCTION get_chiral_space_group

    SUBROUTINE equivalent_origins ( sp_group )
!
!       Purpose:
!       =======
!       Finds alternative/equivalent origins
!
!       Date:         Programmer:            History of changes:
!       ====          ==========             ==================
!       Dec 2003      B.Strokopytov          Lifted from A.A.Vagin
!       Jan 2003      B.Strokopytov          Modified and adapted for f95 compiler     
!
        TYPE (space_group), INTENT(INOUT)    :: sp_group
!       Local variables:
        TYPE (matrix_int)                    :: symd
        TYPE (vector_int)                    :: vecd
        TYPE (vector_int)                    :: twelve    ! this is Vagin's tsym_factor
        INTEGER, DIMENSION(6)                :: id
        INTEGER, DIMENSION(3)                :: is
        INTEGER, DIMENSION(3)                :: denom
        INTEGER, DIMENSION(3)                :: fract
        CHARACTER(len=60)                    :: line
!       Counters:
        INTEGER                              :: run
        INTEGER                              :: i
        INTEGER                              :: j 
        INTEGER                              :: m
        INTEGER                              :: n
        INTEGER                              :: k1
        INTEGER                              :: k2
        INTEGER                              :: k3
        INTEGER                              :: norig

        norig           =  1000000           ! suitable for P1
        sp_group%lpaxis = .TRUE.             ! set all axes as polar

!       Too many origins in P1 - no need to continue:
        IF ( sp_group%space_group_number /= 1 ) THEN

!           Non-P1 space group detected:

!           Possible origins along xyz axes according to Int. Tables -> 0, 1/2, 1/3, 2/3, 1/4, 3/4
            id     = (/  0, 6, 4, 8, 3, 9 /) 
            twelve = (/ 12, 12, 12 /)
            CALL get_polar_axes ( sp_group, sp_group%lpaxis )

!           Have to run twice to figure out how much memory should be allocated:
            DO run = 1, 2

!               Set origin counter:
                norig = 1
                IF ( run == 2 ) THEN

!                   Set first origin is at [0,0,0]:
                    sp_group%sym_int_cheshire(norig) = 1   ! Just unity operator

                    WRITE(*,'('' EQUIVALENT ORIGINS> '', ''['',A,'']'',3(1X,I2,''/'',I1))') TRIM(int_to_c(norig)),&
                    (0, 0, i = 1, 3)
                ENDIF

!               Loop along x, y and z directions:
                DO k1 = 1, 6
                    DO k2 = 1, 6
                        z:DO k3 = 1, 6

                            IF ( k1 * k2 * k3 == 1 ) CYCLE z

!                           Check this origin
                            is = (/ id(k1), id(k2), id(k3) /)

!                           No need to search for possible new origin along polar axis:
                            DO i = 1, 3
                                IF ( sp_group%lpaxis(i) .AND. is(i) /= 0 ) CYCLE z
                            ENDDO
     
                            DO m = 1, sp_group%number_of_primitive_symops - 1
                                DO n = m + 1, sp_group%number_of_primitive_symops
!                                   Conversion to matrix via matrix_part_int (symd is a matrix):
                                    symd = sp_group%sym_int(m) - sp_group%sym_int(n)
                                    vecd = symd * is         ! using matrix_int_times_array_int
                                    vecd = vecd .MMOD. twelve

!                                   Check for non-zero component presence:
                                    IF ( ( vecd .DOT. vecd ) /= 0 ) CYCLE z      ! This is not a possible origin
                                ENDDO
                            ENDDO

                            norig = norig + 1

                            IF ( run == 2 ) THEN

!                               Set up primitive/unit matrix ('cause we are looking for translation vector):
                                sp_group%sym_int_cheshire(norig) = 1
!                               Rescale integer translation vector if tsym_factor is different from 12:
                                sp_group%sym_int_cheshire(norig) = sp_group%sym_int_cheshire(norig) &
                                                                   .SYMV. ( is * ( tsym_factor / 12 ) )
!                               Got final cheshire operator here...

!                               *****************************************
!                               This loop is for needed just for pretty printing of possible origins
!                               Perhaps could be thrown out later Oct 2005 BVS
!                               *****************************************

                                fract = is        ! Note absence of any scaling
                                denom = 12

!                               Reduce 6/12 to 1/2, etc.:
                                DO i = 1, 3
                                    DO j = 8, 2, -1
                                        IF ( MOD ( fract(i), j ) == 0 .AND. MOD ( denom(i), j ) == 0 ) THEN
                                            fract(i) = fract(i)/j
                                            denom(i) = denom(i)/j
                                        ENDIF
                                    ENDDO
                                    IF ( fract(i) == 0 ) denom(i) = 0
                                ENDDO

!                               Would be nice to avoid printing of 0/0 later...
                                WRITE(*,'('' EQUIVALENT_ORIGINS> '', ''['',A,'']'',3(1X,I2,''/'',I1))') & 
                                TRIM(int_to_c(norig)), (fract(i), denom(i), i = 1, 3)
                            ENDIF

                        ENDDO z
                    ENDDO
                ENDDO

                IF ( run == 1 ) THEN
                    WRITE(*,'('' EQUIVALENT_ORIGINS> '', ''Number of equivalent origins for '', A,&
                             &'' space group= '',A)')  TRIM(sp_group%space_group_name), TRIM(int_to_c(norig) )
                    CALL allocate_array ( sp_group%sym_int_cheshire, norig, 'sym_int_cheshire' )
                ENDIF

            ENDDO

!           CALL set_tsym_factor ( sp_group%sym_int_cheshire, sp_group%sym_int(1)%tsym_factor )
            CALL set_tsym_factor ( sp_group%sym_int_cheshire, tsym_factor )

!           Print final result:
            DO i = 1, norig
                line = sp_group%sym_int_cheshire(i)
                WRITE(*,'('' EQUIVALENT_ORIGINS> '', ''Possible origin ['',A''] in form of operator='',&
                    &''('',A,'')'')') TRIM(int_to_c(i)), TRIM(line)
            ENDDO

!           Reallocate non_polar_axes array:
            CALL deallocate_array( sp_group%non_polar_axes )
 
!           This seems elegant (F90 PACK function):
            CALL allocate_array( sp_group%non_polar_axes, COUNT (.NOT.sp_group%lpaxis) )
            sp_group%non_polar_axes = PACK( (/1,2,3/), .NOT. sp_group%lpaxis )
            WRITE(*,'('' EQUIVALENT_ORIGINS> '', ''Non polar axes array = '', 3I3)') sp_group%non_polar_axes
        ELSE
!           P1 case:
            CALL allocate_array(sp_group%non_polar_axes, 0)
            CALL messag('All axes are polar. Setting non_polar_axes array to zero size.', 'equivalent_origins')
        ENDIF

!       Final message:
        IF ( sp_group%lpaxis(1) .AND. sp_group%lpaxis(2) .AND. sp_group%lpaxis(3) ) THEN
            CALL messag('This polar++ space group has infinite number of origins in A B C volume.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(1) .AND. sp_group%lpaxis(2) ) THEN
            CALL messag('This polar+ space group has origin anywhere in A B plane.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(1) .AND. sp_group%lpaxis(3) ) THEN
            CALL messag('This polar+ space group has origin anywhere in A C plane.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(2) .AND. sp_group%lpaxis(3) ) THEN
            CALL messag('This polar+ space group has origin anywhere in B C plane.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(1) ) THEN
            CALL messag('This is a polar space group : origin is not fixed along A axis.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(2) ) THEN
            CALL messag('This is a polar space group : origin is not fixed along B axis.', 'equivalent_origins')
        ELSE IF ( sp_group%lpaxis(3) ) THEN
            CALL messag('This is a polar space group : origin is not fixed along C axis.', 'equivalent_origins')
        ELSE
            CALL messag('This is non-polar space group.', 'equivalent_origins')
        ENDIF
        
    END SUBROUTINE equivalent_origins

    SUBROUTINE centric_and_epsilon_tests ( h, centre, sysabs, epsilon_factor, phi_centric, sp_group, multiplicity )
!
!
!       Purpose:
!       =======
!       Estimates various paparameters for reciprocal vector h.
!
!       Date:            Programmer:               Description of changes:
!       ====             ==========                ======================
!       Sep 2003         B.Strokopytov             F95 code based on K.Cowtan's code
!       Oct 2004         B.Strokopytov             Multiplicity of reflections added
!       Oct 2005         B.Strokopytov             Cosmetic changes
!       Mar 2007         B.Strokopytov             Minor cosmetic changes
!
        TYPE(vector_int), INTENT(IN)            :: h
        LOGICAL,          INTENT(OUT)           :: centre
        LOGICAL,          INTENT(OUT)           :: sysabs
        REAL(KIND=wp),    INTENT(OUT)           :: epsilon_factor
        REAL(KIND=wp),    INTENT(OUT)           :: phi_centric
        TYPE(space_group)                       :: sp_group
        REAL(KIND=wp),    INTENT(OUT), OPTIONAL :: multiplicity
!       Local varibales:
        TYPE(vector_int)                        :: symmetry_related_h
        REAL(KIND=wp)                           :: phi_int
!       Counters:
        INTEGER                                 :: i

!       Initialize:        
        centre         = .FALSE.            ! centric refl
        sysabs         = .FALSE.            ! systematically absent
        epsilon_factor = 1.0_wp             
        phi_centric    = 0.0_wp             ! phase

!       This was absent. BUG CORRECTED NOV 2007 BVS:
        IF ( PRESENT ( multiplicity ) ) THEN
            multiplicity   = 1.0_wp
        ENDIF

!       In P1 sp. group we are done already:
        IF ( sp_group%number_of_symops == 1 ) RETURN
!
!       epsilon_factor would be probably better computed using `number_of_primitive_symops'
!       but then problems with sysabs arises... BVS Dec 12 2003:
        DO i = 2, sp_group%number_of_symops
            symmetry_related_h = sp_group%SYM_HKL(i) * h
            IF ( symmetry_related_h == h ) THEN
                epsilon_factor = epsilon_factor + 1.0_wp

!               Bug below corrected Nov 2003 BVS. Kevin was right:
                IF ( COS ( twopi * ( h .DOT. sp_group%sym(i) ) ) <= 0.999_wp ) sysabs = .TRUE.
            ENDIF

            IF ( -symmetry_related_h == h ) THEN
                centre      = .TRUE.
                phi_centric =  h .DOT. sp_group%sym(i) 
                phi_int     =  REAL( INT(phi_centric), KIND=wp )
                IF ( phi_centric > phi_int ) THEN
                    phi_centric =  pi * ( phi_centric - phi_int )
                ELSE
                    phi_centric =  pi * ( phi_centric - phi_int + 1.0_wp )
                ENDIF
            ENDIF
        ENDDO

!       Small addition to Kevin's code to calculate multiplicity factors:
        multi:IF ( PRESENT ( multiplicity ) ) THEN

            multiplicity = 1.0_wp
            DO i = 2, SIZE(sp_group%SYM_HKL) 
                symmetry_related_h = sp_group%SYM_HKL(i) * h

                IF ( symmetry_related_h == h .OR. -symmetry_related_h == h  ) THEN

!                   This is true multiplicity of a reflection:
                    multiplicity = multiplicity + 1.0_wp
                ENDIF

            ENDDO

!           Duplicates of a reflection in a hemisphere:
            multiplicity = SIZE ( sp_group%SYM_HKL ) / multiplicity

        ENDIF multi

    END SUBROUTINE centric_and_epsilon_tests

    SUBROUTINE centric_phase(hkl, sp_group, centric_phase_value, iprint)
!
!       Purpose:
!       =======
!       Determines value of centric phase modulo 180.0 using simple test
!
!       Note:
!       ====
!       Centric phase value always positive and less than 180.0 degrees
!       Another possibility is centric_phase_value +(or -) 180.0
!
!
!       Date:         Programmer:        Description of changes:
!       ====          ==========         ======================
!       Oct 2005      Strokopytov B.     Original code
!
        TYPE(vector_int),  INTENT(IN)           :: hkl
        TYPE(space_group), INTENT(IN)           :: sp_group
        REAL(KIND=wp),     INTENT(OUT)          :: centric_phase_value
        INTEGER,           INTENT(IN), OPTIONAL :: iprint 
!       Local
        INTEGER                                 :: m
        TYPE(vector)                            :: xyz
        REAL(KIND=wp)                           :: asf
        REAL(KIND=wp)                           :: bsf
        TYPE(vector)                            :: xyz_test
        INTEGER, DIMENSION(3)                   :: hkl_print

!       Checkz:
        IF ( .NOT. ALLOCATED(sp_group%SYM) ) THEN
            CALL die('Programming error... Array sp_group%SYM has not been allocated.', 'centric_phase')
        ENDIF

!       Figure out centrosymmetric phase:
        xyz_test = (/SQRT(2.0_wp), SQRT(3.0_wp), SQRT(5.0_wp)/)
        asf      = 0.0_wp
        bsf      = 0.0_wp
        DO m = 1, SIZE(sp_group%SYM)
            xyz = sp_group%SYM(m) * xyz_test
            asf = asf + COS(twopi * ( hkl .DOT. xyz ) )
            bsf = bsf + SIN(twopi * ( hkl .DOT. xyz ) )
        ENDDO

        centric_phase_value = MOD( ((180.0_wp/pi) * ATAN2(Bsf,Asf) + 360.0_wp), 180.0_wp)

        IF ( PRESENT ( iprint ) ) THEN
            hkl_print = hkl
            WRITE(*,"(' CENTRIC_PHASE> ', 'Centric phase possibilities ', 3I4, 2X, F8.1, ' or ', F8.1)") &
            hkl_print, centric_phase_value, centric_phase_value - 180.0_wp
        ENDIF

    END SUBROUTINE centric_phase

    FUNCTION in_ASU ( hkl, nlaue )
!
!       Purpose:
!       =======
!       Determines whether reflection belongs to ASU or not
!
!       Date:             Programmer:                     History of changes:
!       ====              ==========================      ==================
!       Sep 2003          Strokopytov B. & K. Cowtan      Almost original code
!
!
        INTEGER                        :: in_ASU
        TYPE(vector_int),   INTENT(IN) :: hkl
        INTEGER,            INTENT(IN) :: nlaue
!       Local
        INTEGER, DIMENSION(3)          :: my_hkl                        
        INTEGER                        :: h
        INTEGER                        :: k
        INTEGER                        :: l
        INTEGER                        :: jsign
        INTEGER                        :: inASU
!
        jsign = 1
        inASU = 0

        DO jsign = 1, -1, -2
             my_hkl = jsign * hkl
             h = my_hkl(1)
             k = my_hkl(2)
             l = my_hkl(3)

            SELECT CASE ( nlaue )
!           1 bar ( alternate 1 )
            CASE (1)
                IF ( h > 0 .OR. ( h == 0 .AND. k > 0 ) .OR. ( h == 0 .AND. k == 0 .AND. l >= 0 ) ) inASU = jsign
!           1 bar ( alternate 2 )
            CASE (2)
                IF ( k > 0 .OR. ( k == 0 .AND. l > 0 ) .OR. ( k == 0 .AND. l == 0 .AND. h >= 0 ) ) inASU = jsign 
!           1 bar ( alternate 3 )
            CASE (3)
                IF ( l > 0 .OR. ( l == 0 .AND. h > 0 ) .OR. ( h == 0 .AND. l == 0 .AND. k >= 0 ) ) inASU = jsign 
!           2/m ( alternate 1 ) 
            CASE(4)
                IF ( ( k >= 0 .AND. l > 0 ) .OR. ( k >= 0 .AND. l == 0 .AND. h >= 0 ) ) inASU = jsign
!           2/m  ( alternate 2 )
            CASE(5)
                IF ( ( l >= 0 .AND. k > 0 ) .OR. ( l >= 0 .AND. k == 0 .AND. h >= 0 ) ) inASU = jsign 
!           pg222
            CASE(6)
                IF ( h >= 0 .AND. k >=0 .AND. l >= 0 ) inASU = jsign          
!           pg4 -> 4/m
            CASE(7)
                IF ( (l >= 0 .AND. h == 0 .AND. k >= 0) .OR. ( l >= 0 .AND. h >  0 .AND. k >  0 ) ) inASU = jsign 
!           pg422 -> 4/mmm
            CASE(8)
                IF ( h >= 0 .AND. k >=0 .AND. l >= 0  .AND. h >= k ) inASU = jsign
!           pg3 -> 3bar
            CASE(9)
                IF ( ( h >= 0 .AND. k > 0 ) .OR. ( h ==0 .AND. k == 0 .AND. l > 0 ) ) inASU = jsign
!           pg312 -> 3barM
            CASE(10)
                IF ( ( h >= 0 .AND. k == 0 .AND. l >= 0 ) .OR. &
                     ( h >=0 .AND. k > 0 .AND. h >= k ) ) inASU = jsign
!           pg321 -> 3bar1M
            CASE(11)
                IF ( ( h >= 0 .AND. k >= 0 .AND. h > k ) .OR. &
                     ( h >= 0 .AND. k >= 0 .AND. h == k .AND. l >= 0 ) ) inASU = jsign
!           pg6 -> 6/m 
            CASE(12)
                IF ( ( l >=0 .AND. h > 0 .AND. k > 0 ) .OR. ( l >= 0 .AND. h == 0 .AND. k >= 0 ) ) inASU = jsign
!           pg622 -> 6/mmm
            CASE(13)
                IF (  h >= 0 .AND. k >= 0 .AND. h >= k .AND. l >= 0 ) inASU = jsign
!           pg23 -> m3
            CASE(14)
                IF ( h >=0 .AND. k >= 0 .AND. l >= 0 ) THEN 
                    IF ( ( l == h .AND. k >= h ) .OR. ( l >  h .AND. k > h ) )  inASU = jsign
                ENDIF

!           pg432 -> m3m
            CASE(15)
                IF ( h >= 0 .AND. k >=0 .AND. l >= 0 .AND. k >= l .AND. l >= h ) inASU = jsign
          
            CASE DEFAULT
                WRITE(*,'('' IN_ASU> '',''Laue group number on input is '',I4)')  nlaue
                CALL die('No such laue group whatsoever.','in_ASU')  
            END SELECT

            IF ( inASU /= 0 ) THEN
                in_ASU = inASU
                RETURN
            ENDIF
        ENDDO
        in_ASU = 0
    END FUNCTION in_ASU

    SUBROUTINE copy_space_group ( sp_1, sp_2 )
!
!       Purpose:
!       =======
!       Makes deep copy of space_group structure
!
!       Date:       Programmer:        Description of changes:
!       ====        ==========         ======================
!       Oct 2003    Strokopytov B.     Original code
!       Oct 2005    Strokopytov B.     INTENT(OUT) changed to INTENT(INOUT) 
!                                      Otherwise test for allocation of sp_1 makes no sense
!
!       Oct 2005    Strokopytov B.            
!
        TYPE (space_group), INTENT(INOUT) :: sp_1
        TYPE (space_group), INTENT(IN)    :: sp_2
!       Local variables:
        CHARACTER(LEN=1)                  :: s

        sp_1%space_group_number               = sp_2%space_group_number
        sp_1%number_of_symops                 = sp_2%number_of_symops
        sp_1%number_of_primitive_symops       = sp_2%number_of_primitive_symops
        sp_1%space_group_name                 = sp_2%space_group_name
        sp_1%laue_group_number                = sp_2%laue_group_number
        sp_1%laue_group_name                  = sp_2%laue_group_name
        sp_1%point_group_name                 = sp_2%point_group_name
        sp_1%syngony                          = sp_2%syngony
        sp_1%cif_space_group_name             = sp_2%cif_space_group_name
        sp_1%lattice_type                     = sp_2%lattice_type
        sp_1%chiral                           = sp_2%chiral
        sp_1%lpaxis                           = sp_2%lpaxis
        sp_1%primitive_cheshire_symops_number = sp_2%primitive_cheshire_symops_number    
 
        s = 's'
        IF ( sp_1%number_of_symops == 1 ) s = ' '
        WRITE(*,"(' COPY_SPACE_GROUP> Allocating space for ', A,' symmetry operator', A)") &
        TRIM ( int_to_c ( sp_1%number_of_symops ) ), s

        IF ( .NOT. ALLOCATED ( sp_2%sym_int ) ) THEN
            CALL die ('Programming error: sym_int array has not been initialized.',&
                      'copy_space_group' )
        ELSE IF ( .NOT. ALLOCATED ( sp_2%sym_hkl ) ) THEN
            CALL die ('Programming error: sym_hkl array has not been initialized.',&
                      'copy_space_group' )
        ELSE IF ( .NOT. ALLOCATED ( sp_2%sym ) ) THEN
            CALL die ('Programming error: sym array has not been initialized.',&
                      'copy_space_group' )
        ENDIF

!       This checks below are not absolutely necessary but more safe,
!       simple allocation will do:
        IF ( .NOT. ALLOCATED ( sp_1%sym_int) ) THEN
            CALL allocate_array( sp_1%sym_int, sp_1%number_of_symops )
        ELSE IF ( SIZE ( sp_1%sym_int ) /= SIZE ( sp_2%sym_int ) ) THEN
            CALL allocate_array( sp_1%sym_int, sp_1%number_of_symops )
        ENDIF
       
        IF ( .NOT. ALLOCATED ( sp_1%sym_hkl ) ) THEN
            CALL allocate_array ( sp_1%sym_hkl, sp_1%number_of_symops )
        ELSE IF ( SIZE(sp_1%sym_hkl) /= SIZE(sp_2%sym_hkl) ) THEN
            CALL allocate_array ( sp_1%sym_hkl, sp_1%number_of_symops )
        ENDIF

        IF ( .NOT. ALLOCATED ( sp_1%sym ) ) THEN
            CALL allocate_array ( sp_1%sym, sp_1%number_of_symops )
        ELSE IF ( SIZE(sp_1%sym) /= SIZE(sp_2%sym) ) THEN
            CALL allocate_array ( sp_1%sym, sp_1%number_of_symops )
        ENDIF

        sp_1%sym_int = sp_2%sym_int
        sp_1%sym_hkl = sp_2%sym_hkl
        sp_1%sym     = sp_2%sym  

        IF ( ALLOCATED ( sp_2%sym_int_cheshire ) ) THEN
            CALL allocate_array(sp_1%sym_int_cheshire, SIZE(sp_2%sym_int_cheshire))
            sp_1%sym_int_cheshire = sp_2%sym_int_cheshire
        ENDIF
      
        IF ( ALLOCATED ( sp_2%sym_cheshire ) ) THEN
            CALL allocate_array(sp_1%sym_cheshire, SIZE(sp_2%sym_cheshire))
            sp_1%sym_cheshire = sp_2%sym_cheshire
        ENDIF

        IF ( ALLOCATED ( sp_2%non_polar_axes ) ) THEN
            CALL allocate_array(sp_1%non_polar_axes, SIZE(sp_2%non_polar_axes))
            sp_1%non_polar_axes = sp_2%non_polar_axes
        ENDIF

    END SUBROUTINE copy_space_group

    SUBROUTINE get_laue_group ( sp_gr )
!
!       Purpose:
!       =======
!       Choose Laue group from PG name.
!
!       On entry:
!       ========
!       sp_gr%point_group_name   (Int Tab A or `CCP4' (from PGDEFN))
!
!       On exit:
!       =======
!       sp_gr%laue_group_number     Laue group number
!       sp_gr%laue_group_name       LAUNAM    Laue group name (Int Tab A)
!
!
!
!       Date:         Programmer:         Description of changes:
!       ====          ==========          ======================
!       Oct 2003      Strokopytov B.      Original code
!       Oct 2005      Strokopytov B.      PACK intrinsic introduced
!       Oct 2005      Strokopytov B.      Cosmetic changes
!
!
!
        TYPE(space_group), INTENT(INOUT)           :: sp_gr
!       Local parameters:
        INTEGER, PARAMETER                         :: laue_total = 42
        CHARACTER(LEN=LEN(sp_gr%point_group_name)) :: local_name
        CHARACTER(LEN=6), DIMENSION(laue_total)    :: point_group_name
        CHARACTER(LEN=6), DIMENSION(laue_total)    :: laue_group_name
        INTEGER,          DIMENSION(laue_total)    :: laue_group_number
!       For PACK intrinsic:
        INTEGER,          DIMENSION(1)             :: vector
!       Counters:
        INTEGER                                    :: i
!
        point_group_name   =                                                              &
                             (/                                                           &
                               '23    ', 'M3BAR ',                                        &
                               '432   ', '4BAR3M', 'M3BARM',                              &
                               '422   ', '4/MMM ', '4BAR2M', '4BARM2',                    &
                               '321   ', '32    ', '3M1   ', '3BARM1', '3BARM ','3M    ', &
                               '622   ', '6MM   ', '6BAR2M', '6BARM2', '6/MMM ',          &
                               '312   ', '31M   ', '3BAR1M',                              &
                               '222   ', 'MMM   ', 'MM2   ', '2MM   ', 'M2M   ',          &
                               '2C    ',                                                  &
                               '2     ', 'M     ', '2/M   ',                              &
                               '1     ', '1BAR  ',                                        &
                               '3     ', '3BAR  ',                                        &
                               '4     ', '4/M   ', '4BAR  ',                              &
                               '6     ', '6/M   ', '6BAR  '                               &
                            /)

        laue_group_number  = (/                                                           &
                                 14,       14,                                            &
                                 15,       15,       15,                                  &
                                  8,        8,        8,        8,                        &
                                 11,       11,       11,       11,       11,       11,    &
                                 13,       13,       13,       13,       13,              &
                                 10,       10,       10,                                  &
                                  6,        6,        6,        6,        6,              &
                                  5,                                                      &
                                  4,        4,        4,                                  &
                                  3,        3,                                            &
                                  9,        9,                                            &
                                  7,        7,        7,                                  &
                                 12,       12,       12                                   &
                             /)    

        laue_group_name    = (/                                                           &
                              'm3bar ', 'm3bar ',                                         &
                              'm3barm', 'm3barm', 'm3barm',                               &
                              '4/mmm ', '4/mmm ', '4/mmm ', '4/mmm ',                     &
                              '3barm ', '3barm ', '3barm ', '3barm ', '3barm ', '3barm ', &
                              '6/mmm ', '6/mmm ', '6/mmm ', '6/mmm ', '6/mmm ',           &
                              '3bar1m', '3bar1m', '3bar1m',                               &
                              'mmm   ', 'mmm   ', 'mmm   ', 'mmm   ', 'mmm   ',           &
                              '2/m   ',                                                   &
                              '2/m   ', '2/m   ', '2/m   ',                               &
                              '-1    ', '-1    ',                                         &
                              '-3    ', '-3    ',                                         &
                              '4/m   ', '4/m   ', '4/m   ',                               &
                              '6/m   ', '6/m   ', '6/m   '                                &
                             /)

!
      sp_gr%laue_group_number = 0
      local_name              = sp_gr%point_group_name
      CALL ucase ( local_name )

!     Strip off 'PG' if present:
      IF ( local_name(1:2) == 'PG' ) local_name = local_name(3:)

      IF ( COUNT ( local_name == point_group_name )  == 0 ) THEN
!         Have not found anything useful:
          WRITE(*,*) ' Incorrect point group name :', local_name
          WRITE(*,*) ' Possible names:'
          WRITE(*,'(8(1X,A6))') point_group_name
          CALL die('Can not find suitable laue group. Check point group name.', 'get_laue_group')
      ELSE IF ( COUNT ( local_name == point_group_name )  /= 1 ) THEN
!         Number of matches > 1:
          CALL die('Programming error. Duplicates in point_group_name array.',  'get_laue_group')
      ELSE
!         Single unambiguous match:
          vector = PACK( (/(i, i = 1, laue_total)/), local_name == point_group_name )
          sp_gr%laue_group_number = laue_group_number(vector(1))
          sp_gr%laue_group_name   = laue_group_name  (vector(1))
      ENDIF
    END SUBROUTINE get_laue_group

    SUBROUTINE find_symop_in_library(sp_gr, file_1)
!
!       Purpose:
!       =======
!       Reads space group from the library and runs various tests
!
!       Input:
!       =====
!       Either sp_gr%space_group_number or sp_gr%space_group_name 
!       must be defined before calling this routine
!
!       Output:
!       ======
!       Almost all thinkable attributes of sp_gr
!
!       
!       Date:         Programmer:         Description of changes:
!       ====          ==========          ======================
!       Sep 2003      Strokopytov B.      Original code
!       Sep 2005      Strokopytov B.      Some checking routines moved to module `basic_symmetry_manip'
!       Oct 2005      Strokopytov B.      Cosmetic changes (e.g., euclid-> cheshire)
!
!
        TYPE (space_group), INTENT(INOUT)              :: sp_gr
        CHARACTER(LEN=*),   INTENT(IN)                 :: file_1
!       Local variables:
        INTEGER                                        :: space_group_number
        CHARACTER(LEN=LEN(sp_gr%space_group_name))     :: space_group_name
        INTEGER                                        :: istat
!       Counters:
        INTEGER                                        :: i
        INTEGER                                        :: j
        INTEGER                                        :: k
!       Needed to read lines:
        CHARACTER(LEN=80)                              :: line
        INTEGER                                        :: nwords
        CHARACTER(LEN=wlen), DIMENSION(:), ALLOCATABLE :: words
!
        LOGICAL                                        :: ok
        CHARACTER(len=1)                               :: s
        CHARACTER(len=60)                              :: sym_i
        TYPE(symop_int),     DIMENSION(:), ALLOCATABLE :: new_int_cheshire
        INTEGER                                        :: isymop

        CALL messag(' ', 'find_symop_in_library')

!       Check subroutine input:
        IF ( sp_gr%space_group_number <= 0 .AND. LEN_TRIM(sp_gr%space_group_name) <= 0 ) THEN

            WRITE(*,"(' FIND_SYMOP_IN_LIBRARY> ', 'space group number= ', I4)") sp_gr%space_group_number
            WRITE(*,"(' FIND_SYMOP_IN_LIBRARY> ', 'space group name  = ', A)")  sp_gr%space_group_name
            CALL die ('Incorrect usage: either space group name or number must be present.','find_symop_in_library')

        ELSE IF ( sp_gr%space_group_number <= 0 .AND. sp_gr%space_group_name == '') THEN

            WRITE(*,"(' FIND_SYMOP_IN_LIBRARY> ', 'space group number= ',I4)") sp_gr%space_group_number
            WRITE(*,"(' FIND_SYMOP_IN_LIBRARY> ', 'space group name  = ',A)")  sp_gr%space_group_name
            CALL die ('Incorrect usage: either space group name or number must be present.',&
                      'find_symop_in_library')

        ELSE IF ( sp_gr%space_group_number > 0 ) THEN

            sp_gr%space_group_name=' '

        ENDIF

!       No need to invent numbers:
        isymop = get_next_io_unit()

!       Open symop library:
        OPEN(UNIT=isymop, FILE=file_1, ACCESS='SEQUENTIAL', FORM='FORMATTED', &
             ACTION='READ', STATUS='OLD', IOSTAT = istat )
        IF ( istat /= 0 ) CALL die('Failed to open symmetry operators library. Check names please:'//file_1, &
                                   'find_symop_in_library')

!       Initialize necessary constants:
        ok    = .FALSE.
        istat = 0
        DO WHILE( istat==0 )
             READ(UNIT=isymop, FMT='(A)', IOSTAT = istat ) line
             IF ( DEBUG > 30 ) WRITE(*,*) line
             CALL mysplit( line, words)
             nwords = SIZE(words)
             IF ( nwords >= 7 ) THEN

                 IF (debug > 20 ) WRITE(*,*) line
                 READ(words(1), * ) space_group_number
                 READ(words(2), * ) sp_gr%number_of_symops

                 IF ( sp_gr%number_of_symops > 192 ) THEN
                     CALL die ('Too many symmetry operators in this sp.group: '//TRIM(line), 'find_symop_in_library')
                 ELSE
                     line = words(3)
                     READ( line, * ) sp_gr%number_of_primitive_symops
                     space_group_name = words(4)
                 ENDIF

                 IF ( LEN_TRIM(sp_gr%space_group_name) > 0 ) THEN
                     ok =  TRIM(words(4)) ==  TRIM(sp_gr%space_group_name) 
                 ELSE IF ( sp_gr%space_group_number > 0 ) THEN
                     ok = space_group_number == sp_gr%space_group_number 
                 ENDIF

                 IF ( ok ) THEN
!                    That is it:
                     sp_gr%space_group_name     = space_group_name
                     sp_gr%space_group_number   = space_group_number
                     sp_gr%point_group_name     = words(5)
                     sp_gr%syngony              = words(6)
                     sp_gr%cif_space_group_name = words(7)
                     sp_gr%lattice_type         = sp_gr%space_group_name(1:1)                    

!                    Lets suddenly allocate memory for symm operators:
                     s = 's'
                     IF ( sp_gr%number_of_symops == 1 ) s=' '
                     WRITE(*,"(' FIND_SYMOP_IN_LIBRARY> ','Allocating space for ', A, &
                              &' symmetry operator', A)") TRIM ( int_to_c ( sp_gr%number_of_symops ) ), s

                     CALL allocate_array ( sp_gr%sym_int, sp_gr%number_of_symops, 'sym_int' )
                     CALL allocate_array ( sp_gr%sym,     sp_gr%number_of_symops, 'sym' )
                     CALL allocate_array ( sp_gr%sym_hkl, sp_gr%number_of_symops, 'sym_hkl' )

                     DO i = 1, sp_gr%number_of_symops
                         READ (UNIT=isymop, FMT='(A)', IOSTAT=istat) line
                         IF ( istat /= 0 ) THEN
                             CALL die ('Failed to read line containing symmetry operator.', &
                                       'find_symop_in_library')
                         ENDIF
                         sp_gr%sym_int(i) = line
                     ENDDO

                     EXIT

                 ENDIF
             ENDIF
        ENDDO

!       Deallocate words array:
        IF ( ALLOCATED ( words ) ) DEALLOCATE ( words )

!       Print result:
        IF ( ok ) THEN
            CALL messag ('Successfully finished reading the library of symmetry operators.',&
                         'find_symop_in_library' )
            CLOSE ( 11, IOSTAT = istat )
            IF ( istat /= 0 ) CALL die ('Failed to close SYMOP file','find_symop_in_library')
        ELSE
            IF ( LEN_TRIM(sp_gr%space_group_name) > 0 ) WRITE(*,*) ' Space group name   : ', sp_gr%space_group_name
            IF ( sp_gr%space_group_number > 0  )     WRITE(*,*) ' Space group number : ', sp_gr%space_group_number
            CALL die('Failed to find anything in library of symmetry operators.','find_symop_in_library')
        ENDIF

!       Check if the first symop is identity:
        IF ( sp_gr%sym_int(1) /= 1 ) THEN
            CALL die (' First symmetry operator is not unity.','find_symop_in_library')
        ENDIF

!       Check for duplicates:
        CALL check_duplicate_symops ( sp_gr%sym_int )

!       Major test for closed group:
        CALL check_symop_times_symop_matrix ( sp_gr%sym_int )

        DO i = 1, sp_gr%number_of_symops

!           Get matrix part of the operator via overloaded matrix_part_int:
            sp_gr%sym_hkl(i) = sp_gr%sym_int(i) 

!           Check determinant:
            IF ( .DET.( sp_gr%sym_hkl(i) ) == -1 ) THEN
                sym_i = sp_gr%sym_int(i)
                CALL warn('CENTRIC space group operator :'//TRIM(sym_i), 'find_symop_in_library')
            ENDIF
!
!           get real operator from integer: see MODULE basic_symmetry_operations
            sp_gr%sym(i)     = sp_gr%sym_int(i)
!
!           get transpose of the inverse matrix: this will be the symmetry operation for hkl
!           This produces problems with phase shifts and should be thoroughly avoided. BVS FEB 2004
!           sp_gr%sym_hkl(i) =  sp_gr%sym_hkl(i) **(-1)

!           whether we need an inverse should be thorougly checked in e.g.hexagonal sp. group!!!
            sp_gr%SYM_HKL(i) =  TRANSPOSE ( sp_gr%sym_hkl(i) )
        ENDDO

        CALL get_laue_group ( sp_gr )

!       Let's continue our games for T3 translation function
        CALL equivalent_origins ( sp_gr )
        IF ( .NOT. ALLOCATED ( sp_gr%sym_int_cheshire ) .OR. sp_gr%space_group_number == 1 ) THEN
           CALL messag('P1 space group encountered: returning to main program.', 'find_symop_in_library')
           RETURN
        ENDIF

!       Check first symmetry operator: it must be X,Y,Z
        IF ( sp_gr%sym_int_cheshire(1) /= 1 ) THEN
            CALL die('Programming error: 1st Cheshire symmetry operator is not identity operator.', &
                     'find_symop_in_library')
        ENDIF

!       Form a new Cheshire space group:
        sp_gr%primitive_cheshire_symops_number = SIZE(sp_gr%sym_int_cheshire)
        CALL allocate_array ( new_int_cheshire, sp_gr%number_of_symops * SIZE(sp_gr%sym_int_cheshire) ) 
!        CALL allocate_array(sp_gr%sym_int_cheshire, sp_gr%number_of_symops * SIZE(sp_gr%sym_int_cheshire) , 'sym_int_cheshire')
        k = 0
        DO i = 1, sp_gr%number_of_symops
            DO j = 1, SIZE( sp_gr%sym_int_cheshire)
                k = k + 1
                new_int_cheshire(k) = sp_gr%sym_int(i) * sp_gr%sym_int_cheshire(j)
            ENDDO
        ENDDO
        CALL eliminate_duplicate_symops ( new_int_cheshire )
        CALL allocate_array ( sp_gr%sym_int_cheshire, SIZE ( new_int_cheshire ) )
        sp_gr%sym_int_cheshire = new_int_cheshire
!       Free memory:
        CALL deallocate_array ( new_int_cheshire )

        CALL messag(' ', 'find_symop_in_library')

        IF ( debug > 25 ) THEN

!           A small bug corrected OCT 2007 BVS -> "k" must be reinitialized:
            k = SIZE ( new_int_cheshire )
            CALL messag('Resulting '//TRIM ( int_to_c ( k ) )//' operators in cheshirean space group:',&
                        'find_symop_in_library')

            DO i = 1, SIZE ( sp_gr%sym_int_cheshire )
                sym_i = sp_gr%sym_int_cheshire(i)
                CALL messag('Operator #'//TRIM(int_to_c(i))//'=('//TRIM(sym_i)//')', 'find_symop_in_library')
            ENDDO
        ENDIF

!       Perform checks again:      
        CALL check_duplicate_symops(sp_gr%sym_int_cheshire)
        CALL check_symop_times_symop_matrix(sp_gr%sym_int_cheshire)

!       Create set of real(wp) operators:
        CALL allocate_array(sp_gr%sym_cheshire, SIZE(sp_gr%sym_int_cheshire), 'sym_cheshire')

!       Get wp symmetry operators:  
        DO i = 1, SIZE(sp_gr%sym_int_cheshire)
            sp_gr%sym_cheshire(i) = sp_gr%sym_int_cheshire(i)
        ENDDO

        CALL messag('Done.', 'find_symop_in_library')
        CALL messag(' ', 'find_symop_in_library')

    END SUBROUTINE find_symop_in_library

    SUBROUTINE print_space_group_contents(sp_gr)
        TYPE(space_group), INTENT(IN) :: sp_gr
!       Local variables:
        CHARACTER(LEN=60)             :: sym_xyz
        CHARACTER(LEN=60)             :: sym_hkl
        CHARACTER(LEN=256)            :: line
!       Counters:
        INTEGER                       :: i
!
        CALL messag(' ', 'print_space_group_contents')        
        CALL messag('Space group number=                      '//TRIM(int_to_c(sp_gr%space_group_number)),&
                    'print_space_group_contents')
        CALL messag('Total number of symmetry operations=     '//TRIM(int_to_c(sp_gr%number_of_symops)),&
                    'print_space_group_contents')
        CALL messag('Number of primitive symmetry operations= '//TRIM(int_to_c(sp_gr%number_of_primitive_symops)),&
                    'print_space_group_contents')

        CALL messag('Space group name=     '//TRIM(sp_gr%space_group_name),     'print_space_group_contents')
        CALL messag('Point group name=     '//TRIM(sp_gr%point_group_name),     'print_space_group_contents')
        CALL messag('Syngony=              '//TRIM(sp_gr%syngony),              'print_space_group_contents')
        CALL messag('CIF space group name= '//TRIM(sp_gr%cif_space_group_name), 'print_space_group_contents')
        CALL messag('Lattice_type=         '//TRIM(sp_gr%lattice_type),         'print_space_group_contents')
        CALL messag('Laue group number=    '//TRIM(int_to_c(sp_gr%laue_group_number)),&
                    'print_space_group_contents')
        CALL messag('Laue group name=      '//TRIM(sp_gr%laue_group_name),      'print_space_group_contents')
        CALL messag (' ','print_space_group_contents')
        CALL messag (' Real space operations                  Reciprocal space operations', &
                     'print_space_group_contents')

        DO i = 1, sp_gr%number_of_symops
            sym_xyz = sp_gr%sym_int(i)
            sym_hkl = sp_gr%sym_hkl(i)
            WRITE( line, '(7X, A40, 4X, A40)' ) sym_xyz, sym_hkl
            CALL messag( TRIM ( line ), 'print_space_group_contents')
        ENDDO

        IF ( sp_gr%chiral ) THEN 
            WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'This is chiral space group.')")
        ELSE
            WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'This space group is not chiral.')")
        ENDIF 

        IF ( ALLOCATED ( sp_gr%non_polar_axes ) ) THEN
            i = COUNT(sp_gr%lpaxis)
            IF ( i == 1 ) THEN
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'It has ', I1, ' polar axis.')") i
            ELSE IF ( i > 1 ) THEN
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'It has ', I1, ' polar axes.')") i
            ELSE
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'It has no polar axes.')")
            ENDIF
        ELSE
            WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'It has 3 polar axes.')") 
        ENDIF

        IF ( ALLOCATED ( sp_gr%sym_int_cheshire ) ) THEN

            IF ( ALLOCATED ( sp_gr%sym_int ) ) THEN
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'It has ', A, ' possible origins.')") &
                TRIM ( int_to_c ( SIZE ( sp_gr%sym_int_cheshire ) / SIZE ( sp_gr%sym_int ) ) )
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', '... and ', A, ' primitive cheshire operators.')") &
                TRIM ( int_to_c ( sp_gr%primitive_cheshire_symops_number ) )
            ENDIF

            WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', ' It has ', A, ' cheshirean operators.')") &
            TRIM ( int_to_c ( SIZE ( sp_gr%sym_int_cheshire ) ) )

            DO i = 1, SIZE ( sp_gr%sym_int_cheshire )
                line = sp_gr%sym_int_cheshire(i)
                WRITE(*,"(' PRINT_SPACE_GROUP_CONTENTS> ', 'cheshire operator= (', A, ')')")&
                TRIM ( line )
            ENDDO
        ENDIF

    END SUBROUTINE print_space_group_contents
   
    SUBROUTINE deallocate_space_group ( sp_gr )
!
!       Purpose:
!       =======
!       Frees memory for sp_gr object
!
!       Date:          Programmer:           History of changes:
!       ====           ==========            ==================
!       Oct 2003       Strokopytov B.        Original code
!       Oct 2005          -"-                Less printing plans
!
        TYPE(space_group), INTENT(INOUT) :: sp_gr

        IF ( ALLOCATED(sp_gr%sym_int)        )   CALL deallocate_array ( sp_gr%sym_int,        'sym_int' )
        IF ( ALLOCATED(sp_gr%sym)            )   CALL deallocate_array ( sp_gr%sym,            'sym'     )
        IF ( ALLOCATED(sp_gr%sym_hkl)        )   CALL deallocate_array ( sp_gr%sym_hkl,        'sym_hkl' )

!        IF ( ALLOCATED(sp_gr%non_polar_axes)   ) CALL deallocate_array ( sp_gr%non_polar_axes,   'non_polar_axes'   )
        IF ( ALLOCATED(sp_gr%non_polar_axes)   ) CALL deallocate_array ( sp_gr%non_polar_axes  )
        IF ( ALLOCATED(sp_gr%sym_int_cheshire) ) CALL deallocate_array ( sp_gr%sym_int_cheshire, 'sym_int_cheshire' )
        IF ( ALLOCATED(sp_gr%sym_cheshire)     ) CALL deallocate_array ( sp_gr%sym_cheshire,     'sym_cheshire'     )

    END SUBROUTINE deallocate_space_group

    SUBROUTINE get_polar_axes_from_sp_group ( sp_group, lpaxis )
!
!       Purpose:
!       =======
!       Recalculates polar axes info from sp_group data type
!
        TYPE(space_group),     INTENT(IN)  :: sp_group
        LOGICAL, DIMENSION(:), INTENT(OUT) :: lpaxis
!       Local variables:
        TYPE(vector)                       :: xyz
        TYPE(vector)                       :: xyzm
        REAL(KIND=wp), DIMENSION(3)        :: xyz_test
!       Counters:
        INTEGER                            :: i
        INTEGER                            :: m

!       Set all axes as polar:
        lpaxis = .TRUE.

!       Not a P1 space group:
        IF ( sp_group%space_group_number /= 1 ) THEN
            xyz = (/ 0.13_wp, 0.17_wp, 0.19_wp /)
            DO m = 2, sp_group%number_of_symops
                xyzm     = ( .SYMA. sp_group%SYM(m) ) * xyz
                xyz_test = xyzm - xyz
                DO i = 1, 3
                    IF ( ABS ( xyz_test(i) ) > 0.001_wp ) lpaxis(i) = .FALSE.
                ENDDO
            ENDDO
        ENDIF
    END SUBROUTINE get_polar_axes_from_sp_group

    SUBROUTINE put_in_ASU ( xyz, sp_group, iprint )   
!
!       Purpose:
!       =======
!       Moves xyz vector into quasi-asymmetric unit trying to get Z co-ordinate as small as possible
!       
!       Date:                Programmer:            History of changes:
!       ====                 ==========             ==================
!       Nov 2004             Strokopytov B.         Original code
!       Sep 2005                -"-                 Changed mtz_1 to sp_group
!       Oct 2005                -"-                 iprint added
!       Dec 2006                -"-                 SUBTLE BUG corrected (see below)
!
        TYPE(vector),      INTENT(INOUT)        :: xyz
        TYPE(space_group), INTENT(IN)           :: sp_group
        INTEGER,           INTENT(IN), OPTIONAL :: iprint
!       Local varibables:
        INTEGER, PARAMETER                      :: trmax = 2
        REAL(KIND=wp), DIMENSION(3)             :: xyz_sym
        REAL(KIND=wp), DIMENSION(3)             :: xyz_temp
        REAL(KIND=wp), DIMENSION(3)             :: xyz_ASU
        REAL(KIND=wp), DIMENSION(3)             :: translation
!       Counters:
        INTEGER                                 :: iu
        INTEGER                                 :: iv
        INTEGER                                 :: iw
        INTEGER                                 :: m

        xyz_ASU = (/1.0_wp, 1.0_wp, 1.0_wp/)
        symmetry:DO m = 1, sp_group%number_of_symops
            xyz_sym = sp_group%SYM(m) * xyz
            along_z:DO iw = -trmax, trmax
                translation(3) = REAL(iw, KIND=wp)
                along_y:DO iv = -trmax, trmax
                    translation(2) = REAL(iv, KIND=wp)
                    along_x:DO iu = -trmax, trmax
                        translation(1) = REAL(iu, KIND=wp)
                        xyz_temp = xyz_sym + translation
!                       WRITE(*,"(' xyz_temp=',3ES14.6, ' xyz_ASU=',3ES14.6)") xyz_temp, xyz_ASU
                        IF ( xyz_temp(1) >= 0.0_wp .AND. xyz_temp(2) >= 0.0_wp .AND. xyz_temp(3) >= 0.0_wp ) THEN
!
!                           Very subtle bug detected in some rare circumstances when x-coordinate is
!                           extremely small (less than machine precision ) we may get -eps but 1 - eps = 1.0_wp. 
!                           Solution is to allow equality to unity... This is better then missing
!                           appropriate coordinate arrangement:
                            IF ( xyz_temp(1) <= 1.0_wp .AND. xyz_temp(2) <= 1.0_wp .AND. xyz_temp(3) <= 1.0_wp ) THEN
!                                WRITE(*,"(' xyz_temp=',3F7.2, ' xyz_ASU=',3F7.2)") xyz_temp, xyz_ASU

!                               We are inside unit cell 0-1,0-1,0-1 here:
                                IF ( xyz_temp(3) <= xyz_ASU(3) .AND. xyz_temp(2) <= xyz_ASU(2) &
                                                               .AND. xyz_temp(1) <  xyz_ASU(1) ) THEN
                                    xyz_ASU = xyz_temp
                                ELSE IF ( xyz_temp(3) <= xyz_ASU(3) .AND. xyz_temp(2) < xyz_ASU(2) ) THEN
                                    xyz_ASU = xyz_temp
                                ELSE IF ( xyz_temp(3) < xyz_ASU(3) ) THEN
                                    xyz_ASU = xyz_temp 
                                ENDIF

                                IF ( PRESENT ( iprint ) ) THEN
                                    WRITE(*,"(' xyz_temp=',3F7.2, ' xyz_ASU=',3F7.2)") xyz_temp, xyz_ASU
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO along_x
                ENDDO along_y
            ENDDO along_z
        ENDDO symmetry
        xyz = xyz_ASU
        IF ( PRESENT ( iprint ) ) WRITE(*,"(' resulting xyz=', 3F7.2)") xyz_ASU
    END SUBROUTINE put_in_ASU

    FUNCTION gridsize ( grid_size0, sp_group )
!
!       Purpose:
!       =======
!       Calculates a valid grid size for FFT
!
!       Date:       Programmer:      Description of changes:
!       ====        ==========       ======================
!       Oct 2005    Strokopytov B.   Lifted from `dm' by K.Cowtan
!
        TYPE(vector_int)              :: gridsize
        TYPE(vector_int),  INTENT(IN) :: grid_size0
        TYPE(space_group), INTENT(IN) :: sp_group
!       Local variables:
        INTEGER, DIMENSION(3)         ::  uvw
        INTEGER, DIMENSION(3)         ::  v
!       Counters:
        INTEGER                       ::  i
        INTEGER                       ::  k
!
!       Set up factors required by symmetry and fft:
        uvw = 2
        DO i = 1, sp_group%number_of_symops
            v = .SYMV.  sp_group%sym_int(i)
            DO k = 1, SIZE ( v )
                IF ( v(k) /= 0 ) THEN
                    uvw(k)  = MAX( uvw(k), tsym_factor/v(k) )
                ENDIF
            ENDDO
        ENDDO

!       *****************************************
!       No odd numbers
!       might be useful for huge map calculations
!       uvw = 8
!       *****************************************
        v = grid_size0
        DO k = 1, 3
            IF ( MOD (uvw(k), 2) == 1 )  uvw(k) = 2 * uvw(k)
            v(k) = fftfac( v(k), uvw(k) )
        ENDDO

!       Print the result:
        WRITE (*,110) v, uvw
 110    FORMAT( ' GRIDSIZE> Grid dimensions ', 3I4,'  MUST contain the following prime factors',/, &
                ' GRIDSIZE> for agreement with symmetry restrictions - ', 3I4)

        gridsize = v
    END FUNCTION gridsize

!   Set of allocation/deallocation routines:
    SUBROUTINE allocate_1d_symop_int_array ( array, n, array_name )
        TYPE(symop_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                    INTENT(IN)           :: n
        CHARACTER(LEN=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local variables:
        INTEGER                                                          :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_symop_int_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_symop_int_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_symop_int_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_symop_int_array')
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_symop_int_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_symop_int_array

    SUBROUTINE deallocate_1d_symop_int_array ( array, array_name )
        TYPE(symop_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(LEN=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local variables:
        INTEGER                                                          :: istat

        IF ( ALLOCATED ( array ) ) THEN
           DEALLOCATE ( array, STAT = istat )
           IF ( istat /= 0 ) THEN
                IF ( PRESENT( array_name ) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_symop_int_array')
                ELSE
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_symop_int_array')
                ENDIF
           ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_symop_int_array')
            ENDIF
        ENDIF

    END SUBROUTINE deallocate_1d_symop_int_array

    SUBROUTINE allocate_1d_symop_array ( array, n, array_name )
        TYPE(symop), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                INTENT(IN)           :: n
        CHARACTER(LEN=*),                       INTENT(IN), OPTIONAL :: array_name
!       Local variables:
        INTEGER                                                      :: istat
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_symop_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_symop_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT=istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_symop_array')
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_symop_array')
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_symop_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_symop_array

    SUBROUTINE deallocate_1d_symop_array ( array, array_name )
        TYPE(symop),      DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(LEN=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local variables:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_symop_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_symop_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_symop_array')
            ENDIF
        ENDIF

    END SUBROUTINE deallocate_1d_symop_array

END MODULE symmetry_manip
