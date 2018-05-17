MODULE basic_pdb_manip
USE aniso_manip
USE aniso_symmetry_manip
USE atom_images
USE basic_rot_operations
USE basic_symmetry_operations
!USE mr_list_manip
USE select_kinds
USE symmetry_manip
USE sparse_basic
USE sorting_facilities
USE ucards
USE util
USE vectors
IMPLICIT NONE
TYPE :: atomsf
!   First line:
    CHARACTER(LEN=4)                         :: ID
    INTEGER                                  :: ngauss != 5 ! if nothing than 5
!   Second line:
    INTEGER                                  :: atomic_charge
    INTEGER                                  :: number_of_electrons
!   REAL(KIND=wp)                            :: c
!   Third line:
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: a
!   Forth line:
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: b
!   Fifth line:
    REAL(KIND=wp), DIMENSION(2)              :: CU
    REAL(KIND=wp), DIMENSION(2)              :: MO
END TYPE atomsf

TYPE :: sf_database
    TYPE(atomsf), DIMENSION(:), ALLOCATABLE  :: record
END TYPE sf_database

TYPE :: pdb
!   vector of rigid body parameteres
    REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: p
!
!   connectivity? from pdb? SSBOND? LINK?
    TYPE(atom_image),  DIMENSION(:),   ALLOCATABLE :: all_atom_images
    TYPE(coo_matrix)                               :: COO
    TYPE(csr_matrix)                               :: CSR
!   Normally array of pointers to diagonal blocks:
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: pblock   ! size = np
    TYPE(sf_database)                              :: atomsf_table 
    REAL(KIND=sp),     DIMENSION(:,:), ALLOCATABLE :: su
    REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: w        ! overall weight on various geometric/X-ray terms
    REAL(KIND=wp)                                  :: trace
    REAL(KIND=wp)                                  :: trace_xyz
    REAL(KIND=wp)                                  :: avg_delta_r
    REAL(KIND=wp)                                  :: sum_sigsq_xyz
    REAL(KIND=wp)                                  :: sum_sigsq_all
    REAL(KIND=wp)                                  :: f0
!   Old PDB params:
    CHARACTER(len=80), DIMENSION(:),   ALLOCATABLE :: remark                          
    CHARACTER(len=80), DIMENSION(:),   ALLOCATABLE :: scale_cards                     
    CHARACTER(len=6),  DIMENSION(:),   ALLOCATABLE :: label
!   ANISOU label:                          
    CHARACTER(LEN=6),  DIMENSION(:),   ALLOCATABLE :: ulabel
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: atom_number                     
    CHARACTER(LEN=4),  DIMENSION(:),   ALLOCATABLE :: atom_name                       
    CHARACTER(LEN=1),  DIMENSION(:),   ALLOCATABLE :: altloc
    CHARACTER(len=3),  DIMENSION(:),   ALLOCATABLE :: residue_name                    
    CHARACTER(len=1),  DIMENSION(:),   ALLOCATABLE :: chain_name                      
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: residue_number                  
    CHARACTER(LEN=1),  DIMENSION(:),   ALLOCATABLE :: insertion_code
    TYPE (vector),     DIMENSION(:),   ALLOCATABLE :: xyz                             
    REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: occ                             
    REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: biso
    CHARACTER(LEN=4),  DIMENSION(:),   ALLOCATABLE :: segid
    CHARACTER(LEN=2),  DIMENSION(:),   ALLOCATABLE :: elsymb      ! element symbol
    CHARACTER(LEN=2),  DIMENSION(:),   ALLOCATABLE :: charge      ! atom charge 
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: indxyz      ! pointer to atom block            
!   Scattering factors from Int. Tables:
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: atom_type   ! atom_sf_table entry number
    LOGICAL(KIND=1),   DIMENSION(:,:), ALLOCATABLE :: refinable
!   Anisotropic factors for individual atoms:
    TYPE(aniso),       DIMENSION(:),   ALLOCATABLE :: u
!   Additional arrays:
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: residue_start                   
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: chain_start                     
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: number_of_atoms_in_this_residue 
!   Unit cell parameters:
    REAL(KIND=wp),     DIMENSION(6)                :: cell
!   Conversion matrices:
    TYPE(matrix)                                   :: ORT
    TYPE(matrix)                                   :: DEORT
!   Alternative location identifier introduced Feb 2009 BVS:
    CHARACTER(len=80)                              :: fmt='(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,A4,A2,A2)'
    CHARACTER(len=80)                              :: fmt_aniso='(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,1X,6I7,2X,A4,A2,A2)'
    LOGICAL                                        :: aniso_in_pdb
END TYPE

PRIVATE
PUBLIC :: ALLOCATED
PUBLIC :: allocate_pdb
PUBLIC :: allocate_pdb_array
PUBLIC :: apply_shifts
PUBLIC :: atomsf
PUBLIC :: ASSIGNMENT ( = )
PUBLIC :: centroid
PUBLIC :: count_atom_types
PUBLIC :: deallocate_pdb
PUBLIC :: deallocate_pdb_array
PUBLIC :: ellipsoid_axes
PUBLIC :: get_Bmin
PUBLIC :: get_xold
PUBLIC :: max_dist_in_pdb
PUBLIC :: matrix_pointers
PUBLIC :: min_ellipsoid_half_axis
PUBLIC :: OPERATOR ( - )
PUBLIC :: OPERATOR ( * )
PUBLIC :: OPERATOR ( + )
PUBLIC :: OPERATOR ( .RMSD. )
PUBLIC :: pdb
PUBLIC :: pdb_radius
PUBLIC :: prepare_vector_for_pairs
PUBLIC :: prepare_vector_of_same_type_atoms
PUBLIC :: read_pdb
PUBLIC :: read_pdb_array
PUBLIC :: refinable_occupancies
PUBLIC :: rms
PUBLIC :: set_xold !restore parameters
PUBLIC :: sf_database
PUBLIC :: SIZE
PUBLIC :: subtract_centroid_from_pdb_array
PUBLIC :: write_pdb
INTERFACE ALLOCATED
    MODULE PROCEDURE allocated_atomsf
    MODULE PROCEDURE allocated_sf_database
    MODULE PROCEDURE allocated_pdb
END INTERFACE

INTERFACE SIZE
    MODULE PROCEDURE size_of_atomsf
    MODULE PROCEDURE size_of_sf_database
    MODULE PROCEDURE pdb_size
END INTERFACE

! Some useful interfaces:
INTERFACE allocate_pdb
    MODULE PROCEDURE allocate_pdb_using_template
    MODULE PROCEDURE allocate_pdb_using_sizes
END INTERFACE allocate_pdb

! Declare interface operators:
INTERFACE ASSIGNMENT ( = )
    MODULE PROCEDURE copy_atomsf
    MODULE PROCEDURE copy_sf_database
    MODULE PROCEDURE copy_pdb 
    MODULE PROCEDURE quadratic_tensor
END INTERFACE

INTERFACE OPERATOR ( + )
    MODULE PROCEDURE add_vector_to_pdb
END INTERFACE

INTERFACE OPERATOR ( - )
    MODULE PROCEDURE subtract_vector_from_pdb
END INTERFACE

INTERFACE OPERATOR ( * )
    MODULE PROCEDURE matrix_times_pdb
    MODULE PROCEDURE pdb_times_matrix
    MODULE PROCEDURE pdb_times_symop
    MODULE PROCEDURE symop_times_pdb
END INTERFACE

INTERFACE OPERATOR( .RMSD. )
    MODULE PROCEDURE rms
END INTERFACE

CONTAINS 
    FUNCTION allocated_pdb ( pdb_1 )
        LOGICAL               :: allocated_pdb
        TYPE(pdb), INTENT(IN) :: pdb_1

        allocated_pdb = ALLOCATED ( pdb_1%xyz )
    END FUNCTION allocated_pdb
 
    FUNCTION pdb_size ( pdb_1 )
        INTEGER               :: pdb_size
        TYPE(pdb), INTENT(IN) :: pdb_1

        pdb_size = SIZE ( pdb_1%xyz )
    END FUNCTION pdb_size

    SUBROUTINE allocate_pdb_array ( pdb_array, n, cards_1 )
        TYPE(pdb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: pdb_array
        INTEGER,                              INTENT(IN), OPTIONAL :: n
        TYPE(cards),                          INTENT(IN), OPTIONAL :: cards_1
!       Local variables:
        INTEGER                                             :: istat

        IF ( ALLOCATED ( pdb_array ) ) THEN
            CALL deallocate_pdb_array ( pdb_array )
        ENDIF

!       Note how pdb_array is allocated:
        IF ( PRESENT ( cards_1 ) ) THEN
            ALLOCATE ( pdb_array(cards_1%start_model:cards_1%number_of_models), STAT=istat )
        ELSE IF ( PRESENT ( n ) ) THEN
            ALLOCATE ( pdb_array(1:n), STAT=istat )
        ELSE
            CALL die ( 'Programming error. Either N or CARDS_1 MUST BE PRESENT', 'allocate_pdb_array' )
        ENDIF
        IF ( istat /= 0 ) THEN
            CALL die ( 'Failed to allocate pdb array.', 'allocate_pdb_array' )
        ENDIF

    END SUBROUTINE allocate_pdb_array

    SUBROUTINE deallocate_pdb_array ( pdb_array )
        TYPE(pdb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: pdb_array
!       Local variables:
        INTEGER                                             :: istat
        INTEGER,   DIMENSION(1)                             :: lb
        INTEGER,   DIMENSION(1)                             :: ub
!       Counters:
        INTEGER                                             :: i

        
        IF ( ALLOCATED ( pdb_array ) ) THEN
            lb = LBOUND ( pdb_array )
            ub = UBOUND ( pdb_array )
            DO i = lb(1), ub(1)
                IF ( ALLOCATED ( pdb_array(i) ) ) CALL deallocate_pdb ( pdb_array(i) )
            ENDDO

            DEALLOCATE ( pdb_array, STAT=istat )
            IF ( istat /= 0 ) THEN
                CALL messag( 'Failed to deallocate pdb array.', 'deallocate_pdb_array' )
            ENDIF
        ELSE
            CALL warn( 'pdb array has been deallocated already.', 'deallocate_pdb_array' )
        ENDIF

    END SUBROUTINE deallocate_pdb_array

    SUBROUTINE copy_pdb ( pdb_1, pdb_2 )
        TYPE(pdb), INTENT(INOUT) :: pdb_1
        TYPE(pdb), INTENT(IN)    :: pdb_2
!       Local variables:
        INTEGER                  :: n

!       Deal with remarks and their memory:
        IF ( ALLOCATED ( pdb_2%remark) ) THEN
            CALL allocate_array(pdb_1%remark, SIZE ( pdb_2%remark ) )
            pdb_1%remark       = pdb_2%remark    
        ENDIF

!       Deal with scale_cards:
        IF ( ALLOCATED (pdb_2%scale_cards ) ) THEN
            CALL allocate_array(pdb_1%scale_cards, SIZE ( pdb_2%scale_cards ) )
            pdb_1%scale_cards  = pdb_2%scale_cards
        ENDIF

        IF ( ALLOCATED ( pdb_2%xyz ) ) THEN
            n = SIZE(pdb_2%xyz)
        ELSE
            CALL die('Programming error. pdb_2%xyz has not been allocated.', 'copy_pdb')
        ENDIF

!       PBLOCK added. MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_2%pblock ) ) THEN
            CALL allocate_array ( pdb_1%pblock,  SIZE ( pdb_2%pblock ) )
            pdb_1%pblock = pdb_2%pblock
        ENDIF
    
!       Start copying xyz data:
        IF ( ALLOCATED ( pdb_2%label ) ) THEN
            CALL allocate_array(pdb_1%label, n)
            pdb_1%label = pdb_2%label
        ENDIF

        IF ( ALLOCATED ( pdb_2%ulabel ) ) THEN
            CALL allocate_array(pdb_1%ulabel, n)
            pdb_1%ulabel = pdb_2%ulabel
        ENDIF

        IF ( ALLOCATED ( pdb_2%atom_number ) ) THEN
            CALL allocate_array(pdb_1%atom_number, n)
            pdb_1%atom_number = pdb_2%atom_number
        ENDIF

        IF ( ALLOCATED ( pdb_2%atom_name ) ) THEN
            CALL allocate_array(pdb_1%atom_name, n)
            pdb_1%atom_name = pdb_2%atom_name
        ENDIF

!       MAR 2009 ALTLOC ADDED:
        IF ( ALLOCATED ( pdb_2%altloc ) ) THEN
            CALL allocate_array ( pdb_1%altloc, n )
            pdb_1%altloc = pdb_2%altloc
        ENDIF

        IF ( ALLOCATED(pdb_2%residue_name) ) THEN
            CALL allocate_array(pdb_1%residue_name, n)
            pdb_1%residue_name   = pdb_2%residue_name
        ENDIF

        IF ( ALLOCATED(pdb_2%chain_name) ) THEN
            CALL allocate_array(pdb_1%chain_name, n)
            pdb_1%chain_name     = pdb_2%chain_name
        ENDIF

        IF ( ALLOCATED ( pdb_2%residue_number ) ) THEN
            CALL allocate_array(pdb_1%residue_number, n)
            pdb_1%residue_number = pdb_2%residue_number
        ENDIF

!       Rarely needed but necessary:
        IF ( ALLOCATED ( pdb_2%insertion_code ) ) THEN
            CALL allocate_array ( pdb_1%insertion_code, n )
            pdb_1%insertion_code = pdb_2%insertion_code
        ENDIF

        IF ( ALLOCATED(pdb_2%xyz) ) THEN
            CALL allocate_array(pdb_1%xyz, n)
            pdb_1%xyz = pdb_2%xyz
        ENDIF

        IF ( ALLOCATED(pdb_2%occ) ) THEN
            CALL allocate_array(pdb_1%occ, n)
            pdb_1%occ = pdb_2%occ
        ENDIF

        IF ( ALLOCATED ( pdb_2%biso ) ) THEN
            CALL allocate_array(pdb_1%biso, n)
            pdb_1%biso = pdb_2%biso
        ENDIF

!       JUL 2007 SEGID ADDED:
        IF ( ALLOCATED ( pdb_2%segid ) ) THEN
            CALL allocate_array ( pdb_1%segid, n )
            pdb_1%segid = pdb_2%segid
        ENDIF

!       MAR 2009 ELSYMB & CHARGE ADDED:
        IF ( ALLOCATED ( pdb_2%elsymb ) ) THEN
            CALL allocate_array ( pdb_1%elsymb, n )
            pdb_1%elsymb = pdb_2%elsymb
        ENDIF

        IF ( ALLOCATED ( pdb_2%charge ) ) THEN
            CALL allocate_array ( pdb_1%charge, n )
            pdb_1%charge = pdb_2%charge
        ENDIF

!       MAR 2009 INDXYZ ADDED:
        IF ( ALLOCATED ( pdb_2%indxyz ) ) THEN
            CALL allocate_array ( pdb_1%indxyz, n )
            pdb_1%indxyz = pdb_2%indxyz
        ENDIF

        IF ( ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL allocate_array ( pdb_1%refinable, &
                                  SIZE ( pdb_2%refinable, DIM=1 ), SIZE ( pdb_2%refinable, DIM=2 ) )
            pdb_1%refinable = pdb_2%refinable
        ENDIF

!       Copy table of scattering factors and atom_type:
        IF ( ALLOCATED ( pdb_2%atom_type ) ) THEN

!           Copy sf_database:
            pdb_1%atomsf_table = pdb_2%atomsf_table

            CALL allocate_array ( pdb_1%atom_type, n )
            pdb_1%atom_type = pdb_2%atom_type
        ENDIF

        IF ( ALLOCATED ( pdb_2%atomsf_table ) ) THEN
            pdb_1%atomsf_table = pdb_2%atomsf_table
        ENDIF

!       Copy individual anisothermal factors:
        IF ( ALLOCATED ( pdb_2%u ) ) THEN
            CALL allocate_array( pdb_1%u, n )
            pdb_1%u = pdb_2%u
        ENDIF

!       Copy other data:
        IF ( ALLOCATED ( pdb_2%residue_start ) ) THEN
            CALL allocate_array(pdb_1%residue_start, SIZE ( pdb_2%residue_start ))
            pdb_1%residue_start = pdb_2%residue_start
        ENDIF

        IF ( ALLOCATED ( pdb_2%chain_start ) ) THEN
            CALL allocate_array(pdb_1%chain_start, SIZE ( pdb_2%chain_start ) )
            pdb_1%chain_start   = pdb_2%chain_start
        ENDIF

        IF ( ALLOCATED ( pdb_2%number_of_atoms_in_this_residue ) ) THEN
            CALL allocate_array(pdb_1%number_of_atoms_in_this_residue, &
                                SIZE ( pdb_2%number_of_atoms_in_this_residue ) )
            pdb_1%number_of_atoms_in_this_residue = &
            pdb_2%number_of_atoms_in_this_residue
        ENDIF

        IF ( ALLOCATED ( pdb_2%p ) ) THEN
            CALL allocate_array ( pdb_1%p, SIZE ( pdb_2%p ) )
            pdb_1%p = pdb_2%p
        ENDIF

!  FIXME
        IF ( ALLOCATED ( pdb_2%COO ) ) THEN
        ENDIF    
        
!       Stereochemical data: FIXME
!       GROUPS
!       LINK, SSBOND etc

!       Copy non-pointer data:
        pdb_1%cell  = pdb_2%cell
        pdb_1%ORT   = pdb_2%ORT
        pdb_1%DEORT = pdb_2%DEORT
        pdb_1%fmt   = pdb_2%fmt
        pdb_1%fmt_aniso    = pdb_2%fmt_aniso
        pdb_1%aniso_in_pdb = pdb_2%aniso_in_pdb

!       Hopefully we are done here...

    END SUBROUTINE copy_pdb
   
    FUNCTION rms ( pdb_1, pdb_2 )
        REAL(KIND=wp)         :: rms
        TYPE(pdb), INTENT(IN) :: pdb_1
        TYPE(pdb), INTENT(IN) :: pdb_2
!       Local variables:
        INTEGER               :: n
        REAL(KIND=wp)         :: rmsd
        TYPE(vector)          :: d
!       Counters:
        INTEGER               :: i

!       Checkz:
        IF ( .NOT. ALLOCATED ( pdb_1%xyz ) ) THEN
            CALL die ('pdb_1 has not been properly initialized.','rms')
        ELSE IF ( .NOT. ALLOCATED ( pdb_2%xyz ) ) THEN
            CALL die('pdb_2 has not been properly initialized.','rms')
        ELSE
            IF ( SIZE ( pdb_1 ) /= SIZE ( pdb_2 ) ) THEN
                WRITE(*,"(' RMS> ',' n1 =',I6, ' n2 =',I6)")&
                SIZE(pdb_1%xyz), SIZE(pdb_2%xyz)
                CALL die(' pdb_1 and pdb_2 are not of equal size','rms')
            ENDIF
        ENDIF

        n = SIZE ( pdb_1 )
        rmsd = 0.0_wp
        DO i = 1, n
            d    = pdb_1%xyz(i) - pdb_2%xyz(i)
            rmsd = rmsd + ( d .DOT. d )
        ENDDO
        rms = SQRT ( rmsd / n )
    END FUNCTION rms

    FUNCTION max_dist_in_pdb ( pdb_1 )
        REAL(KIND=wp)         :: max_dist_in_pdb
        TYPE(pdb), INTENT(IN) :: pdb_1
!       Local variables:
        INTEGER               :: n
        REAL(KIND=wp)         :: r2
        REAL(KIND=wp)         :: r2_max
        TYPE(vector)          :: d
!       Counters:
        INTEGER               :: i
        INTEGER               :: j

        IF ( ALLOCATED ( pdb_1%xyz ) ) THEN
            r2_max = 0.0_wp
            n = SIZE ( pdb_1%xyz )
            DO i = 1, n - 1
                DO j = i + 1, n
                    d  = pdb_1%xyz(j) - pdb_1%xyz(i)
                    r2 = d .DOT. d
                    IF ( r2 > r2_max ) THEN
                        r2_max = r2
                    ENDIF
                ENDDO
            ENDDO
            max_dist_in_pdb = SQRT ( r2_max )
        ELSE
            CALL die ('pdb_1 has not been properly initialized', 'max_dist_in_pdb')
        ENDIF
    END FUNCTION max_dist_in_pdb

    FUNCTION pdb_radius ( pdb_1 )
        REAL(KIND=wp)         :: pdb_radius
        TYPE(pdb), INTENT(IN) :: pdb_1
!       Local variables:
        TYPE(pdb)             :: pdb_temp
        INTEGER               :: n
        REAL(KIND=wp)         :: r2
        REAL(KIND=wp)         :: r2_max
        TYPE(vector)          :: d
!       Counters:
        INTEGER               :: i

        n = SIZE ( pdb_1%xyz )
        pdb_temp = pdb_1 - centroid(pdb_1)
        r2_max   = 0.0_wp
        DO i = 1, n
            d  = pdb_temp%xyz(i)
            r2 = d .DOT. d
            IF ( r2 > r2_max ) THEN
                r2_max = r2
            ENDIF
        ENDDO
        pdb_radius = SQRT ( r2_max )
        CALL deallocate_pdb(pdb_temp)
    END FUNCTION pdb_radius

    FUNCTION symop_times_pdb ( sym_1, pdb_2 ) RESULT ( pdb_temp )
!
!       Purpose:
!       =======
!       Applies symmetry operator to pdb coordinates
!       Note that sym_1 has to be orthogonalized first (unless it was orthogonal already)
!       using function convert_symop_to_ort_system from symmetry_manip module
!
        TYPE(pdb)               :: pdb_temp
        TYPE(symop), INTENT(IN) :: sym_1 
        TYPE(pdb),   INTENT(IN) :: pdb_2
!       Counters:
        INTEGER                 :: i

        IF ( ALLOCATED ( pdb_2 ) ) THEN
            pdb_temp = pdb_2
            DO i = 1, SIZE ( pdb_2%xyz )
                pdb_temp%xyz(i) = sym_1 * pdb_2%xyz(i)
            ENDDO
        ELSE
            CALL die( 'pdb_2 has not been properly initialised', 'symop_times_pdb' )
        ENDIF

    END FUNCTION symop_times_pdb

    FUNCTION pdb_times_symop ( pdb_1, sym_2 ) RESULT ( pdb_temp )
!
!       Purpose:
!       =======
!       Applies symmetry operator to pdb co-ordinates
!       Note that sym_1 has to be orthogonalized first (unless it was orthogonal already)
!       using function convert_symop_to_ort_system from symmetry_manip module
!
        TYPE(pdb)               :: pdb_temp
        TYPE(pdb),   INTENT(IN) :: pdb_1
        TYPE(symop), INTENT(IN) :: sym_2
!       Counters:
        INTEGER                 :: i

        IF ( ALLOCATED ( pdb_1 ) ) THEN
            pdb_temp = pdb_1
            DO i = 1, SIZE ( pdb_1 )
                pdb_temp%xyz(i) = sym_2 * pdb_1%xyz(i)
            ENDDO
        ELSE
            CALL die ( 'pdb_1 has not been properly initialised', 'pdb_times_symop' )
        ENDIF

    END FUNCTION pdb_times_symop

    FUNCTION matrix_times_pdb ( ROT, pdb_1 ) RESULT ( pdb_temp )
!
!       Purpose:
!       =======
!       Rotates pdb via given matrix ROT.
!
!       Date:          Programmer:           History of changes:
!       ====           ==========            ==================
!       Oct 2003       B.Strokopytov         Original code.
!       Oct 2007       B.Strokopytov         Cosmetic changes. Comments.
!       
!       Note:
!       ====
!       Any type of matrix can be used. Including deorthogonalisation matrices.
!       THEREFORE NO MATRIX DETERMINANT checking.
!
        TYPE(pdb)                :: pdb_temp
        TYPE(matrix), INTENT(IN) :: ROT 
        TYPE(pdb),    INTENT(IN) :: pdb_1
!       Counters:
        INTEGER                  :: i

        IF ( ALLOCATED ( pdb_1 ) ) THEN

            pdb_temp = pdb_1 

            DO i = 1, SIZE ( pdb_1 ) 
                pdb_temp%xyz(i) = ROT * pdb_1%xyz(i)
            ENDDO

        ELSE

            CALL die ( 'pdb_1 structure has not been initialized', 'matrix_times_pdb' )

        ENDIF

    END FUNCTION matrix_times_pdb

    FUNCTION pdb_times_matrix ( pdb_1, ROT ) RESULT ( pdb_temp )
        TYPE(pdb)                :: pdb_temp
        TYPE(pdb),    INTENT(IN) :: pdb_1
        TYPE(matrix), INTENT(IN) :: ROT
!       Counters:
        INTEGER                  :: i

        IF ( ALLOCATED ( pdb_1 ) ) THEN

            pdb_temp = pdb_1

            DO i = 1, SIZE ( pdb_1 )
                pdb_temp%xyz(i) = ROT * pdb_1%xyz(i)
            ENDDO
        ELSE
            CALL die ( 'pdb_1 structure has not been initialized', 'pdb_times_pdb' )
        ENDIF

    END FUNCTION pdb_times_matrix

    FUNCTION centroid ( pdb_1 )
        TYPE(vector)                :: centroid
        TYPE(pdb), INTENT(IN)       :: pdb_1
!       Local variables:
        INTEGER                     :: n
        TYPE(vector)                :: geom_centre_of_mass
!       For printing:
        REAL(KIND=wp), DIMENSION(3) :: xyz
!       Counters:
        INTEGER                     :: i
!
        IF ( ALLOCATED ( pdb_1%xyz ) ) THEN
            n = SIZE ( pdb_1 )
            IF ( n <= 0 ) CALL die('Number of atoms on input equals zero.', 'centroid')
            IF ( n == 1 ) CALL warn('Just one atom on input?', 'centroid')

!           Initialize:
            geom_centre_of_mass = 0.0_wp
            DO i = 1, n
                geom_centre_of_mass = geom_centre_of_mass + pdb_1%xyz(i)
            ENDDO
            centroid = geom_centre_of_mass / REAL ( n, KIND=wp )

!           Global debugging controls:
            IF ( debug > 5 ) THEN
                xyz = centroid
                WRITE(*,"(' CENTROID> ',' centroid=', 3F8.2)") xyz
            ENDIF

            RETURN
        ELSE
            CALL die ('pdb file has not been initialized.','centroid')
        ENDIF     
    END FUNCTION centroid

    SUBROUTINE quadratic_tensor ( mat_1, pdb_2 )
        TYPE(matrix), INTENT(OUT)     :: mat_1
        TYPE(pdb), INTENT(IN)         :: pdb_2
!       Local variables:
        INTEGER                       :: na
        REAL(KIND=wp), DIMENSION(3,3) :: A
        REAL(KIND=wp), DIMENSION(3)   :: v
!       Counters:
        INTEGER                       :: iat
        INTEGER                       :: i
        INTEGER                       :: j
!
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            na = SIZE ( pdb_2 )
            IF ( na <= 0 ) CALL die ('Number of atoms on input equals zero.','quadratic_tenzor')

!           Set 3x3 matrix to zero:
            A = 0.0_wp
            DO iat = 1, na
                v = pdb_2%xyz(iat)
                DO i = 1, 3
                    DO j = 1, 3
                        A(i,j) = A(i,j) + v(i) * v(j)
                    ENDDO
                ENDDO
            ENDDO

!           Normalize matrix:
            A = A / REAL (na, KIND=wp )

!           get position of centre of mass:
            v = centroid  ( pdb_2 )

!           extract centroid tensor:
            DO i = 1, 3
                DO j = 1, 3
                    A(i,j) = A(i,j) - v(i) * v(j)
                ENDDO
            ENDDO
        ELSE
            CALL die ('Programming error: pdb_2 has not been allocated.','quadratic_tensor')
        ENDIF

        mat_1 = A

    END SUBROUTINE quadratic_tensor
 
    FUNCTION squared_ellipsoid_half_axes ( pdb_1 )
        TYPE (vector)          :: squared_ellipsoid_half_axes
        TYPE (pdb), INTENT(IN) :: pdb_1
!       Local variables:
        TYPE (matrix)          :: mat_1

        mat_1 = pdb_1
        squared_ellipsoid_half_axes = eigen_values ( mat_1 )
    END FUNCTION squared_ellipsoid_half_axes 
!
    FUNCTION ellipsoid_axes ( pdb_1 )
        TYPE(vector)          :: ellipsoid_axes
        TYPE(pdb), INTENT(IN) :: pdb_1
!       Local variables:
        TYPE(vector)          :: ellipsoid_half_axes

        ellipsoid_half_axes   = squared_ellipsoid_half_axes ( pdb_1 )
        ellipsoid_axes        = 2.0_wp * SQRT ( ellipsoid_half_axes )
    END FUNCTION ellipsoid_axes

    FUNCTION min_ellipsoid_half_axis ( pdb_1, iprint )
        REAL(KIND=wp)                       :: min_ellipsoid_half_axis
        TYPE(pdb),     INTENT(IN)           :: pdb_1
        INTEGER,       INTENT(IN), OPTIONAL :: iprint
!       Local variables:
        REAL(KIND=wp), DIMENSION(3) :: ellipsoid_half_axes

        ellipsoid_half_axes = squared_ellipsoid_half_axes ( pdb_1 )
        ellipsoid_half_axes = SQRT ( ellipsoid_half_axes )

        IF ( PRESENT ( iprint ) ) THEN
            IF ( iprint > 0 ) THEN
                WRITE ( *, "(' MIN_ELLIPSOID_HALF_AXIS> ', 'Ellipsoid half axes     =', 3F8.2)" ) &
                ellipsoid_half_axes
                WRITE ( *, "(' MIN_ELLIPSOID_HALF_AXIS> ', 'Min ellipsoid half axis =', 3F8.2)" ) &
                MINVAL( ellipsoid_half_axes )
            ENDIF
        ENDIF

        min_ellipsoid_half_axis = MINVAL ( ellipsoid_half_axes )
    END FUNCTION min_ellipsoid_half_axis

    FUNCTION add_vector_to_pdb ( pdb_1, vec_1 ) RESULT ( pdb_temp )
        TYPE(pdb)                :: pdb_temp 
        TYPE(pdb),    INTENT(IN) :: pdb_1
        TYPE(vector), INTENT(IN) :: vec_1
!       Counters:
        INTEGER                  :: i

        IF ( ALLOCATED ( pdb_1%xyz ) ) THEN
            pdb_temp = pdb_1
            DO i = 1, SIZE ( pdb_1 )
                pdb_temp%xyz(i) = pdb_1%xyz(i) + vec_1
            ENDDO
        ELSE
            CALL die('pdb_1 has not been properly initialised','add_vector_to_pdb')
        ENDIF
    
    END FUNCTION add_vector_to_pdb

    FUNCTION subtract_vector_from_pdb ( pdb_1, vec_1 ) RESULT ( pdb_temp )
        TYPE (pdb)                :: pdb_temp
        TYPE (pdb), INTENT(IN)    :: pdb_1        
        TYPE (vector), INTENT(IN) :: vec_1
!       Counters:
        INTEGER                   :: i

        IF ( ALLOCATED ( pdb_1%xyz ) ) THEN
            pdb_temp = pdb_1
            DO i = 1, SIZE ( pdb_1 )
                pdb_temp%xyz(i) = pdb_1%xyz(i) - vec_1
            ENDDO       
        ELSE
            CALL die ( 'pdb_1 has not been properly initialised', 'subtract_vector_to_pdb' )
        ENDIF

    END FUNCTION subtract_vector_from_pdb

    SUBROUTINE read_pdb ( pdb_1, file_2, cards_1 )
!
!       Purpose:
!       =======
!       Estimates number of ATOM/HETATM records in FILE_2
!       Allocates arrays of appropriate length for PDB_1 object
!       Estimates number of chain, segids, residues etc.
!      
!       Date:           Programmer:            History of changes:
!       ====            ==========             ==================
!       Nov 2003        B.Strokopytov          Original code
!       Sep 2005        B.Strokopytov          ALLOCATABLE arrays intoduced
!       Jul 2007        B.Strokopytov          Minor bugs corrected. Start a new residue
!                                              as soon as new chain detected.
!       Jul 2007        B.Strokopytov          ACTION='READ' added
!       Mar 2008        B.Strokopytov          Skipping of hydrogens added(temporarily).
!
        TYPE(pdb),        INTENT(OUT)  :: pdb_1
        CHARACTER(len=*), INTENT(IN)   :: file_2
        TYPE(cards)                    :: cards_1
!       Local variables:
        CHARACTER(LEN=flen)            :: CCP4_ATOMSF_location
        CHARACTER(LEN=80)              :: line
        CHARACTER(LEN=120)             :: fmt
        CHARACTER(LEN=120)             :: fmt_aniso
        REAL(KIND=wp), DIMENSION(3)    :: xyz
        INTEGER                        :: incoor
        INTEGER                        :: istat
        LOGICAL                        :: exists
        CHARACTER(LEN=1)               :: altloc
        INTEGER, DIMENSION(6)          :: U
        TYPE(matrix)                   :: UMAT
        REAL(KIND=wp), DIMENSION(3)    :: ueigen
        LOGICAL                        :: refine_Biso
        LOGICAL                        :: refine_Us
        LOGICAL                        :: mixed
        LOGICAL                        :: no_BU_ref
        LOGICAL                        :: need_aniso
!       Counters:
        INTEGER                        :: remark_cards_number
        INTEGER                        :: scale_cards_number
        INTEGER                        :: atom_records_number
        INTEGER                        :: iat
        INTEGER                        :: natom
        INTEGER                        :: i
        INTEGER                        :: irem
        INTEGER                        :: iscal
        INTEGER                        :: main_chain_counter
        INTEGER                        :: nmod
        INTEGER                        :: naniso
        INTEGER                        :: nmessag
        CHARACTER(LEN=32),        SAVE :: srname='read_pdb'

!       Find next available unit and open the file:
        incoor = get_next_io_unit()
        OPEN ( UNIT=incoor, FILE=TRIM(file_2), ACCESS='SEQUENTIAL', FORM='FORMATTED', &
               ACTION='READ', STATUS='OLD', IOSTAT=istat )

        IF ( istat /= 0 ) THEN
            CALL die(' Failed to open .pdb file: '//file_2, 'read_pdb')
        ENDIF

!       Don't try to make it better:
        fmt = '('' READ_PDB> '',' // TRIM ( pdb_1%fmt ) // ')'
        fmt_aniso='('' READ_PDB> '',' // TRIM ( pdb_1%fmt_aniso ) // ')'

!       Check whether we want aniso refinement or/and read aniso records:
        refine_Us   = INDEX ( cards_1%refinement_mode, 'U' ) > 0
        refine_Biso = INDEX ( cards_1%refinement_mode, 'B' ) > 0
        no_BU_ref = .NOT. refine_Us .AND. .NOT. refine_Biso 
        mixed = refine_Biso .AND. refine_Us

        IF ( mixed ) THEN
            CALL messag('MIxed thermal factor refinement will be carried out.', 'read_pdb')
        ENDIF
      
!       Initialize:
        remark_cards_number  = 0
        scale_cards_number   = 0
        atom_records_number  = 0
        main_chain_counter   = 0
        naniso               = 0

!       Make first "counting" run:
        DO WHILE (.TRUE.)
            READ ( INCOOR, FMT='(A)', IOSTAT=istat ) line
            IF ( istat /= 0 ) THEN
                EXIT
            ENDIF

            IF ( line(1:6) == 'REMARK' ) remark_cards_number  = remark_cards_number  + 1

!           Do not skip hydrogens:
            IF ( line(1:4) == 'ATOM' .OR. line(1:6) == 'HETATM' ) THEN
                atom_records_number  = atom_records_number  + 1
            ENDIF

            IF ( line(1:4) == 'SCAL'   ) scale_cards_number   = scale_cards_number   + 1
            IF ( line(1:6) == 'ANISOU' ) naniso = naniso + 1
        ENDDO

!       Set to .TRUE. if anisotropic thermal factors found in the PDB:
        pdb_1%aniso_in_pdb = naniso > 0
 
        IF ( debug > 10 ) THEN
            WRITE(*,*) ' aniso_in_pdb=',  pdb_1%aniso_in_pdb
            WRITE(*,*) ' mixed=',         mixed
        ENDIF

!       Check for impossible combination of parameters:       
        IF ( .NOT. pdb_1%aniso_in_pdb .AND. mixed ) THEN
            WRITE(*,*) pdb_1%aniso_in_pdb, mixed
!           CALL deallocate_pdb ( pdb_1 )
            CALL messag('Please decide which type of refinement (B or U) your prefer.',&
                        'read_pdb')
            CALL die('User error. You can''t carry out mixed type (B and U) refinement &
                     &having no ANISOU records in the PDB file.', 'read_pdb' )
        ENDIF

!       Now we are ready to allocate arrays:
        IF ( atom_records_number > 0 ) THEN
            CALL messag('Allocating memory for pdb file...', 'read_pdb')

!           This requires a comment. We need to read aniso values in case:
!           a)  Obviously when they present in the file with condition c)
!           b)  We want to refine them even if they not in the file, again with condition c)
!           c)  No Biso refinement unless we have mixed refinement or no temperature refinement whatsoever.
!           In all this cases we need to read (and allocate) U values
!           In case aniso records present but Biso refinement wanted -- we simply
!           ignore aniso records. Obviously we can't refine both type of parameters for a single atom.
!           We could still refine aniso in this case but aniso ref. was not required by the user!!
!
!      
            need_aniso = (pdb_1%aniso_in_pdb .OR. refine_Us ) .AND.  (.NOT. refine_Biso .OR. mixed .OR. no_bu_ref )  
!            WRITE(*,*) ' refine_Us=', refine_Us
!            WRITE(*,*) ' long expression 2=', (.NOT. refine_Biso .OR. mixed .OR. no_bu_ref )
!            WRITE(*,*) .NOT. refine_Biso, mixed, no_bu_ref 
!            STOP 'logicals'
            CALL allocate_pdb ( pdb_1, natom=atom_records_number, nremarks=remark_cards_number, &
                                nscale_cards=scale_cards_number, anom=cards_1%anom, &
                                aniso=need_aniso )

!           Need this for mixed U/B(iso) refinement:
            IF ( ALLOCATED ( pdb_1%U ) ) THEN
                DO iat = 1, SIZE ( pdb_1 )
                    pdb_1%U(iat) = 0.0_wp   ! set_to_const ROUTINE
                ENDDO
            ENDIF

            CALL messag ( ' - - - Final .pdb file statistics - - - ', 'read_pdb' )
            WRITE(*, "(' READ_PDB> ', 'remark_cards_number=', A)" ) &
            TRIM ( int_to_c ( remark_cards_number ) )
            WRITE(*, "(' READ_PDB> ', 'scale_cards_number= ', A)" ) &
            TRIM ( int_to_c ( scale_cards_number ) )
            WRITE(*, "(' READ_PDB> ', 'atom_records_number=', A)" ) &
            TRIM ( int_to_c ( atom_records_number ) )
            WRITE(*, "(' READ_PBD> ', 'anisotropic atoms  =', A)" ) &
            TRIM ( int_to_c ( naniso ) )
!           Printing frequency:
            nmod = MAX ( atom_records_number / 10, 1 )
        ELSE
            WRITE(*,*)  atom_records_number
            CALL die ('No valid atomic records has been found in input file.', 'read_pdb')
        ENDIF

!       Prepare to read again:
        REWIND ( UNIT = incoor, IOSTAT = istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' incoor=', incoor
            CALL die('Failed to REWIND incoor file.', 'read_pdb')
        ENDIF

        IF ( debug > 2 ) THEN
            WRITE(*,"(' READ_PDB> ', 'Memory allocated. Reading the .pdb file...')")
        ENDIF

        nmessag = 0
        iat   = 0
        iscal = 0
        irem  = 0
        DO WHILE (.TRUE. )

            READ ( UNIT=incoor, FMT='(A)', IOSTAT=istat ) line
            IF ( istat /= 0 ) THEN
                EXIT
            ENDIF

            IF ( line(1:4) == 'ATOM' .OR. line(1:4) == 'HETA' ) THEN
                CALL ucase (line)
                iat = iat + 1
                READ(line, pdb_1%fmt) pdb_1%label(iat),          &
                                      pdb_1%atom_number(iat),    &
                                      pdb_1%atom_name(iat),      &
                                      pdb_1%altloc(iat),         &
                                      pdb_1%residue_name(iat),   & 
                                      pdb_1%chain_name(iat),     &
                                      pdb_1%residue_number(iat), &
                                      pdb_1%insertion_code(iat), &
                                      xyz,                       &
                                      pdb_1%occ(iat),            &
                                      pdb_1%biso(iat),           &
                                      pdb_1%segid(iat),          &
                                      pdb_1%elsymb(iat),         &
                                      pdb_1%charge(iat)
!
                pdb_1%xyz(iat) = xyz
                
                IF ( pdb_1%occ(iat) == 0.0_wp ) THEN 
                    IF (  pdb_1%elsymb(iat)(1:2) /= '  ' ) THEN
                        IF ( ADJUSTL ( pdb_1%elsymb(iat)(1:2) ) /= 'H ' ) THEN
                            WRITE(*,"('|',A2,'|',A2)") ADJUSTL ( pdb_1%elsymb(iat)(1:2) ), 'H '
                            WRITE(*,*) ' iat=', iat
                            WRITE(*,"(1X,A)") line
                            CALL die('Zero occupancies for non H atoms are not allowed. Please correct this...', srname)
                        ENDIF
                    ELSE
!                       Blank elsymb case, check atom name since ELSYMB is empty:
                        IF ( INDEX ( pdb_1%atom_name(iat)(1:2), 'H' ) == 0 ) THEN
                            WRITE(*,*) pdb_1%atom_name(iat)(1:2), INDEX ( pdb_1%atom_name(iat)(1:2), 'H' )
                            CALL warn('Elsymb is empty...', srname)
                            WRITE(*,*) ' iat=', iat
                            WRITE(*,"(1X,A)") line
!                           This is not a hydrogen (only hydrogens may have zero occupnacies):
                            CALL die('Zero occupancies for non H atoms are not allowed. Please correct this.', srname)
                        ENDIF              
                    ENDIF
                ENDIF

                IF ( refine_Us .AND. .NOT. mixed ) THEN

!                   No aniso in file but aniso refinement has been declared => convert Biso to aniso:
                    pdb_1%U(iat) = .CONVERT. ( pdb_1%biso(iat) / pi_squared_8 )

!                   Take care of ANISOU label (otherwise it will be blank/undefined): BUG corrected JULY 2009 BVS:
                    pdb_1%ulabel(iat) = 'ANISOU'
                ENDIF

!               Checking labels:
                IF ( pdb_1%label(iat) /= 'ATOM  ' .AND. pdb_1%label(iat) /= 'HETATM' ) THEN
                    CALL warn('Problem reading labels: '//pdb_1%label(iat), 'read_pdb')
                ENDIF

!               Here we have enough info to determine scattering factors:
                IF ( MOD ( iat, nmod ) == 1 .OR. iat == atom_records_number .OR. iat == 1  ) THEN

                    xyz = pdb_1%xyz(iat)

                    WRITE(*, fmt) pdb_1%label(iat),          &
                                  pdb_1%atom_number(iat),    &
                                  pdb_1%atom_name(iat),      &
                                  pdb_1%altloc(iat),         &
                                  pdb_1%residue_name(iat),   &
                                  pdb_1%chain_name(iat),     &
                                  pdb_1%residue_number(iat), &
                                  pdb_1%insertion_code(iat), &
                                  xyz,                       &
                                  pdb_1%occ(iat),            &
                                  pdb_1%biso(iat),           &
                                  pdb_1%segid(iat),          &
                                  pdb_1%elsymb(iat),         &
                                  pdb_1%charge(iat)

                    IF ( debug > 5 ) THEN
                        IF ( refine_Us .AND. .NOT. pdb_1%aniso_in_pdb ) THEN
                            WRITE(*,"(' READ_PDB> ', 'aniso thermal factors=', 6F10.6)") &
                            pdb_1%U(iat)%u
                        ENDIF
                    ENDIF

                ENDIF

!               Make sure atom name & residue name do not start with blank character:
                pdb_1%atom_name(iat)    = ADJUSTL ( pdb_1%atom_name(iat) )
                pdb_1%residue_name(iat) = ADJUSTL ( pdb_1%residue_name(iat) )
                pdb_1%elsymb(iat)       = ADJUSTL ( pdb_1%elsymb(iat) )

!           Take care of Anisotropic temperature factors:
            ELSE IF ( line(1:6) == 'ANISOU' .AND. need_aniso ) THEN
                CALL ucase(line)
!               FIXME
                READ(line, pdb_1%fmt_aniso) pdb_1%ulabel(iat),         &
                                            pdb_1%atom_number(iat),    &
                                            pdb_1%atom_name(iat),      &              
                                            altloc,                    &
                                            pdb_1%residue_name(iat),   &
                                            pdb_1%chain_name(iat),     &
                                            pdb_1%residue_number(iat), &
                                            pdb_1%insertion_code(iat), &
                                            U,                         &
                                            pdb_1%segid(iat),          &
                                            pdb_1%elsymb(iat),         &
                                            pdb_1%charge(iat)          


!               Copy and scale anisotropic termal factors:
                pdb_1%U(iat) = REAL (U, KIND=wp) / aniso_pdb_factor

                IF ( MOD ( iat, nmod ) == 1 .OR. iat == atom_records_number .OR. iat == 1  ) THEN

!                   Need to print something here: BVS
                    WRITE(*, fmt_aniso) pdb_1%ulabel(iat),         &
                                        pdb_1%atom_number(iat),    &
                                        pdb_1%atom_name(iat),      &
                                        altloc,                    &
                                        pdb_1%residue_name(iat),   &
                                        pdb_1%chain_name(iat),     &
                                        pdb_1%residue_number(iat), &
                                        pdb_1%insertion_code(iat), &
                                        U,                         &
                                        pdb_1%segid(iat),          &
                                        pdb_1%elsymb(iat),         &
                                        pdb_1%charge(iat)
                                   

                ENDIF

!               Make sure atom name & residue name do not start with blank character:
                pdb_1%atom_name(iat)    = ADJUSTL ( pdb_1%atom_name(iat) )
                pdb_1%residue_name(iat) = ADJUSTL ( pdb_1%residue_name(iat) )
                pdb_1%elsymb(iat)       = ADJUSTL ( pdb_1%elsymb(iat) )

!               It's not quite clear how to do this:
                IF ( pdb_1%elsymb(iat)(1:1) == ' ' ) THEN
!                   Use atom name to  recreate elsymb:
                ENDIF                       
                

!               Avoid bad aniso tensors:
                UMAT = pdb_1%U(iat)
                ueigen = eigen_values ( UMAT )

                IF ( ANY ( ueigen(1:2)/ueigen(3) < 0.1_wp ) .OR. ALL ( ueigen < 0.0_wp ) ) THEN
                    WRITE(*,*) ' iat=', iat
                    WRITE(*,*) ' eigenvalues=', ueigen
                    CALL warn('Resseting ellipsoid parameters according to MIN_ADP_RATIO...', srname)
                    pdb_1%U(iat) = pdb_1%U(iat) .CORRANI. 0.1_wp
                   
!                   Recheck:
                    UMAT = pdb_1%U(iat)
                    ueigen = eigen_values ( UMAT )
                    IF ( ANY ( ueigen < 0.0_wp ) ) THEN
                        WRITE(*,*) ' iat=', iat
                        WRITE(*,"( ' eigenvalues=', 3ES12.4)") ueigen 
                        CALL die('O-o-o-p-s-s. Negative eigenvalue(s) detected...CORRANI does not work properly.', srname)
                    ENDIF

                ENDIF

            ELSE IF ( line(1:4) == 'REMA' ) THEN

                irem  = irem + 1
                pdb_1%remark(irem) = line
                IF ( irem <= 10 ) WRITE(*,1010) pdb_1%remark(irem)

            ELSE IF ( line(1:4) == 'SCAL' ) THEN

                iscal = iscal + 1

!               No need to read more than 3 SCALE cards:
                IF ( iscal > 3 ) CYCLE
                pdb_1%scale_cards(iscal) = line
                WRITE(*,1010) pdb_1%scale_cards(iscal)

            ENDIF
!           BUG corrected: removed ADJUSTL functions from here since IAT maybe equal to zero...

        ENDDO

        natom = iat
        WRITE(*,"(' READ_PDB> ', A, ' atoms have been read in.')") & 
        TRIM ( int_to_c ( natom ) )

!        Kostya's test:
!        pdb_1%elsymb = 'C'

!       Figure out start atom number of residues, etc.:
        CALL residue_location ( pdb_1 )
        CALL ugtenv ( 'ATOMSF', CCP4_ATOMSF_location )        

!       Massive checking:
        IF ( LEN_TRIM ( CCP4_ATOMSF_location ) == 0 ) THEN
            WRITE(*,*) ' CCP4_ATOMSF_location=', CCP4_ATOMSF_location
            CALL die ( 'Failed to find ATOMSF database...', 'read_pdb' )
        ENDIF

        INQUIRE ( file =  CCP4_ATOMSF_location, exist = exists )
        IF ( .NOT. exists ) THEN
            WRITE(*,*) ' CCP4_ATOMSF_location=', CCP4_ATOMSF_location
            CALL die ('Failed to locate ATOMSF file...', 'read_pdb' )
        ENDIF

!       More general approach is to make a separate SUBROUTINE for this: FIXME?
        CALL read_sf_database ( pdb_1%atomsf_table, cards_1, CCP4_ATOMSF_location )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' READ_PDB> ', 20A4)") (pdb_1%atomsf_table%record(i)%id, i=1, SIZE ( pdb_1%atomsf_table ))
            WRITE(*,"(' READ_PDB> ', 16I6)") pdb_1%residue_start
            WRITE(*,"(' READ_PDB> ', 'Number of atoms in each residue:')")
            WRITE(*,"(' READ_PDB> ', 16I5)") pdb_1%number_of_atoms_in_this_residue
        ENDIF

        CALL assign_atom_types ( pdb_1, cards_1)

!       Deallocate scaterring factors database:
!        CALL deallocate_sf_database ( sf_data )

!       Close PDB unit:
        CLOSE ( UNIT=incoor, IOSTAT=istat )
        IF ( istat /= 0 ) CALL die ( 'Failed to close pdb file...', 'read_pdb' )
  1010  FORMAT (1X, A)
     
    END SUBROUTINE read_pdb

    SUBROUTINE read_pdb_array ( my_pdb, my_cards )
        TYPE(pdb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: my_pdb
        TYPE(cards),                          INTENT(IN)    :: my_cards
!       Local vars:
        CHARACTER(LEN=132)                                  :: file_1
!       Counters:
        INTEGER                                             :: pdb_model_number
        
        CALL allocate_pdb_array ( my_pdb, cards_1=my_cards )

!       Static pdb must go to my_pdb(0):
        CALL ugtenv ( 'STATIC', file_1 )
        IF ( TRIM ( file_1 ) /= '' ) THEN
            CALL messag ( 'Reading static pdb...', 'bruteptf' )
            CALL read_pdb ( my_pdb(0), file_1, my_cards )
!           No subtraction - static model:
        ENDIF

        DO pdb_model_number = 1, my_cards%number_of_models
            CALL ugtenv ( 'PDB'//TRIM ( int_to_c( pdb_model_number ) ), file_1 )
            CALL messag ( 'pdb file name= '//TRIM ( file_1 ), 'bruteptf' )
            CALL read_pdb ( my_pdb(pdb_model_number), file_1, my_cards )
        ENDDO

    END SUBROUTINE read_pdb_array

    SUBROUTINE subtract_centroid_from_pdb_array ( my_pdb )
!
!       Purpose:
!       =======
!       Subtracts centroids from a set of PDBs
!       except static model PDB(0)
!
!
        TYPE(pdb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: my_pdb
!       Local variables:
        INTEGER,   DIMENSION(1)                             :: lb
        INTEGER,   DIMENSION(1)                             :: ub
!       Counters:
        INTEGER                                             :: pdb_model_number

        lb = LBOUND ( my_pdb )
        ub = UBOUND ( my_pdb )
        DO pdb_model_number = lb(1), ub(1)
            my_pdb(pdb_model_number) = my_pdb(pdb_model_number) - centroid ( my_pdb(pdb_model_number) )
        ENDDO

    END SUBROUTINE subtract_centroid_from_pdb_array

    SUBROUTINE write_pdb ( pdb_1, cards_1, remarks, chain_name, outfile, suffix )
!
!       Need to add fofc_space_cc
!       Chain business
!
        TYPE(pdb),        INTENT(INOUT)          :: pdb_1
        TYPE(cards),      INTENT(IN)             :: cards_1
        CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: remarks
        CHARACTER(LEN=1), INTENT(IN),   OPTIONAL :: chain_name
        CHARACTER(LEN=*), INTENT(IN),   OPTIONAL :: outfile
        CHARACTER(LEN=*), INTENT(IN),   OPTIONAL :: suffix
!       Local variables:
        CHARACTER(len=132)                       :: file_1
        REAL(KIND=wp), DIMENSION(3)              :: xyz
        REAL(KIND=wp), DIMENSION(6)              :: RU
        CHARACTER(len=80)                        :: line
        INTEGER                                  :: outcoor
        INTEGER                                  :: istat
        LOGICAL                                  :: refine_Us
        LOGICAL                                  :: aniso_active
        LOGICAL                                  :: su_present
!       Counters:
        INTEGER                                  :: iat

        IF ( .NOT. PRESENT ( outfile ) ) THEN
            CALL ugtenv ( 'PDBOUT', file_1 )
            IF ( LEN_TRIM ( file_1 ) == 0 ) THEN
                CALL ugtenv ( '-O', file_1 )
                IF ( LEN_TRIM ( file_1 ) == 0 ) THEN
                    CALL die ( 'Failed to detect name of the file in the environment.', 'write_pdb' )
                ENDIF
            ENDIF
        ELSE
            file_1 = outfile
        ENDIF

        outcoor = get_next_io_unit()
        IF ( PRESENT ( suffix ) ) THEN
            file_1 = TRIM ( file_1 ) // TRIM ( suffix )
        ENDIF

        WRITE(*,"(' WRITE_PDB> ', 'Output PDB file= ', A, ' on unit ',A)")&
        TRIM( file_1 ), TRIM( int_to_c ( outcoor ) )

        OPEN ( UNIT=outcoor, FILE=TRIM(file_1), ACCESS='SEQUENTIAL', FORM='FORMATTED',&
               STATUS='UNKNOWN', IOSTAT=istat )

        IF ( istat /= 0 ) THEN
            CALL die ( 'Programming error. Failed to open file :'//TRIM(file_1), 'write_pdb' )
        ENDIF

        IF ( PRESENT ( remarks ) ) THEN
            WRITE ( outcoor, "(A)" ) remarks
        ELSE
            IF ( ALLOCATED ( pdb_1%remark ) ) THEN
                DO iat = 1, SIZE ( pdb_1%remark )
                    IF ( LEN_TRIM ( pdb_1%remark(iat) ) /= 0 ) THEN
                        WRITE ( outcoor, "(A)" ) pdb_1%remark(iat)
                    ENDIF
                ENDDO
                WRITE ( outcoor, "(A)" ) 'REMARK'
            ENDIF
        ENDIF
   
        IF ( ALLOCATED ( pdb_1%scale_cards ) ) THEN
            WRITE ( outcoor, "(A)" ) pdb_1%scale_cards
        ENDIF

        DO iat = 1, SIZE(pdb_1%label)
            IF ( pdb_1%label(iat) /= 'ATOM  ' .AND. pdb_1%label(iat) /= 'HETATM' ) THEN
                WRITE ( *, "(' WRITE_PDB> ', '******* Changing bad label in write pdb: ', A)" ) pdb_1%label(iat)
                pdb_1%label(iat) = 'ATOM  '
            ENDIF
        ENDDO

        refine_Us = INDEX ( cards_1%refinement_mode, 'U' ) > 0
        WRITE(*,"(' WRITE_PDB> ', 'ANISO_IN_PDB=', L1)") pdb_1%aniso_in_pdb
        aniso_active = ALLOCATED ( pdb_1%U )
        su_present   = ALLOCATED ( pdb_1%su )
        WRITE(*,"(' WRITE_PDB> ', 'ANISO_ACTIVE=  ', L1)") aniso_active
        WRITE(*,"(' WRITE_PDB> ', 'SU_PRESENT=  ', L1)")   su_present

        DO iat = 1, SIZE ( pdb_1%xyz )

            IF ( PRESENT ( chain_name ) ) THEN
                pdb_1%chain_name(iat) = chain_name
            ENDIF

!           Need this to overcome privacy:
            xyz = pdb_1%xyz(iat)

!           Move atom names to original position(s):
            IF ( LEN_TRIM ( pdb_1%atom_name(iat) ) /= LEN ( pdb_1%atom_name(iat) ) ) THEN
                pdb_1%atom_name(iat)(2:4) = pdb_1%atom_name(iat)(1:3)
                pdb_1%atom_name(iat)(1:1) = ' ' 
            ENDIF

            IF ( ALLOCATED ( pdb_1%U ) .AND. refine_Us ) THEN

!               Before conversion make sure that atom has been refined anisotropically:
                IF ( COUNT ( pdb_1%refinable(6:11,iat) ) > 0 ) THEN
                    pdb_1%biso(iat) = .CONVERT. ( pi_squared_8 * pdb_1%U(iat) )
                ENDIF

            ENDIF

!           Insertion code and SEGID added OCT 2007; ALTLOC, ELSYMB, CHARGE --MAR 2009:
            WRITE ( line, pdb_1%fmt, IOSTAT=istat) pdb_1%label(iat),               &
                                                   pdb_1%atom_number(iat),         &
                                                   pdb_1%atom_name(iat),           &
                                                   pdb_1%altloc(iat),              &
                                                   pdb_1%residue_name(iat),        &
                                                   pdb_1%chain_name(iat),          &
                                                   pdb_1%residue_number(iat),      &
                                                   pdb_1%insertion_code(iat),      &
                                                   xyz,                            &
                                                   pdb_1%occ(iat),                 &
                                                   pdb_1%biso(iat),                &
                                                   pdb_1%segid(iat),               &
                                                   ADJUSTR(pdb_1%elsymb(iat)),     &
                                                   pdb_1%charge(iat)               

            IF ( istat /= 0 ) THEN
                CALL warn ( 'Minor problems when writing to line: '//TRIM(line),&
                            'write_pdb' )
            ENDIF

            WRITE ( outcoor, '(A)', IOSTAT=istat )  line
            IF ( istat /= 0 ) THEN
                WRITE(*,*) ' istat=', istat
                CALL warn ( 'Problems with line: '//TRIM(line(7:)), 'write_pdb' )
            ENDIF
            
            IF ( su_present ) THEN

                 IF ( ANY ( pdb_1%su(1:11,iat) /= 0.0_sp ) ) THEN

                     IF ( aniso_active ) THEN

                         IF ( ANY ( pdb_1%su(6:8,iat) /= 0.0_sp ) ) THEN
                             pdb_1%su(4,iat) = 8.0_wp * pi_squared * SQRT ( SUM ( pdb_1%su(6:8,iat) ** 2 ) ) / 3.0_wp
                         ENDIF

                     ENDIF

                     WRITE ( line, pdb_1%fmt, IOSTAT=istat) 'SIGATM',                   &
                                                            pdb_1%atom_number(iat),     &
                                                            pdb_1%atom_name(iat),       &
                                                            pdb_1%altloc(iat),          &
                                                            pdb_1%residue_name(iat),    &
                                                            pdb_1%chain_name(iat),      &
                                                            pdb_1%residue_number(iat),  &
                                                            pdb_1%insertion_code(iat),  &
                                                            pdb_1%su(1:3,iat),          &
                                                            pdb_1%su(5,iat),            &
                                                            pdb_1%su(4,iat),            &
                                                            pdb_1%segid(iat),           &
                                                            ADJUSTR(pdb_1%elsymb(iat)), &
                                                            pdb_1%charge(iat)


                    IF ( istat /= 0 ) THEN
                        CALL warn ( 'Minor problems when writing to line: '//TRIM(line),&
                                    'write_pdb' )
                    ENDIF

                    WRITE ( outcoor, '(A)', IOSTAT=istat )  line
                    IF ( istat /= 0 ) THEN
                        WRITE(*,*) ' istat=', istat
                        CALL warn ( 'Problems with line: '//TRIM(line(7:)), 'write_pdb' )
                    ENDIF

                ENDIF

            ENDIF    

            output_aniso:IF ( aniso_active ) THEN

!               Convert back to PDB format:
                RU = pdb_1%U(iat) * aniso_pdb_factor 
                WRITE( line, pdb_1%fmt_aniso, IOSTAT=istat) pdb_1%ulabel(iat),              &
                                                            pdb_1%atom_number(iat),         &
                                                            pdb_1%atom_name(iat),           &
                                                            pdb_1%altloc(iat),              &
                                                            pdb_1%residue_name(iat),        &
                                                            pdb_1%chain_name(iat),          &
                                                            pdb_1%residue_number(iat),      &
                                                            pdb_1%insertion_code(iat),      &
                                                            NINT ( RU ),                    &
                                                            pdb_1%segid(iat),               &
                                                            ADJUSTR ( pdb_1%elsymb(iat) ),  &
                                                            pdb_1%charge(iat)
                IF ( istat /= 0 ) THEN
                    CALL warn ( 'Minor problems when writing to line: '//TRIM(line),&
                                'write_pdb' )
                ENDIF

                WRITE ( outcoor, '(A)', IOSTAT=istat )  line

                IF ( istat /= 0 ) THEN
                    WRITE(*,*) ' istat=', istat
                    CALL warn ( 'Problems with line: '//TRIM(line(7:)), &
                                'write_pdb' )
                ENDIF

                IF ( su_present ) THEN

                    IF ( ANY ( pdb_1%su(6:11,iat) /= 0.0_sp ) ) THEN

                        RU = pdb_1%su(6:11,iat) * aniso_pdb_factor  

                        WRITE( line, pdb_1%fmt_aniso, IOSTAT=istat) 'SIGUIJ',                      &
                                                                    pdb_1%atom_number(iat),        &
                                                                    pdb_1%atom_name(iat),          &
                                                                    pdb_1%altloc(iat),             &
                                                                    pdb_1%residue_name(iat),       &
                                                                    pdb_1%chain_name(iat),         &
                                                                    pdb_1%residue_number(iat),     &
                                                                    pdb_1%insertion_code(iat),     &
                                                                    NINT ( RU ),                   &
                                                                    pdb_1%segid(iat),              &
                                                                    ADJUSTR ( pdb_1%elsymb(iat) ), &
                                                                    pdb_1%charge(iat)

                        IF ( istat /= 0 ) THEN
                            CALL warn ( 'Minor problems when writing to line: '//TRIM(line),&
                                        'write_pdb' )
                        ENDIF

                        WRITE ( outcoor, '(A)', IOSTAT=istat )  line
                        IF ( istat /= 0 ) THEN
                            WRITE(*,*) ' istat=', istat
                            CALL warn ( 'Problems with line: '//TRIM(line(7:)), 'write_pdb' )
                        ENDIF

                    ENDIF
                ENDIF

            ENDIF output_aniso


        ENDDO

        WRITE ( outcoor, '(A)' ) 'END'

!       Close OUTCOOR and check closure:
        CLOSE ( UNIT=outcoor, IOSTAT=istat )

        IF ( istat /= 0 ) THEN
            CALL die('Failed to close UNIT='//TRIM ( int_to_c ( outcoor ) )//' properly',&
                     'write_pdb')
        ENDIF

    END SUBROUTINE write_pdb

    SUBROUTINE allocate_pdb_using_template ( pdb_1, pdb_2 )
!
!       Purpose:
!       =======
!       Allocate pdb structure (pdb_1) using known pdb structure (pdb_2)
!       This is the most straighforward method.
!
        TYPE(pdb), INTENT(INOUT) :: pdb_1
        TYPE(pdb), INTENT(IN)    :: pdb_2

        CALL messag ( ' ', 'allocate_pdb_using_template' )
        
        IF ( ALLOCATED ( pdb_2%p ) ) THEN
            CALL allocate_array(pdb_1%p, SIZE ( pdb_2%p ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%remark ) ) THEN
            CALL allocate_array(pdb_1%remark, SIZE ( pdb_2%remark ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%scale_cards ) ) THEN
            CALL allocate_array(pdb_1%scale_cards, SIZE ( pdb_2%scale_cards ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%label ) ) THEN
            CALL allocate_array(pdb_1%label, SIZE ( pdb_2%label ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%ulabel ) ) THEN
            CALL allocate_array(pdb_1%ulabel, SIZE ( pdb_2%ulabel ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%atom_number ) ) THEN
            CALL allocate_array(pdb_1%atom_number, SIZE ( pdb_2%atom_number))
        ENDIF

        IF ( ALLOCATED ( pdb_2%atom_name ) ) THEN
            CALL allocate_array(pdb_1%atom_name, SIZE ( pdb_2%atom_name ))
        ENDIF

!       ALTLOC added MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_2%altloc ) ) THEN
            CALL allocate_array(pdb_1%altloc, SIZE ( pdb_2%altloc ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%residue_name  ) ) THEN
            CALL allocate_array(pdb_1%residue_name, SIZE ( pdb_2%residue_name ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%chain_name ) ) THEN
            CALL allocate_array(pdb_1%chain_name, SIZE ( pdb_2%chain_name ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%residue_number ) ) THEN
            CALL allocate_array(pdb_1%residue_number, SIZE ( pdb_2%residue_number ))
        ENDIF

!       Insertion code added OCT 2007:
        IF ( ALLOCATED ( pdb_2%insertion_code ) ) THEN
            CALL allocate_array(pdb_1%insertion_code, SIZE ( pdb_2%insertion_code ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%xyz ) ) THEN
            CALL allocate_array(pdb_1%xyz, SIZE ( pdb_2%xyz ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%occ ) ) THEN
            CALL allocate_array(pdb_1%occ, SIZE ( pdb_2%occ ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%biso ) ) THEN
            CALL allocate_array(pdb_1%biso, SIZE ( pdb_2%biso )) 
        ENDIF

!       SEGID added OCT 2007 BVS:
        IF ( ALLOCATED ( pdb_2%segid ) ) THEN
            CALL allocate_array(pdb_1%segid, SIZE ( pdb_2%segid ))
        ENDIF

!       ELSYMB & CHARGE added MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_2%elsymb ) ) THEN
            CALL allocate_array(pdb_1%elsymb, SIZE ( pdb_2%elsymb ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%charge ) ) THEN
            CALL allocate_array(pdb_1%charge, SIZE ( pdb_2%charge ))
        ENDIF

!       INDXYZ added MAR 2009 BVS
        IF ( ALLOCATED ( pdb_2%indxyz ) ) THEN
            CALL allocate_array(pdb_1%indxyz, SIZE ( pdb_2%indxyz ))
        ENDIF

!       Scattering data:
        IF ( ALLOCATED ( pdb_2%atom_type ) ) THEN
            CALL allocate_array(pdb_1%atom_type, SIZE ( pdb_2%atom_type ))
        ENDIF

!       Refinable params:
        IF ( ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL allocate_array(pdb_1%refinable, &
                                SIZE ( pdb_2%refinable, DIM=1), SIZE ( pdb_2%refinable, DIM=2 ))
        ENDIF
 
!       Special treatment for ATOMSF_TABLE: 
        IF ( ALLOCATED ( pdb_2%atomsf_table ) ) THEN
            CALL allocate_sf_database ( pdb_1%atomsf_table, SIZE ( pdb_2%atomsf_table ) )
        ENDIF


!       Anisothermal parameters:
        IF ( ALLOCATED ( pdb_2%u ) ) THEN
            CALL allocate_array(pdb_1%u, SIZE ( pdb_2%u ) )
        ENDIF

!       Residue/chain parameters:
        IF ( ALLOCATED ( pdb_2%residue_start ) ) THEN
            CALL allocate_array(pdb_1%residue_start, SIZE ( pdb_2%residue_start ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%chain_start ) ) THEN
            CALL allocate_array(pdb_1%chain_start, SIZE ( pdb_2%chain_start ))
        ENDIF

        IF ( ALLOCATED ( pdb_2%number_of_atoms_in_this_residue ) ) THEN
            CALL allocate_array ( pdb_1%number_of_atoms_in_this_residue,&
                                  SIZE ( pdb_2%number_of_atoms_in_this_residue ))
        ENDIF

        pdb_1%cell  = pdb_2%cell
        pdb_1%ORT   = pdb_2%ORT
        pdb_1%DEORT = pdb_2%DEORT
        pdb_1%fmt   = pdb_2%fmt
        pdb_1%fmt_aniso = pdb_2%fmt_aniso
        pdb_1%aniso_in_pdb = pdb_2%aniso_in_pdb

    END SUBROUTINE allocate_pdb_using_template

    SUBROUTINE allocate_pdb_using_sizes ( pdb_1, natom, nremarks, nscale_cards, aniso, anom, &
                                          number_of_residues, number_of_chains )
!      
!       Purpose:
!       =======
!       Allocates pdb using basic parameters
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Sep 2003       B.Strokopytov          Original code
!       Oct 2007       B.Strokopytov          SEGID, ATOM_TYPE, INSERTION_CODE arrays added
!       Oct 2007       B.Strokopytov          REFINABLE_OCC array added
!

!       More safe to use INOUT to keep certain parameters:
        TYPE(pdb), INTENT(INOUT)        :: pdb_1
        INTEGER,   INTENT(IN)           :: natom
        INTEGER,   INTENT(IN), OPTIONAL :: nremarks
        INTEGER,   INTENT(IN), OPTIONAL :: nscale_cards
        LOGICAL,   INTENT(IN), OPTIONAL :: aniso
        LOGICAL,   INTENT(IN), OPTIONAL :: anom
        INTEGER,   INTENT(IN), OPTIONAL :: number_of_residues
        INTEGER,   INTENT(IN), OPTIONAL :: number_of_chains
        CHARACTER(LEN=32),         SAVE :: srname='allocate_pdb_using_sizes'
!
        CALL messag ( ' ', srname )
        IF ( PRESENT ( nremarks ) ) THEN
            IF ( nremarks > 0 ) THEN
                CALL allocate_array ( pdb_1%remark, MAX ( 1, nremarks ) )
            ENDIF
        ENDIF

        IF ( PRESENT ( nscale_cards ) ) THEN
            IF ( nscale_cards > 0 ) THEN
                CALL allocate_array(pdb_1%scale_cards, 3)
            ENDIF
        ENDIF
        
        positive_natom:IF ( natom >  0 ) THEN
            CALL allocate_array(pdb_1%label,          natom)
            CALL allocate_array(pdb_1%atom_number,    natom)
            CALL allocate_array(pdb_1%atom_name,      natom)
            CALL allocate_array(pdb_1%altloc,         natom) ! added MAR 2009
            CALL allocate_array(pdb_1%residue_name,   natom)
            CALL allocate_array(pdb_1%chain_name,     natom)
            CALL allocate_array(pdb_1%residue_number, natom)
            CALL allocate_array(pdb_1%insertion_code, natom) ! added
            CALL allocate_array(pdb_1%xyz,            natom)
            CALL allocate_array(pdb_1%occ,            natom)
            CALL allocate_array(pdb_1%biso,           natom)
            CALL allocate_array(pdb_1%segid,          natom) ! added OCT 2007
            CALL allocate_array(pdb_1%elsymb,         natom) ! added MAR 2009
            CALL allocate_array(pdb_1%charge,         natom) ! added MAR 2009
            CALL allocate_array(pdb_1%indxyz,         natom) ! added MAR 2009
            CALL allocate_array(pdb_1%atom_type,      natom) ! This record was missing. BUG corrected 11 OCT 2007 BVS 
            CALL allocate_array(pdb_1%refinable,  11, natom)

            IF ( PRESENT ( aniso ) ) THEN
                IF ( aniso ) THEN
                    CALL allocate_array ( pdb_1%ulabel, natom )
!                   BUG CORRECTED ANISOU was not set. APR/JUL 2009 BVS:
                    pdb_1%ulabel = 'ANISOU'
                    CALL allocate_array ( pdb_1%u,   natom )
                ENDIF
            ENDIF

        ENDIF positive_natom

        IF ( PRESENT ( number_of_residues ) ) THEN
            IF ( number_of_residues > 0 ) THEN
                CALL allocate_array ( pdb_1%residue_start,  number_of_residues + 1 )
                CALL allocate_array ( pdb_1%number_of_atoms_in_this_residue, natom )
            ENDIF
        ENDIF

        IF ( PRESENT ( number_of_chains ) ) THEN
            IF ( number_of_chains > 0 ) THEN
                CALL allocate_array ( pdb_1%chain_start, number_of_chains + 1 )
            ENDIF
        ENDIF

        IF ( PRESENT ( anom ) ) THEN
            IF ( anom ) CALL messag('Program is not ready to use anomalous scattering.', srname)
        ENDIF

!       No allocation of atomsf_table at the moment since
!       its size is generally not known in advance:
        CALL messag('Done.', srname)
    END SUBROUTINE allocate_pdb_using_sizes

    SUBROUTINE deallocate_pdb ( pdb_1 )
        TYPE (pdb), INTENT(INOUT) :: pdb_1
        
!       Parameter record added OCT 2007 BVS:
        IF ( ALLOCATED ( pdb_1%p ) ) THEN
            CALL deallocate_array(pdb_1%p)
        ENDIF

!       PBLOCK added APR 2009 BVS:
        IF ( ALLOCATED ( pdb_1%pblock ) ) THEN
            CALL deallocate_array(pdb_1%pblock)
        ENDIF

        IF ( ALLOCATED ( pdb_1%remark ) ) THEN
            CALL deallocate_array(pdb_1%remark)
        ENDIF

        IF ( ALLOCATED ( pdb_1%scale_cards ) ) THEN
            CALL deallocate_array(pdb_1%scale_cards)
        ENDIF

        IF ( ALLOCATED ( pdb_1%label ) ) THEN
            CALL deallocate_array(pdb_1%label)
        ENDIF

!       ULABEL added MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_1%ulabel ) ) THEN
            CALL deallocate_array(pdb_1%ulabel)
        ENDIF
 
        IF ( ALLOCATED ( pdb_1%atom_number ) ) THEN
            CALL deallocate_array(pdb_1%atom_number)
        ENDIF

        IF ( ALLOCATED ( pdb_1%atom_name ) ) THEN
            CALL deallocate_array(pdb_1%atom_name)
        ENDIF

!       ALTLOC added MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_1%altloc ) ) THEN
            CALL deallocate_array(pdb_1%altloc)
        ENDIF

        IF ( ALLOCATED ( pdb_1%residue_name) ) THEN
            CALL deallocate_array(pdb_1%residue_name)
        ENDIF

        IF ( ALLOCATED( pdb_1%chain_name ) ) THEN
            CALL deallocate_array(pdb_1%chain_name)
        ENDIF

        IF ( ALLOCATED ( pdb_1%residue_number) ) THEN
            CALL deallocate_array(pdb_1%residue_number)
        ENDIF

!       Insertion code added OCT 2007:
        IF ( ALLOCATED ( pdb_1%insertion_code ) ) THEN
            CALL deallocate_array(pdb_1%insertion_code)
        ENDIF

        IF ( ALLOCATED ( pdb_1%xyz ) ) THEN
            CALL deallocate_array(pdb_1%xyz)
        ENDIF

        IF ( ALLOCATED ( pdb_1%occ ) ) THEN
            CALL deallocate_array(pdb_1%occ)
        ENDIF

        IF ( ALLOCATED ( pdb_1%biso ) ) THEN
            CALL deallocate_array(pdb_1%biso)
        ENDIF

!       Segid added OCT 2007:
        IF ( ALLOCATED ( pdb_1%segid ) ) THEN
            CALL deallocate_array(pdb_1%segid)
        ENDIF

!       ELSYMB & CHARGE added MAR 2009:
        IF ( ALLOCATED ( pdb_1%elsymb ) ) THEN
            CALL deallocate_array(pdb_1%elsymb)
        ENDIF

        IF ( ALLOCATED ( pdb_1%charge ) ) THEN
            CALL deallocate_array(pdb_1%charge)
        ENDIF

        IF ( ALLOCATED ( pdb_1%refinable ) ) THEN
            CALL deallocate_array(pdb_1%refinable)
        ENDIF

!       INDXYZ added MAR 2009:
        IF ( ALLOCATED ( pdb_1%indxyz ) ) THEN
            CALL deallocate_array(pdb_1%indxyz)
        ENDIF

!       ATOMSF table part:
        IF ( ALLOCATED ( pdb_1%atom_type ) ) THEN
            CALL deallocate_array(pdb_1%atom_type)
        ENDIF

        IF ( ALLOCATED ( pdb_1%atomsf_table ) ) THEN
            CALL deallocate_sf_database ( pdb_1%atomsf_table )
        ENDIF

!       Anisothermal params rearranged in MAR 2009 BVS:
        IF ( ALLOCATED ( pdb_1%u ) ) THEN
            CALL deallocate_array(pdb_1%u)
        ENDIF
 
!       other... 
        IF ( ALLOCATED ( pdb_1%residue_start ) ) THEN
            CALL deallocate_array(pdb_1%residue_start)
        ENDIF

        IF ( ALLOCATED ( pdb_1%chain_start ) ) THEN
            CALL deallocate_array(pdb_1%chain_start)
        ENDIF

        IF ( ALLOCATED ( pdb_1%number_of_atoms_in_this_residue ) ) THEN
            CALL deallocate_array(pdb_1%number_of_atoms_in_this_residue)
        ENDIF

        IF ( ALLOCATED ( pdb_1%COO ) ) THEN
            CALL deallocate_matrix ( pdb_1%COO )
        ENDIF

        IF ( ALLOCATED ( pdb_1%CSR ) ) THEN
            CALL deallocate_matrix ( pdb_1%CSR )
        ENDIF

!       ALL_ATOM_IMAGES added APR 2009 BVS:
        IF ( ALLOCATED ( pdb_1%all_atom_images ) ) THEN
            CALL deallocate_array(pdb_1%all_atom_images)
        ENDIF

        IF ( ALLOCATED ( pdb_1%su ) ) THEN
            CALL deallocate_array(pdb_1%su)
        ENDIF

        IF ( ALLOCATED ( pdb_1%w ) ) THEN
            CALL deallocate_array(pdb_1%w)
        ENDIF

    END SUBROUTINE deallocate_pdb

    SUBROUTINE residue_location ( pdb_1 )
!
!       Purpose:
!       =======
!       a)Locates start and end of chains
!       b)Locates start and end of residues
! 
!       Date:           Programmer:          History of changes:
!       ====            ==========           ==================
!       Oct 2003        B.Strokopytov        Original code
!       Jul 2007        B.Strokopytov        Corrected residue start definition:
!                                            new chain means the start of new residue
!                                            
        TYPE(pdb), INTENT(INOUT) :: pdb_1
!       Local variables:
        INTEGER                  :: istart
        INTEGER                  :: start_chain
        INTEGER                  :: natom
        INTEGER                  :: previous_residue_number
!       Char variables:
        CHARACTER(KIND=1)        :: previous_insertion_code
        CHARACTER(KIND=1)        :: previous_chain_name
!       Counters:
        INTEGER                  :: i
        IF ( .NOT. ALLOCATED ( pdb_1%residue_number ) ) THEN
            CALL die('Programming error. PDB file has not been properly initialized.',&
                     'residue_location')
        ELSE
            natom  = SIZE ( pdb_1%residue_number )
        ENDIF

        previous_chain_name     = '-'
        previous_residue_number = -100
        previous_insertion_code = ' '
        istart                  =  0
        start_chain             =  0

        DO i = 1, natom

            IF ( pdb_1%residue_number(i) /= previous_residue_number &
            .OR. pdb_1%insertion_code(i) /= previous_insertion_code &
            .OR. pdb_1%chain_name(i) /= previous_chain_name ) THEN
                istart = istart + 1
                previous_residue_number = pdb_1%residue_number(i)
                previous_insertion_code = pdb_1%insertion_code(i)
            ENDIF

            IF ( pdb_1%chain_name(i) /= previous_chain_name ) THEN
                                                                       
                start_chain = start_chain + 1
                previous_chain_name = pdb_1%chain_name(i)
            ENDIF

        ENDDO

!       Allocate memory for residue_start pointers:
        ALLOCATE(pdb_1%residue_start(istart + 1))
        ALLOCATE(pdb_1%number_of_atoms_in_this_residue(natom))
        ALLOCATE(pdb_1%chain_start(start_chain + 1) )

!       Initialize again:
        istart      = 0
        start_chain = 0
        previous_chain_name = '-'
        previous_insertion_code = ' '
!        WRITE(*,*) ' pdb_1%chain_name=', pdb_1%chain_name
        DO i = 1, natom
            IF ( pdb_1%residue_number(i) /= previous_residue_number &
            .OR. pdb_1%insertion_code(i) /= previous_insertion_code &
            .OR. pdb_1%chain_name(i) /= previous_chain_name ) THEN
                istart                      = istart + 1
!               Long standing correction ("-1") added MAR 2009 BVS:
                pdb_1%residue_start(istart) = i - 1
                previous_residue_number     = pdb_1%residue_number(i)
                previous_insertion_code     = pdb_1%insertion_code(i)
            ENDIF

            IF ( pdb_1%chain_name(i) /= previous_chain_name ) THEN
                start_chain = start_chain + 1
                pdb_1%chain_start(start_chain) = i - 1
                previous_chain_name = pdb_1%chain_name(i)
            ENDIF
        ENDDO

!       number_of_atoms_in_residue = start(i+1) - start(i). Hopefully this is obvious:
        pdb_1%residue_start(istart + 1)     = natom
        pdb_1%chain_start(start_chain + 1)  = natom

!       Loop over number of residues:
        DO i = 1, istart
!           Corrected due to new scheme for RESIDUE_START:
            pdb_1%number_of_atoms_in_this_residue(pdb_1%residue_start(i)+1:pdb_1%residue_start(i + 1)) = &
            pdb_1%residue_start(i + 1) - pdb_1%residue_start(i)
        ENDDO

    END SUBROUTINE residue_location

    SUBROUTINE refinable_occupancies(pdb_1)
!
!       Purpose:
!       =======
!       Tries to figure out which atoms have refinable occupnacies.
!       
!       Date:            Programmer:            History of changes:
!       ====             ==========             ==================
!       Oct 2007         B.Strokopytov          Original code
!
!       Note:
!       ====
!       This is rather preliminary version.
!       Need to have provisions from PDB file itself and possibly from user cards.
!       But `HETATM' mechanism seems very universal except for multiple conformations.
!
        TYPE(pdb), INTENT(INOUT) :: pdb_1
!       Counters:
        INTEGER                  :: iat

        DO iat = 1, SIZE ( pdb_1 )

!           All atoms having non-unity occupancies are supposed to be refinable:
            IF ( pdb_1%occ(iat) /= 1.0_wp ) THEN
                pdb_1%refinable(5,iat) = .TRUE.
            ENDIF

!           All atoms having special labels are supposed to have refinable occupancies:
!           FIXME in this case we must restrain occupancies to have sum equal unity, e.g.
!           q(OG A) + q(OG B) + q(OG C) = 1:
            IF ( INDEX ('ABC',  pdb_1%atom_name(iat)(4:4) )  > 0 ) THEN
                 pdb_1%refinable(5,iat) = .TRUE.
                 CYCLE
            ENDIF

!           All single atoms (most likely solvent or heavy metals):            
            IF ( pdb_1%number_of_atoms_in_this_residue(iat) == 1 ) THEN
                pdb_1%refinable(5,iat) = .TRUE. 
                CYCLE
            ENDIF

!           All hetero atoms / residues (most useful for bound molecules):
            IF ( pdb_1%label(iat)(1:4) == 'HETA' ) THEN
                pdb_1%refinable(5,iat) = .TRUE.
                CYCLE
            ENDIF
!
            IF ( pdb_1%atom_name(iat)(1:3) == 'HOH' .OR.  pdb_1%atom_name(iat)(1:3) == 'WAT' ) THEN
                pdb_1%refinable(5,iat) = .TRUE.
                CYCLE
            ENDIF


        ENDDO

!       Global debugging if necessary:
        IF ( debug > 2 ) THEN

            DO iat = 1, SIZE ( pdb_1 )
                WRITE(*,"(' REFINABLE_OCCUPANCIES> ', 'label= ', A, ' atom name= ', A, ' refinable occ= ', L1)")&
                pdb_1%label(iat), pdb_1%atom_name(iat), pdb_1%refinable(5,iat)
            ENDDO

        ENDIF

    END SUBROUTINE refinable_occupancies

    FUNCTION get_bmin ( pdb_1 )
!
!       Purpose:
!       =======
!       Calculates minimal Biso taking into account gaussian B's for scattering factors:
!
!       Date:            Programmer:           History of changes:
!       ====             ==========            ==================
!       Oct 2007         B.Strokopytov         Original code
!
!       Note:
!       ====
!       Could be almost one-liner but would look awful.
!       FIXME: need some provision for anisotropic u11, etc. values
!
        REAL(KIND=wp)               :: get_bmin
        TYPE(pdb),       INTENT(IN) :: pdb_1
!       Local variables:
        INTEGER                     :: ityp
        INTEGER                     :: ngauss 
        REAL(KIND=wp)               :: effective_biso
        REAL(KIND=wp)               :: min_biso
        LOGICAL                     :: aniso_active
        TYPE(matrix)                :: UMAT
        REAL(KIND=wp), DIMENSION(3) :: ueigen
        REAL(KIND=wp)               :: bmin_gauss
        REAL(KIND=wp)               :: ratio
        REAL(KIND=wp)               :: min_ratio
!       Counters:
        INTEGER                     :: iat
        CHARACTER(LEN=32)           :: srname
        
        srname = 'get_bmin'

!       Basic checks:
        IF ( .NOT. ALLOCATED ( pdb_1 ) ) THEN
            CALL die ('Programming error. PDB_1 has not been allocated.', srname )
        ELSE IF ( .NOT. ALLOCATED ( pdb_1%atomsf_table ) ) THEN
            CALL die ('Programming error. PDB_1%ATOMSF_TABLE  has not been allocated.',& 
                       srname )
        ENDIF
        
        min_biso  = 10.0_wp ** 20
        min_ratio = 10.0_wp ** 20

!       Check whether aniso thermal factors have been allocated:
        aniso_active = ALLOCATED ( pdb_1%U )

        DO iat = 1, SIZE ( pdb_1 )

            ityp   =  pdb_1%atom_type(iat)
            ngauss =  pdb_1%atomsf_table%record(ityp)%ngauss

            IF ( ALLOCATED ( pdb_1%atomsf_table%record(ityp) ) ) THEN

!               This is likely to be equal to zero since normally we use 5-gaussian approximation:
                bmin_gauss = MINVAL ( pdb_1%atomsf_table%record(ityp)%b(1:ngauss) )

!               Check whether we deal with anisotropic thermal factors:
                IF ( aniso_active ) THEN
!                   Check U array (just first three values):
                    IF ( ANY ( pdb_1%U(iat) > 0.0_wp ) ) THEN
                         UMAT = pdb_1%U(iat)
                         ueigen = eigen_values ( UMAT )

!                        Assuming that the third eigenvalue is the smallest:
                         ratio  = ueigen(1) / ueigen(3)

!                        Just a check that CORRANI function works correctly:
                         IF ( ratio < MIN_ADP_RATIO - EPSILON (1.0_wp) ) THEN
                             WRITE(*,*) ' Problems with atom=', iat
                             CALL print_mat(UMAT)
                             WRITE(*,*) ' eigvals=',ueigen
                             WRITE(*,*) ' sum=', SUM(ueigen)
!                             pdb_1%U(iat) = pdb_1%U(iat) .CORRANI. 0.25_wp
                         ENDIF

!                        Monitor ratios:
                         min_ratio = MIN ( ratio, min_ratio )

                         effective_Biso = 8.0_wp * pi_squared * MINVAL ( ueigen ) + bmin_gauss
                   
                    ELSE

!                       Mixed refinement --> fall back to calculations using Biso:
                        effective_Biso = pdb_1%biso(iat) + bmin_gauss

                    ENDIF

                ELSE
                    effective_Biso = pdb_1%biso(iat) + bmin_gauss
                ENDIF

            ELSE

                CALL die ('Programming error. PDB_1%ATOMSF_TABLE%RECORD has not been allocated for atom #'&
                          //TRIM ( int_to_c ( iat ) ), srname )

            ENDIF

            min_biso = MIN ( min_biso, effective_Biso )

        ENDDO

!       Monitor minimal ratio:
        IF ( aniso_active ) THEN
            WRITE(*,"(' GET_BMIN> ', 'Smallest eigenvalue ratio=', ES9.2)") min_ratio
            IF ( min_ratio < 0.25_wp ) CALL warn ('This small ratio requires some attention.', srname)
        ENDIF

        get_bmin = min_biso   

    END FUNCTION get_bmin

    SUBROUTINE uniquefy_atom_types
    END SUBROUTINE uniquefy_atom_types
    
    SUBROUTINE matrix_pointers ( pdb_2, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes )
        TYPE(pdb),             INTENT(INOUT) :: pdb_2
        LOGICAL,               INTENT(IN)    :: refine_xyz
        LOGICAL,               INTENT(IN)    :: refine_Biso
        LOGICAL,               INTENT(IN)    :: refine_occ
        LOGICAL,               INTENT(IN)    :: refine_Us
        LOGICAL, DIMENSION(:), INTENT(IN)    :: polar_axes
!       For various tests:
        INTEGER                              :: np
        LOGICAL                              :: mixed
        LOGICAL                              :: aniso_active
        LOGICAL                              :: Us_present
!       Counters:
        INTEGER                              :: l
        INTEGER                              :: indxyz
        INTEGER                              :: iat
        CHARACTER(LEN=32),       SAVE        :: srname = 'matrix_pointers'

!       Biso test:

!       Allocate and initialize:
        CALL allocate_array ( pdb_2%refinable, 11, SIZE ( pdb_2 ) )
!       Initialize:
        pdb_2%refinable = .FALSE.

!       Figure out difficulties with thermal factors refinement:
        mixed = refine_Biso .AND. refine_Us
        aniso_active = ALLOCATED ( pdb_2%U )

!       Refinable occupancies should be determined outside this routine:
        IF ( refine_occ ) THEN
            CALL refinable_occupancies ( pdb_2 )
        ENDIF


        DO iat = 1, SIZE ( pdb_2 ) 

!           Currently no refinement for hydrogens:
            IF ( pdb_2%elsymb(iat) == 'H ' ) THEN

!                Skipping hydrogens altogether: FIXME
                 pdb_2%occ(iat) = 0.0_wp
            ENDIF

            IF ( refine_xyz ) THEN
                pdb_2%refinable(1:3, iat) = .TRUE.
            ENDIF



!           We need (at least) an allocated U array to proceed with aniso refinement:
            aniso_in_play:IF ( aniso_active ) THEN

!               Check whether Us have been read (appropriate record present):
                Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )

                IF ( mixed ) THEN
                    IF ( refine_Us .AND. Us_present ) THEN
                        pdb_2%refinable(6:11, iat) = .TRUE.
                    ENDIF
                ELSE
                    IF ( refine_Us ) THEN
                        pdb_2%refinable(6:11, iat) = .TRUE.
                    ENDIF
                ENDIF

!               Mixed refinement (B and U refinement) declared and we have missing U record:
                IF ( mixed ) THEN

                    IF ( refine_biso .AND. .NOT. Us_present ) THEN
                        pdb_2%refinable(4, iat) = .TRUE.
                    ENDIF

                ENDIF

            ELSE aniso_in_play

!               Normal B-refinement:
                IF ( refine_biso ) THEN
                    pdb_2%refinable(4, iat) = .TRUE.

!                   Check for impossible situation:
                    IF ( COUNT ( pdb_2%refinable(6:11,iat) ) > 0 ) THEN
                        WRITE(*,"(1X,L1,2X,6L1)") pdb_2%refinable(4, iat), pdb_2%refinable(6:11,iat)
                        CALL die('Programming error: can''t refine both isotropic and anisotropic &
                                 &values for the same atom.', 'matrix_pointers')
                    ENDIF 
                ENDIF

            ENDIF aniso_in_play

            IF ( pdb_2%elsymb(iat) == 'H ' ) THEN
            ENDIF

!           Atoms with zero occupancies do not deserve to be refined ==> set all to false:
            IF ( pdb_2%occ(iat) == 0.0_wp ) THEN
                pdb_2%refinable(1:11,iat) = .FALSE.
            ENDIF


        ENDDO

!       Temporary measure in the absence of restraints:
        IF ( SIZE ( pdb_2 )  > 3  ) THEN

!           In case of small number of atoms for tests ignore this:
            IF ( refine_xyz ) THEN
                pdb_2%refinable(1:3,1) = .NOT. polar_axes
            ENDIF

        ENDIF

!       Check:
        IF ( debug > 1000 ) THEN

            WRITE(*,"(' MATRIX_POINTERS> ', ' refine_xyz= ', L1, ' refine_biso= ', L1, &
                      ' refine_occ= ', L1, ' refune_Us=', L1, ' mixed=',L1, ' aniso=', L1 )") &
            refine_xyz, refine_biso, refine_occ, refine_Us, mixed, aniso_active
        ENDIF

        l = 0
        DO iat = 1, SIZE ( pdb_2 )
            indxyz = l
            pdb_2%indxyz(iat) = l
!            IF ( debug > 1000 ) THEN
                WRITE(*,"(' MATRIX_POINTERS> ', ' atom=', I6, ' pdb atom no.=',I6, ' atom name=', A4,&
                                               &' indxyz=', I8, ' refinable= ', 11L1, 5X, 'elsymb=',A4 )") &
                iat, pdb_2%atom_number(iat), pdb_2%atom_name(iat), indxyz, pdb_2%refinable(1:11, iat), &
                pdb_2%elsymb(iat)//pdb_2%charge(iat)
!            ENDIF
            l = l + COUNT ( pdb_2%refinable(1:11, iat) )
!            IF ( l > 3491 + 50 ) THEN 
!                WRITE(*,*) iat, pdb_2%atom_name(iat); STOP 'bad atom'
!            ENDIF
        ENDDO

        np = COUNT ( pdb_2%refinable )
        
        WRITE(*,"(' MATRIX_POINTERS> ', 'Total number of parameters= ', A)") &
        TRIM ( int_to_c ( np ) )

        IF ( np <= 0 ) THEN
            CALL die ('O-o-o-p-s. Nothing to refine.', srname)
        ENDIF

    END SUBROUTINE matrix_pointers

    SUBROUTINE apply_shifts(p, pdb_2, damp, final)
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: p
        TYPE(pdb),                   INTENT(INOUT) :: pdb_2
        REAL(KIND=wp),               INTENT(IN)    :: damp
        LOGICAL,                     INTENT(IN)    :: final
!       Local variables:
        REAL(KIND=wp), DIMENSION(3)                :: shift_xyz
        TYPE(vector)                               :: dxyz
!       Statistics:
        REAL(KIND=wp)                              :: max_shift_biso
        REAL(KIND=wp)                              :: max_shift_occ
        REAL(KIND=wp)                              :: max_shift_xyz
        REAL(KIND=wp)                              :: max_shift_U
!       Uiso shift = ( d U11 + d U22 + d U33 ) / 3:
        REAL(KIND=wp)                              :: Uiso_shift 
        REAL(KIND=wp)                              :: b_crit 
!       Counters:
        INTEGER                                    :: iat
        INTEGER                                    :: indxyz
        INTEGER                                    :: k
        INTEGER                                    :: l
        INTEGER                                    :: ll
        CHARACTER(LEN=32),                    SAVE :: srname = 'apply_shifts'

        CALL messag(' ', srname)
        WRITE(*,"(' APPLY_SHIFTS> ', 'damp=', F6.2)") damp
 
        IF ( MINVAL ( pdb_2%cell(1:3) ) < 16.0_wp ) THEN
            b_crit = 100.0_wp
        ELSE
            b_crit = 200.0_wp
        ENDIF

!       Statistics:
        max_shift_xyz  = 0.0_wp
        max_shift_biso = 0.0_wp
        max_shift_occ  = 0.0_wp
        max_shift_U    = 0.0_wp

!       Initialize pointer:
        indxyz = 0

!       Apply shifts:
        DO iat = 1, SIZE ( pdb_2 )

!           Initialise xyz shift:
            shift_xyz = 0.0_wp

            l = 0
            DO k = 1, 3
                IF ( pdb_2%refinable(k,iat) ) THEN
                    l = l + 1
                    max_shift_xyz = MAX ( max_shift_xyz,  MIN ( ABS ( damp * p(indxyz+l) ), MAX_ABS_XYZ_SHIFT ) )
                    shift_xyz(k) = p(indxyz+l)
                ENDIF
            ENDDO

!           Max abs shift is normally limited to 0.5 Angstroms:
            DO k = 1, 3
                IF ( ABS ( shift_xyz(k) * damp ) > MAX_ABS_XYZ_SHIFT ) THEN
!                   Cannot get in here if shift_xyz(k) == 0:
                    shift_xyz(k) = MAX_ABS_XYZ_SHIFT * shift_xyz(k) / ABS ( damp * shift_xyz(k) )
                ENDIF
            ENDDO

!           Convert to vector:
            dxyz = shift_xyz

!           Apply shift:
            pdb_2%xyz(iat) = pdb_2%xyz(iat) + damp * dxyz

            IF ( pdb_2%refinable(4,iat) ) THEN

                l = l + 1
                max_shift_biso = MAX ( max_shift_biso,  MIN ( ABS ( damp * p(indxyz+l) ), MAX_ABS_BISO_SHIFT ) )
                IF ( ABS ( damp * p(indxyz+l) ) > MAX_ABS_BISO_SHIFT ) THEN
                    p(indxyz+l) = MAX_ABS_BISO_SHIFT * p(indxyz+l) / ABS ( damp * p(indxyz+l) )
                ENDIF

                pdb_2%biso(iat) = pdb_2%biso(iat) + damp * p(indxyz+l)

                IF ( final ) THEN

!                   Need to avoid negative Biso values:
                    pdb_2%biso(iat) = MAX ( 1.0_wp, pdb_2%biso(iat) )

!                   Need to avoid too high B-values:
                    pdb_2%biso(iat) = MIN ( b_crit, pdb_2%biso(iat) )  

                ENDIF
            ENDIF

            IF ( pdb_2%refinable(5,iat) ) THEN
                l = l + 1
                pdb_2%occ(iat) = pdb_2%occ(iat) + damp * p(indxyz+l)
                max_shift_occ = MAX ( max_shift_occ,  ABS ( damp * p(indxyz+l) ) )

!               Avoid negative or zero occupancies and too positive occupancies:
                IF ( final ) THEN
                    pdb_2%occ(iat) = MAX ( 0.01_wp, pdb_2%occ(iat) )
                    pdb_2%occ(iat) = MIN ( 1.0_wp, pdb_2%occ(iat) )
                ENDIF
            ENDIF

            aniso6:IF ( ANY ( pdb_2%refinable(6:11,iat ) ) ) THEN
                Uiso_shift = 0.0_wp

!               Keep L index intact:
                ll = l 
                DO k = 1, 3
                    IF ( pdb_2%refinable(k+5,iat) ) THEN
                        ll = ll + 1
                        Uiso_shift = Uiso_shift + p(indxyz+ll)
                    ENDIF
                ENDDO

                Uiso_shift = Uiso_shift / 3.0_wp

!               scale = 20.0_wp / (8.0_wp * pi ** 2) * p(indxyz+l) / ABS ( damp * Uiso_shift )
                DO k = 1, 6
                    IF ( pdb_2%refinable(k+5,iat) ) THEN
                        l = l + 1
                        IF ( ABS ( damp * Uiso_shift ) > MAX_ABS_BISO_SHIFT / ( 8.0_wp * pi ** 2 ) ) THEN
                            p(indxyz+l) =  MAX_ABS_BISO_SHIFT / ( 8.0_wp * pi ** 2) &
                                        /  ABS ( damp * Uiso_shift ) * p(indxyz+l ) 
                        ENDIF

                        max_shift_U = MAX ( max_shift_U,  MIN ( ABS ( damp * p(indxyz+l) ), &
                                            MAX_ABS_BISO_SHIFT / ( 8.0_wp * pi ** 2 ) ) )

                        pdb_2%U(iat)%U(k) = pdb_2%U(iat)%U(k) + damp * p(indxyz+l)
                    ENDIF
                ENDDO

                IF ( final ) THEN

!                   This is not entirely correct but effective protection against negative shifts:
                    pdb_2%U(iat)%U(1:3) = MAX ( 1.0/ ( 8.0 * pi ** 2), pdb_2%U(iat)%U(1:3) )

!                   Check and correct (if necessary) ADP ellipsoid :
                    pdb_2%U(iat) = pdb_2%U(iat) .CORRANI. MIN_ADP_RATIO

                ENDIF
            ENDIF aniso6

!           Increment shift/gradient/matrix pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )
        ENDDO

        CALL messag('Shifts have been applied...', srname)
        IF ( max_shift_xyz > 0.1_wp ** 5 .OR.  &
             max_shift_biso > 0.1_wp ** 5 .OR. &
             max_shift_occ  > 0.1_wp ** 5 .OR. &
             max_shift_U > 0.1_wp ** 5 ) THEN
            WRITE(*,"(' APPLY_SHIFTS> ', 'max xyz  shift=', F10.5)") max_shift_xyz
            WRITE(*,"(' APPLY_SHIFTS> ', 'max biso shift=', F10.5)") max_shift_biso
            WRITE(*,"(' APPLY_SHIFTS> ', 'max occ  shift=', F10.5)") max_shift_occ
            WRITE(*,"(' APPLY_SHIFTS> ', 'max U    shift=', F10.5)") max_shift_U
        ELSE
            WRITE(*,"(' APPLY_SHIFTS> ', 'max xyz  shift=', ES9.2)") max_shift_xyz
            WRITE(*,"(' APPLY_SHIFTS> ', 'max biso shift=', ES9.2)") max_shift_biso
            WRITE(*,"(' APPLY_SHIFTS> ', 'max occ  shift=', ES9.2)") max_shift_occ
            WRITE(*,"(' APPLY_SHIFTS> ', 'max U    shift=', ES9.2)") max_shift_U
        ENDIF

    END SUBROUTINE apply_shifts

    SUBROUTINE get_xold(p, pdb_2)
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: p
        TYPE(pdb),                   INTENT(IN)    :: pdb_2
!       Local variables:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Counters:
        INTEGER                                    :: iat
        INTEGER                                    :: indxyz
        INTEGER                                    :: k
        INTEGER                                    :: l
        CHARACTER(LEN=32),                    SAVE :: srname = 'get_xold'

        CALL messag(' ', srname)

!       Checkz:
        IF ( COUNT ( pdb_2%refinable ) /= SIZE ( p ) ) THEN
            CALL die('Programming error. Inconsistent size of P array.', srname)
        ENDIF

!       Initialize parameter pointer:
        indxyz = 0

!       Apply shifts:
        DO iat = 1, SIZE ( pdb_2 )

!           Convert to array:
            xyz = pdb_2%xyz(iat)

            l = 0
            DO k = 1, 3
                IF ( pdb_2%refinable(k,iat) ) THEN
                    l = l + 1
                    p(indxyz+l) = xyz(k)
                ENDIF
            ENDDO

            IF ( pdb_2%refinable(4,iat) ) THEN
                l = l + 1
                p(indxyz+l) = pdb_2%biso(iat)
            ENDIF

            IF ( pdb_2%refinable(5,iat) ) THEN
                l = l + 1
                p(indxyz+l) = pdb_2%occ(iat)
            ENDIF

            IF ( ANY ( pdb_2%refinable(6:11,iat) ) ) THEN
                DO k = 1, 6
                    IF ( pdb_2%refinable(k+5,iat) ) THEN
                        l = l + 1
                        p(indxyz+l) = pdb_2%U(iat)%u(k)   ! FIXME
                    ENDIF
                ENDDO 
            ENDIF

!           Increment shift/gradient/matrix pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )

        ENDDO

        CALL messag('Old model parmeters have been collected...', srname)

    END SUBROUTINE get_xold

    SUBROUTINE set_xold(p, pdb_2)
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: p
        TYPE(pdb),                   INTENT(INOUT) :: pdb_2
!       Local variables:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Counters:
        INTEGER                                    :: iat
        INTEGER                                    :: indxyz
        INTEGER                                    :: k
        INTEGER                                    :: l
        CHARACTER(LEN=32),                    SAVE :: srname = 'set_xold'

        CALL messag(' ', srname)

!       Checkz:
        IF ( COUNT ( pdb_2%refinable ) /= SIZE ( p ) ) THEN
            CALL die('Programming error. Inconsistent size of P array.', srname)
        ENDIF

!       Initialize parameter pointer:
        indxyz = 0

!       Apply shifts:
        DO iat = 1, SIZE ( pdb_2 )

!           Initialize xyz:
            xyz = pdb_2%xyz(iat)

            l = 0
            DO k = 1, 3

!               Refinable coordinates will be replaced with new values:
                IF ( pdb_2%refinable(k,iat) ) THEN
                    l = l + 1
                    xyz(k) =  p(indxyz+l) 
                ENDIF
            ENDDO

!           Set final params:
            pdb_2%xyz(iat) = xyz

!           Do not change B-value if it is unrefinable:
            IF ( pdb_2%refinable(4,iat) ) THEN
                l = l + 1
                pdb_2%biso(iat) = p(indxyz+l)
            ENDIF

            IF ( pdb_2%refinable(5,iat) ) THEN
                l = l + 1
                pdb_2%occ(iat) =  p(indxyz+l)
            ENDIF

            IF ( ANY ( pdb_2%refinable(6:11,iat) ) ) THEN
                DO k = 1, 6
                    IF ( pdb_2%refinable(k+5,iat) ) THEN
                        l = l + 1
                        pdb_2%U(iat)%u(k) = p(indxyz+l)
                    ENDIF
                ENDDO
            ENDIF

!           Increment shift/gradient/matrix pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )
        ENDDO

        CALL messag('PDB parameters have been reset...', srname)
    END SUBROUTINE set_xold

    SUBROUTINE count_atom_types ( pdb_2, list_of_atom_types )
        TYPE(pdb),                          INTENT(IN)    :: pdb_2
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: list_of_atom_types
!       Local variables:
        INTEGER                                           :: ityp
        INTEGER                                           :: number_of_sf_types
        INTEGER, DIMENSION(:), ALLOCATABLE                :: bins
!       Counters:
        INTEGER                                           :: i
        INTEGER                                           :: iat
        INTEGER                                           :: l

        CALL allocate_array ( bins, SIZE ( pdb_2%atomsf_table ) )
        bins = 0
        DO iat = 1, SIZE ( pdb_2 )
            ityp = pdb_2%atom_type(iat)
            bins(ityp) = bins(ityp) + 1
        ENDDO        
        number_of_sf_types = COUNT ( bins > 0 )
        CALL allocate_array ( list_of_atom_types, number_of_sf_types )
        l = 0
        DO i = 1, SIZE ( bins )
            IF ( bins(i) > 0 ) THEN
                l = l + 1
                list_of_atom_types(l) = i
            ENDIF
        ENDDO
        IF ( debug > 100 ) THEN
            WRITE(*,*) ' list_of_atom_types=', list_of_atom_types
        ENDIF

        CALL deallocate_array(bins)
    END SUBROUTINE count_atom_types

    SUBROUTINE prepare_vector_for_pairs ( vector, pdb_2, pairs, ityp, jtyp )
        INTEGER(KIND=eb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector
        TYPE(pdb),                                   INTENT(IN)    :: pdb_2
        INTEGER, DIMENSION(:,:),                     INTENT(IN)    :: pairs
        INTEGER,                                     INTENT(IN)    :: ityp
        INTEGER,                                     INTENT(IN)    :: jtyp
!       Local variables:
        INTEGER                                                    :: iat
        INTEGER                                                    :: jat
        LOGICAL                                                    :: ok
!       Counters:
        INTEGER(KIND=eb)                                           :: pair
        INTEGER(KIND=eb)                                           :: n
        INTEGER                                                    :: run

        DO run = 0, 1
            n = 0
            DO pair = 1, SIZE ( pairs, DIM=2 )
                iat = pairs(1,pair)
                jat = pairs(2,pair)
!               It's important to take swap in atom types into account, i.e. 'CN' or 'NC' pair is indistinguishable:
                ok =  ( pdb_2%atom_type(iat) == ityp .AND. pdb_2%atom_type(jat) == jtyp ) .OR. &
                      ( pdb_2%atom_type(iat) == jtyp .AND. pdb_2%atom_type(jat) == ityp )

                IF ( ok ) THEN
                    n = n + 1
                    IF ( run == 1 ) THEN
!                       Store location of suitable pair:
                        vector(n) = pair
                    ENDIF
                ENDIF
            ENDDO
!           Counted number of suitable pairs, time to allocate:
            IF ( run == 0 ) THEN
                CALL allocate_array(vector, n)
            ENDIF
        ENDDO

    END SUBROUTINE prepare_vector_for_pairs

    SUBROUTINE prepare_vector_of_same_type_atoms( vector, pdb_2, ityp )
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: vector
        TYPE(pdb),                          INTENT(IN)    :: pdb_2
        INTEGER,                            INTENT(IN)    :: ityp
!       Local variables:
        INTEGER                                           :: iat
        LOGICAL                                           :: ok
!       Counters:
        INTEGER                                           :: n
        INTEGER                                           :: run

        DO run = 0, 1
            n = 0
            DO iat = 1, SIZE ( pdb_2 )
                ok =  ( pdb_2%atom_type(iat) == ityp ) 
                IF ( ok ) THEN
                    n = n + 1
                    IF ( run == 1 ) THEN
!                       Store the location of suitable atom:
                        vector(n) = iat 
                    ENDIF
                ENDIF
            ENDDO

!           Counted number of atoms having certain ITYP, time to allocate:
            IF ( run == 0 ) THEN
                CALL allocate_array(vector, n)
            ENDIF
        ENDDO

    END SUBROUTINE prepare_vector_of_same_type_atoms

!  Purpose:
!  ========
!  Reads, creates and deallocate  derived type sf_database
!  (database of scattering factors)
! 
!  Date:          Programmer:            History of changes:
!  ====           ==========             ==================
!  Oct 2003       B.Strokopytov          Original code
!  Oct 2007       B.Strokopytov          ATOMSF, SF_DATABASE types have been implemented 
!  Mar 2009       B.Strokopytov          Cosmetic changes. Need to become independent of CCP4 format. GLobal task, isn't it?
!
    FUNCTION allocated_atomsf ( atomsf_1 )
        LOGICAL                  :: allocated_atomsf
        TYPE(atomsf), INTENT(IN) :: atomsf_1 
        allocated_atomsf = ALLOCATED ( atomsf_1%a ) .AND. ALLOCATED ( atomsf_1%b )
    END FUNCTION allocated_atomsf

    SUBROUTINE allocate_atomsf ( atomsf_1, ngauss )
        TYPE(atomsf), INTENT(INOUT) :: atomsf_1
        INTEGER,      INTENT(IN)    :: ngauss
        IF ( ALLOCATED (  atomsf_1) ) THEN
            CALL warn( 'Programming error. THis record has been allocated already.', 'allocate_atomsf')
        ELSE
            CALL allocate_array ( atomsf_1%a, ngauss )
            CALL allocate_array ( atomsf_1%b, ngauss )
        ENDIF
    END SUBROUTINE

    SUBROUTINE deallocate_atomsf ( atomsf_1 )
        TYPE(atomsf), INTENT(INOUT) :: atomsf_1
        IF ( ALLOCATED ( atomsf_1 ) ) THEN
            CALL deallocate_array ( atomsf_1%a )
            CALL deallocate_array ( atomsf_1%b )
        ELSE
            CALL warn ('Programming error. Nothing to deallocate.', 'allocate_atomsf')
        ENDIF
    END SUBROUTINE deallocate_atomsf

    FUNCTION size_of_atomsf ( atomsf_1 )
        INTEGER                  :: size_of_atomsf
        TYPE(atomsf), INTENT(IN) :: atomsf_1
        size_of_atomsf = SIZE ( atomsf_1%a )
    END FUNCTION size_of_atomsf

    SUBROUTINE copy_atomsf ( atomsf_1, atomsf_2 )
        TYPE(atomsf), INTENT(INOUT) :: atomsf_1
        TYPE(atomsf), INTENT(IN)    :: atomsf_2
!       Checkz:
        IF ( .NOT. ALLOCATED ( atomsf_2 ) ) THEN
            CALL die ('Programming error. ATOMSF_2 has not been allocated...', 'copy_atomsf' )
        ENDIF
        atomsf_1%id = atomsf_2%id
        atomsf_1%ngauss = atomsf_2%ngauss
        atomsf_1%atomic_charge = atomsf_2%atomic_charge ! Bug corrected MAR 2009 BVS
        atomsf_1%number_of_electrons = atomsf_2%number_of_electrons
!       Allocatable part:
        CALL allocate_array ( atomsf_1%a, SIZE ( atomsf_2 ) )
        CALL allocate_array ( atomsf_1%b, SIZE ( atomsf_2 ) )
        atomsf_1%a = atomsf_2%a
        atomsf_1%b = atomsf_2%b
!       Anomalous part:
        atomsf_1%CU  = atomsf_2%CU
        atomsf_1%MO =  atomsf_2%MO
    END SUBROUTINE copy_atomsf

    FUNCTION allocated_sf_database ( sf_database_1 )
        LOGICAL                       :: allocated_sf_database
        TYPE(sf_database), INTENT(IN) :: sf_database_1
        allocated_sf_database = ALLOCATED ( sf_database_1%record )
    END FUNCTION allocated_sf_database
 
    FUNCTION size_of_sf_database ( sf_database_1 )
        INTEGER                       :: size_of_sf_database
        TYPE(sf_database), INTENT(IN) :: sf_database_1
        size_of_sf_database = SIZE ( sf_database_1%record )
    END FUNCTION size_of_sf_database

    SUBROUTINE allocate_sf_database ( sf_database_1, number_of_db_entries )
        TYPE(sf_database), INTENT(INOUT) :: sf_database_1
        INTEGER,           INTENT(IN)    :: number_of_db_entries
!       Local variables:
        INTEGER                          :: istat

        IF ( ALLOCATED ( sf_database_1 ) ) THEN
            CALL deallocate_sf_database ( sf_database_1 )
        ENDIF

        ALLOCATE(sf_database_1%record(number_of_db_entries), STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL die('Failed to allocate SF_DATABASE_1', 'allocate_sf_database')
        ENDIF

    END SUBROUTINE allocate_sf_database
  
    SUBROUTINE deallocate_sf_database ( sf_database_1 )
        TYPE(sf_database), INTENT(INOUT) :: sf_database_1
!       Local variables:
        INTEGER                          :: istat
!       Counters:
        INTEGER                          :: i

!       Deallocate entries in records:
        DO i = 1, SIZE ( sf_database_1 )
            IF ( ALLOCATED ( sf_database_1%record(i) ) ) THEN
                CALL deallocate_array( sf_database_1%record(i)%a )
                CALL deallocate_array( sf_database_1%record(i)%b )
            ENDIF
        ENDDO

!       Deallocate DB itself:
        IF( ALLOCATED ( sf_database_1 ) ) THEN
           DEALLOCATE ( sf_database_1%record, STAT=istat )
        ENDIF
    END SUBROUTINE deallocate_sf_database

    SUBROUTINE copy_sf_database ( sf_database_1, sf_database_2 )
        TYPE(sf_database), INTENT(INOUT) :: sf_database_1
        TYPE(sf_database), INTENT(IN)    :: sf_database_2
!       Counters:
        INTEGER                          :: i

        CALL allocate_sf_database ( sf_database_1, SIZE ( sf_database_2 ) )

        DO i = 1, SIZE ( sf_database_2 )
            sf_database_1%record(i) = sf_database_2%record(i)            
        ENDDO

    END SUBROUTINE copy_sf_database

    SUBROUTINE assign_atom_types ( pdb_2, cards_1 )
!       may work incorrectly if END is added in PDB file (with one atom)
!
!       Purpose:
!       =======
!       Tries to assign sf atom type to all atoms
!
!       MAR 2009. First position in the atom name may be a digit. CORRECTED 28 MAR 2009. BVS.
!
        TYPE(pdb),                      INTENT(INOUT) :: pdb_2
        TYPE(cards),                    INTENT(IN)    :: cards_1
!       Local variables:
        INTEGER                                       :: pos
!       Counters:
        INTEGER                                       :: db_entry
        INTEGER                                       :: iat
        INTEGER                                       :: no_of_entries_in_sf_database
        INTEGER                                       :: nwarn
        INTEGER                                       :: k
        CHARACTER(LEN=32),                       SAVE :: srname = 'assign_atom_types'

        no_of_entries_in_sf_database = SIZE ( pdb_2%atomsf_table%record )

!       Necessary for further checks:
        pdb_2%atom_type = 0
        nwarn = 0

        main_loop:DO iat = 1, SIZE ( pdb_2%atom_name )

            db_loop:DO db_entry =  1, no_of_entries_in_sf_database

                elsymb_business:IF ( ANY ( pdb_2%elsymb /= '  ' ) ) THEN

                    elsymb_not_blank:IF ( pdb_2%elsymb(iat) /= '  ' ) THEN

!                       Try to find exact match, e.g. Cf+3 or He, 
!                       note that 2-gaussian approximation may not be available for such an atom:
                        IF (( TRIM ( pdb_2%elsymb(iat) ) // TRIM ( pdb_2%charge(iat) )) == pdb_2%atomsf_table%record(db_entry)%id & 
                            .AND. pdb_2%atomsf_table%record(db_entry)%id(2:2) /= '*' ) THEN

                            pdb_2%atom_type(iat) = db_entry                      

!                           Atom type has been assigned, assign atom type for next atom:
                            CYCLE main_loop

                        ENDIF

!                       Take care of the most popular atoms:
                        IF ( TRIM ( pdb_2%elsymb(iat) ) // '*' == pdb_2%atomsf_table%record(db_entry)%id(1:2) .AND. &
                             pdb_2%atomsf_table%record(db_entry)%ngauss == cards_1%ngauss ) THEN

                             pdb_2%atom_type(iat) = db_entry
                        
!                            Atom type has been assigned, assign atom type for next atom:
                            CYCLE main_loop
                        ENDIF

                    ENDIF elsymb_not_blank

                ENDIF elsymb_business

            ENDDO db_loop

!           ==== Switch to use of atom names in case elsymb and charge records are absent ====
            nwarn = nwarn + 1
            IF ( nwarn <= 10 ) THEN
                WRITE(*,"('ASSIGN_ATOM_TYPES> ', 'atom name=', A, ' iat=', I8, ' residue=', A3,1X,I4)") &
                          pdb_2%atom_name(iat), iat, ' residue=', pdb_2%residue_name(iat), pdb_2%residue_number(iat)
                CALL warn('Switching to the use of atom names to identify scattering type... Check your PDB file.', &
                          srname)
            ENDIF

!           Try to find an exact match with any number of gaussians for heavy metals, e.g. Cf+3:
            DO db_entry = 1, no_of_entries_in_sf_database

                IF ( pdb_2%atom_name(iat) == pdb_2%atomsf_table%record(db_entry)%id .AND.&
                    pdb_2%number_of_atoms_in_this_residue(iat) == 1 .AND.     &
                    pdb_2%atomsf_table%record(db_entry)%id(2:2) /= '*' ) THEN   ! this is probably unnecessary since, e.g.`C*'
                    pdb_2%atom_type(iat) = db_entry                             ! atom name is unlikely
                    CYCLE main_loop
                ENDIF
            ENDDO

!           If this is heavy metal with long name it's done already:
            IF ( cards_1%ngauss == 2 ) THEN

!               Try to find inexact match for 2-gaussian (H,C,N,O,P,S) atoms:
                DO db_entry = 1, no_of_entries_in_sf_database
                    pos = 1
                    IF ( ch_is_digit (pdb_2%atom_name(iat)(1:1) ) ) pos = 2 
                    IF ( pdb_2%atom_name(iat)(pos:pos)//'*' == pdb_2%atomsf_table%record(db_entry)%id(1:2) .AND. &
                        pdb_2%atomsf_table%record(db_entry)%ngauss == cards_1%ngauss .AND.                 &
                        INDEX ( TWO_GAUSSIAN_ATOMS, pdb_2%atom_name(iat)(pos:pos) ) > 0 ) THEN

                        IF ( debug >= 20 ) THEN
                            WRITE(*,*) pdb_2%atom_name(iat),' ', INDEX ( TWO_GAUSSIAN_ATOMS, pdb_2%atom_name(iat)(1:1)  )
                        ENDIF

                        pdb_2%atom_type(iat) = db_entry
                        CYCLE main_loop
                    ENDIF
                ENDDO            
            ENDIF

!           Try to find inexact match for atoms with any number of gaussians (like FE):
            DO k = 2, 0, -1
                DO db_entry = 1, no_of_entries_in_sf_database
                    pos = 1
                    IF ( ch_is_digit (pdb_2%atom_name(iat)(1:1) ) ) pos = 2

!                   '*' requires 1:2+k index: 
                    IF ( pdb_2%atom_name(iat)(pos:pos+k)//'*' == pdb_2%atomsf_table%record(db_entry)%id(1:2+k) .OR. &
                         pdb_2%atom_name(iat)(pos:pos+k)      == pdb_2%atomsf_table%record(db_entry)%id(1:1+k) ) THEN

                        pdb_2%atom_type(iat) = db_entry
                        CYCLE main_loop
                    ENDIF
                ENDDO
            ENDDO
        ENDDO main_loop

        DO iat = 1, SIZE ( pdb_2%atom_name )
            IF ( debug > 25 ) THEN
                IF ( pdb_2%atom_type(iat) /= 0 ) THEN
                    WRITE(*,"(' ASSIGN_ATOM_TYPES> ', 'atom= ',A, ' assigned type=', I4, ' ngauss=',I2)")&
                               pdb_2%atom_name(iat), pdb_2%atom_type(iat), &
                               pdb_2%atomsf_table%record(pdb_2%atom_type(iat))%ngauss
                ELSE
                    WRITE(*,"(' ASSIGN_ATOM_TYPES> ', 'atom= ',A, ' assigned type=', I4)")&
                               pdb_2%atom_name(iat), pdb_2%atom_type(iat)
                ENDIF
            ENDIF

            IF ( pdb_2%atom_type(iat) == 0 ) THEN

!               Total failure:
                CALL messag('We have only the following atom names in the database : ',&
                            srname)

                WRITE(*,'(16A5)') pdb_2%atomsf_table%record%id
                WRITE(*,"(' ASSIGN_ATOM_TYPES> ', 'Failed to assign ATOM TYPE to atom: ',A, ' iat=', I8)") &
                          pdb_2%atom_name(iat), iat

!               This has to be removed:
!               WRITE(*,*) ' number of atoms in this residue=', number_of_atoms_in_this_residue(iat)
                CALL die('User error. Please rename you atoms accordingly...', srname)

            ENDIF 

        ENDDO

    END SUBROUTINE assign_atom_types

    SUBROUTINE read_sf_database ( sf_data, cards_1, x_scatt_file)
!
!
!        Date:            Programmer:           History of changes:
!        ====             ==========            ==================
!        Sep 2003         B.Strokopytov         Original code
!        Oct 2007         B.Strokopytov         Substanially re-worked
!                                               a) atomsf type introduced
!                                               b) sf_database type added
!                                               c) The only atoms for which we can have
!                                                  2-gaussian approximation is now: "H","C","N" & "O"
!                                                  The rest suspects like "S" or "P" will have always 5-gaussian
!                                                  approximation plus heavy metals and etc. as usual
!                                                  Those elements occur very infrequently in PDB files to make any
!                                                  significant speed improvements
!                                               d) "2-gaussians" are hard-coded into the program at the moment
!                                                  but in principle can be done to be changed by user from input cards
!
!
        TYPE(sf_database), INTENT(INOUT)                :: sf_data
        TYPE(cards),       INTENT(IN)                   :: cards_1
        CHARACTER(len=*),  INTENT(IN),      OPTIONAL    :: x_scatt_file
!       Local params to read a line:
        CHARACTER(LEN=flen)                             :: scatt_file
        CHARACTER(LEN=llen)                             :: line
        CHARACTER(LEN=flen)                             :: CCP4_ATOMSF_file 
        INTEGER                                         :: nwords
        CHARACTER(LEN=wlen), DIMENSION(:), ALLOCATABLE  :: words
        INTEGER                                         :: ngauss
!       Unit for read:
        INTEGER                                         :: iscatt
!       Counters:
        INTEGER                                         :: i
        INTEGER                                         :: db_entry
!       Test file open/close:
        INTEGER                                         :: istat
        
!       try to open a file:
        IF ( PRESENT ( x_scatt_file ) ) THEN
            scatt_file = x_scatt_file
        ELSE
            CALL ugtenv('ATOMSF', CCP4_ATOMSF_file)
            scatt_file = TRIM ( CCP4_ATOMSF_file )
        ENDIF

!       No need to invent numbers with get_next_io_unit:        
        iscatt = get_next_io_unit()

        OPEN ( UNIT=iscatt, file=scatt_file, ACCESS='SEQUENTIAL', &
               ACTION='READ', FORM='FORMATTED', STATUS='OLD', IOSTAT=istat)

        IF ( istat /= 0 ) CALL die ( 'Problems opening sf_database_file : '//scatt_file, &
                                     'read_sf_database')

!       Common macromolecular atoms:
        CALL messag('Common atoms marked with a ''*'' are: '//common_atoms, &
                    'read_sf_database' )

!       Initialize counter:
        db_entry = 0
   
!       Loop through all records in the ATOMSF database:
        DO WHILE (.TRUE.)
            READ( UNIT = iscatt, FMT='(A)', IOSTAT = istat ) line

!           uppercase here to avoid problems... MAR 2009 BVS:
            CALL ucase (line)

!           Check whether everything has been read:
            IF ( istat /= 0 ) EXIT

            ad: IF ( line(1:2) /= 'AD' ) THEN
                 CALL mysplit ( line, words )
                 nwords = SIZE(words)

!                Count number of entries in the database:
                 IF ( line(1:1) /= ' ' .AND. SIZE ( words ) <= 2 ) THEN
                       db_entry = db_entry + 1

!                      skip 4 records:
                       DO i = 1, 4
                           READ ( iscatt, '(A)', IOSTAT = istat ) line
                           IF ( istat /= 0 ) EXIT
                       ENDDO
                 ELSE
                     WRITE(*,*) ' nwords=',nwords
                     WRITE(*,'(A)') line
                     CALL die ( 'Inconsistent lines in ATOMSF database :'//line, 'read_sf_database' )
                 ENDIF 

            ELSE ad
                IF ( debug > 2 ) THEN
                    WRITE(*,"(' READ_CCP4_SF_DATABASE> ', A)") TRIM ( line(1:90) )
                ENDIF
            ENDIF ad
        ENDDO

        WRITE(*,"(' READ_SF_DATABASE> ', 'number of atom entries in ATOMSF database= ', A)")&
        TRIM ( int_to_c ( db_entry ) )

!       Rewind to read again:
        REWIND ( UNIT=iscatt, IOSTAT=istat )
        IF ( istat /= 0 ) THEN
            CALL die ('Failed to rewind unit ISCATT','read_sf_database')
        ENDIF

!       Allocate DB itself (not records):
        IF ( db_entry > 0 ) THEN
            CALL messag ('Allocating SF_DATA records...', 'read_sf_database' )
            CALL allocate_sf_database ( sf_data, db_entry )
        ELSE
            CALL die ('Errors when reading sf_database from file:'//scatt_file, &
                      'read_sf_database')
        ENDIF

!       --------------------
!       Start all over again:
!       --------------------

        db_entry = 0
        DO WHILE( .TRUE. )
            READ ( iscatt, '(A)', IOSTAT = istat) line 
            CALL ucase(line)

!           Check whether everything has been read:
            IF ( istat /= 0 ) THEN
                WRITE(*,"(' READ_SF_DATABASE> ', 'status on exit=', I4)") istat
                EXIT
            ENDIF

            IF ( line(1:2) /= 'AD' ) THEN
                CALL mysplit ( line, words )
                nwords = SIZE ( words )
                IF ( nwords <= 2 .AND. line(1:1) /= ' ' ) THEN
                    
                    db_entry = db_entry + 1

!                   Atom ID (character data):
                    sf_data%record(db_entry)%id = words(1)(1:4)

!                   Number of gaussians to read:
                    ngauss = 5
                    IF ( SIZE ( words ) == 2 ) THEN
                        READ(words(2),*) ngauss
                    ENDIF


!                   Mark with a star all potential macromolecule atoms having name length= 1:
                    IF ( LEN_TRIM ( sf_data%record(db_entry)%id) == 1 ) THEN
                        IF ( INDEX ( common_atoms, TRIM ( sf_data%record(db_entry)%id ) ) > 0 ) THEN
                            sf_data%record(db_entry)%id(2:2) = '*'
                        ENDIF
                    ENDIF

                    IF ( cards_1%debug > 25 ) THEN
                        WRITE(*,"(' READ_SF_DATABASE> ','Reading atom_name=',A,' ngauss=',I2)")&
                        sf_data%record(db_entry)%id, ngauss
                    ELSE IF ( cards_1%debug > 2 ) THEN
                        IF ( INDEX ( sf_data%record(db_entry)%id, '*' ) > 0 ) THEN
                            WRITE(*,"(' READ_SF_DATABASE> ','Reading atom_name=',A,' ngauss=',I2)")&
                            sf_data%record(db_entry)%id, ngauss
                        ENDIF
                    ENDIF

!                   Allocate arrays according to ngauss:
                    sf_data%record(db_entry)%ngauss = ngauss
                    CALL allocate_atomsf(sf_data%record(db_entry), ngauss )

                    IF ( SIZE ( sf_data%record(db_entry) ) /= ngauss ) THEN
                        WRITE(*,*) SIZE ( sf_data%record(db_entry) ), ngauss
                        CALL die ('Programming error. Strange allocation...', 'read_sf_database' )
                    ENDIF

!                   Read next line:
                    READ ( iscatt, '(A)', IOSTAT = istat) line
                    CALL mysplit ( line, words )
                    nwords = SIZE ( words )

!                   Check:
                    IF ( nwords /= 3 ) THEN
                         WRITE(*,"(' READ_SF_DATABASE> ', 'atom_name= ',A)")&
                         sf_data%record(db_entry)%id
                         CALL die (' Syntax error in sf_database :'//line,'read_sf_database')
                    ELSE
                        READ (words(1),*)  sf_data%record(db_entry)%atomic_charge
                        READ (words(2),*)  sf_data%record(db_entry)%number_of_electrons

!                       Read constant "c" in case of 5 gaussians only:
                        IF ( ngauss == 5 ) THEN
                            READ (words(3),*)  sf_data%record(db_entry)%a(ngauss)
                            sf_data%record(db_entry)%b(ngauss) = 0.0_wp
                        ENDIF
                    ENDIF

                    READ ( iscatt, '(A)', IOSTAT = istat) line 
                    CALL mysplit ( line, words)
                    nwords = SIZE ( words )

!                   Always check what we are reading:
                    IF ( nwords /= 4 ) THEN
                        CALL die (' Syntax error in sf_databse :'//line,'read_sf_database')
                    ELSE
!                       Read "A" line with care otherwise it won't fit into allocated arrays: 
                        READ(line,*) sf_data%record(db_entry)%a(1:MIN(4,ngauss))
                    ENDIF

                    READ ( iscatt, '(A)', IOSTAT = istat) line
                    CALL mysplit ( line, words )
                    nwords = SIZE(words)
                    IF ( nwords /= 4 ) THEN
                        CALL die (' Syntax error in sf_databse :'//line,'read_sf_database')
                    ELSE
!                       Read "B" line with care otherwise it won't fit into allocated arrays:
                        READ(line,*) sf_data%record(db_entry)%b(1:MIN(4,ngauss))
                    ENDIF

!                   Read anomalous data:
                    READ ( iscatt, '(A)', IOSTAT=istat) line
                    CALL mysplit(line, words)
                    nwords = SIZE(words)
                    IF ( nwords /=  4 ) THEN
                        WRITE(*, '(A)') line
                        CALL die ( 'Syntax error in sf_databse :'//line, 'read_sf_database' )
                    ELSE
                        READ(line,*) sf_data%record(db_entry)%CU(1:2), sf_data%record(db_entry)%MO(1:2)
                    ENDIF

                ELSE
!                   In no circumstances we arrive here:
                    WRITE(*,"(' Ignoring: ', A, 1X, A)") words(1:2)
                    CALL die('Syntax error when reading SF_DATA:'//line, 'read_sf_database' )
                ENDIF
            ENDIF    
        ENDDO

!       Database in memory now:
        CLOSE ( UNIT = iscatt, IOSTAT = istat )
        IF ( istat /= 0 ) THEN
             CALL warn ('Failed to close sf database','read_sf_database')
        ELSE
             WRITE(*,"(' READ_SF_DATABASE> ','sf database has been successfully closed.')")
        ENDIF

!       Free memory:
        IF ( ALLOCATED(words) ) CALL deallocate_array ( words )

    END SUBROUTINE read_sf_database

END MODULE basic_pdb_manip
