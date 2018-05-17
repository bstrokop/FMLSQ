MODULE block_diagonal_util
USE anymatrix_operations
USE constants
USE basic_pdb_manip
USE fail
USE sparse_basic
USE util
IMPLICIT NONE
CONTAINS
    SUBROUTINE residue_block_diagonal_pairs(pair, pdb_2, residues_in_block, nout)
!
!       Purpose:
!       =======
!       Prepares list of atomic pairs for each residue.
!       This is necessary to create block diagonal matrix
!       which does not suffer from major drawbacks of arbitrary 
!       sparse matrix (indefineteness).
!
!       Block diagonal matrix created this way is guaranteed
!       to have positive eigenvalues provided that number
!       of observations exceeds number of refinable parameters.
!
!       Date:               Programmer:                 History of changes:
!       ====                ==========                  ==================
!       Mar 2008            B.Strokopytov               Original code
!
        INTEGER,              DIMENSION(:,:),   ALLOCATABLE, INTENT(INOUT) :: pair
        TYPE(pdb),                                           INTENT(INOUT) :: pdb_2
        INTEGER,                                             INTENT(IN)    :: residues_in_block
        INTEGER,                                             INTENT(IN)    :: nout
        REAL(KIND=sp)                                                      :: average_number_of_contacts
!       Local variables:
        INTEGER                                                            :: modcon
        INTEGER                                                            :: nres
        INTEGER                                                            :: natom
        INTEGER                                                            :: np_block
!       Counters:
        INTEGER                                                            :: i
        INTEGER                                                            :: ires
        INTEGER                                                            :: j
        INTEGER(KIND=eb)                                                   :: number_of_saved_pairs
        INTEGER                                                            :: run
        INTEGER                                                            :: nnz
        INTEGER                                                            :: number_of_diagonal_blocks
!       SUB name:
        CHARACTER(LEN=32)                                                  :: srname='residue_block_diagonal_pairs'

!       Checkz:
        nres = SIZE ( pdb_2%residue_start ) - 1
        IF ( nres < 1 ) THEN
            WRITE(*,*) nres
            CALL die('Nothing to process.', srname)
        ENDIF

        outer_loop:DO run = 0, 1

            number_of_saved_pairs = 0
            number_of_diagonal_blocks = 0
            nnz = 0
            IF ( run == 0 ) THEN
                WRITE(*,"(' RESIDUE_BLOCK_DIAGONAL_PAIRS> ', 'Going to process ', A, ' residues (including solvent)')") &
                TRIM ( int_to_c ( nres ) )
            ENDIF

!           Create blocks. Obviously enough number of residues in block cannot excced NRES:
            res_loop:DO ires = 1, nres, MIN ( residues_in_block, nres )

!               Number of atoms in this residue block, ires+residue_block_size cannot exceed nrez+1:
                natom = pdb_2%residue_start(MIN(ires+residues_in_block,nres+1)) - pdb_2%residue_start(ires)

                IF ( natom > 0 ) THEN

                    number_of_diagonal_blocks = number_of_diagonal_blocks + 1

                    IF ( run == 1 ) THEN

!                       Counter number of refinable parameters in this block:
                        np_block = COUNT ( pdb_2%refinable(1:11,pdb_2%residue_start(ires)+1:pdb_2%residue_start(ires)+natom) )

!                       This will help us to identify number of blocks and their positions in the normal matrix:
                        pdb_2%pblock(nnz+1:nnz+np_block) = number_of_diagonal_blocks
                        nnz = nnz + np_block

                    ENDIF

                ELSE IF ( natom == 0 ) THEN

!                   Just warn:
                    WRITE(*,*) pdb_2%residue_start(ires+1), pdb_2%residue_start(ires)
                    WRITE(*,*) ' natom=', pdb_2%residue_start(ires+1), pdb_2%residue_start(ires)
                    CALL die('Number of atoms in the block equal to zero... Skipping the block...', srname)

                ENDIF

!               Need to run this cycles to create a list of (saved) pairs:
                DO i = 1, natom 
                    DO j = i, natom 
                        
                        number_of_saved_pairs = number_of_saved_pairs + 1 
!                       PAIR array has been allocated:
                        IF ( run == 1 ) THEN
                            pair(1,number_of_saved_pairs) = pdb_2%residue_start(ires) + i
                            pair(2,number_of_saved_pairs) = pdb_2%residue_start(ires) + j
!                           Monitor atom pairs:
                            IF ( nout > 0 ) THEN
                                IF ( MOD ( number_of_saved_pairs, modcon ) == 0 )  THEN
                                    WRITE(*,"(' RESIDUE_BLOCK_DIAGONAL_PAIRS> ', 'entry number=', I9, ' atomic pair=', 2I7)") &
                                    number_of_saved_pairs, pair(1:2,number_of_saved_pairs)
                                ENDIF
                            ENDIF

                        ENDIF

                    ENDDO
                ENDDO

!               Allocate PAIR array:
                IF ( run == 0 ) THEN
                    CALL allocate_array ( pair, 2_eb, number_of_saved_pairs )
                    CALL allocate_array ( pdb_2%pblock, COUNT ( pdb_2%refinable ) )
                    pdb_2%pblock = 0
                    modcon = MAX ( number_of_saved_pairs / 10, 1)
                ENDIF
            ENDDO res_loop
        ENDDO outer_loop

!       Simple check:
        IF ( number_of_saved_pairs > 0 ) THEN
            CALL messag('Finished pair list preparation...', srname)
        ELSE
            WRITE(*,*) number_of_saved_pairs
            CALL messag('Failed to create list of atom pairs...', srname)
        ENDIF

    END SUBROUTINE residue_block_diagonal_pairs

    SUBROUTINE atom_block_diagonal_pairs(pair, pdb_2, atoms_in_block, nout)
!
!       Purpose:
!       =======
!       Prepares list of atomic pairs for simplest block diagonal
!       approximation of sparse precoditioner.
!
!       Date:               Programmer:                 History of changes:
!       ====                ==========                  ==================
!       Mar 2008            B.Strokopytov               Original code
!
        INTEGER,              DIMENSION(:,:),   ALLOCATABLE, INTENT(INOUT) :: pair
        TYPE(pdb),                                           INTENT(INOUT) :: pdb_2
        INTEGER,                                             INTENT(IN)    :: atoms_in_block     ! step                                         
        INTEGER,                                             INTENT(IN)    :: nout
!       Local variables:
        INTEGER                                                            :: modcon
        INTEGER                                                            :: natom
        INTEGER                                                            :: np_block
!       Counters:
        INTEGER                                                            :: i
        INTEGER                                                            :: number_of_diagonal_blocks
        INTEGER                                                            :: nnz
!       SUB name:
        CHARACTER(LEN=32)                                                  :: srname='atom_block_diagonal_pairs'

!       Checkz:
        natom = SIZE ( pdb_2 )
        IF ( natom < 1 ) THEN
            WRITE(*,*) natom
            CALL die('Nothing to process.', srname)
        ENDIF

        modcon = MAX ( natom / 10, 1 )
        CALL allocate_array ( pair, 2, natom )

!       Allocate and intialize pointers for sparse matrix:
        CALL allocate_array ( pdb_2%pblock, COUNT ( pdb_2%refinable ) )
        pdb_2%pblock = 0
        nnz = 0
        number_of_diagonal_blocks = 0

        DO i = 1, natom, atoms_in_block
            pair(1:2,i) = i
            np_block = COUNT ( pdb_2%refinable(1:11,i:MIN(i+atoms_in_block-1,natom)) )
            number_of_diagonal_blocks = number_of_diagonal_blocks + 1
            pdb_2%pblock(nnz+1:nnz+np_block) = number_of_diagonal_blocks
            nnz = nnz + np_block
            IF ( nout > 0 ) THEN
                IF ( MOD ( i, modcon ) == 0 ) THEN
                    WRITE(*,"(' ATOM_BLOCK_DIAGONAL_PAIRS> ', ' iat=', I6, ' jat=', I6)") &
                    pair(1:2,i)
                ENDIF
            ENDIF
        ENDDO
        IF ( debug > 100 ) THEN
            WRITE(*,"(30I6)") pdb_2%pblock
        ENDIF
        CALL messag('Finished pair list preparation...', srname)

    END SUBROUTINE atom_block_diagonal_pairs

    FUNCTION calc_nnz ( pdb_2, atom_pair ) RESULT ( nnz )
!       
!       Purpose:
!       =======
!       Predicts number of non-zero elements in sparse normal matrix.
!       Necessary for proper allocation of sparse matrices.
!
!       FIXME: Probably we need a function, not a subroutine.
!
        INTEGER                                           :: nnz
        TYPE(pdb),                            INTENT(IN)  :: pdb_2
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)  :: atom_pair
!       Local variables:
        INTEGER                                           :: block_size 
!       Pointers:
        INTEGER                                           :: iat
        INTEGER                                           :: jat
!       Counters:
        INTEGER                                           :: i
        INTEGER                                           :: matrix_elements
!       Test:
        INTEGER                                           :: iatold
        INTEGER                                           :: duplicates 
        CHARACTER(LEN=1)                                  :: warning
!       Name it:
        CHARACTER(LEN=32)                                 :: srname='calc_nnz'

!       Start with basic checks:
        IF ( COUNT ( pdb_2%refinable ) == 0 ) THEN
            CALL die('Programming error. Nothing to refine.', 'calc_nnz')
        ENDIF

!       Initialize counter of non-zero elements:
        matrix_elements = 0
        iatold     = -1
        duplicates = 0

        atomic_pairs:DO i = 1, SIZE ( atom_pair, DIM=2 )

            iat = atom_pair(1,i)
            jat = atom_pair(2,i)

            IF ( iat > jat ) THEN
                WRITE(*,*) ' iat=', iat, ' jat=', jat
                CALL die('Programming error. Incorrectly prepared list. IAT > JAT.', srname)
            ENDIF

!           Count number of matrix elements:
            block_size = COUNT ( pdb_2%refinable(1:MAX_BLOCK_SIZE,iat) ) &
                       * COUNT ( pdb_2%refinable(1:MAX_BLOCK_SIZE,jat) )

!           Add all atomic pair block sizes up:
            IF ( iat == jat ) THEN
                matrix_elements = matrix_elements + block_size
            ELSE
!               Scale up due to matrix symmetry:
                matrix_elements = matrix_elements + 2 * block_size
            ENDIF

            IF ( debug > 100 ) THEN
                WRITE(*,*)  ' iat=', iat, ' jat= ', jat, ' block_size=', block_size           
                WRITE(*,*)  ' curr num matrix elem=', matrix_elements
            ENDIF

!           Simple test for main diagonal duplicate entries:
            IF ( iat == jat ) THEN

                warning = ' '
                IF ( iat == iatold ) THEN 
                    warning = '*'
                    WRITE(*,"(' CALC_NNZ> ', 'iat/jat=', 2I6, ' entry number=', I10, A)") &
                    iat, jat, matrix_elements, warning
                    duplicates = duplicates + 1
                ENDIF

                iatold = iat                            
                
            ENDIF

        ENDDO atomic_pairs

        IF ( duplicates >  0 ) THEN
            WRITE(*,"(' CALC_NNZ> ', 'Oops... ', A, ' diagonal duplicates found...')") &
            TRIM ( int_to_c ( duplicates ) )
            CALL die('Programming error. Duplicates found.', srname)
        ENDIF

        nnz = matrix_elements

    END FUNCTION calc_nnz

    SUBROUTINE sparse_matrix_atom_block_pointers ( pdb_2, atom_pair, atom_block )
        TYPE(pdb),                            INTENT(IN)    :: pdb_2
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)    :: atom_pair
        INTEGER,  DIMENSION(:),  ALLOCATABLE, INTENT(INOUT) :: atom_block
!       Pointers:
        INTEGER                                             :: iat
        INTEGER                                             :: jat
!       Counters:
        INTEGER                                             :: ip
        INTEGER                                             :: n
        INTEGER                                             :: nnz

!       Need testing:
        CALL allocate_array ( atom_block, SIZE ( atom_pair, DIM=2 ) )
        nnz = 0
        atom_block(1) = 0
        DO ip = 1, SIZE ( atom_pair, DIM=2 ) - 1
            iat = atom_pair(1,ip)
            jat = atom_pair(2,ip)
            IF ( jat > iat ) THEN
                nnz = nnz + COUNT ( pdb_2%refinable(1:max_block_size,iat) ) &
                          * COUNT ( pdb_2%refinable(1:max_block_size,jat) )
            ELSE IF ( jat == iat ) THEN
                n =  COUNT ( pdb_2%refinable(1:max_block_size,iat) )
                nnz = nnz + ( n * (n + 1) ) / 2 
            ENDIF
            atom_block(ip+1) = nnz
        ENDDO

    END SUBROUTINE sparse_matrix_atom_block_pointers

    SUBROUTINE copy_coo_matrix_to_array_of_anymatrices ( H, pdb_2, ANYMAT )
        TYPE(coo_matrix),                            INTENT(IN)    :: H
        TYPE(pdb),                                   INTENT(IN)    :: pdb_2
        TYPE(anymatrix),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ANYMAT
!       Local variables:
        INTEGER                                                    :: number_of_diagonal_blocks
        INTEGER,          DIMENSION(:), ALLOCATABLE                :: block_sizes
        INTEGER,          DIMENSION(:), ALLOCATABLE                :: pblock_sizes
        INTEGER                                                    :: diagonal_block_number
        INTEGER                                                    :: shift
        INTEGER                                                    :: ia
        INTEGER                                                    :: ja
        INTEGER                                                    :: n
        INTEGER                                                    :: m 
        REAL(KIND=wp)                                              :: val
!       Counters:
        INTEGER                                                    :: i
        INTEGER                                                    :: l
        CHARACTER(LEN=64),                                 SAVE    :: srname='copy_coo_matrix_to_array_of_anymatrices'
  
!       Calculate nessary block sizes and pointers:
        CALL calculate_diagonal_block_sizes(block_sizes, pblock_sizes, number_of_diagonal_blocks, pdb_2)

        IF ( .NOT. ALLOCATED ( ANYMAT ) ) THEN
            CALL allocate_array ( ANYMAT, number_of_diagonal_blocks )
            DO i = 1, number_of_diagonal_blocks
                CALL allocate_anymatrix ( ANYMAT(i), block_sizes(i), block_sizes(i) )
                ANYMAT(i) = 0.0_wp
            ENDDO
        ENDIF

!       Check parameters of sparse coordinate matrix:
        IF ( H%nnz == 0 ) THEN
            WRITE(*,*) ' Number of non-zero elements=', H%nnz
            CALL die('What kind of coordinate matrix is this?', srname)
        ELSE IF ( H%nnz /= SIZE ( H%val ) ) THEN
            WRITE(*,*) ' H%nnz=', H%nnz, ' size=', SIZE ( H%val )
            CALL die('What kind of sparse coordinate matrix is this?', srname)
        ENDIF

!       Loop over all available elements:
        DO l= 1, H%nnz

!           Zero indices can happen for unexpanded matrices:
            IF ( H%ia(l) == 0 .OR. H%ja(l) == 0 ) CYCLE

            diagonal_block_number = pdb_2%pblock(H%ja(l))
            shift = pblock_sizes(diagonal_block_number)

!           For block diagonal matrices both IA and JA should always give the same block number due to symmetry:
            IF ( pdb_2%pblock(H%ia(l)) /= diagonal_block_number ) THEN
                WRITE(*,*) ' JA=', H%ja(l), ' IA=', H%ia(l)
                WRITE(*,*) ' block No.=', diagonal_block_number, ' block size=', block_sizes(diagonal_block_number)
                WRITE(*,*) ' shift=', shift
                CALL die('Programming error. This does not look like block diagonal matrix.', srname)
            ENDIF


!           Confirm size of current matrix:
            n = SIZE ( ANYMAT(diagonal_block_number), DIM=1 )
            m = SIZE ( ANYMAT(diagonal_block_number), DIM=2 )

            ia = H%ia(l) - shift
            ja = H%ja(l) - shift

            IF ( ia <= 0 .OR. ja <= 0 ) THEN
                WRITE(*,*) ia, ja
                WRITE(*,*) H%ia(l), H%ja(l), shift
                CALL die('Programming error. Unreasonable indices detected.', srname)
            ENDIF


!           Make checkz:
            IF ( ia > n .OR. ja > m ) THEN
                WRITE(*,*) ' ia=', ia, ' ja=', ja, n, m 
                CALL die('Programming error. Matrix indices exceed allocated size of matrix.', srname)
            ENDIF

            CALL set_anymatrix_elem ( ANYMAT(diagonal_block_number), H%val(l), ia, ja ) 

        ENDDO
 
        my_diagonal_blocks:DO i = 1, number_of_diagonal_blocks
            shift = pblock_sizes(i)

!           Debugging print:
            IF ( debug > 100 ) THEN
                WRITE(*,*) ' block No.=', i, ' shift=', shift
                CALL print_anymatrix(ANYMAT(i))
            ENDIF

            n = SIZE ( ANYMAT(i), DIM=1 )
            m = SIZE ( ANYMAT(i), DIM=2 )

!           Nothing to do if matrix has zero size:
            IF ( n == 0 ) THEN
                CYCLE
            ENDIF

!           Everything should be in the upper triangular part: 
            DO ia = 1, n 
                DO ja = ia, m
                    CALL get_anymatrix_elem(val, ANYMAT(i), ia, ja)
                    IF ( val == 0.0_wp ) THEN
                        WRITE(*,*) ' anymat_sizes=', n, m 
                        WRITE(*,*) ' ia=', ia, ' ja=', ja, ' val=', val
                        CALL warn('Matrix element is equal exactly to zero.', srname)
                    ENDIF

!                   Check for non-zeroes lower part:
!                    IF ( ia /= ja ) THEN
!                        CALL get_anymatrix_elem(val, ANYMAT(i), ja, ia)
!                        IF ( val == 0.0_wp ) THEN
!                            WRITE(*,*) ' anymat_sizes=', n, m
!                            WRITE(*,*) ' ia=', ja, ' ja=', ia, ' val=', val
!                            CALL warn('Matrix element is equal exactly to zero.', srname)
!                        ENDIF
!                    ENDIF

                ENDDO
            ENDDO

!           Won't hurt to check main diagonal again:
            DO ia = 1, n
                IF ( ANYMAT(i)%A(ia,ia) <= 0.0_wp ) THEN
                    WRITE(*,*) ' diagonal block=',i, ' shift=', shift
                    WRITE(*,*) ' block diagonal matrix element', ia,ia, 'is zero or negative.'
                    CALL die('Negative matrix element on the diagonal...', srname)
                ENDIF
            ENDDO

            IF ( i == 1 .AND. debug > 1000 ) THEN
                WRITE(*,*) ' First block:'
                CALL print_anymatrix ( ANYMAT(i) )
            ENDIF

!           Calculate eigenvalues and eigenvectors:
            CALL store_anymatrix_as_eigenvectors(ANYMAT(i))

!           Check for zero or negative eigenvalues:
            m = 0
            DO ia = 1, n
                IF ( ANYMAT(i)%e(ia) <= 0.0_wp ) THEN
                    m = m + 1
                    WRITE(*,"(' COPY_COO_MATRIX_TO_ARRAY_OF_ANYMATRICES> ', 'shift=', I8)") shift

!                   It's important to find contributions largest in magnitude (therefore ABS function used):
                    WRITE(*,"(' COPY_COO_MATRIX_TO_ARRAY_OF_ANYMATRICES> ', 'Lambda No.=', I5, ' value=', ES12.5,&
                              & ' offending param= ', A,' its value=',ES9.2 )") ia, ANYMAT(i)%e(ia),    &
                              TRIM ( int_to_c ( shift+MAXLOC ( ABS ( ANYMAT(i)%a(1:n,ia) ), DIM=1) ) ), &
                              MAXVAL ( ABS ( ANYMAT(i)%a(1:n,ia) ) )
                ENDIF
            ENDDO

!           Print excessive matrix block condition number:
            IF ( ANYMAT(i)%e(n) / ANYMAT(i)%e(1) > 10_wp ** 6 ) THEN
                WRITE(*,"(' COPY_COO_MATRIX_TO_ARRAY_OF_ANYMATRICES> ', 'Block=',I5,&
                         &' Matrix condition number for this block= ',ES9.2)") i, ANYMAT(i)%e(n) / ANYMAT(i)%e(1)
            ENDIF

            IF ( m /= 0 ) THEN
                WRITE(*,"(' COPY_COO_MATRIX_TO_ARRAY_OF_ANYMATRICES> ', 'Number of negative lamdbas detected=', I4)") m
                WRITE(*,*) ' Eigenvalues:'
                WRITE(*,"(2000ES9.2)") ANYMAT(i)%e
                WRITE(*,*) ' Eigenvectors:'
                DO l = 1, m
                    WRITE(*,"(1X,I8,')',1X,2000ES9.2)") l, ANYMAT(i)%a(1:n,l)
                ENDDO

!               That's an overkill to print all this (subject for removal):
!               CALL print_anymatrix ( ANYMAT(i) )

!               If ATOMS unusually close to each other a lambda close to zero may occur 
!               and even become negative due to errors in calcns:
                IF ( m == 1 .AND. n > 1 ) THEN

!                   Correct offending negative lambda:
                    ANYMAT(i)%e(1) = 0.1_wp * ANYMAT(i)%e(2)
                    CALL messag('Resetting one negative lambda...', srname)
                ELSE
!                   Too many negative lambdas. Need to investigate. Don't continue calcns:
                    CALL die ( 'Too many negative lambdas detected.', srname )
                ENDIF
            ENDIF

        ENDDO my_diagonal_blocks

!       Free memory:
        CALL deallocate_array(block_sizes)
        CALL deallocate_array(pblock_sizes)

    END SUBROUTINE copy_coo_matrix_to_array_of_anymatrices

    SUBROUTINE copy_array_of_inverted_anymatrices_to_coo_matrix(H, pdb_2, ANYMAT, xpower)
        TYPE(coo_matrix),                            INTENT(INOUT)        :: H
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        TYPE(anymatrix),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: ANYMAT
        REAL(KIND=wp),                               INTENT(IN), OPTIONAL :: xpower
!       Local variables:
        REAL(KIND=wp)                                                     :: power
        INTEGER                                                           :: number_of_diagonal_blocks
        INTEGER,          DIMENSION(:), ALLOCATABLE                       :: block_sizes
        INTEGER,          DIMENSION(:), ALLOCATABLE                       :: pblock_sizes
        INTEGER                                                           :: diagonal_block_number
        INTEGER                                                           :: shift
        INTEGER                                                           :: ia
        INTEGER                                                           :: ja
        INTEGER                                                           :: n
        INTEGER                                                           :: m
        REAL(KIND=wp)                                                     :: val
        REAL(KIND=wp),    DIMENSION(:,:), ALLOCATABLE                     :: AINV
!       Counters:
        INTEGER                                                           :: i
        INTEGER                                                           :: l
        INTEGER                                                           :: ii
        INTEGER                                                           :: jj
        CHARACTER(LEN=64),                                           SAVE :: srname

        srname = 'copy_array_of_inverted_anymatrices_to_coo_matrix'
        power = -1.0_wp
        IF ( PRESENT ( xpower ) ) THEN
            power = xpower
        ENDIF 

!       Calculate nessary block sizes and pointers:
        CALL calculate_diagonal_block_sizes(block_sizes, pblock_sizes, number_of_diagonal_blocks, pdb_2)
        H%val = 0
        H%ia = 0
        H%ja = 0
        l = 0
        DO i = 1, number_of_diagonal_blocks
            n = SIZE ( ANYMAT(i), DIM=1 )
            m = SIZE ( ANYMAT(i), DIM=2 )
            CALL allocate_array(AINV, n, m)
            shift = pblock_sizes(i)

!           Invert matrices one by one:
            IF ( power == -1.0_wp ) THEN
!               Too many routines FIXME:
                CALL get_inverted_anymatrix(AINV, ANYMAT(i))
            ELSE IF ( power == 0.5_wp ) THEN
                CALL get_sqrt_of_anymatrix(AINV, ANYMAT(i))
            ELSE IF ( power == -0.5_wp ) THEN
                CALL get_sqrt_Of_inverted_anymatrix(AINV, ANYMAT(i))
            ELSE
                WRITE(*,*) ' power=', power
                CALL die('Programming error. Unreasonable (X)POWER option.', srname)
            ENDIF

            IF ( debug > 100 ) THEN
                DO ia = 1, n
                    WRITE(*,"(' COPY_ARRAY_OF_INVERTED_ANYMATRICES_TO_COO_MATRIX> ', 1000ES9.2)") &
                                (AINV(ia,ja),ja=1,m)
                ENDDO
            ENDIF
!           Create new inverted (and sorted) H matrix:
            DO ia = 1, n
                DO ja = 1, n
                    l = l + 1
                    H%ia(l) = ia + shift
                    H%ja(l) = ja + shift
                    H%val(l) = AINV(ia,ja)
                    IF ( ia == ja ) THEN
                        IF ( H%val(l) <= 0.0_wp ) THEN
                            WRITE(*,*) ' diagonal block=', i,  ' shift=', shift
                            WRITE(*,*) ' H%val(l)=', H%val(l), ' ia=', ia, ' ja=', ja
                            IF ( debug > 3 ) THEN
                                IF ( n < 1000 ) THEN
                                    WRITE(*,"(51X,1000I9)") (ii,ii=1,m)
                                    DO ii = 1, n
                                        WRITE(*,"(' COPY_ARRAY_OF_INVERTED_ANYMATRICES_TO_COO_MATRIX> ', 1000ES9.2)") &
                                                 (AINV(ii,jj),jj=1,m)
                                    ENDDO
                                ENDIF
                            ENDIF
                            CALL die('Negative matrix element on the main diagonal.', srname)
   
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            CALL deallocate_array(AINV)
        ENDDO

        IF ( l /= H%nnz ) THEN
            WRITE(*,*) ' l=', l, ' H%nnz=', H%nnz
            CALL die('Programming error. Inconsistent number of elements...', srname)
        ELSE IF (  ANY ( H%val == 0.0_wp ) ) THEN
            CALL warn('Possible programming error. Zero values found in H%val', srname) 
        ELSE IF ( ANY (H%ia == 0 ) ) THEN
            CALL die('Programming error. Zero indices found.', srname)
        ENDIF

!       Free memory:
        CALL deallocate_array(block_sizes)
        CALL deallocate_array(pblock_sizes)
!       No need to keep this array:
        CALL deallocate_array(ANYMAT) 
    END SUBROUTINE copy_array_of_inverted_anymatrices_to_coo_matrix

    SUBROUTINE calculate_diagonal_block_sizes ( block_sizes, pblock_sizes, number_of_diagonal_blocks, pdb_2 )
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: block_sizes
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: pblock_sizes
        INTEGER,                            INTENT(INOUT) :: number_of_diagonal_blocks
        TYPE(pdb),                          INTENT(IN)    :: pdb_2
!       Counters:
        INTEGER                                           :: i
        INTEGER                                           :: l
        CHARACTER(LEN=32),                        SAVE    :: srname = 'calculate_diagonal_block_sizes'
 
!       Last elements of PBLOCK contains the much needed number:
        number_of_diagonal_blocks = pdb_2%pblock( SIZE ( pdb_2%pblock ) )

!       Checkz:
        IF ( number_of_diagonal_blocks ==  0 ) THEN
            WRITE(*,*) number_of_diagonal_blocks
            CALL die('NUMBER OF DIAGONAL BLOCKS equals zero.', srname)
        ENDIF

        CALL allocate_array(block_sizes, number_of_diagonal_blocks)
!        DO i = 1, number_of_diagonal_blocks
!            block_sizes(i) = COUNT ( pdb_2%pblock == i )
!        ENDDO

        block_sizes = 0
        DO i = 1, COUNT ( pdb_2%refinable )
            block_sizes(pdb_2%pblock(i)) = block_sizes(pdb_2%pblock(i)) + 1
        ENDDO

        IF ( debug > 1000 ) THEN
            WRITE(*,*) 'block_sizes:'
            WRITE(*,"(10I5)") block_sizes
        ENDIF

!       Rearrange BLOCK_SIZES into pointers:
        CALL allocate_array( pblock_sizes, SIZE ( block_sizes ) )
        l = 0
        DO i = 1, number_of_diagonal_blocks
            pblock_sizes(i) = l
            l = l + block_sizes(i)
        ENDDO

    END SUBROUTINE calculate_diagonal_block_sizes 


END MODULE block_diagonal_util
