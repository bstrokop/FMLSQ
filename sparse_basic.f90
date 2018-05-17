MODULE sparse_basic
USE sparse
USE util
IMPLICIT NONE
TYPE sparse_matrix
    INTEGER                                  :: ncol
    INTEGER                                  :: nrow
    INTEGER(KIND=eb)                         :: nnz
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ia
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ja
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: val 
END TYPE
TYPE coo_matrix
    INTEGER                                  :: ncol
    INTEGER                                  :: nrow
    INTEGER(KIND=eb)                         :: nnz
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ia
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ja
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: val
END TYPE
TYPE csr_matrix
    INTEGER                                  :: ncol
    INTEGER                                  :: nrow
    INTEGER(KIND=eb)                         :: nnz
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ia
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ja
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: val
END TYPE

TYPE diag_matrix
    INTEGER                                  :: ncol
    INTEGER                                  :: nrow
    INTEGER(KIND=eb)                         :: nnz
!   FIXME: do we need IA, JA for diag matrix?:
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ia
    INTEGER,       DIMENSION(:), ALLOCATABLE :: ja
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: val
END TYPE

TYPE dense_matrix
    INTEGER                                    :: ncol
    INTEGER                                    :: nrow
    INTEGER(KIND=eb)                           :: nnz
    REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: val
END TYPE

PRIVATE
PUBLIC :: coo_matrix
PUBLIC :: csr_matrix
PUBLIC :: ASSIGNMENT( = )
PUBLIC :: OPERATOR ( * )
PUBLIC :: ALLOCATED
PUBLIC :: ALLOCATE_MATRIX
PUBLIC :: DEALLOCATE_MATRIX
PUBLIC :: condition_coo_matrix
PUBLIC :: expand_coo_matrix
PUBLIC :: matvec
INTERFACE ALLOCATED
    MODULE PROCEDURE coo_matrix_allocated
    MODULE PROCEDURE csr_matrix_allocated
    MODULE PROCEDURE diag_matrix_allocated
END INTERFACE

INTERFACE ASSIGNMENT( = )
    MODULE PROCEDURE copy_coo_to_csr_matrix
    MODULE PROCEDURE copy_coo_to_dense_matrix
END INTERFACE

INTERFACE OPERATOR ( * )
    MODULE PROCEDURE matrix_times_vector
END INTERFACE

INTERFACE allocate_matrix
    MODULE PROCEDURE allocate_coo_matrix
    MODULE PROCEDURE allocate_csr_matrix
END INTERFACE


INTERFACE deallocate_matrix
    MODULE PROCEDURE deallocate_coo_matrix
    MODULE PROCEDURE deallocate_csr_matrix
END INTERFACE


CONTAINS
    FUNCTION coo_matrix_allocated ( A )
        LOGICAL                      :: coo_matrix_allocated
        TYPE(coo_matrix), INTENT(IN) :: A
        coo_matrix_allocated = ALLOCATED ( A%val )
    END FUNCTION coo_matrix_allocated

    FUNCTION csr_matrix_allocated ( B ) 
        LOGICAL                      :: csr_matrix_allocated
        TYPE(csr_matrix), INTENT(IN) :: B 
        csr_matrix_allocated = ALLOCATED ( B%val )
    END FUNCTION csr_matrix_allocated

    FUNCTION diag_matrix_allocated ( C )
        LOGICAL                       :: diag_matrix_allocated
        TYPE(diag_matrix), INTENT(IN) :: C
        diag_matrix_allocated = ALLOCATED ( C%val )
    END FUNCTION diag_matrix_allocated

    SUBROUTINE allocate_coo_matrix(A, nrow, ncol, nnz)
        TYPE(coo_matrix), INTENT(INOUT) :: A
        INTEGER,          INTENT(IN)    :: nrow
        INTEGER,          INTENT(IN)    :: ncol 
        INTEGER(KIND=eb), INTENT(IN)    :: nnz
        A%nrow = nrow
        A%ncol = ncol
        A%nnz  = nnz
        CALL allocate_array ( A%val, nnz )
        CALL allocate_array ( A%ia, nnz )
        CALL allocate_array ( A%ja, nnz )
    END SUBROUTINE allocate_coo_matrix

    SUBROUTINE deallocate_coo_matrix(A)
        TYPE(coo_matrix) :: A

        IF ( ALLOCATED ( A%val ) ) CALL deallocate_array ( A%val )
        IF ( ALLOCATED ( A%ia ) )  CALL deallocate_array ( A%ia )
        IF ( ALLOCATED ( A%ja ) )  CALL deallocate_array ( A%ja )
        
    END SUBROUTINE deallocate_coo_matrix

    SUBROUTINE allocate_csr_matrix ( A, nrow, ncol, nnz )
        TYPE(csr_matrix), INTENT(INOUT) :: A
        INTEGER,          INTENT(IN)    :: nrow
        INTEGER,          INTENT(IN)    :: ncol
        INTEGER(KIND=eb), INTENT(IN)    :: nnz

        A%nrow = nrow
        A%ncol = ncol
        A%nnz  = nnz
        CALL allocate_array(A%val, nnz)
        CALL allocate_array(A%ia, nrow+1)
        CALL allocate_array(A%ja, nnz)
    END SUBROUTINE allocate_csr_matrix

    SUBROUTINE deallocate_csr_matrix(A)
        TYPE(csr_matrix) :: A

        IF ( ALLOCATED ( A%val ) ) CALL deallocate_array ( A%val )
        IF ( ALLOCATED ( A%ia ) )  CALL deallocate_array ( A%ia )
        IF ( ALLOCATED ( A%ja ) )  CALL deallocate_array ( A%ja )
    END SUBROUTINE deallocate_csr_matrix

    SUBROUTINE copy_coo_to_csr_matrix(Acsr, Acoo)
        TYPE(csr_matrix), INTENT(INOUT) :: Acsr
        TYPE(coo_matrix), INTENT(IN)    :: Acoo

        Acsr%nrow = Acoo%nrow
        Acsr%ncol = Acoo%ncol
        Acsr%nnz  = Acoo%nnz
        IF ( .NOT. ALLOCATED ( Acsr ) ) THEN
             CALL die('Programming error. Acsr has not been allocated.', 'copy_coo_to_csr_matrix')
!            CALL deallocate_matrix ( Acsr )
        ELSE IF ( SIZE ( Acsr%val ) /= SIZE ( Acoo%val ) ) THEN
            WRITE(*,*)  SIZE ( Acsr%val ), SIZE ( Acoo%val )
            CALL die('Sizes of ACSR and ACOO differ.', 'copy_coo_to_csr_matrix')
        ENDIF
!        CALL allocate_csr_matrix ( Acsr, Acsr%nrow, Acsr%ncol, Acsr%nnz )

        CALL coocsr ( Acsr%nrow, Acsr%nnz, Acoo%val, Acoo%ia, Acoo%ja, &
                                           Acsr%val, Acsr%ja, Acsr%ia )

    END SUBROUTINE copy_coo_to_csr_matrix

    SUBROUTINE copy_coo_to_dense_matrix ( D, A )
!
!       Purpose:
!       =======
!       Mainly for debugging purposes.
!
!       Date:            Programmer:           History of changes:
!       ====             ==========            ==================
!       Feb 2008         B.Strokopytov         Original code
!
        REAL(KIND=wp),    DIMENSION(:,:), ALLOCATABLE , INTENT(INOUT) :: D
        TYPE(coo_matrix),                 INTENT(IN)    :: A
!       Local variables:
        INTEGER                                         :: i
        CHARACTER(LEN=32)                               :: srname='convert_coo_to_dense_matrix'

!       Check proper allocation of dense matrix:
        IF ( .NOT. ALLOCATED ( D ) ) THEN
            CALL die('Programming error. Dense matrix D has not been allocated properly.', srname)
        ENDIF

!       Checkz:
        IF ( A%nrow /= SIZE (D, DIM=1) ) THEN

            WRITE(*,*) A%nrow, SIZE (D, DIM=1)
            CALL die('Programming error. Matrix sizes do not match.', srname)

        ELSE IF ( A%nnz /= SIZE ( A%val ) ) THEN

            WRITE(*,*) A%nnz, SIZE ( A%val )
            CALL die('Programming error. Inconsistent matrix parameters.', srname)

        ELSE IF ( MAXVAL ( A%ia ) > SIZE ( D, DIM=1 ) .OR. &
                  MAXVAL ( A%ja ) > SIZE ( D, DIM=2 ) ) THEN

            WRITE(*,*) MAXVAL ( A%ia ), SIZE ( D, DIM=1 )
            WRITE(*,*) MAXVAL ( A%ja ), SIZE ( D, DIM=2 )
            CALL die('Programming error. Matrix sizes do not match.', srname)

        ENDIF

!       Initialize D matrix for error-checking:
        D = 0.0_wp

!       Just copy:
        DO i = 1, A%nnz
            IF ( A%ia(i) > 0 .AND.  A%ja(i) > 0 ) THEN
                D(A%ia(i),A%ja(i)) = A%val(i)
            ELSE
                WRITE(*,*) A%ia(i), A%ja(i)
                CALL warn('Programming error. Bad indices.', srname)
            ENDIF
        ENDDO

    END SUBROUTINE copy_coo_to_dense_matrix

    SUBROUTINE condition_coo_matrix ( A, nnz, sk )
!
!       Purpose:
!       =======
!       Conditions sparse matrix in coordinate format
!
!       Date:             Programmer:              History of changes:
!       ====              ==========               ==================
!       Mar 2008          B.Strokopytov            Original code
!
!       Discussion:
!       ==========
!       On INPUT nnz is current number of non-zero elements.
!       This number may differ from A%nnz since we want some CPU economy during scaling.

        TYPE(coo_matrix),               INTENT(INOUT) :: A
        INTEGER(KIND=eb),               INTENT(IN)    :: nnz
        REAL(KIND=wp),    DIMENSION(:), INTENT(INOUT) :: sk
!       Local variables:
        INTEGER                                       :: ia
        INTEGER                                       :: ja        
!       Counters:
        INTEGER(KIND=eb)                              :: i
!       Name:
        CHARACTER(LEN=32)                             :: srname = 'condition_coo_matrix'
        
!       Checkz:
        IF ( SIZE ( A%val ) /= A%nnz ) THEN

            WRITE(*,*) SIZE ( A%val ), '/= ', A%nnz
            CALL die('Programming error. Inconsistent matrix parameters.', srname)

        ELSE IF ( nnz > A%nnz ) THEN

            WRITE(*,*) nnz, '> ', A%nnz
            CALL die('Programming error. NNZ exceeds available storage.', srname)

        ENDIF

        DO i = 1, nnz
            ia = A%ia(i)
            ja = A%ja(i)
            A%val(i) = A%val(i) * sk(ia) * sk(ja)            
        ENDDO

        CALL messag('Done.', srname)

    END SUBROUTINE condition_coo_matrix

    SUBROUTINE expand_coo_matrix ( H, np, nnz )
        TYPE(coo_matrix), INTENT(INOUT) :: H
        INTEGER,          INTENT(IN)    :: np
        INTEGER(KIND=eb), INTENT(IN)    :: nnz
!       Local variables:
        INTEGER                         :: i
        INTEGER                         :: j
!       Counters:
        INTEGER(KIND=eb)                :: elem
        INTEGER(KIND=eb)                :: l
        CHARACTER(LEN=32)               :: srname='expand_coo_matrix'
        
!       Initialize counter:
        elem = 0

!       Loop over asymmetric part:
        DO l = 1, nnz
            i = H%ia(l)
            j = H%ja(l)

!           Check for zero indices:
            IF ( i == 0 .OR. j== 0 ) THEN
                WRITE(*,*) i,j, ' from', l
                CALL die('Programming error. Zero indices for sparse matrix MUST be avoided.',&
                          srname)
            ENDIF

!           Do expansion for non-diagonal elements only:
            IF ( i /= j ) THEN

                elem = elem + 1

                matrix_size_check:IF ( nnz+elem <= H%nnz ) THEN

!                   It's very unlikely that in reality we can obtain zero matrix element,
!                   most probable explanation is programming error:
                    IF ( H%val(l) == 0.0_wp ) THEN
                        WRITE(*,*) i, j, H%val(l)
                        CALL warn('Possible programming error. Zero matrix element.', &
                                  srname)
                    ENDIF

!                   Check that current value is zero before expansion:
                    IF ( H%val(nnz+elem) == 0.0_wp ) THEN
                        H%val(nnz+elem) = H%val(l)
                        H%ia(nnz+elem) = j
                        H%ja(nnz+elem) = i
                    ELSE
!                       Non-zero value means we have already calculated this element:
                        WRITE(*,"(2I8)") nnz, elem
                        WRITE(*," (' SPARSE_NORMAL_MATRIX> ', 2I8, ' Unexpected matrix element=', ES9.2, &
                             & 'Current matrix element= ', ES9.2, ' index= ', I10)") &
                              j,i, H%val(nnz+elem), H%val(l), l
                        CALL die('Non-diagonal duplicate matrix element has been found.', &
                                  srname)
                    ENDIF

                ELSE

                    WRITE(*,*) H%nnz, nnz+elem
                    WRITE(*,*) j, i, H%val(l)
                    CALL die('Programming error. Matrix size is not sufficient to keep all the data.',&
                              srname)

                ENDIF matrix_size_check

            ENDIF
        ENDDO

        H%nnz = nnz + elem
        WRITE(*,"(' EXPAND_COO_MATRIX> ', A, ' elements has been added to ', A, ' elements')") &
        TRIM ( int_to_c ( elem ) ), TRIM ( int_to_c ( nnz ) )
        WRITE(*,"(' EXPAND_COO_MATRIX> ', 'Total number of matrix elements= ', A)") &
        TRIM ( int_to_c ( H%nnz ) )
        IF ( H%nnz /= SIZE ( H%val ) ) THEN
            WRITE(*,*) H%nnz, SIZE ( H%val )
            CALL die('Programming error. Total number of elemenets disagrees with the allocated size of H%VAL.',&
                      srname)
        ENDIF
    END SUBROUTINE expand_coo_matrix

    FUNCTION matrix_times_vector ( A, x )  RESULT ( y )
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE                :: y
        TYPE(csr_matrix),                            INTENT(IN)    :: A
        REAL(KIND=wp),    DIMENSION(:),              INTENT(IN)    :: x

        CALL allocate_array ( y, A%nrow )
        CALL amux(A%nrow, x, y, A%val, A%ja, A%ia)

    END FUNCTION matrix_times_vector 

    SUBROUTINE matvec ( n, y, x,  A)                  ! y := A*x
!
!       Purpose:
!       =======
!       Multiplies given matrix A by arbitrary vector x:

        INTEGER,                     INTENT(IN)  :: n
        REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: y
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: x
        TYPE(csr_matrix),            INTENT(IN)  :: A
!       Local variables:
        CHARACTER(LEN=32)                        :: srname = 'matvec'

!       Checkz:
        IF ( n /= A%nrow ) THEN
            WRITE(*,*) n, A%nrow
            CALL die('Programming error. Inconsistent parameters.', srname)
        ENDIF

!       Multiply:
        CALL amux(A%nrow, x, y, A%val,  A%ja, A%ia)
    END SUBROUTINE matvec

END MODULE sparse_basic
