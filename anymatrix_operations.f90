MODULE anymatrix_operations
USE constants
USE fail
USE mkl95_lapack, ONLY : syevd
USE select_kinds
USE util
IMPLICIT NONE
TYPE :: anymatrix
    REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: e
END TYPE

INTERFACE allocate_anymatrix
    MODULE PROCEDURE allocate_anymatrix_using_dims
    MODULE PROCEDURE allocate_anymatrix_using_array
END INTERFACE

INTERFACE ALLOCATED
    MODULE PROCEDURE allocated_anymatrix
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_of_anymatrices
END INTERFACE

INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE array_to_anymatrix
    MODULE PROCEDURE anymatrix_to_array
    MODULE PROCEDURE set_anymatrix_to_const
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_array_of_anymatrices
END INTERFACE

INTERFACE SIZE
    MODULE PROCEDURE size_of_anymatrix
    MODULE PROCEDURE size_of_anymatrix_with_dim
END INTERFACE

CONTAINS
    PURE FUNCTION allocated_anymatrix ( A )
        LOGICAL                     :: allocated_anymatrix
        TYPE(anymatrix), INTENT(IN) :: A
        allocated_anymatrix = ALLOCATED ( A%a )
    END FUNCTION allocated_anymatrix

    PURE FUNCTION size_of_anymatrix ( A )
        INTEGER                     :: size_of_anymatrix
        TYPE(anymatrix), INTENT(IN) :: A
        size_of_anymatrix = SIZE ( A%a )
    END FUNCTION size_of_anymatrix

    PURE FUNCTION size_of_anymatrix_with_dim ( A, DIM )
        INTEGER                     :: size_of_anymatrix_with_dim
        TYPE(anymatrix), INTENT(IN) :: A
        INTEGER,         INTENT(IN) :: DIM
        size_of_anymatrix_with_dim = SIZE ( A%a, DIM=DIM )
    END FUNCTION size_of_anymatrix_with_dim

    SUBROUTINE array_to_anymatrix ( A, arr )
        TYPE(anymatrix),               INTENT(INOUT) :: A
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)    :: arr

        IF ( .NOT. ALLOCATED ( A ) ) THEN
            CALL allocate_array(A%a, SIZE(arr,DIM=1), SIZE(arr,DIM=2) )
        ELSE IF ( SIZE ( A, DIM=1 ) /= SIZE ( arr, DIM=1 ) ) THEN
            CALL allocate_array(A%a, SIZE(arr,DIM=1), SIZE(arr,DIM=2) )
        ELSE IF ( SIZE ( A, DIM=2 ) /= SIZE ( arr, DIM=2 ) ) THEN
            CALL allocate_array(A%a, SIZE(arr,DIM=1), SIZE(arr,DIM=2) )
        ENDIF

        A%a = arr
    END SUBROUTINE array_to_anymatrix

    SUBROUTINE anymatrix_to_array ( arr, A )
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)    :: arr
        TYPE(anymatrix),                            INTENT(IN)       :: A
        IF ( .NOT. ALLOCATED ( arr ) ) THEN
            CALL allocate_array(arr, SIZE (A%a, DIM=1), SIZE(A%a, DIM=2))
        ENDIF
        arr = A%a
    END SUBROUTINE anymatrix_to_array

    SUBROUTINE allocate_anymatrix_using_dims ( A, n, m)
        TYPE(anymatrix),               INTENT(INOUT) :: A
        INTEGER,                       INTENT(IN)    :: n
        INTEGER,                       INTENT(IN)    :: m

        CALL allocate_array ( A%a, n, m )
        IF ( n == m ) THEN
            CALL allocate_array ( A%e, n )
            A%e = 0.0_wp
        ENDIF
    END SUBROUTINE allocate_anymatrix_using_dims

    SUBROUTINE allocate_anymatrix_using_array ( A, sizes )
        TYPE(anymatrix),       INTENT(INOUT) :: A
        INTEGER, DIMENSION(:), INTENT(IN)    :: sizes
        CALL allocate_anymatrix ( A, sizes(1), sizes(2) )
    END SUBROUTINE allocate_anymatrix_using_array

    SUBROUTINE deallocate_anymatrix ( A )
        TYPE(anymatrix),  INTENT(INOUT) :: A
        IF ( ALLOCATED ( A%a ) ) THEN
            CALL deallocate_array ( A%a )
        ENDIF
        IF ( ALLOCATED ( A%e ) ) THEN
            CALL deallocate_array ( A%e )
        ENDIF
    END SUBROUTINE deallocate_anymatrix

    SUBROUTINE allocate_array_of_anymatrices ( A, n )
        TYPE(anymatrix), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A
        INTEGER,                                    INTENT(IN)    :: n
!       Local variables:
        INTEGER                                                   :: istat

        ALLOCATE(A(n), STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL die('O-o-o-p-s-s. Failed to allocate array of matrices.', 'allocate_array_of_anymatrices')
        ENDIF
    END SUBROUTINE allocate_array_of_anymatrices

    SUBROUTINE deallocate_array_of_anymatrices ( A )
        TYPE(anymatrix), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A
!       Local variables:
        INTEGER                                      :: istat
!       Counters:
        INTEGER                                      :: i
        IF ( ALLOCATED ( A ) ) THEN
            DO i = 1, SIZE ( A )
                IF ( ALLOCATED ( A(i) ) ) CALL deallocate_anymatrix ( A(i) )
            ENDDO
        ENDIF
        DEALLOCATE ( A, STAT = istat )
        IF ( istat /= 0 ) THEN
            CALL die('O-o-o-p-s-s. Failed to deallocate array of matrices.', 'deallocate_array_of_anymatrices')
        ENDIF
    END SUBROUTINE deallocate_array_of_anymatrices

    SUBROUTINE set_anymatrix_to_const ( A, c )
        TYPE(anymatrix), INTENT(INOUT) :: A
        REAL(KIND=wp),   INTENT(IN)    :: c
        A%a = c
    END SUBROUTINE set_anymatrix_to_const

    SUBROUTINE set_anymatrix_elem (A, val, i, j)
        TYPE(anymatrix), INTENT(INOUT) :: A
        REAL(KIND=wp),   INTENT(IN)    :: val
        INTEGER,         INTENT(IN)    :: i
        INTEGER,         INTENT(IN)    :: j
        A%a(i,j) = val
    END SUBROUTINE set_anymatrix_elem

    SUBROUTINE get_anymatrix_elem (val, A, i, j)
        REAL(KIND=wp),   INTENT(INOUT) :: val
        TYPE(anymatrix), INTENT(IN)    :: A
        INTEGER,         INTENT(IN)    :: i
        INTEGER,         INTENT(IN)    :: j
        val = A%a(i,j)
    END SUBROUTINE get_anymatrix_elem

    SUBROUTINE print_anymatrix ( A )
        TYPE(anymatrix), INTENT(INOUT) :: A
!       Counters:
        INTEGER                        :: i
        INTEGER                        :: j

        DO i = 1, SIZE ( A%a, DIM=1 ) - 1
            DO j = 1,  SIZE ( A%a, DIM=2 )
                A%a(j,i) = A%a(i,j)
            ENDDO
        ENDDO
        WRITE(*,"(1X,'ANYMATRIX PRINT (extended by symmetry):')")
        WRITE(*,"(6X, 1000I9)") (i, i=1, SIZE ( A%a,DIM=1 ) )
        DO i = 1, SIZE ( A%a, DIM=1 )
            WRITE(*,"(1X,I5,1000ES9.2)") i, (A%a(i,j), j=1, SIZE ( A%a, DIM=2 ))
        ENDDO

    END SUBROUTINE print_anymatrix

    SUBROUTINE store_anymatrix_as_eigenvectors ( A )
        TYPE(anymatrix),             INTENT(INOUT) :: A
!       Local variables:
        INTEGER                                    :: i
        INTEGER                                    :: n
        INTEGER                                    :: m
        INTEGER                                    :: info
        CHARACTER(LEN=32),                    SAVE :: srname='store_anymatrix_as_eigenvectors'
        INTEGER(KIND=eb)                           :: mkl_memstat
        INTEGER(KIND=fb)                           :: allocated_buffers

        n = SIZE ( A, DIM=1 )
        m = SIZE ( A, DIM=2 )

        IF ( m /= n ) THEN
            WRITE(*,*) n, m
            CALL die('Programming error n /= m. What eigenvectors you are talking about?', srname )
        ENDIF

        IF ( .NOT. ALLOCATED ( A%e ) ) THEN
            CALL allocate_array ( A%e, SIZE ( A%a, DIM=1 ) )
        ENDIF

        CALL syevd(A%a, A%e, jobz='V', UPLO='U', info=info)
        
        IF ( debug > 10000 ) THEN
            WRITE(*,"(1X, A, ' MB of memory in ', A, ' allocated buffers')") &
            TRIM(int_to_c(MKL_MemStat(Allocated_Buffers)/2**20)), TRIM(int_to_c(Allocated_Buffers))
        ENDIF

        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('Illegal info value. Possible programming error.', srname)
        ENDIF
        
    END SUBROUTINE store_anymatrix_as_eigenvectors

    SUBROUTINE get_original_anymatrix ( O, A )
        REAL(KIND=wp),  DIMENSION(:,:), ALLOCATABLE,          INTENT(INOUT) :: O
        TYPE(anymatrix),                                      INTENT(IN)    :: A
!       Local automatic array:
        REAL(KIND=wp), DIMENSION(SIZE(A%a,DIM=1),SIZE(A%a,DIM=2))           :: W
!       Local variables:
        INTEGER                                                             :: n
        INTEGER                                                             :: m
        INTEGER                                                             :: ja
        CHARACTER(LEN=32),                                           SAVE   :: srname='get_original_anymatrix'


!       Figure out matrix dimensions:
        n = SIZE ( A, DIM=1 )
        m = SIZE ( A, DIM=2 )
        IF ( n /= m ) THEN
            WRITE(*,*) n, m
            CALL die('What eigenvectors you are talking about?', srname)
        ENDIF

        CALL allocate_array ( O, n,m )

        IF ( ANY ( ABS ( A%e ) == 0.0_wp ) ) THEN
            WRITE(*,*) MINLOC ( ABS ( A%e ) ), A%e(  MINLOC ( ABS ( A%e ) ) )
            CALL die('Inversion is not possible. Matrix contains zero eigenvalues.', srname)
        ENDIF

!       Scale eigenvector along columns:
        DO ja = 1, n
            W(1:n,ja) = A%a(1:n,ja) * A%e(ja)
        ENDDO

!       BTW FIXME: need overloaded transpose function:
        O = MATMUL ( W, TRANSPOSE ( A%a ) )

    END SUBROUTINE get_original_anymatrix

    SUBROUTINE get_inverted_anymatrix ( AINV, A)
        REAL(KIND=wp),  DIMENSION(:,:), ALLOCATABLE,          INTENT(INOUT) :: AINV
        TYPE(anymatrix),                                      INTENT(IN)    :: A
!       Local automatic array:
        REAL(KIND=wp), DIMENSION(SIZE(A%a,DIM=1),SIZE(A%a,DIM=2))           :: W
!       Local variables:
        INTEGER                                                             :: n
        INTEGER                                                             :: m
        INTEGER                                                             :: ja
        CHARACTER(LEN=32),                                           SAVE   :: srname='get_inverted_anymatrix'


!       Figure out matrix dimensions:
        n = SIZE ( A, DIM=1 )
        m = SIZE ( A, DIM=2 )
        IF ( n /= m ) THEN
            WRITE(*,*) n, m
            CALL die('What eigenvectors you are talking about?', srname)
        ENDIF

        CALL allocate_array ( AINV, n,m )

        IF ( ANY ( ABS ( A%e ) == 0.0_wp ) ) THEN
            WRITE(*,*) MINLOC ( ABS ( A%e ) ), A%e(  MINLOC ( ABS ( A%e ) ) )
            CALL die('Inversion is not possible. Matrix contains zero eigenvalues.', srname)
        ENDIF

!       Scale eigenvector along columns:
        DO ja = 1, n
            W(1:n,ja) = A%a(1:n,ja) / A%e(ja)
        ENDDO

!       BTW FIXME: need overloaded transpose function:
        AINV = MATMUL ( W, TRANSPOSE ( A%a ) )
    
    END SUBROUTINE get_inverted_anymatrix

    SUBROUTINE get_sqrt_of_inverted_anymatrix ( SQRT_AINV, A )
        REAL(KIND=wp),  DIMENSION(:,:), ALLOCATABLE,          INTENT(INOUT) :: SQRT_AINV
        TYPE(anymatrix),                                      INTENT(IN)    :: A
!       Local automatic array:
        REAL(KIND=wp), DIMENSION(SIZE(A,DIM=1),SIZE(A,DIM=2))               :: W
!       Local variables:
        INTEGER                                                             :: n
        INTEGER                                                             :: m
        INTEGER                                                             :: ja
        CHARACTER(LEN=32),                                           SAVE   :: srname='get_sqrt_of_inverted_anymatrix'
        
        
!       Figure out matrix dimensions:
        n = SIZE ( A, DIM=1 )
        m = SIZE ( A, DIM=2 )
        IF ( n /= m ) THEN
            WRITE(*,*) n, m
            CALL die('What eigenvectors you are talking about?', srname)
        ENDIF

        CALL allocate_array ( SQRT_AINV, n, m )

        IF ( ANY ( ABS ( A%e ) == 0.0_wp ) ) THEN
            WRITE(*,*) MINLOC ( ABS ( A%e ) ), A%e(  MINLOC ( ABS ( A%e ) ) )
            CALL die('Inversion is not possible. Matrix contains zero eigenvalues.', srname)
        ENDIF

!       Scale eigenvector along columns:
        DO ja = 1, n
            W(1:n,ja) = A%a(1:n,ja) / SQRT ( A%e(ja) )
        ENDDO

!       BTW FIXME: need overloaded transpose function:
        SQRT_AINV = MATMUL ( W, TRANSPOSE ( A%a ) )

    END SUBROUTINE get_sqrt_of_inverted_anymatrix

    SUBROUTINE get_sqrt_of_anymatrix ( SQRT_A, A )
        REAL(KIND=wp),  DIMENSION(:,:), ALLOCATABLE,          INTENT(INOUT) :: SQRT_A
        TYPE(anymatrix),                                      INTENT(IN)    :: A
!       Local automatic array:
        REAL(KIND=wp), DIMENSION(SIZE(A,DIM=1),SIZE(A,DIM=2))               :: W
!       Local variables:
        INTEGER                                                             :: n
        INTEGER                                                             :: m
        INTEGER                                                             :: ja
        CHARACTER(LEN=32),                                           SAVE   :: srname='store_sqrt_of_anymatrix'


!       Figure out matrix dimensions:
        n = SIZE ( A, DIM=1 )
        m = SIZE ( A, DIM=2 )
        IF ( n /= m ) THEN
            WRITE(*,*) n, m
            CALL die('What eigenvectors you are talking about?', srname)
        ENDIF

        IF ( ANY ( ABS ( A%e ) == 0.0_wp ) ) THEN
            WRITE(*,*) MINLOC ( ABS ( A%e ) ), A%e(  MINLOC ( ABS ( A%e ) ) )
            CALL die('Inversion is not possible. Matrix contains zero eigenvalues.', srname)
        ENDIF

!       Scale eigenvector along columns:
        DO ja = 1, n
            W(1:n,ja) = A%a(1:n,ja) * SQRT ( A%e(ja) )
        ENDDO

!       BTW FIXME: need overloaded transpose function:
        SQRT_A = MATMUL ( W, TRANSPOSE ( A%a ) )

    END SUBROUTINE get_sqrt_of_anymatrix


END MODULE anymatrix_operations
