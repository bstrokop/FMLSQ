MODULE dense_matrix_util
USE select_kinds
USE fail
IMPLICIT NONE
CONTAINS
    SUBROUTINE correct_potri(A, uplo)
    REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: A
    CHARACTER(LEN=*),                           INTENT(IN)    :: uplo
    INTEGER                                                   :: i
    INTEGER                                                   :: j
    INTEGER                                                   :: n

    n = SIZE ( A, DIM=1 )

    IF ( n /= SIZE ( A, DIM=1) ) THEN
        WRITE(*,*) n, SIZE ( A, DIM=1)
        CALL die('Matrix must be symmetrical.', 'correct_potri')
    ENDIF

    DO i = 1, n - 1
       DO j = i+1, n
           IF ( uplo == 'U' ) THEN
               A(j,i) = A(i,j)
           ELSE
               A(i,j) = A(j,i)
           ENDIF
       ENDDO
    ENDDO 

END SUBROUTINE correct_potri

    SUBROUTINE check_matrix_symmetry ( A )
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: A
    END SUBROUTINE check_matrix_symmetry

    SUBROUTINE average_symmetrical_matrix
    END SUBROUTINE average_symmetrical_matrix

END MODULE dense_matrix_util
