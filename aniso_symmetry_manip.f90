MODULE aniso_symmetry_manip
USE aniso_manip
USE basic_symmetry_operations
IMPLICIT NONE
TYPE :: aniso_symop
    PRIVATE
    REAL(KIND=wp), DIMENSION(6,6) :: a
END TYPE
PRIVATE
PUBLIC :: aniso_symop
PUBLIC :: allocate_array
PUBLIC :: deallocate_array
PUBLIC :: ASSIGNMENT (=)
PUBLIC :: OPERATOR (.SIXBYSIX.)
PUBLIC :: OPERATOR (*)
PUBLIC :: TRANSPOSE
PUBLIC :: aniso_symop_matmul
PRIVATE :: mult_aniso_symop_by_2d_array
PRIVATE :: mult_2d_array_by_aniso_symop
INTERFACE allocate_array
    MODULE PROCEDURE allocate_aniso_symop_array
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_aniso_symop_array
END INTERFACE

INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE convert_6x6_array_to_matrix
    MODULE PROCEDURE convert_aniso_symop_to_array
END INTERFACE

INTERFACE OPERATOR (.SIXBYSIX.)
     MODULE PROCEDURE convert_3x3_array_to_6x6_array
     MODULE PROCEDURE convert_symop_to_aniso_symop
     MODULE PROCEDURE convert_matrix_to_aniso_symop
END INTERFACE

INTERFACE OPERATOR (*)    
    MODULE PROCEDURE mult_aniso_symop_by_aniso
    MODULE PROCEDURE mult_aniso_symop_by_array
    MODULE PROCEDURE mult_array_by_aniso_symop
END INTERFACE

INTERFACE aniso_symop_matmul
!   Don't touch this:
    MODULE PROCEDURE mult_2d_array_by_aniso_symop
    MODULE PROCEDURE mult_aniso_symop_by_2d_array 
END INTERFACE 

INTERFACE TRANSPOSE
    MODULE PROCEDURE transpose_aniso_symop_matrix
END INTERFACE

CONTAINS
    SUBROUTINE allocate_aniso_symop_array (RU, n)
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: RU
        INTEGER,                                      INTENT(IN)    :: n
!       Local variables:
        INTEGER                                                     :: istat

!       Usual strategy (deallocate if allocated):
        IF ( ALLOCATED ( RU ) ) THEN
            CALL deallocate_array ( RU )
        ENDIF

        ALLOCATE ( RU(n), STAT=istat )

!       Test whether array has been allocated:
        IF ( istat /= 0 ) THEN
            CALL die('O-o-opps. Programming error. Failed to allocate RU array.',&
                     'allocate_aniso_symop_array')
        ENDIF

    END SUBROUTINE allocate_aniso_symop_array

    SUBROUTINE deallocate_aniso_symop_array(RU)
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: RU
!       Local variables:
        INTEGER                                                     :: istat
        
        IF ( ALLOCATED ( RU ) ) THEN
            DEALLOCATE ( RU, STAT=istat )
            IF ( istat /= 0 ) THEN
                CALL die('O-o-o-p-s-s... Programming error. Failed to deallocate RU array.', &
                         'deallocate_aniso_symop_array')
            ENDIF
        ENDIF
      
    END SUBROUTINE deallocate_aniso_symop_array

    SUBROUTINE convert_6x6_array_to_matrix (RU, R)
        TYPE(aniso_symop),              INTENT(OUT) :: RU
        REAL(KIND=wp),  DIMENSION(:,:), INTENT(IN)  :: R
        IF ( SIZE ( R, DIM=1 ) /= 6 .OR. SIZE ( R, DIM=2) /= 6)  THEN
            WRITE(*,*) SIZE ( R, DIM=1 ), SIZE ( R, DIM=2 )
            CALL die('Incorrect dims.', 'convert_6x6_array_to_matrix')
        ELSE
            RU%a = R
        ENDIF
    END SUBROUTINE convert_6x6_array_to_matrix

    SUBROUTINE convert_aniso_symop_to_array(R, RU)
        REAL(KIND=wp),  DIMENSION(6,6), INTENT(OUT) :: R
        TYPE(aniso_symop),              INTENT(IN)  :: RU
        R = RU%a
    END SUBROUTINE convert_aniso_symop_to_array

    FUNCTION convert_3x3_array_to_6x6_array (R) RESULT (RU)
! 
!       Purpose:
!       =======
!       Conversion of 3x3 array to 6x6 array to simplify
!       function differentiation according to Murshudov et al. (1999).
!
        REAL(KIND=wp), DIMENSION(6,6)             :: RU
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: R
!       Counters:
        INTEGER                                   :: i
        INTEGER                                   :: j

        IF ( SIZE ( R, DIM=1) /= 3 .OR. SIZE ( R, DIM=2 ) /= 3 ) THEN
             WRITE(*,*)  SIZE ( R, DIM=1), SIZE ( R, DIM=2 )
             CALL die('Programming error. Conversion is not possible.', &
                     'convert_3x3_array_to_6x6_array')
        ENDIF

        DO i = 1, 3
            DO j = 1, 3
!               In Garib's paper there is a printing error but this line below is correct BVS MAR 2009:
                RU(i,j) = R(i,j) * R(i,j)
            ENDDO

            RU(i,4) = 2.0_wp * R(i,1) * R(i,2)
            RU(i,5) = 2.0_wp * R(i,1) * R(i,3)
            RU(i,6) = 2.0_wp * R(i,2) * R(i,3)

            RU(4,i) = R(1,i) * R(2,i)
            RU(5,i) = R(1,i) * R(3,i)
            RU(6,i) = R(2,i) * R(3,i)    
        ENDDO

        RU(4,4) = R(1,1) * R(2,2) + R(1,2) * R(2,1)
        RU(5,4) = R(1,1) * R(3,2) + R(1,2) * R(3,1)
        RU(6,4) = R(2,1) * R(3,2) + R(2,2) * R(3,1)

        RU(4,5) = R(1,1) * R(2,3) + R(1,3) * R(2,1)
        RU(5,5) = R(1,1) * R(3,3) + R(1,3) * R(3,1)
        RU(6,5) = R(2,1) * R(3,3) + R(2,3) * R(3,1)

        RU(4,6) = R(1,2) * R(2,3) + R(1,3) * R(2,2)
        RU(5,6) = R(1,2) * R(3,3) + R(1,3) * R(3,2)
        RU(6,6) = R(2,2) * R(3,3) + R(2,3) * R(3,2)

    END FUNCTION convert_3x3_array_to_6x6_array

    FUNCTION convert_symop_to_aniso_symop (sym_1) RESULT (RU)
        TYPE(aniso_symop)                        :: RU
        TYPE(symop),                  INTENT(IN) :: sym_1
!       Local arrays:
        REAL(KIND=wp), DIMENSION(3,3)            :: temp

!       Extract rotational part converting to array type:
        temp = .SYMA. sym_1
        RU%a = convert_3x3_array_to_6x6_array (temp)                   
    END FUNCTION convert_symop_to_aniso_symop

    FUNCTION convert_matrix_to_aniso_symop (A) RESULT (RU)
        TYPE(aniso_symop)                        :: RU
        TYPE(matrix),                 INTENT(IN) :: A
!       Local arrays:
        REAL(KIND=wp), DIMENSION(3,3)            :: temp
        temp = A
        RU%a = convert_3x3_array_to_6x6_array (temp)        
    END FUNCTION convert_matrix_to_aniso_symop

    FUNCTION mult_aniso_symop_by_aniso(RU, U) RESULT ( RUU )
!
!       This routine is quite suitable for parallel computations
!       since no allocation of any arrays takes place...
        TYPE(aniso)                   :: RUU
        TYPE(aniso_symop), INTENT(IN) :: RU
        TYPE(aniso),       INTENT(IN) :: U
        RUU%u = MATMUL (RU%a, U%u)   
    END FUNCTION mult_aniso_symop_by_aniso

    FUNCTION mult_aniso_symop_by_array(RU, arr) RESULT ( RUU )
!
!       At present this routine is not quite suitable for parallel computations due
!       to allocation of array and requires "$OMP CRITICAL" instruction.
!
        REAL(KIND=wp),     DIMENSION(6)              :: RUU
        TYPE(aniso_symop),               INTENT(IN)  :: RU
        REAL(KIND=wp),     DIMENSION(:), INTENT(IN)  :: arr
        CHARACTER(LEN=32),                    SAVE   :: srname = 'mult_aniso_symop_by_array'
        IF ( SIZE ( arr ) == 6 ) THEN
!            CALL allocate_array ( RUU, 6 )
            RUU = MATMUL (RU%a, arr)
        ELSE
            CALL die('O-o-o-o-p-s-s. This operation is only permitted for arrays of size 6.', srname)
        ENDIF
    END FUNCTION mult_aniso_symop_by_array

    FUNCTION mult_array_by_aniso_symop (arr, RU) RESULT ( RUU )
        REAL(KIND=wp),     DIMENSION(6)              :: RUU
        REAL(KIND=wp),     DIMENSION(:), INTENT(IN)  :: arr
        TYPE(aniso_symop),               INTENT(IN)  :: RU
        CHARACTER(LEN=32),                    SAVE   :: srname = 'mult_aniso_symop_by_array'

        IF ( SIZE ( arr ) == 6 ) THEN
!           MATMUL easily allows that:
            RUU = MATMUL (arr, RU%a)
        ELSE
            CALL die('O-o-o-o-p-s-s. This operation is only permitted for arrays of size 6.', srname)
        ENDIF

    END FUNCTION mult_array_by_aniso_symop

    SUBROUTINE mult_2d_array_by_aniso_symop(RUU, arr, RU)
!
!       Main problem:
!       ============
!       RUU has variable size here but cannot be allocatable here
!       because of the problems with parallel computing.
!
!       We are lucky enough that MATMUL is thread-safe routine.
!
        REAL(KIND=wp),      DIMENSION(:,:), INTENT(INOUT) :: RUU
        REAL(KIND=wp),      DIMENSION(:,:), INTENT(IN)    :: arr
        TYPE(aniso_symop),                  INTENT(IN)    :: RU

        RUU = MATMUL ( arr, RU%a )
    END SUBROUTINE mult_2d_array_by_aniso_symop

    SUBROUTINE mult_aniso_symop_by_2d_array(RUU, RU, arr)
!
!       Main problem:
!       ============
!       RUU has variable size here but cannot be allocatable here
!       because of the problems with parallel computing.
!
!       We are lucky enough that MATMUL is thread-safe routine.
!
        REAL(KIND=wp),      DIMENSION(:,:), INTENT(INOUT) :: RUU
        TYPE(aniso_symop),                  INTENT(IN)    :: RU
        REAL(KIND=wp),      DIMENSION(:,:), INTENT(IN)    :: arr

        RUU = MATMUL ( RU%a, arr )
    END SUBROUTINE mult_aniso_symop_by_2d_array

    FUNCTION transpose_aniso_symop_matrix ( RU ) RESULT ( RUT )
        TYPE(aniso_symop)             :: RUT
        TYPE(aniso_symop)             :: RU
        RUT%a = TRANSPOSE ( RU%a )
    END FUNCTION transpose_aniso_symop_matrix                   
       
END MODULE aniso_symmetry_manip
