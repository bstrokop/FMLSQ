MODULE vectors
!
!
!   Purpose:
!   =======
!   To define a derived type called vector and matrix and the
!   operations that can be defined on them. The module
!   defines eight operations that can be performed on vectors and/or
!   matrices:
!
!                   Operation                                Operator
!                   =========                                ========
!    1. Creation from a wp precision array                  =
!    2. Conversion to a wp precision array                  =
!    3. Vector/matrix addition                                  +
!    4. Vector/matrix subtraction                               -
!    5. Vector/matrix/scalar multiplication                     *
!    6. Vector/matrix/scalar division                           /
!    7. Dot product                                           .DOT.
!    8. Cross product                                           *
!    9. Matrix transposition                                    No operator
!   10. Matrix determinant                                      No operator
!   11. Matrix inversion                                        /
!   12. Matrix inversion                                        **
!   13. Any matrix element > const                              >
!
!
!   It contains a total of >120 procedures to implement those operations:
!   wp_array_to_vector 
!   real_array_to_vector
!   ......................
!   inverse
!   inverse_int
!
!   Additional non-operator proceudres are:
!   det
!   det_int
!   transpose_of_the_inverse
!   transpose_of_the_inverse_int
!   print_mat
!   print_mat_int
!
!    Record of revisions:
!               Date           Programmer              Description of change
!               ====           ==========              =====================
!             01/05/96        S. J. Chapman            Original code
!             09/17/03        B. V. Strokopytov        Some wp precision procedures added
!                                                      for completeness
!             09/18/03        B. V. Strokopytov        "PRIVATE" line added
!             09/28/03        B. V. Strokopytov        Numerous additions of various wp/int functions
!             09/29/03        B. V. Strokopytov        Printing routines added
!             09/29/03        B. V. Strokopytov        108 MODULE PROCEDURES +det +inverse_of_the_transpose
!                                                                            +matrix printing routines
!                                                      112 functions/subroutines in total
!           
!             Oct 2005        B. V. Strokopytov        .DET. interface has been added
!             Oct 2005        B. V. Strokopytov        NINT interface added           
!             Oct 2005        B. V. Strokopytov        MOD  interface added
!             Oct 2005        B. V. Strokopytov        matrix_int_to_array added
!             Oct 2005        B. V. Strokopytov        vector, vector_int, matrix, matrix_int have PRIVATE components now
!             Oct 2005        B. V. Strokopytov        `double' in function names replaced with `wp'
!             Nov 2005        B. V. Strokopytov        ABS interface added
!             Nov 2005        B. V. Strokopytov        OPERATOR(>) rearranged
!             Nov 2005        B. V. Strokopytov        OPERATOR(<) added
!             Nov 2005        B. V. Strokopytov        MINVAL, MAXVAL interfaces added 
!             Nov 2005        B. V. Strokopytov        OPERATOR(.x.) added
!
!             Nov 2005        B. V. Strokopytov        vector_int_over_vector_int 
!                                                      vector_int_over_vector
!                                                      vector_over_vector_int
!                                                      vector_over_vector
!                                                      added to OPERATOR(/)
!
!             Nov 2005        B. V. Strokopytov        REAL interface added 
!             Nov 2005        B. V. Strokopytov        min_vec removed - MINVAL will do 
! ----------------------------------------------------------------------------------
!                                                      
USE constants
USE mkl95_lapack, ONLY : geevx, syevd
USE select_kinds
USE fail 
IMPLICIT NONE

! Declare vector data types:
TYPE :: vector
    PRIVATE
    REAL(KIND=wp), DIMENSION(3) :: v
END TYPE

! Declare vector_int data types:
TYPE :: vector_int
      PRIVATE
      INTEGER, DIMENSION(3)       :: v
END TYPE

! Declare matrix data types:
TYPE :: matrix_int
    PRIVATE
    INTEGER, DIMENSION(3,3)       :: a
END TYPE

TYPE :: matrix
    PRIVATE
    REAL(KIND=wp), DIMENSION(3,3) :: a
END TYPE

TYPE :: tensor
    PRIVATE
    REAL(KIND=wp)                 :: r11
    REAL(KIND=wp)                 :: r12
    REAL(KIND=wp)                 :: r13
    REAL(KIND=wp)                 :: r22
    REAL(KIND=wp)                 :: r23
    REAl(KIND=wp)                 :: r33
END TYPE


! Declare all items to be private except for type vector and
! the operators defined for it
PRIVATE
PUBLIC :: ABS
PUBLIC :: ASSIGNMENT (=)
PUBLIC :: eigen_values
PUBLIC :: extract_diagonal
PUBLIC :: INT
PUBLIC :: inverse
PUBLIC :: inverse_int
PUBLIC :: lmod
PUBLIC :: matrix
PUBLIC :: matrix_int
PUBLIC :: MAXVAL
PUBLIC :: MINVAL
PUBLIC :: MOD
PUBLIC :: NINT
PUBLIC :: OPERATOR ( == )
PUBLIC :: OPERATOR ( /= )
PUBLIC :: OPERATOR ( / )
PUBLIC :: OPERATOR ( * )
PUBLIC :: OPERATOR ( ** )
PUBLIC :: OPERATOR ( + )
PUBLIC :: OPERATOR ( - )
PUBLIC :: OPERATOR ( > )
PUBLIC :: OPERATOR ( < )
PUBLIC :: OPERATOR ( .DET. )
PUBLIC :: OPERATOR ( .DIAG. )
PUBLIC :: OPERATOR ( .DOT. )
PUBLIC :: OPERATOR ( .HMOD. )
PUBLIC :: OPERATOR ( .HTUMOD. )
PUBLIC :: OPERATOR ( .INV. )
PUBLIC :: OPERATOR ( .MMOD. )
PUBLIC :: print_mat
PUBLIC :: print_mat_int
PUBLIC :: print_tensor
PUBLIC :: REAL
PUBLIC :: rmod
PUBLIC :: SQRT
PUBLIC :: trace
PUBLIC :: TRANSPOSE
PUBLIC :: tensor
PUBLIC :: vector
PUBLIC :: vector_int
PUBLIC :: OPERATOR ( .X. )

! Declare interface operators:
INTERFACE ABS
    MODULE PROCEDURE abs_matrix
    MODULE PROCEDURE abs_vector
END INTERFACE

INTERFACE ASSIGNMENT ( = )
!   Vector part
!   Double/real operations:
    MODULE PROCEDURE wp_array_to_vector
    MODULE PROCEDURE real_array_to_vector
    MODULE PROCEDURE vector_to_array
    MODULE PROCEDURE set_vector_to_const
!   Integer operations:
    MODULE PROCEDURE array_int_to_vector_int
    MODULE PROCEDURE vector_int_to_array_int
    MODULE PROCEDURE set_vector_int_to_const
    MODULE PROCEDURE vector_int_to_vector_wp
    MODULE PROCEDURE vector_int_to_array_wp
!
!   Matrix part
!   Double/real operations:
    MODULE PROCEDURE wp_array_to_matrix
    MODULE PROCEDURE matrix_to_array
    MODULE PROCEDURE set_matrix_to_const
!   Integer operations:
    MODULE PROCEDURE array_int_to_matrix_int
    MODULE PROCEDURE matrix_int_to_array_int
    MODULE PROCEDURE matrix_int_to_array
    MODULE PROCEDURE set_matrix_int_to_const
    MODULE PROCEDURE matrix_int_to_matrix_wp
!   Tensor operations:
    MODULE PROCEDURE wp_array_to_tensor
END INTERFACE

INTERFACE OPERATOR ( .DIAG. )
    MODULE PROCEDURE set_matrix_diag_arr
    MODULE PROCEDURE set_matrix_diag_vec
    MODULE PROCEDURE set_matrix_diag_const
    MODULE PROCEDURE set_matrix_int_diag
    MODULE PROCEDURE set_matrix_int_diag_const
END INTERFACE

INTERFACE OPERATOR ( == )
!   Vector part:
    MODULE PROCEDURE vector_eq_zero
    MODULE PROCEDURE vector_eq_vector
!   Integer operations:
    MODULE PROCEDURE vector_int_eq_zero
    MODULE PROCEDURE vector_int_eq_vector_int
!   Matrix operations:
    MODULE PROCEDURE matrix_eq_matrix
    MODULE PROCEDURE matrix_eq_unity
!   Integer operations:
    MODULE PROCEDURE matrix_int_eq_matrix_int
    MODULE PROCEDURE matrix_int_eq_unity
END INTERFACE

INTERFACE OPERATOR ( /= )
!   Vector operations:
    MODULE PROCEDURE vector_ne_vector
    MODULE PROCEDURE vector_ne_zero
!   Integer operations:
    MODULE PROCEDURE vector_int_ne_vector_int
    MODULE PROCEDURE vector_int_ne_zero
!   Matrix operations:
    MODULE PROCEDURE matrix_ne_matrix
    MODULE PROCEDURE matrix_ne_unity
!   Integer operations:
    MODULE PROCEDURE matrix_int_ne_matrix_int
    MODULE PROCEDURE matrix_int_ne_unity
END INTERFACE

INTERFACE OPERATOR ( + )
!   Vector operations:
    MODULE PROCEDURE vector_add_vector
    MODULE PROCEDURE vector_int_add_vector_int
    MODULE PROCEDURE vector_int_add_vector
    MODULE PROCEDURE vector_add_vector_int
    MODULE PROCEDURE vector_add_const
    MODULE PROCEDURE vector_int_add_const
    MODULE PROCEDURE const_add_vector
    MODULE PROCEDURE const_add_vector_int
!   Matrix operations:
    MODULE PROCEDURE matrix_add_matrix
    MODULE PROCEDURE matrix_int_add_matrix_int
!   New 2009:
    MODULE PROCEDURE matrix_add_array
    MODULE PROCEDURE array_add_matrix
END INTERFACE

INTERFACE OPERATOR ( - )
!   Vector operations:
    MODULE PROCEDURE vector_minus
    MODULE PROCEDURE vector_subtract_vector
!   Integer operations:
    MODULE PROCEDURE vector_int_minus               ! one argument
    MODULE PROCEDURE vector_int_subtract_vector_int
    MODULE PROCEDURE vector_int_subtract_vector
    MODULE PROCEDURE vector_subtract_vector_int
!   Matrix operations:
    MODULE PROCEDURE matrix_minus
    MODULE PROCEDURE matrix_subtract_matrix
!   New 2009:
    MODULE PROCEDURE matrix_subtract_array
    MODULE PROCEDURE array_subtract_matrix
!   Integer operations:
    MODULE PROCEDURE matrix_int_minus
    MODULE PROCEDURE matrix_int_subtract_matrix_int
END INTERFACE

INTERFACE OPERATOR ( * )
!   Vector part:
    MODULE PROCEDURE vector_times_real
    MODULE PROCEDURE real_times_vector
    MODULE PROCEDURE wp_times_vector
    MODULE PROCEDURE vector_times_wp
    MODULE PROCEDURE int_times_vector
    MODULE PROCEDURE vector_times_int
    MODULE PROCEDURE cross_product
!   Integer vector part:
    MODULE PROCEDURE vector_int_times_int
    MODULE PROCEDURE int_times_vector_int
    MODULE PROCEDURE vector_int_times_real
    MODULE PROCEDURE real_times_vector_int
    MODULE PROCEDURE vector_int_times_wp
    MODULE PROCEDURE wp_times_vector_int
!   Matrix part:
    MODULE PROCEDURE matrix_times_matrix
    MODULE PROCEDURE matrix_times_vector
    MODULE PROCEDURE matrix_times_int
    MODULE PROCEDURE matrix_times_real
    MODULE PROCEDURE matrix_times_wp
    MODULE PROCEDURE int_times_matrix
    MODULE PROCEDURE real_times_matrix
    MODULE PROCEDURE wp_times_matrix
    MODULE PROCEDURE matrix_times_array_wp
!   Integer matrix part:
    MODULE PROCEDURE matrix_int_times_matrix_int
    MODULE PROCEDURE matrix_int_times_array_int
    MODULE PROCEDURE matrix_int_times_array_wp
    MODULE PROCEDURE matrix_int_times_vector_int
    MODULE PROCEDURE matrix_int_times_matrix
    MODULE PROCEDURE int_times_matrix_int
    MODULE PROCEDURE real_times_matrix_int
    MODULE PROCEDURE wp_times_matrix_int
    MODULE PROCEDURE matrix_int_times_int
    MODULE PROCEDURE matrix_int_times_real
    MODULE PROCEDURE matrix_int_times_wp
!   Mixed matrix part:
    MODULE PROCEDURE matrix_times_matrix_int
    MODULE PROCEDURE matrix_times_vector_int
    MODULE PROCEDURE matrix_int_times_vector
!   Tensor operations:
    MODULE PROCEDURE tensor_times_vector
    MODULE PROCEDURE tensor_times_vector_int
    MODULE PROCEDURE tensor_times_array
    MODULE PROCEDURE tensor_times_array_int
END INTERFACE

INTERFACE OPERATOR ( .x. )
    MODULE PROCEDURE vector_times_vector
    MODULE PROCEDURE vector_int_times_vector
    MODULE PROCEDURE vector_times_vector_int
    MODULE PROCEDURE vector_int_times_vector_int
END INTERFACE

INTERFACE OPERATOR ( / )
!   Vector part:
    MODULE PROCEDURE vector_div_wp
    MODULE PROCEDURE vector_div_real
    MODULE PROCEDURE vector_div_int
    MODULE PROCEDURE wp_div_vector
    MODULE PROCEDURE real_div_vector
!   Integer vector part:
    MODULE PROCEDURE vector_int_div_int
    MODULE PROCEDURE vector_int_div_real
    MODULE PROCEDURE vector_int_div_wp
    MODULE PROCEDURE wp_div_vector_int
    MODULE PROCEDURE real_div_vector_int
    MODULE PROCEDURE vector_over_vector
    MODULE PROCEDURE vector_over_vector_int
    MODULE PROCEDURE vector_int_over_vector
    MODULE PROCEDURE vector_int_over_vector_int
!   Matrix part:
    MODULE PROCEDURE matrix_div_wp
    MODULE PROCEDURE matrix_div_real
    MODULE PROCEDURE matrix_div_int
    MODULE PROCEDURE int_div_matrix
    MODULE PROCEDURE real_div_matrix
    MODULE PROCEDURE wp_div_matrix
    MODULE PROCEDURE matrix_div_matrix
!   Integer matrix part:
    MODULE PROCEDURE matrix_int_div_matrix_int
    MODULE PROCEDURE matrix_int_div_int
    MODULE PROCEDURE matrix_int_div_real
    MODULE PROCEDURE matrix_int_div_wp
    MODULE PROCEDURE int_div_matrix_int
    MODULE PROCEDURE real_div_matrix_int
    MODULE PROCEDURE wp_div_matrix_int
END INTERFACE

INTERFACE OPERATOR ( > )
    MODULE PROCEDURE matrix_gt_wp_const
    MODULE PROCEDURE vector_gt_wp_const
END INTERFACE

INTERFACE OPERATOR ( < )
    MODULE PROCEDURE matrix_lt_wp_const
    MODULE PROCEDURE vector_lt_wp_const
END INTERFACE

INTERFACE OPERATOR ( .DOT. )
!   Vector part only:
    MODULE PROCEDURE dot_product_wp
    MODULE PROCEDURE dot_product_int
!   Mixed operations:
    MODULE PROCEDURE dot_product_int_wp
    MODULE PROCEDURE dot_product_wp_int
!   Array operation:
    MODULE PROCEDURE vec_int_2d_dot_product_vec
END INTERFACE

INTERFACE OPERATOR ( ** )
    MODULE PROCEDURE vector_inverse
END INTERFACE

INTERFACE OPERATOR ( .INV. )
    MODULE PROCEDURE inverse_int
    MODULE PROCEDURE inverse
END INTERFACE 

INTERFACE MINVAL
    MODULE PROCEDURE minval_vector
    MODULE PROCEDURE minval_matrix
END INTERFACE

INTERFACE MAXVAL
    MODULE PROCEDURE maxval_vector
    MODULE PROCEDURE maxval_matrix
END INTERFACE

INTERFACE TRANSPOSE
    MODULE PROCEDURE transpose_wp
    MODULE PROCEDURE transpose_int
END INTERFACE

INTERFACE OPERATOR ( .DET. )
    MODULE PROCEDURE det
    MODULE PROCEDURE det_int
END INTERFACE 

INTERFACE OPERATOR ( .MMOD. )
    MODULE PROCEDURE vec_int_mod_vec_int
    MODULE PROCEDURE vec_int_mod_arr_int
    MODULE PROCEDURE arr_int_mod_arr_int
    MODULE PROCEDURE arr_int_mod_vec_int
END INTERFACE

INTERFACE OPERATOR ( .HTUMOD. )
    MODULE PROCEDURE h2mod
END INTERFACE

INTERFACE OPERATOR ( .HMOD. )
    MODULE PROCEDURE hmod
END INTERFACE

INTERFACE INT
    MODULE PROCEDURE int_matrix
    MODULE PROCEDURE int_vector
END INTERFACE

INTERFACE NINT
    MODULE PROCEDURE nint_matrix
    MODULE PROCEDURE nint_vector
END INTERFACE

INTERFACE MOD
    MODULE PROCEDURE mod_vector
    MODULE PROCEDURE mod_vector_int
END INTERFACE

INTERFACE REAL
    MODULE PROCEDURE real_wp_vector
    MODULE PROCEDURE real_wp_matrix
END INTERFACE 

INTERFACE SQRT
    MODULE PROCEDURE sqrt_vector
END INTERFACE

! Now define the implementing functions:
CONTAINS
    SUBROUTINE wp_array_to_vector(vec_result, array)
        TYPE(vector),  INTENT(OUT)              :: vec_result
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: array
!       Local variables:
        INTEGER, DIMENSION(1)                   :: extent

        extent = UBOUND ( array ) - LBOUND ( array ) + 1
        IF ( extent(1) < 3 ) THEN
            CALL die('Extent of the input array < 3','wp_array_to_vector')
        ELSE
            vec_result%v = array(1:3)
        ENDIF
    END SUBROUTINE wp_array_to_vector
   
    SUBROUTINE real_array_to_vector(vec_result, array)
        TYPE(vector),                 INTENT(OUT) :: vec_result
        REAL(KIND=sp) , DIMENSION(:), INTENT(IN)  :: array
!       Local variables:
        INTEGER,        DIMENSION(1)              :: extent

        extent =  UBOUND ( array ) - LBOUND ( array ) + 1
        IF ( extent(1) < 3 ) THEN
            CALL die ( 'Extent of the input array < 3','real_array_to_vector')
        ELSE
            vec_result%v = array(1:3)
        ENDIF
    END SUBROUTINE real_array_to_vector

    PURE SUBROUTINE array_int_to_vector_int ( vec_result_int, array_int )
        TYPE(vector_int),      INTENT(OUT) :: vec_result_int
        INTEGER, DIMENSION(:), INTENT(IN)  :: array_int

        vec_result_int%v = array_int(1:3)
    END SUBROUTINE array_int_to_vector_int

    PURE SUBROUTINE vector_int_to_array_int ( array_result_int, vec_1_int )
        TYPE(vector_int),      INTENT(IN)  :: vec_1_int
        INTEGER, DIMENSION(:), INTENT(OUT) :: array_result_int

        array_result_int(1:3) = vec_1_int%v
    END SUBROUTINE vector_int_to_array_int

    SUBROUTINE vector_to_array(array_result, vec_1)
        REAL(KIND=wp), DIMENSION(:), INTENT(OUT) :: array_result
        TYPE(vector),                INTENT(IN)  :: vec_1

        array_result(1:3) = vec_1%v
    END SUBROUTINE vector_to_array

    SUBROUTINE wp_array_to_matrix ( mat_result, array )
        TYPE(matrix),                  INTENT(OUT) :: mat_result
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)  :: array
!       Local variables:
        INTEGER,       DIMENSION(2)                :: extent

        extent = UBOUND ( array ) - LBOUND ( array ) + 1
        IF ( extent(1) /= 3 .OR. extent(2) /= 3 ) THEN
            CALL die ('Input array dims must be 3x3','wp_array_to_matrix')
        ELSE
            mat_result%a = array
        ENDIF
    END SUBROUTINE wp_array_to_matrix

    SUBROUTINE array_int_to_matrix_int( mat_result, array_int)
         TYPE(matrix_int),        INTENT(OUT) :: mat_result
         INTEGER, DIMENSION(:,:), INTENT(IN)  :: array_int
!        Local variables:
         INTEGER, DIMENSION(2)                :: extent

         extent = UBOUND(array_int) - LBOUND(array_int) + 1
         IF ( extent(1) /= 3 .OR. extent(2) /= 3 ) THEN
             CALL die ('Input array dims must be 3x3','array_int_to_matrix_int')
         ELSE
             mat_result%a = array_int
         ENDIF     
    END SUBROUTINE array_int_to_matrix_int

    SUBROUTINE set_matrix_int_to_const ( mat_result, constant )
         TYPE(matrix_int), INTENT(OUT) :: mat_result
         INTEGER, INTENT(IN)           :: constant

         mat_result%a = constant
    END SUBROUTINE set_matrix_int_to_const

    SUBROUTINE set_vector_int_to_const ( vec_result, constant )
        TYPE(vector_int), INTENT(OUT) :: vec_result
        INTEGER,          INTENT(IN)  :: constant

        vec_result%v = constant
    END SUBROUTINE set_vector_int_to_const

    SUBROUTINE set_vector_to_const ( vec_result, constant )
        TYPE(vector),  INTENT(OUT) :: vec_result
        REAL(KIND=wp), INTENT(IN)  :: constant

        vec_result%v = constant
    END SUBROUTINE set_vector_to_const

    SUBROUTINE set_matrix_to_const ( mat_result, constant )
         TYPE(matrix), INTENT(OUT) :: mat_result
         REAL(KIND=wp), INTENT(IN) :: constant

         mat_result%a = constant 
    END SUBROUTINE set_matrix_to_const

!   OPERATOR ( .DIAG. )
    FUNCTION set_matrix_diag_arr ( arr ) RESULT ( mat_result )
        TYPE(matrix)                             :: mat_result
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: arr
!       Local variables:
        INTEGER,       DIMENSION(1)              :: extent
        INTEGER                                  :: i
!
        extent = UBOUND ( arr ) - LBOUND ( arr ) + 1
        IF ( extent(1) < 3 ) CALL die ( 'Too small array to set a diagonal', 'set_matrix_diag' )

!       Set matrix to zero:  
        mat_result%a = 0.0_wp

!       Set matrix main diagonal:
        DO i = 1, 3
            mat_result%a(i,i) = arr(i)
        ENDDO
    END FUNCTION set_matrix_diag_arr

    FUNCTION set_matrix_diag_vec ( vec_1 ) RESULT ( mat_result )
        TYPE(matrix)             :: mat_result
        TYPE(vector), INTENT(IN) :: vec_1

        mat_result = set_matrix_diag_arr ( vec_1%v )         
    END FUNCTION set_matrix_diag_vec

    FUNCTION set_matrix_diag_const ( const ) RESULT ( mat_result )
        TYPE(matrix)               :: mat_result
        REAL(KIND=wp), INTENT(IN)  :: const
!       Counters:
        INTEGER                    :: i

!       Set matrix to zero:
        mat_result%a = 0.0_wp

!       Set diag:
        DO i = 1, 3
            mat_result%a(i,i) = const
        ENDDO
    END FUNCTION set_matrix_diag_const

    FUNCTION set_matrix_int_diag ( arr ) RESULT ( mat_result )
        TYPE(matrix_int)                  :: mat_result
        INTEGER, DIMENSION(:), INTENT(IN) :: arr
!       Local variales:
        INTEGER, DIMENSION(1)             :: extent
!       Counters:
        INTEGER                           :: i

!       Checks:
        extent = UBOUND ( arr ) - LBOUND ( arr ) + 1
        IF ( extent(1) < 3 ) CALL die ( 'Too small array to set a diagonal', 'set_matrix_int_diag' )

        mat_result%a = 0
        DO i = 1, 3
            mat_result%a(i,i) = arr(i)
        ENDDO
    END FUNCTION set_matrix_int_diag

    FUNCTION set_matrix_int_diag_const ( const ) RESULT ( mat_result )
        TYPE(matrix_int)                  :: mat_result
        INTEGER,               INTENT(IN) :: const
!       Counters:
        INTEGER                           :: i

        mat_result%a = 0
        DO i = 1, 3
            mat_result%a(i,i) = const
        ENDDO
    END FUNCTION set_matrix_int_diag_const

    SUBROUTINE matrix_int_to_array_int ( array_int_result, mat_int )
         INTEGER, DIMENSION(3,3), INTENT(OUT) :: array_int_result
         TYPE(matrix_int),        INTENT(IN)  :: mat_int

         array_int_result = mat_int%a
    END SUBROUTINE matrix_int_to_array_int

    SUBROUTINE matrix_int_to_array( array_int_result, mat_int )
         REAL(KIND=wp), DIMENSION(3,3), INTENT(OUT) :: array_int_result
         TYPE(matrix_int),              INTENT(IN)  :: mat_int

         array_int_result = mat_int%a
    END SUBROUTINE matrix_int_to_array

    SUBROUTINE matrix_to_array ( array_result, mat )
         REAL(KIND=wp), DIMENSION(3,3), INTENT(OUT) :: array_result
         TYPE(matrix),                  INTENT(IN)  :: mat

         array_result = mat%a
    END SUBROUTINE matrix_to_array

    SUBROUTINE matrix_int_to_matrix_wp ( mat_result, mat )
        TYPE(matrix),     INTENT(OUT) :: mat_result
        TYPE(matrix_int), INTENT(IN)  :: mat

        mat_result%a = mat%a
    END SUBROUTINE matrix_int_to_matrix_wp

    SUBROUTINE vector_int_to_vector_wp ( vec_result, vec )
        TYPE(vector),     INTENT(OUT) :: vec_result
        TYPE(vector_int), INTENT(IN)  :: vec

        vec_result%v = vec%v 
    END SUBROUTINE vector_int_to_vector_wp

!   Small addition for completeness of '=' interface (OCT 2007):
    SUBROUTINE vector_int_to_array_wp ( vec_result, vec )
        REAL(KIND=wp),    DIMENSION(3), INTENT(INOUT) :: vec_result
        TYPE(vector_int),               INTENT(IN)    :: vec

        vec_result = vec%v
    END SUBROUTINE vector_int_to_array_wp

    FUNCTION matrix_int_eq_matrix_int ( mat_1, mat_2 )
        LOGICAL                      :: matrix_int_eq_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1, mat_2

        matrix_int_eq_matrix_int = ALL( mat_1%a == mat_2%a)
    END FUNCTION matrix_int_eq_matrix_int

    FUNCTION matrix_int_ne_matrix_int ( mat_1, mat_2 )
        LOGICAL                      :: matrix_int_ne_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_ne_matrix_int = ANY( mat_1%a /= mat_2%a)
    END FUNCTION matrix_int_ne_matrix_int


    FUNCTION matrix_eq_matrix ( mat_1, mat_2 )
        LOGICAL                  :: matrix_eq_matrix
        TYPE(matrix), INTENT(IN) :: mat_1
        TYPE(matrix), INTENT(IN) :: mat_2

        matrix_eq_matrix = ALL ( ABS(mat_1%a - mat_2%a) < eps )
    END FUNCTION matrix_eq_matrix

    FUNCTION matrix_ne_matrix ( mat_1, mat_2 )
        LOGICAL                  :: matrix_ne_matrix
        TYPE(matrix), INTENT(IN) :: mat_1
        TYPE(matrix), INTENT(IN) :: mat_2

        matrix_ne_matrix = ANY ( ABS(mat_1%a - mat_2%a) > eps )
    END FUNCTION matrix_ne_matrix

    FUNCTION vector_int_eq_vector_int ( vec_1, vec_2 )
        LOGICAL                      :: vector_int_eq_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        vector_int_eq_vector_int = ALL( vec_1%v == vec_2%v )
    END FUNCTION vector_int_eq_vector_int

    FUNCTION vector_int_ne_vector_int (vec_1, vec_2 )
        LOGICAL                      :: vector_int_ne_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        vector_int_ne_vector_int = ANY( vec_1%v /= vec_2%v )
    END FUNCTION vector_int_ne_vector_int

    FUNCTION vector_eq_vector ( vec_1, vec_2 )
        LOGICAL                  :: vector_eq_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_eq_vector = ALL ( ABS ( vec_1%v - vec_2%v ) < eps )
    END FUNCTION vector_eq_vector

    FUNCTION vector_ne_vector ( vec_1, vec_2 )
        LOGICAL                  :: vector_ne_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_ne_vector = ANY ( ABS ( vec_1%v - vec_2%v ) > eps )
    END FUNCTION vector_ne_vector

    FUNCTION matrix_eq_unity (mat_1, wp_2)
        LOGICAL                   :: matrix_eq_unity
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=wp), INTENT(IN) :: wp_2
!       Counters:
        INTEGER                   :: i
        INTEGER                   :: j

        DO i = 1, 3
            DO j = 1, 3
                IF ( i == j ) THEN
                    IF ( mat_1%a(i,j) /= wp_2 ) THEN
                        matrix_eq_unity = .FALSE.
                        RETURN
                    ENDIF                  
                ELSE
                    IF ( mat_1%a(i,j) /= 0.0_wp ) THEN
                        matrix_eq_unity =.FALSE.
                        RETURN
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        matrix_eq_unity = .TRUE.
    END FUNCTION matrix_eq_unity       

    FUNCTION matrix_ne_unity ( mat_1, wp_2 )
        LOGICAL                   :: matrix_ne_unity
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

        matrix_ne_unity = .NOT. matrix_eq_unity ( mat_1, wp_2 )
    END FUNCTION matrix_ne_unity

    FUNCTION matrix_int_eq_unity (mat_1, int_2)
        LOGICAL                      :: matrix_int_eq_unity
        TYPE(matrix_int), INTENT(IN) :: mat_1
        INTEGER,          INTENT(IN) :: int_2
!       Counters:
        INTEGER                      :: i
        INTEGER                      :: j

        DO i = 1, 3
            DO j = 1, 3
                IF ( i == j ) THEN
                    IF ( mat_1%a(i,j) /= int_2 ) THEN
                        matrix_int_eq_unity = .FALSE.
                        RETURN
                    ENDIF
                ELSE
                    IF ( mat_1%a(i,j) /= 0) THEN
                        matrix_int_eq_unity = .FALSE.
                        RETURN
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        matrix_int_eq_unity = .TRUE.
    END FUNCTION matrix_int_eq_unity

    FUNCTION matrix_int_ne_unity ( mat_1, int_2 )
        LOGICAL                       :: matrix_int_ne_unity
        TYPE (matrix_int), INTENT(IN) :: mat_1
        INTEGER, INTENT(IN)           :: int_2

        matrix_int_ne_unity = .NOT. matrix_int_eq_unity ( mat_1, int_2 )
    END FUNCTION matrix_int_ne_unity

    FUNCTION vector_eq_zero ( vec_1, wp_2 )
        LOGICAL                   :: vector_eq_zero
        TYPE (vector), INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

        vector_eq_zero = ALL ( vec_1%v == wp_2 )
    END FUNCTION vector_eq_zero

    FUNCTION vector_ne_zero ( vec_1, wp_2 )
        LOGICAL                   :: vector_ne_zero
        TYPE (vector), INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

        vector_ne_zero = ANY ( vec_1%v /= wp_2 ) 
    END FUNCTION vector_ne_zero

    FUNCTION vector_int_eq_zero ( vec_1, int_2 )
        LOGICAL                       :: vector_int_eq_zero
        TYPE(vector_int), INTENT(IN) :: vec_1
        INTEGER, INTENT(IN)           :: int_2

        vector_int_eq_zero = ALL ( vec_1%v == int_2 )
    END FUNCTION vector_int_eq_zero

    FUNCTION vector_int_ne_zero ( vec_1, int_2 )
        LOGICAL                       :: vector_int_ne_zero
        TYPE(vector_int), INTENT(IN) :: vec_1
        INTEGER, INTENT(IN)           :: int_2

        vector_int_ne_zero = ANY ( vec_1%v /= int_2 ) 
    END FUNCTION vector_int_ne_zero
!
!   (+) operations
    FUNCTION vector_add_vector(vec_1, vec_2)
        TYPE(vector)             :: vector_add_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_add_vector%v = vec_1%v + vec_2%v
    END FUNCTION vector_add_vector

    FUNCTION vector_int_add_vector_int(vec_1_int, vec_2_int)
        TYPE(vector_int)             :: vector_int_add_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_int_add_vector_int%v = vec_1_int%v + vec_2_int%v
    END FUNCTION vector_int_add_vector_int

    FUNCTION matrix_add_matrix ( mat_1, mat_2 )
        TYPE(matrix)             :: matrix_add_matrix
        TYPE(matrix), INTENT(IN) :: mat_1
        TYPE(matrix), INTENT(IN) :: mat_2

        matrix_add_matrix%a = mat_1%a + mat_2%a
    END FUNCTION matrix_add_matrix

    FUNCTION matrix_int_add_matrix_int ( mat_1, mat_2 )
        TYPE(matrix_int)             :: matrix_int_add_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_add_matrix_int%a = mat_1%a + mat_2%a
    END FUNCTION matrix_int_add_matrix_int

    FUNCTION matrix_add_array ( mat_1, arr ) RESULT(a)
        REAL(KIND=wp), DIMENSION(3,3) :: a
        TYPE(matrix),                  INTENT(IN) :: mat_1
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: arr
        a = mat_1%a + arr
    END FUNCTION matrix_add_array

    FUNCTION array_add_matrix ( arr, mat_1 ) RESULT(a)
        REAL(KIND=wp), DIMENSION(3,3) :: a
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: arr
        TYPE(matrix),                  INTENT(IN) :: mat_1
        a = arr + mat_1%a
    END FUNCTION

    FUNCTION matrix_int_subtract_matrix_int ( mat_1, mat_2 )
        TYPE(matrix_int)             :: matrix_int_subtract_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_subtract_matrix_int%a = mat_1%a - mat_2%a
    END FUNCTION matrix_int_subtract_matrix_int

    FUNCTION matrix_subtract_matrix ( mat_1, mat_2 )
        TYPE(matrix)             :: matrix_subtract_matrix
        TYPE(matrix), INTENT(IN) :: mat_1
        TYPE(matrix), INTENT(IN) :: mat_2

        matrix_subtract_matrix%a = mat_1%a - mat_2%a
    END FUNCTION matrix_subtract_matrix

    FUNCTION matrix_int_minus ( mat_2 )
        TYPE(matrix_int)             :: matrix_int_minus
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_minus%a = -mat_2%a
    END FUNCTION matrix_int_minus

    FUNCTION matrix_subtract_array ( mat_1, arr ) RESULT(a)
        REAL(KIND=wp), DIMENSION(3,3)             :: a
        TYPE(matrix),                  INTENT(IN) :: mat_1
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: arr
        a = mat_1%a - arr
    END FUNCTION matrix_subtract_array

    FUNCTION array_subtract_matrix ( arr, mat_1 ) RESULT(a)
        REAL(KIND=wp), DIMENSION(3,3)             :: a
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: arr
        TYPE(matrix),                  INTENT(IN) :: mat_1
        a = arr - mat_1%a
    END FUNCTION array_subtract_matrix

    FUNCTION matrix_minus ( mat_2 )
        TYPE(matrix)             :: matrix_minus
        TYPE(matrix), INTENT(IN) :: mat_2

        matrix_minus%a = -mat_2%a
    END FUNCTION matrix_minus

    FUNCTION vector_subtract_vector ( vec_1, vec_2 )
        TYPE(vector)             :: vector_subtract_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_subtract_vector%v = vec_1%v - vec_2%v
    END FUNCTION vector_subtract_vector

    FUNCTION vector_int_subtract_vector_int ( vec_1_int, vec_2_int )
        TYPE(vector_int)             :: vector_int_subtract_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_int_subtract_vector_int%v = vec_1_int%v - vec_2_int%v
    END FUNCTION vector_int_subtract_vector_int

    FUNCTION vector_int_subtract_vector ( vec_1_int, vec_2 )
        TYPE(vector)                 :: vector_int_subtract_vector
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector),     INTENT(IN) :: vec_2

        vector_int_subtract_vector%v = vec_1_int%v - vec_2%v 
    END FUNCTION vector_int_subtract_vector

    FUNCTION vector_subtract_vector_int ( vec_1, vec_2_int )
        TYPE(vector)                 :: vector_subtract_vector_int
        TYPE(vector), INTENT(IN)     :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_subtract_vector_int%v = vec_1%v - vec_2_int%v
    END FUNCTION vector_subtract_vector_int

    FUNCTION vector_int_add_vector ( vec_1_int, vec_2)
        TYPE(vector)                 :: vector_int_add_vector
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector),     INTENT(IN) :: vec_2

        vector_int_add_vector%v = REAL(vec_1_int%v, KIND=wp) + vec_2%v
    END FUNCTION vector_int_add_vector

    FUNCTION vector_add_vector_int ( vec_1, vec_2_int)
        TYPE(vector)                 :: vector_add_vector_int
        TYPE(vector),     INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_add_vector_int%v = vec_1%v + REAL(vec_2_int%v, KIND=wp)
    END FUNCTION vector_add_vector_int

!   End of vector mixed operations

    FUNCTION vector_add_const ( vec_1, const_2 )
        TYPE(vector)              :: vector_add_const
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: const_2

        vector_add_const%v = vec_1%v + const_2
    END FUNCTION vector_add_const

    FUNCTION vector_int_add_const ( vec_1, const_2 )
        TYPE(vector_int)             :: vector_int_add_const
        TYPE(vector_int), INTENT(IN) :: vec_1
        INTEGER,          INTENT(IN) :: const_2

        vector_int_add_const%v = vec_1%v + const_2
    END FUNCTION vector_int_add_const

    FUNCTION const_add_vector ( const_1, vec_2 )
        TYPE(vector)              :: const_add_vector
        REAL(KIND=wp), INTENT(IN) :: const_1
        TYPE(vector),  INTENT(IN) :: vec_2

        const_add_vector%v = const_1 + vec_2%v
    END FUNCTION const_add_vector

    FUNCTION const_add_vector_int ( const_1, vec_2 )
        TYPE(vector_int)             :: const_add_vector_int
        REAL(KIND=wp),    INTENT(IN) :: const_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        const_add_vector_int%v = const_1 + vec_2%v
    END FUNCTION const_add_vector_int
 
!   operator(-)
    FUNCTION vector_int_minus ( vec_2_int )
        TYPE(vector_int)             :: vector_int_minus
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_int_minus%v = -vec_2_int%v
    END FUNCTION vector_int_minus

    FUNCTION vector_minus ( vec_2 )
        TYPE(vector)             :: vector_minus
        TYPE(vector), INTENT(IN) :: vec_2

        vector_minus%v = -vec_2%v
    END FUNCTION vector_minus

!   (*) operations
    FUNCTION vector_times_real ( vec_1, real_2 )
        TYPE(vector)              :: vector_times_real
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=sp), INTENT(IN) :: real_2

        vector_times_real%v = vec_1%v * real_2
    END FUNCTION vector_times_real

    FUNCTION vector_int_times_int ( vec_1_int, int_2 )
        TYPE(vector_int)             :: vector_int_times_int
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        INTEGER,          INTENT(IN) :: int_2

        vector_int_times_int%v = vec_1_int%v * int_2
    END FUNCTION vector_int_times_int

    FUNCTION vector_int_times_real ( vec_1_int, real_2 )
        TYPE(vector)                  :: vector_int_times_real
        TYPE(vector_int),  INTENT(IN) :: vec_1_int
        REAL(KIND=sp),     INTENT(IN) :: real_2

        vector_int_times_real%v = vec_1_int%v * real_2
    END FUNCTION vector_int_times_real

    FUNCTION real_times_vector_int ( real_1, vec_2_int )
        TYPE(vector)                 :: real_times_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int
        REAL(KIND=sp),    INTENT(IN) :: real_1

        real_times_vector_int%v = real_1 * vec_2_int%v
    END FUNCTION real_times_vector_int

    FUNCTION vector_int_times_wp ( vec_1_int, wp_2 )
        TYPE(vector)                 :: vector_int_times_wp
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        REAL(KIND=wp),    INTENT(IN) :: wp_2

        vector_int_times_wp%v = vec_1_int%v * wp_2
    END FUNCTION vector_int_times_wp

    FUNCTION wp_times_vector_int ( wp_1, vec_2_int )
        TYPE(vector)                 :: wp_times_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int
        REAL(KIND=wp),    INTENT(IN) :: wp_1

        wp_times_vector_int%v = wp_1 * vec_2_int%v
    END FUNCTION wp_times_vector_int
    
    FUNCTION real_times_vector ( real_1, vec_2 )
        TYPE(vector)              :: real_times_vector
        TYPE(vector),  INTENT(IN) :: vec_2
        REAL(KIND=sp), INTENT(IN) :: real_1

        real_times_vector%v = real_1 * vec_2%v
    END FUNCTION real_times_vector

    FUNCTION int_times_vector_int ( int_1, vec_2_int )
        TYPE(vector_int)             :: int_times_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int
        INTEGER,          INTENT(IN) :: int_1

        int_times_vector_int%v = int_1 * vec_2_int%v
    END FUNCTION int_times_vector_int

    FUNCTION vector_times_wp ( vec_1, wp_2 )
        TYPE(vector)              :: vector_times_wp
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

        vector_times_wp%v = vec_1%v * wp_2
    END FUNCTION vector_times_wp

    FUNCTION wp_times_vector ( wp_1, vec_2 )
        TYPE(vector)              :: wp_times_vector
        TYPE(vector),  INTENT(IN) :: vec_2
        REAL(KIND=wp), INTENT(IN) :: wp_1

        wp_times_vector%v = wp_1 * vec_2%v
    END FUNCTION wp_times_vector

    FUNCTION vector_times_int ( vec_1, int_2 )
        TYPE(vector)             :: vector_times_int
        TYPE(vector), INTENT(IN) :: vec_1
        INTEGER,      INTENT(IN) :: int_2

        vector_times_int%v = vec_1%v * int_2
    END FUNCTION vector_times_int

    FUNCTION int_times_vector ( int_1, vec_2 )
        TYPE(vector)             :: int_times_vector
        TYPE(vector), INTENT(IN) :: vec_2
        INTEGER,      INTENT(IN) :: int_1

        int_times_vector%v = int_1 * vec_2%v
    END FUNCTION int_times_vector

!   OPERATOR ( / ):    
    FUNCTION vector_div_real ( vec_1, real_2 )
        TYPE(vector)              :: vector_div_real
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=sp), INTENT(IN) :: real_2

        vector_div_real%v = vec_1%v / real_2
    END FUNCTION vector_div_real

    FUNCTION vector_div_wp ( vec_1, wp_2 )
        TYPE(vector)              :: vector_div_wp
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

        vector_div_wp%v = vec_1%v / wp_2
    END FUNCTION vector_div_wp

    FUNCTION vector_div_int ( vec_1, int_2 )
        TYPE(vector)             :: vector_div_int
        TYPE(vector), INTENT(IN) :: vec_1
        INTEGER,      INTENT(IN) :: int_2

        vector_div_int%v = vec_1%v / int_2
    END FUNCTION vector_div_int

    FUNCTION vector_int_div_int ( vec_1_int, int_2 )
        TYPE(vector)                 :: vector_int_div_int
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        INTEGER,          INTENT(IN) :: int_2

        vector_int_div_int%v = vec_1_int%v / REAL(int_2, KIND=wp)
    END FUNCTION vector_int_div_int

    FUNCTION vector_int_div_real ( vec_1_int, real_2 )
        TYPE(vector)                 :: vector_int_div_real
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        REAL(KIND=sp),    INTENT(IN) :: real_2

        vector_int_div_real%v = vec_1_int%v / real_2
    END FUNCTION vector_int_div_real

    FUNCTION vector_int_div_wp ( vec_1_int, wp_2 )
        TYPE(vector)                 :: vector_int_div_wp
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        REAL(KIND=wp),    INTENT(IN) :: wp_2

        IF ( wp_2 /= 0.0_wp ) THEN
            vector_int_div_wp%v = vec_1_int%v / wp_2
        ELSE
            vector_int_div_wp%v = 0.0_wp
            CALL die('Programming error. Division by zero wp_2 = 0.', 'vector_int_div_wp')
        ENDIF
    END FUNCTION vector_int_div_wp

    FUNCTION wp_div_vector_int ( wp_1, vec_2_int )
        TYPE(vector)                 :: wp_div_vector_int
        REAL(KIND=wp),    INTENT(IN) :: wp_1
        TYPE(vector_int), INTENT(IN) :: vec_2_int

!       No checks for zero denominator yet:
        wp_div_vector_int%v = wp_1 / vec_2_int%v
    END FUNCTION wp_div_vector_int

    FUNCTION real_div_vector_int ( real_1, vec_2_int )
        TYPE(vector)                 :: real_div_vector_int
        REAL(KIND=sp),    INTENT(IN) :: real_1
        TYPE(vector_int), INTENT(IN) :: vec_2_int

!       No checks for zero denominator yet:
        real_div_vector_int%v = real_1 / REAL(vec_2_int%v, KIND=wp)
    END FUNCTION real_div_vector_int

!   Continue with OPERATOR (/):
    FUNCTION vector_over_vector ( vec_1, vec_2 )
        TYPE(vector)             :: vector_over_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_over_vector%v = vec_1%v / vec_2%v
    END FUNCTION vector_over_vector

    FUNCTION vector_int_over_vector ( vec_1_int, vec_2 )
        TYPE(vector)                 :: vector_int_over_vector
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector),     INTENT(IN) :: vec_2

        vector_int_over_vector%v = vec_1_int%v / vec_2%v
    END FUNCTION vector_int_over_vector

    FUNCTION vector_over_vector_int ( vec_1, vec_2_int )
        TYPE(vector)                 :: vector_over_vector_int
        TYPE(vector),     INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_over_vector_int%v = vec_1%v / vec_2_int%v
    END FUNCTION vector_over_vector_int


!   This is special and probably not very useful:
    FUNCTION vector_int_over_vector_int ( vec_1_int, vec_2_int )
!
!       Purpose:
!       =======
!       Divides integer vector by integer vector
!
!       Date:            Programmer:                Description of changes:
!       ====             ==========                 ======================
!       Nov 2005         Strokopytov B.             Original code
!
!
        TYPE(vector_int)             :: vector_int_over_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1_int
        TYPE(vector_int), INTENT(IN) :: vec_2_int

        vector_int_over_vector_int%v =  vec_1_int%v /  vec_2_int%v
    END FUNCTION vector_int_over_vector_int

    FUNCTION wp_div_vector ( wp_1, vec_2 )
        TYPE(vector)              :: wp_div_vector
        REAL(KIND=wp), INTENT(IN) :: wp_1
        TYPE(vector),  INTENT(IN) :: vec_2

!       No checks for zero denominator yet:
        wp_div_vector%v = wp_1 / vec_2%v
    END FUNCTION wp_div_vector

    FUNCTION real_div_vector ( real_1, vec_2 )
        TYPE(vector)              :: real_div_vector
        REAL(KIND=sp), INTENT(IN) :: real_1
        TYPE(vector),  INTENT(IN) :: vec_2

!       No checks for zero denominator yet:
        real_div_vector%v = real_1 / vec_2%v
    END FUNCTION real_div_vector

    FUNCTION matrix_div_int ( mat_1, int_2 )
        TYPE(matrix)             :: matrix_div_int
        TYPE(matrix), INTENT(IN) :: mat_1
        INTEGER,      INTENT(IN) :: int_2

!       No checks for zero denominator yet:
        matrix_div_int%a = mat_1%a / int_2
    END FUNCTION matrix_div_int

    FUNCTION matrix_div_real ( mat_1, real_2 )
        TYPE(matrix)              :: matrix_div_real
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=sp), INTENT(IN) :: real_2

!       No checks for zero denominator yet:
        matrix_div_real%a = mat_1%a / real_2
    END FUNCTION matrix_div_real

    FUNCTION matrix_div_wp ( mat_1, wp_2 )
        TYPE(matrix)              :: matrix_div_wp
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=wp), INTENT(IN) :: wp_2

!       No checks for zero denominator yet:
        matrix_div_wp%a = mat_1%a / wp_2
    END FUNCTION matrix_div_wp

    FUNCTION matrix_div_matrix ( mat_1, mat_2 )
         TYPE(matrix)             :: matrix_div_matrix
         TYPE(matrix), INTENT(IN) :: mat_1
         TYPE(matrix), INTENT(IN) :: mat_2

!        No checks whether matrix inversion exists yet:
         matrix_div_matrix = mat_1 * .INV. mat_2
    END FUNCTION matrix_div_matrix

    FUNCTION matrix_int_div_int ( mat_1_int, int_2 )
        TYPE(matrix)                 :: matrix_int_div_int
        TYPE(matrix_int), INTENT(IN) :: mat_1_int
        INTEGER,          INTENT(IN) :: int_2

!       No checks for zero denominator yet:
        matrix_int_div_int%a = mat_1_int%a / REAL ( int_2, KIND=wp )       
    END FUNCTION matrix_int_div_int

    FUNCTION matrix_int_div_real ( mat_1_int, real_2 )
        TYPE(matrix)                 :: matrix_int_div_real
        TYPE(matrix_int), INTENT(IN) :: mat_1_int
        REAL(KIND=sp),    INTENT(IN) :: real_2

!       No checks for zero denominator yet:
        matrix_int_div_real%a = mat_1_int%a / real_2
    END FUNCTION matrix_int_div_real

    FUNCTION matrix_int_div_wp ( mat_1_int, wp_2 )
        TYPE(matrix)                 :: matrix_int_div_wp
        TYPE(matrix_int), INTENT(IN) :: mat_1_int
        REAL(KIND=wp),    INTENT(IN) :: wp_2

!       No checks for zero denominator yet:
        matrix_int_div_wp%a = mat_1_int%a / wp_2
    END FUNCTION matrix_int_div_wp

    FUNCTION int_div_matrix_int ( int_1, mat_2 )
        TYPE(matrix)                 :: int_div_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_2
        INTEGER,          INTENT(IN) :: int_1

        int_div_matrix_int = int_1 * .INV. mat_2
    END FUNCTION int_div_matrix_int

    FUNCTION real_div_matrix_int ( real_1, mat_2 )
        TYPE(matrix)                 :: real_div_matrix_int
        REAL(KIND=sp),    INTENT(IN) :: real_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        real_div_matrix_int = real_1 * .INV. mat_2
    END FUNCTION real_div_matrix_int

    FUNCTION wp_div_matrix_int ( wp_1, mat_2 )
        TYPE(matrix)                 :: wp_div_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_2
        REAL(KIND=wp),    INTENT(IN) :: wp_1

        wp_div_matrix_int =  wp_1 * .INV. mat_2
    END FUNCTION wp_div_matrix_int

    FUNCTION matrix_int_div_matrix_int ( mat_1, mat_2 )
        TYPE(matrix_int)             :: matrix_int_div_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_div_matrix_int = mat_1 * .INV. mat_2
        
    END FUNCTION matrix_int_div_matrix_int

    FUNCTION int_div_matrix ( int_1, mat_2 )
        TYPE(matrix)              :: int_div_matrix
        INTEGER,       INTENT(IN) :: int_1
        TYPE (matrix), INTENT(IN) :: mat_2

        int_div_matrix = int_1 * .INV. mat_2
    END FUNCTION int_div_matrix
 
    FUNCTION real_div_matrix ( real_1, mat_2 )
        TYPE(matrix)              :: real_div_matrix
        REAL(KIND=sp), INTENT(IN) :: real_1
        TYPE (matrix), INTENT(IN) :: mat_2

        real_div_matrix =  real_1 * .INV. mat_2
    END FUNCTION real_div_matrix

    FUNCTION wp_div_matrix ( wp_1, mat_2 )
         TYPE(matrix)              :: wp_div_matrix
         REAL(KIND=wp), INTENT(IN) :: wp_1
         TYPE(matrix),  INTENT(IN) :: mat_2

         wp_div_matrix = wp_1 * .INV. mat_2
    END FUNCTION wp_div_matrix

    FUNCTION dot_product_wp ( vec_1, vec_2 ) 
        REAL(KIND=wp)            :: dot_product_wp
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        dot_product_wp = DOT_PRODUCT ( vec_1%v, vec_2%v )
    END FUNCTION dot_product_wp

    FUNCTION dot_product_int ( vec_1, vec_2 )
        INTEGER                      :: dot_product_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        dot_product_int = DOT_PRODUCT ( vec_1%v, vec_2%v )
    END FUNCTION dot_product_int

    PURE FUNCTION dot_product_int_wp ( vec_1, vec_2 )
        REAL(KIND=wp)                :: dot_product_int_wp
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector),     INTENT(IN) :: vec_2

        dot_product_int_wp =  DOT_PRODUCT( vec_1%v, vec_2%v )
    END FUNCTION dot_product_int_wp

    PURE FUNCTION dot_product_wp_int ( vec_1, vec_2 )
        REAL(KIND=wp)                :: dot_product_wp_int
        TYPE(vector),     INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        dot_product_wp_int = DOT_PRODUCT ( vec_1%v, vec_2%v ) 
    END FUNCTION dot_product_wp_int

!   Special array operation(s):
    FUNCTION vec_int_2d_dot_product_vec ( vec_1, vec_2 ) RESULT ( new )
        REAL(KIND=wp),    DIMENSION(:,:), ALLOCATABLE             :: new
        TYPE(vector_int), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: vec_1
        TYPE(vector),                                  INTENT(IN) :: vec_2
!       Local variables:
        INTEGER                                                   :: i
        INTEGER                                                   :: j

        ALLOCATE ( new(SIZE(vec_1,DIM=1),SIZE(vec_1,DIM=2)) )
        DO i = 1, SIZE ( vec_1, DIM=1 )
            DO j = 1, SIZE ( vec_1, DIM=2 )
                new(i,j) = DOT_PRODUCT ( vec_1(i,j)%v, vec_2%v )
            ENDDO
        ENDDO
     END FUNCTION vec_int_2d_dot_product_vec

    FUNCTION cross_product ( vec_1, vec_2 )
        TYPE(vector)             :: cross_product
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        cross_product%v(1) = vec_1%v(2)*vec_2%v(3) - vec_1%v(3)*vec_2%v(2)
        cross_product%v(2) = vec_1%v(3)*vec_2%v(1) - vec_1%v(1)*vec_2%v(3)
        cross_product%v(3) = vec_1%v(1)*vec_2%v(2) - vec_1%v(2)*vec_2%v(1)
    END FUNCTION cross_product

    FUNCTION matrix_int_times_int ( mat_1, int_2 )
        TYPE(matrix_int)             :: matrix_int_times_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        INTEGER,          INTENT(IN) :: int_2

        matrix_int_times_int%a = mat_1%a * int_2
    END FUNCTION matrix_int_times_int

    FUNCTION matrix_int_times_real ( mat_1 , real_2 )
        TYPE(matrix)                 :: matrix_int_times_real
        TYPE(matrix_int), INTENT(IN) :: mat_1
        REAL(KIND=sp),    INTENT(IN) :: real_2

        matrix_int_times_real%a = mat_1%a * real_2
    END FUNCTION matrix_int_times_real

    FUNCTION int_times_matrix_int ( int_1, mat_2 )
        TYPE(matrix_int)             :: int_times_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_2
        INTEGER,          INTENT(IN) :: int_1

        int_times_matrix_int%a = int_1 * mat_2%a
    END FUNCTION int_times_matrix_int

    FUNCTION real_times_matrix_int ( real_1, mat_2 )
        TYPE(matrix)                 :: real_times_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_2
        REAL(KIND=sp),    INTENT(IN) :: real_1

        real_times_matrix_int%a = real_1 * mat_2%a
    END FUNCTION real_times_matrix_int

    FUNCTION wp_times_matrix_int ( wp_1, mat_2 )
        TYPE(matrix)                 :: wp_times_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_2
        REAL(KIND=wp),    INTENT(IN) :: wp_1

        wp_times_matrix_int%a = wp_1 * mat_2%a
    END FUNCTION wp_times_matrix_int

    FUNCTION matrix_int_times_wp ( mat_1 , wp_2 )
        TYPE(matrix)                  :: matrix_int_times_wp
        TYPE(matrix_int),  INTENT(IN) :: mat_1
        REAL(KIND=wp),     INTENT(IN) :: wp_2

        matrix_int_times_wp%a = mat_1%a * wp_2
    END FUNCTION matrix_int_times_wp

    FUNCTION matrix_int_times_matrix_int ( mat_1, mat_2 )
        TYPE(matrix_int)             :: matrix_int_times_matrix_int
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_int_times_matrix_int%a = MATMUL(mat_1%a, mat_2%a)
    END FUNCTION matrix_int_times_matrix_int     

    FUNCTION matrix_int_times_vector_int ( mat_1, vec_2 )
        TYPE(vector_int)              :: matrix_int_times_vector_int
        TYPE(matrix_int), INTENT(IN)  :: mat_1
        TYPE(vector_int), INTENT(IN)  :: vec_2

        matrix_int_times_vector_int%v = MATMUL ( mat_1%a, vec_2%v )
    END FUNCTION matrix_int_times_vector_int

    FUNCTION matrix_times_real ( mat_1, real_2 )
        TYPE (matrix)               :: matrix_times_real
        TYPE (matrix),   INTENT(IN) :: mat_1
        REAL (KIND=sp),  INTENT(IN) :: real_2

        matrix_times_real%a = mat_1%a * real_2
    END FUNCTION matrix_times_real

    FUNCTION matrix_times_wp ( mat_1, wp_2 )
        TYPE(matrix)               :: matrix_times_wp
        TYPE(matrix),   INTENT(IN) :: mat_1
        REAL(KIND=wp),  INTENT(IN) :: wp_2

        matrix_times_wp%a = mat_1%a * wp_2
    END FUNCTION matrix_times_wp

    FUNCTION matrix_times_int ( mat_1, int_2 )
        TYPE(matrix)             :: matrix_times_int
        TYPE(matrix), INTENT(IN) :: mat_1
        INTEGER,      INTENT(IN) :: int_2

        matrix_times_int%a = mat_1%a * int_2
    END FUNCTION matrix_times_int

    FUNCTION real_times_matrix ( real_1, mat_2 )
        TYPE(matrix)              :: real_times_matrix
        REAL(KIND=sp), INTENT(IN) :: real_1
        TYPE(matrix),  INTENT(IN) :: mat_2

        real_times_matrix%a = real_1 * mat_2%a
    END FUNCTION real_times_matrix    

!   This should be renamed:
    FUNCTION wp_times_matrix ( wp_1, mat_2 )
        TYPE(matrix)               :: wp_times_matrix
        REAL(KIND=wp),  INTENT(IN) :: wp_1
        TYPE(matrix),   INTENT(IN) :: mat_2

        wp_times_matrix%a = wp_1 * mat_2%a
    END FUNCTION wp_times_matrix

    FUNCTION int_times_matrix ( int_1, mat_2 )
        TYPE(matrix)             :: int_times_matrix
        INTEGER,      INTENT(IN) :: int_1
        TYPE(matrix), INTENT(IN) :: mat_2

        int_times_matrix%a = int_1 * mat_2%a
    END FUNCTION int_times_matrix

    FUNCTION matrix_times_matrix ( mat_1, mat_2 )
        TYPE(matrix)             :: matrix_times_matrix
        TYPE(matrix), INTENT(IN) :: mat_1, mat_2

        matrix_times_matrix%a = MATMUL ( mat_1%a, mat_2%a )
    END FUNCTION matrix_times_matrix

    FUNCTION matrix_int_times_matrix ( mat_1, mat_2 )
        TYPE(matrix)                 :: matrix_int_times_matrix
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(matrix),     INTENT(IN) :: mat_2
        matrix_int_times_matrix%a = MATMUL ( mat_1%a, mat_2%a )
    END FUNCTION matrix_int_times_matrix

    FUNCTION matrix_times_matrix_int ( mat_1, mat_2 )
        TYPE(matrix)                 :: matrix_times_matrix_int
        TYPE(matrix),     INTENT(IN) :: mat_1
        TYPE(matrix_int), INTENT(IN) :: mat_2

        matrix_times_matrix_int%a = MATMUL ( mat_1%a, mat_2%a )
    END FUNCTION matrix_times_matrix_int

    FUNCTION matrix_times_vector ( mat_1, vec_2 )
        TYPE(vector)             :: matrix_times_vector
        TYPE(matrix), INTENT(IN) :: mat_1
        TYPE(vector), INTENT(IN) :: vec_2

        matrix_times_vector%v = MATMUL ( mat_1%a, vec_2%v )
    END FUNCTION matrix_times_vector
   
    FUNCTION matrix_times_vector_int ( mat_1, vec_2 )
        TYPE(vector)                 :: matrix_times_vector_int
        TYPE(matrix),     INTENT(IN) :: mat_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        matrix_times_vector_int%v= MATMUL ( mat_1%a, vec_2%v )
    END FUNCTION matrix_times_vector_int


    FUNCTION matrix_int_times_vector ( mat_1, vec_2 )
        TYPE(vector)                 :: matrix_int_times_vector
        TYPE(matrix_int), INTENT(IN) :: mat_1
        TYPE(vector),     INTENT(IN) :: vec_2

        matrix_int_times_vector%v = MATMUL ( mat_1%a, vec_2%v )
    END FUNCTION matrix_int_times_vector

    FUNCTION matrix_int_times_array_int ( mat_1, arr_2 )
        TYPE(vector_int)                  :: matrix_int_times_array_int
        TYPE(matrix_int),      INTENT(IN) :: mat_1
        INTEGER, DIMENSION(:), INTENT(IN) :: arr_2

        matrix_int_times_array_int%v = MATMUL ( mat_1%a, arr_2 )
    END FUNCTION matrix_int_times_array_int

    FUNCTION matrix_int_times_array_wp ( mat_1, arr_2 )
        TYPE(vector)                            :: matrix_int_times_array_wp
        TYPE(matrix_int),            INTENT(IN) :: mat_1
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: arr_2

        matrix_int_times_array_wp%v = MATMUL ( mat_1%a, arr_2 )
    END FUNCTION matrix_int_times_array_wp

    FUNCTION matrix_times_array_wp ( mat_1, arr_2 )
        TYPE(vector)                            :: matrix_times_array_wp
        TYPE(matrix),                INTENT(IN) :: mat_1
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: arr_2

        matrix_times_array_wp%v = MATMUL ( mat_1%a, arr_2 ) 
    END FUNCTION matrix_times_array_wp

    FUNCTION transpose_int ( mat )
        TYPE(matrix_int)             :: transpose_int
        TYPE(matrix_int), INTENT(IN) :: mat

        transpose_int%a = TRANSPOSE ( mat%a )
    END FUNCTION transpose_int

    FUNCTION transpose_wp ( mat )
        TYPE(matrix)             :: transpose_wp
        TYPE(matrix), INTENT(IN) :: mat

        transpose_wp%a = TRANSPOSE ( mat%a )
    END FUNCTION transpose_wp

    FUNCTION det_int ( mat )
        INTEGER                      :: det_int
        TYPE(matrix_int), INTENT(IN) :: mat

        det_int = mat%a(1,1)*( mat%a(2,2)*mat%a(3,3) - mat%a(2,3)*mat%a(3,2) ) &
                + mat%a(1,2)*( mat%a(2,3)*mat%a(3,1) - mat%a(2,1)*mat%a(3,3) ) &
                + mat%a(1,3)*( mat%a(2,1)*mat%a(3,2) - mat%a(2,2)*mat%a(3,1) )
    END FUNCTION det_int

    FUNCTION det ( mat )
        REAL(KIND=wp)            :: det
        TYPE(matrix), INTENT(IN) :: mat

        det     = mat%a(1,1)*( mat%a(2,2)*mat%a(3,3) - mat%a(2,3)*mat%a(3,2) ) &
                + mat%a(1,2)*( mat%a(2,3)*mat%a(3,1) - mat%a(2,1)*mat%a(3,3) ) &
                + mat%a(1,3)*( mat%a(2,1)*mat%a(3,2) - mat%a(2,2)*mat%a(3,1) )
    END FUNCTION det

    FUNCTION transpose_of_the_inverse_int ( mat )

        TYPE(matrix_int)             :: transpose_of_the_inverse_int
        TYPE(matrix_int), INTENT(IN) :: mat
!       Local variables:
        INTEGER                      :: det_i
!
        det_i = det_int ( mat )

!       This has to do with unitary matrices used for symmetry operators:
        IF ( ABS ( det_i ) /= 1 ) THEN
            CALL die (' Invalid abs(determinant) /=1',' transpose_of_the_inverse...')
        ENDIF

!       Matrix inversion is impossible:
        IF ( det_i == 0 ) THEN
            CALL print_mat_int(mat)
            CALL die( 'zero determinant','transpose_of_the_inverse_int')
        ENDIF

!       Calculate the inverse:
        transpose_of_the_inverse_int%a(1,1) = ( mat%a(2,2)*mat%a(3,3) - mat%a(2,3)*mat%a(3,2) ) / det_i
        transpose_of_the_inverse_int%a(2,1) = ( mat%a(3,2)*mat%a(1,3) - mat%a(3,3)*mat%a(1,2) ) / det_i
        transpose_of_the_inverse_int%a(3,1) = ( mat%a(1,2)*mat%a(2,3) - mat%a(1,3)*mat%a(2,2) ) / det_i

        transpose_of_the_inverse_int%a(1,2) = ( mat%a(2,3)*mat%a(3,1) - mat%a(2,1)*mat%a(3,3) ) / det_i
        transpose_of_the_inverse_int%a(2,2) = ( mat%a(3,3)*mat%a(1,1) - mat%a(3,1)*mat%a(1,3) ) / det_i
        transpose_of_the_inverse_int%a(3,2) = ( mat%a(1,3)*mat%a(2,1) - mat%a(1,1)*mat%a(2,3) ) / det_i
       
        transpose_of_the_inverse_int%a(1,3) = ( mat%a(2,1)*mat%a(3,2) - mat%a(2,2)*mat%a(3,1) ) / det_i
        transpose_of_the_inverse_int%a(2,3) = ( mat%a(3,1)*mat%a(1,2) - mat%a(3,2)*mat%a(1,1) ) / det_i
        transpose_of_the_inverse_int%a(3,3) = ( mat%a(1,1)*mat%a(2,2) - mat%a(1,2)*mat%a(2,1) ) / det_i
    END FUNCTION transpose_of_the_inverse_int

    FUNCTION transpose_of_the_inverse ( mat )
        TYPE(matrix)             :: transpose_of_the_inverse
        TYPE(matrix), INTENT(IN) :: mat
!       Local variables:
        REAL(KIND=wp)            :: det_wp

!       Check whether inverse is possible:
        det_wp = det( mat )
        IF ( det_wp == 0.0_wp) CALL die ( ' Zero determinant','transpose_of_the_inverse')

        transpose_of_the_inverse%a(1,1) = ( mat%a(2,2)*mat%a(3,3) - mat%a(2,3)*mat%a(3,2) ) / det_wp
        transpose_of_the_inverse%a(2,1) = ( mat%a(3,2)*mat%a(1,3) - mat%a(3,3)*mat%a(1,2) ) / det_wp
        transpose_of_the_inverse%a(3,1) = ( mat%a(1,2)*mat%a(2,3) - mat%a(1,3)*mat%a(2,2) ) / det_wp

        transpose_of_the_inverse%a(1,2) = ( mat%a(2,3)*mat%a(3,1) - mat%a(2,1)*mat%a(3,3) ) / det_wp
        transpose_of_the_inverse%a(2,2) = ( mat%a(3,3)*mat%a(1,1) - mat%a(3,1)*mat%a(1,3) ) / det_wp
        transpose_of_the_inverse%a(3,2) = ( mat%a(1,3)*mat%a(2,1) - mat%a(1,1)*mat%a(2,3) ) / det_wp

        transpose_of_the_inverse%a(1,3) = ( mat%a(2,1)*mat%a(3,2) - mat%a(2,2)*mat%a(3,1) ) / det_wp
        transpose_of_the_inverse%a(2,3) = ( mat%a(3,1)*mat%a(1,2) - mat%a(3,2)*mat%a(1,1) ) / det_wp
        transpose_of_the_inverse%a(3,3) = ( mat%a(1,1)*mat%a(2,2) - mat%a(1,2)*mat%a(2,1) ) / det_wp

    END FUNCTION transpose_of_the_inverse

    FUNCTION inverse_int ( mat )
        TYPE(matrix_int)             :: inverse_int
        TYPE(matrix_int), INTENT(IN) :: mat

        inverse_int = transpose_int ( transpose_of_the_inverse_int ( mat ) )
    END FUNCTION inverse_int

    FUNCTION inverse ( mat )
        TYPE(matrix)             :: inverse
        TYPE(matrix), INTENT(IN) :: mat

        inverse = transpose_wp ( transpose_of_the_inverse ( mat ) )
    END FUNCTION inverse

    FUNCTION vector_inverse ( vec, pow )
        TYPE(vector)             :: vector_inverse
        TYPE(vector), INTENT(IN) :: vec
        INTEGER,      INTENT(IN) :: pow

        IF (pow >= 0) THEN
            CALL die('Programming error. pow is supposed to be negative in function inverse...',&
                     'vector_inverse')
        ENDIF

        vector_inverse%v(1) = 1.0_wp / ( vec%v(1) ** ABS( pow ) )
        vector_inverse%v(2) = 1.0_wp / ( vec%v(2) ** ABS( pow ) )
        vector_inverse%v(3) = 1.0_wp / ( vec%v(3) ** ABS( pow ) )
    END FUNCTION vector_inverse

!   Printing utilities:
    SUBROUTINE print_mat ( m1, m2, m3, fmtin )
        TYPE(matrix),     INTENT(IN), OPTIONAL :: m1
        TYPE(matrix),     INTENT(IN), OPTIONAL :: m2
        TYPE(matrix),     INTENT(IN), OPTIONAL :: m3
!       Local variables:
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: fmtin
        CHARACTER(len=60)                      :: fmt
        INTEGER                                :: i
        INTEGER                                :: j

        IF ( LEN( fmtin ) > 60 ) THEN
            CALL die ( 'Length of format string exceeds 60 characters.', 'print_mat' )
        ENDIF

        IF ( PRESENT ( fmtin ) ) THEN
            fmt = fmtin
        ELSE
            IF ( m1%a(1,1) > 1.0_wp ) THEN
                fmt = '(" PRINT_MAT> ", 3F9.3, 4X, 3F10.6, 4X, 3F9.5)'
            ELSE
                fmt = '(" PRINT_MAT> ", 3(4X,3F12.6) )'
            ENDIF
        ENDIF

        IF ( PRESENT ( m1 ) .AND. PRESENT ( m2 ) .AND. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*, fmt) (m1%a(i, j), j = 1, 3), (m2%a(i, j), j = 1, 3), (m3%a(i, j), j = 1, 3)
            ENDDO
        ELSE IF ( PRESENT ( m1 ) .AND. PRESENT ( m2 ) .AND. .NOT. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*, fmt) (m1%a(i, j), j = 1, 3), (m2%a(i, j), j =1, 3)
            ENDDO
        ELSE IF ( PRESENT ( m1 ) .AND. .NOT. PRESENT ( m2 ) .AND. .NOT. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*, fmt) (m1%a(i, j), j = 1, 3)
            ENDDO
        ELSE
            CALL warn('Nothing to print.', 'print_mat')
        ENDIF

    END SUBROUTINE print_mat

    SUBROUTINE print_mat_int ( m1, m2, m3, fmtin )
        TYPE(matrix_int), INTENT(IN), OPTIONAL :: m1
        TYPE(matrix_int), INTENT(IN), OPTIONAL :: m2
        TYPE(matrix_int), INTENT(IN), OPTIONAL :: m3
!       Local variables:
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: fmtin
        CHARACTER(len=60)                      :: fmt
        INTEGER                                :: i
        INTEGER                                :: j

        IF ( LEN ( fmtin ) > 60 ) STOP ' Error in print_mat_int...'

        IF ( PRESENT ( fmtin ) ) THEN
            fmt = fmtin
        ELSE
            fmt = '(3(4X,3I4))'
        ENDIF

        IF ( PRESENT ( m1 ) .AND. PRESENT ( m2 ) .AND. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*,fmt) (m1%a(i,j),j=1,3), (m2%a(i,j),j=1,3), (m3%a(i,j),j=1,3)
            ENDDO
        ELSE IF ( PRESENT ( m1 ) .AND. PRESENT ( m2 ) .AND. .NOT. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*,fmt) (m1%a(i,j),j=1,3), (m2%a(i,j),j=1,3)
            ENDDO
        ELSE IF ( PRESENT ( m1 ) .AND. .NOT. PRESENT ( m2 ) .AND. .NOT. PRESENT ( m3 ) ) THEN
            DO i = 1, 3
                WRITE(*, fmt) (m1%a(i,j),j=1,3)
            ENDDO
        ENDIF
    END SUBROUTINE print_mat_int
   
    SUBROUTINE wp_array_to_tensor ( tensor_out, arr_in )
        TYPE(tensor),                  INTENT(OUT) :: tensor_out
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)  :: arr_in

        tensor_out%r11 = arr_in(1,1)*arr_in(1,1) + arr_in(2,1)*arr_in(2,1) + arr_in(3,1)*arr_in(3,1)
        tensor_out%r22 = arr_in(1,2)*arr_in(1,2) + arr_in(2,2)*arr_in(2,2) + arr_in(3,2)*arr_in(3,2)
        tensor_out%r33 = arr_in(1,3)*arr_in(1,3) + arr_in(2,3)*arr_in(2,3) + arr_in(3,3)*arr_in(3,3)

        tensor_out%r12 = 2.0 * (arr_in(1,1)*arr_in(1,2) + arr_in(2,1)*arr_in(2,2) + arr_in(3,1)*arr_in(3,2))
        tensor_out%r13 = 2.0 * (arr_in(1,1)*arr_in(1,3) + arr_in(2,1)*arr_in(2,3) + arr_in(3,1)*arr_in(3,3))
        tensor_out%r23 = 2.0 * (arr_in(1,2)*arr_in(1,3) + arr_in(2,2)*arr_in(2,3) + arr_in(3,2)*arr_in(3,3))
    END SUBROUTINE wp_array_to_tensor

!   .MMOD. interface:
    FUNCTION vec_int_mod_vec_int ( vec_1, vec_2 )
        TYPE(vector_int)             :: vec_int_mod_vec_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        vec_int_mod_vec_int%v(1) = MOD ( vec_1%v(1) + 10000*vec_2%v(1), vec_2%v(1) )
        vec_int_mod_vec_int%v(2) = MOD ( vec_1%v(2) + 10000*vec_2%v(2), vec_2%v(2) ) 
        vec_int_mod_vec_int%v(3) = MOD ( vec_1%v(3) + 10000*vec_2%v(3), vec_2%v(3) ) 
    END FUNCTION  vec_int_mod_vec_int

    FUNCTION vec_int_mod_arr_int ( vec_1, arr_2 )
        TYPE(vector_int)                  :: vec_int_mod_arr_int
        TYPE(vector_int),      INTENT(IN) :: vec_1
        INTEGER, DIMENSION(3), INTENT(IN) :: arr_2

        vec_int_mod_arr_int%v(1) = MOD ( vec_1%v(1) + 10000*arr_2(1), arr_2(1) )
        vec_int_mod_arr_int%v(2) = MOD ( vec_1%v(2) + 10000*arr_2(2), arr_2(2) )
        vec_int_mod_arr_int%v(3) = MOD ( vec_1%v(3) + 10000*arr_2(3), arr_2(3) )
    END FUNCTION  vec_int_mod_arr_int

!   These 2 routines require some simplification
    FUNCTION arr_int_mod_arr_int ( arr_1, arr_2 ) RESULT ( temp )
        INTEGER, DIMENSION(3)             :: temp
        INTEGER, DIMENSION(3), INTENT(IN) :: arr_1
        INTEGER, DIMENSION(3), INTENT(IN) :: arr_2

        temp(1) = MOD (arr_1(1) + 10000 * arr_2(1), arr_2(1))
        temp(2) = MOD (arr_1(2) + 10000 * arr_2(2), arr_2(2))
        temp(3) = MOD (arr_1(3) + 10000 * arr_2(3), arr_2(3))
    END FUNCTION  arr_int_mod_arr_int

    FUNCTION arr_int_mod_vec_int ( arr_1, vec_2 ) RESULT ( temp )
        INTEGER, DIMENSION(3)             :: temp
        INTEGER, DIMENSION(:), INTENT(IN) :: arr_1
        TYPE(vector_int),      INTENT(IN) :: vec_2

        temp(1) = MOD(arr_1(1) + 10000*vec_2%v(1), vec_2%v(1))
        temp(2) = MOD(arr_1(2) + 10000*vec_2%v(2), vec_2%v(2))
        temp(3) = MOD(arr_1(3) + 10000*vec_2%v(3), vec_2%v(3))
    END FUNCTION arr_int_mod_vec_int 

!   OPERATOR ( * ) :
    FUNCTION tensor_times_vector ( tensor_1, vec_2 )
        REAL(KIND=wp)              :: tensor_times_vector
        TYPE(tensor),   INTENT(IN) :: tensor_1
        TYPE(vector),   INTENT(IN) :: vec_2

        tensor_times_vector = tensor_1%r11 * vec_2%v(1) ** 2 +         &
                              tensor_1%r22 * vec_2%v(2) ** 2 +         &
                              tensor_1%r33 * vec_2%v(3) ** 2 +         &
                              tensor_1%r12 * vec_2%v(1) * vec_2%v(2) + &
                              tensor_1%r13 * vec_2%v(1) * vec_2%v(3) + &
                              tensor_1%r23 * vec_2%v(2) * vec_2%v(3)              
    END FUNCTION tensor_times_vector

    FUNCTION tensor_times_array ( tensor_1, arr_2 )
        REAL(KIND=wp)                           :: tensor_times_array
        TYPE(tensor),                INTENT(IN) :: tensor_1
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: arr_2

        tensor_times_array = tensor_1%r11 * arr_2(1) ** 2 +       &
                             tensor_1%r22 * arr_2(2) ** 2 +       &
                             tensor_1%r33 * arr_2(3) ** 2 +       &
                             tensor_1%r12 * arr_2(1) * arr_2(2) + &
                             tensor_1%r13 * arr_2(1) * arr_2(3) + &
                             tensor_1%r23 * arr_2(2) * arr_2(3)
    END FUNCTION tensor_times_array

    FUNCTION tensor_times_vector_int ( tensor_1, vec_2 )
        REAL(KIND=wp)                :: tensor_times_vector_int
        TYPE(tensor),     INTENT(IN) :: tensor_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        tensor_times_vector_int = tensor_1%r11 * vec_2%v(1) ** 2 +         &
                                  tensor_1%r22 * vec_2%v(2) ** 2 +         &
                                  tensor_1%r33 * vec_2%v(3) ** 2 +         &
                                  tensor_1%r12 * vec_2%v(1) * vec_2%v(2) + &
                                  tensor_1%r13 * vec_2%v(1) * vec_2%v(3) + &
                                  tensor_1%r23 * vec_2%v(2) * vec_2%v(3)
    END FUNCTION tensor_times_vector_int

    FUNCTION tensor_times_array_int ( tensor_1, arr_2 )
        REAL(KIND=wp)                     :: tensor_times_array_int
        TYPE(tensor),          INTENT(IN) :: tensor_1
        INTEGER, DIMENSION(:), INTENT(IN) :: arr_2

        tensor_times_array_int = tensor_1%r11 * arr_2(1) ** 2 +       &
                                 tensor_1%r22 * arr_2(2) ** 2 +       &
                                 tensor_1%r33 * arr_2(3) ** 2 +       &
                                 tensor_1%r12 * arr_2(1) * arr_2(2) + &
                                 tensor_1%r13 * arr_2(1) * arr_2(3) + &
                                 tensor_1%r23 * arr_2(2) * arr_2(3)
    END FUNCTION tensor_times_array_int

!   Print private components of tensor:    
    SUBROUTINE print_tensor ( t_1, t_2 ) 
        TYPE(tensor), OPTIONAL :: t_1
        TYPE(tensor), OPTIONAL :: t_2
!       Local vars:
        CHARACTER(len=120)     :: fmt_1, fmt_2, fmt

!       FIXME this is not general:
!       Need to redo all tensor stuff... BVS OCT 2007

        fmt_1 = '('' PRINT_TENSOR> Real space tensor:'',6X, 6F12.1,/,&
                 &'' PRINT_TENSOR> Reciprocal space tensor:'', 6F12.9)'
        fmt_2 = '('' PRINT_TENSOR> Reciprocal space tensor:'', 6F12.8,/,&
                 &'' PRINT_TENSOR> Real space tensor:'', 6X, 6F12.1)'

        WRITE(*,"(' PRINT_TENSOR> ',30X,'   r11         r22         r33         r12         r13         r23')")           

        IF ( PRESENT ( t_1 ) .AND. PRESENT ( t_2 ) ) THEN

            IF ( t_1%r11 > 1.0_wp ) THEN
                fmt = fmt_1
            ELSE
                fmt = fmt_2
            ENDIF
            
            WRITE(*, fmt) t_1%r11, t_1%r22, t_1%r33, t_1%r12, t_1%r13, t_1%r23,  &
                          t_2%r11, t_2%r22, t_2%r33, t_2%r12, t_2%r13, t_2%r23

        ELSE IF ( PRESENT ( t_1 ) .AND. .NOT. PRESENT ( t_2 ) ) THEN

            WRITE(*, fmt ) t_1%r11, t_1%r22, t_1%r33, t_1%r12, t_1%r13, t_1%r23

        ENDIF 

    END SUBROUTINE print_tensor

    FUNCTION extract_diagonal ( mat_1 )
        TYPE(vector)             :: extract_diagonal
        TYPE(matrix), INTENT(IN) :: mat_1

        extract_diagonal%v(1) = mat_1%a(1,1)
        extract_diagonal%v(2) = mat_1%a(2,2)
        extract_diagonal%v(3) = mat_1%a(3,3)       
    END FUNCTION extract_diagonal

    FUNCTION trace ( mat_1 )
        REAL(KIND=wp)            :: trace
        TYPE(matrix), INTENT(IN) :: mat_1

        trace = mat_1%a(1,1) + mat_1%a(2,2) + mat_1%a(3,3)
    END FUNCTION trace 

    FUNCTION h2mod ( vec_1, vec_2 ) 
        TYPE(vector_int)             :: h2mod
        TYPE(vector_int), INTENT(IN) :: vec_1 
        TYPE(vector_int), INTENT(IN) :: vec_2

        h2mod%v(1) = 2 * lmod ( vec_1%v(1), vec_2%v(1) )
        h2mod%v(2) =     lmod ( vec_1%v(2), vec_2%v(2) )
        h2mod%v(3) =     lmod ( vec_1%v(3), vec_2%v(3) )
    END FUNCTION h2mod

    FUNCTION hmod ( vec_1, vec_2 )
        TYPE(vector_int)             :: hmod
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2

        hmod%v(1) = lmod ( vec_1%v(1), vec_2%v(1) )
        hmod%v(2) = lmod ( vec_1%v(2), vec_2%v(2) )
        hmod%v(3) = lmod ( vec_1%v(3), vec_2%v(3) )
    END FUNCTION hmod

    FUNCTION lmod ( i, j )
!
!       Purpose:
!       =======
!       fixed integer mod (optimised for 0<=i<j)
!
!       Programmer:                      Date:
!       ==========                       ====
!       Lifted from K.Cowtan             Sep 2003
!       Adapted for F95 by B.Strokopytov Sep 2003
!
        INTEGER             :: lmod
        INTEGER, INTENT(IN) :: i
        INTEGER, INTENT(IN) :: j

        IF ( j == 0 ) THEN
            CALL die('Programming error. j = 0.', 'lmod')
        ENDIF

        IF ( i < 0 ) THEN
            lmod = MOD(i + 1, j) + j - 1
        ELSE IF ( i >= j ) THEN
            lmod = MOD(i, j)
        ELSE
            lmod = i
        ENDIF
    END FUNCTION lmod
  
    FUNCTION rmod ( r, s )
        REAL(KIND=wp)             :: rmod
        REAL(KIND=wp), INTENT(IN) :: r, s
!       Local
        REAL(KIND=wp)             :: t      
      
        t = MOD(r,s)
        IF ( t >= 0.0_wp) THEN
            rmod = t
        ELSE
            rmod = t + s
        ENDIF
    END FUNCTION rmod

    FUNCTION eigen_values ( mat_1 )
        TYPE(vector) :: eigen_values
        TYPE(matrix) :: mat_1
!       CALL mkl lapack95 interface (symmetric matrices only):
        CALL syevd  (mat_1%A, eigen_values%v, 'N', UPLO='U')
    END FUNCTION eigen_values

!   NINT interface:
    FUNCTION nint_vector ( vec_1 )
        TYPE(vector_int)         :: nint_vector
        TYPE(vector), INTENT(IN) :: vec_1

        nint_vector = NINT ( vec_1%v )
    END FUNCTION nint_vector

    FUNCTION nint_matrix ( mat_1 )
        TYPE(matrix_int)         :: nint_matrix
        TYPE(matrix), INTENT(IN) :: mat_1

        nint_matrix = NINT ( mat_1%a )
    END FUNCTION nint_matrix

!   INT interface:
    FUNCTION int_vector ( vec_1 )
        TYPE(vector_int)         :: int_vector
        TYPE(vector), INTENT(IN) :: vec_1

        int_vector = INT ( vec_1%v )
    END FUNCTION int_vector

    FUNCTION int_matrix ( mat_1 )
        TYPE(matrix_int)         :: int_matrix
        TYPE(matrix), INTENT(IN) :: mat_1

        int_matrix = INT ( mat_1%a )
    END FUNCTION int_matrix

!   REAL interface:
    FUNCTION real_wp_vector ( vec_1_int )
        TYPE(vector)     :: real_wp_vector
        TYPE(vector_int) :: vec_1_int

        real_wp_vector%v = REAL ( vec_1_int%v, KIND=wp )
    END FUNCTION real_wp_vector

    FUNCTION real_wp_matrix ( mat_1_int )
        TYPE(matrix)     :: real_wp_matrix
        TYPE(matrix_int) :: mat_1_int

        real_wp_matrix%a = REAL ( mat_1_int%a, KIND=wp )
    END FUNCTION real_wp_matrix

!   MOD interface:
    FUNCTION mod_vector ( vec_1, factor )
        TYPE(vector)              :: mod_vector
        TYPE(vector),  INTENT(IN) :: vec_1
        REAL(KIND=wp), INTENT(IN) :: factor
        mod_vector%v = MOD ( vec_1%v, factor )
    END FUNCTION mod_vector

    FUNCTION mod_vector_int ( vec_1, factor )
        TYPE(vector_int)             :: mod_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        INTEGER,          INTENT(IN) :: factor

        mod_vector_int%v = MOD ( vec_1%v, factor )
    END FUNCTION mod_vector_int

!   ABS interface:
    FUNCTION abs_matrix ( mat_1 )
        TYPE(matrix)             :: abs_matrix
        TYPE(matrix), INTENT(IN) :: mat_1
        
        abs_matrix = ABS ( mat_1%a )
    END FUNCTION abs_matrix

    FUNCTION abs_vector ( vec_1 )
        TYPE(vector)             :: abs_vector
        TYPE(vector), INTENT(IN) :: vec_1

        abs_vector = ABS ( vec_1%v )
    END FUNCTION abs_vector 

    FUNCTION matrix_gt_wp_const ( mat_1, wp_const )
        LOGICAL, DIMENSION(3,3)   :: matrix_gt_wp_const
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=wp), INTENT(IN) :: wp_const
      
        matrix_gt_wp_const = mat_1%a > wp_const
    END FUNCTION matrix_gt_wp_const

    FUNCTION vector_gt_wp_const ( vec_1, wp_const )
        LOGICAL,       DIMENSION(3) :: vector_gt_wp_const
        TYPE(vector),  INTENT(IN)   :: vec_1
        REAL(KIND=wp), INTENT(IN)   :: wp_const

        vector_gt_wp_const = vec_1%v > wp_const
    END FUNCTION vector_gt_wp_const

    FUNCTION matrix_lt_wp_const ( mat_1, wp_const )
        LOGICAL, DIMENSION(3,3)   :: matrix_lt_wp_const
        TYPE(matrix),  INTENT(IN) :: mat_1
        REAL(KIND=wp), INTENT(IN) :: wp_const

        matrix_lt_wp_const = mat_1%a < wp_const
    END FUNCTION matrix_lt_wp_const

    FUNCTION vector_lt_wp_const ( vec_1, wp_const )
        LOGICAL,       DIMENSION(3) :: vector_lt_wp_const
        TYPE(vector),  INTENT(IN)   :: vec_1
        REAL(KIND=wp), INTENT(IN)   :: wp_const

        vector_lt_wp_const = vec_1%v < wp_const
    END FUNCTION vector_lt_wp_const

    FUNCTION sqrt_vector ( vec_1 )
        TYPE(vector)     :: sqrt_vector
        TYPE(vector)     :: vec_1

        sqrt_vector = SQRT ( vec_1%v )
    END FUNCTION sqrt_vector

!   MINVAL interface:
    FUNCTION minval_vector ( vec_1 )
        REAL(KIND=wp)            :: minval_vector
        TYPE(vector), INTENT(IN) :: vec_1

        minval_vector = MINVAL (vec_1%v, DIM=1 )
    END FUNCTION minval_vector

    FUNCTION minval_matrix ( mat_1 )
        REAL(KIND=wp)            :: minval_matrix
        TYPE(matrix), INTENT(IN) :: mat_1 

        minval_matrix = MINVAL ( mat_1%a )
    END FUNCTION minval_matrix

!   MAXVAL interface:
    FUNCTION maxval_vector ( vec_1 )
        REAL(KIND=wp) :: maxval_vector
        TYPE(vector)  :: vec_1
    
        maxval_vector = MAXVAL ( vec_1%v, DIM=1 )
    END FUNCTION maxval_vector

    FUNCTION maxval_matrix ( mat_1 )
        REAL(KIND=wp) :: maxval_matrix
        TYPE(matrix)  :: mat_1

        maxval_matrix = MAXVAL ( mat_1%a )
    END FUNCTION maxval_matrix

!   OPERATOR (.x.)   
    FUNCTION vector_times_vector ( vec_1, vec_2 )
        TYPE(vector)             :: vector_times_vector
        TYPE(vector), INTENT(IN) :: vec_1
        TYPE(vector), INTENT(IN) :: vec_2

        vector_times_vector%v = vec_1%v * vec_2%v
    END FUNCTION vector_times_vector

    FUNCTION vector_int_times_vector ( vec_1, vec_2 )
        TYPE(vector)                 :: vector_int_times_vector
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector),     INTENT(IN) :: vec_2
        
        vector_int_times_vector%v = vec_1%v * vec_2%v
    END FUNCTION vector_int_times_vector

    FUNCTION vector_times_vector_int ( vec_1, vec_2 )
        TYPE(vector)                 :: vector_times_vector_int
        TYPE(vector),     INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2
        
        vector_times_vector_int%v = vec_1%v * vec_2%v
    END FUNCTION vector_times_vector_int

    FUNCTION vector_int_times_vector_int ( vec_1, vec_2 )
        TYPE(vector_int)             :: vector_int_times_vector_int
        TYPE(vector_int), INTENT(IN) :: vec_1
        TYPE(vector_int), INTENT(IN) :: vec_2
    
        vector_int_times_vector_int%v = vec_1%v * vec_2%v
    END FUNCTION vector_int_times_vector_int

END MODULE vectors
