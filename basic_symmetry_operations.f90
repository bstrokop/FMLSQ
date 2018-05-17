MODULE basic_symmetry_operations 
!
!   Purpose:
!   =======
!   Define a derived type called symop/symop_int and the
!   operations that can be defined on it. The module
!   defines eight operations that can be performed on symmetry
!   operators:
!
!
!
!    Operation                                           Operator
!    =========                                           ========
!
!    1. Creation from a double precision array                  =
!    2. Conversion to a double precision array                  =
!    3. Logical operators                                      ==,/=
!    3. Symmetry operator addition                              +
!    4. Symmetry operator subtraction                           -
!    5. Symmetry operator multiplications                       *
!
!
!
!    Record of revisions:
!    ===================
!
!    Date           Programmer              Description of change
!    ====           ==========              =====================
!    Sep 2003       Strokopytov B.          Original code
!    Sep 2005       Strokopytov B.          "PRIVATE" line added
!    Oct 2007       B.Strokopytov           eliminate_duplicate_symops added 
!
!
USE constants
!USE iso_varying_string
USE fail
USE parser_library 
USE select_kinds
USE string_manip
USE util
USE vectors

IMPLICIT NONE

! Declare real symop data types:
TYPE :: symop
    PRIVATE
    TYPE(matrix)                  :: a
    TYPE(vector)                  :: v
    INTEGER                       :: tsym_factor = 0
END TYPE

! Declare integer symop data types:
TYPE :: symop_int
    PRIVATE
    TYPE (matrix_int)             :: a
    TYPE (vector_int)             :: v
    INTEGER                       :: tsym_factor = 0
END TYPE

! Declare all items to be private
! except for type vector, matrix and
! the operators defined for it:

PRIVATE
PUBLIC :: ASSIGNMENT ( = )
PUBLIC :: check_duplicate_symops
PUBLIC :: check_symop_times_symop_matrix
PUBLIC :: eliminate_duplicate_symops
PUBLIC :: OPERATOR ( + )
PUBLIC :: OPERATOR ( - )
PUBLIC :: OPERATOR ( * )
PUBLIC :: OPERATOR ( == )
PUBLIC :: OPERATOR ( /= )
PUBLIC :: OPERATOR ( .DOT. )
PUBLIC :: OPERATOR ( .MODU. )
PUBLIC :: OPERATOR ( .SYMA. )
PUBLIC :: OPERATOR ( .SYMV. )
PUBLIC :: OPERATOR ( .INV. )
PUBLIC :: print_symop
PUBLIC :: set_tsym_factor
PUBLIC :: symop
PUBLIC :: symop_int

! Declare interface operators:


INTERFACE ASSIGNMENT ( = )

!   Char/int operations:
    MODULE PROCEDURE char_string_to_symop_int
    MODULE PROCEDURE symop_int_to_char_string
    MODULE PROCEDURE matrix_int_to_char_string

!   Integer operations:
    MODULE PROCEDURE symop_int_to_array_int
    MODULE PROCEDURE array_int_to_symop_int
    MODULE PROCEDURE set_symop_int
    MODULE PROCEDURE matrix_part_int

!   Double precision operations:
    MODULE PROCEDURE set_symop
    MODULE PROCEDURE matrix_part
    MODULE PROCEDURE vector_part 
    MODULE PROCEDURE wp_array_to_symop

!   Mixed operations:
    MODULE PROCEDURE symop_int_to_symop_double
    MODULE PROCEDURE symop_double_to_symop_int
    MODULE PROCEDURE symop_int_to_array_real
    MODULE PROCEDURE vector_part_int
END INTERFACE

INTERFACE OPERATOR ( .SYMA. )
    MODULE PROCEDURE set_symop_matrix_part
    MODULE PROCEDURE set_symop_int_matrix_part
    MODULE PROCEDURE set_symop_matrix_part_mat

    MODULE PROCEDURE extract_symop_int_matrix_part
    MODULE PROCEDURE extract_symop_matrix_part
END INTERFACE

INTERFACE OPERATOR ( .SYMV. )
    MODULE PROCEDURE set_symop_vector_part_vec
    MODULE PROCEDURE set_symop_vector_part_arr
    MODULE PROCEDURE set_symop_int_vector_part_vec
    MODULE PROCEDURE set_symop_int_vector_part_arr
   
    MODULE PROCEDURE extract_symop_int_vector_part
    MODULE PROCEDURE extract_symop_vector_part
END INTERFACE

INTERFACE OPERATOR ( == )
    MODULE PROCEDURE symop_int_eq_symop_int
    MODULE PROCEDURE symop_eq_symop
    MODULE PROCEDURE symop_int_eq_unity
    MODULE PROCEDURE symop_eq_unity
END INTERFACE

INTERFACE OPERATOR ( /= )
    MODULE PROCEDURE symop_int_ne_symop_int
    MODULE PROCEDURE symop_ne_symop
    MODULE PROCEDURE symop_int_ne_unity
    MODULE PROCEDURE symop_ne_unity
END INTERFACE

INTERFACE OPERATOR ( + )
    MODULE PROCEDURE symop_int_add_symop_int
    MODULE PROCEDURE symop_add_symop
    MODULE PROCEDURE symop_add_vector_int
    MODULE PROCEDURE symop_add_vector
END INTERFACE

INTERFACE OPERATOR ( - )
    MODULE PROCEDURE symop_int_subtract_symop_int
    MODULE PROCEDURE symop_subtract_symop
    MODULE PROCEDURE symop_subtract_vector_int
    MODULE PROCEDURE symop_subtract_vector
END INTERFACE

INTERFACE OPERATOR ( * )
    MODULE PROCEDURE symop_int_times_symop_int
    MODULE PROCEDURE symop_times_symop
    MODULE PROCEDURE matrix_int_times_symop_int
    MODULE PROCEDURE matrix_int_times_symop
    MODULE PROCEDURE matrix_times_symop
    MODULE PROCEDURE symop_times_matrix
    MODULE PROCEDURE symop_int_times_vector_int
    MODULE PROCEDURE symop_int_times_vector
    MODULE PROCEDURE symop_int_times_array_int
    MODULE PROCEDURE symop_times_vector
    MODULE PROCEDURE symop_int_times_array_double
    MODULE PROCEDURE string_times_string
END INTERFACE

INTERFACE OPERATOR ( .INV. )
    MODULE PROCEDURE invert_symop
    MODULE PROCEDURE invert_symop_int
END INTERFACE

INTERFACE OPERATOR ( .MODU. )
    MODULE PROCEDURE mod_tsym_int
    MODULE PROCEDURE mod_tsym
END INTERFACE

INTERFACE OPERATOR ( .DOT. )
    MODULE PROCEDURE vector_int_dot_symop
    MODULE PROCEDURE symop_dot_vector_int
END INTERFACE

INTERFACE print_symop
    MODULE PROCEDURE print_symop_int
    MODULE PROCEDURE print_symop_double
END INTERFACE

! Now define the implementing functions:
CONTAINS
    SUBROUTINE convert_symop ( line, sym )
!

!       List of parameters:
        CHARACTER(LEN=*),                     INTENT(INOUT) :: line
        INTEGER,              DIMENSION(3,4), INTENT(OUT)   :: sym
!       Declare local variables:
        CHARACTER(LEN=wlen),  DIMENSION(:),   ALLOCATABLE   :: words                                
        CHARACTER(LEN=LEN(line))                            :: string 
        INTEGER, DIMENSION(4)                               :: row
!       Counters:
        INTEGER                                             :: k

!       Make copy:
        string = line
!
!       Eliminate blanks:
        CALL elimbl ( string )
        IF ( debug > 25 ) WRITE(*,*) 'convert_symop>> ', string
!
!       Convert line to uppercase:
        CALL ucase ( string )
        IF ( debug > 25 ) WRITE(*,*) 'convert_symop>> ', string
!
        string = replace_commas ( string )
        CALL mysplit ( string, words )
        IF ( SIZE ( words ) /= 3 ) THEN
            CALL die (' Syntax error or not 3D operator.', 'convert_symop' )
        ENDIF
!
!       Get final matrix:
        DO k = 1, 3
            CALL getrow ( words(k), row, tsym_factor )
            sym(k, 1:4) = row(1:4)
        ENDDO

!       Deallocate storage:
        IF ( ALLOCATED ( words ) ) CALL deallocate_array ( words )
    END SUBROUTINE convert_symop

    FUNCTION replace_commas ( line ) RESULT ( new )
        CHARACTER(LEN=*)              :: line
        CHARACTER(LEN=LEN_TRIM(line)) :: new
!       Counters
        INTEGER                       :: i
        INTEGER                       :: commas

        commas = 0
        DO i = 1, LEN_TRIM ( line )
            IF ( line(i:i) /= ',' ) THEN
                new(i:i) = line(i:i)
            ELSE
                new(i:i) = ' '
                commas   = commas + 1
            ENDIF
        ENDDO
        IF ( commas < 2 ) THEN
            WRITE(*,*) line
            WRITE(*,*) commas
            CALL die ( 'Syntax error. Number of commas < 2', 'replace_commas' )
        ENDIF
    END FUNCTION replace_commas

    SUBROUTINE getrow ( string, sym, tsym_factor )
!
!       Purpose:
!       =======
!       converts Internatiol Tables style symmetry
!       operator to a row of integer matrix.
!       Factorizes translational part.
!
!       Date      Programmer       Description of change
!       ====      =========        =====================
!       08/08/03  Strokopytov B.   Original code
!
!

!       Declare calling parameters:
        CHARACTER(LEN=*),      INTENT(IN)       :: string
        INTEGER, DIMENSION(4), INTENT(OUT)      :: sym
        INTEGER,               INTENT(IN)       :: tsym_factor
!
!       Declare local variables:
        INTEGER                                 :: factor
        INTEGER                                 :: is
        CHARACTER(len=LEN(string))              :: temp
        CHARACTER(len=2), DIMENSION(6), SAVE    :: symbols  = (/'-X', '-Y', '-Z', '+X', '+Y', '+Z'/)
        CHARACTER(len=1), DIMENSION(3), SAVE    :: symbols1 = (/'X','Y','Z'/)
        INTEGER,          DIMENSION(6), SAVE    :: vector   = (/-1, -1, -1, 1, 1, 1/)
!       Counters:
        INTEGER                                 :: i
        INTEGER                                 :: j

!       Initialize sym:
        sym  = (/0, 0, 0, 0/)
!       Copy string:
        temp = string
!
!       Translate row from matrix part and clean up interpreted chars:
        row: IF ( LEN(temp) > 1 ) THEN
            loop: DO i = 1, 6
                j = INDEX( temp, symbols(i) )

                inner_if: IF ( j > 0 ) THEN
                    is = i
                    if ( is > 3 ) is = is - 3
                    sym(is)     = vector(i)
                    temp(j:j+1) = ' '
                END IF inner_if

            END DO loop
        ENDIF row

        new_loop: DO i = 1, 3
            symchk:  IF ( sym(i) == 0 ) THEN
                j = INDEX ( temp, symbols1(i) )
                IF ( j > 0 ) THEN
                    sym(i)    = 1
                    temp(j:j) = ' '
                END IF
            ENDIF symchk
        END DO new_loop

!       Do translational part:
        j = INDEX( temp, '/')

        trans: IF ( j > 0 ) THEN

!           Checks:
            IF ( j < 2 ) THEN
                CALL messag(string, 'getrow')
                CALL die('Syntax error 1 in translational part of a symmetry operator.', 'getrow')
            ENDIF
            
            IF ( j == LEN( temp ) ) THEN
                CALL messag(string, 'getrow')
                CALL die ('Syntax error 2 in translation part of a symmetry operator.',  'getrow')
            ENDIF

!           Get numerator and clean up:
            IF ( j == 2 ) THEN
                READ ( temp(j-1:j-1), * ) sym(4)
                temp(j-1:j-1) = ' '
            ELSE IF ( j > 2 .AND. INDEX( '-+', temp(j-2:j-2) ) > 0 ) THEN
                READ ( temp(j-2:j-1), * ) sym(4)
                temp(j-2:j-1) = ' '
            ELSE
                CALL messag(string, 'getrow')
                CALL die ('Syntax error 3 in translation part of a symmetry operator.','getrow')
            END IF

!           Get denumerator:
            READ ( temp(j+1:j+1), * ) factor

!           Clean up denumerator:
            temp(j+1:j+1) = ' '

!           Clear "/":
            temp(j:j) = ' '

!           Factorize to form x/24 (e.g. 3/4 becomes 9/12):
            factor = tsym_factor/factor
            sym(4) = factor*sym(4)
            sym(4) = MOD( sym(4) + 100*tsym_factor, tsym_factor )
        ENDIF trans

!       Check whether any non-blank characters left:
        DO i = 1, LEN_TRIM ( temp )
            nonblank: IF ( temp(i:i) /= ' ' ) THEN
                CALL die ('The following uninterpreted characters left: '//temp,'getrow')
            END IF nonblank
        END DO

    END SUBROUTINE getrow

    SUBROUTINE char_string_to_symop_int ( symop_result, line )
!
!
!       Purpose:
!       =======
!       Converts character string, e.g. 'X,Z,Y'
!       to integer symmetry operator (type symop_int)
!
!       Date:           Programmer:              Description of changes:
!       ====            ==========               ======================
!       Oct 2003        Strokopytov B.           Original code
!
!
!
!
        CHARACTER(len=*), INTENT(IN)  :: line
        TYPE (symop_int), INTENT(OUT) :: symop_result
!       Local variables: 
        CHARACTER(LEN=LEN(line))      :: temp
        INTEGER, DIMENSION(3,4)       :: sym

        temp = line

!       Convert temp to 3x4 matrix (translation part will be factorized):
        CALL convert_symop ( temp, sym )

!       Store 3x4 matrix in symop_result:
        symop_result%a = sym(1:3,1:3)
        symop_result%v = sym(1:3,4)

!       Main entry point for tsym_factor. Defined here only once:
        symop_result%tsym_factor = tsym_factor
    END SUBROUTINE char_string_to_symop_int

    SUBROUTINE symop_int_to_char_string (line_1, symop_int_2)
!
!
!       Purpose:
!       =======
!       Converts integer symop type (symop_int) to character string, e.g. 'X,Z,Y+1/2'
!
!       Date:            Programmer:                Description of changes:
!       ====             ==========                 ======================
!       Sep 2003         Strokopytov B.             Original code
!       Oct 2005         Strokopytov B.             SAVE statement added
!       Oct 2005         Strokopytov B.             Cosmetic changes 
!
        CHARACTER(len=60), INTENT(OUT)       :: line_1
        TYPE(symop_int),   INTENT(IN)        :: symop_int_2
!       Local variables:
        INTEGER,          DIMENSION(3,3)     :: a
        CHARACTER(len=1), DIMENSION(3), SAVE :: xyz = (/'x','y','z'/)
        CHARACTER(len=6)                     :: chartrn
        INTEGER                              :: first
        INTEGER, DIMENSION(3)                :: fract
        INTEGER                              :: denom
!       Counters:
        INTEGER                              :: i
        INTEGER                              :: j

        IF ( symop_int_2%tsym_factor <= 0 ) THEN
            WRITE(*,*) symop_int_2%tsym_factor
            CALL die('Programming error. tsym_factor has not been assigned.', 'symop_int_to_char_string')
        ENDIF
        line_1 = ''
        fract  = symop_int_2%v
        a      = symop_int_2%a
        rows: DO i = 1, 3
            
            translation: IF ( fract(i) /= 0 ) THEN
                denom = symop_int_2%tsym_factor
                DO j = 8, 2, -1
                    IF ( MOD( ABS( fract(i)), j ) == 0 .AND. MOD ( denom, j ) == 0 ) THEN
                        fract(i) = fract(i) / j
                        denom    = denom    / j
                    ENDIF
                ENDDO
                WRITE ( chartrn, '(I2,A1,I2)' ) fract(i), '/', denom

!               Eliminate blanks:
                CALL elimbl ( chartrn )
                line_1 = TRIM(line_1) // TRIM(chartrn)
            ENDIF translation

            first = 0

            ! No need to check if symop_int_2%v%v(i) /=0 
            ! since translation part goes first and + or - sign has to be used
            IF ( fract(i) == 0 ) THEN
                prerun: DO j = 1, 3
                    IF( a(i,j) == 1 ) THEN
                        first = j
                        line_1 = TRIM(line_1) // xyz(j)
                        EXIT
                    ENDIF
                ENDDO prerun
            ENDIF

            final: DO j = 1, 3
                noneed: IF ( j /= first ) THEN

                            IF ( a(i,j) == 1 ) THEN
                                line_1 = TRIM ( line_1 ) // '+' // xyz(j)
                            ELSE IF ( a(i,j) == -1 ) THEN
                                line_1 = TRIM ( line_1 ) // '-' // xyz(j)
                            ENDIF

                       check: IF (LEN_TRIM(line_1) > 59 ) THEN
                           CALL die ('Possible character length exceeded', 'symop_int_to_char_string')
                       ENDIF check
 
               ENDIF noneed
           ENDDO final

          IF ( i < 3 ) line_1 = TRIM(line_1)//','
        ENDDO rows
    END SUBROUTINE symop_int_to_char_string

    SUBROUTINE matrix_int_to_char_string( line_1, matrix_int_2 )
        CHARACTER(len=60), INTENT(OUT)       :: line_1
        TYPE(matrix_int),  INTENT(IN)        :: matrix_int_2
!       Local variables:
        INTEGER, DIMENSION(3,3)              :: a
        CHARACTER(len=1), DIMENSION(3), SAVE :: hkl = (/'h', 'k', 'l'/)
        INTEGER                              :: first
!       Counters:
        INTEGER                              :: i
        INTEGER                              :: j

        line_1 = ''
        a      = matrix_int_2
        rows: DO i = 1, 3
            first = 0
            prerun: DO j = 1, 3
                IF( a(i,j) == 1 ) THEN
                    first = j
                    line_1 = TRIM ( line_1 ) // hkl(j)
                    EXIT
                ENDIF
            ENDDO prerun

            final: DO j = 1, 3
                noneed: IF ( j /= first ) THEN

                   IF ( a(i,j) == 1 ) THEN
                        line_1 = TRIM ( line_1 ) // '+' // hkl(j)
                   ELSE IF ( a(i,j) == -1 ) THEN
                        line_1 = TRIM ( line_1 ) // '-' // hkl(j)
                   ENDIF
 
                   check: IF (LEN_TRIM ( line_1 ) > 59 ) THEN
                      CALL die ('Possible character length exceeded','matrix_int_to_char_string')
                   ENDIF check

                ENDIF noneed
            ENDDO final

            IF ( i < 3 ) line_1 = TRIM ( line_1 ) // ','
        ENDDO rows
    END SUBROUTINE matrix_int_to_char_string

    SUBROUTINE symop_int_to_array_int ( array_result, symop_1 )
        TYPE (symop_int),                 INTENT(IN)  :: symop_1
        INTEGER,          DIMENSION(:,:), INTENT(OUT) :: array_result

!       This requires some checks on upper bounds:
        array_result(1:3, 1:3) = symop_1%a
        array_result(1:3, 4  ) = symop_1%v
     END SUBROUTINE symop_int_to_array_int

     SUBROUTINE array_int_to_symop_int ( symop_result, sym )
!
!        Purpose:
!        =======
!        Converts integer matrix 3x4 into symop_int object
!
!
!        Date:                 Programmer:                History of changes:
!        ====                  ==========                 ==================
!        Sep 2003              Strokopytov B.             Original code
!        Oct 2005              Strokopytov B.             Array slices introduced (last 2 lines)

!
         TYPE(symop_int),                 INTENT(OUT) :: symop_result
         INTEGER,         DIMENSION(:,:), INTENT(IN)  :: sym
!        Local array:
         INTEGER,         DIMENSION(2)                :: extent

!        Checks for precaution:
         extent = UBOUND ( sym ) - LBOUND ( sym ) + 1
         IF ( extent(1) /= 3 .OR. extent(2) /= 4 ) THEN
             WRITE(*,*) ' Input array dims : ', extent
             CALL die ( 'Incompatible array size.', 'array_int_to_symop_int' )
         ENDIF
         symop_result%a = sym(1:3,1:3)
         symop_result%v = sym(1:3,4  )
     END SUBROUTINE array_int_to_symop_int

     SUBROUTINE symop_int_to_symop_double ( symop_result, sym )
         TYPE (symop)    , INTENT(OUT) :: symop_result
         TYPE (symop_int), INTENT(IN)  :: sym

         symop_result%a           = sym%a
         symop_result%v           = sym%v / REAL ( sym%tsym_factor, KIND=wp )
         symop_result%tsym_factor = sym%tsym_factor
     END SUBROUTINE symop_int_to_symop_double

     SUBROUTINE symop_double_to_symop_int ( symop_result_int, sym )
         TYPE(symop_int), INTENT(OUT) :: symop_result_int
         TYPE(symop),     INTENT(IN)  :: sym

         symop_result_int%a           = NINT ( sym%a )
         symop_result_int%v           = MOD ( NINT ( sym%v * REAL ( sym%tsym_factor, KIND=wp ) ), sym%tsym_factor )
         symop_result_int%tsym_factor = sym%tsym_factor
     END SUBROUTINE symop_double_to_symop_int

     SUBROUTINE symop_int_to_array_real ( array_result, sym )
!
!        Purpose:
!        =======
!        Conversion to CCP4-style operator 4x4 in single precision
!        Barely needed and used
!
!        Date:                 Programmer:                History of changes:
!        ====                  ==========                 ==================
!        Sep 2003              Strokopytov B.             Original code
!
!
         REAL(KIND=wp), DIMENSION(4,4), INTENT(OUT) :: array_result
         TYPE(symop_int),               INTENT(IN)  :: sym

         array_result          = 0.0_wp
         array_result(4,4)     = 1.0_wp
         array_result(1:3,1:3) = sym%a
         array_result(1:3,4)   = sym%v / REAL ( sym%tsym_factor, KIND=wp ) 
     END SUBROUTINE symop_int_to_array_real

     SUBROUTINE wp_array_to_symop ( sym_1, arr_2 )
         TYPE(symop),                   INTENT(OUT) :: sym_1
         REAL(KIND=wp), DIMENSION(3,4), INTENT(IN)  :: arr_2
 
         sym_1%a = arr_2(1:3,1:3)
         sym_1%v = arr_2(1:3,4)
     END SUBROUTINE wp_array_to_symop 

     SUBROUTINE set_symop_int ( symop_result, constant )
!
!        Purpose:
!        =======
!        Normally used to setup unity operator
!
         TYPE(symop_int), INTENT(OUT) :: symop_result
         INTEGER,          INTENT(IN) :: constant

!        Matrix part:
         symop_result%a = 0
         symop_result%a = .DIAG. constant

!        Vector part:
         symop_result%v = 0

!        Entry point for tsym_factor:
         symop_result%tsym_factor = tsym_factor
     END SUBROUTINE set_symop_int

     SUBROUTINE set_symop ( symop_result, constant )
!
!        Purpose:
!        =======
!        Normally used to setup unity operator
!
         TYPE (symop),   INTENT(OUT) :: symop_result
         REAL (KIND=wp), INTENT(IN)  :: constant
! 
         symop_result%a = .DIAG. constant
         symop_result%v = 0.0_wp
         symop_result%tsym_factor = tsym_factor
     END SUBROUTINE set_symop

     FUNCTION set_symop_vector_part_vec ( symop_1, vec_2 ) RESULT ( symop_result )
!
!        Purpose:
!        =======
!        Sets vector part of symmetry operator
!
         TYPE(symop)              :: symop_result
         TYPE(symop),  INTENT(IN) :: symop_1
         TYPE(vector), INTENT(IN) :: vec_2

         symop_result%a = symop_1%a
!        Set translation part:
         symop_result%v = vec_2
         
     END FUNCTION set_symop_vector_part_vec

     FUNCTION set_symop_vector_part_arr ( symop_1, arr_2 ) RESULT ( symop_result )
!
!        Purpose:
!        =======
!        Sets vector part of symmetry operator
!
         TYPE(symop)                             :: symop_result
         TYPE(symop),                 INTENT(IN) :: symop_1
         REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: arr_2

         symop_result%a = symop_1%a

!        Set translation part:
         symop_result%v = arr_2

     END FUNCTION set_symop_vector_part_arr

     FUNCTION set_symop_int_vector_part_vec ( symop_1, vec_2 ) RESULT ( symop_result )
!
!        Purpose:
!        =======
!        Sets vector part of symmetry operator (integer version) 
!
         TYPE(symop_int)              :: symop_result
         TYPE(symop_int),  INTENT(IN) :: symop_1
         TYPE(vector_int), INTENT(IN) :: vec_2

         symop_result%a = symop_1%a

!        Set translation part:
         symop_result%v = vec_2
     END FUNCTION set_symop_int_vector_part_vec

     FUNCTION set_symop_int_vector_part_arr ( symop_1, arr_2 ) RESULT ( symop_result )
!
!        Purpose:
!        =======
!        Sets vector part of symmetry operator (integer version)
!
         TYPE(symop_int)                   :: symop_result
         TYPE(symop_int),       INTENT(IN) :: symop_1
         INTEGER, DIMENSION(:), INTENT(IN) :: arr_2

         symop_result%a = symop_1%a

!        Set translation part:
         symop_result%v = arr_2
     END FUNCTION set_symop_int_vector_part_arr

     FUNCTION extract_symop_vector_part ( symop_1 ) RESULT ( vec_1 )
         TYPE(vector)            :: vec_1
         TYPE(symop), INTENT(IN) :: symop_1
     
         vec_1 = symop_1%v
     END FUNCTION extract_symop_vector_part

     FUNCTION extract_symop_int_vector_part ( symop_1 ) RESULT ( vec_1 )
         TYPE(vector_int)            :: vec_1
         TYPE(symop_int), INTENT(IN) :: symop_1

         vec_1 = symop_1%v
     END FUNCTION extract_symop_int_vector_part

     FUNCTION set_symop_int_matrix_part ( symop_1, int_2 ) RESULT ( symop_result )
         TYPE(symop_int)             :: symop_result
         TYPE(symop_int), INTENT(IN) :: symop_1
         INTEGER,         INTENT(IN) :: int_2

         symop_result%a = .DIAG. int_2
         symop_result%v = symop_1%v
     END FUNCTION set_symop_int_matrix_part

     FUNCTION set_symop_matrix_part ( symop_1, const ) RESULT ( symop_result )
         TYPE(symop)               :: symop_result
         TYPE(symop),   INTENT(IN) :: symop_1
         REAL(KIND=wp), INTENT(IN) :: const 

         symop_result%a = .DIAG. const
         symop_result%v = symop_1%v
     END FUNCTION set_symop_matrix_part
 
     FUNCTION set_symop_matrix_part_mat ( symop_1, mat_2 ) RESULT ( symop_result )
         TYPE(symop)             :: symop_result
         TYPE(symop), INTENT(IN) :: symop_1
         TYPE(matrix), INTENT(IN) :: mat_2

         symop_result%a = mat_2
         symop_result%v = symop_1%v
     END FUNCTION set_symop_matrix_part_mat

     FUNCTION extract_symop_int_matrix_part ( symop_1 ) RESULT ( mat_1 )
         TYPE(matrix_int)            :: mat_1
         TYPE(symop_int), INTENT(IN) :: symop_1

         mat_1 = symop_1%a 
     END FUNCTION extract_symop_int_matrix_part

     FUNCTION extract_symop_matrix_part ( symop_1 ) RESULT ( mat_1 )
         TYPE(matrix)            :: mat_1
         TYPE(symop), INTENT(IN) :: symop_1

         mat_1 = symop_1%a
     END FUNCTION extract_symop_matrix_part
       
     FUNCTION mod_tsym_int ( vec_1_int, factor )
!
!        Purpose:
!        =======
!        Make integer vector positive MODULO(factor)
!        Needed for operators comparison
!        It might be that we have to compare say, -1/4 and +3/4 (which is the same thing!)
!
!        Note: 
!        =====
!        This routine and the next one probably should be moved to `vectors' module
!
         TYPE (vector_int)             :: mod_tsym_int
         TYPE (vector_int), INTENT(IN) :: vec_1_int
         INTEGER,           INTENT(IN) :: factor

         mod_tsym_int = MOD ( vec_1_int + factor * 10000, factor )

     END FUNCTION mod_tsym_int

     FUNCTION mod_tsym ( vec_1, factor )
!
!        Purpose:
!        =======
!        Make real vector positive MODULO(factor)
!        Basically this is for operator comparison (translation part)
!        Clearly it's more convenient manipulating integer symops than real ones.
!
         TYPE (vector)             :: mod_tsym
         TYPE (vector), INTENT(IN) :: vec_1
         REAL(KIND=wp), INTENT(IN) :: factor

!        Array to array conversion is faster:
         mod_tsym = MOD ( vec_1 + factor * 100.0_wp, factor )
     END FUNCTION mod_tsym

     FUNCTION symop_int_eq_symop_int ( sym_1, sym_2 )
         LOGICAL                     :: symop_int_eq_symop_int
         TYPE(symop_int), INTENT(IN) :: sym_1
         TYPE(symop_int), INTENT(IN) :: sym_2

         symop_int_eq_symop_int = ( sym_1%a == sym_2%a ) .AND. & 
         mod_tsym_int ( sym_1%v, sym_1%tsym_factor ) == mod_tsym_int ( sym_2%v, sym_1%tsym_factor )
     END FUNCTION symop_int_eq_symop_int

     FUNCTION symop_int_ne_symop_int ( sym_1, sym_2 )
         LOGICAL                     :: symop_int_ne_symop_int
         TYPE(symop_int), INTENT(IN) :: sym_1
         TYPE(symop_int), INTENT(IN) :: sym_2

         symop_int_ne_symop_int = .NOT. symop_int_eq_symop_int ( sym_1, sym_2 ) 
     END FUNCTION symop_int_ne_symop_int
   
     FUNCTION symop_eq_symop ( sym_1, sym_2 )
         LOGICAL                 :: symop_eq_symop
         TYPE(symop), INTENT(IN) :: sym_1
         TYPE(symop), INTENT(IN) :: sym_2
!        Local variables:
         TYPE(symop_int)         :: sym_int_1
         TYPE(symop_int)         :: sym_int_2

!        Safest way to compare real operators:
 
         sym_int_1 = sym_1
         sym_int_2 = sym_2

!        After conversion to integer operators it should be a snap
         symop_eq_symop = ( sym_int_1 == sym_int_2 )
     END FUNCTION symop_eq_symop

     FUNCTION symop_ne_symop ( sym_1, sym_2 )
         LOGICAL                 :: symop_ne_symop
         TYPE(symop), INTENT(IN) :: sym_1
         TYPE(symop), INTENT(IN) :: sym_2

         symop_ne_symop = .NOT. symop_eq_symop ( sym_1, sym_2 )
     END FUNCTION symop_ne_symop

     FUNCTION symop_int_eq_unity ( sym_1, int_2 )
         LOGICAL                     :: symop_int_eq_unity
         TYPE(symop_int), INTENT(IN) :: sym_1
         INTEGER,         INTENT(IN) :: int_2

         symop_int_eq_unity = (sym_1%a == int_2 ) .AND. ( sym_1%v == 0 )
     END FUNCTION symop_int_eq_unity

     FUNCTION symop_int_ne_unity ( sym_1, int_2 )
!
!        Purpose:
!        =======
!        Normally used to check whether given symmetry operator is 'X,Y,Z'
!
         LOGICAL                     :: symop_int_ne_unity
         TYPE(symop_int), INTENT(IN) :: sym_1
         INTEGER, INTENT(IN)         :: int_2

!        See vectors module for details:
         symop_int_ne_unity = .NOT. ( (sym_1%a == int_2 ) .AND. ( sym_1%v == 0 ) )
     END FUNCTION symop_int_ne_unity

     FUNCTION symop_eq_unity ( sym_1, double_2 )
!
!        Purpose:
!        =======
!        Normally used to check whether given symmetry operator is 'X,Y,Z'
!
         LOGICAL                   :: symop_eq_unity
         TYPE(symop),   INTENT(IN) :: sym_1
         REAL(KIND=wp), INTENT(IN) :: double_2

!        See vectors module for details:
         symop_eq_unity = ( sym_1%a == double_2 ) .AND. ( sym_1%v == 0.0_wp )
     END FUNCTION symop_eq_unity

     FUNCTION symop_ne_unity ( sym_1, double_2 )
         LOGICAL                   :: symop_ne_unity
         TYPE(symop),   INTENT(IN) :: sym_1
         REAL(KIND=wp), INTENT(IN) :: double_2

         symop_ne_unity = .NOT. ( ( sym_1%a == double_2 ) .AND. ( sym_1%v == 0.0_wp ) )
     END FUNCTION symop_ne_unity

     FUNCTION symop_int_add_symop_int ( sym_1, sym_2 )
         TYPE(symop_int)             :: symop_int_add_symop_int
         TYPE(symop_int), INTENT(IN) :: sym_1
         TYPE(symop_int), INTENT(IN) :: sym_2

         symop_int_add_symop_int%a           = sym_1%a + sym_2%a
         symop_int_add_symop_int%v           = MOD ( sym_1%v + sym_2%v + 100*sym_2%tsym_factor, sym_2%tsym_factor )
         symop_int_add_symop_int%tsym_factor = sym_2%tsym_factor
     END FUNCTION symop_int_add_symop_int

     FUNCTION symop_add_symop ( sym_1, sym_2 )
         TYPE (symop)             :: symop_add_symop
         TYPE (symop), INTENT(IN) :: sym_1
         TYPE (symop), INTENT(IN) :: sym_2

         symop_add_symop%a           = sym_1%a + sym_2%a
         symop_add_symop%v           = sym_1%v + sym_2%v
         symop_add_symop%tsym_factor = sym_2%tsym_factor 
     END FUNCTION symop_add_symop

     FUNCTION symop_add_vector_int ( sym_1, vec_2_int )
         TYPE ( symop )                  :: symop_add_vector_int
         TYPE ( symop ), INTENT(IN)      :: sym_1
         TYPE ( vector_int ), INTENT(IN) :: vec_2_int

         symop_add_vector_int%a           = sym_1%a
         symop_add_vector_int%v           = sym_1%v + vec_2_int
         symop_add_vector_int%tsym_factor = sym_1%tsym_factor
     END FUNCTION symop_add_vector_int

     FUNCTION symop_add_vector ( sym_1, vec_2 )
         TYPE ( symop )              :: symop_add_vector
         TYPE ( symop ),  INTENT(IN) :: sym_1
         TYPE ( vector ), INTENT(IN) :: vec_2

         symop_add_vector%a           = sym_1%a
         symop_add_vector%v           = sym_1%v + vec_2
         symop_add_vector%tsym_factor = sym_1%tsym_factor
     END FUNCTION symop_add_vector

     FUNCTION symop_int_subtract_symop_int ( sym_1, sym_2 )
         TYPE(symop_int)             :: symop_int_subtract_symop_int
         TYPE(symop_int), INTENT(IN) :: sym_1
         TYPE(symop_int), INTENT(IN) :: sym_2

         symop_int_subtract_symop_int%a = sym_1%a - sym_2%a
         symop_int_subtract_symop_int%v = MOD ( sym_1%v - sym_2%v + 100 * sym_1%tsym_factor, sym_1%tsym_factor )
         symop_int_subtract_symop_int%tsym_factor = sym_2%tsym_factor
     END FUNCTION symop_int_subtract_symop_int
 
     FUNCTION symop_subtract_symop ( sym_1, sym_2)
         TYPE(symop)             :: symop_subtract_symop
         TYPE(symop), INTENT(IN) :: sym_1
         TYPE(symop), INTENT(IN) :: sym_2

         symop_subtract_symop%a = sym_1%a - sym_2%a
         symop_subtract_symop%v = sym_1%v - sym_2%v
         symop_subtract_symop%tsym_factor = sym_2%tsym_factor
     END FUNCTION symop_subtract_symop

     FUNCTION symop_subtract_vector_int ( sym_1, vec_2_int )
         TYPE ( symop )                  :: symop_subtract_vector_int
         TYPE ( symop ), INTENT(IN)      :: sym_1
         TYPE ( vector_int ), INTENT(IN) :: vec_2_int

         symop_subtract_vector_int%a = sym_1%a
         symop_subtract_vector_int%v = sym_1%v - vec_2_int
         symop_subtract_vector_int%tsym_factor = sym_1%tsym_factor
     END FUNCTION symop_subtract_vector_int

     FUNCTION symop_subtract_vector ( sym_1, vec_2 )
         TYPE ( symop )              :: symop_subtract_vector
         TYPE ( symop ),  INTENT(IN) :: sym_1
         TYPE ( vector ), INTENT(IN) :: vec_2
         symop_subtract_vector%a = sym_1%a
         symop_subtract_vector%v = sym_1%v - vec_2
         symop_subtract_vector%tsym_factor = sym_1%tsym_factor
     END FUNCTION symop_subtract_vector

     FUNCTION symop_int_times_symop_int ( sym_1, sym_2 )
         TYPE(symop_int)             :: symop_int_times_symop_int
         TYPE(symop_int), INTENT(IN) :: sym_1
         TYPE(symop_int), INTENT(IN) :: sym_2
         symop_int_times_symop_int%a = sym_1%a * sym_2%a
         symop_int_times_symop_int%v = mod_tsym_int ( sym_1%a * sym_2%v + sym_1%v, sym_2%tsym_factor )
         symop_int_times_symop_int%tsym_factor = sym_2%tsym_factor
     END FUNCTION symop_int_times_symop_int

     FUNCTION symop_times_symop ( sym_1, sym_2 )
         TYPE (symop)             :: symop_times_symop
         TYPE (symop), INTENT(IN) :: sym_1
         TYPE (symop), INTENT(IN) :: sym_2

         symop_times_symop%a = sym_1%a * sym_2%a
         symop_times_symop%v = ( sym_1%a * sym_2%v ) + sym_1%v
         symop_times_symop%tsym_factor = sym_2%tsym_factor
     END FUNCTION symop_times_symop

     FUNCTION symop_int_times_vector_int ( sym_1, vec_2 )
         TYPE(vector)                 :: symop_int_times_vector_int
         TYPE(symop_int),  INTENT(IN) :: sym_1
         TYPE(vector_int), INTENT(IN) :: vec_2

         symop_int_times_vector_int = sym_1%a * vec_2 + ( sym_1%v / REAL ( sym_1%tsym_factor, KIND=wp ) )
     END FUNCTION symop_int_times_vector_int
     
     FUNCTION symop_int_times_vector ( sym_1, vec_2 )
         TYPE (vector)                :: symop_int_times_vector
         TYPE (symop_int), INTENT(IN) :: sym_1
         TYPE (vector),    INTENT(IN) :: vec_2

         symop_int_times_vector = sym_1%a * vec_2 + ( sym_1%v / REAL ( sym_1%tsym_factor, KIND=wp ) )
     END FUNCTION symop_int_times_vector

     FUNCTION symop_times_vector ( sym_1, vec_2 )
         TYPE(vector)             :: symop_times_vector
         TYPE(symop),  INTENT(IN) :: sym_1
         TYPE(vector), INTENT(IN) :: vec_2

         symop_times_vector = sym_1%a * vec_2 + sym_1%v
     END FUNCTION symop_times_vector

     FUNCTION matrix_int_times_symop_int ( mat_1, sym_2 )
         TYPE(symop_int)              :: matrix_int_times_symop_int
         TYPE(matrix_int), INTENT(IN) :: mat_1
         TYPE(symop_int),  INTENT(IN) :: sym_2

         matrix_int_times_symop_int%a = mat_1 * sym_2%a
         matrix_int_times_symop_int%v = mat_1 * sym_2%v
         matrix_int_times_symop_int%tsym_factor = sym_2%tsym_factor
     END FUNCTION matrix_int_times_symop_int

     FUNCTION matrix_int_times_symop ( mat_1, sym_2 )
         TYPE(symop)                   :: matrix_int_times_symop
         TYPE(matrix_int ), INTENT(IN) :: mat_1
         TYPE(symop),       INTENT(IN) :: sym_2

         matrix_int_times_symop%a = mat_1 * sym_2%a
         matrix_int_times_symop%v = mat_1 * sym_2%v
         matrix_int_times_symop%tsym_factor = sym_2%tsym_factor
     END FUNCTION matrix_int_times_symop

     FUNCTION matrix_times_symop ( mat_1, sym_2 )
         TYPE(symop)              :: matrix_times_symop
         TYPE(matrix), INTENT(IN) :: mat_1
         TYPE(symop),  INTENT(IN) :: sym_2

         matrix_times_symop%a = mat_1 * sym_2%a
         matrix_times_symop%v = mat_1 * sym_2%v
         matrix_times_symop%tsym_factor = sym_2%tsym_factor
     END FUNCTION matrix_times_symop

     FUNCTION symop_times_matrix ( sym_1, mat_2 )
         TYPE(symop)              :: symop_times_matrix
         TYPE(symop),  INTENT(IN) :: sym_1
         TYPE(matrix), INTENT(IN) :: mat_2

         symop_times_matrix%a = sym_1%a * mat_2
         symop_times_matrix%v = sym_1%v
         symop_times_matrix%tsym_factor = sym_1%tsym_factor
     END FUNCTION symop_times_matrix

     FUNCTION symop_int_times_array_int ( sym_1, arr_2 )
         TYPE(vector)                      :: symop_int_times_array_int
         TYPE(symop_int), INTENT(IN)       :: sym_1
         INTEGER, DIMENSION(:), INTENT(IN) :: arr_2
!        Local variables:
         TYPE(vector_int)                  :: vec_2

         vec_2                     = arr_2
         symop_int_times_array_int = sym_1 * vec_2     
     END FUNCTION symop_int_times_array_int

     FUNCTION symop_int_times_array_double ( sym_1, arr_2 )
         TYPE(vector)                            :: symop_int_times_array_double
         TYPE(symop_int),             INTENT(IN) :: sym_1
         REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: arr_2
!        Local variables:
         TYPE(vector)                            :: vec_2
         vec_2                        = arr_2
         symop_int_times_array_double = sym_1 * vec_2
     END FUNCTION symop_int_times_array_double

     FUNCTION string_times_string ( string_1, string_2 )
         CHARACTER(len=60)            :: string_times_string
         CHARACTER(len=*), INTENT(IN) :: string_1, string_2
!        Local variables:
         TYPE(symop_int)              :: sym_1
         TYPE(symop_int)              :: sym_2

         sym_1               = string_1
         sym_2               = string_2
         string_times_string = sym_1 * sym_2
     END FUNCTION string_times_string

     SUBROUTINE matrix_part_int ( mat_1_int, sym_1_int )
         TYPE(matrix_int), INTENT(OUT) :: mat_1_int
         TYPE(symop_int),  INTENT(IN)  :: sym_1_int

         mat_1_int = sym_1_int%a
     END SUBROUTINE matrix_part_int

     SUBROUTINE matrix_part ( mat_1, sym_1 )
         TYPE(matrix), INTENT(OUT) :: mat_1
         TYPE(symop),  INTENT(IN)  :: sym_1

         mat_1 = sym_1%a
     END SUBROUTINE matrix_part

     SUBROUTINE vector_part_int ( vec_1, sym_1_int )
         TYPE(vector),    INTENT(OUT) :: vec_1
         TYPE(symop_int), INTENT(IN)  :: sym_1_int 

         vec_1 = sym_1_int%v / REAL ( tsym_factor, KIND=wp )
     END SUBROUTINE vector_part_int

     SUBROUTINE vector_part ( vec_1, sym_1 )
         TYPE(vector), INTENT(OUT) :: vec_1
         TYPE(symop),  INTENT(IN)  :: sym_1
         vec_1 = sym_1%v
     END SUBROUTINE vector_part

     FUNCTION vector_int_dot_symop ( vec_1, symop_2 )
         REAL(KIND=wp)                :: vector_int_dot_symop
         TYPE(vector_int), INTENT(IN) :: vec_1
         TYPE(symop),      INTENT(IN) :: symop_2

         vector_int_dot_symop = vec_1 .DOT. symop_2%v
     END FUNCTION vector_int_dot_symop

     FUNCTION symop_dot_vector_int ( symop_1, vec_2 )
         REAL(KIND=wp) :: symop_dot_vector_int
         TYPE(symop),      INTENT(IN) :: symop_1
         TYPE(vector_int), INTENT(IN) :: vec_2

         symop_dot_vector_int = symop_1%v .DOT. vec_2
     END FUNCTION symop_dot_vector_int

!    ADDED OCT 2007:
     FUNCTION invert_symop ( symop_1 )
!
!        Note:
!        ====
!        Because of translation modification we can't use SYMOP TYPE as general rotational operator:
         TYPE(symop)             :: invert_symop
         TYPE(symop), INTENT(IN) :: symop_1

         invert_symop%a = .INV. symop_1%a
         invert_symop%tsym_factor = 24

!        Some mudification to get rid of whole translations, like 1.0, 2.0 etc.:
         invert_symop%v = MOD ( NINT ( -symop_1%v * invert_symop%tsym_factor ) + 24 * 100, 24 ) / 24
     END FUNCTION invert_symop

!    ADDED OCT 2007:
     FUNCTION invert_symop_int ( symop_1 )
         TYPE(symop_int)             :: invert_symop_int
         TYPE(symop_int), INTENT(IN) :: symop_1

         invert_symop_int%a = .INV. symop_1%a
         invert_symop_int%v = -symop_1%v
         invert_symop_int%tsym_factor = 24
     END FUNCTION invert_symop_int

     SUBROUTINE check_duplicate_symops ( sym_int )
         TYPE(symop_int), DIMENSION(:), INTENT(INOUT) :: sym_int
!        Counters:
         INTEGER                                      :: i
         INTEGER                                      :: j
!
         DO i = 1, SIZE ( sym_int ) - 1
             DO j = i + 1, SIZE ( sym_int )
                 IF ( sym_int(i) == sym_int(j) ) THEN
                     WRITE(*, '( '' FIND_SYMOP_IN_LIBRARY> '',&
                              &  '' Duplicate symmetry operators : '', 2I5)' ) i, j
                     CALL die ( 'Duplicate symops found...', 'check_duplicate_symops' )
                 ENDIF
            ENDDO
        ENDDO
        CALL messag ( 'No duplicates found...', 'check_duplicate_symops' )
     END SUBROUTINE check_duplicate_symops

     SUBROUTINE eliminate_duplicate_symops ( sym_int )
!
!        Purpose:
!        =======
!        Removes duplicate symops. Basically needed for cheshirian operators.
!
!        Date:                    Programmer:                History of changes:
!        ====                     ==========                 ==================
!        Jul 2006 - Mar 2007      B.Strokopytov              Original code
!        Oct 2007                 B.Strokopytov              Cosmetic changes.                   
!
!
         TYPE(symop_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: sym_int
!        Local variables:
         INTEGER                                                   :: istat
         TYPE(symop_int), DIMENSION(:), ALLOCATABLE                :: sym_temp
         LOGICAL(KIND=1), DIMENSION(:), ALLOCATABLE                :: deleted
!        Counters:
         INTEGER                                                   :: i
         INTEGER                                                   :: j
!        S/R name:
         CHARACTER(LEN=32)                                         :: srname

         srname = 'eliminate_duplicate_symops'

         CALL allocate_array ( deleted, SIZE ( sym_int ) )
         deleted = .FALSE.

!        Check for duplicates:
         DO i = 1, SIZE ( sym_int ) - 1
             IF ( deleted(i) ) CYCLE
             DO j = i + 1, SIZE ( sym_int )
                 IF ( sym_int(i) == sym_int(j) ) THEN
                     IF ( debug > 5 ) THEN
                         WRITE(*, "( ' ELIMINATE_DUPLICATE_SYMOPS> ',&
                                  &  ' Duplicate symmetry operators : ', 2I5)" ) i, j
                     ENDIF
                     deleted(j) = .TRUE.
                 ENDIF
            ENDDO
        ENDDO

        WRITE(*,"(' ELIMINATE_DUPLICATE_SYMOPS> ', A, ' duplicates found...')") &
        TRIM ( int_to_c ( COUNT ( deleted ) ) )

        IF ( COUNT ( deleted ) == 0 ) RETURN

!       Rearrange sym_int using temp array:
        ALLOCATE ( sym_temp ( COUNT ( .NOT. deleted ) ), STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat
            CALL die ('Failed to allocate SYM_TEMP.', srname )
        ENDIF

!       Remove operators marked as deleted:
        j = 0
        DO i = 1, SIZE ( sym_int )
            IF ( .NOT. deleted(i) ) THEN
                j = j + 1
                sym_temp(j) = sym_int(i)
            ENDIF
        ENDDO

!       Reallocate SYM_INT with correct size:
        DEALLOCATE ( sym_int, STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat
            CALL die ('Failed to deallocate SYM_INT.', srname )
        ENDIF

        ALLOCATE  ( sym_int(SIZE ( sym_temp ) ), STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat
            CALL die ('Failed to allocate SYM_INT.', srname )
        ENDIF

!       Save final result:
        sym_int = sym_temp

!       Free memory:
        DEALLOCATE ( sym_temp, STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat
            CALL die ('Failed to deallocate SYM_TEMP.', srname )
        ENDIF

        CALL deallocate_array ( deleted )

     END SUBROUTINE eliminate_duplicate_symops
       
     SUBROUTINE check_symop_times_symop_matrix ( sym_int )
         TYPE(symop_int),  DIMENSION(:),    INTENT(IN) :: sym_int
!        Local variables:
         INTEGER(KIND=fb), DIMENSION(:,:), ALLOCATABLE :: sym_times_sym
         CHARACTER(len=60)                             :: sym_i
         CHARACTER(len=60)                             :: sym_j
         CHARACTER(len=60)                             :: sym_k
         LOGICAL                                       :: closed_group      
!        Counters:
         INTEGER                                       :: i
         INTEGER                                       :: j
         INTEGER                                       :: k

         CALL messag ( 'Testing whether operators form a closed group...',&
                      'check_symop_times_symop_matrix' )
         CALL allocate_array ( sym_times_sym, SIZE ( sym_int ), SIZE ( sym_int ) )

         sym_times_sym = 0
         DO i = 1, SIZE ( sym_int )
             DO j = 1, SIZE ( sym_int )
                 closed_group = .FALSE.
                 DO k = 1, SIZE ( sym_int )
                     IF ( sym_int(i) * sym_int(j) == sym_int(k) ) THEN
                         closed_group        = .TRUE.
                         sym_times_sym(i, j) = k
                         EXIT
                     ENDIF
                 ENDDO
                 IF ( .NOT. closed_group ) THEN
                     sym_i = sym_int(i)
                     sym_j = sym_int(j)
                     sym_k = sym_int(i) * sym_int(j)
                     WRITE(*,*) ' OP_1       : ', TRIM ( sym_i ), ' OP_2: ', TRIM ( sym_j )
                     WRITE(*,*) ' OP_1 x OP_2: ', TRIM ( sym_k )
                     CALL messag &
                     (           &
                      'Operators #'//TRIM(int_to_c(i))//' and #'//TRIM(int_to_c(j))//&
                      ' produce operator=('//TRIM(sym_i*sym_j)//') outside this group.', 'find_symop_in_library'&
                     )
                     CALL die ( 'This space group is not closed...', 'check_symop_times_symop_matrix' )
                 ENDIF
             ENDDO
         ENDDO

!        print out matrix to confirm formation of closed group:
         IF ( debug > 15 ) THEN
             DO i = 1, SIZE ( sym_int )
                 WRITE(*,"(' SYMOP*SYMOP> ', 192I4)") &
                 (sym_times_sym(i,j), j = 1, SIZE ( sym_int ))
             ENDDO
        ENDIF
        CALL deallocate_array ( sym_times_sym )
        CALL messag('Looks good...', 'check_symop_times_symop_matrix')
     END SUBROUTINE check_symop_times_symop_matrix

     SUBROUTINE print_symop_int ( symop_1, tsym_factor )
         TYPE(symop_int), INTENT(IN)           :: symop_1
         INTEGER,         INTENT(IN), OPTIONAL :: tsym_factor
!        Local vars:
         CHARACTER(LEN=80)                     :: sym_int_1

!        Convert to character string first:
         sym_int_1 = symop_1
         IF ( PRESENT ( tsym_factor ) ) THEN
             WRITE(*,"(' PRINT_SYMOP_INT> ', A, ' tsym_factor: ',A)") TRIM ( ADJUSTL ( sym_int_1 ) ), &
                                                                      TRIM ( int_to_c ( symop_1%tsym_factor ) )
         ELSE
            WRITE(*,"(' PRINT_SYMOP_INT> ', A)") TRIM ( ADJUSTL ( sym_int_1 ) )
         ENDIF
    END SUBROUTINE print_symop_int

    SUBROUTINE print_symop_double ( symop_1, tsym_factor )
         TYPE(symop), INTENT(IN)           :: symop_1
         INTEGER,     INTENT(IN), OPTIONAL :: tsym_factor
!        Local vars:
         TYPE(symop_int)                   :: symop_int_1
         CHARACTER(LEN=60)                 :: sym_1

!        Convert to symop_int first otherwise cannot convert to character string:
         symop_int_1 = symop_1

!         Convert to character string:
          sym_1 = symop_int_1

!         How to print:
          IF ( PRESENT ( tsym_factor ) ) THEN
              WRITE(*,"(' PRINT_SYMOP> ', A, ' tsym_factor= ',A)") TRIM ( ADJUSTL ( sym_1 ) ), &
                                                                  TRIM ( int_to_c ( symop_1%tsym_factor ) )
          ELSE
              WRITE(*,"(' PRINT_SYMOP> ', A)") TRIM ( ADJUSTL ( sym_1 ) )
          ENDIF
 
    END SUBROUTINE print_symop_double

    SUBROUTINE set_tsym_factor ( symop_int_arr, tsym_factor_2 )
         TYPE(symop_int), DIMENSION(:), INTENT(INOUT) :: symop_int_arr
         INTEGER,                       INTENT(IN)    :: tsym_factor_2

         symop_int_arr%tsym_factor = tsym_factor_2
    END SUBROUTINE set_tsym_factor


END MODULE basic_symmetry_operations
