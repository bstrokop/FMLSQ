MODULE aniso_manip
USE constants
USE fail
USE mkl95_lapack , ONLY : syevd
USE vectors
USE select_kinds
USE util
IMPLICIT NONE
TYPE :: aniso
    REAL(KIND=wp), DIMENSION(6) :: u
END TYPE
PUBLIC 
PRIVATE array_to_2d_array 
INTERFACE allocate_array
    MODULE PROCEDURE allocate_aniso_array
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_aniso_array
END INTERFACE

INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE set_to_const
    MODULE PROCEDURE aniso_to_array
    MODULE PROCEDURE array_to_aniso
    MODULE PROCEDURE aniso_to_aniso
    MODULE PROCEDURE aniso_to_mat
END INTERFACE

INTERFACE OPERATOR ( .CONVERT. )
    MODULE PROCEDURE b_to_u
    MODULE PROCEDURE u_to_b
    MODULE PROCEDURE array_to_2d_array
    MODULE PROCEDURE array_2d_to_1d_array
!    MODULE PROCEDURE array_to_mat  ! redundant
END INTERFACE

INTERFACE OPERATOR(+)
    MODULE PROCEDURE add_u_to_u
    MODULE PROCEDURE add_b_to_u
    MODULE PROCEDURE add_u_To_b
END INTERFACE

INTERFACE OPERATOR(-) 
    MODULE PROCEDURE u_minus_u
END INTERFACE

INTERFACE OPERATOR(*)
    MODULE PROCEDURE scale_u_by_const
    MODULE PROCEDURE scale_const_by_u
END INTERFACE

INTERFACE OPERATOR(/)
    MODULE PROCEDURE divide_u
END INTERFACE

INTERFACE OPERATOR ( > )
    MODULE PROCEDURE gt_const
END INTERFACE

INTERFACE OPERATOR ( >= ) 
    MODULE PROCEDURE ge_const
END INTERFACE

INTERFACE OPERATOR ( == )
    MODULE PROCEDURE eq_const
END INTERFACE

INTERFACE SUM
    MODULE PROCEDURE usum3
END INTERFACE

INTERFACE MAXVAL
    MODULE PROCEDURE max3val
END INTERFACE

INTERFACE MINVAL
    MODULE PROCEDURE min3val
END INTERFACE

INTERFACE OPERATOR ( .CORRANI. )
    MODULE PROCEDURE correct_aniso_tensor
END INTERFACE
CONTAINS
    SUBROUTINE allocate_aniso_array(U, n)
        TYPE(aniso), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: U
        INTEGER,                                INTENT(IN)    :: n
!       Local variables:
        INTEGER                                               :: istat

        IF ( ALLOCATED ( U ) ) THEN
            CALL deallocate_array ( U )
        ENDIF

        ALLOCATE(U(n), STAT=istat)

        IF ( istat /= 0 ) THEN
            CALL die('Failed to allocate U array.', &
                     'allocate_aniso_array')
        ENDIF 

    END SUBROUTINE allocate_aniso_array

    SUBROUTINE deallocate_aniso_array(U)
        TYPE(aniso), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: U
!       Local variables:
        INTEGER                                               :: istat

        IF ( ALLOCATED ( U ) ) THEN

            DEALLOCATE(U, STAT=istat)

            IF ( istat /= 0 ) THEN
                CALL die('Failed to deallocate U array.', &
                         'deallocate_aniso_array')
            ENDIF

        ENDIF

    END SUBROUTINE deallocate_aniso_array

    FUNCTION indmat() RESULT ( a )
        INTEGER, DIMENSION(2,9) :: a
        a(1,1) = 1; a(2,1) = 1
        a(1,2) = 2; a(2,2) = 2
        a(1,3) = 3; a(2,3) = 3

        a(1,4) = 1; a(2,4) = 2
        a(1,5) = 1; a(1,5) = 3
        a(1,6) = 2; a(2,6) = 3

        a(1,7) = a(2,4); a(2,7) = a(1,4)
        a(1,8) = a(2,5); a(2,8) = a(1,5)
        a(1,9) = a(2,6); a(2,9) = a(1,6)        
    END FUNCTION indmat

    SUBROUTINE set_to_const (u, const)
        TYPE(aniso),   INTENT(OUT) :: u
        REAL(KIND=wp), INTENT(IN)  :: const
        u%u = const
    END SUBROUTINE set_to_const

    FUNCTION array_to_2d_array(b) RESULT(a)
        REAL(KIND=wp), DIMENSION(3,3)            :: a
        REAL(KIND=wp), DIMENSION(:),  INTENT(IN) :: b

        IF ( SIZE ( b ) /= 6 ) THEN
            WRITE(*,*) ' SIZE(b)=', SIZE(b)
            CALL die ('Unreasonable array size.', 'array_to_2d_array')
        ENDIF

        a(1,1) = b(1)
        a(2,2) = b(2)
        a(3,3) = b(3)
        a(1,2) = b(4)
        a(1,3) = b(5)
        a(2,3) = b(6)
        a(2,1) = b(4)
        a(3,1) = b(5)
        a(3,2) = b(6)        

    END FUNCTION array_to_2d_array

!   This is redundant since we have array to matrix conversion in VECTORS.
    FUNCTION array_to_mat(b) RESULT (A)
        TYPE(matrix)                             :: A
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: b
!       Local array:
        REAL(KIND=wp), DIMENSION(3,3)            :: temp
        temp = .CONVERT. b
        A = temp
    END FUNCTION array_to_mat

    FUNCTION array_2d_to_1d_array(a) RESULT ( b ) 
        REAL(KIND=wp), DIMENSION(6)                :: b
         REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)  :: a

        IF ( SIZE ( a, DIM=1 ) /= 3 .OR. SIZE ( a, DIM=2 ) /= 3 ) THEN
            WRITE(*,*) SIZE ( a, DIM=1 ), SIZE ( a, DIM=2 )  
            CALL die('Unreasonable array size.', 'array_2d_to_1d_array')
        ENDIF

        b(1) = a(1,1)
        b(2) = a(2,2)
        b(3) = a(3,3)
        b(4) = a(1,2)
        b(5) = a(1,3)
        b(6) = a(2,3)
    END FUNCTION array_2d_to_1d_array

    SUBROUTINE aniso_to_array (arr, aniso_1)
        REAL(KIND=wp), DIMENSION(6), INTENT(OUT) :: arr
        TYPE(aniso),                 INTENT(IN)  :: aniso_1
        arr = aniso_1%u
    END SUBROUTINE aniso_to_array

    SUBROUTINE array_to_aniso(aniso_1, arr)
        TYPE(aniso),                 INTENT(OUT) :: aniso_1
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: arr
        aniso_1%u = arr
    END SUBROUTINE array_to_aniso

    SUBROUTINE aniso_to_aniso(aniso_1, aniso_2)
        TYPE(aniso),                 INTENT(OUT) :: aniso_1
        TYPE(aniso),                 INTENT(IN)  :: aniso_2
        aniso_1%u = aniso_2%u
    END SUBROUTINE aniso_to_aniso

    FUNCTION b_to_u( b )
        TYPE(aniso)                :: b_to_u 
        REAL(KIND=wp), INTENT(IN)  :: b
        b_to_u%u(1:3) = b
        b_to_u%u(4:6) = 0.0_wp
    END FUNCTION b_to_u

    FUNCTION u_to_b ( u )
        REAL(KIND=wp)              :: u_to_b
        TYPE(aniso),   INTENT(IN)  :: u
        u_to_b = SUM ( u%u(1:3) ) / 3.0_wp
    END FUNCTION u_to_b

    SUBROUTINE aniso_to_mat(U, aniso_1)
        TYPE(matrix), INTENT(OUT) :: U
        TYPE(aniso),  INTENT(IN)  :: aniso_1
!       Local array:
        REAL(KIND=wp), DIMENSION(3,3) :: a

        a(1,1) = aniso_1%u(1)
        a(2,2) = aniso_1%u(2)
        a(3,3) = aniso_1%u(3)

        a(1,2) = aniso_1%u(4)
        a(1,3) = aniso_1%u(5)
        a(2,3) = aniso_1%u(6)
!
        a(2,1) = aniso_1%u(4)
        a(3,1) = aniso_1%u(5)
        a(3,2) = aniso_1%u(6)

!       Use VECTORS module:
        U = a

    END SUBROUTINE aniso_to_mat
  
    FUNCTION add_u_To_u(u1, u2)
        TYPE(aniso)             :: add_u_to_u
        TYPE(aniso), INTENT(IN) :: u1
        TYPE(aniso), INTENT(IN) :: u2
        add_u_to_u%u = u1%u + u2%u
   END FUNCTION add_u_to_u
   
   FUNCTION add_u_to_b(u1, b)
       TYPE(aniso)               :: add_u_to_b
       TYPE(aniso),   INTENT(IN) :: u1
       REAL(KIND=wp), INTENT(IN) :: b
       add_u_to_b%u(1:3) = U1%u(1:3) + b 
       add_u_to_b%u(4:6) = U1%u(4:6)
   END FUNCTION add_u_to_b

   FUNCTION add_b_to_u(b, u1)
       TYPE(aniso)               :: add_b_to_u
       TYPE(anisO),   INTENT(IN) :: u1
       REAL(KIND=wp), INTENT(IN) :: b
       add_b_to_u%u(1:3) = U1%u(1:3) + b
       add_b_to_u%u(4:6) = U1%u(4:6)
   END FUNCTION add_b_to_u

   FUNCTION u_minus_u ( u1, u2 )
       TYPE(aniso)             :: u_minus_u
       TYPE(aniso), INTENT(IN) :: u1
       TYPE(aniso), INTENT(IN) :: u2
       u_minus_u%u = u1%u - u2%u
   END FUNCTION u_minus_u

   FUNCTION scale_u_by_const(u, const)
       TYPE(aniso)                :: scale_u_by_const
       TYPE(aniso),    INTENT(IN) :: u
       REAL(KIND=wp),  INTENT(IN) :: const

       scale_u_by_const%u = u%u * const

   END FUNCTION scale_u_by_const

   FUNCTION scale_const_by_u(const, u)
       TYPE(aniso)                :: scale_const_by_u
       REAL(KIND=wp),  INTENT(IN) :: const
       TYPE(aniso),    INTENT(IN) :: u

       scale_const_by_u%u = const * u%u

   END FUNCTION scale_const_by_u

   FUNCTION divide_u(u, const)
       TYPE(aniso)               :: divide_u
       TYPE(aniso),   INTENT(IN) :: u
       REAL(KIND=wp), INTENT(IN) :: const
       divide_u%u = u%u / const
   END FUNCTION divide_u

   FUNCTION gt_const(u, c)  RESULT(arr) 
       LOGICAL,     DIMENSION(3)            :: arr
       TYPE(aniso),              INTENT(IN) :: U
       REAL(KIND=wp),            INTENT(IN) :: c
         
       arr = U%u(1:3) > c
   END FUNCTION gt_const

   FUNCTION ge_const(u, c)  RESULT(arr)
       LOGICAL,     DIMENSION(3)            :: arr
       TYPE(aniso),              INTENT(IN) :: U
       REAL(KIND=wp),            INTENT(IN) :: c

       arr = U%u(1:3) >= c
   END FUNCTION ge_const

   FUNCTION eq_const(u, c)  RESULT(arr)
       LOGICAL,     DIMENSION(3)            :: arr
       TYPE(aniso),              INTENT(IN) :: U
       REAL(KIND=wp),            INTENT(IN) :: c

       arr = U%u(1:3) == c
   END FUNCTION eq_const

   FUNCTION usum3 ( u )
       REAL(KIND=wp)           :: usum3
       TYPE(aniso), INTENT(IN) :: U
       usum3 = SUM (U%u(1:3))
   END FUNCTION usum3

   FUNCTION max3val ( u )
       REAL(KIND=wp)           :: max3val
       TYPE(aniso), INTENT(IN) :: U
       max3val = MAXVAL(U%u(1:3))
   END FUNCTION max3val

   FUNCTION min3val ( u )
       REAL(KIND=wp)           :: min3val
       TYPE(aniso), INTENT(IN) :: U
       min3val = MINVAL(U%u(1:3))
   END FUNCTION min3val

   FUNCTION udet_aniso(U)
       REAL(KIND=wp)            :: udet_aniso
       TYPE(aniso),  INTENT(IN) :: U 

       udet_aniso = U%u(1) * U%u(2) * U%u(3)          &
                  - U%u(1) * U%u(6) ** 2              &
                  - U%u(2) * U%u(5) ** 2              &
                  - U%u(3) * U%u(4) ** 2              &
                  + 2.0_wp * U%u(4) * U%u(5) * U%u(6)
    END FUNCTION udet_aniso

    FUNCTION udet_aniso_using_minors ( U, minors )
        REAL(KIND=wp)                           :: udet_aniso_using_minors
        TYPE(aniso),  INTENT(IN) :: U 
        REAL(KIND=wp), DIMENSION(6), INTENT(IN) :: minors
!       Local:
        REAL(KIND=wp)                           :: temp

!       More economical version of UDET_ANISO above (if minors are available):
!        temp = U%u(1) * minors(1) - U%u(4) * minors(4) + U%u(5) * minors(5)
        temp = U%u(1) * minors(1)  + U%u(4) * minors(4) + U%u(5) * minors(5)

        IF ( temp <= 0.0_wp ) THEN
            WRITE(*,*) ' det=', temp
            WRITE(*,"(' U=', 6ES9.2)") U%u
            WRITE(*,"(' minors=', 6ES9.2)") minors
            CALL die('Possible programming error. Negative or zero determinant.', 'udet_aniso_using_minors')
        ENDIF

        udet_aniso_using_minors = temp 

    END FUNCTION udet_aniso_using_minors

    FUNCTION det_first_deriv_aniso (U) RESULT (d)
!
!       Purpose:
!       =======
!       Derivative of determinant with respect to matrix elements:
!       det|U|/dUij
!       Massive derivative calculations.
!
        REAL(KIND=wp), DIMENSION(6)   :: d
        TYPE(aniso),   INTENT(IN)     :: U
!       Local arrays:
        REAL(KIND=wp), DIMENSION(3,3) :: a

        d(1) = U%u(2) * U%u(3) - U%u(6) ** 2
        d(2) = U%u(1) * U%u(3) - U%u(5) ** 2
        d(3) = U%u(1) * U%u(2) - U%u(4) ** 2

        d(4) = -2.0_wp * ( U%u(3) * U%u(4) - U%u(5) * U%u(6) )
        d(5) = -2.0_wp * ( U%u(2) * U%u(5) - U%u(4) * U%u(6) )
        d(6) = -2.0_wp * ( U%u(1) * U%u(6) - U%u(4) * U%u(5) )

    END FUNCTION det_first_deriv_aniso

    FUNCTION det_first_deriv_selected (U, ij) RESULT (d)
!
!       Purpose:
!       =======
!       Derivative of determinant with respect to matrix elements:
!       det|U|/dUij
!       Massive derivative calculations.
!
        REAL(KIND=wp)                 :: d
        TYPE(aniso),   INTENT(IN)     :: U
        INTEGER,       INTENT(IN)     :: ij
!       Local arrays:
        REAL(KIND=wp), DIMENSION(3,3) :: a

        SELECT CASE (ij)

        CASE(1)

        d = U%u(2) * U%u(3) - U%u(6) ** 2

        CASE(2)

        d = U%u(1) * U%u(3) - U%u(5) ** 2

        CASE(3)

        d = U%u(1) * U%u(2) - U%u(4) ** 2

        CASE(4)

        d = -2.0_wp * ( U%u(3) * U%u(4) - U%u(5) * U%u(6) )

        CASE(5)

        d = -2.0_wp * ( U%u(2) * U%u(5) - U%u(4) * U%u(6) )

        CASE(6)

        d = -2.0_wp * ( U%u(1) * U%u(6) - U%u(4) * U%u(5) )

        CASE DEFAULT

        WRITE(*,*) ' ij=', ij
        CALL die('Incorrect value of IJ > 6.', 'det_first_deriv_select')

        END SELECT

    END FUNCTION det_first_deriv_selected

    FUNCTION det_second_deriv_mat(U) RESULT (dd)
        REAL(KIND=wp), DIMENSION(6,6) :: dd
        TYPE(aniso),   INTENT(IN)     :: U
!       Counters:
        INTEGER                       :: i
        INTEGER                       :: j

        dd(1,1) = 0.0_wp
        dd(1,2) = U%u(3)
        dd(1,3) = U%u(2)
        dd(1,4) = 0.0_wp
        dd(1,5) = 0.0_wp
        dd(1,6) = -2.0_wp * U%u(6)

        dd(2,2) = 0.0_wp
        dd(2,3) = U%u(1)
        dd(2,4) = 0.0_wp
        dd(2,5) = -2.0_wp * U%u(5)
        dd(2,6) = 0.0_wp

        dd(3,3) = 0.0_wp
        dd(3,4) = -2.0_wp * U%u(4)
        dd(3,5) = 0.0_wp
        dd(3,6) = 0.0_wp

        dd(4,4) = -2.0_wp * U%u(3)
        dd(4,5) =  2.0_wp * U%u(6)
        dd(4,6) =  2.0_wp * U%u(5)

        dd(5,5) = -2.0_wp * U%u(2)
        dd(5,6) =  2.0_wp * U%u(4)

        dd(6,6) = -2.0_wp * U%u(1)

!       Apply symmetry around main diagonal:
        DO i = 1, 5
            DO j = i+1, 6
                dd(j,i) = dd(i,j)
            ENDDO
        ENDDO
!        write(*,*) ' Second deriv of the det:'
!        DO i = 1,6
!            WRITE(*,"(6F10.5)") (dd(i,j),j=1,6)
!        ENDDO

    END FUNCTION det_second_deriv_mat

    FUNCTION minors_aniso(U)  RESULT(a)
        REAL(KIND=wp), DIMENSION(6)   :: a
        TYPE(aniso),  INTENT(IN)      :: U

        a(1) = U%u(2) * U%u(3) - U%u(6) ** 2
        a(2) = U%u(1) * U%u(3) - U%u(5) ** 2
        a(3) = U%u(1) * U%u(2) - U%u(4) ** 2

        a(4) = U%u(5) * U%u(6) - U%u(3) * U%u(4)
        a(5) = U%u(4) * U%u(6) - U%u(2) * U%u(5)
        a(6) = U%u(4) * U%u(5) - U%u(1) * U%u(6)

    END FUNCTION minors_aniso

    FUNCTION minors_first_deriv_test(U,ts) RESULT (d)
!
!       Purpose:
!       =======
!       Calculates first derivatives for all minors.
!       
        REAL(KIND=wp), DIMENSION(6)  :: d
        TYPE(aniso),   INTENT(IN)    :: U
        INTEGER,       INTENT(IN)    :: ts

        SELECT CASE ( ts )

        CASE(1)

!           First minor derivatives:
            d(1) =  0.0_wp
            d(2) =  U%u(3)
            d(3) =  U%u(2)
            d(4) =  0.0_wp
            d(5) =  0.0_wp
            d(6) = -2.0_wp * U%u(6)

        CASE(2)

!           Second minor derivatives:
            d(1) =  U%u(3)
            d(2) =  0.0_wp
            d(3) =  U%u(1)
            d(4) =  0.0_wp
            d(5) = -2.0_wp * U%u(5)
            d(6) =  0.0_wp

        CASE(3)

!           Third minor derivatives:
            d(1) =  U%u(2)
            d(2) =  U%u(1)
            d(3) =  0.0_wp
            d(4) = -2.0_wp * U%u(4)
            d(5) =  0.0_wp
            d(6) =  0.0_wp

        CASE(4)

!           Fourth minor derivatives:
            d(1) =  0.0_wp
            d(2) =  0.0_wp
            d(3) = -U%u(4)
            d(4) = -U%u(3)
            d(5) =  U%u(6)
            d(6) =  U%u(5)

        CASE(5)

!           Fifth minor derivatives:
            d(1) =  0.0_wp
            d(2) = -U%u(5); 
            d(3) =  0.0_wp; 
            d(4) =  U%u(6); 
            d(5) = -U%u(2); 
            d(6) =  U%u(4);

        CASE(6)

!           Sixth minor derivatives:
            d(1) = -U%u(6); 
            d(2) =  0.0_wp; 
            d(3) =  0.0_wp; 
            d(4) =  U%u(5); 
            d(5) =  U%u(4); 
            d(6) = -U%u(1)

        CASE DEFAULT

            WRITE(*,*) ' TS=', ts
            CALL die ( 'Unreasonable TS argument on input.', 'minors_first_derivatives' )

        END SELECT

    END FUNCTION minors_first_deriv_test

    FUNCTION all_minors_first_deriv_wrt_U(U,ij) RESULT (d)
!
!       Purpose:
!       =======
!       Calculates first derivatives for all minors with respect to U(ij).
!       
        REAL(KIND=wp), DIMENSION(6)  :: d
        TYPE(aniso),   INTENT(IN)    :: U
        INTEGER,       INTENT(IN)    :: ij

        SELECT CASE ( ij )

        CASE(1)

!           Derivatives of all minors with respect to U11:
            d(1) =  0.0_wp
            d(2) =  U%u(3)
            d(3) =  U%u(2)
            d(4) =  0.0_wp
            d(5) =  0.0_wp
            d(6) = -U%u(6); 


        CASE(2)

!           Derivatives of all minors with respect to U22:
            d(1) =  U%u(3)
            d(2) =  0.0_wp
            d(3) =  U%u(1)
            d(4) =  0.0_wp
            d(5) = -U%u(5); 
            d(6) =  0.0_wp; 

        CASE(3)

!           Derivatives of all minors with respect to U33:
            d(1) =  U%u(2)
            d(2) =  U%u(1)
            d(3) =  0.0_wp
            d(4) = -U%u(4)
            d(5) =  0.0_wp; 
            d(6) =  0.0_wp; 

        CASE(4)

!           Derivatives of all minors with respect to U12:
            d(1) =  0.0_wp
            d(2) =  0.0_wp
            d(3) = -2.0_wp * U%u(4)
            d(4) = -U%u(3)
            d(5) =  U%u(6); 
            d(6) =  U%u(5); 

        CASE(5)

!           Derivatives of all minors with respect to U13:
            d(1) =  0.0_wp
            d(2) = -2.0_wp * U%u(5)
            d(3) =  0.0_wp
            d(4) =  U%u(6)
            d(5) = -U%u(2); 
            d(6) =  U%u(4); 

        CASE(6)

!           Derivatives of all minors with respect to U23:
            d(1) = -2.0_wp * U%u(6)
            d(2) =  0.0_wp
            d(3) =  0.0_wp
            d(4) =  U%u(5)
            d(5) =  U%u(4);
            d(6) = -U%u(1)

        CASE DEFAULT

            WRITE(*,*) ' ij=', ij
            CALL die ( 'Unreasonable ij argument on input.', 'minors_first_derivatives' )

        END SELECT

    END FUNCTION all_minors_first_deriv_wrt_u

    FUNCTION minors_second_deriv(ts) RESULT( dd )
        REAL(KIND=wp), DIMENSION(6,6) :: dd 
        INTEGER,       INTENT(IN)     :: ts

!       Avoiding assignment of too many zeroes:
        dd = 0.0_wp

        SELECT CASE ( ts )

        CASE(1)

!       First minor second derivatives:
        dd(2,3) =  1
        dd(3,2) =  1
        dd(6,6) = -2

        CASE(2)

!       Second minor second derivatives:
        dd(1,3) =  1
        dd(3,1) =  1
        dd(5,5) = -2
        
        CASE(3)

!       Thrid minor second derivatives:
        dd(1,2) =  1
        dd(2,1) =  1
        dd(4,4) = -2

        CASE(4)

!       Fourth minor second derivatives:
        dd(3,4) = -1
        dd(4,3) = -1
        dd(5,6) =  1
        dd(6,5) =  1

        CASE(5)
   
!       Fifth minor second derivatives:
        dd(2,5) = -1
        dd(4,6) =  1
        dd(5,2) = -1
        dd(6,4) =  1

        CASE(6)

!       Sixth minor second derivatives:
        dd(1,6) = -1
        dd(4,5) =  1
        dd(5,4) =  1
        dd(6,1) = -1

        CASE DEFAULT    
        
            WRITE(*,*) ' TS=', ts
            CALL die('Programming error. Unreasonable TS argument on input.', 'minors_second_deriv')

        END SELECT

!   dd(6,6,6)
!   First index ts:
!   Second index ij:
!   Third index kl:

    END FUNCTION minors_second_deriv

    FUNCTION UINV(U) RESULT (VTS)
        TYPE(matrix)                  :: VTS
        TYPE(aniso),   INTENT(IN)     :: U
!       Local arrays:            
        REAL(KIND=wp), DIMENSION(6)   :: a
        REAL(KIND=wp), DIMENSION(3,3) :: temp
!
        a = minors_aniso (U) / udet_aniso (U)
        temp = .CONVERT. a
        
        VTS = temp 
    END FUNCTION UINV

    FUNCTION VTS_first_deriv(U, ij) RESULT(dvts)
!       should be Uinv_first_deriv
!       FIXME dvts -> 6x6 for all ij
!
!       Calculates derivatives of TSth minor
!       with respect to all aniso tensor values.
!
        REAL(KIND=wp), DIMENSION(6) :: dvts
        TYPE(aniso),   INTENT(IN)   :: U
        INTEGER,       INTENT(IN)   :: ij
!       Local arrays:
        REAL(KIND=wp), DIMENSION(6) :: MINORS
        REAL(KIND=wp), DIMENSION(6) :: VTS
        REAL(KIND=wp)               :: udet
        REAL(KIND=wp)               :: ddetUdU
        REAL(KIND=wp), DIMENSION(6) :: dumtsdu

!       Calculate DET:
        udet   = udet_aniso(U)

!       Calculate Minors:
        MINORS = minors_aniso(U)

!       Calculate inverted matrix:
        VTS    = MINORS / Udet

!       Calculate determinant first derivatives. We need just Uijth derivative...
!       FIXME
        ddetUdU = det_first_deriv_selected(U,ij)

!       All minors first derivatives against Uij:
        dUmtsdU  = all_minors_first_deriv_wrt_u(U, ij)

!       Need to introduce 6 x 6 matrix:
!       FIXME:
        dvts = dUmtsdU - VTS * ddetUdU
!
!       Final scaling according to (26):
        dvts = dvts / Udet

    END FUNCTION vts_first_deriv

    FUNCTION correct_aniso_tensor(U, ratio) RESULT ( UC )
        TYPE(aniso)                   :: UC
        TYPE(aniso),     INTENT(IN)   :: U
        REAL(KIND=wp),   INTENT(IN)   :: ratio
!       Local variables:
        REAL(KIND=wp), DIMENSION(3,3) :: UMAT
        TYPE(aniso)                   :: UCTEMP
        INTEGER                       :: info
        REAL(KIND=wp), DIMENSION(3)   :: e
        REAL(KIND=wp)                 :: old_trace
        REAL(KIND=wp)                 :: scale
        REAL(KIND=wp)                 :: UMAX 
        REAL(KIND=wp), DIMENSION(3,3) :: W
!       Counters:
        INTEGER                       :: itns
        INTEGER                       :: ja
        CHARACTER(LEN=32),      SAVE  :: srname='correct_aniso_tensor'

!       Convert to 2d array to 1d array (don't touch this please):
        UMAT = .CONVERT. U%u ! BUG corrected

!       UMAT will be destroyed and will contain eigenvectors after this call:
        CALL syevd(UMAT, e, jobz='V', UPLO='U', info=info)

!       Total reset of aniso parameters if negative ellipsoid axis is detected:
        IF ( ALL ( e <= 0.0_wp ) ) THEN
            CALL messag('All three ellipsoid axes are negative. Resetting...', srname)
            UCTEMP = .CONVERT. ( SUM ( e ) / 3.0_wp )

!           FIXME: need more special functions:
            UCTEMP%U(1:3) = 1.0_wp / ( 8.0 * pi ** 2)  ! setting Uiso equal to 1.0 A ** 2
            UC = UCTEMP
            RETURN
        ENDIF

!       Make a copy of initial ADP:
        UCTEMP = U

!       --- Protect against too large Uiso ---:
        UMAX = MAXVAL ( UCTEMP )


!       This corresponds to 3.0 * 8 * pi ** 2 ~ 240 A**2:
        IF ( UMAX > 3.0_wp ) THEN
            UCTEMP = UCTEMP * ( 3.0_wp / UMAX )
            e  = e * ( 3.0_wp / UMAX )
        ENDIF

        IF ( e(1) / e (3) > ratio ) THEN
            UC = UCTEMP
!           Ratio looks good -> return:   
            RETURN
        ENDIF

!       --- End of protection against too large Uiso ---


!       Iterate to correct ration between ellipsoid axes:
        itns = 0
        DO WHILE ( ANY (  e(1:2)  / e(3)  < ratio ) )
            IF ( itns > 3 ) THEN
                WRITE(*,"(' itns=', I3, ' current ADP=', 6F8.4, ' ratios=', 2F6.3)") &
                            itns, UCTEMP%u, e(1:2)/e(3)
                CALL warn('Too many iterations.', srname)
                EXIT
            ENDIF
            IF ( debug > 1000 ) THEN
                IF ( itns == 0 ) WRITE(*,*)
                WRITE(*,"(' itns=', I3, ' current ADP=', 6F8.4, ' ratios=', 2F6.3)") &
                            itns, UCTEMP%u, e(1:2)/e(3)
            ENDIF

!           2d array to 1d conversion:
            UMAT = .CONVERT. UCTEMP%u

!           Eigenvalue analysis, UMAT will be destroyed:
            CALL syevd(UMAT, e, jobz='V', UPLO='U', info=info)
            IF ( info /= 0 ) THEN
                WRITE(*,*) ' info=', info
                CALL warn ('O-o-o-op-s-s... Problems during eigenvector analysis...', srname)
            ENDIF

            old_trace = SUM ( e )

!           Cannot have too small trace:
            old_trace = MAX ( old_trace, 3.0_wp / (8.0_wp * pi ** 2) )  

!           Add EPS to avoid useless cycling:
            e(1) = e(3) * ( ratio + EPSILON ( 1.0_wp) )

            scale = old_trace / SUM ( e )

!           Keep trace intact:
            e = e * scale

!           Scale by columns:
            DO ja = 1, 3
                W(1:3,ja) = UMAT(1:3,ja) * e(ja)
            ENDDO
            UMAT = MATMUL ( W, TRANSPOSE ( UMAT ) ) ! UMAT recreated

!           Convert matrix to 6 values of UC:
            UCTEMP%u = .CONVERT. UMAT

            itns = itns + 1
        ENDDO

        IF ( debug > 1000 ) THEN
            WRITE(*,"(' itns=', I3, ' current ADP=', 6F8.4, ' ratios=', 2F6.3, ' RESULT')")&
            itns, UCTEMP%u, e(1:2)/e(3)
        ENDIF

!       Save final result:
        UC = UCTEMP 

    END FUNCTION correct_aniso_tensor

    SUBROUTINE check_xyz_deriv(U, x)

!       Purpose:
!       =======
!       Checks xyz first derivatives (without usage of central differences).
!
        TYPE(aniso),                 INTENT(IN)    :: U
        REAL(KIND=wp), DIMENSION(3), INTENT(INOUT) :: x
!       Local:
        REAL(KIND=wp)                              :: scale_aniso
        REAL(KIND=wp)                              :: qform
        REAL(KIND=wp)                              :: udet
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(3)                :: f1
        REAL(KIND=wp), DIMENSION(3)                :: drodx
        REAL(KIND=wp)                              :: delta
        REAL(KIND=wp), DIMENSION(6)                :: v
        REAL(KIND=wp), DIMENSION(3,3)              :: VTS
        TYPE(matrix)                               :: temp_mat
        INTEGER                                    :: i
        INTEGER                                    :: j

        udet = udet_aniso(U)
        scale_aniso = 1.0_wp / ( twopi ** 1.5 * SQRT ( udet ) )
        v    = minors_aniso(U) / Udet
        VTS  = .CONVERT. v
        temp_mat = VTS
        CALL print_mat(temp_mat)
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
        f0 = scale_aniso * EXP (-0.5_wp * qform )
        WRITE(*,*) ' f0=', f0
!       Calc dro/dx:
        drodx = -scale_aniso * MATMUL(VTS, x) *  EXP (-0.5_wp * qform )
        WRITE(*,*) ' drodx=', drodx
        WRITE(*,"(6(ES17.9))") drodx

        delta= 0.1_wp ** 9 
        DO i = 1, 3
            x(i) = x(i) + delta
            qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
            f1(i) = scale_aniso * EXP (-0.5_wp * qform )
            x(i) = x(i) - delta
        ENDDO
        WRITE(*,*) ' Please compare with:'
        WRITE(*,"(6(ES17.9))") (f1 - f0)/delta
    END SUBROUTINE check_xyz_deriv

    SUBROUTINE xyz_deriv(dRodX, U, x)
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: Drodx 
        TYPE(aniso),                 INTENT(IN)    :: U
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: x
!       Local:
        REAL(KIND=wp)                              :: scale_aniso
        REAL(KIND=wp)                              :: qform
        REAL(KIND=wp)                              :: udet
        REAL(KIND=wp), DIMENSION(6)                :: v
        REAL(KIND=wp), DIMENSION(3,3)              :: VTS

        udet = udet_aniso(U)
        scale_aniso = 1.0_wp / SQRT ( twopi ** 3 * udet ) 
        v    = minors_aniso(U) / udet
        VTS  = .CONVERT. v
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
!       Calc dro/dx:
        drodx = -scale_aniso * MATMUL(VTS, x) *  EXP (-0.5_wp * qform )
!       WRITE(*,*) ' drodx=', drodx
!       WRITE(*,"(6(ES17.9))") drodx
    END SUBROUTINE xyz_deriv

    SUBROUTINE second_xyz_deriv(dRodX, U, x)
        REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: Drodx
        TYPE(aniso),                   INTENT(IN)    :: U
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN)    :: x
!       Local:
        REAL(KIND=wp)                                :: scale_aniso
        REAL(KIND=wp)                                :: qform
        REAL(KIND=wp)                                :: udet
        REAL(KIND=wp), DIMENSION(6)                  :: v
        REAL(KIND=wp), DIMENSION(3,3)                :: VTS
        REAL(KIND=wp), DIMENSION(3,3)                :: XXT
        INTEGER                                      :: i,j
        TYPE(matrix)                                 :: temp_mat
        udet = udet_aniso(U)
        scale_aniso = 1.0_wp / SQRT ( twopi ** 3 * udet )
        v    = minors_aniso(U) / udet
        VTS  = .CONVERT. v
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
!       Calc dro/dx:
        v(1:3) = MATMUL(VTS, x)
        DO i = 1, 3
            DO j = 1,3
                XXT(i,j) = v(i) * v(j)
            ENDDO
        ENDDO
!       No minus sign here:
        drodx = scale_aniso * (-VTS + XXT) *  EXP (-0.5_wp * qform )
        WRITE(*,*) ' Second derivatives against xyz(analytical)='
!        WRITE(*,"(9(ES17.9))") drodx
        temp_mat = drodx
        CALL print_mat(temp_mat)
    END SUBROUTINE second_xyz_deriv
 
    SUBROUTINE second_U_deriv(dROdu, U, x)
        REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: dRodU
        TYPE(aniso),                   INTENT(IN)    :: U
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN)    :: x
!       Local:
        REAL(KIND=wp)                                :: scale_aniso
        REAL(KIND=wp)                                :: qform
        REAL(KIND=wp)                                :: udet
        REAL(KIND=wp), DIMENSION(6)                  :: v
        REAL(KIND=wp), DIMENSION(6)                  :: minors
!       Special arrays for transpose:
        REAL(KIND=wp), DIMENSION(6,1)                :: vdet
        REAL(KIND=wp), DIMENSION(6,1)                :: vm
!       Inverted U matrix:
        REAL(KIND=wp), DIMENSION(3,3)                :: VTS
        REAL(KIND=wp), DIMENSION(6)                  :: XXT
        REAL(KIND=wp), DIMENSION(6,6)                :: DZERO
        REAL(KIND=wp), DIMENSION(6,6,6)              :: CTS
        REAL(KIND=wp), DIMENSION(6,6)                :: D2U
        REAL(KIND=wp), DIMENSION(6,6)                :: DDT
        REAL(KIND=wp), DIMENSION(6,6)                :: VVT
        REAL(KIND=wp), DIMENSION(6,6)                :: M2D
        INTEGER                                      :: ts
        INTEGER                                      :: ij
       
        udet = udet_aniso(U)
        scale_aniso = 0.5_wp / SQRT ( twopi ** 3 * udet )
        minors = minors_aniso(U) / udet
        VTS  = .CONVERT. minors
        
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )

!       Zero order coefs:
        vdet(1:6,1) = det_first_deriv_aniso ( U ) / UDET
        VVT   = MATMUL ( vdet, TRANSPOSE(vdet) ) 
        D2U   = det_second_deriv_mat(U) / UDET
        DZERO = 1.5_wp * VVT - d2U

!       Second_order coefs:
        xxt(1) = x(1) ** 2
        xxt(2) = x(2) ** 2
        xxt(3) = x(3) ** 2
        xxt(4) = 2.0_wp * x(1) * x(2)
        xxt(5) = 2.0_wp * x(1) * x(3)
        xxt(6) = 2.0_wp * x(2) * x(3)

        CTS = 0.0_wp
        DO ts = 1, 6        
!           Derivative of all Uij against single Ts minor:
            vm(1:6,1) = minors_first_deriv_test(U, ts)
            DDT = MATMUL ( vm, TRANSPOSE(vdet) ) 
            M2D = minors_second_deriv ( ts )
            CTS(1:6,1:6,ts) =  CTS(1:6,1:6,ts) + ( minors(ts) * (-3.0_wp * VVT + d2U) & 
            + (1.0_wp / UDET) *  ( 1.5_wp * ( DDT + TRANSPOSE (DDT)  ) &
                                 - M2D ) ) 
        ENDDO

!       xxt dependent section (secoond order coefs):
        D2U = 0
        DO ts = 1, 6
            D2U = D2U + CTS(1:6,1:6,ts) * xxt(ts)
        ENDDO

!       Loop over dV(ts)/dU(ij). xxt dependent (preparing for fourth order):
        vdet = 0
        DO ij = 1, 6
            v  = VTS_first_deriv(U, ij)
            vdet(ij,1) = vdet(ij,1) + SUM ( v * xxt )
        ENDDO

!       Fourth order coefs:
        DDT = 0.5_wp * MATMUL ( vdet, TRANSPOSE ( vdet ) )

!       Final sum:
        dROdU = scale_aniso * EXP ( -0.5 * qform ) * (DZERO + D2U + DDT) 
        WRITE(*,"(' xxt=', 6F8.2)")  xxt
        WRITE(*,*) ' Final probably correct result (UU analytical):'
        DO ts = 1, 6
            WRITE(*,"(1X,6F12.6)") (dRodU(ts,ij),ij=1,6)
        ENDDO
    END SUBROUTINE second_u_deriv

    FUNCTION funani(U, x, occ)
        REAL(KIND=wp)                           :: funani
        TYPE(aniso),                 INTENT(IN) :: U
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: x
        REAL(KIND=wp),               INTENT(IN) :: occ
!       Local:
        REAL(KIND=wp)                           :: scale_aniso
        REAL(KIND=wp)                           :: qform
        REAL(KIND=wp)                           :: udet
        REAL(KIND=wp), DIMENSION(6)             :: v
        REAL(KIND=wp), DIMENSION(3,3)           :: VTS

        udet = udet_aniso(U)
        scale_aniso = 1.0_wp / SQRT ( twopi ** 3 * udet )
        v    = minors_aniso(U) / udet
        VTS  = .CONVERT. v
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
        funani = occ * scale_aniso * EXP ( -0.5_wp * qform )
    END FUNCTION funani

    FUNCTION finite_diff_coord2(U,x,occ,i,j,h,k)
        REAL(KIND=wp)                           :: finite_diff_coord2
        TYPE(aniso),                 INTENT(IN) :: U
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: x
        REAL(KIND=wp),               INTENT(IN) :: occ
        INTEGER,                     INTENT(IN) :: i
        INTEGER,                     INTENT(IN) :: j
        REAL(KIND=wp),               INTENT(IN) :: h
        REAL(KIND=wp),               INTENT(IN) :: k
!       Local:
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: y
        TYPE(aniso)                              :: Utemp
        REAL(KIND=wp)                            :: occ_temp
        REAL(KIND=wp), DIMENSION(4)              :: f
!       Local variables:
        REAL(KIND=wp), DIMENSION(2)              :: sign_h
        REAL(KIND=wp), DIMENSION(2)              :: sign_k
!       Counters:
        INTEGER                                  :: n 
        INTEGER                                  :: l 
        INTEGER                                  :: m
        DATA sign_h/1.0_wp,-1.0_wp/
        DATA sign_k/1.0_wp,-1.0_wp/

        CALL allocate_array(y, SIZE(x))

        IF ( i /= j ) THEN
            m = 0
            DO n = 1, 2
                DO l = 1, 2
                    m = m + 1
!                   Make sure we have fresh copy of y array all the time:
                    y = x
                    Utemp%u = U%u

!                   Assume the following order: x,y,z, q, U1, U2, etc.
                    IF ( i <= 3 ) THEN
                        y(i) = x(i) + sign_h(n) * h
                    ELSE IF ( i == 4 ) THEN
                        occ_temp = occ + sign_h(n) * h
                    ELSE    
                        Utemp%u(i-4) = U%u(i-4) + sign_h(n) * h
                    ENDIF

                    IF ( j <=3 ) THEN
                        y(j) = x(j) + sign_k(l) * k
                    ELSE IF ( j == 4 ) THEN
                        occ_temp = occ + sign_k(l) * k
                    ELSE
                        Utemp%u(j-4) = U%u(j-4) + sign_k(l) * k
                    ENDIF

!                   Adapt sign for summation:
                    f(m) =  sign_h(n) * sign_k(l) * funani(Utemp, y, occ_temp)

                ENDDO
            ENDDO
            CALL deallocate_array(y)
!           finite_diff_coord2 = (f(1) - f(2) - f(3) + f(4) )  / (4.0_wp * h * k):
            finite_diff_coord2 = SUM ( f ) / (4.0_wp * h * k)
        ELSE
            y = x
            Utemp = U
            occ_temp = occ

!           Increment ith parameter:
            IF ( i <= 3 ) THEN
                y(i) = x(i) + h
            ELSE IF ( i == 4 ) THEN
                occ_temp =  (occ + h) 
            ELSE 
                Utemp%u(i-4) = U%u(i-4) + h
            ENDIF

            f(1) = funani(Utemp, y, occ_temp)
!            IF ( i == 4 ) f(1) = f(1) / 2
!           Calculation function with initial params:        
            y = x
            Utemp = U
            occ_temp = occ

            IF ( occ_temp == 0.0_wp ) STOP 'ERROR in finite_diff_coord2'

            f(2) = funani(Utemp, y, occ_temp)
!            IF ( I == 4 ) f(2) = f(2)/2 
!           Refresh and decrement ith parameter:
            y = x
            Utemp = U
            occ_temp = occ

            IF ( i <=3 ) THEN
                y(i) = x(i) - h
            ELSE IF ( i == 4 ) THEN
                occ_temp =  occ - h
            ELSE 
                Utemp%u(i-4) = U%u(i-4) - h
            ENDIF

            f(3) = funani(Utemp, y, occ_temp)
!            IF ( i == 4 ) f(3) = f(3) / 2
            finite_diff_coord2 = (f(1) - 2.0_wp * f(2) + f(3)) / h ** 2            
        ENDIF
        CALL deallocate_array(y)
    END FUNCTION finite_diff_coord2

    FUNCTION finite_diff_coord(U,x,i,j,h,k)
        REAL(KIND=wp)                           :: finite_diff_coord
        TYPE(aniso),                 INTENT(IN) :: U
        REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: x
        INTEGER,                     INTENT(IN) :: i
        INTEGER,                     INTENT(IN) :: j
        REAL(KIND=wp),               INTENT(IN) :: h
        REAL(KIND=wp),               INTENT(IN) :: k
!       Local:
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: y
!       Local variables:
        REAL(KIND=wp)                            :: fpp
        REAL(KIND=wp)                            :: fpm
        REAL(KIND=wp)                            :: fmp
        REAL(KIND=wp)                            :: fmm

        CALL allocate_array(y, SIZE(x))

!       Always refresh y array:
        y = x
        y(i) = x(i) + h
        y(j) = x(j) + k
        fpp = funani(U, y, 1.0_wp)

        y = x
        y(i) = x(i) + h
        y(j) = x(j) - k
        fpm = funani(U, y, 1.0_wp)

        y = x
        y(i) = x(i) - h
        y(j) = x(j) + k
        fmp = funani(U, y, 1.0_wp)
        y = x
        y(i) = x(i) - h
        y(j) = x(j) - k
        fmm = funani(U, y,1.0_wp)

        CALL deallocate_array(y)

        finite_diff_coord = (fpp - fpm - fmp + fmm) / (4.0_wp * h * k)

    END FUNCTION finite_diff_coord

    SUBROUTINE check_u_deriv(U, x)
        TYPE(aniso),                 INTENT(INOUT) :: U
        REAL(KIND=wp), DIMENSION(3), INTENT(INOUT) :: x
!       Local:
        REAL(KIND=wp)                           :: scale_aniso
        REAL(KIND=wp)                           :: qform
        REAL(KIND=wp)                           :: udet
        REAL(KIND=wp)                           :: f0
        REAL(KIND=wp), DIMENSION(6)             :: f1
        REAL(KIND=wp), DIMENSION(6)             :: drodU
        REAL(KIND=wp)                           :: delta
        REAL(KIND=wp), DIMENSION(6)             :: v
        REAL(KIND=wp), DIMENSION(6)             :: vtsd
        REAL(KIND=wp), DIMENSION(3,3)           :: VTS
        REAL(KIND=wp), DIMENSION(6)             :: xxt
        INTEGER                                 :: i
        INTEGER                                 :: ij 
        udet = udet_aniso(U)
        scale_aniso = 1.0_wp / ( twopi ** 1.5 * SQRT ( udet ) )
        v    = minors_aniso(U) / Udet
        VTS  = .CONVERT. v
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
        f0 = scale_aniso * EXP (-0.5_wp * qform )
        WRITE(*,*) ' f0=', f0
!       Calc dro/dU:
!       Scaled derivatives of the determinant :
        v = det_first_deriv_aniso ( U ) / UDET
        xxt(1) = x(1) ** 2 
        xxt(2) = x(2) ** 2
        xxt(3) = x(3) ** 2
        xxt(4) = 2.0_wp * x(1) * x(2)
        xxt(5) = 2.0_wp * x(1) * x(3)
        xxt(6) = 2.0_wp * x(2) * x(3)
 
!       Loop over dV(ts)/dU(ij)
        DO ij = 1, 6
            vtsd  = VTS_first_deriv(U, ij)
            v(ij) = v(ij) + SUM ( vtsd * xxt )
        ENDDO
        dRodU = -0.5_wp * scale_aniso * v *  EXP (-0.5_wp * qform ) 
        WRITE(*,*) ' drodU='
        WRITE(*,"(6(ES17.9))") drodU

        delta= 0.1_wp ** 8 
        
        DO i = 1, 6
            U%u(i) = U%u(i) + delta
            udet = udet_aniso(U)
            scale_aniso = 1.0_wp / ( twopi ** 1.5 * SQRT ( udet ) )
            v    = minors_aniso(U) / Udet
!           Convert 1d array to 3x3 2d array:
            VTS  = .CONVERT. v
            qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
            f1(i) = scale_aniso * EXP (-0.5_wp * qform )
            U%u(i) = U%u(i) - delta
        ENDDO
        WRITE(*,*) ' Please compare with:'
        WRITE(*,"(6(ES17.9))") (f1 - f0)/delta
        WRITE(*,"(6(ES17.9))") ((f1 - f0)/delta) / drodu
    END SUBROUTINE check_u_deriv

    SUBROUTINE u_deriv(dRoDu, U, x)
        REAL(KIND=wp), DIMENSION(6)                :: drodU
        TYPE(aniso),                 INTENT(INOUT) :: U
        REAL(KIND=wp), DIMENSION(3), INTENT(INOUT) :: x
!       Local:
        REAL(KIND=wp)                              :: scale_aniso
        REAL(KIND=wp)                              :: qform
        REAL(KIND=wp)                              :: udet
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(6)                :: f1
        REAL(KIND=wp), DIMENSION(6)                :: v
        REAL(KIND=wp), DIMENSION(6)                :: vtsd
        REAL(KIND=wp), DIMENSION(3,3)              :: VTS
        REAL(KIND=wp), DIMENSION(6)                :: xxt
        INTEGER                                    :: ij 
!
!       It's pretty obvious now -- we should start with minors calcns:
        v    = minors_aniso(U)
!        udet = udet_aniso(U)
!       Not very aesthetically pleasing but economical: 
        udet = U%u(1)*v(1) - U%u(4) * v(4) + U%u(5) * v(5)
        v    = v / Udet
        scale_aniso = 1.0_wp / ( twopi ** 1.5 * SQRT ( udet ) )
!       COnvert to matrix form (3x3 array):
        VTS  = .CONVERT. v
!       Ready to calculate eigenvalues:
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
        f0 = scale_aniso * EXP (-0.5_wp * qform )
!       Calc dro/dU:
!       Scaled derivatives of the determinant :
        v = det_first_deriv_aniso ( U ) / UDET
        xxt(1) = x(1) ** 2 
        xxt(2) = x(2) ** 2
        xxt(3) = x(3) ** 2
        xxt(4) = 2.0_wp * x(1) * x(2)
        xxt(5) = 2.0_wp * x(1) * x(3)
        xxt(6) = 2.0_wp * x(2) * x(3)
 
!       Loop over dV(ts)/dU(ij)
        DO ij = 1, 6
            vtsd  = VTS_first_deriv(U, ij)
            v(ij) = v(ij) + SUM ( vtsd * xxt )
        ENDDO
        dRodU = -0.5_wp * scale_aniso * v *  EXP (-0.5_wp * qform ) 
!        WRITE(*,*) ' drodU='
!        WRITE(*,"(6(ES17.9))") drodU
    END SUBROUTINE u_deriv

    SUBROUTINE xu_cross_deriv(dRoDu, U, x)
        REAL(KIND=wp), DIMENSION(3,6)              :: drodU
        TYPE(aniso),                 INTENT(INOUT) :: U
        REAL(KIND=wp), DIMENSION(3), INTENT(INOUT) :: x
!       Local:
        REAL(KIND=wp)                              :: scale_aniso
        REAL(KIND=wp)                              :: qform
        REAL(KIND=wp)                              :: udet
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(6)                :: f1
        REAL(KIND=wp), DIMENSION(6)                :: v
        REAL(KIND=wp), DIMENSION(6,1)              :: vk
        REAL(KIND=wp), DIMENSION(3,1)              :: vm
        REAL(KIND=wp), DIMENSION(6)                :: vtsd
        REAL(KIND=wp), DIMENSION(3,3)              :: VTS
        REAL(KIND=wp), DIMENSION(6)                :: xxt
        INTEGER                                    :: ij 
        INTEGER                                    :: ts
!
!       It's pretty obvious now -- we should start with minors calcns:
        v    = minors_aniso(U)
!        udet = udet_aniso(U)
!       Not very aesthetically pleasing but economical: 
        udet = U%u(1)*v(1) - U%u(4) * v(4) + U%u(5) * v(5)
        v    = v / Udet
        scale_aniso = 1.0_wp / ( twopi ** 1.5 * SQRT ( udet ) )
!       COnvert to matrix form (3x3 array):
        VTS  = .CONVERT. v
!       Ready to calculate eigenvalues:
        qform = DOT_PRODUCT( x, MATMUL ( VTS, x) )
        f0 = scale_aniso * EXP (-0.5_wp * qform )

        xxt(1) = x(1) ** 2 
        xxt(2) = x(2) ** 2
        xxt(3) = x(3) ** 2
        xxt(4) = 2.0_wp * x(1) * x(2)
        xxt(5) = 2.0_wp * x(1) * x(3)
        xxt(6) = 2.0_wp * x(2) * x(3)
 
!       Calc dro/dU:
!       Scaled derivatives of the determinant :
        v = det_first_deriv_aniso ( U ) / UDET
!       Loop over dV(ts)/dU(ij), don't attempt to set v=0 here, please:
        DO ij = 1, 6
            vtsd  = VTS_first_deriv(U, ij)
            v(ij) = v(ij) + SUM ( vtsd * xxt )
        ENDDO
        vk(1:6,1) = v
        vm(1:3,1) = MATMUL ( VTS, x )
        dRodU = 0.5_wp * MATMUL ( vm, TRANSPOSE ( vk ) )
        
        DO ij = 1,6
            vtsd =  VTS_first_deriv(U, ij)
            VTS  = .CONVERT. vtsd
            dRodU(1:3,ij) = dRodU(1:3,ij) - MATMUL ( VTS, x ) 
        ENDDO
        dRodU = scale_aniso * dRodU *  EXP (-0.5_wp * qform ) 
        WRITE(*,*) ' drodU='
        DO ts = 1,3
            WRITE(*,"(6(F10.6))") (drodU(ts,ij),ij=1,6)
        ENDDO
    END SUBROUTINE xu_cross_deriv

    SUBROUTINE all_second_deriv ( dRodP, U, duvwort, occ_i, occ_j)
!
!       Purpose:
!       =======
!       Calculates various second derivatives of electron denstiy of anisotropically moving atom.
!     
!       Garib made a small error claiming he was calculating minors of the matrix
!       In fact, those are cofactors since sign has (-1)**(i+j) has been included in
!       the calculation of minors.
!
        REAL(KIND=wp), DIMENSION(10,10), INTENT(INOUT) :: dRodP
        TYPE(aniso),                     INTENT(IN)    :: U
        REAL(KIND=wp), DIMENSION(3),     INTENT(IN)    :: duvwort
        REAL(KIND=wp),                   INTENT(IN)    :: occ_i
        REAL(KIND=wp),                   INTENT(IN)    :: occ_j
!       Local variables and arrays:
        REAL(KIND=wp)                                  :: scale_aniso
        REAL(KIND=wp)                                  :: qform
        REAL(KIND=wp)                                  :: udet
        REAL(KIND=wp)                                  :: udet_minors
        REAL(KIND=wp), DIMENSION(6)                    :: minors
!       Scaled minors:
        REAL(KIND=wp), DIMENSION(6)                    :: minors_divided_by_det
!       Special arrays for transpose:
        REAL(KIND=wp), DIMENSION(6)                    :: udet_deriv
        REAL(KIND=wp), DIMENSION(6)                    :: det_scaled_deriv
        REAL(KIND=wp), DIMENSION(6,1)                  :: minors_deriv
!       Inverted U matrix:
        REAL(KIND=wp), DIMENSION(3,3)                  :: VTS
!       Deriv. wrt aniso U:
        REAL(KIND=wp), DIMENSION(6,6)                  :: VTSD
        REAL(KIND=wp), DIMENSION(3,3)                  :: VTSD_ROW
        REAL(KIND=wp), DIMENSION(3)                    :: VTS_duvwort
        REAL(KIND=wp), DIMENSION(3,1)                  :: vtemp
        REAL(KIND=wp), DIMENSION(6,1)                  :: utemp
        REAL(KIND=wp), DIMENSION(1,6)                  :: wtemp
!       Normal VXXT matrix:
        REAL(KIND=wp), DIMENSION(3,3)                  :: VXXT        
!       Mixed xxt vector with last three components doubled:
        REAL(KIND=wp), DIMENSION(6)                    :: wxxt                               
!       Zero order matrix:
        REAL(KIND=wp), DIMENSION(6,6)                  :: DZERO 
!       Garibs 2nd order coefficient matirx:
        REAL(KIND=wp), DIMENSION(6,6,6)                :: CTS       
        REAL(KIND=wp), DIMENSION(6,6)                  :: DSECOND
        REAL(KIND=wp), DIMENSION(6,6)                  :: DSECOND_DERIV
!       Fourth order:
        REAL(KIND=wp), DIMENSION(6,6)                  :: DFOURTH
        REAL(KIND=wp), DIMENSION(6,6)                  :: MINORS_SEC_DER
!       Additional matrices:
        REAL(KIND=wp), DIMENSION(6,6)                  :: DDT
        REAL(KIND=wp), DIMENSION(6,6)                  :: VVT
!       Temp:
        REAL(KIND=wp), DIMENSION(3,6)                  :: DRODU
        INTEGER                                        :: ts
        INTEGER                                        :: ij
        INTEGER                                        :: k
        INTEGER                                        :: l
        INTEGER                                        :: ind
        TYPE(matrix)                                   :: temp_matrix

!       General outer loop section:
        minors = minors_aniso ( U )
!        udet_minors   = udet_aniso_using_minors ( U, minors )
        udet = udet_aniso ( u )
        WRITE(*,*) ' udet =', udet
        scale_aniso = 1.0_wp / SQRT ( twopi ** 3 * udet )
        minors_divided_by_det = minors / udet
        VTS = .CONVERT. minors_divided_by_det
!       U derivative section (maybe not needed in some cases):
        DO k = 1, 6
            VTSD (k,1:6) = VTS_first_deriv ( U, k )
        ENDDO
        
!       Former V vector (in some cases):        
        det_scaled_deriv = det_first_deriv_aniso ( U ) / udet ! vdet
        utemp(1:6,1) = det_scaled_deriv
        wtemp(1,1:6) = det_scaled_deriv
        VVT = MATMUL ( utemp, TRANSPOSE ( utemp ) )           ! VVT
        DSECOND_DERIV = det_second_deriv_mat (U) / udet       ! D2U
        DZERO = 1.5_wp * VVT - DSECOND_DERIV

!       Formula (24):
        CTS = 0.0_wp
        DO ts = 1, 6
!           Derivative of all U(ij) against single TS minor:
            utemp(1:6,1) = minors_first_deriv_test(U, ts)
            DDT          = MATMUL ( utemp, wtemp )
            MINORS_SEC_DER = minors_second_deriv ( ts )
            CTS(1:6,1:6,ts) =  CTS(1:6,1:6,ts) + ( minors_divided_by_det(ts)           &
                            * (-3.0_wp * VVT + DSECOND_DERIV)                          &
                            + (1.0_wp / UDET) *  ( 1.5_wp * ( DDT + TRANSPOSE (DDT)  ) &
                            - MINORS_SEC_DER ) )
        ENDDO

!       Distance depending XYZ section. This is usefull in inner loop for other derivs:
        VTS_duvwort = MATMUL ( VTS, duvwort ) 
        qform = DOT_PRODUCT ( duvwort, VTS_duvwort )

!       -------------------------------------------------------------------------      
!       XX area
!       VXXT matrix:
        vtemp(1:3,1) = VTS_duvwort
        VXXT = MATMUL ( vtemp, TRANSPOSE ( vtemp ) )
        dRodP(1:3,1:3) = occ_i * occ_j * scale_aniso * (-VTS + VXXT) * EXP ( -0.5_wp * qform )

!       -------------------------------------------------------------------------      
!       Omitting XYZB section for a while... XYZ deriv. based on aniso value while B is B.
!       Generally we can calculate this using chain rule.
!       Convert B to U
!       We can calculate XYZU (see below) as normal, and get (1:3,1:6) matrix
!       Then convert it to (1:3,1) matrix applying chain rule.
!       -------------------------------------------------------------------------      
!       XYZB area

!       -------------------------------------------------------------------------      
!       XO area
!
!       XYZ against Oj section (basically looks like a XYZ gradient):
!       Calc dro/dx: This should be 5 after all:
        ind = 4
        DRODP(1:3,ind) = -occ_i * scale_aniso * VTS_duvwort * EXP ( -0.5_wp * qform )
        DRODP(ind,1:3) = -occ_j * scale_aniso * VTS_duvwort * EXP ( -0.5_wp * qform )
!       -------------------------------------------------------------------------      
!       XYZ against U section (XYZU) XU aread:
        wxxt(1) = duvwort(1) ** 2
        wxxt(2) = duvwort(2) ** 2
        wxxt(3) = duvwort(3) ** 2
        wxxt(4) = 2.0_wp * duvwort(1) * duvwort(2)
        wxxt(5) = 2.0_wp * duvwort(1) * duvwort(3)
        wxxt(6) = 2.0_wp * duvwort(2) * duvwort(3)

        utemp(1:6,1) = det_scaled_deriv(1:6) + MATMUL ( VTSD(1:6,1:6), wxxt )
        vtemp(1:3,1) = VTS_duvwort

!       Initialise:
        dRoDU =  0.5 * MATMUL ( vtemp, TRANSPOSE ( utemp ) )
!       Accumulate: 
        DO ij = 1, 6
!           Convert to squared matrix:
            VTSD_ROW = .CONVERT. VTSD(ij,1:6)
            DRODU(1:3,ij) = DRODU(1:3,ij) - MATMUL ( VTSD_ROW, duvwort )
        ENDDO

        ind = 5
!       No need to apply 0.5 factor (see xu_cross_deriv)
        DRODP(1:3,ind:ind+5) = occ_i * occ_j * scale_aniso * EXP ( -0.5_wp * qform ) * DRODU
        DRODP(ind:ind+5,1:3) = occ_i * occ_j * scale_aniso * EXP ( -0.5_wp * qform ) * TRANSPOSE ( DRODU )

!       -------------------------------------------------------------------------
!       OX area has been computed in XO area... Skipping...
!       -------------------------------------------------------------------------
!       OB area... Skipping...
!       -------------------------------------------------------------------------
!       OO area...

        DRODP(4,4) =  scale_aniso * EXP ( -0.5_wp * qform )
!       
!       -------------------------------------------------------------------------
!       Ocnncy Oi agaist U ( QU area):
        
!       Loop over dV(ts)/dU(ij)
!        DO ij = 1, 6
!            vtsd  = VTS_first_deriv(U, ij)
!            v(ij) = v(ij) + SUM ( vtsd * xxt )
!        ENDDO
        utemp(1:6,1) = det_scaled_deriv + MATMUL ( VTSD, wxxt )
        
        DRODP(4,5:10) = -0.5_wp * occ_i * scale_aniso * utemp(1:6,1) *  EXP (-0.5_wp * qform )
        DRODP(5:10,4) = -0.5_wp * occ_j * scale_aniso * utemp(1:6,1) * EXP (-0.5_wp * qform )
        
!       -------------------------------------------------------------------------
!       Skipping Ux area... has been computed in xU section.
!       Skipping UB area... No clue yet... Need two atoms...
!       Skippin  QU area... has been computed in OU section. Just above
!       -------------------------------------------------------------------------
!       Infamous UU area:
        
!       xxt dependent section (second order coefs):
        DSECOND = 0.0_wp
        DO ts = 1, 6
            DSECOND = DSECOND + CTS(1:6,1:6,ts) * wxxt(ts)
        ENDDO 
        utemp(1:6,1) =  MATMUL ( VTSD, wxxt )

!       Fourth order coefs:
        DFOURTH = 0.5_wp * MATMUL ( utemp, TRANSPOSE ( utemp ) )

!       Final sum, factor of 0.5 is necessay to correct common scale aniso:
        dROdp(5:10,5:10) = ( 0.5_wp * occ_i * occ_j * scale_aniso * EXP ( -0.5_wp * qform ) ) &
                         * ( DZERO + DSECOND + DFOURTH )

        WRITE(*,*) ' all_second_deriv:'
        DO k = 1, 10
            WRITE(*,"(10F10.6)") (drodp(k,l),l=1,10)
        ENDDO

END SUBROUTINE all_second_deriv

END MODULE aniso_manip
