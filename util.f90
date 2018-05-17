MODULE util
USE constants
USE fail
!USE iso_varying_string
USE select_kinds
USE string_manip
USE vectors
IMPLICIT NONE
PRIVATE
PUBLIC :: ASSIGNMENT ( = )
PUBLIC :: allocate_array
PUBLIC :: ch_is_digit
PUBLIC :: check_file_name
PUBLIC :: check_directory_name
PUBLIC :: deallocate_array
PUBLIC :: get_next_io_unit
PUBLIC :: init_random_number
PUBLIC :: int_to_c
PUBLIC :: OPERATOR ( .NORM. )
PUBLIC :: OPERATOR ( .NUMERIC. )
PUBLIC :: OPERATOR ( .EQS. )
PUBLIC :: is_numeric
PUBLIC :: ugtenv2
PUBLIC :: swap
INTERFACE OPERATOR ( .NUMERIC. )
    MODULE PROCEDURE get_numeric
END INTERFACE

INTERFACE OPERATOR ( .EQS. )
    MODULE PROCEDURE word_equal_char
END INTERFACE

INTERFACE int_to_c
    MODULE PROCEDURE int4_to_c
    MODULE PROCEDURE int8_to_c
END INTERFACE

INTERFACE swap
    MODULE PROCEDURE swap_int
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_1d_wp_array
    MODULE PROCEDURE allocate_1d_wp_array_eb
    MODULE PROCEDURE allocate_2d_wp_array
    MODULE PROCEDURE allocate_3d_wp_array
    MODULE PROCEDURE allocate_2d_sp_array
    MODULE PROCEDURE allocate_1d_complex_wp_array
    MODULE PROCEDURE allocate_2d_complex_wp_array
    MODULE PROCEDURE allocate_1d_vector_int_array
    MODULE PROCEDURE allocate_2d_vector_int_array
    MODULE PROCEDURE allocate_1d_vector_array
    MODULE PROCEDURE allocate_2d_vector_array
    MODULE PROCEDURE allocate_1d_integer_array
    MODULE PROCEDURE allocate_1d_integer_array_eb
    MODULE PROCEDURE allocate_1d_integer_eb_array
    MODULE PROCEDURE allocate_2d_integer_array
    MODULE PROCEDURE allocate_2d_integer_array_eb
    MODULE PROCEDURE allocate_2d_int_wb_array
    MODULE PROCEDURE allocate_3d_int_wb_array
    MODULE PROCEDURE allocate_3d_int_fb_array
    MODULE PROCEDURE allocate_3d_int_grid_array
    MODULE PROCEDURE allocate_1d_matrix_int_array
    MODULE PROCEDURE allocate_1d_matrix_array
    MODULE PROCEDURE allocate_1d_sp_array
    MODULE PROCEDURE allocate_1d_complex_sp_array
    MODULE PROCEDURE allocate_1d_logical_bt_array
    MODULE PROCEDURE allocate_2d_logical_bt_array
    MODULE PROCEDURE allocate_1d_character_array
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_1d_wp_array
    MODULE PROCEDURE deallocate_2d_wp_array
    MODULE PROCEDURE deallocate_3d_wp_array
    MODULE PROCEDURE deallocate_2d_sp_array
    MODULE PROCEDURE deallocate_1d_complex_wp_array
    MODULE PROCEDURE deallocate_2d_complex_wp_array
    MODULE PROCEDURE deallocate_1d_vector_int_array
    MODULE PROCEDURE deallocate_2d_vector_int_array
    MODULE PROCEDURE deallocate_1d_vector_array
    MODULE PROCEDURE deallocate_2d_vector_array
    MODULE PROCEDURE deallocate_1d_integer_array
    MODULE PROCEDURE deallocate_1d_integer_eb_array
    MODULE PROCEDURE deallocate_2d_integer_array
    MODULE PROCEDURE deallocate_2d_int_wb_array
    MODULE PROCEDURE deallocate_3d_int_wb_array
    MODULE PROCEDURE deallocate_3d_int_fb_array
    MODULE PROCEDURE deallocate_1d_matrix_int_array
    MODULE PROCEDURE deallocate_1d_matrix_array
    MODULE PROCEDURE deallocate_1d_sp_array
    MODULE PROCEDURE deallocate_1d_complex_sp_array
    MODULE PROCEDURE deallocate_1d_logical_bt_array
    MODULE PROCEDURE deallocate_2d_logical_bt_array
    MODULE PROCEDURE deallocate_1d_character_array
END INTERFACE

INTERFACE OPERATOR ( .NORM. )
    MODULE PROCEDURE normalize_wp_array
END INTERFACE

CONTAINS
    SUBROUTINE allocate_1d_wp_array(array, n, array_name)
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                  INTENT(IN)           :: n
        CHARACTER(len=*),                         INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                        :: istat

        IF ( ALLOCATED(array) ) THEN
            DEALLOCATE(array, STAT=istat)
            IF (istat /= 0) THEN
                CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_wp_array' )
            ENDIF
        ENDIF

        ALLOCATE(array(n), STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_wp_array' )
        ELSE
            IF ( PRESENT(array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_wp_array

    SUBROUTINE allocate_1d_wp_array_eb(array, n, array_name)
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER(KIND=eb),                         INTENT(IN)           :: n
        CHARACTER(len=*),                         INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                        :: istat

        IF ( ALLOCATED(array) ) THEN
            DEALLOCATE(array, STAT=istat)
            IF (istat /= 0) THEN
                CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_wp_array_eb' )
            ENDIF
        ENDIF

        ALLOCATE(array(n), STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_wp_array_eb' )
        ELSE
            IF ( PRESENT(array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_wp_array_eb' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_wp_array_eb


    SUBROUTINE allocate_1d_sp_array ( array, n, array_name )
        REAL(KIND=sp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                  INTENT(IN)           :: n
        CHARACTER(len=*),                         INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                        :: istat

        IF ( ALLOCATED(array) ) THEN
            DEALLOCATE(array, STAT=istat)
            IF ( istat /= 0 ) THEN
                CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_sp_array' )
            ENDIF
        ENDIF

        ALLOCATE(array(n), STAT=istat)
        IF ( istat /= 0 ) THEN
            CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_sp_array' )
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_sp_array

    SUBROUTINE deallocate_1d_wp_array ( array, array_name )
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                         INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                        :: istat
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_wp_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_1d_wp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_wp_array

    SUBROUTINE deallocate_1d_sp_array ( array, array_name )
        REAL(KIND=sp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                         INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                        :: istat
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_sp_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_1d_sp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_sp_array

    SUBROUTINE allocate_2d_wp_array ( array, n, m, array_name )
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                    INTENT(IN)           :: n
        INTEGER,                                    INTENT(IN)           :: m
        CHARACTER(len=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                          :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT=istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_wp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_wp_array' )
                ENDIF
            ENDIF
        ENDIF
         
        ALLOCATE ( array(n, m), STAT=istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_wp_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_wp_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_wp_array

    SUBROUTINE deallocate_2d_wp_array ( array, array_name )
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                          :: istat
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_wp_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_wp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_wp_array

    SUBROUTINE allocate_2d_sp_array ( array, n, m, array_name )
        REAL(KIND=sp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                    INTENT(IN)           :: n
        INTEGER,                                    INTENT(IN)           :: m
        CHARACTER(len=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                          :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_sp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_sp_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_sp_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_sp_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_sp_array

    SUBROUTINE deallocate_2d_sp_array ( array, array_name )
        REAL(KIND=sp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                           INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                          :: istat
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_sp_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_sp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_sp_array
    
    SUBROUTINE allocate_3d_wp_array ( array, n, m, p, array_name )
        REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)          :: array
        INTEGER,                                      INTENT(IN)             :: n
        INTEGER,                                      INTENT(IN)             :: m
        INTEGER,                                      INTENT(IN)             :: p
        CHARACTER(len=*),                             INTENT(IN), OPTIONAL   :: array_name
!       Local vars:
        INTEGER                                                              :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_3d_wp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'3-D'//' array.', 'allocate_3d_wp_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(0:n-1, 0:m-1, 0:p-1), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
               CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_3d_wp_array' )
            ELSE
               CALL die('Failed to allocate '//'3-D'//' array.', 'allocate_3d_wp_array' )
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_3d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_3d_wp_array

    SUBROUTINE deallocate_3d_wp_array ( array, array_name )
        REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                             INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                            :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_3d_wp_array')
                ELSE
                    CALL die('Failed to deallocate '//'3-D'//' array.', 'deallocate_3d_wp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_3d_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_3d_wp_array

!   vector_int routines
    SUBROUTINE allocate_1d_vector_int_array ( array, n, array_name )
        TYPE(vector_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_vector_int_array' )
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_vector_int_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name ) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_vector_int_array' )
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_vector_int_array' )
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_vector_int_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_vector_int_array

    SUBROUTINE deallocate_1d_vector_int_array ( array, array_name )
        TYPE(vector_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)         :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL  ::  array_name
!       Local
        INTEGER                                                            :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_vector_int_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_vector_int_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_vector_int_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_vector_int_array

    SUBROUTINE allocate_2d_vector_int_array ( array, n, m, array_name )
        TYPE(vector_int), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_vector_int_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_vector_int_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_vector_int_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_vector_int_array' )
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_vector_int_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_vector_int_array

    SUBROUTINE deallocate_2d_vector_int_array ( array, array_name )
        TYPE(vector_int), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_vector_int_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_vector_int_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_vector_int_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_vector_int_array

    SUBROUTINE allocate_1d_vector_array ( array, n, array_name )
        TYPE(vector),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_vector_array' )
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_vector_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_vector_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_1d_vector_array' )
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_vector_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_vector_array

    SUBROUTINE allocate_2d_vector_array ( array, n, m, array_name )
        TYPE(vector),     DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m 
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

!       Start with checks
        IF ( n <= 0 .OR. m <= 0 ) THEN
            CALL die('Programming error. n or m is equal to zero. Allocation halted.', 'allocate_2d_vector_array')
        ENDIF

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_vector_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_vector_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n,m), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_vector_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_vector_array' )
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_vector_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_vector_array

    SUBROUTINE deallocate_1d_vector_array ( array, array_name )
        TYPE(vector),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_vector_int_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_1d_vector_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_vector_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_vector_array

    SUBROUTINE deallocate_2d_vector_array ( array, array_name )
        TYPE(vector),     DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_vector_int_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_vector_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_vector_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_vector_array

    SUBROUTINE allocate_1d_complex_wp_array ( array, n, array_name )
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_complex_wp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_complex_wp_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_complex_wp_array' )
            ELSE
                CALL die('Failed to allocate '//'1-D'//' array.', 'allocate_1d_complex_wp_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_complex_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_complex_wp_array

    SUBROUTINE allocate_1d_complex_sp_array ( array, n, array_name )
        COMPLEX(KIND=sp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_complex_sp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_complex_sp_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_complex_sp_array' )
            ELSE
                CALL die('Failed to allocate '//'1-D'//' array.', 'allocate_1d_complex_sp_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_complex_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_complex_sp_array

    SUBROUTINE deallocate_1d_complex_wp_array ( array, array_name )
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_complex_wp_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_complex_wp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_complex_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_complex_wp_array

    SUBROUTINE deallocate_1d_complex_sp_array ( array, array_name )
        COMPLEX(KIND=sp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_complex_sp_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_complex_sp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_complex_sp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_complex_sp_array

    SUBROUTINE allocate_2d_complex_wp_array ( array, n, m, array_name )
        COMPLEX(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_complex_wp_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_complex_wp_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_complex_wp_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_complex_wp_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_vector_int_array' )
            ENDIF
        ENDIF

    END SUBROUTINE allocate_2d_complex_wp_array

    SUBROUTINE deallocate_2d_complex_wp_array ( array, array_name )
        COMPLEX(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_complex_wp_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_complex_wp_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_complex_wp_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_complex_wp_array

    SUBROUTINE allocate_1d_logical_bt_array ( array, n, array_name )
        LOGICAL(KIND=1),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_logical_bt_array' )
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_logical_bt_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_1d_logical_bt_array' )
            ELSE
                CALL die('Failed to allocate '//'1-D'//' array.', 'allocate_1d_logical_bt_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_logcial_bt_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_logical_bt_array

    SUBROUTINE deallocate_1d_logical_bt_array ( array, array_name )
        LOGICAL(KIND=1),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_logical_bt_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_logical_bt_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_logical_bt_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_logical_bt_array

    SUBROUTINE allocate_2d_logical_bt_array ( array, n, m, array_name )
        LOGICAL(KIND=1),  DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m 
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_logical_bt_array' )
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_logical_bt_array' )
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate '//TRIM(array_name)//' array.', 'allocate_2d_logical_bt_array' )
            ELSE
                CALL die('Failed to allocate '//'2-D'//' array.', 'allocate_2d_logical_bt_array' )
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_2d_logcial_bt_array' )
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_logical_bt_array

    SUBROUTINE deallocate_2d_logical_bt_array ( array, array_name )
        LOGICAL(KIND=1),  DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_logical_bt_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_logical_bt_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_logical_bt_array' )
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_logical_bt_array

    SUBROUTINE allocate_1d_integer_array ( array, n, array_name )
        INTEGER(KIND=fb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat
 
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_integer_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_integer_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_1d_integer_array')
            ENDIF
        ELSE 
            CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_integer_array')     
        ENDIF
    END SUBROUTINE allocate_1d_integer_array

    SUBROUTINE deallocate_1d_integer_array ( array, array_name )
        INTEGER(KIND=fb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:  
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            stat:IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_integer_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_integer_array')
                ENDIF
            ELSE stat
                IF ( PRESENT(array_name) ) THEN
                    CALL messag('Array '//TRIM(array_name)//' has been successfully deallocated.', &
                    'deallocate_1d_integer_array')
                ENDIF
            ENDIF stat
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_integer_array')
            ENDIF
        ENDIF
    
    END SUBROUTINE deallocate_1d_integer_array

    SUBROUTINE allocate_1d_integer_array_eb ( array, n, array_name )
        INTEGER(KIND=fb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER(KIND=eb),                            INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat
 
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_integer_array_eb')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_integer_array_eb')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_1d_integer_array_eb')
            ENDIF
        ELSE 
            CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_integer_array_eb')     
        ENDIF
    END SUBROUTINE allocate_1d_integer_array_eb

    SUBROUTINE allocate_1d_integer_eb_array ( array, n, array_name )
        INTEGER(KIND=eb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER(KIND=eb),                            INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat
 
        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_integer_eb_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_integer_eb_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_1d_integer_eb_array')
            ENDIF
        ELSE 
            CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_integer_eb_array')     
        ENDIF
    END SUBROUTINE allocate_1d_integer_eb_array

    SUBROUTINE deallocate_1d_integer_eb_array ( array, array_name )
        INTEGER(KIND=eb), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            stat:IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_integer_eb_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_integer_eb_array')
                ENDIF
            ELSE stat
                IF ( PRESENT(array_name) ) THEN
                    CALL messag('Array '//TRIM(array_name)//' has been successfully deallocated.', &
                    'deallocate_1d_integer_eb_array')
                ENDIF
            ENDIF stat
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_integer_eb_array')
            ENDIF
        ENDIF
    
    END SUBROUTINE deallocate_1d_integer_eb_array

    SUBROUTINE allocate_2d_integer_array ( array, n, m, array_name )
        INTEGER(KIND=fb), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m 
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_integer_array')
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_integer_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_2d_integer_array')
            ENDIF
        ELSE
            CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_2d_integer_array')
        ENDIF
    END SUBROUTINE allocate_2d_integer_array

    SUBROUTINE allocate_2d_integer_array_eb ( array, n, m, array_name )
        INTEGER(KIND=fb), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER(KIND=eb),                              INTENT(IN)           :: n
        INTEGER(KIND=eb),                              INTENT(IN)           :: m
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_integer_array_eb')
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_integer_array_eb')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_2d_integer_array_eb')
            ENDIF
        ELSE
            CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_2d_integer_array_eb')
        ENDIF
    END SUBROUTINE allocate_2d_integer_array_eb

    SUBROUTINE deallocate_2d_integer_array ( array, array_name )
        INTEGER(KIND=fb), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_integer_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_integer_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_integer_array')
            ENDIF
        ENDIF

    END SUBROUTINE deallocate_2d_integer_array

    SUBROUTINE allocate_2d_int_wb_array ( array, n, m, array_name )
        INTEGER(KIND=wb), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                       INTENT(IN)           :: n
        INTEGER,                                       INTENT(IN)           :: m
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_2d_int_wb_array')
                ELSE
                    CALL die('Failed to reallocate '//'2-D'//' array.', 'allocate_2d_int_wb_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m), STAT = istat )
        IF ( istat == 0 ) THEN
            IF( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_2d_int_wb_array')
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_2d_int_wb_array')
            ELSE
                CALL die('Failed to allocate memory for '//'2-D'//' array.', 'allocate_2d_int_wb_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_2d_int_wb_array

    SUBROUTINE deallocate_2d_int_wb_array ( array, array_name )
        INTEGER(KIND=wb), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                              INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                             :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_2d_int_wb_array')
                ELSE
                    CALL die('Failed to deallocate '//'2-D'//' array.', 'deallocate_2d_int_wb_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_2d_int_wb_array')
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_2d_int_wb_array

    SUBROUTINE allocate_3d_int_wb_array ( array, n, m, p, array_name )
        INTEGER(KIND=wb), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                         INTENT(IN)           :: n
        INTEGER,                                         INTENT(IN)           :: m
        INTEGER,                                         INTENT(IN)           :: p 
        CHARACTER(len=*),                                INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                               :: istat

        IF ( n < 1 .OR. m < 1 .OR. p < 1 ) THEN
            CALL die('Programming error. Some of the array dimension are zero.', 'allocate_3d_int_wb_array')
        ENDIF

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_3d_int_wb_array')
                ELSE
                    CALL die('Failed to reallocate '//'3-D'//' array.', 'allocate_3d_int_wb_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n, m, p), STAT = istat )
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_3d_int_wb_array')
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_3d_int_wb_array')
            ELSE
                CALL die('Failed to allocate memory for '//'3-D'//' array.', 'allocate_3d_int_wb_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_3d_int_wb_array

    SUBROUTINE allocate_3d_int_fb_array ( array, n, m, p, array_name, start )
        INTEGER(KIND=fb), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                         INTENT(IN)           :: n
        INTEGER,                                         INTENT(IN)           :: m
        INTEGER,                                         INTENT(IN)           :: p
        CHARACTER(len=*),                                INTENT(IN), OPTIONAL :: array_name
        INTEGER,                                         INTENT(IN), OPTIONAL :: start                    
!       Local vars:
        INTEGER                                                               :: istat

        IF ( n < 1 .OR. m < 1 .OR. p < 1 ) THEN
            CALL die('Programming error. Some of the array dimension are zero.', 'allocate_3d_int_wb_array')
        ENDIF

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_3d_int_wb_array')
                ELSE
                    CALL die('Failed to reallocate '//'3-D'//' array.', 'allocate_3d_int_wb_array')
                ENDIF
            ENDIF
        ENDIF

        IF ( PRESENT(start) ) THEN
            IF ( start == 1 ) THEN
                ALLOCATE ( array(n, m, p), STAT = istat )
            ELSE
                ALLOCATE ( array(0:n-1, 0:m-1, 0:p-1), STAT = istat )
            ENDIF
        ELSE
            ALLOCATE ( array(0:n-1, 0:m-1, 0:p-1), STAT = istat )
        ENDIF
        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_3d_int_wb_array')
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_3d_int_wb_array')
            ELSE
                CALL die('Failed to allocate memory for '//'3-D'//' array.', 'allocate_3d_int_wb_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_3d_int_fb_array

    SUBROUTINE allocate_3d_int_grid_array ( array, grid, array_name, start )
        INTEGER(KIND=fb), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,          DIMENSION(:),                  INTENT(IN)           :: grid 
        CHARACTER(len=*),                                INTENT(IN), OPTIONAL :: array_name
        INTEGER,                                         INTENT(IN), OPTIONAL :: start
!       Local vars:
        INTEGER                                                               :: n
        INTEGER                                                               :: m 
        INTEGER                                                               :: p
        INTEGER                                                               :: istat
       
        n = grid(1)
        m = grid(2)
        p = grid(3)

        IF ( n < 1 .OR. m < 1 .OR. p < 1 ) THEN
            CALL die('Programming error. Some of the array dimension are zero.', 'allocate_3d_int_grid_array')
        ENDIF

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_3d_int_grid_array')
                ELSE
                    CALL die('Failed to reallocate '//'3-D'//' array.', 'allocate_3d_int_grid_array')
                ENDIF
            ENDIF
        ENDIF

        IF ( PRESENT(start) ) THEN
            IF ( start == 1 ) THEN
                ALLOCATE ( array(n, m, p), STAT = istat )
            ELSE
                ALLOCATE ( array(0:n-1, 0:m-1, 0:p-1), STAT = istat )
            ENDIF
        ELSE
            ALLOCATE ( array(0:n-1, 0:m-1, 0:p-1), STAT = istat )
        ENDIF

        IF ( istat == 0 ) THEN
            IF ( PRESENT(array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been successfully allocated.', 'allocate_3d_int_grid_array')
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_3d_int_grid_array')
            ELSE
                CALL die('Failed to allocate memory for '//'3-D'//' array.', 'allocate_3d_int_grid_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_3d_int_grid_array

    SUBROUTINE deallocate_3d_int_wb_array ( array, array_name )
        INTEGER(KIND=wb), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                                INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                               :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_3d_int_wb_array')
                ELSE
                    CALL die('Failed to deallocate '//'3-D'//' array.', 'deallocate_3d_int_wb_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_3d_int_wb_array')
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_3d_int_wb_array

    SUBROUTINE deallocate_3d_int_fb_array ( array, array_name )
        INTEGER(KIND=fb), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                                INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                               :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_3d_int_fb_array')
                ELSE
                    CALL die('Failed to deallocate '//'3-D'//' array.', 'deallocate_3d_int_fb_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_3d_int_fb_array')
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_3d_int_fb_array

    SUBROUTINE allocate_1d_matrix_int_array ( array, n, array_name )
        TYPE(matrix_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT ( array_name ) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_matrix_int_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_matrix_int_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT( array_name ) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_matrix_int_array')
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_matrix_int_array')
            ENDIF
        ELSE
            IF ( PRESENT (array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_matrix_int_array')
            ENDIF
        ENDIF
    END SUBROUTINE allocate_1d_matrix_int_array

    SUBROUTINE deallocate_1d_matrix_int_array ( array, array_name )
        TYPE(matrix_int), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat
 
        IF ( ALLOCATED(array) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_matrix_int_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_matrix_int_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_matrix_int_array')
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_matrix_int_array     

    SUBROUTINE allocate_1d_matrix_array ( array, n, array_name )
        TYPE(matrix),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT ( array_name ) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_matrix_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_matrix_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT( array_name ) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_matrix_array')
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_matrix_array')
            ENDIF
        ELSE
            IF ( PRESENT (array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_matrix_array')
            ENDIF
        ENDIF

    END SUBROUTINE allocate_1d_matrix_array

    SUBROUTINE deallocate_1d_matrix_array ( array, array_name )
        TYPE(matrix),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED(array) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT(array_name) ) THEN
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_matrix_array')
                ELSE
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_matrix_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT(array_name) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_matrix_array')
            ENDIF
        ENDIF
    END SUBROUTINE deallocate_1d_matrix_array

    SUBROUTINE allocate_1d_character_array ( array, n, array_name )
        CHARACTER(len=*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        INTEGER,                                     INTENT(IN)           :: n
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT ( array_name ) ) THEN
                    CALL die('Failed to reallocate '//TRIM(array_name)//' array.', 'allocate_1d_character_array')
                ELSE
                    CALL die('Failed to reallocate '//'1-D'//' array.', 'allocate_1d_character_array')
                ENDIF
            ENDIF
        ENDIF

        ALLOCATE ( array(n), STAT = istat )
        IF ( istat /= 0 ) THEN
            IF ( PRESENT( array_name ) ) THEN
                CALL die('Failed to allocate memory for '//TRIM(array_name)//' array.', 'allocate_1d_character_array')
            ELSE
                CALL die('Failed to allocate memory for '//'1-D'//' array.', 'allocate_1d_character_array')
            ENDIF
        ELSE
            IF ( PRESENT (array_name) ) THEN
                CALL messag('Array '//TRIM(array_name)//' has been allocated successfully.', 'allocate_1d_character_array')
            ENDIF
        ENDIF

    END SUBROUTINE allocate_1d_character_array

    SUBROUTINE deallocate_1d_character_array ( array, array_name )
        CHARACTER(len=*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: array
        CHARACTER(len=*),                            INTENT(IN), OPTIONAL :: array_name
!       Local vars:
        INTEGER                                                           :: istat

        IF ( ALLOCATED ( array ) ) THEN
            DEALLOCATE ( array, STAT = istat )
            IF ( istat /= 0 ) THEN
                IF ( PRESENT ( array_name ) ) THEN
                    WRITE(*,*) SIZE(array)
                    CALL die('Failed to deallocate '//TRIM(array_name)//' array.', 'deallocate_1d_character_array')
                ELSE
                    WRITE(*,*) SIZE(array), array
                    CALL die('Failed to deallocate '//'1-D'//' array.', 'deallocate_1d_character_array')
                ENDIF
            ENDIF
        ELSE
            IF ( PRESENT ( array_name ) ) THEN
                CALL warn('Array '//TRIM(array_name)//' has been deallocated already.', 'deallocate_1d_character_array')
            ENDIF
        ENDIF

    END SUBROUTINE deallocate_1d_character_array

    FUNCTION int4_to_c ( i )
        CHARACTER(len=12)            :: int4_to_c 
        INTEGER,          INTENT(IN) :: i
!       Local vars:
        CHARACTER(len=12)            :: char12

        WRITE(char12, '(I12)') i
        int4_to_c = ADJUSTL ( char12 )
    END FUNCTION int4_to_c

    FUNCTION int8_to_c ( i )
        CHARACTER(len=24)            :: int8_to_c
        INTEGER(KIND=eb), INTENT(IN) :: i
!       Local vars:
        CHARACTER(len=24)            :: char12

        WRITE(char12, '(I24)') i
        int8_to_c = ADJUSTL ( char12 )
    END FUNCTION int8_to_c

    FUNCTION normalize_wp_array( name ) RESULT(new)
        REAL(KIND=wp), DIMENSION(:),           INTENT(IN) :: name
        REAL(KIND=wp), DIMENSION(1:SIZE(name))            :: new

        new = name / SQRT ( DOT_PRODUCT ( name, name ) )
    END FUNCTION normalize_wp_array

    FUNCTION get_next_io_unit () RESULT (next)
!
!       Purpose:
!       =======
!       Finds next available file input/output unit
!
!       Date:         Programmer:        Description of changes:
!       ====          ==========         ======================
!       Oct 2005      Strokopytov B.     Lifted from E. Akin
!                                        "Object-Oriented Programming
!                                         via Fortran 90/95"
!      
!       Oct 2005          -"-            Cosmetic changes
!
        INTEGER            :: next
        INTEGER, PARAMETER :: min_unit  = 10
        INTEGER, PARAMETER :: max_unit  = 999
        INTEGER, SAVE      :: last_unit = 0     ! initialize
        INTEGER            :: count             ! number of failures
        LOGICAL            :: open              ! file status

        count = 0 ; next = min_unit - 1
!       Check next in line:
        my_last_unit:IF ( last_unit > 0 ) THEN
            next = last_unit + 1
            INQUIRE(UNIT=next, OPENED=open)
!           Found it:
            IF ( .NOT. open ) last_unit = next
            RETURN
        ELSE
!           Loop through allowed units:
            unit_numbers:DO
                next =next + 1
                INQUIRE (UNIT=next, OPENED=open)
                IF ( .NOT. open ) THEN
!                   Found it:
                    last_unit = next        
                    EXIT unit_numbers
                ENDIF

!               Attempt reset 3 times:
                reset:IF ( next == max_unit ) THEN
                    last_unit = 0
                    count     = count + 1
                    IF ( count <= 3 ) next = min_unit - 1
                ENDIF reset
              
                abort:IF ( next > max_unit ) THEN
                    CALL die('max unit exceeded.', 'get_next_io_unit')
                ENDIF abort
            ENDDO unit_numbers
        ENDIF my_last_unit

    END FUNCTION get_next_io_unit

    FUNCTION word_equal_char ( word, my_char )
        LOGICAL                      :: word_equal_char
        CHARACTER(LEN=*), INTENT(IN) :: word
        CHARACTER(LEN=*), INTENT(IN) :: my_char
!       Local vars:
        CHARACTER(LEN=LEN(word))     :: temp

        temp = ADJUSTL ( word )

!       If word is shorter than my_char they cannot be equal:
        IF ( LEN ( temp ) < LEN ( my_char ) ) THEN
            word_equal_char = .FALSE.
        ELSE
!           ALL first LEN(my_char) characters must be equal:
            word_equal_char = ( temp(1:LEN ( my_char )) == my_char )
        ENDIF
    END FUNCTION word_equal_char

    FUNCTION is_numeric ( word ) 
!
!       Purpose:
!       =======
!       Checks whether word can be converted to a double precision number.
!
!       FIXME: this function does not work reliably yet.
!       Need CORELS like stuff.
!       The number may contain only digits + '.','+','-','E','D'
!
!       Date:          Programmer:           History of changes:
!       ====           ==========            ==================
!       Aug 2005       B.Strokopytov         Original code
!       Jan 2008       B.Strokopytov         Added additional protection
!                                            against alphabetic characters.
!
!        Note:
!        ====
!        No explicit precaution against occurence of several dots, pluses or minuses.
!        Ideally this problem could be solved with regular expressions package.
!
!
        LOGICAL                      :: is_numeric
        CHARACTER(LEN=*), INTENT(IN) :: word
!       Local vars:
        CHARACTER(LEN=LEN(word))     :: temp
        REAL(KIND=wp)                :: dummy
        INTEGER                      :: istat
        CHARACTER(LEN=1)             :: ch
        CHARACTER(LEN=15),      SAVE :: char15 = '0123456789-+DE.'
!       Counters:
        INTEGER                      :: i
        INTEGER                      :: kexp
        

        temp = ADJUSTL ( word )
        is_numeric = .TRUE.

!       Simple check first:
        IF ( ch_is_alpha ( temp(1:1) ) ) THEN
            is_numeric = .FALSE.
            RETURN
        ENDIF

!       Check string for non-numeric characters:
        kexp = 0
        DO i = 1, LEN_TRIM ( temp )
            ch = temp(i:i)
            IF ( INDEX ( char15, ch ) == 0 ) THEN
                is_numeric = .FALSE.
                RETURN
            ENDIF

            IF ( ch_is_alpha ( ch ) ) THEN
                kexp = kexp + 1
                IF ( kexp > 1 ) THEN
                    is_numeric = .FALSE.
                    RETURN
                ENDIF
            ENDIF
        ENDDO

!       This our last resort to check:
        READ ( temp, *, IOSTAT=istat ) dummy
!        WRITE(*,*) ' istat=',istat, ' dummy=',dummy
        is_numeric = ( istat == 0 )
    END FUNCTION is_numeric

    FUNCTION get_numeric ( word ) RESULT ( numeric )
        REAL(KIND=wp)                    :: numeric
        CHARACTER(LEN=*),     INTENT(IN) :: word
!       Local vars:
        CHARACTER(LEN=LEN(word))         :: temp
        REAL(KIND=wp)                    :: dummy
        INTEGER                          :: istat
        
        temp = ADJUSTL ( word )
        READ ( temp, *, IOSTAT = istat ) dummy
        IF ( istat == 0 ) THEN
            numeric = dummy
        ELSE
            CALL die ( 'Syntax error when converting '//TRIM(temp)//' to double precision', 'get_numeric')
        ENDIF
    END FUNCTION get_numeric
  
    FUNCTION ch_is_alpha ( c )
!
!       Purpose:
!       =======
!       CH_IS_ALPHA returns TRUE if C is an alphabetic character.
!
!
!        Date:                Author:              History of changes:
!        ====                 ======               ==================
!        05 August  1999      John Burkardt        Original code
!        23 January 2008      Boris Strokopytov    Heavily modified
!
!        Parameters:
!        ==========
!        Input, character C, a character to check.
!
!        Output, logical CH_IS_ALPHA is TRUE if C is an alphabetic character.
!

        CHARACTER(LEN=1), INTENT(IN) :: c
        LOGICAL                      :: ch_is_alpha

        ch_is_alpha = ( LLE ( 'a', c ) .AND. LLE ( c, 'z' ) ) .OR. &
                      ( LLE ( 'A', c ) .AND. LLE ( c, 'Z' ) )

    END FUNCTION ch_is_alpha

    FUNCTION ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!

    character c
    logical ch_is_digit

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
    else
        ch_is_digit = .false.
    end if

    return
    END FUNCTION ch_is_digit

    SUBROUTINE init_random_number ( user_seed )
        INTEGER, DIMENSION(:), OPTIONAL    :: user_seed
!       Local variables:
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        INTEGER                            :: seed_size
        INTEGER                            :: istat
        INTEGER, DIMENSION(4)              :: default_seed

!       Larger numbers are probably not recommended:
        default_seed = (/123456789, 987654321, 123456789, 987654321/)
        IF ( debug > 100 ) THEN
            WRITE(*,*) ' default_seed=', default_seed
        ENDIF

!       Figure out seed size:
        CALL RANDOM_SEED ( SIZE=seed_size )

!       Allocate array of appropriate size:
        ALLOCATE ( seed(seed_size), STAT=istat )
        IF ( istat /= 0 ) THEN
            CALL die ('Failed to allocate SEED array.', 'init_random_number' )
        ENDIF

        CALL messag ( 'Initializing random number generator with seed size= '//TRIM(int_to_c(seed_size)),&
                      'init_random_number')

!       Check seed size:
        IF ( seed_size > SIZE ( default_seed ) ) THEN
            CALL die ( 'Unexpected seed size... Are you on 128-bit machine?',&
                       'init_random_number')
        ENDIF
        
        IF ( .NOT. PRESENT ( user_seed ) ) THEN
            seed(1:seed_size) = default_seed(1:seed_size)
        ELSE

!           Check user seed:
            IF ( SIZE ( user_seed ) < seed_size ) THEN
                WRITE(*,*) ' user_seed_size=', SIZE ( user_seed )
                CALL die ('USER_SEED array is too small.', 'init_random_number')
            ENDIF

            IF ( ANY ( user_seed(1:seed_size) == 0 ) ) THEN
                WRITE(*,*) ' user_seed=', user_seed
                CALL die ('USER_SEED array contains zero entries...', 'init_random_number')
            ENDIF

!           It should be OK now to initialise:
            seed(1:seed_size) = user_seed(1:seed_size)
        ENDIF

        CALL RANDOM_SEED ( put=seed(1:seed_size) )
        WRITE (*,"(' INIT_RANDOM_NUMBER> ', 'Random number generator has been &
       &initialised using the following array:', 8I12)" ) seed(1:seed_size)

    END SUBROUTINE init_random_number

    SUBROUTINE check_file_name ( file_name, exists )
        CHARACTER(LEN=*), INTENT(IN)    :: file_name
        LOGICAL                         :: exists
!       Local variables:
        CHARACTER(LEN=32),         SAVE :: srname='check_file_name'

        exists = .TRUE.

!       Check than length of file name is reasonable:
        IF ( LEN_TRIM ( file_name ) == 0 ) THEN
            exists = .FALSE.
            RETURN
        ELSE IF ( LEN_TRIM ( file_name ) > 255 ) THEN
            CALL die('Too long a file name: '//TRIM ( file_name ), srname)
        ENDIF
        
!       Name has been read but we need to check whether such file really exists:
        INQUIRE ( file = file_name, exist = exists )

    END SUBROUTINE check_file_name

    SUBROUTINE check_directory_name ( file_name, exists )
        CHARACTER(LEN=*), INTENT(IN)  :: file_name
        LOGICAL,          INTENT(OUT) :: exists
!       Local variables:
        CHARACTER(LEN=32),       SAVE :: srname='check_directory_name'
        INTEGER                       :: j
        CHARACTER(LEN=256)            :: file_dir
        
        exists = .TRUE.

!       Get directory name in Unix searching backwards:
        j = INDEX ( file_name, '/', .TRUE.)

        IF ( j > 0 ) THEN

            file_dir = file_name(1:j)//'.'

!           Need to figure out whether this suitable for other compilers than Intel:
            INQUIRE ( directory = TRIM ( file_dir ), exist = exists )
            IF ( .NOT. exists ) THEN
                CALL messag('Make sure to check directory name for the output file: '//TRIM ( file_name ), srname)
!               Hard stop:
                CALL die('No such directory name: '//TRIM ( file_dir )//' Please correct this...', srname)
            ENDIF
        ENDIF

    END SUBROUTINE  check_directory_name

    SUBROUTINE ugtenv2 ( keyword, char132 )
        CHARACTER(LEN=*), INTENT(IN)    :: keyword
        CHARACTER(LEN=*), INTENT(INOUT) :: char132
!       Local variables:
        CHARACTER(LEN=512)              :: arg
        INTEGER                         :: n
!       Counter:
        INTEGER                         :: i
        INTEGER                         :: iargc

        n = iargc()
        char132 = ' '
        DO i = 1, n
            CALL getarg ( i, arg )
            CALL ucase  ( arg )
            IF ( TRIM ( keyword ) == TRIM ( arg ) ) THEN
                IF ( i + 1 <= n ) THEN
                    CALL getarg ( i + 1, arg )
                    char132 = TRIM ( ADJUSTL ( arg ) )
                    EXIT
                ENDIF
            ENDIF
        ENDDO

    END SUBROUTINE ugtenv2
 
    SUBROUTINE swap_int (i, j)
        INTEGER, INTENT(INOUT) :: i
        INTEGER, INTENT(INOUT) :: j
!       Local variables:
        INTEGER                :: temp
        temp = i
        i = j
        j = temp
    END SUBROUTINE swap_int

END MODULE util
