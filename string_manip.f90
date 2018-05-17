MODULE string_manip 
IMPLICIT NONE
PUBLIC :: ucase
CONTAINS
    SUBROUTINE ucase(string)
!
!         Purpose:
!         =======
!         Convert string to uppercase in machine independent way.
!         Should give no problems whether it's ASCII or EBCDIC
!         sequence.
!
!         Record of revisions:
!             Date       Programmer
!             ====       =========
!            09/08/03    Lifted from Chapman
!
!

!       Declare calling parameters:
        CHARACTER(len=*), INTENT(INOUT) :: string

!       Declare local variables:
        INTEGER :: i
        INTEGER :: length
!
        length = LEN ( string )

!       Now shift lowercase letters to uppercase in machine-independent way.
        DO i = 1, length
            IF ( LGE(string(i:i),'a') .AND. LLE(string(i:i),'z') ) THEN
                string(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
            ENDIF
        END DO

    END SUBROUTINE ucase

END MODULE string_manip

