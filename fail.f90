MODULE fail 
USE string_manip
IMPLICIT NONE
PRIVATE
PUBLIC :: bruteptf_failure
PUBLIC :: die
PUBLIC :: warn
PUBLIC :: messag
CONTAINS
    SUBROUTINE die(message, sub)
        CHARACTER(len=*), INTENT(IN) :: message
        CHARACTER(len=*), INTENT(IN) :: sub
!       Local
        CHARACTER(len=LEN_TRIM(sub)) :: temp

        WRITE(*,"(' DIE> ', '################################')")
        WRITE(*,"(' DIE> ', '################################')")
        WRITE(*,"(' DIE> ', '###          ERROR           ###')")
        WRITE(*,"(' DIE> ', '################################')")
        WRITE(*,"(' DIE> ', '################################')")
        WRITE(*,"(' . . . . . . . . . . . . . . . . . . . . .')")
        temp = sub(1:LEN_TRIM(sub))
        CALL ucase( temp )
        WRITE(*,1000) temp, TRIM(message)
        CALL EXIT(1)
1000    FORMAT (1X,A,'> ERROR:',1X,A,/,1X,'DIE> Bailing out...')
    END SUBROUTINE die

    SUBROUTINE warn (message, sub)
        CHARACTER(len=*), INTENT(IN) :: message
        CHARACTER(len=*), INTENT(IN) :: sub
!
        CHARACTER(len=LEN_TRIM(sub)) :: temp
        temp = sub
        CALL ucase( temp )
        WRITE(*,1000) temp, TRIM(message)
1000    FORMAT (1X,A,'> WARNING : ', A)
    END SUBROUTINE warn

    SUBROUTINE messag (message, sub)
        CHARACTER(len=*), INTENT(IN) :: message
        CHARACTER(len=*), INTENT(IN) :: sub
!
        CHARACTER(len=LEN_TRIM(sub)) :: temp
        temp = TRIM(sub)
        CALL ucase( temp )
        WRITE(*,1000) temp, TRIM(message)
1000    FORMAT (1X,A,'> ', A)
    END SUBROUTINE messag
    
    SUBROUTINE bruteptf_failure ( sub )
        CHARACTER(len=*), INTENT(IN) :: sub
!
        CHARACTER(len=LEN_TRIM(sub)) :: temp
        temp = TRIM(sub)
        CALL ucase ( temp )
        CALL messag('All solutions have been rejected due to self_symmetry packing considerations.', temp)
        CALL messag('Use less strict packing criteria or try to divide your molecule', temp)
        CALL messag('into separate domains (they may have moved), remove suspicios loops,', temp)
        CALL messag('N-terminus, C-terminus, etc.', temp)
        CALL messag('Another possibility is to get better experimental phases and run again.', temp)
        CALL messag(' ', 'bruteptf_failure')
        CALL die   ('Zero number of solutions left... ', 'bruteptf_failure')

    END SUBROUTINE bruteptf_failure
END MODULE fail 
