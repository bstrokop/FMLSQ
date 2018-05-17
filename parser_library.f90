MODULE parser_library 
USE constants
USE fail 
!USE iso_varying_string
USE util
IMPLICIT NONE
PRIVATE
PUBLIC :: elimbl
PUBLIC :: mysplit
CONTAINS
    SUBROUTINE elimbl(line)
!
!      Purpose:
!      =======
!      Removes blank spaces INSIDE character string
!
!      Record of revisions
!          Date      Programmer
!          ====      ==========
!          2003      B.Strokopytov
!
!       List of parameters:
        CHARACTER(LEN=*)         :: line
!       Local variables:
        CHARACTER(LEN=LEN(line)) :: temp
!       Counters:
        INTEGER                  :: i
        INTEGER                  :: j

        j = 0
        loop: DO i = 1, LEN ( line )
            elim: IF ( line(i:i) /= ' ' ) THEN
                j         = j + 1
                temp(j:j) = line(i:i)
            ENDIF elim
        ENDDO loop

        line(1:j) = temp(1:j)

!       Clean up the end of line
        IF (j < LEN ( line ) ) line(j + 1:LEN ( line ) ) = ' '

    END SUBROUTINE elimbl
 
    SUBROUTINE mysplit ( string, words )
!
!        Purpose:
!        =======
!        Splits line into several arbitrary length strings (type varying_string)
!        This as of now (Oct 2005) probably still not a standard F95
!        Uses blank spaces and quotes as split points
!        (Currently no distinction between quote types)
!
!        Thanks to authors of module `iso_varying string'
!
!
!        Date:                       Programmer:     Description of changes
!        ====                        ==========      ==============================
!        Sep 2003                    Strokopytov B.  initial version
!        Sep 2005                    Strokopytov B.  iso_varying_string module added
!        Oct 2005                    Strokopytov B.  words has now INTENT(INOUT). INTENT(OUT) probably
!                                                    is a bug which has to do with realloaction of
!                                                    `words' array
!
!        Oct 2007                    B.Strokopytov   iso_varying_string has been removed
!                                                    due to problems with XLF, GFORTRAN compilers
!
         CHARACTER(LEN=*),                            INTENT(IN)    :: string
         CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: words
!        Local variables:
         INTEGER,          DIMENSION(:), ALLOCATABLE                :: start
         INTEGER,          DIMENSION(:), ALLOCATABLE                :: finis
         CHARACTER(LEN=LEN(string) )                                :: line
         CHARACTER(LEN=2), PARAMETER                                :: quote_types = '''"'
         LOGICAL                                                    :: quote
!        Counters:
         INTEGER                                                    :: i
         INTEGER                                                    :: nwords
         INTEGER                                                    :: nquotes

!        Check whether string on input is not empty:
         IF ( LEN_TRIM ( string ) == 0 ) THEN

!            If so allocate zero length words array and return:
             CALL allocate_array ( words, 0 )
             RETURN
         ENDIF

!        Make a copy:
         line = string

!        Initialize:
         quote   = .FALSE.
         nwords  = 0
         nquotes = 0

!        Estimate number of words:
         DO i = 1, LEN ( line )
             IF ( INDEX ( quote_types, line(i:i) ) > 0 ) THEN 
                 nquotes = nquotes + 1
                 IF ( MOD ( nquotes, 2 ) == 1 ) THEN
                     nwords = nwords + 1
                     quote  = .TRUE.
                 ELSE
                     quote  = .FALSE.
                 ENDIF
             ENDIF
             IF ( .NOT. quote ) THEN
                 IF ( strt(line, i ) ) nwords = nwords  + 1
             ENDIF
         ENDDO

!        Word count finished:
         CALL allocate_array(start, nwords)
         CALL allocate_array(finis, nwords)

!        Reinitialize:
         quote   = .FALSE.
         nwords  = 0
         nquotes = 0
!
         search: DO i = 1, LEN ( line )

!            Find quotation marks:
             outerif: IF ( INDEX ( quote_types, line(i:i) ) > 0 ) THEN
                 nquotes = nquotes + 1
             
!                If odd quote  - we are "inside"
!                If even quote - we are "outside"
                 IF ( MOD (nquotes, 2) == 1 ) THEN
                     nwords = nwords + 1
                     IF ( nwords > SIZE ( start ) ) THEN
                         WRITE(*,*) nwords, SIZE(start)
                         CALL die('Programming error. nwords > SIZE(start)...', 'mysplit')
                     ENDIF
                        start(nwords) = i + 1
                        quote         = .TRUE.
                    ELSE
                        finis(nwords) = i - 1
                        quote         = .FALSE.
                    ENDIF

               ENDIF outerif

               notquote: IF ( .NOT. quote ) THEN
                   strtword: IF ( strt( line, i ) ) THEN
                       nwords        = nwords + 1
                       start(nwords) = i
                   ENDIF strtword

                   fnshword: IF ( fnsh ( line, i ) ) THEN
                       finis(nwords) = i
                   ENDIF fnshword
               ENDIF notquote
        ENDDO search
      
!       Check consistency of quotes:
        IF  ( MOD  ( nquotes, 2 ) == 1 ) THEN
            CALL die ( 'Syntax error. Number of quotes on input is not even.', 'mysplit' )
        ENDIF

!       Allocate words array:
        CALL allocate_array ( words, nwords )

!       Final count:
        nw: DO i = 1, nwords
            IF ( debug >= 200 ) WRITE(*,"(' start, finis: ',2I8)") start(i),finis(i)
            words(i) = line(start(i):finis(i))

!           Blank out read words for further checks:
            line( start(i):finis(i) ) = ' '
        ENDDO nw

!       Check if any uninterpreted characters left:
        DO i = 1, LEN_TRIM ( line )
            IF ( .NOT. ( line(i:i) == quote_types(1:1) .OR. line(i:i) == quote_types(2:2) .OR. line(i:i)==' ' ) ) THEN
                CALL die ( 'Syntax error. Uninterpreted characters left : '//TRIM(line), 'mysplit' )
            ENDIF   
        ENDDO
!        
        CALL deallocate_array ( start )
        CALL deallocate_array ( finis )

     END SUBROUTINE mysplit

     FUNCTION strt(line, i)
!
!        Purpose:
!        =======
!        Detects start of non-blank sequence of characters
!
!        Date:                Programmer:         Version:
!        ====                 ==========          =======  
!        Sep 2003             Strokopytov B.      Original code
!
         LOGICAL                      :: strt
         CHARACTER(LEN=*), INTENT(IN) :: line
         INTEGER,          INTENT(IN) :: i
!        Local variables:
         CHARACTER(LEN=3), PARAMETER  :: empty=' ''"'
         INTEGER                      :: i1
!
         strt = .FALSE.
         i1   = i - 1
! 
         IF ( INDEX ( empty, line(i:i) ) <= 0 ) THEN

!            If we are here then we're dealing with non-blank character at position "i":
             IF ( i1 > 0 ) THEN
                 strt = INDEX ( empty, line(i1:i1) ) > 0
             ELSE
                 strt = .TRUE.
             ENDIF

         ENDIF
    END FUNCTION strt

    FUNCTION fnsh(line, i)
!
!       Purpose:
!       =======
!       Detects end of non-blank sequence of characters
!
!       Date:                Programmer:         Version:
!       ====                 ==========          =======
!       Sep 2003             Strokopytov B.      Original code
!       
        LOGICAL                      :: fnsh
        CHARACTER(LEN=*), INTENT(IN) :: line
        INTEGER,          INTENT(IN) :: i
!       Local parameters and variables:
        CHARACTER(LEN=3), PARAMETER  :: empty=' ''"'
        INTEGER                      :: i1

        fnsh = .FALSE.
        i1   = i + 1
!
        IF ( INDEX ( empty, line(i:i) ) <= 0 ) THEN

!           If we are here we're dealing with non-blank character at position "i":
            IF ( i < LEN ( line ) ) THEN
                fnsh = INDEX ( empty, line(i1:i1) ) > 0 
            ELSE
                fnsh = .TRUE.
            ENDIF

        ENDIF
    END FUNCTION fnsh

END MODULE parser_library
