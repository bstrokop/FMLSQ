MODULE ucards
USE constants
USE fail
!USE iso_varying_string
USE parser_library
USE select_kinds
USE string_manip
USE util
USE vectors
IMPLICIT NONE
!
!   Purpose:
!   =======
!   Structure for storing user cards information.
!
!   Note:
!   ====
!   TOO GENERAL AT THIS MOMENT. AUG 2007. Need smaller version
!   for each particular program. On the other hand it is not yet clear
!   what final version of the program will look like...
!
!   Date:                  Programmer:                     History of changes:
!   ====                   ==========                      ==================
!   OCT 2007               B.Strokopytov                   Global debug included
!
!
TYPE cards
    INTEGER                                     :: nmol         ! expected number of molecules
    REAL(KIND=wp), DIMENSION(6)                 :: cell
    INTEGER                                     :: debug        ! debug level
    LOGICAL                                     :: xray         
    LOGICAL                                     :: anom
    LOGICAL                                     :: aniso
    REAL(KIND=wp)                               :: resol_min
    REAL(KIND=wp)                               :: resol_max
    REAL(KIND=wp)                               :: smin
    REAL(KIND=wp)                               :: smax
    REAL(KIND=wp)                               :: fo_low
    REAL(KIND=wp)                               :: fo_high
    REAL(KIND=wp)                               :: sigma_cutoff
!   Refinement:
    INTEGER                                     :: ngauss
    INTEGER                                     :: ncycles
    REAL(KIND=wp)                               :: bw
    LOGICAL                                     :: sigma_weighting
!   
!   Must be at most six chars, e.q. 'xyzbqu'
    CHARACTER(LEN=6)                            :: refinement_mode
!   Rotation:
    CHARACTER(LEN=4)                            :: angle_system
    CHARACTER(LEN=4)                            :: package
!   R-free:
    REAL(KIND=wp)                               :: percentage
!   PDB part:
    INTEGER                                     :: number_of_models
    INTEGER                                     :: number_of_static
    INTEGER                                     :: start_model
!   Phased translation function:
    LOGICAL                                     :: phased_translation_function
    INTEGER,          DIMENSION(:), ALLOCATABLE :: phase_sign
    INTEGER                                     :: mpeaks
    LOGICAL                                     :: interpolation
    REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE :: step
    REAL(KIND=wp)                               :: grid_factor
!   Packing:                                 
    LOGICAL                                     :: packing
    LOGICAL                                     :: real_shape 
    REAL(KIND=wp)                               :: overlap
    REAL(KIND=wp)                               :: atom_radius
    CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: atom_name
    INTEGER                                     :: maxsol
!   qualex:
    CHARACTER(LEN=132)                          :: weights_file_name
    CHARACTER(LEN=132)                          :: dimacs_file_name
    LOGICAL                                     :: greedy
    LOGICAL                                     :: qualex
!   QPB:
    LOGICAL                                     :: qpb
!   SF Table:
    LOGICAL                                     :: table
    REAL(KIND=wp)                               :: shannon_sampling
    REAL(KIND=wp)                               :: shannon_rate
!   Program labels:
    CHARACTER(LEN=30), DIMENSION(200)           :: lsprgi
    CHARACTER(LEN=1),  DIMENSION(200)           :: ctprgi
END TYPE cards
CONTAINS
    SUBROUTINE read_cards ( cards_1 )
!
!       Purpose:
!       =======
!       Reads input cards using module ISO_VARYING_STRING
!       FIXME:
!       Need to eliminate ISO_VARYING_STRING for the purpose of generality.
!       Currently this module is not possible to compile using gfortran.
!       There are also difficulties with xlf 8.1 compiler. Also I find it
!       somewhat difficult to do conversion to integers or real numbers
!       directly from VARYING_STRINGs.
!       
!       Date:           Programmer:              Description of changes:
!       ====            ==========               ======================
!       Oct 2005        B.Strokopytov            Original code
!
!       Nov 2005        B.Strokopytov            Functions .EQS., is_numeric, .NUMERIC. added
!                                                eliminating the need for temporary character
!                                                strings
!
!       Nov 2005        B.Strokopytov            Several bugs corrected: words(i) -> words(1)
!
!       Aug 2007        B.Strokopytov            Serious bug detected using Intel 10.0.025 compiler:
!                                                Unit 5 cannot REWIND properly
!                                                and therefore reading of cards
!                                                by CPP4 parser in MTZ_IO module
!                                                failed.
!
        TYPE(cards), INTENT(INOUT)                      :: cards_1
!       Local variables:
        CHARACTER(LEN=400)                              :: line
        CHARACTER(LEN=132)                              :: temp_file
        REAL(KIND=wp)                                   :: resol_temp
        CHARACTER(LEN=256)                              :: polarrfn_log_file
        INTEGER,              DIMENSION(2)              :: phase_sign
!       Parser variables:       
        CHARACTER(LEN=wlen),  DIMENSION(:), ALLOCATABLE :: words
        INTEGER                                         :: istat
!       Counters:
        INTEGER                                         :: i
        INTEGER                                         :: j
        INTEGER                                         :: pdb_model_number

!       Initialize default parameters:
        cards_1%nmol = 1
        cards_1%number_of_models = 0
        cards_1%number_of_static = 0
        cards_1%start_model = 1
        cards_1%cell         = (/0.0_wp, 0.0_wp, 0.0_wp, 90.0_wp, 90.0_wp, 90.0_wp/)
        cards_1%debug        = 0
        cards_1%anom         = .FALSE.
        cards_1%aniso        = .FALSE.
        cards_1%xray         = .TRUE.
        cards_1%resol_min    = 0.0_wp
        cards_1%resol_max    = 0.0_wp

!       Fobs limits:
        cards_1%fo_low       = 0.001_wp
        cards_1%fo_high      = 10.0_wp ** 7
        cards_1%sigma_cutoff = 0.0_wp

!       Refinement:
        cards_1%ngauss = 5 
        cards_1%refinement_mode = 'XYZBQ'
        cards_1%ncycles         = 3
        cards_1%bw              = 0.0_wp
        cards_1%sigma_weighting = .FALSE.

!       Phased translation:       
        cards_1%phased_translation_function = .FALSE.
        phase_sign = 0
        cards_1%interpolation = .TRUE.
        CALL allocate_array ( cards_1%step, 1 )
        cards_1%step = 9.0_wp
        cards_1%angle_system = 'EULE'
        cards_1%package = 'XPLO'

!       Peak search:
        cards_1%mpeaks  = 5
        cards_1%grid_factor = 4.0_wp
        cards_1%shannon_rate = cards_1%grid_factor / 2.0_wp
        cards_1%percentage  = 0.0_wp

!       Packing:
        cards_1%packing = .TRUE.
        cards_1%real_shape = .FALSE.
        cards_1%atom_radius = 3.0_wp
        cards_1%overlap = 0.06_wp
        CALL allocate_array ( cards_1%atom_name, 5 )
        cards_1%atom_name = (/ 'N   ', 'CA  ', 'CB  ', 'C   ', 'O   ' /)

!       Graph theory methods:
        cards_1%greedy = .FALSE.
        cards_1%qualex = .FALSE.
        cards_1%qpb = .FALSE.

!       Table:
        cards_1%table = .FALSE.
        cards_1%shannon_sampling = 4.0_wp
       
        CALL ccprcs( 6, 'ROTLSQ', '$Date: 2005/10/31 11:39:00 $' )

!       Estimate number of non-static PDB models on input:
        pdb_model_number = 0
        DO WHILE ( .TRUE. )
            pdb_model_number = pdb_model_number + 1
            CALL ugtenv ( 'PDB'//TRIM ( int_to_c ( pdb_model_number ) ), temp_file )
            IF ( LEN_TRIM ( temp_file ) == 0 ) EXIT
        ENDDO

        cards_1%number_of_models = pdb_model_number - 1

!       Check whether static model is present:
        CALL ugtenv ( 'STATIC', temp_file )
        IF ( LEN_TRIM ( temp_file ) /= 0 ) THEN
            cards_1%number_of_static = 1
        ENDIF    

!       Rewind the unit (just in case) other routines already read the cards:
        REWIND(UNIT=5, IOSTAT=istat)

!       Really big loop reading cards:
        DO WHILE ( .TRUE. )
!           Bug corrected AUG 2007 need to use UNIT=5, otherwise nothing to rewind later(?!):
            
            READ (UNIT=5, FMT="(A)", IOSTAT=istat) line

            IF ( istat == 0 ) THEN

                WRITE(*, "(' > ', A)" ) TRIM ( line )

!               Clean up everything after comment:
                j = INDEX ( line, '#' )
                IF ( j > 0 ) line(j:) = ' '

!               Move everyhting to upper case:
                CALL ucase   ( line )

!               Split line into separate words:
                CALL mysplit ( line, words )

!               Empty line detected:
                IF ( SIZE ( words ) == 0 ) CYCLE

            ELSE
                CALL messag ( 'End of cards.', 'read_cards' )
                EXIT
            ENDIF

            IF ( words(1) .EQS. 'CELL' ) THEN
                IF ( SIZE ( words ) == 7 ) THEN
                    DO j = 1, 6
                        cards_1%cell(j) = .NUMERIC. words(1+j)
                    ENDDO
                ELSE IF ( SIZE ( words ) == 4 ) THEN
                    DO j = 1, 3
                        cards_1%cell(j) = .NUMERIC. words(1+j)
                    ENDDO
                ELSE
                   CALL die ( 'Syntax error.', 'read_cards' )
                ENDIF
            ELSE IF ( words(1) .EQS. 'MODE' ) THEN

                IF ( SIZE ( words )  >= 2 )  THEN
                    cards_1%refinement_mode = TRIM ( words(2) )
                ELSE
                    CALL die ('Syntax error. No mode has been specified.', 'read_cards')
                ENDIF

            ELSE IF ( words(1) .EQS. 'NCYC' ) THEN
                IF ( SIZE ( words ) >= 2 ) THEN
                    READ(words(2),*) cards_1%ncycles
                ELSE
                    CALL die ('Syntax error. Number of cycles has not been specified.',&
                              'read_cards')
                ENDIF
            ELSE IF ( words(1) .EQS. 'NGAU' ) THEN
                IF ( SIZE ( words ) >= 2 ) THEN
                    READ(words(2),*) cards_1%ngauss
                    IF ( cards_1%ngauss /= 2 .AND. cards_1%ngauss /= 5 ) THEN
                        CALL die('Incorrect number of gaussians. Must use 2 or 5.','read_cards')
                    ENDIF
                ELSE
                    CALL die ('Syntax error. Number of gaussians has not been specified.',&
                              'read_cards')
                ENDIF
            ELSE IF ( words(1) .EQS. 'WEIG' ) THEN
                IF ( SIZE ( words ) > 1 )   THEN
                    IF ( is_numeric ( words(2) ) ) THEN
                        READ(words(2),*) cards_1%bw
                    ELSE
                        cards_1%sigma_weighting = .TRUE.
                    ENDIF
                ELSE
                    CALL die ('Syntax error. Weighting scheme not specified correctly.',&
                             'read_cards')
                ENDIF

!           Need to deal with it...                
            ELSE IF ( words(1) .EQS. 'NMOL' ) THEN

                IF ( SIZE ( words)  == 2 ) THEN
                    IF ( is_numeric ( words(2) ) ) THEN
                        cards_1%nmol = NINT ( .NUMERIC. words(2) )               
                    ELSE 
                        CALL ccperr ( 1, ' *** Syntax error: '//TRIM(line)//' ***' )
                    ENDIF

                ELSE IF ( SIZE ( words ) == 1 ) THEN
                    CALL ccperr( 1, ' *** Need to supply integer number of molecules ***' )
                ENDIF              
            ELSE IF ( words(1) .EQS. 'DEBU' ) THEN
                IF ( SIZE(words) == 1 ) THEN
                    CALL ccperr( 1, ' *** Need to supply integer value for DEBUG ***' )
                ELSE
                    IF ( is_numeric ( words(2) ) ) cards_1%debug = NINT ( .NUMERIC. words(2) )
                    debug = cards_1%debug
                ENDIF
            ELSE IF ( words(1) .EQS. 'FLIM' ) THEN

                IF ( SIZE ( words ) > 1 ) THEN
                    cards_1%fo_low = .NUMERIC. words(2)
                ENDIF

                IF ( SIZE ( words ) > 2 ) THEN
                    cards_1%fo_high = .NUMERIC. words(3)
                ENDIF

                IF ( SIZE ( words ) > 3 ) THEN
                    IF ( SIZE ( words ) == 4 ) THEN                        
                        CALL ccperr ( 1, ' *** Need to supply SIGMAcut keyword ***' )
                    ELSE
                        IF ( is_numeric ( words(5) ) ) THEN
                             cards_1%sigma_cutoff = .NUMERIC. words(5)
                        ELSE
                            CALL ccperr(1,' *** Need to supply real value for sigma cutoff ***')
                        ENDIF
                    ENDIF
                ENDIF

            ELSE IF ( words(1) .EQS. 'XRAY' ) THEN

                cards_1%xray = .TRUE.

                IF ( SIZE ( words ) > 1 ) THEN
                    IF ( words(2) .EQS. 'ANOM' ) THEN
                        cards_1%anom = .TRUE.
                    ELSE IF ( words(2) .EQS. 'NOAN' ) THEN
                        cards_1%anom = .FALSE.
                    ELSE
                        CALL ccperr(1, ' *** Syntax error: ANOM/NOANOM keyword should be supplied ***')
                    ENDIF
                ENDIF
! grid
                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'FACT' ) THEN
                        IF ( is_numeric ( words(i+1) ) ) THEN
                            cards_1%grid_factor = .NUMERIC. words(i+1)
                        ENDIF
                    ELSE IF ( words(i) .EQS. 'GRID' ) THEN
                        IF ( is_numeric ( words(i+1) ) ) THEN 
                            cards_1%grid_factor = .NUMERIC. words(i+1)
                        ENDIF
                    ENDIF
                ENDDO


            ELSE IF ( words(1) .EQS. 'NOXR' ) THEN
                cards_1%xray = .FALSE.
            ELSE IF ( words(1) .EQS. 'FILT' ) THEN
                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'GREE' ) THEN
                        cards_1%greedy = .TRUE.
                    ELSE IF ( words(i) .EQS. 'QUAL' ) THEN
                        cards_1%qualex = .TRUE.
                        CALL ugtenv ( 'DIMACS', cards_1%dimacs_file_name )
                        CALL ugtenv ( 'WEIGHTS', cards_1%weights_file_name )
                        IF ( LEN_TRIM ( cards_1%dimacs_file_name ) == 0 ) THEN
                            CALL die ('Please supply arbitrary name for DIMACS file...', 'read_cards' )
                        ELSE IF ( LEN_TRIM ( cards_1%weights_file_name ) == 0 ) THEN
                            CALL die ('Please supply arbiatrary name for WEIGHTS file...', 'read_cards' )
                        ENDIF
                    ELSE IF ( words(i) .EQS. 'QPB' ) THEN
                        cards_1%qpb = .TRUE.
                    ENDIF
                ENDDO
            ELSE IF ( words(1) .EQS. 'PACK' ) THEN
                cards_1%packing = .TRUE.
                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'OVER' ) THEN
                        IF ( i + 1 <= SIZE ( words ) ) THEN
                            IF ( is_numeric ( words(i+1) ) ) THEN
                                cards_1%overlap = .NUMERIC. words(i+1)
                            ELSE
                                CALL die ('Syntax error. Real value should follow OVERLAP keyword.', 'read_cards' )
                            ENDIF
                        ELSE
                            CALL die ('Syntax error. Real value should follow OVERLAP keyword.', 'read_cards' )
                        ENDIF
                    ENDIF
                
                    IF ( words(i) .EQS. 'RADI' ) THEN
                        IF ( i + 1 <= SIZE ( words ) ) THEN
                            IF ( is_numeric ( words(i+1) ) ) cards_1%atom_radius = .NUMERIC. words(i+1)
                            cards_1%atom_radius = MIN ( 3.0_wp, cards_1%atom_radius )
                            cards_1%atom_radius = MAX ( 2.2_wp, cards_1%atom_radius )
                        ENDIF    
                    ENDIF

                    IF ( ( words(i) .EQS. 'ATOM' ) .OR. ( words(i) .EQS. 'NAME' ) ) THEN
                        IF ( i + 1 <= SIZE ( words ) ) THEN
                            IF ( is_numeric ( words(i+1) ) ) j = NINT ( .NUMERIC. words(i+1) )
                            ioerr:IF ( is_numeric ( words(i+1) ) ) THEN
                                CALL allocate_array ( cards_1%atom_name, j )
                                IF (  i + 2 + j - 1  <= SIZE ( words ) ) THEN
                                    cards_1%atom_name = words(i+2:i+2+j-1)
                                ELSE
                                    WRITE(*,*) i + 2 + j - 1, SIZE ( words )
                                    CALL die ( 'Syntax error. Number of atoms should be '//TRIM(int_to_c(j)), 'read_cards' )
                                ENDIF
                            ELSE ioerr
                                IF ( words(i+1) .EQS. 'ALL' ) THEN
                                    CALL allocate_array ( cards_1%atom_name, 1 )
                                    cards_1%atom_name = 'ALL'
                                ELSE
                                    DO j = i + 1, SIZE ( words )
                                        IF ( words(j) .EQS. 'END' ) EXIT
                                    ENDDO
                                    j = j - 1
                                    CALL allocate_array ( cards_1%atom_name, j - (i+1) + 1 )
                                    cards_1%atom_name = words(i+1:j)
                                ENDIF
                            ENDIF ioerr
                        ENDIF
                    ENDIF
                    
                    IF ( ( words(i) .EQS. 'REAL' ) .OR. ( words(i) .EQS. 'SHAP' ) ) THEN
                        cards_1%real_shape = .TRUE.
                    ENDIF
 
                    IF ( ( words(i) .EQS. 'MAXS' ) .OR. ( words(i) .EQS. 'SOLU' ) ) THEN
                        IF ( i + 1 <= SIZE ( words ) ) THEN
                            IF ( is_numeric ( words(i+1) ) ) THEN
                                cards_1%maxsol = .NUMERIC. words(i+1) 
                            ELSE
                                WRITE ( *,* ) TRIM ( line )
                                CALL die ( 'Syntax error when reading maxsol.', 'read_cards')
                            ENDIF
                        ENDIF
                    ENDIF

                ENDDO
            ELSE IF ( words(1) .EQS. 'RESO' ) THEN
                IF ( is_numeric ( words(2) ) .AND. is_numeric ( words(3) ) ) THEN
                    cards_1%resol_min = .NUMERIC. words(2)
                    cards_1%resol_max = .NUMERIC. words(3)
                ELSE
                    WRITE ( *,* ) TRIM ( line )
                    CALL die ( 'Syntax error when reading resol_min/resol_max.', 'read_cards')
                ENDIF

!               Swap if necessary:
                IF ( cards_1%resol_max > cards_1%resol_min ) THEN
                    resol_temp        = cards_1%resol_min
                    cards_1%resol_min = cards_1%resol_max
                    cards_1%resol_max = resol_temp
                ENDIF
                cards_1%smin = 1.0_wp / cards_1%resol_min ** 2
                cards_1%smax = 1.0_wp / cards_1%resol_max ** 2

            ELSE IF ( words(1) .EQS. 'PHAS' ) THEN

                cards_1%phased_translation_function = .TRUE.

                DO i = 2, SIZE ( words )
                    IF ( words(i) == '+' ) THEN
                        phase_sign(1) = +1
                    ELSE IF ( words(i) == '-' ) THEN
                        phase_sign(2) = -1
                    ELSE IF ( words(i) == '+-' ) THEN
                        phase_sign = (/1,-1/)
                    ELSE IF ( words(i) == '-+' ) THEN
                        phase_sign = (/1,-1/)
                    ENDIF  
                ENDDO

                IF ( COUNT ( phase_sign /= 0 ) == 0 ) THEN
                    CALL allocate_array ( cards_1%phase_sign, 1 )
                    cards_1%phase_sign(1) = +1
                ELSE
                    CALL allocate_array ( cards_1%phase_sign, COUNT ( phase_sign /= 0 ) )
                    cards_1%phase_sign = PACK ( (/1, -1/), phase_sign /= 0 )
                ENDIF

                DO i = 2, SIZE ( words )
                    IF ( ( words(i) .EQS. 'PEAK' ) .OR. ( words(i) .EQS. 'NPEA' ) ) THEN
                        IF ( SIZE ( words ) < i + 1 ) THEN
                            WRITE ( *, "(A)" ) TRIM ( line )
                            CALL die ( 'Syntax error. Integer value for PEAKS must be supplied.', 'read_cards' )
                        ENDIF
                        IF ( is_numeric ( words(i+1) ) ) THEN
                            cards_1%mpeaks = .NUMERIC. words(i+1) 
                        ELSE
                            WRITE ( *, "(A)") TRIM ( line )
                            CALL die ( 'Syntax error. Integer value for PEAKS must be supplied.', 'read_cards' )
                        ENDIF
                    ENDIF    
                ENDDO

!               Check whether we want interpolation or not:
                DO i = 2, SIZE( words )
                    IF ( words(i) .EQS. 'INTE' ) THEN
                        cards_1%interpolation = .TRUE.
                    ENDIF
                ENDDO
                
                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'NOINTE' ) THEN
                        cards_1%interpolation = .FALSE.
                    ENDIF

                    IF ( SIZE(words) >= i + 1 ) THEN
                        IF ( ( words(i) .EQS. 'NO' ) .AND. ( words(i+1) .EQS. 'INTE' ) ) THEN
                            cards_1%interpolation = .FALSE.
                        ENDIF
                    ENDIF
                ENDDO

!               Check whether angle step has been supplied:
                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'STEP' ) THEN
                        IF ( SIZE ( words ) < i + cards_1%number_of_models ) THEN
                            WRITE ( *, "(A)" ) TRIM ( line )
                            CALL die ( 'Syntax error. Number of angle step entries must be '&
                                       //TRIM(int_to_c(cards_1%number_of_models)), 'read_cards' )
                        ENDIF

                        CALL allocate_array ( cards_1%step, cards_1%number_of_models )

                        DO j = 1, cards_1%number_of_models
                            IF ( is_numeric ( words(i+j) ) ) THEN
                                cards_1%step(j) = .NUMERIC. words(i+j)
                            ELSE
                                WRITE ( *, "( '> ', A )" ) TRIM ( line )
                                CALL die ( 'Syntax error. Number of angle step entries must be equal&
                                           & to the number of non-static pdb models on input', 'read_cards' )
                            ENDIF
                        ENDDO                        
                        EXIT
                    ENDIF
                ENDDO

                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'FACT' ) THEN
                        cards_1%grid_factor = .NUMERIC. words(i+1)
                    ELSE IF ( words(i) .EQS. 'GRID' ) THEN
                        cards_1%grid_factor = .NUMERIC. words(i+1)
                    ENDIF 
                ENDDO

                DO i = 2, SIZE ( words )
                    IF ( words(i) .EQS. 'LATT' ) THEN
                        cards_1%angle_system = 'LATT' 
                    ELSE IF ( words(i) .EQS. 'EULE' ) THEN
                        cards_1%angle_system = 'EULE'
                    ELSE IF ( words(i) .EQS. 'XPLO' ) THEN
                        cards_1%package = 'XPLO'
                    ELSE IF ( words(i) .EQS. 'CCP4' ) THEN
                        cards_1%package = 'CCP4'
                    ENDIF
                ENDDO

            ELSE IF ( words(1) .EQS. 'LABI' ) THEN
                cards_1%xray = .TRUE.
            ELSE IF ( words(1) .EQS. 'LABO' ) THEN
                cards_1%xray = .TRUE.
            ELSE IF ( words(1) .EQS. 'TABL' ) THEN
                cards_1%table = .TRUE.
                DO i = 1, SIZE ( words )
                    IF ( SIZE (words) >= i + 1 ) THEN
                        IF ( words(i+1) .EQS. 'SHAN' ) THEN
                            IF ( SIZE (words) >= i + 2 ) THEN
                                IF ( is_numeric ( words(i+2) ) ) THEN
                                    cards_1%shannon_sampling = .NUMERIC. words(i+2)
                                ELSE
                                    WRITE(*,*) TRIM ( line )
                                    CALL die ( 'Failed to read SHANON SAMPLING value...', 'read_cards' )
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ELSE IF ( words(1) .EQS. 'NOTA' ) THEN
                cards_1%table = .FALSE.
            ELSE IF ( words(1) .EQS. 'END' ) THEN
                CALL messag ( 'Finished reading cards.', 'read_cards' )
                EXIT
            ELSE
                IF ( ALLOCATED ( words ) ) THEN
                    WRITE( *, "(' Unrecognized keyword: ')" ) TRIM(words(1))
                ENDIF
            ENDIF
        ENDDO

        IF ( ALLOCATED ( words ) ) CALL deallocate_array ( words )

        cards_1%maxsol = MAX ( cards_1%maxsol, 200  )
        cards_1%maxsol = MIN ( cards_1%maxsol, 15000 )
        IF ( ALL ( cards_1%cell(1:3) > 0.0_wp ) ) THEN
            WRITE ( *, "(' READ_CARDS> ', 'CELL=', 6F10.3 )" ) cards_1%cell
        ENDIF
        WRITE(*,"(' READ_CARDS> ', 'RESOLUTION MIN=', F10.2, ' MAX=',F10.2)") cards_1%resol_min, cards_1%resol_max
        WRITE(*,"(' READ_CARDS> ', 'S**2       MIN=', F10.7, ' MAX=',F10.7)") cards_1%smin, cards_1%smax
        WRITE(*,"(' READ_CARDS> ', 'F_LOW=', ES9.2, ' F_HIGH= ', ES9.2, ' SIGMA CUT=', F5.1)") &
        cards_1%fo_low, cards_1%fo_high, cards_1%sigma_cutoff
        WRITE(*,"(' READ_CARDS> ', 'XRAY=', L1, ' ANOM=',L1, ' ANISO=',L1)") &
        cards_1%xray, cards_1%anom, cards_1%aniso
        WRITE(*,"(' READ_CARDS> ', 'NMOL=', I3)") cards_1%nmol

!       Start counting from zero if static model available:
        IF ( cards_1%number_of_static > 0 ) THEN
            cards_1%start_model = 0
        ENDIF
        WRITE(*,"(' READ_CARDS> ', 'NUMBER OF MODELS=', I3, ' NUMBER OF STATIC MODELS=', I2, ' START MODEL=',I2)")&
        cards_1%number_of_models, cards_1%number_of_static, cards_1%start_model

!       Check that we have some models:
!        IF ( cards_1%number_of_models == 0 ) THEN
!            CALL die ('Number of models to test equals zero...', 'read_cards' )
!        ENDIF

        IF ( cards_1%phased_translation_function ) THEN
            WRITE(*,"(' READ_CARDS> ', 'PHASED_TRANSLATION_FUNCTION ACTIVE= ',L1, ' PHASE SIGNS=', 2I3)")&
            cards_1%phased_translation_function, cards_1%phase_sign
        ELSE
            WRITE(*,"(' READ_CARDS> ', 'PHASED_TRANSLATION_FUNCTION ACTIVE= ',L1)") &
            cards_1%phased_translation_function
        ENDIF
        WRITE(*,"(' READ_CARDS> ', 'COMPLEX STRUCTURE FACTOR TABLE=', L1, ' SHANNON SAMPLING=', F5.1)") &
        cards_1%table, cards_1%shannon_sampling
        WRITE(*,"(' READ_CARDS> ', 'ANGLE SYSTEM=', A, ' FROM PACKAGE=', A, ' MPEAKS=', I3, ' ANGLE STEP=',200F5.1 )" ) &
        cards_1%angle_system, cards_1%package, cards_1%mpeaks, cards_1%step
        WRITE(*,"(' READ_CARDS> ', 'FFT GRID FACTOR=', F5.1, ' 3D INTERPOLATION=',L1)" ) &
        cards_1%grid_factor, cards_1%interpolation

        IF ( cards_1%packing ) THEN
            IF ( cards_1%atom_name(1) == 'ALL' ) THEN

!               If all atoms included in packing reduce atom radius:
                IF ( cards_1%atom_radius > 2.6_wp ) THEN
                    CALL messag ( 'Reducing atom radius since all atoms used for logical map calculations.',&
                                  'read_cards')
                    cards_1%atom_radius = 2.6_wp 
                ENDIF
            ENDIF

            IF ( cards_1%overlap > 0.12_wp ) THEN
                WRITE ( *, "(' READ_CARDS> ', 'max overlap=', F6.1, '%')") 100*cards_1%overlap
                CALL warn ( 'Trying to use very loose value for overlap', 'read_cards' )
                CALL messag ( 'Correcting overlap...', 'read_cards' )
                cards_1%overlap = 0.12_wp
            ELSE IF ( cards_1%overlap < 0.06_wp ) THEN
                WRITE ( *, "(' REAL_CARDS> ', 'max_overlap=', F6.1, '%')") 100*cards_1%overlap
                CALL warn ( 'Trying to use very tight value for overlap', 'read_cards' )
                CALL messag ( 'Correcting overlap...', 'read_cards' )
                cards_1%overlap = 0.06_wp
            ENDIF


            WRITE ( *, "(' READ_CARDS> ', 'MAX OVERLAP=', F6.3, ' ATOM RADIUS=', F5.1, ' REAL SHAPE=', L1)") &
            cards_1%overlap, cards_1%atom_radius, cards_1%real_shape
            WRITE ( *, "(' READ_CARDS> ', 'NUMBER OF SOLUTION TO ANALYZE=', I7, ' ATOM NAMES=', 100(1X,A))") &
            cards_1%maxsol, cards_1%atom_name
            
        ENDIF

        IF ( .NOT. cards_1%greedy .AND. .NOT. cards_1%qualex ) THEN
            cards_1%greedy = .TRUE.
        ENDIF

        IF ( cards_1%greedy ) THEN
            CALL messag ( 'Going to use greedy searches...', 'read_cards')
        ENDIF

        IF ( cards_1%bw == 0.0_wp .AND. .NOT. cards_1%sigma_weighting ) THEN
            CALL messag( 'Going to use unit weighting...', 'read_cards' )
        ELSE IF ( cards_1%bw /= 0.0_wp .AND. .NOT. cards_1%sigma_weighting ) THEN
            CALL messag( 'Going to use exponential B weighting...', 'read_cards' )
        ELSE
            CALL messag( 'Going to use SIGMA based weighting scheme...', 'read_cards')
        ENDIF
        IF ( cards_1%qualex ) THEN
            CALL messag ( 'Going to use qualex_ms algorithm...', 'read_cards' )
            WRITE ( *, "(' READ_CARDS> ', 'DIMACS FILE=', A )" )  TRIM ( cards_1%dimacs_file_name )
            WRITE ( *, "(' READ_CARDS> ', 'WEIGHTS FILE=', A )" ) TRIM ( cards_1%weights_file_name )
        ENDIF

        IF ( cards_1%qpb ) THEN
            CALL messag ( 'Going to use QPB algorithm...', 'read_cards' )
        ENDIF
        

!       Rewind default input unit for further reading by CCP4 parsers or others:
        REWIND ( UNIT=5, IOSTAT=istat )
        IF ( istat /= 0 ) THEN
            CALL die ( 'Failed to rewind UNIT 5', 'read_cards' )
        ENDIF

    END SUBROUTINE read_cards
 
    SUBROUTINE change_2_euler ( my_cards )
!
!       Purpose:
!       =======
!       Switches definition of angle system to eulerian.
!
!       Date:             Programmer:           History of changes:
!       ====              ==========            ==================
!       Aug 2005          B.Strokopytov         Original code
!
        TYPE(cards), INTENT(INOUT) :: my_cards
        my_cards%angle_system = 'EULE'
    END SUBROUTINE

    SUBROUTINE deallocate_cards ( my_cards )
!
!       Purpose:
!       =======
!       Deallocates CARDS type arrays. Should be used near the main program
!       completion.
!
!       Date:              Programmer:            History of changes:
!       ====               ==========             ==================
!       Aug 2005           B.Strokopytov          Original code
!       Aug 2007           B.Strokopytov          ISTAT check added
! 
        TYPE(cards), INTENT(INOUT) :: my_cards
!       Local variables:
        INTEGER                    :: istat

        IF ( ALLOCATED ( my_cards%phase_sign ) ) THEN
            DEALLOCATE ( my_cards%phase_sign, STAT=istat )
            IF ( istat /= 0 ) THEN
                WRITE(*,*) ' istat=', istat
                CALL die ('Failed to deallocate PHASE SIGN array.', 'deallocate_cards' )
            ENDIF
        ENDIF

        IF ( ALLOCATED ( my_cards%step ) ) THEN
            DEALLOCATE ( my_cards%step, STAT=istat )
            IF ( istat /= 0 ) THEN
                WRITE(*,*) ' istat=', istat
                CALL die ('Failed to deallocate STEP array.', 'deallocate_cards' )
            ENDIF
        ENDIF

        IF ( ALLOCATED ( my_cards%atom_name ) ) THEN
            DEALLOCATE ( my_cards%atom_name, STAT=istat )
            IF ( istat /= 0 ) THEN
                WRITE(*,*) ' istat=', istat
                CALL die ('Failed to deallocate ATOM_NAME array.', 'deallocate_cards' )
            ENDIF
        ENDIF

        CALL messag ('Cards have been successfully deallocated...', 'deallocate_cards' )

    END SUBROUTINE deallocate_cards

END MODULE ucards
