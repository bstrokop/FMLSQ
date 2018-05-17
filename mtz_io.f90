MODULE mtz_io
USE constants
USE cell_manip
USE sorting_facilities
USE fail
USE fft_util
USE ucards
USE util
IMPLICIT NONE
!
!
!   Date:              Programmer:         Description of changes:
!   ====               ==========          ======================
!   
!   Oct 2005           Strokopytov B.      Complete rewrite
!   Nov 2005           Strokopytov B.      Bug (fc wp real) corrected
!   Aug 2007           Strokopytov B.      Problem when reading cards
!                                          through CCP4 PARSER detected
!                                          It is better consistently use
!                                          unit=5 throughout the program
!                                          Actually we need to read input
!                                          twice. Once in UCARDS module,
!                                          and once in MTZ_IO module to get
!                                          labels.
!                                          This needs further attention
!                                          FIXME
! 
!
!   Aug 2007           Strokopytov B.      Cosmetic changes
!   Aug 2007                -"-            LABOUT still needs attention
!                                          (has not been implemeted yet)
!
TYPE :: mtz
     INTEGER,          DIMENSION(:),    ALLOCATABLE :: lookup 
     INTEGER                                        :: number_of_columns_on_input
     REAL(KIND=wp)                                  :: grid_factor
     REAL(KIND=wp)                                  :: resol_min
     REAL(KIND=wp)                                  :: resol_max
     REAL(KIND=wp)                                  :: s_min_squared
     REAL(KIND=wp)                                  :: s_max_squared
     TYPE(space_group)                              :: sp_group
     REAL(KIND=wp),    DIMENSION(6)                 :: cell
     REAL(KIND=wp)                                  :: volume
     TYPE(matrix)                                   :: ORT
     TYPE(matrix)                                   :: DEORT
     INTEGER                                        :: hmin
     INTEGER                                        :: hmax
     INTEGER                                        :: kmin
     INTEGER                                        :: kmax
     INTEGER                                        :: lmin
     INTEGER                                        :: lmax
     TYPE(tensor)                                   :: real_tensor
     TYPE(tensor)                                   :: recip_tensor
     INTEGER                                        :: estimated_number_of_reflections
     INTEGER                                        :: current_number_of_reflections
!    Default percentage of "TEST" reflections:
     REAL(KIND=wp)                                  :: percentage = 5.0_wp
     REAL(KIND=wp)                                  :: fft_radius
     REAL(KIND=wp)                                  :: badd                ! Navaza's anitalias
     REAL(KIND=wp)                                  :: b_scale             ! Current anitalias
!    Refinement:
     REAL(KIND=wp)                                  :: fofc_scale
     REAL(KIND=wp)                                  :: fofc_biso
     REAL(KIND=wp)                                  :: bulk_scale
     REAL(KIND=wp)                                  :: bulk_biso   
!    Size of ndatasets:
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: pname
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: xname
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: dname
!    Size of nlprgo:
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: pname_out
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: xname_out
     CHARACTER(LEN=64), DIMENSION(:),   ALLOCATABLE :: dname_out
!    Size of NDATASETS
     REAL(KIND=sp),     DIMENSION(:,:), ALLOCATABLE :: datcell
     REAL(KIND=sp),     DIMENSION(:),   ALLOCATABLE :: datwave
     INTEGER                                        :: ndatasets
     TYPE(vector_int),  DIMENSION(:),   ALLOCATABLE :: hkl
     TYPE(vector_int),  DIMENSION(:),   ALLOCATABLE :: hkl_in_P1
!    Indices for refinement in orthogonal system:
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: HO
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: KO
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: LO
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: s
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: s_in_P1
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: eps
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: fo
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: sigfo
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: fomo
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: mult
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: phio
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: hla
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: hlb
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: hlc
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: hld
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: a
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: b
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fc
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: phic
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fc_in_P1           ! in P1 hemisphere
     COMPLEX(KIND=wp),  DIMENSION(:,:), ALLOCATABLE :: fcq_in_P1          ! in P1 hemisphere
     COMPLEX(KIND=wp),  DIMENSION(:,:), ALLOCATABLE :: fcq_sym            ! in P1 hemisphere
     COMPLEX(KIND=wp),  DIMENSION(:,:), ALLOCATABLE :: fc_temp            ! in P1 hemisphere
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: wphases            ! in P1 hemisphere
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fobs               ! phased fo in P1 hemisphere
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: fo_in_P1           ! - " -
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: weight             ! - " -
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: weight_in_P1       ! - " -
     INTEGER,           DIMENSION(:),   ALLOCATABLE :: Rfree_flag
     INTEGER,           DIMENSION(:),   ALLOCATABLE :: Rfree_flag_in_P1
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: sigfo_in_P1
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fc_in_P1_init
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fpart
     COMPLEX(KIND=wp),  DIMENSION(:),   ALLOCATABLE :: fpart_in_P1
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: sc
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: sc_in_P1
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: rres
     INTEGER,           DIMENSION(:),   ALLOCATABLE :: bin
     TYPE(vector_int),  DIMENSION(:),   ALLOCATABLE :: ideal_hkl
     INTEGER,           DIMENSION(:),   ALLOCATABLE :: ind
!     INTEGER,           DIMENSION(:),   ALLOCATABLE :: tag
     INTEGER,           DIMENSION(:),   ALLOCATABLE :: isort
!    Bulk solvent:
     REAL(KIND=wp),     DIMENSION(:),   ALLOCATABLE :: pbulk
     REAL(KIND=wp)                                  :: radius     
!    Need this for output:(need to deallocate) FIXME:
     CHARACTER(LEN=1),  DIMENSION(:),   ALLOCATABLE :: ctprgo
     CHARACTER(LEN=30), DIMENSION(:),   ALLOCATABLE :: lsprgo
!     INTEGER,           DIMENSION(:),   ALLOCATABLE :: csetout
END TYPE

INTERFACE allocated
    MODULE PROCEDURE allocated_mtz
END INTERFACE allocated

! Set Global CCP4 Parameters for mtz i/o operations..
INTEGER, PARAMETER :: MAXPAR  = 200
INTEGER, PARAMETER :: MCOLS   = 200
INTEGER, PARAMETER :: MAXSETS =  20
INTEGER, PARAMETER :: MAXSYM  = 192 
INTEGER, PARAMETER :: MCOLS_OUT = 16

CONTAINS
     
!   Bits cut out of wilson.f
!   5/3/2004: Updated to show reading of crystals and datasets

!   This program reads a basic MTZ file. Nothing is done
!   with the data read in.

!   It should be run as:
!
!   mtz_read HKLIN foo.mtz << eof
!   LABI F=filelabel SIGF=filelabel
!   END
!   eof

    SUBROUTINE read_mtz ( my_mtz, cards_1 )
        TYPE(mtz),                    INTENT(INOUT) :: my_mtz
!       Need to change grid_factor sometimes:
        TYPE(cards),                  INTENT(INOUT) :: cards_1
!       Local variables:
        INTEGER                                     :: lentit
        INTEGER                                     :: mtzerr
        INTEGER                                     :: mtzin
        INTEGER                                     :: mtzprt
        INTEGER                                     :: nlprgi
        INTEGER                                     :: numcol
        INTEGER                                     :: nreflx
        INTEGER                                     :: nspgrp
!        INTEGER                                     :: nsym
        INTEGER                                     :: nsymp
        INTEGER                                     :: ndatasets
        INTEGER, DIMENSION(3)                       :: inhkl
        INTEGER, DIMENSION(MCOLS)                   :: mtzlok
        INTEGER, DIMENSION(MAXSETS)                 :: isets
        INTEGER                                     :: ncol
        INTEGER, DIMENSION(MCOLS)                   :: csetid
        INTEGER, DIMENSION(MCOLS)                   :: csetout
        INTEGER                                     :: setid
        REAL(KIND=sp)                               :: s
        TYPE(vector_int)                            :: hkl
        REAL(KIND=wp)                               :: fo
        REAL(KIND=wp)                               :: sfo
        REAL(KIND=wp)                               :: phio
        REAL(KIND=wp)                               :: fomo
        REAL(KIND=wp)                               :: eps
        REAL(KIND=wp)                               :: phicen
        REAL(KIND=wp)                               :: mult
!       Changed to allocatable AUG 2007 BVS:
!        INTEGER, DIMENSION(:), ALLOCATABLE          :: seed
!        INTEGER                                     :: seed_size
!        REAL(KIND=wp)                               :: xseed
        LOGICAL                                     :: centre
        LOGICAL                                     :: sysabs

!        REAL(KIND=sp),     DIMENSION(4,4,MAXSYM)    :: rsym
        REAL(KIND=sp),     DIMENSION(6)             :: celmtz
        REAL(KIND=sp),     DIMENSION(2,MCOLS)       :: rngmtz
        REAL(KIND=sp),     DIMENSION(MCOLS)         :: recin
        REAL(KIND=sp)                               :: sminwl
        REAL(KIND=sp)                               :: smaxwl
        REAL(KIND=sp),     DIMENSION(6,MAXSETS)     :: datcell
        REAL(KIND=sp),     DIMENSION(MAXSETS)       :: datwave
        LOGICAL                                     :: mtzeof
        LOGICAL,           DIMENSION(MCOLS)         :: logmss
        CHARACTER(LEN=1)                            :: lattyp
        CHARACTER(LEN=10)                           :: namspg
        CHARACTER(LEN=10)                           :: pgname
        CHARACTER(LEN=10)                           :: versnx
        CHARACTER(LEN=1),  DIMENSION(MCOLS)         :: ctprgi
        CHARACTER(LEN=30), DIMENSION(MCOLS)         :: lsprgi
        CHARACTER(LEN=1),  DIMENSION(MCOLS)         :: ctyps
        CHARACTER(LEN=30), DIMENSION(MCOLS)         :: clabs
        CHARACTER(LEN=80)                           :: titin
        CHARACTER(LEN=64), DIMENSION(MAXSETS)       :: pname
        CHARACTER(LEN=64), DIMENSION(MAXSETS)       :: xname
        CHARACTER(LEN=64), DIMENSION(MAXSETS)       :: dname

!       Parser variables:
        REAL(KIND=sp),     DIMENSION(MAXPAR)        :: fvalue
        INTEGER,           DIMENSION(MAXPAR)        :: ibeg
        INTEGER,           DIMENSION(MAXPAR)        :: idec
        INTEGER,           DIMENSION(MAXPAR)        :: iend
        INTEGER,           DIMENSION(MAXPAR)        :: ityp
        INTEGER                                     :: ntok
        LOGICAL                                     :: lend
        CHARACTER(LEN=4)                            :: key
        CHARACTER(LEN=400)                          :: line
        CHARACTER(LEN=4),  DIMENSION(MAXPAR)        :: cvalue
!
!       Counters:
        INTEGER                                     :: i
        INTEGER                                     :: j
        INTEGER                                     :: run 
!       Statistics:
        INTEGER                                     :: num_refl_in_file
        INTEGER                                     :: num_hi_resol_reflns
        INTEGER                                     :: num_sigma_cut_off
        INTEGER                                     :: num_not_in_asu
        INTEGER                                     :: num_centric_reflns
        INTEGER                                     :: num_centric_err_reflns
        INTEGER                                     :: num_unobs
        INTEGER                                     :: num_sysabs
        INTEGER                                     :: num_passed
        INTEGER                                     :: num_eps
        INTEGER                                     :: num_fomo_zero
        INTEGER                                     :: num_fo_low
        INTEGER                                     :: num_fo_high
        INTEGER                                     :: num_fomo_very_low
        INTEGER                                     :: num_isNAN
        CHARACTER(LEN=120)                          :: file_1
!       Integration sphere fit check:
        REAL(KIND=wp), DIMENSION(3)                 :: ort_diag
        INTEGER                                     :: Rfree_flag
        INTEGER                                     :: nlprgo
        INTEGER                                     :: iset
!       Printing:
        INTEGER                                     :: nmod
        DATA (lsprgi(j),j=1,200) /'H','K','L','FP','SIGFP','PHIO','FOMO', 193*' '/
        DATA (ctprgi(j),j=1,200) /'H','H','H','F', 'Q',    'P',   'W',    193*' '/
        DATA  mtzlok             /7*-1,                                   193*0/

!       Sort order mainly for output HKL file:
        CALL allocate_array(my_mtz%isort,5)

!       This seems to be different in REFMAC (see subvag.f, rcard_tor1.f):
        my_mtz%isort=(/1,2,3,0,0/)

!       To start with set number of active program labels:
        nlprgi = 5
        IF ( cards_1%phased_translation_function ) THEN

!           Need 7 labels to run phased translation function:
            nlprgi = 7

        ELSE

!           No need to read phased mtz -> clear program labels:
            DO i = 6, 7
                lsprgi(i) = ' '
                ctprgi(i) = ' '
                mtzlok(i) = 0
            ENDDO
        ENDIF

!       CCP4 initialisations:
        CALL ccpfyp

!       MTZ-specific initialisations:
        CALL mtzini

!       Open input file:
        mtzin  = 1
        mtzprt = 0
        mtzerr = 0
        CALL lropen ( mtzin, 'HKLIN', mtzprt, mtzerr )

!       Parse keyworded input:
 10     CONTINUE
        line = '   '
        key  = ' '
        ntok = maxpar
        lend = .FALSE.
        
!
!       ***********************************************************
        CALL parser ( key, line, ibeg, iend, ityp, fvalue, cvalue, idec, ntok, lend, .TRUE. )
!       ***********************************************************

!       End of file?:
        IF ( lend ) THEN
            CALL warn ( 'EOF detected before reading END card...', 'read_mtz' )
            GO TO 50
        ENDIF

        IF ( key == 'LABI' ) THEN

!           Parse LABIN line:
            CALL lkyin ( mtzin, lsprgi, nlprgi, ntok, line, ibeg, iend )

        ELSE IF ( key(1:3) == 'END') THEN

            CALL messag ( 'End of cards detected...', 'read_mtz' )
            GO TO 50

        ELSE
!            WRITE (6,FMT=6022) LINE(1:LENSTR(LINE))
! 6022       FORMAT (' Error: Key_Word line NOT Understood:-',/' ',A)
            GO TO 10
        ENDIF

        GO TO 10

!       End of keyworded input:
     50 CONTINUE
        
!       Get header info:
        numcol = mcols
        CALL lrinfo ( mtzin, versnx, numcol, nreflx, rngmtz )
        CALL lrtitl ( mtzin, titin, lentit )
        CALL lrcell ( mtzin, celmtz )

!       This is MTZ cell:
        my_mtz%cell = REAL ( celmtz, KIND=wp )
        WRITE(*, "(' READ_MTZ> ', 'MTZ cell=', 6F10.4)") my_mtz%cell


        CALL lrsymi ( mtzin, nsymp, lattyp, nspgrp, namspg, pgname )
        IF ( debug >= 10 ) THEN
            WRITE(*,*) ' namspg=', '"', TRIM(namspg), '"'
        ENDIF

!       Taking just the most unique identifier:
        my_mtz%sp_group%space_group_number = nspgrp

!       Prepare all symmetry operators:
        CALL ugtenv ( 'SYMOP', file_1 )
        CALL find_symop_in_library ( my_mtz%sp_group, file_1 )
        CALL print_space_group_contents ( my_mtz%sp_group )


        IF ( debug >= 20 ) THEN
            WRITE(*,*) ' lsprgi->', lsprgi(1:10)(1:10), SIZE ( lsprgi )
            WRITE(*,*) ' ctprgi->', ctprgi(1:10), SIZE ( ctprgi )
            WRITE(*,*) ' mtzlok->', mtzlok(1:10), SIZE ( mtzlok )
            WRITE(*,*) ' nlprgi->', nlprgi
        ENDIF

!       Set up column assignments based on LABIN:
        CALL lrassn ( mtzin, lsprgi, nlprgi, mtzlok, ctprgi )
        CALL lrclab(mtzin,clabs,ctyps,ncol)
        ndatasets = MAXSETS
!       Check integer location array:
        IF ( debug >= 20 ) THEN        
            WRITE(*,*) ' mtzlok-> ', mtzlok(1:10)
        ENDIF

        IF ( .NOT. cards_1%phased_translation_function ) THEN

! FIXME
!           This should be 6 in case of Rfree:
            IF ( ANY ( mtzlok(1:5) <= 0 ) ) THEN
                CALL die ( 'Failed to get all column pointers', 'read_mtz' )
            ENDIF

! FIXME
!           This should be 6 in case of Rfree:
            CALL allocate_array ( my_mtz%lookup, 5 )
            my_mtz%lookup = mtzlok(1:5)

        ELSE

!           This may change to 8 in case of Rfree(FIXME):
            IF ( ANY ( mtzlok(1:7) <= 0 ) ) THEN
                CALL messag ( 'Check your input labels', 'read_mtz' )
                j = MAXVAL ( LEN_TRIM( lsprgi(1:7) ) )
                WRITE(*, "(' READ_MTZ> ',1X,A,1X,I3)") ( lsprgi(i)(1:j), mtzlok(i), i=1,7 )
                CALL die ( 'Failed to get all column pointers', 'read_mtz' )
            ENDIF
            CALL allocate_array ( my_mtz%lookup, 7 )
            my_mtz%lookup = mtzlok(1:7)

        ENDIF

!       *** For CCP4 5.0 libraries
!       Read in the crystal, project, dataset and set id
!       information, plus the cell parameters and wavelengths,
!       associated with each dataset

        ndatasets = maxsets
        CALL lridx ( mtzin, pname, xname, dname, isets, &
                     datcell, datwave, ndatasets )
        my_mtz%ndatasets = ndatasets

!       MCOLS_OUT columns on output would be sufficient:
        nlprgo = MCOLS_OUT
        CALL allocate_array(my_mtz%lsprgo,  MCOLS_OUT)
        CALL allocate_array(my_mtz%ctprgo,  MCOLS_OUT)
!        CALL allocate_array(my_mtz%csetout, MCOLS_OUT)

!       Copy FP SIGFP information:
        my_mtz%lsprgo(4:5) = clabs(my_mtz%lookup(4:5))
        my_mtz%ctprgo(4:5) = ctyps(my_mtz%lookup(4:5))

!       Allocate various arrays:
        CALL allocate_array(my_mtz%pname, ndatasets)
        CALL allocate_array(my_mtz%xname, ndatasets)
        CALL allocate_array(my_mtz%dname, ndatasets)
        CALL allocate_array(my_mtz%datwave, ndatasets)
        CALL allocate_array(my_mtz%datcell, 6, ndatasets)

        IF ( ndatasets > 0 ) THEN
            CALL lrclid(mtzin, csetid, ncol)
            my_mtz%datcell(1:6,1:ndatasets) = datcell(1:6,1:ndatasets)
            my_mtz%datwave(1:ndatasets) = datwave(1:ndatasets)
        ENDIF

!       Need these labels for output:
        CALL allocate_array(my_mtz%pname_out, nlprgo)
        CALL allocate_array(my_mtz%xname_out, nlprgo)
        CALL allocate_array(my_mtz%dname_out, nlprgo)

        my_mtz%cell = 0.0_wp

        IF ( my_mtz%ndatasets > 0 ) THEN

            my_mtz%pname = pname(1:ndatasets)
            my_mtz%xname = xname(1:ndatasets)
            my_mtz%dname = dname(1:ndatasets)

            csetout(1:5) = csetid(my_mtz%lookup(1:5))
!           All new columns inherit dataset of FP:
            csetout(6:nlprgo) = csetid(my_mtz%lookup(4))

            IF ( debug > 10000 ) THEN
                WRITE(*,*) ' CSETID=', CSETID(my_mtz%lookup(1:5))
                WRITE(*,*) ' INHERIT CSETID=', CSETID(my_mtz%lookup(4))
                WRITE(*,*) ' isets=',          isets(1:5)
            ENDIF
          
            DO i = 1, nlprgo
                setid = csetout(i)
                DO iset=1, my_mtz%ndatasets
                    IF ( isets(iset) == setid ) THEN
                        my_mtz%pname_out(i) = pname(iset)
                        my_mtz%xname_out(i) = xname(iset)
                        my_mtz%dname_out(i) = dname(iset)
                        my_mtz%cell(1:6) = datcell(1:6,iset)
                    ENDIF
                ENDDO
            ENDDO

        ENDIF

!       If NDATASETiS equals zero then we need to call LRCELL:        
        IF ( my_mtz%cell(1) == 0.0_wp ) THEN
            CALL lrcell(mtzin, my_mtz%cell)
        ENDIF

!       Check N datasets on input:
        IF ( ndatasets > 0 ) THEN 
            DO i = 1, ndatasets 
                IF ( pname(i) /= ' ') THEN 
                    WRITE(*,"(' READ_MTZ> ', 'project=',A, ' crystal=',A,' dname=',A,    &
             &                ' datwave=',F10.5, ' set=',I3)")                           & 
                                TRIM ( pname(i) ), TRIM ( xname(i) ), TRIM ( dname(i) ), &
                                datwave(i), i
                 ENDIF
            ENDDO
        ENDIF
        CALL lrrsol ( mtzin, sminwl, smaxwl )
        WRITE(*,"(' READ_MTZ_limits> ', 'smin(whole)= ',F8.5,' smax(whole)=',F8.5)")&
        sminwl, smaxwl

!       resol max not higher than 0.1 Ang resolution allowed:
        my_mtz%resol_max = 1.0_wp / SQRT ( MIN( 100.0_wp, smaxwl ) )

!       calculate hmin, hmax limits:
        IF ( cards_1%resol_min == 0.0_wp ) THEN


!           resol min not lower than 10000.0 Ang resolution allowed:
            my_mtz%resol_min = 1.0_wp / SQRT ( MAX( 0.1_wp ** 8, sminwl ) )

        ELSE

            my_mtz%resol_min = cards_1%resol_min

!           Resol max on cards will be ignored if it exceeds resolution in file:
            my_mtz%resol_max = MAX ( my_mtz%resol_max, cards_1%resol_max )

        ENDIF

!       Calculate hmin, hmax, etc. limits:
        my_mtz%hmax = INT ( my_mtz%cell(1) / my_mtz%resol_max )
        my_mtz%kmax = INT ( my_mtz%cell(2) / my_mtz%resol_max )
        my_mtz%lmax = INT ( my_mtz%cell(3) / my_mtz%resol_max )

        my_mtz%hmin = 0
        my_mtz%kmin = -my_mtz%kmax
        my_mtz%lmin = -my_mtz%lmax

        my_mtz%s_min_squared = 1.0_wp / my_mtz%resol_min ** 2
        my_mtz%s_max_squared = 1.0_wp / my_mtz%resol_max ** 2


        CALL calcmetric ( my_mtz%cell, my_mtz%ORT, my_mtz%DEORT, &
                          my_mtz%REAL_TENSOR, my_mtz%RECIP_TENSOR, my_mtz%volume )

        CALL check_sym_cell_consistency(my_mtz%REAL_TENSOR, my_mtz%sp_group)



!        Alternatively:
!        *** For CCP4 4.2 (and earlier)
!        Read in the project, dataset and set id
!        information, plus the cell parameters and wavelengths,
!        associated with each dataset
!        NDATASETS = MaxSets
!        CALL LRIDC(MTZIN,PNAME,DNAME,ISETS,
!     +       DATCELL,DATWAVE,NDATASETS)

!       Initialise mtz data:
        double_run:DO run = 0, 1

!           Double run is needed to estimate number of reflections in mtz file:
            run_1:IF ( run == 1 ) THEN
                nmod = MAX ( my_mtz%estimated_number_of_reflections / 20, 1 )

!               Allocate_mtz arrays:
                CALL allocate_array ( my_mtz%hkl, my_mtz%estimated_number_of_reflections, 'hkl')
                CALL allocate_array ( my_mtz%s,   my_mtz%estimated_number_of_reflections, 's' )
                CALL allocate_array ( my_mtz%eps, my_mtz%estimated_number_of_reflections, 'eps' )  
                CALL allocate_array ( my_mtz%mult, my_mtz%estimated_number_of_reflections, 'mult' )  
                CALL allocate_array ( my_mtz%fo,  my_mtz%estimated_number_of_reflections, 'fo' )

                IF ( cards_1%percentage > 0.0_wp ) THEN
                    CALL allocate_array ( my_mtz%Rfree_flag, my_mtz%estimated_number_of_reflections, &
                                         'Rfree_flag')
                ENDIF

                IF ( my_mtz%lookup(5) > 0 ) THEN
                    CALL allocate_array ( my_mtz%sigfo, my_mtz%estimated_number_of_reflections, 'sigfo' )
                ENDIF

                IF ( cards_1%phased_translation_function ) THEN
                    CALL allocate_array ( my_mtz%phio, my_mtz%estimated_number_of_reflections, 'phio' )
                    CALL allocate_array ( my_mtz%fomo, my_mtz%estimated_number_of_reflections, 'fomo' )
                ENDIF

!               This does not take into account Rfree flag (FIXME):
                IF ( SIZE ( my_mtz%lookup ) == 11 ) THEN
                    IF ( cards_1%phased_translation_function .AND. my_mtz%lookup(8) > 0 ) THEN
                        CALL warn ('Pointer to FOMO column does not exist.', 'read_mtz')
                        CALL allocate_array ( my_mtz%hla, my_mtz%estimated_number_of_reflections, 'hla' )
                        CALL allocate_array ( my_mtz%hlb, my_mtz%estimated_number_of_reflections, 'hlb' )
                        CALL allocate_array ( my_mtz%hlc, my_mtz%estimated_number_of_reflections, 'hlc' )
                        CALL allocate_array ( my_mtz%hld, my_mtz%estimated_number_of_reflections, 'hld' )
                    ELSE
                        CALL messag('Hendrickson-Lattmann coefs. are absent on input at the moment.', 'read_mtz')
                    ENDIF
                ENDIF

                IF ( cards_1%phased_translation_function ) THEN
                    CALL messag('    H   K   L      FO      SIG(FO)     PHIO      FOMO    MULT', 'read_mtz' )
                ELSE
                    CALL messag('    H   K   L      FO      SIG(FO)     EPS       MULT', 'read_mtz' )
                ENDIF
            ENDIF run_1

!           Initialize counters:
            num_hi_resol_reflns    = 0
            num_unobs              = 0
            num_sigma_cut_off      = 0
            num_not_in_ASU         = 0
            num_sysabs             = 0
            num_centric_reflns     = 0
            num_centric_err_reflns = 0
            num_eps                = 0
            num_fomo_zero          = 0
            num_fo_low             = 0
            num_fo_high            = 0
            num_fomo_very_low      = 0
            num_passed             = 0
            num_refl_in_file       = 0
            num_isNaN              = 0

!           Loop over reflections.
!           lrrefl returns columns in file order:
            read_hkl:DO WHILE ( .TRUE. )       
                CALL lrrefl ( mtzin, s, recin, mtzeof )
                IF ( mtzeof ) EXIT
                num_refl_in_file = num_refl_in_file + 1

!               Check for Missing Number Flags:
                CALL lrrefm ( mtzin, logmss )

!               IF( logmss(mtzlok(4)) .OR. logmss(mtzlok(5))) GO TO 110
                DO i = 1, SIZE ( my_mtz%lookup )
                    IF ( logmss(i) ) recin(i) = 0.0
                ENDDO

!               s is now slightly more accurate:
                s = my_mtz%RECIP_TENSOR * NINT ( recin(1:3) )

!               Check resolution limits:
!               IF ( s < sminwl .OR. s > smaxwl ) GO TO 110
        
                IF ( s < my_mtz%s_min_squared .OR. s > my_mtz%s_max_squared ) THEN
                    num_hi_resol_reflns  = num_hi_resol_reflns + 1
                    CYCLE 
                ENDIF

                IF ( 1.0_wp / SQRT ( s ) < my_mtz%resol_max ) THEN
                    WRITE ( *, "('passed',3I5,A,F12.9)" ) NINT ( recin(1:3) ),' d=', 1.0_wp/SQRT(s)
                    CALL warn ( 'Programming error. Incorrect calculation of s.', 'read_mtz' )
                ENDIF

!               Read data:
                inhkl = NINT ( recin(1:3) )
                fo    = recin(mtzlok(4))
                IF ( isNAN(fo) ) THEN
                    num_isNAN = num_isNAN + 1
                    CYCLE
                ENDIF
                sfo   = recin(mtzlok(5))
                IF ( isNAN(sfo) ) THEN
                    num_isNAN = num_isNAN + 1
                   CYCLE
                ENDIF

!                Rfree_flag = recin(mtzlog(8))
                Rfree_flag = 1
 
!               Reject reflection if its amplitude out of range:
                IF ( fo < cards_1%fo_low ) THEN
                    num_fo_low = num_fo_low + 1
                    CYCLE
                ENDIF

                IF ( fo > cards_1%fo_high ) THEN
                    num_fo_high = num_fo_high + 1
                    CYCLE
                ENDIF

!               Reject reflection on the basis of sigma:
                IF ( sfo == 0.0 ) THEN
                    num_unobs = num_unobs + 1
                ELSE IF ( fo < cards_1%sigma_cutoff * sfo ) THEN
                    num_sigma_cut_off = num_sigma_cut_off + 1
                    CYCLE
                ENDIF    

!               Indices predefined by lsprgi array:
                IF ( cards_1%phased_translation_function ) THEN
                    IF ( ANY ( mtzlok(1:7) == 0 ) ) THEN
                        WRITE(*,"(' mtzlok=',7I4)") mtzlok(1:7)
                        CALL die('Zero pointers to columns detected.', 'read_mtz')
                    ENDIF
                    phio = recin(mtzlok(6))
                    fomo = recin(mtzlok(7))
                   
                ENDIF
!               WRITE(*,"(3I4, 2F10.2)") inhkl, fo, sfo

!               Assign hkl vector:
                hkl = NINT ( recin(1:3) )
                IF ( in_ASU( hkl, my_mtz%sp_group%laue_group_number ) /= 1 ) THEN
                    num_not_in_asu = num_not_in_asu + 1
!                   In general we don't care if the reflection in ASU or not since we gonna expand to P1 anyway
                ENDIF

!               phio:
                IF ( cards_1%phased_translation_function ) THEN
                    phio = ( pi / 180.0_wp ) * recin(my_mtz%lookup(6))
                ENDIF
                CALL centric_and_epsilon_tests ( hkl, centre, sysabs, eps, phicen, my_mtz%sp_group, mult )

                IF ( sysabs ) THEN
                    num_sysabs = num_sysabs + 1
                    CYCLE
                ENDIF

!               centric reflection:
                IF ( centre )  THEN
                    eps = -eps
                    num_centric_reflns = num_centric_reflns + 1
                    IF ( cards_1%phased_translation_function ) THEN

                        IF ( COS(2.0_wp * (phio - phicen)) < 1.0_wp - 0.1_wp**13 &
                             .AND. .NOT. logmss(my_mtz%lookup(6)) ) THEN

                            num_centric_err_reflns = num_centric_err_reflns + 1

                            IF ( MOD ( num_centric_err_reflns, 100 ) == 1 ) THEN
                                WRITE(*,*) ' phio_old=', phio
                                WRITE(*,"(' READ_MTZ> ', 'Correcting phase ', F8.5, ' for centric reflection:', 3I4, &
                                         &' with phicen=', F8.5)") phio, NINT ( recin(1:3) ), phicen
                            ENDIF

                            phio = NINT( (phio - phicen) / pi ) * pi + phicen

                            IF ( MOD ( num_centric_err_reflns, 100 ) == 1 ) THEN
                                WRITE(*,*) ' phio_new=', phio
                            ENDIF

                        ENDIF
                    ENDIF
                ENDIF

!               fom section:
                IF ( cards_1%phased_translation_function ) THEN
                    fomo = REAL ( recin(my_mtz%lookup(7)), KIND=wp )

!                   Correct fomo:
                    IF ( fomo > 1.0_wp .OR. fomo < 0.0_wp ) THEN
                        WRITE ( *,"(' READ_MTZ> ', 'Error in input mtz FOM must be in 0..1 range for reflection :',&
                                  &3I4, F5.3)" )  NINT(recin(1:3)), fomo
                        CALL messag('Correcting fom...','read_mtz')
                        fomo = MIN(fomo, 1.0_wp)
                        fomo = MAX(fomo, 0.0_wp)
                        IF ( fomo == 0.0_wp ) THEN
                            num_fomo_zero = num_fomo_zero + 1
                            CYCLE
                        ENDIF
                    ENDIF

!                   Too low resulting amplitude:
                    IF ( fomo * fo < cards_1%fo_low ) THEN
                        num_fomo_very_low = num_fomo_very_low + 1
                        CYCLE
                    ENDIF
                ENDIF

                num_passed = num_passed + 1

                final_ok:IF ( run == 1 ) THEN
                    my_mtz%hkl(num_passed)   = hkl
                    my_mtz%s(num_passed)     = my_mtz%RECIP_TENSOR * hkl ! This is better
                    my_mtz%eps(num_passed)   = eps
                    my_mtz%fo(num_passed)    = fo
                    my_mtz%sigfo(num_passed) = ABS(sfo)
                    my_mtz%mult(num_passed)   = mult

!                   For R-free, CC-free:
                    IF ( ALLOCATED ( my_mtz%Rfree_flag ) ) THEN
                        my_mtz%Rfree_flag(num_passed) = Rfree_flag 
                    ENDIF

                    IF ( ALLOCATED(my_mtz%phio) ) THEN
                        my_mtz%phio  ( num_passed ) = phio
                    ENDIF

                    IF ( ALLOCATED(my_mtz%fomo ) ) THEN
                        my_mtz%fomo  ( num_passed ) = fomo
                    ENDIF
!
                    IF ( ABS(eps) > 1.01_wp ) THEN
                        num_eps = num_eps + 1
                    ENDIF

!                   Nmod is definied after 0 run only:
                    IF ( MOD ( num_passed, nmod ) == 0 .OR. num_passed == 1 ) THEN
                        inhkl = my_mtz%hkl(num_passed)
                        IF ( cards_1%phased_translation_function ) THEN
                            
                            WRITE ( *,"(' READ_MTZ> ', 1X, 3I4, 8F10.2)") inhkl, my_mtz%fo(num_passed),&
                                                                          my_mtz%sigfo(num_passed), &
                                                                          180.0/pi*my_mtz%phio(num_passed),&
                                                                          my_mtz%fomo(num_passed), &
                                                                          my_mtz%mult(num_passed)
                        ELSE

                            WRITE ( *,"(' READ_MTZ> ', 1X, 3I4, 8F10.2)") inhkl, &
                                                                          my_mtz%fo(num_passed),&
                                                                          my_mtz%sigfo(num_passed), &
                                                                          my_mtz%eps(num_passed), &
                                                                          my_mtz%mult(num_passed)
                        ENDIF
                    ENDIF
                ENDIF final_ok

            ENDDO read_hkl  
            my_mtz%estimated_number_of_reflections = num_passed
    
            IF ( run == 0 ) THEN
                CALL messag('Estimated number of reflections in the mtz file= '//TRIM(int_to_c(num_passed)), 'read_mtz')
            ENDIF

!           Rewind mtz file:
            CALL lrrewd(mtzin)

!       End of file:
        ENDDO double_run

!       Close file:
        CALL lrclos ( mtzin )

!       Print statistics:
        WRITE(*,"(' READ_MTZ> ', 'final number of reflections          =', I8)") num_passed
        WRITE(*,"(' READ_MTZ> ', 'number of centrosymmetric reflections=', I8)") num_centric_reflns
        WRITE(*,"(' READ_MTZ> ', 'reflections not in ASU               =', I8)") num_not_in_ASU
        WRITE(*,"(' READ_MTZ> ', 'systematically absent                =', I8)") num_sysabs
        WRITE(*,"(' READ_MTZ> ', 'hi resol reflections                 =', I8)") num_hi_resol_reflns
        WRITE(*,"(' READ_MTZ> ', 'reflections in epsilon zone          =', I8)") num_eps
        WRITE(*,"(' READ_MTZ> ', 'reflections with zero fom            =', I8)") num_fomo_zero
        WRITE(*,"(' READ_MTZ> ', 'reflections with low amplitude       =', I8)") num_fo_low
        WRITE(*,"(' READ_MTZ> ', 'reflections with high amplitude      =', I8)") num_fo_high
        WRITE(*,"(' READ_MTZ> ', 'rejected on low sigma basis          =', I8)") num_sigma_cut_off
        WRITE(*,"(' READ_MTZ> ', 'NaN records                          =', I8)") num_isNaN

        IF ( ALLOCATED ( my_mtz%fomo ) ) THEN
            WRITE(*,"(' READ_MTZ> ', 'rejected on low fom basis            = ',I6)") &
            num_fomo_very_low
        ENDIF

        IF ( ALLOCATED ( my_mtz%phio ) ) THEN
            WRITE(*,"(' READ_MTZ> ', 'phases corrected for centric refl    = ',I6)")&
            num_centric_err_reflns
        ENDIF

        my_mtz%current_number_of_reflections = num_passed

!       Last check:
        IF ( SIZE ( my_mtz%hkl ) /= my_mtz%current_number_of_reflections ) THEN
            WRITE(*,*) SIZE ( my_mtz%hkl ), my_mtz%current_number_of_reflections
            CALL die('Programming error: length of mtz hkl array does not match current_number_of_reflections',&
                     'read_mtz')
        ENDIF

        IF ( SIZE ( my_mtz%hkl ) == 0 ) THEN
            WRITE(*,*) SIZE ( my_mtz%hkl )
            CALL die('User error. No reflections found in this resolution range.', 'read_mtz')
        ENDIF

!       Calculate optimal parameters for accurate structure factor calculation:
        CALL messag('We are using an approximation of Navaza''s algorithm...', 'read_mtz')

        ort_diag = extract_diagonal ( my_mtz%ORT )

!       Test various grid factors using Navaza's method until suitable Shannon rate found:
        DO i = 6, 48
            my_mtz%grid_factor = i / 2.0_wp
            WRITE(*,"(' READ_MTZ> ', 60('-'))")
            WRITE(*,"(' READ_MTZ> ', 'Testing grid factor=', F6.1, ' Shannon rate=', F6.2)")&
            my_mtz%grid_factor, my_mtz%grid_factor / 2.0_wp

!           Optimal sigma' calculcation a la NAVAZA:
            CALL optimal_sigma(my_mtz, my_mtz%grid_factor, my_mtz%fft_radius, my_mtz%badd)

!           Report results:
            WRITE (*,"(' READ_MTZ> ', 'Setting optimal integration radius= ', F8.2)") my_mtz%fft_radius
            WRITE (*,"(' READ_MTZ> ', 'Effective cell params in orthogonal system= ', 3F8.2)" ) ort_diag 

!           Radius will assume values close or less than MY_MTZ%FFT_RADIUS (no need to increase or reduce it here):

!           A factor of 2 is necessary if we going to calculate diagonal matrix elements using old approaches:
            IF ( 2.0_wp  * my_mtz%fft_radius > 0.5 * MINVAL ( ort_diag ) ) THEN

!               We do not want a sphere occupying more that half of unit cell in any direction:
!               This may be the case when calculating maps containing (2 pi i h) coefficients, then
!               periodicity of the cell is effectively reduced by a factor of 2: f(x) = -f(x+T/2)
                WRITE(*,"(' READ_MTZ> ', 'scaled fft_radius of ', F6.2, ' should be less than ', F6.2)") &
                2.0_wp * my_mtz%fft_radius,  0.5 * MINVAL ( ort_diag )
                CALL messag ('This sphere won''t fit into this cell. Changing Shannon factor.','read_mtz')
                IF ( debug > 20 ) THEN
                    CALL messag('Going to sleep a little...', 'read_mtz')
                    CALL sleep(2)
                ENDIF
            ELSE
                EXIT
            ENDIF
        ENDDO

!       Set b_scale as well for safety:
        my_mtz%b_scale = my_mtz%badd
        WRITE ( *, "(' READ_MTZ> ', 'Setting anti-alias b-factor=         ', F8.2)") my_mtz%badd
        CALL messag ('Still need to reduce it by value of min Biso', 'read_mtz')
        WRITE ( *, "(' READ_MTZ> ', 'Shannon rate=                        ', F8.2)") my_mtz%grid_factor / 2
        CALL messag ( '.mtz file has been read successfully.', 'read_mtz' )
        cards_1%grid_factor = my_mtz%grid_factor

    END SUBROUTINE read_mtz

!   Adapted loosely from f2mtz.f
!   5/3/2004: Updated to show writing of crystals and datasets

!   This program writes a basic MTZ file, using information
!   in the DATA statements. Obviously, you will want to replace
!   this with your own data.

!   It should be run as "mtz_write HKLOUT foo.mtz"

    SUBROUTINE write_mtz(mtz_1)
!       Should be used as a reference to certain params:
        TYPE(mtz),                INTENT(INOUT) :: mtz_1
!       Max number of labels
        INTEGER, parameter                   :: maxlab = 200
        INTEGER, PARAMETER                   :: MAXSYM = 192
        INTEGER                              :: mtzout
        INTEGER                              :: NSYM
        INTEGER                              :: NSYMP
!        INTEGER                              :: NLPRGI
        INTEGER, PARAMETER                   :: nlprgo = 11
        INTEGER                              :: NREFW
!        REAL(KIND=sp)                        :: wavelength = 0.0_sp
        REAL(KIND=sp), DIMENSION(6)          :: cell
        REAL(KIND=sp), DIMENSION(4,4,MAXSYM) :: rsym
        REAL(KIND=sp), DIMENSION(MAXLAB)     :: adata
!        CHARACTER(LEN=30), DIMENSION(MAXLAB) :: LSPRGI
!        CHARACTER(LEN=1),  DIMENSION(MAXLAB) :: CTPRGI
        CHARACTER(LEN=30), DIMENSION(MAXLAB) :: LSPRGO
        CHARACTER(LEN=1),  DIMENSION(MAXLAB) :: CTPRGO
        CHARACTER(LEN=10)                    :: PGNAM
        INTEGER                              :: ndatasets
        INTEGER                              :: isymlib
        INTEGER,           DIMENSION(3)      :: ihkl
!       Counters:
        INTEGER                              :: i
        INTEGER                              :: iappnd
        INTEGER                              :: iset

!       Example data from distributed toxd example
!        DATA NLPRGI,LSPRGI /34,'H','K','L','F','SIGF',195*' '/

!       Nice to have bulk solvent separately:
        DATA LSPRGO /'H','K','L','FP','SIGFP','FC','PHIC','FWT','PHWT','DELFWT','PHDELWT','FOM', 'PHCOMB', 187*' '/
        DATA CTPRGO /'H','H','H','F', 'Q',    'F', 'P',   'F',  'P',   'F',     'P',      'W',   'P',      187*' '/

        isymlib = get_next_io_unit()
        iappnd  = 0

!       Open mtz file:
        mtzout = 2 
        call lwopen(mtzout,'HKLOUT')

!       Write title to MTZ header
        call lwtitl (MTZOUT,' Output mtz file from FMLSQ',1)

!       Copy cell parameters:
        cell = mtz_1%cell
        IF ( ALL ( cell(4:6) < 5.0_wp ) ) THEN
            cell(4:6) = cell(4:6) * 180.0_wp / pi
        ENDIF

        IF ( debug > 10000 ) THEN
            WRITE(*,*) 'cell=',cell
        ENDIF

        CALL lwcell(mtzout, cell)

!       Store sort order to be used:
        CALL lwsort(mtzout, mtz_1%isort)

!       Store the symmetry in the mtz header:
        CALL MSYMLB(isymlib, mtz_1%sp_group%space_group_number, &
                             mtz_1%sp_group%space_group_name,   &
                             PGNAM, NSYMP, NSYM, RSYM)
        call lwsymm (MTZOUT,NSYM,NSYMP,RSYM,mtz_1%sp_group%lattice_type,mtz_1%sp_group%space_group_number,&
                     mtz_1%sp_group%space_group_name,PGNAM)
!
         
        CALL lwclab(mtzout, lsprgo, nlprgo, ctprgo,0)

! FIXME data sets
        ndatasets = mtz_1%ndatasets
        IF ( debug > 10000 ) THEN
            WRITE(*,*) ' ndatasets=',ndatasets
            WRITE(*,*) ' lookup=', mtz_1%lookup
        ENDIF

        IF ( ndatasets > 0 ) THEN
            DO iset = 1, ndatasets
                IF ( debug > 10000 ) THEN
                    WRITE(*,*) ' iset=',iset, mtz_1%pname(iset), ' ', mtz_1%xname(iset)
                    WRITE(*,*)  mtz_1%datcell(1:6,iset)
                ENDIF
                CALL lwidx(mtzout, mtz_1%pname(iset), mtz_1%xname(iset), mtz_1%dname(iset),&
                           mtz_1%datcell(1,iset), mtz_1%datwave(iset))
            ENDDO
            IF ( debug > 10000 ) THEN
                WRITE(*,*) mtzout, ' nlprgo=', nlprgo, ' iappnd=', iappnd
            ENDIF

!           Note that how, e.g. XNAME_OUT and XNAME have different dimensionality/size:
            CALL lwidasx(mtzout, nlprgo, mtz_1%xname_out, mtz_1%dname_out, iappnd)
        ENDIF
!   
!       Store program details in history header:
        call lwhstl (MTZOUT, 'FMLSQ is a standard full matrix refinement pgm.')

!       Loop over ouput reflections
!       Data for 2 reflections given here - normally calculated and
!       lots more!
        nrefw = SIZE ( mtz_1%hkl ) 
        do i = 1, nrefw

!           Initialise ADATA with Missing Number Flags
            CALL EQUAL_MAGIC(MTZOUT,ADATA,NLPRGO)

!           Fill in data columns where possible (all in this simple e.g.)
            ihkl =  mtz_1%hkl(i)
            adata(1:3) = ihkl(1:3)
!           Make substitution for 
            IF ( mtz_1%sigfo(i) > 0.0_wp ) THEN
                adata(4) = ABS ( mtz_1%fo(i) )
                adata(5) = ABS ( mtz_1%sigfo(i) )
            ELSE
                adata(4) = ABS ( mtz_1%fc(i) )
                adata(5) = 0.1 * adata(4)
            ENDIF

!           Fcalc (columns 6 and 7):
            adata(6) = ABS ( mtz_1%fc(i) )
            adata(7) = 180.0/pi * phase_of ( mtz_1%fc(i) )
            IF ( adata(7) < -180.0_wp ) adata(7) = adata(7) + 360.0_wp

!           2Fo-Fc difference Fourier:
            adata(8) = ABS ( 2.0_wp * mtz_1%fo(i) ) - ABS ( mtz_1%fc(i) )
            adata(9) = adata(7)
            IF ( adata(8) < 0.0 ) THEN
                adata(8) = -adata(8)
                adata(9) = adata(9) - 180.0
                IF ( adata(9) < -180.0_wp ) adata(9) = adata(9) + 360.0_wp
            ENDIF

!           Fo-Fc difference Fourier:
            adata(10) = ABS ( mtz_1%fo(i) ) - ABS ( mtz_1%fc(i) )
            adata(11) = adata(7)
            IF ( adata(10) < 0.0 ) THEN
                adata(10) = -adata(10)
                adata(11) = adata(11) - 180.0
                IF ( adata(11) < -180.0_wp ) adata(11) = adata(11) + 360.0_wp
            ENDIF

!           Write out reflection
            call lwrefl (MTZOUT, adata)

        ENDDO

!       Close mtz file, print full header:
        call lwclos (MTZOUT, 3)

    END SUBROUTINE write_mtz


!  Adapted from mtzmnf.f



!   This program reads a basic MTZ file, calculates a couple of new
!   columns based on the input data, and writes out a new file.

!   It should be run as:
!
!   mtz_rw HKLIN foo1.mtz HKLOUT foo2.mtz << eof
!   LABI F=filelabel SIGF=filelabel
!   LABO I=filelabel SIGI=filelabel
!   END
!   eof

    SUBROUTINE mtz_rw

!       .. Parameters ..
        INTEGER, PARAMETER :: MAXPAR=200
        INTEGER, PARAMETER :: MCOLS=200
        INTEGER, PARAMETER :: MAXSYM=96

        INTEGER                              :: MINDX
        INTEGER                              :: MTZPRT
        INTEGER                              :: MTZERR
        INTEGER,           DIMENSION(MCOLS)  :: LOOKUP
        INTEGER                              :: NLPRGI
        INTEGER                              :: JDO40
        INTEGER                              :: LENSTR
        INTEGER                              :: NLPRGO
        INTEGER                              :: IAPPND
        INTEGER                              :: NCOLX
        REAL(KIND=sp),     DIMENSION(MCOLS)  :: ADATA
        REAL(KIND=sp),     DIMENSION(MCOLS)  :: BDATA
        REAL(KIND=sp)                        :: S
        LOGICAL                              :: MTZEOF
        LOGICAL                              :: LOGMSS(MCOLS)
        CHARACTER(LEN=1),  DIMENSION(MCOLS)  :: CTPRGI
        CHARACTER(LEN=1),  DIMENSION(MCOLS)  :: CTPRGO
        CHARACTER(LEN=30), DIMENSION(MCOLS)  :: LSPRGI
        CHARACTER(LEN=30), DIMENSION(MCOLS)  :: LSPRGO

!       Parser variables:
        REAL(KIND=sp),     DIMENSION(MAXPAR) :: FVALUE
        INTEGER,           DIMENSION(MAXPAR) :: IBEG
        INTEGER,           DIMENSION(MAXPAR) :: IDEC
        INTEGER,           DIMENSION(MAXPAR) :: IEND
        INTEGER,           DIMENSION(MAXPAR) :: ITYP
        INTEGER                              :: NTOK
        LOGICAL                              :: LEND
        CHARACTER(LEN=4)                     :: KEY
        CHARACTER(LEN=400)                   :: LINE
        CHARACTER(LEN=4),  DIMENSION(MAXPAR) :: CVALUE

        DATA NLPRGI/5/
        DATA LSPRGI/'H','K','L','F','SIGF',195*' '/
        DATA CTPRGI/'H','H','H','F','Q',195*' '/
        DATA LOOKUP/-1,-1,-1,-1,-1,195*0/
        DATA NLPRGO/2/
        DATA LSPRGO/'I','SIGI',198*' '/
        DATA CTPRGO/'J','Q',198*' '/

!       CCP4 initialisations:
        call ccpfyp
!       MTZ-specific initialisations:
        call mtzini

!       Open input and output files on same index:
        MINDX = 1
        MTZprt = 1
        MTZERR = 0
        CALL LROPEN(MINDX,'HKLIN',MTZPRT,MTZERR)
        CALL LWOPEN(MINDX,'HKLOUT')

!       Parse keyworded input:
 10     CONTINUE
        LINE = '   '
        KEY = ' '
        NTOK = MAXPAR
!     
!       ***********************************************************
        CALL PARSER(KEY,LINE,IBEG,IEND,ITYP,FVALUE,CVALUE,IDEC,NTOK,LEND, .TRUE.)
!       ***********************************************************
   
!       End of file?
        IF (LEND) GO TO 50
      
        IF (KEY == 'LABI') THEN

!           Parse LABIN line
            CALL LKYIN(MINDX,LSPRGI,NLPRGI,NTOK,LINE,IBEG,IEND)

        ELSEIF (KEY == 'LABO') THEN

!           Parse LABOUT line
            CALL LKYOUT(MINDX,LSPRGO,NLPRGO,NTOK,LINE,IBEG,IEND)

        ELSE IF (KEY == 'END') THEN
            GO TO 50
        ELSE
            WRITE (6,FMT=6022) LINE(1:LENSTR(LINE))
 6022       FORMAT (' Error: Key_Word line NOT Understood:-',/' ',A)
            GO TO 10
        END IF
      
        GO TO 10

!       End of keyworded input:
   50   CONTINUE

!       Return no. of columns in input file:
        CALL LRNCOL(MINDX,NCOLX)
!       Assign input labels:
        CALL LRASSN(MINDX,LSPRGI,NLPRGI,LOOKUP,CTPRGI)
!       Assign output labels, appended to input labels:
        IAPPND = 1
        CALL LWASSN(MINDX,LSPRGO,NLPRGO,CTPRGO,IAPPND)

!       ---- Loop for reflections-----------
!
 100    CONTINUE

!       Read reflection data in file order:
        CALL LRREFL(MINDX,S,ADATA,MTZEOF)
        IF (MTZEOF) GO TO 200

!       Check for Missing Number Flags
        CALL LRREFM(MINDX,LOGMSS)

!       default is to leave columns unchanged:
        DO 40 JDO40 = 1,MCOLS
            BDATA(JDO40) = ADATA(JDO40)
 40     CONTINUE   

!       If either F or SIGF columns missing, mark new columns as missing:
        IF (LOGMSS(LOOKUP(4)) .OR. LOGMSS(LOOKUP(5))) THEN
            CALL EQUAL_MAGIC(MINDX,BDATA(NCOLX+1),2)
!           Else, write out new columns:
        ELSE
            BDATA(NCOLX+1) = ADATA(LOOKUP(4))**2
            BDATA(NCOLX+2) = (2.*ADATA(LOOKUP(5))*ADATA(LOOKUP(4)) + &
                                 ADATA(LOOKUP(5))*ADATA(LOOKUP(5)))
        ENDIF
 
        CALL LWREFL(MINDX,BDATA)

!       Return for next reflection:
        GOTO 100

 200    CONTINUE
        CALL LRCLOS(MINDX)
        CALL LWCLOS(MINDX,MTZPRT)

    END SUBROUTINE mtz_rw

    SUBROUTINE addopcol(label, ctype, lsprgo, ctprgo, nlprgo, index)
        CHARACTER(len= *),               INTENT(IN)    :: label
        CHARACTER(len= 1),               INTENT(IN)    :: ctype
        CHARACTER(len=30), DIMENSION(:), INTENT(INOUT) :: lsprgo
        CHARACTER(len= 1), DIMENSION(:), INTENT(INOUT) :: ctprgo
        INTEGER,                         INTENT(INOUT) :: nlprgo
        INTEGER,                         INTENT(OUT)   :: index

        nlprgo         = nlprgo+1
        index          = nlprgo
        lsprgo(nlprgo) = label
        ctprgo(nlprgo) = ctype
    END SUBROUTINE addopcol

!   FIXME Figure out whether we really need this:
    FUNCTION convert_symop_to_ort_system_mtz ( sym_1, mtz_2 )
        TYPE(symop)             :: convert_symop_to_ort_system_mtz
        TYPE(symop), INTENT(IN) :: sym_1
        TYPE(mtz), INTENT(IN)   :: mtz_2
!
        convert_symop_to_ort_system_mtz = ( mtz_2%ORT * sym_1 ) * mtz_2%DEORT
    END FUNCTION convert_symop_to_ort_system_mtz

    FUNCTION allocated_mtz ( mtz_1 )
        LOGICAL               :: allocated_mtz
        TYPE(mtz), INTENT(IN) :: mtz_1
        allocated_mtz = ALLOCATED ( mtz_1%hkl )
    END FUNCTION allocated_mtz

    SUBROUTINE deallocate_mtz ( mtz_1 )
!
!       Purpose:
!       =======
!       Deallocates MTZ type/structure
!
!       Date:                 Programmer:              History of changes:
!       ====                  ==========               ==================
!       Oct 2003              B.Strokopytov            Original code
!       Oct 2005              B.Strokopytov            ALLOCATABLE arrays
!       Oct 2007              B.Strokopytov            Rfree_flag arrays
!       Oct 2007              B.Strokopytov            Orthogonal Ho, Ko, Lo arrays     
!
!
        TYPE(mtz), INTENT(INOUT) :: mtz_1
!
        CALL messag(' ', 'deallocate_mtz')
        IF ( ALLOCATED(mtz_1%lookup) ) THEN
            CALL deallocate_array( mtz_1%lookup, 'lookup')
        ENDIF

!       Deallocate space group contents:
        IF ( ALLOCATED ( mtz_1%sp_group ) ) THEN
            CALL deallocate_space_group ( mtz_1%sp_group )
        ENDIF

        IF ( ALLOCATED ( mtz_1%datwave ) ) THEN
            CALL deallocate_array(mtz_1%datwave, 'datwave')
        ENDIF

        IF ( ALLOCATED ( mtz_1%datcell ) ) THEN
            CALL deallocate_array(mtz_1%datcell, 'datcell')
        ENDIF

        IF ( ALLOCATED ( mtz_1%xname ) ) THEN
            CALL deallocate_array(mtz_1%xname, 'xname')
        ENDIF
        
        IF ( ALLOCATED ( mtz_1%xname_out ) ) THEN
            CALL deallocate_array(mtz_1%xname_out, 'xname_out')
        ENDIF

        IF ( ALLOCATED ( mtz_1%pname ) ) THEN
            CALL deallocate_array(mtz_1%pname, 'pname')
        ENDIF

        IF ( ALLOCATED ( mtz_1%pname_out ) ) THEN
            CALL deallocate_array(mtz_1%pname_out, 'pname_out')
        ENDIF
        
        IF ( ALLOCATED ( mtz_1%dname ) ) THEN
            CALL deallocate_array(mtz_1%dname, 'dname')
        ENDIF

        IF ( ALLOCATED ( mtz_1%dname_out ) ) THEN
            CALL deallocate_array(mtz_1%dname_out, 'dname_out')
        ENDIF

        IF ( ALLOCATED ( mtz_1%hkl ) ) THEN
            CALL deallocate_array(mtz_1%hkl, 'hkl')
        ENDIF

        IF ( ALLOCATED( mtz_1%hkl_in_P1) ) THEN
            CALL deallocate_array(mtz_1%hkl_in_P1, 'hkl_in_P1')
        ENDIF

        IF ( ALLOCATED ( mtz_1%HO ) ) THEN
            CALL deallocate_array(mtz_1%HO, 'H0')
        ENDIF

        IF ( ALLOCATED ( mtz_1%KO ) ) THEN
            CALL deallocate_array(mtz_1%KO, 'K0')
        ENDIF

        IF ( ALLOCATED ( mtz_1%LO ) ) THEN
            CALL deallocate_array(mtz_1%LO, 'L0')
        ENDIF

        IF ( ALLOCATED( mtz_1%s ) ) THEN
            CALL deallocate_array(mtz_1%s, 's')
        ENDIF

        IF ( ALLOCATED( mtz_1%s_in_P1 ) ) THEN
            CALL deallocate_array(mtz_1%s_in_P1, 's_in_P1')
        ENDIF

        IF ( ALLOCATED( mtz_1%eps ) ) THEN
            CALL deallocate_array(mtz_1%eps, 'eps')
        ENDIF

        IF ( ALLOCATED( mtz_1%fo ) ) THEN
            CALL deallocate_array(mtz_1%fo, 'fo')
        ENDIF

        IF ( ALLOCATED( mtz_1%sigfo ) ) THEN
            CALL deallocate_array(mtz_1%sigfo, 'sigfo')
        ENDIF

        IF ( ALLOCATED( mtz_1%fomo ) ) THEN
            CALL deallocate_array(mtz_1%fomo, 'fomo')
        ENDIF

        IF ( ALLOCATED( mtz_1%mult) ) THEN
            CALL deallocate_array(mtz_1%mult, 'mult')
        ENDIF

        IF ( ALLOCATED( mtz_1%phio ) ) THEN
            CALL deallocate_array(mtz_1%phio, 'phio')
        ENDIF

        IF ( ALLOCATED( mtz_1%hla ) ) THEN
            CALL deallocate_array(mtz_1%hla, 'hla')
        ENDIF
        IF ( ALLOCATED( mtz_1%hlb ) ) THEN
            CALL deallocate_array(mtz_1%hla, 'hlb')
        ENDIF

        IF ( ALLOCATED ( mtz_1%hlc ) ) THEN
            CALL deallocate_array(mtz_1%hla, 'hlc')
        ENDIF
         
        IF ( ALLOCATED ( mtz_1%hld ) ) THEN
            CALL deallocate_array(mtz_1%hla, 'hld')
        ENDIF

        IF ( ALLOCATED ( mtz_1%a ) ) THEN
            CALL deallocate_array(mtz_1%a, 'a')
        ENDIF

        IF ( ALLOCATED ( mtz_1%b ) ) THEN
            CALL deallocate_array(mtz_1%b, 'b')
        ENDIF

        IF ( ALLOCATED ( mtz_1%fc ) ) THEN
            CALL deallocate_array(mtz_1%fc, 'fc')
        ENDIF

        IF ( ALLOCATED ( mtz_1%phic ) ) THEN
            CALL deallocate_array(mtz_1%phic, 'phic')
        ENDIF

        IF ( ALLOCATED ( mtz_1%fc_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%fc_in_P1, 'fc_in_P1')
        ENDIF

        IF ( ALLOCATED ( mtz_1%fcq_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%fcq_in_P1, 'fcq_in_P1')
        ENDIF

        IF ( ALLOCATED ( mtz_1%fcq_sym ) ) THEN
            CALL deallocate_array ( mtz_1%fcq_sym, 'fcq_sym' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%fc_temp ) ) THEN
            CALL deallocate_array ( mtz_1%fc_temp, 'fc_temp' )
        ENDIF
        IF ( ALLOCATED ( mtz_1%fcq_sym ) ) THEN
            CALL deallocate_array ( mtz_1%wphases, 'wphases' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%fo_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%fo_in_P1, 'fo_in_P1')
        ENDIF

        IF ( ALLOCATED ( mtz_1%weight ) ) THEN
            CALL deallocate_array ( mtz_1%weight, 'weight')
        ENDIF
      
        IF ( ALLOCATED ( mtz_1%weight_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%weight_in_P1, 'weight_in_P1' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%Rfree_flag ) ) THEN
            CALL deallocate_array ( mtz_1%Rfree_flag, 'Rfree_flag')
        ENDIF

        IF ( ALLOCATED ( mtz_1%Rfree_flag_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%Rfree_flag_in_P1, 'Rfree_flag_in_P1' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%sigfo_in_P1 ) ) THEN
            CALL deallocate_array(mtz_1%sigfo_in_P1, 'sigfo_in_P1')
        ENDIF
        IF ( ALLOCATED ( mtz_1%fobs ) ) THEN
            CALL deallocate_array(mtz_1%fobs, 'fobs') 
        ENDIF

        IF ( ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN 
            CALL deallocate_array ( mtz_1%fc_in_P1_init, 'fc_in_P1_init' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%fpart_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%fpart_in_P1, 'fpart_in_P1' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%fpart ) ) THEN
            CALL deallocate_array ( mtz_1%fpart, 'fpart')
        ENDIF

        IF ( ALLOCATED ( mtz_1%sc ) ) THEN
            CALL deallocate_array ( mtz_1%sc, 'sc' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%sc_in_P1 ) ) THEN
            CALL deallocate_array ( mtz_1%sc_in_P1, 'sc_in_P1' )
        ENDIF

        IF ( ALLOCATED ( mtz_1%rres ) ) THEN
            CALL deallocate_array ( mtz_1%rres )
        ENDIF

        IF ( ALLOCATED ( mtz_1%bin ) ) THEN
            CALL deallocate_array ( mtz_1%bin )
        ENDIF

        IF ( ALLOCATED ( mtz_1%ideal_hkl ) ) THEN
            CALL deallocate_array(mtz_1%ideal_hkl)
        ENDIF
 
        IF ( ALLOCATED ( mtz_1%ind ) ) THEN
            CALL deallocate_array(mtz_1%ind)
        ENDIF

!        IF ( ALLOCATED ( mtz_1%tag ) ) THEN
!            CALL deallocate_array(mtz_1%tag)
!        ENDIF

        IF ( ALLOCATED ( mtz_1%isort ) ) THEN
            CALL deallocate_array(mtz_1%isort)
        ENDIF

        IF ( ALLOCATED ( mtz_1%pbulk ) ) THEN
            CALL deallocate_array(mtz_1%pbulk, 'pbulk')
        ENDIF

        IF ( ALLOCATED ( mtz_1%ctprgo ) ) THEN
            CALL deallocate_array(mtz_1%ctprgo, 'ctprgo')
        ENDIF

        IF ( ALLOCATED ( mtz_1%lsprgo ) ) THEN
            CALL deallocate_array(mtz_1%lsprgo, 'lsprgo')
        ENDIF
        CALL messag('Done.', 'deallocate_mtz')

    END SUBROUTINE deallocate_mtz

    FUNCTION num_reflns_in_a_sphere ( mtz_1 )
!
!       Purpose:
!       =======
!       Calculates theoretically possible number of reflections in a sphere with r=1/D(min)
!
!       Note:
!       ====
!       Obsolete. Double read of MTZ is much safer.
!
        INTEGER               :: num_reflns_in_a_sphere
        TYPE(mtz), INTENT(IN) :: mtz_1
!       Local variables:
        REAL(KIND=wp)         :: s_max_squared
        REAL(KIND=wp)         :: s_min_squared
        REAL(KIND=wp)         :: s_squared
        INTEGER               :: il, ik, ih
        TYPE(vector_int)      :: hkl
        LOGICAL               :: centre
        LOGICAL               :: sysabs
        REAL(KIND=wp)         :: eps
        REAL(KIND=wp)         :: phicen
        INTEGER               :: numref
        INTEGER               :: num_sysabs

        CALL messag(' ', 'num_reflns_in_a_sphere')
        CALL messag('Estimating number of reflections in the ASU...', 'num_reflns_in_a_sphere') 

        s_max_squared     = 1.0_wp/(mtz_1%resol_max)**2
        IF ( mtz_1%resol_min > 0.0_wp ) THEN
            s_min_squared = 1.0_wp/(mtz_1%resol_min)**2
        ELSE
            s_min_squared = 0.0_wp
        ENDIF
        num_sysabs = 0
        numref     = 0
!
        DO il = -mtz_1%lmax, mtz_1%lmax
            DO ik = -mtz_1%kmax, mtz_1%kmax
                DO ih = -mtz_1%hmax, mtz_1%hmax

                    hkl       = (/ih, ik, il/)
                    s_squared = mtz_1%RECIP_TENSOR * hkl

                    IF ( s_squared  >= s_min_squared .AND. s_squared <= s_max_squared ) THEN

                        IF ( in_ASU( hkl, mtz_1%sp_group%laue_group_number ) == 1 ) THEN
                            CALL centric_and_epsilon_tests ( hkl, centre, sysabs, eps, phicen, &
                                                             mtz_1%sp_group )
                            IF ( .NOT. sysabs ) THEN
                                numref = numref + 1
                            ELSE
                                num_sysabs = num_sysabs + 1
                            ENDIF

                        ENDIF

                    ENDIF

                ENDDO
            ENDDO
        ENDDO

        WRITE(*, "(' NUM_REFLNS_IN_A_SPHERE> ', 'Number of systematic absences   (theory)= ', A)") &
        TRIM(int_to_c(num_sysabs))
        WRITE(*, "(' NUM_REFLNS_IN_A_SPHERE> ', 'Estimated number of reflections (theory)= ', A)") &
        TRIM(int_to_c(numref))
        CALL messag(' ', 'num_reflns_in_a_sphere')

        num_reflns_in_a_sphere = numref

    END FUNCTION num_reflns_in_a_sphere

    SUBROUTINE generate_ideal_index_list(mtz_1, x_resol_max)
        TYPE(mtz),                                   INTENT(INOUT)          :: mtz_1
        REAL(KIND=wp),                               INTENT(IN),   OPTIONAL :: x_resol_max
!       Local variables:
        TYPE(vector_int)                                                    :: hkl
        INTEGER,          DIMENSION(3)                                      :: ihkl
        REAL(KIND=wp)                                                       :: s_squared
        LOGICAL                                                             :: centre
        LOGICAL                                                             :: sysabs
        REAL(KIND=wp)                                                       :: eps
        REAL(KIND=wp)                                                       :: phicen
        INTEGER                                                             :: numref
        INTEGER                                                             :: num_sysabs
        REAL(KIND=wp)                                                       :: mult
        REAL(KIND=wp)                                                       :: resol_max
        REAL(KIND=wp)                                                       :: s_min_squared
        REAL(KIND=wp)                                                       :: s_max_squared
        INTEGER                                                             :: hmax
        INTEGER                                                             :: kmax
        INTEGER                                                             :: lmax
        INTEGER                                                             :: hmin
        INTEGER                                                             :: kmin
        INTEGER                                                             :: lmin
!       Counters:
        INTEGER                                                             :: i
        INTEGER                                                             :: il
        INTEGER                                                             :: ik
        INTEGER                                                             :: ih
        INTEGER                                                             :: run

        CALL messag(' ', 'generate_ideal_index_list')

!       Reserving the possibility of using alternative resol_max:
        IF ( PRESENT ( x_resol_max ) ) THEN
            resol_max = x_resol_max
        ELSE
            resol_max = mtz_1%resol_max
        ENDIF
        WRITE(*,*) ' Resolution:', resol_max

!       Lower resolution limit set to infinity:
        s_min_squared = 1.0_wp / 10.0_wp ** 12
        s_max_squared = 1.0_wp / resol_max ** 2

!       Calculate hmin, hmax, etc. limits:
        hmax = INT ( mtz_1%cell(1) / resol_max )
        kmax = INT ( mtz_1%cell(2) / resol_max )
        lmax = INT ( mtz_1%cell(3) / resol_max )

        hmin = -hmax
        kmin = -kmax
        lmin = -lmax

        WRITE(*,"(6I5)") hmin, hmax, kmin, kmax, lmin, lmax

        DO run = 0, 1
            num_sysabs = 0
            numref     = 0

!           Implicit sorting (a la CCP4) h is the slowest index:
            DO ih = hmin, hmax
                DO ik = kmin, kmax
                    DO il = lmin, lmax

                        hkl = (/ih, ik, il/)
                        IF ( in_ASU ( hkl, mtz_1%sp_group%laue_group_number ) /= 1 ) CYCLE

                        s_squared = mtz_1%RECIP_TENSOR * hkl

                        resolution:IF ( s_squared > s_min_squared .AND. s_squared <= s_max_squared ) THEN

!                            This is in disagreement with CCP4 ASU and hkl packing routines:
!                            IF ( .NOT. hemi_h ( hkl ) ) hkl = -hkl

                            CALL centric_and_epsilon_tests(hkl, centre, sysabs, eps, phicen, &
                                                           mtz_1%sp_group, mult)

                            IF ( .NOT. sysabs ) THEN

                                numref = numref + 1

                                IF ( run == 1 ) THEN
                                    mtz_1%ideal_hkl(numref) = hkl

                                    IF ( debug > 20 ) THEN
                                        IF ( numref < 1000 ) THEN
                                            WRITE(*,"(' HKL ',3I4, F5.1, F10.2)" ) &
                                                        ih, ik, il, mult, SQRT ( 1.0_wp / s_squared )
                                        ENDIF
                                    ENDIF

                                ENDIF

                            ELSE

                                num_sysabs = num_sysabs + 1

                            ENDIF

                        ENDIF resolution

                    ENDDO
                ENDDO
            ENDDO

            IF ( run == 0 ) THEN
                WRITE(*, "(' GENERATE_IDEAL_INDEX_LIST> ', 'Number of systematic absences (theory)= ', A)" )  &
                TRIM ( int_to_c ( num_sysabs ) )
                WRITE(*, "(' GENERATE_IDEAL_INDEX_LIST> ', 'Estimated number of reflections (theory)= ', A)" ) &
                TRIM ( int_to_c ( numref ) )
                CALL messag ( 'Allocating appropriate arrays for hkl listing...', 'generate_ideal_index_list' )
                CALL allocate_array ( mtz_1%ideal_hkl,  numref, 'ideal_hkl' )
                CALL messag ( ' ', 'generate_ideal_index_list' )
            ENDIF

        ENDDO

        CALL messag ( 'Ideal list of HKLs has been generated successfully.', 'generate_ideal_index_list' )
        IF ( .NOT. PRESENT ( x_resol_max ) ) THEN
            WRITE(*,"(' GENERATE_IDEAL_INDEX_LIST> ', 'Completeness of the data set is:', F7.1,' %')") &
                   100.0_wp * SIZE ( mtz_1%hkl ) / REAL ( numref )
        ENDIF

!        This seems to be an unnecessary sophistication:
!        CALL pack_ideal_hkl_array(mtz_1)
!        CALL sort_and_unpack_ideal_hkl_array(mtz_1)
        IF ( debug > 10000 ) THEN
            DO i = 1, MIN ( 100, SIZE ( mtz_1%ideal_hkl ) )
                ihkl = mtz_1%ideal_hkl(i)
                WRITE(*,"(' HKL ',3I4)" ) ihkl
            ENDDO
        ENDIF
        CALL messag ( ' ', 'generate_ideal_index_list')

    END SUBROUTINE generate_ideal_index_list

    SUBROUTINE pack_ideal_hkl_array(mtz_1)
        TYPE(mtz), INTENT(INOUT)  :: mtz_1
!       Local variables:
        INTEGER, DIMENSION(3)     :: ihkl
        INTEGER                   :: nr
!       Counters:
        INTEGER                   :: i
        CHARACTER(LEN=32),   SAVE :: srname='pack_ideal_hkl_array'

!       Checkz:
        IF ( .NOT. ALLOCATED ( mtz_1%ideal_hkl ) ) THEN
            CALL die('Programming error. Array IDEAL_HKL has not been allocated.', srname)
        ENDIF        

        nr = SIZE ( mtz_1%ideal_hkl )
        IF ( nr <= 0 ) THEN
            CALL die('Programming error. Array IDEAL_HKL has zero size.', srname)
        ENDIF

        IF ( .NOT. ALLOCATED ( mtz_1%ind ) ) THEN
            CALL allocate_array(mtz_1%ind, nr )
        ENDIF

        DO i = 1, nr
            ihkl = mtz_1%ideal_hkl(i)
            IF ( ihkl(3) >= 0 ) THEN
                mtz_1%ind(i) = (ihkl(3)*(2*mtz_1%kmax+1) + (ihkl(2) + mtz_1%kmax))*(2*mtz_1%hmax+1) &
                             + (ihkl(1) + mtz_1%hmax)
            ELSE
                WRITE(*,*) ihkl
                CALL die('Programming error. Negative L index.', srname)
            ENDIF
        ENDDO
    END SUBROUTINE pack_ideal_hkl_array 

    SUBROUTINE sort_and_unpack_ideal_hkl_array(mtz_1)
        TYPE(mtz),                   INTENT(INOUT)  :: mtz_1
!       Local array:
        TYPE(vector_int), DIMENSION(:), ALLOCATABLE :: sorted_hkl
!       Local variables:
        INTEGER,          DIMENSION(3)              :: ihkl
        INTEGER                                     :: nr
        INTEGER                                     :: itemp1
        INTEGER                                     :: itemp2
        INTEGER                                     :: hmax
        INTEGER                                     :: kmax
        INTEGER                                     :: lmax
!       Counters:
        INTEGER                                     :: i
        CHARACTER(LEN=32),                     SAVE :: srname='unpack_ideal_hkl_array'

!       Checkz:
        IF ( .NOT. ALLOCATED ( mtz_1%ideal_hkl ) ) THEN
            CALL die('Programming error. Array IDEAL_HKL has not been allocated.', srname)
        ENDIF        

        nr = SIZE ( mtz_1%ideal_hkl )
        IF ( nr <= 0 ) THEN
            WRITE(*,*) nr
            CALL die('Programming error. Array IDEAL_HKL has zero size.', srname)
        ENDIF

        IF ( .NOT. ALLOCATED ( mtz_1%ind ) ) THEN
            CALL allocate_array(mtz_1%ind, nr )
        ENDIF

        CALL sortag(mtz_1%ind)
        CALL allocate_array(sorted_hkl, nr)

        hmax = mtz_1%hmax
        kmax = mtz_1%kmax
        lmax = mtz_1%lmax
        
        DO i = 1, nr
            ihkl(3) = mtz_1%ind(i)/((2*kmax+1)*(2*hmax+1))
            itemp1  = ihkl(3)*(2*kmax+1)*(2*hmax+1)
            itemp1  = mtz_1%ind(i)-itemp1
            itemp2  = itemp1/(2*hmax+1)
            ihkl(2) = itemp2-kmax
            itemp2  = itemp2*(2*hmax+1)
            itemp2  = itemp1-itemp2
            ihkl(1) = itemp2-hmax

            sorted_hkl(i) = ihkl            
        ENDDO

        mtz_1%ideal_hkl = sorted_hkl

!       Free memory:
        CALL deallocate_array(sorted_hkl)

    END SUBROUTINE sort_and_unpack_ideal_hkl_array 

    SUBROUTINE reset_polar_axes ( mtz_1 )
!
!       Purpose:
!       =======
!       Resets polar axes to .FALSE.
!       Useful when switching from Patterson searches
!       to phased translation function searches
!
!
!       Array integer array of non_polar_axes is now of length three
!
!       Author:                 Date:                        Description of changes
!       ======                  ====                         ======================
!       B.Strokopytov           Mar 2004                     Original code
!       B.Strokopytov           May 2005                     Explicit deallocation added
!
        TYPE(mtz), INTENT(INOUT) :: mtz_1

!       Reallocate non_polar_axes array:
        IF ( ALLOCATED ( mtz_1%sp_group%non_polar_axes ) ) THEN
            CALL deallocate_array(mtz_1%sp_group%non_polar_axes)
        ENDIF

        CALL allocate_array(mtz_1%sp_group%non_polar_axes, 3, 'non_polar_axes')

!       Re-initialize
        mtz_1%sp_group%non_polar_axes = (/1, 2, 3/)
        mtz_1%sp_group%lpaxis         = (/.FALSE., .FALSE., .FALSE./)

    END SUBROUTINE reset_polar_axes

    SUBROUTINE get_polar_axes_from_mtz (mtz_1, lpaxis )
!
!       Purpose:
!       =======
!       recalculate polar axis info from mtz_1%sp_group
!
!       Date:           Programmer:         History of changes:
!       ====            ==========          ==================
!       Oct 2003        B.Strokopytov       Based on Alexei Vagin code.
!       Oct 2007        B.Strokopytov       Comments, cosmetic changes.
!
!
        TYPE(mtz),             INTENT(IN)  :: mtz_1
        LOGICAL, DIMENSION(:), INTENT(OUT) :: lpaxis
!       Local varibles:
        TYPE(vector)                       :: xyz
        TYPE(vector)                       :: xyzm
        REAL(KIND=wp), DIMENSION(3)        :: xyz_test
!       Counters
        INTEGER                            :: i
        INTEGER                            :: m

        lpaxis = .TRUE.

        IF ( mtz_1%sp_group%space_group_number /= 1 ) THEN

            xyz = (/0.13_wp, 0.17_wp, 0.19_wp/)

            DO m = 2, mtz_1%sp_group%number_of_symops

!               Multiply matrix of symm operators by vector:
                xyzm     = ( .SYMA. mtz_1%sp_group%SYM(m) ) * xyz
                xyz_test = xyzm - xyz

                DO i = 1, 3
                    IF ( ABS(xyz_test(i)) > 0.001_wp ) lpaxis(i) = .FALSE.
                ENDDO

            ENDDO

        ENDIF
    END SUBROUTINE get_polar_axes_from_mtz
 
    FUNCTION gridres ( my_mtz, grid_factor )
!
!       Purpose:
!       =======
!       Calculate grid required for a particular data resolution
!       Lifted from K. Cowtan
!       Modified for Fortran 95 by B. Strokopytov Oct 2003
!
        TYPE(vector_int)                      :: gridres
        TYPE(mtz), INTENT(IN)                 :: my_mtz
        REAL(KIND=wp), INTENT(IN), OPTIONAL   :: grid_factor
        REAL(KIND=wp)                         :: grid_fac
!
!       Local variables:
        TYPE(vector_int)                      :: grid_size0
!
!       Counters:
        INTEGER                               :: i
!
        IF (  PRESENT( grid_factor ) ) THEN
            grid_fac = grid_factor
        ELSE
            grid_fac = 3.0_wp
        ENDIF

        grid_size0 = (/ (INT(  grid_fac * my_mtz%cell(i) / my_mtz%resol_max ) + 1, i=1,3) /)
        gridres    = gridsize( grid_size0 , my_mtz%sp_group )

    END FUNCTION gridres

    SUBROUTINE get_vector_of_possible_origins ( txyz, grid_factor, mtz_1 )
!
!       Purpose:
!       =======
!       Returns vector txyz containing possible origins
!       for monoclinic and higher symmetry groups
!
!       Note:
!       ====
!       txyz returned in fractions of unit cell.
!
!       Caution:
!       =======
!       This subroutine is not quite that useful for P1
!       You get a very very long txyz vector
!       (Primitive group contains infinite number
!        of possible origins.)
!
        TYPE(vector),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: txyz
        REAL(KIND=wp),                  INTENT(IN)                 :: grid_factor
        TYPE(mtz),                      INTENT(IN)                 :: mtz_1
!       Local variables:
        INTEGER,          DIMENSION(3)                             :: limits
        TYPE(vector_int)                                           :: map_size
        INTEGER,          DIMENSION(3)                             :: v
        REAL(KIND=wp),    DIMENSION(3)                             :: xyz
        INTEGER                                                    :: mpeaks
!       Counters:
        INTEGER                                                    :: i
        INTEGER                                                    :: j
        INTEGER                                                    :: k
        INTEGER                                                    :: m
        INTEGER                                                    :: n

        IF ( mtz_1%sp_group%space_group_number == 1 ) THEN
            CALL die('This routine is not well suited for P1 sp. group calculations...',&
                     'get_vector_of_possible_origins')
        ENDIF

        map_size = gridres ( mtz_1, grid_factor )
        limits   = 1
        v        = map_size
        DO k = 1, 3
            IF ( mtz_1%sp_group%lpaxis(k) ) THEN
                limits(k) = v(k)
            ENDIF
        ENDDO

        mpeaks = limits(1) * limits(2) * limits(3) * mtz_1%sp_group%primitive_cheshire_symops_number

        IF ( mpeaks < 1 ) THEN
           WRITE(*,"(' GET_VECTOR_OF_POSSIBLE_ORIGINS> ', 'mpeaks= ', A)") TRIM(int_to_c(mpeaks))
           CALL die('Programming error. Unreasonably low mpeaks number.', 'find_T3_cross_vectors')
        ENDIF

!       txyz size should be equal to mpeaks:
        WRITE(*,*) ' mpeaks=', mpeaks
        CALL allocate_array(txyz, mpeaks, 'txyz')

        n = 0
!       BUG CORRECTED DEC 2006:
        DO m = 1, mtz_1%sp_group%primitive_cheshire_symops_number ! don't change this
            DO k = 0, limits(3) - 1
                DO j = 0, limits(2) - 1
                    DO i = 0, limits(1) - 1
                        n       = n + 1
                        txyz(n) = (/REAL(i, KIND=wp)/limits(1), REAL(j, KIND=wp)/limits(2), REAL(k, KIND=wp)/limits(3)/)
                        txyz(n) =  txyz(n) +  .SYMV. mtz_1%sp_group%sym_cheshire(m)

!                       Need this for print out:
                        xyz     =  txyz(n) 
                        WRITE(*, "(' GET_VECTOR_OF_POSSIBLE_ORIGINS> ', 'Possible origin #', A,&
                                 &' corresponding to [',A,'] euclidian operator in fractions of unit cell = ', 3F6.3)")&
                                 &TRIM(int_to_c(n)), TRIM(int_to_c(m)), xyz
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE get_vector_of_possible_origins

    SUBROUTINE optimal_sigma ( mtz_1, grid_factor, radius, badd )
!
!       Attempts to simulate Navaza's algorithm for optimal anti-alias factor
!       Further research is needed...
!
!       Date:              Programmer:                 History of changes:
!       ====               ==========                  ==================
!       Oct 2007           B.Strokopytov               Original code based on Navaza's method
!                                                      Acta Cryst. (2002)
!
!       Note:
!       ====
!       FIXME Need to test formation of my_tensor
!
        TYPE(mtz),     INTENT(IN)   :: mtz_1
        REAL(KIND=wp), INTENT(IN)   :: grid_factor
        REAL(KIND=wp), INTENT(OUT)  :: radius
        REAL(KIND=wp), INTENT(OUT)  :: badd
!       Local variables:
        REAL(KIND=wp),PARAMETER       :: accuracy=3.5_wp
        TYPE(matrix)                  :: my_tensor
        REAL(KIND=wp), DIMENSION(3)   :: uvw
        REAL(KIND=wp), DIMENSION(3)   :: hkl_max
        REAL(KIND=wp), DIMENSION(3)   :: hkl_max_scaled
        REAL(KIND=wp), DIMENSION(3)   :: hkl_temp
!       To choose reflections as close to dmin as possible:
        REAL(KIND=wp), DIMENSION(3)   :: s_limits
        REAL(KIND=wp)                 :: smin
        REAL(KIND=wp)                 :: smax
        REAL(KIND=wp)                 :: s
!       Navaza's s':
        REAL(KIND=wp)                 :: sigma
        REAL(KIND=wp)                 :: max_sigma
!       Counters:
        INTEGER                       :: i 
        INTEGER                       :: i_chosen
        INTEGER                       :: n1
        INTEGER                       :: n2
        INTEGER                       :: n3

!       This is correct definition of reciprocal metric tensor:
        my_tensor =  mtz_1%DEORT * TRANSPOSE ( mtz_1%DEORT )
        IF ( debug > 25 ) THEN
            CALL print_mat ( my_tensor )
        ENDIF

!       This is correct for any syngonies / groups:
        hkl_temp = extract_diagonal ( my_tensor )

!       Correct hkl_max finally:
        hkl_max  = INT ( SQRT ( 1.0_wp / mtz_1%resol_max ** 2 / hkl_temp ) )
        DO i = 1, 3
            hkl_max_scaled = 0.0_wp
            hkl_max_scaled(i:i) = hkl_max(i:i)            
            hkl_temp    = my_tensor * hkl_max_scaled
            s_limits(i) = DOT_PRODUCT ( hkl_max_scaled, hkl_temp )
        ENDDO

        IF ( debug > 35 ) THEN
            WRITE(*,"(' OPTIMAL_SIGMA> ', 's^2 limits=', 3F10.6)") s_limits
        ENDIF

!       Testing will be done using smin, smax:
        smin = MINVAL ( s_limits ) 
        smax = MAXVAL ( s_limits )

!       Increase/reduce limits by machine precision to make sure all vectors participate:
        smin = smin - smin * SQRT ( EPSILON ( 1.0_wp ) )
        smax = smax + smax * SQRT ( EPSILON ( 1.0_wp ) )

!       Report resolution limits:
        IF ( debug > 3 ) THEN
            WRITE (*,"(' OPTIMAL_SIGMA> ', 'Testing largest vectors between d(max)=', F10.6,&
       &               ' and d(min)=', F10.6, ' limits')") &
                        1.0 / SQRT ( smax ), 1.0 / SQRT ( smin )
        ENDIF

!       Initialise:
        i_chosen = 0
        max_sigma = -1.0_wp

!       Start looping in the direction of reciprocal cell parameters:
        DO n1 = -1, 1
            DO n2 = -1, 1
                DO n3 = -1, 1
                    uvw = (/ n1, n2, n3 /)
                    IF ( ALL ( uvw == 0 ) ) CYCLE

!                   Array by array multiplication:
                    hkl_max_scaled = hkl_max * uvw

!                   Prepare for matrix/array product:
                    hkl_temp = my_tensor *  hkl_max_scaled

!                   Estimate resolution:
                    s = DOT_PRODUCT ( hkl_max_scaled, hkl_temp   )
                    IF ( debug > 40 ) THEN
                        WRITE(*,"(' OPTIMAL_SIGMA> ', 'hkl_max_scaled= ', 3F8.1, ' d=',F8.3)")&
                        hkl_max_scaled, 1.0 / SQRT ( s )
                    ENDIF

                    IF ( s >= smin .AND. s <= smax ) THEN
                        i_chosen = i_chosen + 1
                        WRITE(*,"(' OPTIMAL_SIGMA> ', 'Reflection ', 3F8.1, ' with d(min)=', F8.3,' has been chosen...')")&
                        hkl_max_scaled, 1.0_wp/SQRT(s)

!                       Refine now:
                        CALL optimize_sigma(sigma, grid_factor, hkl_max_scaled, mtz_1)

!                       It's more convenient and safe to have max_sigma positive:
                        max_sigma = MAX ( max_sigma, ABS ( sigma ) )
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

!       Minimum acceptable sigma for ALL tested reflections:
        IF ( i_chosen > 0 ) THEN
            WRITE(*,"(' OPTIMAL_SIGMA> ', 'Minimal SIGMA=', F10.5, ' for ', A, ' tested reflections.' )")&
            max_sigma, TRIM ( int_to_c( i_chosen ) )
        ELSE
            WRITE(*,*) ' i_chosen=', i_chosen
            WRITE(*,"(' OPTIMAL_SIGMA> ', ' d(max)=', F6.2, ' d(min)=', F6.2)")& 
            1.0 / SQRT ( smin ), 1.0 / SQRT ( smax )
            CALL die ('Programming error. Failed to choose any reflections between dmax and dmin...',&
                      'optimal_sigma' )
        ENDIF
! Test  FIXME:
!        max_sigma = max_sigma + 0.2
        
         
!       Optimal cut_off radius corresponding to 10^-5 accuracy:
        radius = max_sigma * SQRT ( LOG ( 100.0_wp ) * 5.0_wp )
        badd   = 8.0_wp * pi ** 2 * max_sigma ** 2

    END SUBROUTINE optimal_sigma

    SUBROUTINE optimize_sigma ( sigma, grid_factor, hkl_real, mtz_1 )
!
!       Purpose:
!       ========
!
!       Refines optimal sigma' using LSQ method. From refined sigma' BADD and optimal
!       integration radius can be recovered. See Navaza (2002), Acta Cryst. A57, 451
!       
!       Date:              Programmer:                    History of changes:
!       ====               ==========                     ==================
!       Oct 2007           B.Strokopytov                  Original code
!
!       Note:
!       ====
!       my_tensor formation must be checked for non-orthogonal cells
!
!
        REAL(KIND=wp),               INTENT(OUT) :: sigma
        REAL(KIND=wp),               INTENT(IN)  :: grid_factor
        REAL(KIND=wp), DIMENSION(3), INTENT(IN)  :: hkl_real
        TYPE(mtz),                   INTENT(IN)  :: mtz_1
!       Local variables:
        REAL(KIND=wp)                            :: accuracy
        REAL(KIND=wp)                            :: factor
        TYPE(matrix)                             :: my_tensor
        REAL(KIND=wp), DIMENSION(3)              :: uvw
        REAL(KIND=wp), DIMENSION(3)              :: hkl_max
        REAL(KIND=wp), DIMENSION(3)              :: hkl_max_scaled
        REAL(KIND=wp), DIMENSION(3)              :: hkl_temp
        REAL(KIND=wp)                            :: sigma_init 
        REAL(KIND=wp)                            :: sum
        REAL(KIND=wp)                            :: temp 
!       Refinement:
        REAL(KIND=wp)                            :: e_theory
        REAL(KIND=wp)                            :: delta
        REAL(KIND=wp)                            :: f
        REAL(KIND=wp)                            :: exp_factor
        REAL(KIND=wp)                            :: dsumds
        REAL(KIND=wp)                            :: d2sumd2s
        REAL(KIND=wp)                            :: shift
!       Counters:
        INTEGER                                  :: i
        INTEGER                                  :: n1
        INTEGER                                  :: n2
        INTEGER                                  :: n3
        CHARACTER(LEN=32),               SAVE    :: srname = 'optimize_sigma'

!       Debugging checks:
        IF ( debug > 20 ) THEN
            WRITE(*,"(' OPTIMIZE_SIGMA> ', ' -- Starting refinement for hkl:', 3F8.1,' --')")&
            hkl_real 
        ENDIF

!       This starting point is suitable for all resolution ranges up to 0.04 A,
!       DO NOT ATTEMPT TO CHANGE THIS:
        IF ( mtz_1%resol_max < 0.5_wp ) THEN
            sigma_init = 0.01_wp
        ELSE
            sigma_init = 0.1_wp
        ENDIF

!       Navaza's eps(theory) and initialization of other parameters:
!        IF ( mtz_1%resol_max > 2.0_wp ) THEN
!            accuracy = 4.0_wp
!        ELSE
!            accuracy = 6.0_wp
!        ENDIF

        accuracy = MAX (4.0_wp, 6.0_wp / mtz_1%resol_max )
        e_theory  = 0.1_wp ** accuracy
            
        my_tensor = mtz_1%DEORT * TRANSPOSE ( mtz_1%DEORT )
        hkl_max   = grid_factor * (/mtz_1%hmax, mtz_1%kmax, mtz_1%lmax/)
        IF ( debug > 25 ) THEN
            CALL print_mat ( my_tensor )
        ENDIF

        i = 0
        sigma = sigma_init
        shift = 1.0_wp

!       Be sure to limit number of cycles in case of trouble:
        main_loop:DO WHILE ( ABS ( shift ) >  2.0 * ABS ( sigma ) * EPSILON (1.0_wp) .AND. i < 50 )
            sum      = 0.0_wp
            dsumds   = 0.0_wp
            d2sumd2s = 0.0_wp
            DO n1 = -1, 1
                DO n2 = -1, 1
                    DO n3 = -1, 1
                        uvw = (/ n1, n2, n3 /)
                        IF ( ALL ( uvw == 0 ) ) CYCLE

!                       Array by array multiplication:
                        hkl_max_scaled = hkl_max * uvw

!                       Prepare for matrix/array product:
                        hkl_temp = my_tensor * ( hkl_max_scaled + 2.0_wp * hkl_real )
                        factor   = DOT_PRODUCT ( hkl_max_scaled, hkl_temp )
                        factor   = 0.5_wp * twopi ** 2 *  factor
                        IF ( factor * sigma ** 2 > ABS ( MINEXPONENT ( 1.0_wp ) ) ) THEN
                            WRITE(*,*) ' factor=', factor, ' MINEXPONENT=', ABS ( MINEXPONENT( 1.0_wp ) )
                            CALL die (' Programming error. MINEXPONENT exceeded...',&
                                       srname )
                        ENDIF
                        IF ( debug > 15 ) THEN
                            WRITE(*,"(' hkl_temp=', 3F10.5, ' factor=',ES12.5)") hkl_temp,factor
                        ENDIF
 
                        exp_factor = EXP ( -factor * sigma ** 2 )
                      
                        sum        = sum + exp_factor
                        temp       = 2.0_wp * factor * sigma * exp_factor
                        dsumds     = dsumds - temp

!                       This could be removed since LSQ is much more safer:
                        d2sumd2s   = d2sumd2s - temp / sigma  + 4 * temp ** 2 / exp_factor
                    ENDDO
                ENDDO
            ENDDO

!           Always check value of SUM ( we are dealing with EXP functions here:
            IF ( sum < e_theory * EPSILON ( 1.0_wp ) ) THEN
                WRITE(*,*) ' sum=', sum, ' e_theory * eps=', e_theory * EPSILON(1.0_wp )
                CALL messag('SUM is too small to proceed further...', srname )
                CALL messag ('Correcting SIGMA...', srname )
                CALL messag ('Continue SIGMA refinement...', srname)
                i = 0
                sigma_init = sigma_init / 2
                shift = 1.0_wp
                sigma = sigma_init 
                WRITE(*,"(' OPTIMIZE_SIGMA> ', 'New value for sigma=', F9.5)") sigma
                CYCLE main_loop 
            ENDIF

!           Calculate difference and function value:
            delta = e_theory - sum
            f = delta ** 2

!           LSQ:
            shift = delta / dsumds

!           Second derivatives (unstable) if we are too far:
!           shift = delta * dsumds / (( dsumds ) ** 2 - delta * d2sumd2s )

            IF ( debug > 20 ) THEN

                IF ( i <= 10 .OR. MOD ( i,10 ) == 0 .OR. ABS ( shift ) < 100.0 * EPSILON ( 1.0_wp ) * sigma ) THEN
                    WRITE(*,*) ' cycle=', i
                    WRITE(*,*) ' f=', f, ' sigma=', sigma
                    WRITE(*,*) ' dsumds=', dsumds, ' d2sumd2s=', d2sumd2s
                    WRITE(*,*) ' shift=', shift
                ENDIF

            ENDIF

!           Apply shift to minimize target function:
            sigma = sigma + shift
            i = i + 1

        ENDDO main_loop

        IF ( debug > 20 ) THEN
            WRITE(*,"(' OPTIMIZE_SIGMA> ', 'SIGMA refinement converged after ', A, ' cycles')")&
            TRIM(int_to_c(i))
            WRITE(*,"(' OPTIMIZE_SIGMA> ', 'Final function value=', ES12.5, ' sigma=', F9.5)")&
            f, sigma
            WRITE(*,"(' OPTIMIZE_SIGMA> ', 'Final shift=', ES12.5)") shift
            WRITE(*,"(' OPTIMIZE_SIGMA> ', 'EPS(theory)=', ES12.5, ' sum=', ES12.5)")&
            e_theory, sum
            IF ( ABS ( shift ) < 2 * EPSILON ( 1.0_wp ) * ABS ( sigma ) ) THEN
                CALL messag ('Machine precision has been achieved for value of SIGMA', srname)
            ENDIF
        ENDIF

    END SUBROUTINE optimize_sigma 
END MODULE mtz_io
