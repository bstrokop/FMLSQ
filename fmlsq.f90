! -----------------------------------------------------------------------
!        FMLSQ - FULL MATRIX LEAST SQUARES refinement 
!                  (VERSION 0.1 21-OCT-2008)
!                  (LINUX PC VERSION)
!
!                  BY  BORIS V. STROKOPYTOv
!                      KONSTANTIN POLAYKOV
!                      PAVEL DOROVATOVSKY
!
!           DEPARTMENT OF SYNCHROTRON RADIATION
!                  KURCHATOV INSTITUTE
!                 MOSCOW, RUSSIA
!
!            TELEPHONE:  +7[499]  196-7895
!            TELEFAX:    +7[499]  135-1204
!            E-MAIL:     STROKOP@GMAIL.COM
!
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PROGRAM fmlsq 
USE basic_pdb_manip
USE lsq
USE ucards
USE refinement_util
IMPLICIT NONE
    TYPE(cards)        :: my_cards
    TYPE(pdb)          :: my_pdb
    TYPE(mtz)          :: my_mtz
    CHARACTER(LEN=256) :: pdb_file_in
    CHARACTER(LEN=256) :: pdb_file_out
    CHARACTER(LEN=256) :: pdb_dir
    CHARACTER(LEN=32)  :: srname
    LOGICAL            :: exists
    INTEGER            :: j
!   Refinement:
    LOGICAL            :: refine_xyz
    LOGICAL            :: refine_Biso
    LOGICAL            :: refine_occ
!   For the future:
    LOGICAL            :: refine_Us

    srname = 'fmlsq'

!   Start with mandatory CCP4 call
    CALL ccpfyp
    CALL read_cards ( my_cards ) ! need at least resmin, resmax 
    CALL ugtenv2 ( '-M', pdb_file_in )
    IF ( LEN_TRIM ( pdb_file_in ) == 0 ) THEN
        CALL ugtenv2 ('-PDBIN', pdb_file_in )
        IF ( LEN_TRIM ( pdb_file_in ) == 0 ) THEN
            CALL messag('Please specify "-m" or "-pdbin" argument...', srname)
            CALL die('Failed to read PDBIN file name.', srname)
        ENDIF
    ENDIF

!   Name has been read but we need to check whether such file really exists:
    CALL check_file_name(pdb_file_in, exists) 
    IF ( .NOT. exists ) THEN
        WRITE(*,"(' ROTLSQ> ', 'No such file name: ', A)") pdb_file_in
        CALL die ('Failed to find PDBIN file.', srname)
    ENDIF

    CALL ugtenv2 ( '-O', pdb_file_out )
    IF ( LEN_TRIM ( pdb_file_out ) == 0 ) THEN
        CALL ugtenv2 ('-PDBOUT', pdb_file_out )
        IF ( LEN_TRIM ( pdb_file_out ) == 0 ) THEN
            CALL messag('Please specify "-o" or "-pdbout" argument...', srname)
            CALL die('Failed to read PDBOUT file name.', srname) 
        ENDIF
    ELSE IF ( LEN_TRIM ( pdb_file_out ) > 255 ) THEN
        CALL die('Too long file name: '//pdb_file_out, srname)
    ENDIF

!   Need this ASAP to avoid problems after the heavy calculations:
    CALL check_directory_name(pdb_file_out, exists)

    CALL read_pdb ( my_pdb, pdb_file_in, my_cards )
    CALL which_params ( my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us )

    IF ( refine_occ ) THEN
        CALL refinable_occupancies(my_pdb)
    ENDIF

!   Moved matrix_pointers. Need polar_axes before that.

    IF ( debug > 3 ) THEN
        WRITE(*,"(' ROTLSQ> ', 'bmin=', F6.2)") get_bmin ( my_pdb )
    ENDIF

     CALL symmlq_refinement ( my_cards, my_pdb ) 
!!   CALL schur_complements_matinv ( my_cards, my_pdb )
!    CALL bekas_diagonal_estimate ( my_cards, my_pdb )
!    CALL lanczos_matinv( my_cards, my_pdb )
!    CALL qr_matinv( my_cards, my_pdb )
!    CALL estimate_standard_uncertainties ( my_cards, my_pdb )

    CALL write_pdb ( my_pdb, my_cards )
    CALL deallocate_cards ( my_cards )
    CALL deallocate_pdb ( my_pdb )
    CALL deallocate_mtz ( my_mtz )
    CALL dfftw_cleanup()
    CALL ccperr ( 0, ' *** Normal termination *** ' )

END PROGRAM fmlsq
