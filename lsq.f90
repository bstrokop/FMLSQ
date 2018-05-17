MODULE lsq
USE arpack
USE block_diagonal_util
USE corr_coef_manip
USE estimators
USE fast_distance
USE fast_genden
USE first_derivatives
USE hermitian_map_manip
USE implicit_matrix_inversion
USE linear_search
USE map_fft
USE mkl95_blas
USE mkl95_lapack
USE mtz_io
!USE second_derivatives
!USE statistics_complete
USE testing
USE ucards
IMPLICIT NONE
CONTAINS
    SUBROUTINE symmlq_refinement ( my_cards, my_pdb )
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz         
!       Number of atoms:
        INTEGER                                    :: n
!       Total number of parameters:                   
        INTEGER                                    :: np
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       SYMMLQ:
!       Need this for tests temporarily:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: w
!       For solution check:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: r1
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: x 
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: y
!       MAtrix conditioning:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: sk
!       Linear search:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: pold
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: xnew 
        LOGICAL                                    :: check
        REAL(KIND=wp)                              :: fold
        REAL(KIND=wp)                              :: fnew
!       SYMMLQ:
        INTEGER, PARAMETER                         :: nout = 6
        INTEGER                                    :: itnlim
        INTEGER                                    :: istop
        INTEGER                                    :: itn
        LOGICAL                                    :: checkA
!       Check symmetry of preconditioner:
        LOGICAL                                    :: goodb
        LOGICAL                                    :: precon
!       Slightly different preconditioning scheme:
        LOGICAL                                       :: myprecon
        REAL(KIND=wp)                                 :: shift
        REAL(KIND=wp)                                 :: rtol
        REAL(KIND=wp)                                 :: Anorm
        REAL(KIND=wp)                                 :: Acond
        REAL(KIND=wp)                                 :: rnorm
        REAL(KIND=wp)                                 :: ynorm
        REAL(KIND=wp)                                 :: r1norm
!       LAPACK:
        REAL(KIND=wp), EXTERNAL                       :: dnrm2
!       CPU vars:
        REAL(KIND=sp)                                 :: time0
        REAL(KIND=sp)                                 :: time1
        REAL(KIND=sp)                                 :: t0
        REAL(KIND=sp)                                 :: symmlq_time0
        REAL(KIND=sp)                                 :: symmlq_time1
!       Counters:
        INTEGER                                       :: i
        INTEGER                                       :: main_itns
        CHARACTER(LEN=32), SAVE                       :: srname = 'symmlq_refinement'
!       Sparse matrix:
!        TYPE(coo_matrix)                             :: Acoo
!        TYPE(csr_matrix)                             :: Acsr
        INTEGER(KIND=eb)                              :: nnz
!       Sparse matrix atom list:
        REAL(KIND=wp)                                 :: r_crit
        INTEGER,          DIMENSION(:,:), ALLOCATABLE :: atom_pair
        LOGICAL                                       :: diag_mode
        TYPE(anymatrix),  DIMENSION(:),   ALLOCATABLE :: ANYMAT       
!        TYPE(vector_int), DIMENSION(:),   ALLOCATABLE :: ideal_hkl
        INTEGER(KIND=eb) :: mkl_memstat
        INTEGER(KIND=fb) :: allocated_buffers
!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, &
                                                    refine_occ, refine_Us)

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 23 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

!       FIXME (test):
!        my_pdb%refinable(1:3,1) = .FALSE. 

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st atom/group along Y axis...', srname)
        ENDIF

!       Need to decide which parameters to refine in rigid body using POLAR_AXES calcns:
!       FIXME

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Set number of atoms in MY_PDB:
        n = SIZE ( my_pdb )

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Choosing r_crit within 2.5 - 5.0 A range:
        r_crit = MAX ( 2.0_wp * my_mtz%resol_max, 2.5_wp )
        r_crit = MIN ( r_crit, 5.0_wp )
! FIXME
!        r_crit = 1.85_wp
!        r_crit =  my_mtz%resol_max

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' SYMMLQ_REFINEMENT> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") &
                      xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )

!       Check whether we can do LSQ refinement without restraints:
        IF ( SIZE ( my_mtz%fo ) < np + 1 ) THEN
            WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'Number of reflections=', I6, ' number of params=', I8)") &
            SIZE ( my_mtz%fo ), np
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:  FIXME
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )

!       Need this for tests:
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) ) 
!       
!        my_cards%bw = 0.0_wp
        my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
        my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )

        IF ( my_cards%bw == 0.0_wp ) THEN
            my_mtz%weight_in_P1 = 1.0_wp
        ENDIF

        IF ( my_cards%sigma_weighting ) THEN
            CALL messag('Using bloody SIGMA weighting.', srname)
            my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
            my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       Allocate array for partial structure contibution:
        CALL allocate_array(my_mtz%fpart_in_P1, SIZE ( my_mtz%fc_in_P1 ) )

!       Allocate arrays for scaling:
        CALL allocate_array(my_mtz%sc_in_P1, SIZE ( my_mtz%fc_in_P1 ) )

!       This might be needed for statistics:
        CALL allocate_array(my_mtz%sc, SIZE ( my_mtz%fo ) )      ! original sp. group

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

!       Bulk solvent maps:
        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       ================ END OF COMMON PART FOR ALL MATRIX-VECTOR LSQ ROUTINES ====================

!       SYMMLQ array allocations:
        CALL allocate_array(g,  np)
        CALL allocate_array(r1, np)

!       Need this line below just for test (coould be removed later):
        CALL allocate_array(w,  np)

        CALL allocate_array(x,  np)
        CALL allocate_array(y,  np)

!       Conditioning based on main diagonal:
        CALL allocate_array(sk, np)

!       Linear search:
        CALL allocate_array(pold, np)
        CALL allocate_array(xnew, np)

!       Initialise parameters for SYMMLQ:
        checkA = .TRUE.
        goodb  = .FALSE.
!        precon = .FALSE.
        precon = .TRUE.
!       Another type of conditioning will be used if PRECON is .FALSE.:
        myprecon = .NOT. precon
        shift  = 0.0_wp
        itnlim = np
        rtol = EPSILON ( 1.0_wp ) ** 0.50_wp
!        Better accuracy for the future (rtol_power in user cards?):
!        rtol = EPSILON ( 1.0_wp ) ** 0.66666666_wp 
        diag_mode = .FALSE.
        IF ( .NOT. diag_mode ) THEN
            CALL residue_block_diagonal_pairs ( atom_pair, my_pdb, residues_in_block=1, nout = 6)
!            CALL atom_block_diagonal_pairs ( atom_pair, my_pdb, atoms_in_block=1, nout=6)
            nnz = calc_nnz ( my_pdb, atom_pair )

!           BUG CORRECTED APR 2010 (this was inside loop over refinement cycles without any deallocation): 
            CALL allocate_matrix ( my_pdb%coo, np, np, nnz )
            CALL allocate_matrix ( my_pdb%csr, np, np, nnz )
            IF ( debug > 100 ) THEN
                DO i = 1, np
                    WRITE(*,"(' parameter ', I7, ' belongs to block ', I6 )") i, my_pdb%pblock(i)
                ENDDO
            ENDIF
        ENDIF
!       CPU timing:
        CALL CPU_TIME(time0)

!       ------------ Refinement loops -----------
!        DO main_itns = 0, my_cards%ncycles - 1
        DO main_itns = 1, my_cards%ncycles 
           

!           b_min is necessary to adjust b_scale = badd - bmin:
            bmin = get_bmin ( my_pdb )
            CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!           Adjust bmin in my_mtz as well:
            my_mtz%b_scale = my_mtz%badd - bmin
            WRITE(*,"(//,1X,90('-'))")
            WRITE(*,"(' SYMMLQ_REFINEMENT> ','***** LEAST SQUARES REFINEMENT CYCLE NO. ', A,' RESULTS *****')")&
            TRIM ( int_to_c( main_itns ) )
            WRITE(*,"(1X,90('-'),/)")

!           Prepare atom images at the beginning of each cycle (depends on module FAST_GENDEN):
            CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
            
            IF ( .NOT. ALLOCATED ( my_pdb%all_atom_images ) ) THEN
                CALL die('O-o-o-o-p-s-s. Failed to allocate atom images. Not enough memory?', srname)
            ENDIF

!           This call has to go first before matrix calcns since we need phases:
            CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                               my_pdb, mode=my_cards%refinement_mode, fun=fold ) 
       

!           Calculate sparse matrix using (ntype+1)*ntype/2 FFTs (improved Tronrud method):
            CALL calculate_sparse_matrix(my_pdb%COO, my_mtz, my_maps, my_pdb, atom_pair, sk)

!           First copy to array of (any)matrices to perform explicit inversion of matrix blocks:
            CALL copy_coo_matrix_to_array_of_anymatrices ( my_pdb%COO, my_pdb, ANYMAT )

!           See if the deallocation of ANYMAT does occur: 
            CALL copy_array_of_inverted_anymatrices_to_coo_matrix(my_pdb%COO, my_pdb, ANYMAT) !xpower)

!           Copy to CSR type matrix (makes no sense to deallocate MY_PDB%COO here):
            my_pdb%csr = my_pdb%coo

            CALL CPU_TIME ( symmlq_time0 )     
            t0 = SECNDS(0.0)
            CALL symmlq (np, Ax, RESMIN_INV, g, shift, checkA, precon,               &
                         x, itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm, &
                         my_mtz, my_maps, my_pdb, my_cards%refinement_mode, sk, my_pdb%csr)

!           Report CPU time spent in linear solver:
            CALL CPU_TIME ( symmlq_time1 )
            WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'CPU time in linear solver (SYMMLQ)=', F10.2, ' s',&
                     &'  Elapsed time=',F10.2, ' s')") symmlq_time1 - symmlq_time0, SECNDS(t0)

            IF ( istop == 6 ) THEN
                CALL die('Programming error. Matrix is not symmetric...', srname)
            ENDIF

!           Check solution:
            IF ( debug > 1 ) THEN
                CALL Ax(np, y, x, my_mtz, my_maps, my_pdb, my_cards%refinement_mode)
                r1   = g - y
                r1norm = dnrm2 ( np, r1, 1 )

                DO i = 1, MIN ( 10, np )
                    WRITE(*,"(' b=', ES12.5, ' y=', ES12.5, ' diff= ', ES9.2)") &
                              g(i), y(i), r1(i)
                ENDDO

!           Print residual:
            WRITE(nout, 2000) r1norm
            ENDIF
       2000 FORMAT(/ ' Final residual =', 1p, e8.1)
            IF ( debug > 30 ) THEN
                WRITE(nout, 2100) (i, x(i), i = 1, np)
       2100     FORMAT(/ ' Solution  x' / 1p, 4(i6, e14.6))
            ENDIF

            CALL book_line_search(pold, fold, g, x, xnew, fnew, check, &
                                  my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                  my_pdb, my_cards%refinement_mode)
            IF ( check ) THEN
                CALL messag('Refinement stopped. Convergence has been achieved...', srname)
                EXIT
            ENDIF

        ENDDO

!       MKL dependent part (dummy routines might be needed for GFORTRAN compiler):
        WRiTE(*,"(1X,  'SYMMLQ_REFINEMENT> ', 'INTEL MKL USED ', A, ' MB of memory in ', A, ' allocated buffers')") &
        TRIM ( int_to_c ( mkl_memstat ( allocated_buffers ) / 2**20 ) ), TRIM( int_to_c ( Allocated_Buffers ) )
        CALL MKL_FREEBUFFERS
        IF ( mkl_memstat ( allocated_buffers ) > 0 ) THEN
            CALL warn('MKL library memory leak detected.', srname)
        ENDIF

        WRITE(*,"(//,1X,90('-'))")
        WRITE(*,"(' SYMMLQ_REFINEMENT> ','***** LEAST SQUARES REFINEMENT RESULTS BEFORE CYCLE NO. ', A,' ******')") &
                    TRIM ( int_to_c ( main_itns ) )
        WRITE(*,"(1X,90('-'),/)")

!       FIXME we need only stats here (no gradients??): to be fixed...
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                           my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 

!       FIXME mtz output needed (without holes):
!       Test mtz write:
        CALL generate_ideal_index_list(my_mtz)
        CALL move_fo_to_new_index_list(my_mtz)

        CALL generate_ideal_fcalc_in_true_space_group(my_mtz, my_maps(my_first_allocated_map), &
                                                      binary_map, logical_map, my_pdb)
        CALL write_mtz(my_mtz)

!       MKL dependent part (dummy routines might be needed for GFORTRAN compiler) and need to check environment FIXME:
        WRiTE(*,"(1X,  'SYMMLQ_REFINEMENT> ', 'INTEL MKL USED ', A, ' MB of memory in ', A, ' allocated buffers')") &
        TRIM ( int_to_c ( mkl_memstat ( allocated_buffers ) / 2**20 ) ), TRIM( int_to_c ( Allocated_Buffers ) )
        CALL MKL_FREEBUFFERS
        IF ( mkl_memstat ( allocated_buffers ) > 0 ) THEN
            CALL warn('MKL library memory leak detected.', srname)
        ENDIF
        CALL CPU_TIME(time1)

!       Massive deallocation:

!       Linear search:
        CALL deallocate_array(xnew)
        CALL deallocate_array(pold)

!       Gradient:
        CALL deallocate_array(g)

!       SYMMLQ:
        CALL deallocate_array(r1)
!        CALL deallocate_array(r2)
!        CALL deallocate_array(v)

!       This array is used for tests:
        CALL deallocate_array(w)

        CALL deallocate_array(x)
        CALL deallocate_array(y)

!       Conditioning:
        CALL deallocate_array(sk)

!       Deallocation of Maps, matrices etc.:
        IF ( ALLOCATED ( my_maps ) )                CALL deallocate_array(my_maps)

!       DO not deallocate most of the arrays in PDB structure except matrices and images:
        IF ( ALLOCATED ( my_pdb%COO ) )             CALL deallocate_matrix(my_pdb%COO)
        IF ( ALLOCATED ( my_pdb%all_atom_images ) ) CALL deallocate_array(my_pdb%all_atom_images)
        IF ( ALLOCATED ( ANYMAT ) )                 CALL deallocate_array(ANYMAT)

!       Bulk maps:
        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       my_pdb will be deallocated in MAIN program:         

        WRITE ( *, "(' SYMMLQ_REFINEMENT> Total CPU time for the whole run =', F10.2,' s')") time1 - time0
        WRITE ( *, "(' SYMMLQ_REFINEMENT> Normal termination.')")


    END SUBROUTINE symmlq_refinement 

    SUBROUTINE prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)
       TYPE(mtz),         INTENT(INOUT) :: my_mtz
       TYPE(map),         INTENT(INOUT) :: binary_map
       TYPE(map_logical), INTENT(INOUT) :: logical_map
!      Local variables:
       REAL(KIND=wp)                    :: b_scale
       REAL(KIND=wp)                    :: fine_grid_factor
       REAL(KIND=wp)                    :: grid_factor
       CHARACTER(LEN=32),         SAVE  :: srname='prepare_bulk_solvent_maps'

!       Allocate maps for aniso/bulk scaling but rearrange some parameters first:
        grid_factor = my_mtz%grid_factor
        b_scale = my_mtz%b_scale

!       Need finer grid than usual:
        fine_grid_factor = 4
        binary_map%grid_factor = fine_grid_factor
        my_mtz%grid_factor = fine_grid_factor
        my_mtz%b_scale = 0.0_wp

!       Allocate logical map:
        logical_map%grid_factor = my_mtz%grid_factor
        CALL allocate_map(logical_map, my_mtz)

!       Allocate binary map:
        CALL allocate_map(binary_map, my_mtz, 0)
        IF ( binary_map%map_size /= logical_map%map_size ) THEN
            CALL die('Programming error. Map sizes differ.', srname)
        ENDIF

!       Prepare FFTW plans for binary map:
        CALL messag('Preparing plans for bulk solvent map...', srname)
        CALL real_to_complex_plan_in_place(binary_map)
        CALL complex_to_real_plan_in_place(binary_map)
        CALL messag('Maps for aniso/bulk solvent scaling have been prepared...', srname)

!       Restore some very important MTZ parameters:
        my_mtz%grid_factor = grid_factor        
        my_mtz%b_scale = b_scale

    END SUBROUTINE prepare_bulk_solvent_maps

    SUBROUTINE destroy_bulk_solvent_maps(binary_map, logical_map)
        TYPE(map),         INTENT(INOUT) :: binary_map
        TYPE(map_logical), INTENT(INOUT) :: logical_map

!       Remove bulk maps:
        CALL destroy_real_fft_plan(binary_map)
        CALL destroy_complex_fft_plan(binary_map)
        IF ( ALLOCATED ( binary_map ) )  CALL deallocate_map(binary_map)
        IF ( ALLOCATED ( logical_map ) ) CALL deallocate_map(logical_map)

    END SUBROUTINE destroy_bulk_solvent_maps

    SUBROUTINE schur_complements_matinv ( my_cards, my_pdb )
        USE dense_matrix_util
        IMPLICIT NONE
!
!       Algorithm:
!       =========
!       Due to Isaac Schur (1875-1941)
!
!       T = A - C*B**(-1)*D
!
!       S = B - D*A**(-1)*C
!
!
!       Note that A and B symmetric, while C and D are not.
!       
!
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
!       Bulk:
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
!
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Total number of observations:
        INTEGER                                    :: nobs         
!       Total number of parameters:                   
        INTEGER                                    :: np
!       Type of refinement:
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: x
        REAL(KIND=wp)                              :: delta_r
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       CPU vars:
        REAL(KIND=sp)                              :: t0
        REAL(KIND=sp)                              :: time0
        REAL(KIND=sp)                              :: time1
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: iat
        INTEGER                                    :: j 
        INTEGER                                    :: k
        INTEGER                                    :: l
!       Matrices:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: A
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: B
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: C
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: D 
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: W
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: AINV
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: BINV
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: SS 
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: T 
        INTEGER                                    :: info
        INTEGER                                    :: r
        INTEGER                                    :: s
        INTEGER,      EXTERNAL                     :: mkl_progress
        CHARACTER(LEN=80), SAVE                    :: srname = 'schur_complements_matinv'

!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us)
        polar_axes = my_mtz%sp_group%lpaxis

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 33 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st group along Y axis...', srname)
        ENDIF

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )
        IF ( np < 2 ) THEN
            WRITE(*,*) ' np=', np
            CALL die('Too small a number of refinable parameters.', srname)
        ENDIF

!       Check whether we can do LSQ refinement without restraints:
        nobs = SIZE ( my_mtz%fo )
        WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', 'Number of reflections=', I6, ' number of params=', I8)") &
        nobs, np
        IF ( nobs < np + 1 ) THEN
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) )

        IF ( .NOT. my_cards%sigma_weighting ) THEN
            WRITE(*,*) ' B for weighting=', my_cards%bw
            my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
            my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )
        ELSE
             IF ( ANY ( my_mtz%sigfo == 0.0_wp ) ) THEN
                 CALL die('Some sigmas are equal to zero.', srname)
             ENDIF
             my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
             my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' SYMMLQ_REFINEMENT> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       SYMMLQ array allocations:
        CALL allocate_array(g,  np)  ! matvec product
        CALL allocate_array(x,  np)

!       Matrix array constants:
        r = np / 2
        s = np - r

!       Need A, C, D=R, B matrix allocation:
        CALL allocate_array (A, r, r)
        CALL allocate_array (B, s, s)
        CALL allocate_array (C, r, s)
        CALL allocate_array (D, s, r)

        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0)

!       Having b_min it is necessary to adjust b_scale = badd - bmin:
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)
!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       This call has to go first before matrix calcns since we need phases, scale and (perhaps) other params:
        CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                           my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 
       
!       CALL dv(matrix, g, my_mtz, my_pdb)

!       Initialize x:
        x = 0.0_wp

!       Calculate the whole matrix (this the most lengthy part of calculations):
        DO i = 1, np

            x(i) = 1.0_wp
            IF ( MOD ( i, 10 ) == 0 .OR. i == 1 ) THEN 
                 WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ','Current number of multipications= ', A)") TRIM ( int_to_c ( i ) )
                IF ( MOD ( i, 100 ) == 0 ) WRITE(*,*)
            ENDIF
       

            CALL Ax(np, g, x, my_mtz, my_maps, my_pdb, my_cards%refinement_mode)

!           Accumulate A, B, C, D matrices (need this to save memory):
            IF ( i <= r ) THEN
                A(1:r,i) = g(1:r)
                D(1:s,i) = g(r+1:r+s)
            ELSE
                C(1:r,i-r) = g(1:r)
                B(1:s,i-r) = g(r+1:r+s)
            ENDIF

            x(i) = 0.0_wp

        ENDDO
        CALL CPU_TIME(time1)
        WRITE(*,"(' SCHUR_COMPLEMENTS_MATINV> ', ' CPU time for matrix calcns=', F10.1, ' s Elapsed time=', F9.1, ' s')") &
        time1 - time0, SECNDS(t0)
        t0 = SECNDS(0.0)

!       Reserve G array for future use:
        CALL deallocate_array(x)
        CALL deallocate_mtz(my_mtz)
        CALL deallocate_array(my_maps)
        IF ( ALLOCATED ( my_pdb%COO ) ) CALL deallocate_matrix(my_pdb%COO)
        IF ( ALLOCATED ( my_pdb%all_atom_images ) ) CALL deallocate_array ( my_pdb%all_atom_images )

!       Can do some symmetry stats and improvement here:
!       FIXME
!        debug = 101
        IF ( debug > 100 ) THEN
            WRITE(*,*)
            WRITE(*,*) ' matrix A:'
            DO i = 1, r
                WRITE(*,"(1000ES10.2)") (A(i,j),j=1,r)
            ENDDO
            WRITE(*,*)
            WRITE(*,*) ' matrix B:'
            DO i = 1, s
                WRITE(*,"(1000ES10.2)") (B(i,j),j=1,r)
            ENDDO
            WRITE(*,*)
            WRITE(*,*) ' matrix C:'
            DO i = 1, r
                WRITE(*,"(1000ES10.2)") (C(i,j),j=1,s)
            ENDDO
            WRITE(*,*)
            WRITE(*,*) ' matrix D:'
            DO i = 1, s
                WRITE(*,"(1000ES10.2)") (C(i,j),j=1,r)
            ENDDO
            WRITE(*,*) ' max C-TRANSPOSE(D):', MAXVAL(C-TRANSPOSE(D))
            STOP 'STOP ABCD' 
        ENDIF

!       Calculate inversion of symmetric B:
        CALL allocate_array(BINV, s, s)

!       Preserve B for future use (1.5 x matrix memory will be used):
        BINV = B
        CALL messag('Cholesky decomposition of matrix B', srname)
        CALL potrf(BINV, uplo='U', info=info)  ! Cholesky decomposition
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix B failed.',srname)
        ENDIF

!       Final inverison of matrix B:
        CALL potri(BINV)

!       Make it really symmetric:
        CALL correct_potri(BINV,uplo='U')

        IF ( debug > 10000 ) THEN
            WRITE(*,*) ' Matrix Binv:'
            DO i = 1, 10
                WRITE(*,"(1000ES9.2)") (BINV(i,j),j=1,10)
            ENDDO
            WRITE(*,*) ' r,s=', r,s
        ENDIF

!       Calculate C * BINV: C is r-by-s and BINV is s-by-s (1.5 x matrix memory):
        CALL allocate_array(W, r, s)
        CALL gemm(C, BINV, W, 'N','N', 1.0_wp, 0.0_wp)
!        W=MATMUL(C,BINV)

!       No longer need that:
        CALL deallocate_array(BINV)

!       Memory 1.5 x (matrix memory) :
        CALL allocate_array(T, r, r)

!       W is r-by-s matrix, D is s-by-r matrix => T is r-by-r matrix:
        CALL gemm(W, D, T, 'N', 'N', 1.0_wp, 0.0_wp)
!        T=MATMUL(W,D)

        IF ( debug > 10000 ) THEN
            write(*,*) ' Matrix T:'
            DO i = 1, 10
                WRITE(*,"(1000ES9.2)") (T(i,j),j=1,10)
            ENDDO
        ENDIF
 
!       No longer need that:
        CALL deallocate_array(W)        

        T = A - T

!       Debug:
        IF ( debug > 10000 ) THEN
            WRITE(*,*) ' Matrix T:'
            DO i = 1, 10
                WRITE(*,"(1000ES9.2)") (T(i,j),j=1,10)
            ENDDO
        ENDIF

        CALL potrf(T, uplo='U', info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix T failed.', srname)
        ENDIF

        CALL potri(T)
!       No need to make T symmetric here:

!       Accumulate trace of T**(-1):
        DO i = 1, r
            g(i) = T(i,i)
        ENDDO
        CALL deallocate_array(T)
        
!       ---> Do the same for S matrix staring with A inversion:
        CALL potrf(A, uplo='U', info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix A failed.', 'potrf')
        ENDIF

        CALL potri(A)
        CALL correct_potri(A, uplo='U')

!       Just for the sake of clarity:
        CALL allocate_array(AINV, r, r)
        AINV = A
        IF ( debug > 10000 ) THEN
            WRITE(*,*) ' Checking AINV:'
            DO i = 1, 10
                WRITE(*,"(1000ES9.2)") (AINV(i,j),j=1,10)
            ENDDO
        ENDIF

!       No longer need matrix A:
        CALL deallocate_array(A)

!       D * A**(-1): D is s-by-r matrix, A**(-1) is r-by-r matrix => W is s-by-r matrix:
        CALL allocate_array(W, s, r)
        CALL gemm(D, AINV, W, 'N', 'N', 1.0_wp, 0.0_wp)
        CALL deallocate_array(AINV)
        CALL deallocate_array(D)

!       W * C: W is s-by-r matrix, C is r-by-s matrix => SS is s-by-s matrix:
        CALL allocate_array(SS, s, s)
        CALL gemm(W, C, SS, 'N', 'N', 1.0_wp, 0.0_wp)

!       Free memory:
        CALL deallocate_array(C)
        CALL deallocate_array(W)

        SS = B - SS
        CALL deallocate_array(B)

        IF ( debug > 10000 ) THEN
            WRITE(*,*) ' Matrix SS (must be symmetrical):'
            DO i = 1, 10
                WRITE(*,"(1000ES9.2)") (SS(i,j),j=1,10)
            ENDDO
        ENDIF

!       Time to invert SS:
        CALL potrf(SS, uplo='U', info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix S failed.', 'potrf')
        ENDIF

        CALL potri(SS)

!       We are interested only in diagonal (no need to call CORRECT_POTRI to make it symmetrical):
        DO i = 1, s
            g(i+r) = SS(i,i)
        ENDDO

        CALL deallocate_array(SS)        
        CALL CPU_TIME(time1)

!       Print result:
        CALL messag('Matrix diagonal has been calculated.', srname)
        WRITE(*,*) '  Parameter        Trace of inverse    Standard uncertainty'

        DO i = 1, np
            WRITE(*,*) i, g(i), SQRT ( f0 * g(i) / ( nobs - np ) )
        ENDDO

!       Accumulate statistics:
        my_pdb%trace = SUM ( g )
        my_pdb%sum_sigsq_all =  my_pdb%trace * f0 / ( nobs - np )
        WRITE(*,*) ' SUM of squared sigmas=', my_pdb%sum_sigsq_all

!       Initialize array of standard uncertainties:
        CALL allocate_array(my_pdb%su, 11, SIZE ( my_pdb ) )
        my_pdb%su = 0.0_sp
        my_pdb%avg_delta_r = 0.0_wp
        my_pdb%trace_xyz = 0.0_wp
        l = 0
        j = 0
        DO iat = 1, SIZE ( my_pdb )
            delta_r = 0.0_wp
            DO k = 1, 11            
                IF ( my_pdb%refinable(k,iat) ) THEN
                    l = l + 1
                    my_pdb%su(k,iat) = SQRT ( g(l) * f0 / ( nobs - np ) )
!                    WRITE(*,"(3I6,2ES11.4)") k,iat,l, my_pdb%su(k,iat), g(l)                    

!                   Special statistics for XYZ s.u.:
                    IF ( k <= 3 ) THEN
                        my_pdb%trace_xyz = my_pdb%trace_xyz + g(l)
                        delta_r = delta_r + g(l)
                        IF ( k == 3 ) THEN
                            delta_r = SQRT ( delta_r * f0 / (nobs - np) )
                            IF ( delta_r /= 0.0_wp ) THEN
                                j = j + 1
                                my_pdb%avg_delta_r = my_pdb%avg_delta_r + delta_r
                            ENDIF
                        ENDIF
                    ENDIF

                ENDIF
            ENDDO
            
            WRITE(*,"(' SUs:', 11F8.3,4X,11L1,2X,F8.3)") my_pdb%su(1:11,iat), my_pdb%refinable(1:11,iat), delta_r
        ENDDO

!       Print xyz statistics:
        my_pdb%sum_sigsq_xyz = my_pdb%trace_xyz * f0 / ( nobs - np )
        my_pdb%f0 = f0
        WRITE(*,*) ' SUM of squared sigmas for xyz=', my_pdb%trace_xyz * f0 / ( nobs - np )
        WRITE(*,*) ' average delta r error=',  my_pdb%avg_delta_r / j, ' for ', j, ' atoms'

!       Massive deallocation:

!       Gradient:
        IF ( ALLOCATED ( g ) ) CALL deallocate_array(g)

        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       MY_PDB will be deallocated in MAIN program:
        
        WRITE ( *, "(' SCHUR_COMPLEMENT_MATINV> Total CPU time for the whole run=       ', F10.1,' s')") time1 - time0
        WRITE ( *, "(' SCHUR_COMPLEMENT_MATINV> Total elapsed time for matrix inversion=', F10.2,' s')") SECNDS(t0)
        WRITE ( *, "(' SCHUR_COMPLEMENT_MATINV> Normal termination.')")

    END SUBROUTINE schur_complements_matinv

    SUBROUTINE qr_matinv ( my_cards, my_pdb )
        USE dense_matrix_util
        IMPLICIT NONE
!
!       Purpose:
!       =======
!       Simple matrix diagonalization.
!
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
!       Bulk:
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
!
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Total number of observations:
        INTEGER                                    :: nobs         
!       Total number of parameters:                   
        INTEGER                                    :: np
!       Type of refinement:
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: x
        REAL(KIND=wp)                              :: trace_xyz 
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       CPU vars:
        REAL(KIND=sp)                              :: t0
        REAL(KIND=sp)                              :: time0
        REAL(KIND=sp)                              :: time1
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: j 
        INTEGER                                    :: k
        INTEGER                                    :: l
!       Matrices:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: A
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: e
        INTEGER                                    :: info
        INTEGER,      EXTERNAL                     :: mkl_progress
        CHARACTER(LEN=80), SAVE                    :: srname = 'qr_matinv'

!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us)
        polar_axes = my_mtz%sp_group%lpaxis

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 33 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st group along Y axis...', srname)
        ENDIF

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' QR_MATINV> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )
        IF ( np < 2 ) THEN
            WRITE(*,*) ' np=', np
            CALL die('Too small a number of refinable parameters.', srname)
        ENDIF

!       Check whether we can do LSQ refinement without restraints:
        nobs = SIZE ( my_mtz%fo )
        WRITE(*,"(' QR_MATINV> ', 'Number of reflections=', I6, ' number of params=', I8)") &
        nobs, np
        IF ( nobs < np + 1 ) THEN
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' QR_MATINV> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' QR_MATINV> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' QR_MATINV> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) )

        IF ( .NOT. my_cards%sigma_weighting ) THEN
            WRITE(*,*) ' B for weighting=', my_cards%bw
            my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
            my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )
        ELSE
             IF ( ANY ( my_mtz%sigfo == 0.0_wp ) ) THEN
                 CALL die('Some sigmas are equal to zero.', srname)
             ENDIF
             my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
             my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' QR_MATINV> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

!       Introducing bulk solvent maps:
        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       SYMMLQ array allocations:
        CALL allocate_array(g,  np)  ! matvec product
        CALL allocate_array(x,  np)

        CALL allocate_array (A, np, np)
        CALL allocate_array (e, np)

        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0)

!       Having b_min it is necessary to adjust b_scale = badd - bmin:
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       This call has to go first before matrix calcns since we need phases, scale and (perhaps) other params:
        CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                       my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 
       
!       CALL dv(matrix, g, my_mtz, my_pdb)

!       Initialize x:
        x = 0.0_wp

!       Calculate the whole matrix (this the most lengthy part of calculations):
        DO i = 1, np

            x(i) = 1.0_wp
            IF ( MOD ( i, 10 ) == 0 .OR. i == 1 ) THEN 
                 WRITE(*,"(' QR_MATINV> ',' Multipication No.=', A)") TRIM ( int_to_c ( i ) )
            ENDIF
       

            CALL Ax(np, g, x, my_mtz, my_maps, my_pdb, my_cards%refinement_mode)

!           Accumulate A, B, C, D matrices (need this to save memory):
            A(1:np, i) = g(1:np)
            x(i) = 0.0_wp

        ENDDO
        CALL CPU_TIME(time1)
        WRITE(*,"(' QR_MATINV> ', ' CPU time for matrix calcns=', F10.1, ' s Elapsed time=', F9.1, ' s')") &
        time1 - time0, SECNDS(t0)
        t0 = SECNDS(0.0)

!       Reserve G array for future use:
        CALL deallocate_array(x)
        CALL deallocate_mtz(my_mtz)
        CALL deallocate_array(my_maps)
        IF ( ALLOCATED ( my_pdb%COO ) ) CALL deallocate_matrix ( my_pdb%COO )
        IF ( ALLOCATED ( my_pdb%all_atom_images ) ) CALL deallocate_array ( my_pdb%all_atom_images )

!       Can do some symmetry stats and improvement here:
!       FIXME
!        debug = 101
        IF ( debug > 100 ) THEN
            WRITE(*,*)
            WRITE(*,*) ' matrix A:'
            DO i = 1, MIN(10,np) 
                WRITE(*,"(1000ES10.2)") (A(i,j),j=1,MIN(10,np))
            ENDDO
            STOP 'STOP ABCD' 
        ENDIF

        CALL messag('Start QR decomposition of matrix A', srname)
        CALL syevd(A, e, jobz='V', UPLO='U', info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix A failed.',srname)
        ENDIF
        WRITE(*,"(' QR_MATINV> ', 'lambda (max)= ',  ES9.2, ' lambda (min)=', ES9.2,  ' ratio (condition number)=', ES9.2)")&
        e(np), e(1), e(np)/e(1)

!       Accumulate trace of T**(-1):
         DO  i = 1, np
             g(i) = DOT_PRODUCT ( A(i,1:np) / e(1:np), A(i,1:np) )
         ENDDO
        
!       No longer need matrix A:
        IF ( ALLOCATED ( A ) ) CALL deallocate_array(A)
        IF ( ALLOCATED ( e ) ) CALL deallocate_array(e)        

!       Print result:
        CALL messag('Matrix diagonal has been calculated.', srname)
        WRITE(*,*) '  Parameter        Trace of inverse    Standard uncertainty'

        DO i = 1, np
            WRITE(*,*) i, g(i), SQRT ( f0 * g(i) / ( nobs - np ) )
        ENDDO
        WRITE(*,*) ' SUM of squared sigmas=', SUM ( g ) * f0 / ( nobs - np ) 

        trace_xyz = 0.0_wp
        l = 0
        DO i = 1, SIZE ( my_pdb )
            DO k = 1, 3            
                IF ( my_pdb%refinable(k,i) ) THEN
                    l = l + 1
                    trace_xyz = trace_xyz + g(l)
                ENDIF
            ENDDO
            l = l + COUNT ( my_pdb%refinable(4:11,i) )
        ENDDO

!       Print xyz statistics:
        WRITE(*,*) ' SUM of squared sigmas for xyz=', trace_xyz * f0 / ( nobs - np )
        WRITE(*,*) ' average delta r error=',  SQRT ( trace_xyz * f0 / ( nobs - np ) / COUNT ( my_pdb%occ > 0.0_wp ) )

!       Massive deallocation:

!       Gradient:
        IF ( ALLOCATED ( g ) ) CALL deallocate_array(g)

!       Remove bulk maps:
        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       MY_PDB will be deallocated in MAIN program:
        
        WRITE ( *, "(' QR_MATINV> Total CPU time for the whole run=       ', F10.1,' s')") time1 - time0
        WRITE ( *, "(' QR_MATINV> Total elapsed time for matrix inversion=', F10.2,' s')") SECNDS(t0)
        WRITE ( *, "(' QR_MATINV> Normal termination.')")

    END SUBROUTINE qr_matinv

    SUBROUTINE lanczos_matinv ( my_cards, my_pdb )
        IMPLICIT NONE
!       INTEGER, PARAMETER :: lwp = SELECTED_REAL_KIND(p=23)
!
!       Purpose:
!       =======
!       Simple matrix tridiagonalization using Lanczos method with full orthogonalization.
!       Extremely suitable for crystallography. The major drawback is a need to keep
!       orthogonalized vectors in random access memory. Otherwise it would indispensable.
!       We could simply move orthogonal vectors to disk. But constant
!       orthogonalization spoiles everything.
!
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
!       Bulk solvent maps:
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
!
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Total number of observations:
        INTEGER                                    :: nobs         
!       Total number of parameters:                   
        INTEGER                                    :: np
!       Type of refinement:
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: x
        REAL(KIND=wp)                              :: trace_xyz 
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       CPU vars:
        REAL(KIND=sp)                              :: t0
        REAL(KIND=sp)                              :: t1
        REAL(KIND=sp)                              :: time0
        REAL(KIND=sp)                              :: time1
        REAL(KIND=wp)                              :: orthogonalization_time
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: j 
        INTEGER                                    :: k
        INTEGER                                    :: l
!       Matrices:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: O
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: alpha
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: beta
!       Potential lwp parameters:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: q 
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: r
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: w
        REAL(KIND=wp)                              :: beta1
        INTEGER                                    :: info
        INTEGER,      EXTERNAL                     :: mkl_progress
        CHARACTER(LEN=80), SAVE                    :: srname = 'lanczos_matinv'

!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us)
        polar_axes = my_mtz%sp_group%lpaxis

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 33 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st group along Y axis...', srname)
        ENDIF

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' LANCZOS_MATINV> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )
        IF ( np < 2 ) THEN
            WRITE(*,*) ' np=', np
            CALL die('Too small a number of refinable parameters.', srname)
        ENDIF

!       Check whether we can do LSQ refinement without restraints:
        nobs = SIZE ( my_mtz%fo )
        WRITE(*,"(' QR_MATINV> ', 'Number of reflections=', I6, ' number of params=', I8)") &
        nobs, np
        IF ( nobs < np + 1 ) THEN
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' LANCZOS_MATINV> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' LANCZOS_MATINV> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' LANCZOS_MATINV> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) )

        IF ( .NOT. my_cards%sigma_weighting ) THEN
            WRITE(*,*) ' B for weighting=', my_cards%bw
            my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
            my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )
        ELSE
             IF ( ANY ( my_mtz%sigfo == 0.0_wp ) ) THEN
                 CALL die('Some sigmas are equal to zero.', srname)
             ENDIF
             my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
             my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' QR_MATINV> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

!       Prepare for bulk solvent calculations:
        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       SYMMLQ array allocations:
        CALL allocate_array(g, np)
        CALL allocate_array(x, np)

        CALL allocate_array (O, np, np)
        CALL allocate_array (q, np)
        CALL allocate_array (r, np)
        CALL allocate_array (w, np)
        CALL allocate_array (alpha, np)
        CALL allocate_array (beta, np - 1)

        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0)

!       Having b_min it is necessary to adjust b_scale = badd - bmin:
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       This call has to go first before matrix calcns since we need phases, scale and (perhaps) other params:
        CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                           my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 
       
!       CALL dv(matrix, g, my_mtz, my_pdb)

!       Initialize vectors:
        q = 0.0_wp; q(1) = 1.0_wp; r = 0.0_wp
        O = 0.0_wp
        orthogonalization_time = 0 

!       Calculate the whole matrix (this the most lengthy part of calculations):
        DO i = 1, np

            IF ( MOD ( i, 10 ) == 0 .OR. i == 1 ) THEN 
                 WRITE(*,"(' QR_MATINV> ',' Multipication No.=', A)") TRIM ( int_to_c ( i ) )
            ENDIF
       
            O(1:np,i) = q
            CALL Ax(np, g, q, my_mtz, my_maps, my_pdb, my_cards%refinement_mode)
            r = r + g
            alpha(i) = DOT_PRODUCT (q, r)
            r = r - alpha(i) * q
            t1 = SECNDS(0.0)
            w = r

!           Need array slices here:
            CALL gemv(O(1:np,1:i), w, g(1:i), 1.0_wp, 0.0_wp, 'T')

! ANother possibility:
!            IF ( 1_eb * i * np > 10_eb ** 4 ) THEN
!OMP PARALLEL 
!                DO j = 1, i
!                    g(j) = DOT_PRODUCT (O(1:np,j), w)
!                ENDDO
!OMP END PARALLEL
!            ELSE
!                DO j = 1, i
!                    g(j) = DOT_PRODUCT (O(1:np,j), w)
!                ENDDO
!            ENDIF

            DO j = 1, i
!                r = r - DOT_PRODUCT (O(1:np,j), w) * O(1:np,j)
                r = r - g(j) * O(1:np,j)
            ENDDO

            orthogonalization_time = orthogonalization_time + SECNDS(t1)

            beta1 = SQRT ( DOT_PRODUCT ( r, r ) )

            IF ( i < np ) THEN
                beta(i) = beta1
            ENDIF

            IF ( beta1 == 0.0_wp ) EXIT
            q = r / beta1 
            r = -beta1 * O(1:np,i)
            WRITE(*,*) i, ' alpha=', alpha(i), ' beta=', beta1
        ENDDO

        CALL CPU_TIME(time1)
        WRITE(*,"(' LANCZOS_MATINV> ', ' CPU time for matrix calcns=', F10.1, ' s Elapsed time=', F10.2, ' s')") &
        time1 - time0, SECNDS(t0)
        t0 = SECNDS(0.0)

        CALL pttrf ( alpha, beta, info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('O-o-o-o-p-s-s. Possible programming error. Inversion of matrix T failed.', srname)
        ENDIF

!       Reserve G array for future use:
        CALL deallocate_array(w)
        CALL deallocate_array(q)
        CALL deallocate_array(r)
        CALL deallocate_mtz(my_mtz)
        CALL deallocate_array(my_maps)

!       Can do some symmetry stats and improvement here:
!       FIXME
!        debug = 101
!        IF ( debug > 100 ) THEN
            WRITE(*,*)
            WRITE(*,*) ' matrix O:'
            DO i = 1, MIN(10,np) 
                WRITE(*,"(1000ES10.2)") (O(i,j),j=1,MIN(10,np))
            ENDDO
!            STOP 'STOP ABCD' 
!        ENDIF

!        WRITE(*,"(' QR_MATINV> ', 'lambda (max)= ',  ES9.2, ' lambda (min)=', ES9.2,  ' ratio (condition number)=', ES9.2)")&
!        e(np), e(1), e(np)/e(1)

!       Accumulate trace of T**(-1):
         DO  i = 1, np
             x = O(i,1:np)

!            Linear solver, tridiagonal matrix has been factorized by PTTRF:
             CALL pttrs(alpha(1:np), beta(1:np-1), x(1:np), info=info )
             g(i) = DOT_PRODUCT ( O(i,1:np), x(1:np) )
         ENDDO
        
!       No longer need matrix A:
        IF ( ALLOCATED ( O ) ) CALL deallocate_array(O)
        IF ( ALLOCATED ( x ) ) CALL deallocate_array(x)
        IF ( ALLOCATED ( alpha ) ) CALL deallocate_array(alpha)
        IF ( ALLOCATED ( beta)  ) CALL deallocate_array(beta)        

!       Print result:
        CALL messag('Matrix diagonal has been calculated.', srname)
        WRITE(*,*) '  Parameter        Trace of inverse    Standard uncertainty'

        DO i = 1, np
            WRITE(*,*) i, g(i), SQRT ( f0 * g(i) / ( nobs - np ) )
        ENDDO
        WRITE(*,*) ' SUM of squared sigmas=', SUM ( g ) * f0 / ( nobs - np ) 

        trace_xyz = 0.0_wp
        l = 0
        DO i = 1, SIZE ( my_pdb )
            DO k = 1, 3            
                IF ( my_pdb%refinable(k,i) ) THEN
                    l = l + 1
                    trace_xyz = trace_xyz + g(l)
                ENDIF
            ENDDO
            l = l + COUNT ( my_pdb%refinable(4:11,i) )
        ENDDO

!       Print xyz statistics:
        WRITE(*,*) ' SUM of squared sigmas for xyz=', trace_xyz * f0 / ( nobs - np )

!       Dividing by the number of atoms with non-zero occ:
        WRITE(*,*) ' average delta r error=',  SQRT ( trace_xyz * f0 / ( nobs - np ) / COUNT ( my_pdb%occ > 0.0_wp ) )

!       Massive deallocation:

!       Gradient:
        IF ( ALLOCATED ( g ) ) CALL deallocate_array(g)
        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       MY_PDB will be deallocated in MAIN program:
        
        WRITE ( *, "(' LANCZOS_MATINV> Total CPU time for the whole run=       ', F10.1,' s')") time1 - time0
        WRITE ( *, "(' LANCZOS_MATINV> Total elapsed time for matrix inversion=', F10.2,' s')") SECNDS(t0)
        WRITE ( *, "(' LANCZOS_MATINV> Total orthogonalization time=',            F10.2,' s')") &
        orthogonalization_time
        WRITE ( *, "(' LANCZOS_MATINV> Normal termination.')")

    END SUBROUTINE lanczos_matinv

    SUBROUTINE bekas_diagonal_estimate ( my_cards, my_pdb )
        IMPLICIT NONE
!
!       INTEGER, PARAMETER :: lwp = SELECTED_REAL_KIND(p=23)
!
!       Purpose:
!       =======
!       Simple matrix diagonalization.
!
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Total number of observations:
        INTEGER                                    :: nobs         
!       Total number of parameters:                   
        INTEGER                                    :: np
!       Type of refinement:
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       List:
        INTEGER,       DIMENSION(:,:), ALLOCATABLE :: atom_pair
        INTEGER(KIND=eb)                           :: nnz
        TYPE(anymatrix), DIMENSION(:), ALLOCATABLE :: ANYMAT       
        TYPE(csr_matrix)                           :: csr_2
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: v
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: av
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: r1
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: sk 
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: y 

        REAL(KIND=wp)                              :: trace_xyz 
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       SYMMLQ:
        INTEGER, PARAMETER                         :: nout = 6
        INTEGER                                    :: itnlim
        LOGICAL                                    :: checkA
!       Check symmetry of preconditioner:
        LOGICAL                                    :: goodb
        LOGICAL                                    :: precon
!       Slightly different preconditioning scheme:
        REAL(KIND=wp)                              :: shift
        REAL(KIND=wp)                              :: rtol
!       LAPACK:
        REAL(KIND=wp), EXTERNAL                    :: dnrm2
!       Bekas parameters:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: diag_estimate
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: old_diag_estimate
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: qk
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: tk
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: H ! Hadamard matrix
        LOGICAL                                    :: hadamard_vectors
        INTEGER                                    :: m
        INTEGER                                    :: n
        INTEGER                                    :: non_positive
!       CPU vars:
        REAL(KIND=sp)                              :: t0
        REAL(KIND=sp)                              :: time0
        REAL(KIND=sp)                              :: time1
        REAL(KIND=sp)                              :: symmlq_time0
        REAL(KIND=wp)                              :: orthogonalization_time
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: j 
        INTEGER                                    :: k
        INTEGER                                    :: l
        INTEGER                                    :: ntrials
        CHARACTER(LEN=80), SAVE                    :: srname = 'bekas_diagonal_estimate'

!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us)
        polar_axes = my_mtz%sp_group%lpaxis

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 33 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st group along Y axis...', srname)
        ENDIF

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' LANCZOS_MATINV> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )
        IF ( np < 2 ) THEN
            WRITE(*,*) ' np=', np
            CALL die('Too small a number of refinable parameters.', srname)
        ENDIF

!       Check whether we can do LSQ refinement without restraints:
        nobs = SIZE ( my_mtz%fo )
        WRITE(*,"(' QR_MATINV> ', 'Number of reflections=', I6, ' number of params=', I8)") &
        nobs, np
        IF ( nobs < np + 1 ) THEN
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' LANCZOS_MATINV> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' LANCZOS_MATINV> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' LANCZOS_MATINV> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) )

        IF ( .NOT. my_cards%sigma_weighting ) THEN
            WRITE(*,*) ' B for weighting=', my_cards%bw
            my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
            my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )
        ELSE
             IF ( ANY ( my_mtz%sigfo == 0.0_wp ) ) THEN
                 CALL die('Some sigmas are equal to zero.', srname)
             ENDIF
             my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
             my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' QR_MATINV> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       SYMMLQ array allocations:
        CALL allocate_array(g,  np)  ! 
        CALL allocate_array(v,  np)  ! random vector?
        CALL allocate_array(av,  np)  ! random vector?
        CALL allocate_array(sk,  np) 
        CALL allocate_array(y,  np) 
        CALL allocate_array(r1, np)
      
        CALL allocate_array (diag_estimate, np)
        CALL allocate_array (old_diag_estimate, np)
        CALL allocate_array (tk, np)
        CALL allocate_array (qk, np)

        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0)

!       Having b_min it is necessary to adjust b_scale = badd - bmin:
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       This call has to go first before matrix calcns since we need phases, scale and (perhaps) other params:
        CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                           my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 
       
!       CALL dv(matrix, g, my_mtz, my_pdb)

!       Initialize vectors:
        diag_estimate = 0.0_wp; tk = 0.0_wp; qk = 0.0_wp
        orthogonalization_time = 0 

!       Calculate very thick preconditioner (note that we need to do this only once):
        CALL residue_block_diagonal_pairs ( atom_pair, my_pdb, residues_in_block=5, nout = 6)
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       Create template for sparse matrix:
        nnz = calc_nnz ( my_pdb, atom_pair )
        CALL allocate_matrix ( my_pdb%coo, np, np, nnz )
        CALL calculate_sparse_matrix(my_pdb%COO, my_mtz, my_maps, my_pdb, atom_pair, sk)
!       End of preconditioner calculations --

        CALL copy_coo_matrix_to_array_of_anymatrices ( my_pdb%COO, my_pdb, ANYMAT )
        CALL copy_array_of_inverted_anymatrices_to_coo_matrix(my_pdb%COO, my_pdb, ANYMAT, xpower=-1.0_wp)

        my_pdb%csr = my_pdb%coo
!       Need to allocate additional CSR matrix to keep Q**(1/2):
        CALL copy_array_of_inverted_anymatrices_to_coo_matrix(my_pdb%COO, my_pdb, ANYMAT, xpower=+0.5_wp)
        CSR_2 = my_pdb%coo

!       Extract diagonal:
        j = 0
        sk = 0
        DO i = 1, my_pdb%coo%nnz
            IF ( my_pdb%coo%ia(i) == my_pdb%coo%ja(i) ) THEN
                j = j + 1
                IF ( j > np ) THEN
                    CALL die('Programming error. Too many diagonal elements.', srname)      
                ENDIF
                sk(j) = my_pdb%coo%val(i)
            ENDIF
        ENDDO
        sk = 1.0_wp / SQRT ( sk )

!       No longer need that:
        CALL deallocate_matrix(my_pdb%coo)
        CALL deallocate_array(ANYMAT)


!       Initialize random number generator:
        CALL init_random_number()

!       Intialize SYMMLQ parameters:
        checkA = .TRUE.
        goodb  = .FALSE.
        precon = .TRUE.
        shift  = 0.0_wp
        itnlim = np/4
        rtol = EPSILON ( 1.0_wp ) ** 0.33333_wp

        hadamard_vectors = .TRUE.
        IF ( hadamard_vectors ) THEN
            m = 256 * 2 
!           Figure out n:
            n = 1
            DO k = 1, 30
                n = 2 * n

!               Found suitable n:
                IF ( n > np ) EXIT
            ENDDO

            WRITE(*,*) k, n, ' for Hadamard matrix'
            CALL allocate_array ( H, m, n)
            CALL hadamard(m, n, H)
            WRITE(*,*) ' Hadamard matrix has been generated. Occupying ', m*n*8/1024**2, ' MB.'
        ENDIF


        CALL CPU_TIME ( symmlq_time0 )
        t0 = SECNDS(0.0)
        old_diag_estimate = -1

        IF ( .NOT. hadamard_vectors ) m = np / 4 

!????
        m = np !etc

!       Calculate the inverse matrix estimate using SYMMLQ and MC Method (this is the most lengthy part of calculations):
        DO ntrials = 1, m

            WRITE(*,"(//,1X,80('-'))")
            IF ( .NOT. hadamard_vectors ) THEN
                CALL RANDOM_NUMBER(v)
            
!               Create random number distribution between -1 and +1:
                v = 2.0_wp * v - 1.0_wp
            ELSE
                v(1:np) =  H(ntrials,1:np)
            ENDIF

!           Scaling is VERY important:
!            v(1:np) =  sk(1:np) * v(1:np)

!            CALL matvec(np, r1, v, CSR_2)
!            v = r1

!           Implicit inverted matrix multiplication:
            CALL Ax(np, r1, v,  my_mtz, my_maps, my_pdb, my_cards%refinement_mode )
!            CALL Ax(np, av, r1, my_mtz, my_maps, my_pdb, my_cards%refinement_mode )

!            CALL symmlq (np, Ax, RESMIN_INV, v, shift, checkA, precon,                       &
!                         av, itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm,     &
!                         my_mtz, my_maps, my_pdb, my_cards%refinement_mode, sk, my_pdb%csr)

!           Recheck solution:
!            CALL Ax(np, y, av, my_mtz, my_maps, my_pdb, my_cards%refinement_mode)

!            r1   = v - y
!            r1norm = dnrm2 ( np, r1, 1 )

!            IF ( debug > 25 ) THEN
!                DO i = 1, MIN ( 10, np )
!                    WRITE(*,"(' b=', ES12.5, ' y=', ES12.5, ' diff= ', ES9.2)") &
!                              v(i), y(i), r1(i)
!                ENDDO
!            ENDIF

!           Print residual:
!            WRITE(nout, 2000) r1norm
!       2000 FORMAT(/ ' Final residual =', 1p, e8.1)
!            IF ( debug > 30 ) THEN
!                WRITE(nout, 2100) (i, av(i), i = 1, np)
!       2100     FORMAT(/ ' Solution  x' / 1p, 4(i6, e14.6))
!            ENDIF
                                  

!            CALL matvec(np, r1, av, CSR_2)
!            av = r1

            tk = tk + av * v
            qk = qk + v * v
            diag_estimate = tk / qk

            WRITE(*,*) ' Parameter  Inverse Diagonal  SU Estimate'
            DO i = 1, MIN(10,np)
                WRITE(*,"( I8,4X,ES11.4,2X,ES10.3)") i, diag_estimate(i), SQRT ( f0 * diag_estimate(i) / ( nobs - np ) )
            ENDDO

!           Check how many non-positive elements we have now:
            non_positive = COUNT ( diag_estimate < 0.0_wp )
            IF ( ALL ( old_diag_estimate > 0.0_wp ) .AND. ALL ( diag_estimate > 0.0_wp ) ) THEN
                old_diag_estimate = ABS ( old_diag_estimate - diag_estimate ) / diag_estimate
                WRITE(*,*) ' Max change:', MAXVAL ( old_diag_estimate )
                IF ( MAXVAL ( old_diag_estimate ) < 0.005_wp ) THEN
                    WRITE(*,*) ' Diagonal stabilized. Leaving the loop.'
                    EXIT
                ENDIF
            ENDIF
            WRITE(*,"(' RANDOM_DIAG_ESTIMATE> ', 'Number of solved linear systems=', I6, &
           &          ' trace=', ES11.4, ' non-positive=',I4)") &
            ntrials, SUM ( diag_estimate ), non_positive

            old_diag_estimate = diag_estimate
            

        ENDDO

        CALL CPU_TIME(time1)
        WRITE(*,"(' BEKAS_DIAG_ESTIMATE> ', ' CPU time for matrix calcns=', F10.1, ' s Elapsed time=', F10.2, ' s')") &
        time1 - time0, SECNDS(t0)
        t0 = SECNDS(0.0)


!       Reserve G array for future use:
        g = diag_estimate
        CALL deallocate_array(diag_estimate)
        CALL deallocate_array(qk)
        CALL deallocate_array(tk)
        CALL deallocate_mtz(my_mtz)
        CALL deallocate_array(my_maps)

        IF ( ALLOCATED ( v ) ) CALL deallocate_array(v)
        IF ( ALLOCATED ( av ) ) CALL deallocate_array(av)

!       Print result:
        CALL messag('Matrix diagonal has been calculated.', srname)
        WRITE(*,*) '  Parameter        Trace of inverse    Standard uncertainty'

        DO i = 1, np
            WRITE(*,*) i, g(i), SQRT ( f0 * g(i) / ( nobs - np ) )
        ENDDO
        WRITE(*,*) ' SUM of squared sigmas=', SUM ( g ) * f0 / ( nobs - np ) 

        trace_xyz = 0.0_wp
        l = 0
        DO i = 1, SIZE ( my_pdb )
            DO k = 1, 3            
                IF ( my_pdb%refinable(k,i) ) THEN
                    l = l + 1
                    trace_xyz = trace_xyz + g(l)
                ENDIF
            ENDDO
            l = l + COUNT ( my_pdb%refinable(4:11,i) )
        ENDDO

!       Print xyz statistics:
        WRITE(*,*) ' SUM of squared sigmas for xyz=', trace_xyz * f0 / ( nobs - np )

!       Dividing by the number of atoms with non-zero occ:
        WRITE(*,*) ' average delta r error=',  SQRT ( trace_xyz * f0 / ( nobs - np ) / COUNT ( my_pdb%occ > 0.0_wp ) )

!       Massive deallocation:

!       Gradient:
        IF ( ALLOCATED ( g ) ) CALL deallocate_array(g)
        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       MY_PDB will be deallocated in MAIN program:
        
        WRITE ( *, "(' BEKAS_DIAGONAL_ESTIMATE> Total CPU time for the whole run=       ', F10.1,' s')") time1 - time0
        WRITE ( *, "(' BEKAS_DIAGONAL_ESTIMATE> Total elapsed time for matrix inversion=', F10.2,' s')") SECNDS(t0)
        WRITE ( *, "(' BEKAS_DIAGONAL_ESTIMATE> Total orthogonalization time=',            F10.2,' s')") &
        orthogonalization_time
        WRITE ( *, "(' BEKAS_DIAGONAL_ESTIMATE> Normal termination.')")

    END SUBROUTINE bekas_diagonal_estimate

    SUBROUTINE estimate_standard_uncertainties ( my_cards, my_pdb )
        IMPLICIT NONE
!
!       INTEGER, PARAMETER :: lwp = SELECTED_REAL_KIND(p=23)
!
!       Purpose:
!       =======
!       Simple matrix diagonalization.
!
        TYPE(cards),                 INTENT(INOUT) :: my_cards
!       PDB:
        TYPE(pdb),                   INTENT(INOUT) :: my_pdb
!       Local variables:
        TYPE(mtz)                                  :: my_mtz
        TYPE(map),     DIMENSION(:),   ALLOCATABLE :: my_maps
!       Bulk:
        TYPE(map)                                  :: binary_map
        TYPE(map_logical)                          :: logical_map
!     
        LOGICAL,       DIMENSION(3)                :: polar_axes
!       Test:
        REAL(KIND=wp), DIMENSION(3)                :: xyz
!       Total number of observations:
        INTEGER                                    :: nobs         
!       Total number of parameters:                   
        INTEGER                                    :: np
!       Type of refinement:
        LOGICAL                                    :: refine_xyz
        LOGICAL                                    :: refine_biso
        LOGICAL                                    :: refine_occ
        LOGICAL                                    :: refine_Us
!       Maps:
        INTEGER                                    :: number_of_maps
        INTEGER                                    :: my_first_allocated_map ! USUALLY 1
!       List:
        INTEGER,       DIMENSION(:,:), ALLOCATABLE :: atom_pair
        INTEGER(KIND=eb)                           :: nnz
        TYPE(anymatrix), DIMENSION(:), ALLOCATABLE :: ANYMAT       
        TYPE(csr_matrix)                           :: Q
        TYPE(csr_matrix)                           :: QINV
!       Target function:
        REAL(KIND=wp)                              :: f0
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: g
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: v
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: av
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: sk 
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: y 
!       Lambdas for Robibson-Wathen method:
        REAL(KIND=wp)                              :: alpha
        REAL(KIND=wp)                              :: beta
        REAL(KIND=wp)                              :: sii
        REAL(KIND=wp), DIMENSION(2)                :: extreme_lambda
        REAL(KIND=wp), DIMENSION(2)                :: bounds 
        REAL(KIND=wp)                              :: dummy
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: ub
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: lb
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: diag_estimate
        REAL(KIND=wp)                              :: unorm2
!       Probably unneccesary but seems to work:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: O
        LOGICAL                                    :: orthog
!       Trace:
        REAL(KIND=wp)                              :: trace_xyz 
!       Anti-alias reduction:
        REAL(KIND=wp)                              :: bmin
!       LAPACK:
!        REAL(KIND=wp), EXTERNAL                    :: dnrm2
!       CPU vars:
        REAL(KIND=sp)                              :: t0
        REAL(KIND=sp)                              :: t1
        REAL(KIND=sp)                              :: time0
        REAL(KIND=sp)                              :: time1
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: j 
        INTEGER                                    :: k
        INTEGER                                    :: l
        INTEGER                                    :: mult 
        CHARACTER(LEN=80), SAVE                    :: srname = 'estimate_standard_uncertainties'

!       Which type of refinement we are dealing with:
        CALL which_params(my_cards%refinement_mode, refine_xyz, refine_biso, refine_occ, refine_Us)
        polar_axes = my_mtz%sp_group%lpaxis

!       Read mtz and allocate appropriate arrays:
        CALL read_mtz ( my_mtz, my_cards )

!       Save definitions of polar axes:
        polar_axes = my_mtz%sp_group%lpaxis
        CALL matrix_pointers(my_pdb, refine_xyz, refine_biso, refine_occ, refine_Us, polar_axes)

!       Fix appropriate coordinates in case of polar axis:
        IF ( COUNT ( my_pdb%refinable ) > 33 ) THEN
            IF ( refine_xyz .AND. ANY ( my_pdb%refinable(1:3,1) ) ) THEN
                my_pdb%refinable(1:3,1) = .NOT. polar_axes
            ENDIF
        ENDIF

        IF ( polar_axes(2) ) THEN
            CALL warn('No refinement for the 1st group along Y axis...', srname)
        ENDIF

!       FIXME
!       Maybe we don't need this anymore:
!       Reset definition of non_polar_axes to calculate proper 3D Fouriers:
        IF ( ANY ( my_mtz%sp_group%lpaxis ) ) THEN
            CALL messag ( 'Introducing changes to polar axes definitions...', srname )
            CALL reset_polar_axes ( my_mtz )
        ENDIF

!       Expand fo to hemisphere - ALL calculation will be done in P1 space group general routines:
        CALL messag('Expanding mtz file to hemisphere since symmetries are still unknown...',&
                    srname)

        CALL reduce_fo_to_hemisphere ( my_mtz ) 

!       Copy data to MY_PDB:
        my_pdb%cell  = my_mtz%cell
        my_pdb%ORT   = my_mtz%ORT
        my_pdb%DEORT = my_mtz%DEORT

!       Check first 3 coordinates:
        IF ( debug > 3 ) THEN
            xyz = my_pdb%xyz(1)
            WRITE(*,"(' LANCZOS_MATINV> ', &
                      'Checking first atom co-ordinate (should be orthogonal): xyz(1)= ', 3F7.2)") xyz
        ENDIF

!       Determine number of parameters:
        np = COUNT ( my_pdb%refinable )
        IF ( np < 2 ) THEN
            WRITE(*,*) ' np=', np
            CALL die('Too small a number of refinable parameters.', srname)
        ENDIF

!       Check whether we can do LSQ refinement without restraints:
        nobs = SIZE ( my_mtz%fo )
        WRITE(*,"(' QR_MATINV> ', 'Number of reflections=', I6, ' number of params=', I8)") &
        nobs, np
        IF ( nobs < np + 1 ) THEN
            CALL die ('Too low number of reflections. Please increase resolution. Refinement halted...', srname)
        ENDIF

        WRITE(*,"(' LANCZOS_MATINV> ', 'Refinement mode               = ', A)")&
        my_cards%refinement_mode
        WRITE(*,"(' LANCZOS_MATINV> ', 'Number of refinable parameters= ', A)")&
        TRIM ( int_to_c ( np ) )

!       Choose maps and prepare FFTW plans:
        CALL allocate_array ( my_maps, my_mtz, my_cards%refinement_mode )
        number_of_maps = SIZE ( my_maps )

!       Find first map which is really allocated:
        my_first_allocated_map = first_allocated_map ( my_maps )

!       See that number of maps equals 5 or 11:
        WRITE(*,"(' LANCZOS_MATINV> ', 'number of maps=', I3)")&
        number_of_maps 

!       Allocate weight_in_P1 array and set simple weighting scheme:
        CALL allocate_array ( my_mtz%weight_in_P1, SIZE ( my_mtz%fc_in_P1 ) )
        CALL allocate_array ( my_mtz%weight,       SIZE ( my_mtz%fo ) )

        IF ( .NOT. my_cards%sigma_weighting ) THEN
            WRITE(*,*) ' B for weighting=', my_cards%bw
            my_mtz%weight = EXP ( -0.25_wp * my_cards%bw * my_mtz%s )
            my_mtz%weight_in_P1 = EXP ( -0.25_wp * my_cards%bw * my_mtz%s_in_P1 )
        ELSE
             IF ( ANY ( my_mtz%sigfo == 0.0_wp ) ) THEN
                 CALL die('Some sigmas are equal to zero.', srname)
             ENDIF
             my_mtz%weight = 1.0_wp / my_mtz%sigfo ** 2
             my_mtz%weight_in_P1 = 1.0_wp / my_mtz%sigfo_in_P1 ** 2
        ENDIF

!       ************************************************
!       Allocate special fcQ arrays in MY_MTZ structure:
!       ************************************************

!       Need to think hard how to reduce their number (currently some may be allocated but unused):
        CALL allocate_array ( my_mtz%fcq_in_P1, SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

!       New array added for convolution speed up, NOV 2007 BVS:
        CALL allocate_array ( my_mtz%fcq_sym,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )
        CALL allocate_array ( my_mtz%fc_temp,   SIZE ( my_mtz%fc_in_P1 ), number_of_maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' QR_MATINV> ', 'number of FCQs= ', A)") &
            TRIM ( int_to_c ( SIZE ( my_mtz%fcq_in_P1, DIM=2 ) ) )
        ENDIF

        CALL prepare_bulk_solvent_maps(my_mtz, binary_map, logical_map)

!       Array allocations:
        CALL allocate_array(g,  np) 
        CALL allocate_array(v,  np)  
        CALL allocate_array(y,  np)  
        CALL allocate_array(av, np) 
        CALL allocate_array(sk, np)

        CALL allocate_array(diag_estimate, np)      
        CALL allocate_array(ub, np)
        CALL allocate_array(lb, np)

!       Special re-orthogonalization array:
        orthog = .FALSE.
        IF ( orthog ) THEN
            CALL allocate_array(O, np, 100)
        ENDIF

        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0)

!       Having b_min it is necessary to adjust b_scale = badd - bmin:
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       This call has to go first before matrix calcns since we need phases, scale and (perhaps) other params:
        CALL prepare_atom_images(my_maps(my_first_allocated_map), my_pdb, my_pdb%all_atom_images)
        CALL new_lsq_first_deriv_wrt_all ( g, my_mtz, my_maps(my_first_allocated_map), binary_map, logical_map, &
                                           my_pdb, mode=my_cards%refinement_mode, fun=f0 ) 
       
!       CALL dv(matrix, g, my_mtz, my_pdb)

!       Initialize vectors:
        ub = 0.0_wp; lb = 0.0_wp

!       Calculate very thick preconditioner (note that we need to do this only once):
        CALL residue_block_diagonal_pairs ( atom_pair, my_pdb, residues_in_block=20, nout = 6)
        bmin = get_bmin ( my_pdb )
        CALL apply_bmin_to_array_of_maps(my_maps, bmin)

!       Adjust bmin in my_mtz as well:
        my_mtz%b_scale = my_mtz%badd - bmin

!       Create template for block diagonal matrix:
        nnz = calc_nnz ( my_pdb, atom_pair )
        CALL allocate_matrix ( my_pdb%coo, np, np, nnz )

!       Calculate preconditioner Q:
        CALL calculate_sparse_matrix(my_pdb%COO, my_mtz, my_maps, my_pdb, atom_pair, sk)

!       Extract diagonal immediately:
        j = 0
        DO i = 1, my_pdb%coo%nnz

!           Note that matrix maybe unexpanded and contain some zero indices:
            IF ( my_pdb%coo%ia(i) == my_pdb%coo%ja(i) .AND. my_pdb%coo%ja(i) /= 0 ) THEN
                j = j + 1
                IF ( j > np ) THEN
                    WRITE(*,*)  my_pdb%coo%ia(i), my_pdb%coo%ja(i)                    
                    WRITE(*,*) j, np                    
                    CALL die('Programming error. Too many diagonal elements.', srname)
                ENDIF
                sk(j) = my_pdb%coo%val(i)
            ENDIF
        ENDDO
!       Check the length of main diagonal:
        IF ( j /= np ) THEN
            WRITE(*,*) j, np
            CALL die('Programming error. Incorrect length of main diagonal', srname)         
        ENDIF
        g = sk
        sk = 1.0_wp / SQRT ( sk )

!       Get ANYMAT representation:
        CALL copy_coo_matrix_to_array_of_anymatrices ( my_pdb%COO, my_pdb, ANYMAT )

!       Need to allocate additional CSR matrix to keep QINV=Q**(-1/2):
        CALL copy_array_of_inverted_anymatrices_to_coo_matrix(my_pdb%COO, my_pdb, ANYMAT, xpower=-0.5_wp)
        QINV = my_pdb%COO

!       Need to allocate additional CSR matrix to keep Q**(1/2) which will be used rarely:
        CALL copy_array_of_inverted_anymatrices_to_coo_matrix(my_pdb%COO, my_pdb, ANYMAT, xpower=+0.5_wp)
        CALL deallocate_array(ANYMAT)
        Q = my_pdb%COO
        CALL deallocate_matrix(my_pdb%COO)

        CALL CPU_TIME ( time0 )
        t0 = SECNDS(0.0)

!       Calculcate extreme eigenvalues for Q**(-1/2) * H * Q**(-1/2):
!       Start with smallest lambda:
        CALL implicit_shift_invert(np, 1, 21, 'LM', my_mtz, my_maps, my_pdb, my_cards%refinement_mode, extreme_lambda(1), &
                                   dummy, sk, ACSR=QINV)

!       Largest lambda: CHECK??
        CALL implicit_dsdrv1(np, 1, 5, 'LM', my_mtz, my_maps, my_pdb, my_cards%refinement_mode, dummy, &
                             extreme_lambda(2), sk, Q=QINV)
!        CALL implicit_dsdrv1(np, 2, 100, 'BE', my_mtz, my_maps, my_pdb, my_cards%refinement_mode, extreme_lambda(1), &
!                             extreme_lambda(2), sk, Q=QINV)

        WRITE(*,"(' IMPLICIT_TRACE_ESTIMATION> ', 'Extreme lambdas (ARPACK est.)=', 2ES9.2)") &
        extreme_lambda(1:2)

        CALL sleep(3)

!       Save extreme values:
        alpha = extreme_lambda(1)
        beta  = extreme_lambda(2)

!       Calculate the Robinson-Wathen inverse matrix trace estimate:
        v = 0
        mult = 0
        DO i = 1, np
            t1 = SECNDS(0.0)
            WRITE(*,"(//)")    
            CALL messag(' Estimate value for quadratic form&
            & v('//TRIM(int_to_c(i))//')*f(A)*v('//TRIM(int_to_c(i))//'), where f(A)=A^-1', srname)
            v(i) = 1

!           Calculate Q**(-1/2) * v:
            CALL matvec(np, av, v, QINV)

!           Implicit matrix multiplication y = H * Q**(-1/2) * v (remember no SK):
            CALL Ax(np, y, av,  my_mtz, my_maps, my_pdb, my_cards%refinement_mode )

!           Calculate Q**(-1/2) * H * Q**(-1/2) * v:
            CALL matvec(np, av, y, QINV)

!           Famous Robinson-Wathen bounds for trace of inverse:
            sii = DOT_PRODUCT ( av, av )
            ub(i) = 1.0_wp / alpha + (alpha - av(i)) ** 2  / (alpha * ( alpha * av(i) - sii ) )
            lb(i) = 1.0_wp / beta  - (av(i) - beta)  ** 2  / (beta  * ( sii - beta * av(i) ) )

            extreme_lambda(1:2) = (/lb(i),ub(i)/)
!                                                            T                                             -1
!           Try to figure out product: Q**(1/2)*v to obtain v [Q**(1/2)*Q**(-1/2)* H * Q**(-1/2) * Q**(1/2)] * v
            CALL matvec(np, y, v, QINV)
            unorm2 =  DOT_PRODUCT ( y, y ) 
            WRITE(*,"(' ')")
            WRITE(*,"('                                     Frobenius              Robinson-Wathen     ')")
            WRITE(*,"(' Parameter  Diagonal element         row norm        Lower bound       Upper bound')")
            WRITE(*,"(I10,4(4X,ES14.6))") i, av(i), sii, lb(i)*unorm2, ub(i)*unorm2
            IF ( .NOT. orthog ) THEN
                CALL golub_implicit(np, my_mtz, my_maps, my_pdb, my_cards%refinement_mode, &
                                    extreme_lambda, y, bounds, nout=6, sk=sk, QINV=QINV, itn=k)
            ELSE
                CALL golub_implicit(np, my_mtz, my_maps, my_pdb, my_cards%refinement_mode, &
                                    extreme_lambda, y, bounds, nout=6, sk=sk, QINV=QINV, O=O, itn=k)
            ENDIF
            mult = mult + k
            WRITE(*,"(' Number of iterations=', I5, ' Total number of iterations=', I9, ' ESU=', F8.3)") &
                      k, mult,  SQRT ( bounds(2) * f0 / (nobs-np) )
            g(i) = MAXVAL ( bounds )
            WRITE(*,"(' Elapsed time to calculate this S.U.=', F10.2,' s')") SECNDS(t1)
!           Reset single-entry vector:            
            v(i) = 0
        ENDDO
        extreme_lambda(1:2) = (/alpha,beta/)

        CALL CPU_TIME(time1)
        WRITE(*,"(' ESTIMATE_STANDARD_UNCERTAINTIES> ', ' CPU time for inverse diagonal calcns=', F10.1, ' s Elapsed time=', F10.2, ' s')") &
        time1 - time0, SECNDS(t0)

!       Reserve G array for future use:
        CALL deallocate_array(diag_estimate)
        CALL deallocate_mtz(my_mtz)
        CALL deallocate_array(my_maps)

        IF ( ALLOCATED ( v ) ) CALL deallocate_array(v)
        IF ( ALLOCATED ( av ) ) CALL deallocate_array(av)

!       Print result:
        CALL messag('Matrix diagonal has been calculated.', srname)
        WRITE(*,*) '  Parameter        Trace of inverse    Standard uncertainty'

        DO i = 1, np
            WRITE(*,*) i, g(i), SQRT ( f0 * g(i) / ( nobs - np ) )
        ENDDO
        WRITE(*,*) ' SUM of squared sigmas=', SUM ( g ) * f0 / ( nobs - np ) 

        trace_xyz = 0.0_wp
        l = 0
        DO i = 1, SIZE ( my_pdb )
            DO k = 1, 3            
                IF ( my_pdb%refinable(k,i) ) THEN
                    l = l + 1
                    trace_xyz = trace_xyz + g(l)
                ENDIF
            ENDDO
            l = l + COUNT ( my_pdb%refinable(4:11,i) )
        ENDDO

!       Print xyz statistics:
        WRITE(*,*) ' SUM of squared sigmas for xyz=', trace_xyz * f0 / ( nobs - np )

!       Dividing by the number of atoms with non-zero occ:
        WRITE(*,*) ' average delta r error=',  SQRT ( trace_xyz * f0 / ( nobs - np ) / COUNT ( my_pdb%occ > 0.0_wp ) )

!       Massive deallocation:
        CALL deallocate_matrix(QINV)
        CALL deallocate_matrix(Q)
   

!       Gradient:
        IF ( ALLOCATED ( g ) ) CALL deallocate_array(g)
        CALL destroy_bulk_solvent_maps(binary_map, logical_map)

!       MY_PDB will be deallocated in MAIN program:
        CALL CPU_TIME(time1)        
        WRITE ( *, "(' ESTIMATE_STANDARD_UNCERTAINTIES> Total CPU time for the whole run=       ', F10.1,' s')") time1 - time0
        WRITE ( *, "(' ESTIMATE_STANDARD_UNCERTAINTIES> Total elapsed time for matrix inversion=', F10.2,' s')") SECNDS(t0)
        WRITE ( *, "(' ESTIMATE_STANDARD_UNCERTAINTIES> Normal termination.')")

    END SUBROUTINE estimate_standard_uncertainties

END MODULE lsq

   FUNCTION mkl_progress( thread, step, stage )
       INTEGER                 :: mkl_progress
       INTEGER                 :: thread, step
       CHARACTER(LEN=*)        :: stage
!       CHARACTER(LEN=32), SAVE :: srname='mkl_progress'
       WRITE(*,"(' MKL_PROGRESS> ', 'thread=', I3, ' stage= ', A, ' step=', I8)") thread, stage, step
!       PRINT*,'Thread:',thread,',stage:',stage,',step:',step
       mkl_progress = 0
       RETURN
    END FUNCTION mkl_progress

