MODULE fast_genden
USE aniso_manip
USE aniso_symmetry_manip
USE atom_images
USE basic_pdb_manip
USE map_fft
USE logical_map_manip
USE refinement_util
USE vectors
USE OMP_LIB
IMPLICIT NONE
CONTAINS
    SUBROUTINE prepare_atom_images(map_1, pdb_2, all_atom_images, true_sp_group)
!
!       Purpose:
!       =======
!       Generates atomic images for dramatic speed up of real space operations
!       during implicit matrix-vector product calculations.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2008      B.Strokopytov           Original code based on GENERATE_DENSITY routine.
!                                             Changed JUT to VECTOR_INT type.
!                                             Introduced atom_images module.
!                                             Added RUN variable.
!                                             Added N   counter.
!                                             Added NGP counter (number of grid points).
!
!       N.B.: Usage of symmetry is not recommended option in the context of this routine.
!
!       
        TYPE(map),                                   INTENT(INOUT)        :: map_1
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        TYPE(atom_image), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: all_atom_images
        TYPE(space_group),                           INTENT(IN), OPTIONAL :: true_sp_group
!       Local variables:
        INTEGER                                                           :: natom
        INTEGER                                                           :: nsym
!       e.d. stuff:
        INTEGER                                                           :: ityp      ! atom type 
        INTEGER                                                           :: ngauss   
        REAL(KIND=wp), DIMENSION(5)                                       :: a
        REAL(KIND=wp), DIMENSION(5)                                       :: b 
        REAL(KIND=wp), DIMENSION(5)                                       :: be
        REAL(KIND=wp)                                                     :: b_scale
        REAL(KIND=wp)                                                     :: r_cut_off
        REAL(KIND=wp)                                                     :: r2_cut_off
        REAL(KIND=wp)                                                     :: r2
        TYPE(vector_int)                                                  :: uvw
        TYPE(vector)                                                      :: duvw
!       Atom/box parameters:
        TYPE(vector)                                                      :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                                                  :: closest_grid_point
        TYPE(vector_int)                                                  :: box_size
        TYPE(vector_int)                                                  :: jut
        INTEGER,       DIMENSION(3)                                       :: box_lower_limit
        INTEGER,       DIMENSION(3)                                       :: box_upper_limit
        REAL(KIND=wp), DIMENSION(3)                                       :: ort_diag
        TYPE(vector)                                                      :: deort_diag
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                                     :: r2_cell_max
!       Aniso:
        LOGICAL                                                           :: aniso_active
        LOGICAL                                                           :: Us_present
!       Need to calculate eigenvalues:
        TYPE(matrix)                                                      :: UMAT
!       Need U data for each gaussian:
        TYPE(aniso),   DIMENSION(5)                                       :: U
        TYPE(matrix),  DIMENSION(5)                                       :: VTS
        REAL(KIND=wp), DIMENSION(5)                                       :: udet
        REAL(KIND=wp), DIMENSION(5)                                       :: qform
        REAL(KIND=wp)                                                     :: Biso !Uiso
        REAL(KIND=wp), DIMENSION(3)                                       :: ueigen
        REAL(KIND=wp), DIMENSION(6)                                       :: v
        TYPE(vector)                                                      :: duvwort
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE                      :: RU 
!       Map array dimensions:
        INTEGER,       DIMENSION(3)                                       :: dims
!       Counters:
        INTEGER                                                           :: iat
        INTEGER                                                           :: l 
        INTEGER                                                           :: m
!       Number of grid point for single atom image:
        INTEGER                                                           :: n
        INTEGER                                                           :: run
!       Box counters:
        INTEGER                                                           :: ju
        INTEGER                                                           :: jv
        INTEGER                                                           :: jw
        CHARACTER(LEN=32),                                         SAVE   :: srname = 'prepare_atom_images'
!       Figure out how many density grid points is involved:
        INTEGER(KIND=eb)                                                  :: dgp
        
!       Determine map array bounds (important for 1-d reindexing):
        dims = UBOUND(map_1%array) - LBOUND(map_1%array) + 1

!       Make very basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname) 
        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' GENERATE_DENSITY> ', 'natom=', A)" ) TRIM ( int_to_c ( natom ) )
            CALL die('Zero number of atoms on input.', srname)
        ENDIF

        IF ( ALLOCATED ( pdb_2%atomsf_table ) ) THEN
            IF ( SIZE ( pdb_2%atomsf_table ) < 1 ) THEN
                WRITE(*,*) SIZE ( pdb_2%atomsf_table )
                CALL die('Programming error. ATOMSF_TABLE has not been initialized properly.',&
                         srname )
            ENDIF
        ELSE
            CALL die('ATOMSF_TABLE has not been initialized.', srname)
        ENDIF

        IF ( .NOT. ALLOCATED ( pdb_2%atom_type ) ) THEN
            CALL die('Programming error. ATOM_TYPE has not been initialized properly.',&
                     srname )
        ELSE
            IF ( COUNT ( pdb_2%atom_type > 0 ) /= natom ) THEN
                CALL die('Some atoms have unassigned ATOM_TYPE.', srname)
            ENDIF
        ENDIF

        aniso_active = ALLOCATED ( pdb_2%U )

!       Initialize map (simply for test):
        map_1%array = 0.0_wp

!       Initialize parameters:
        b_scale    = map_1%b_scale
        deort_diag = extract_diagonal ( map_1%DEORT )
        ort_diag   = extract_diagonal ( map_1%ORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1%cell
            CALL die ( 'Programming error. Unit cell has not been initialized for current map.',&
                       srname )
        ENDIF

!       Calculate maximum sphere radius around any atom:
        r2_cell_max = 0.25_wp * MINVAL ( ort_diag(1:3) ** 2 )
        
!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
            CALL allocate_array(RU, nsym)
!           Convert symmetry operators to 6x6 matrix:
            DO m = 1, nsym
!               This needs testing:
                RU(m) = .SIXBYSIX. ( map_1%ORT * (.SYMA. map_1%sp_group%sym(m) * map_1%DEORT ) )
            ENDDO

        ELSE
            nsym = 1
        ENDIF

!       Allocate array of atom images:
        CALL allocate_array(all_atom_images, natom)        

        runs:DO run = 0, 1

!           Estimate total number of grid points to store in memory:
            dgp = 0
            symmetry:DO m = 1, nsym
!$OMP PARALLEL DO DEFAULT (PRIVATE) SHARED ( RUN, M, NSYM, DIMS, NATOM, MAP_1, PDB_2, ALL_ATOM_IMAGES, &
!$OMP                                        B_SCALE, ANISO_ACTIVE, ORT_DIAG, DEORT_DIAG, R2_CELL_MAX, &
!$OMP                                        DGP, RU, DEBUG, SRNAME)
                coord:DO iat = 1, natom

!                   This implicitly assumes that we use standard GENERATE_DENSITY ROUTINE (BVS JUN 2009).
!                   The purpose of excluding atom with all unrefinable parameters is memory economy,
!                   nevertheless if we want to generate density for ALL atoms which is normally the case,
!                   we should use STANDARD ROUTINES and NOT images generated here:
                    IF ( pdb_2%occ(iat) == 0.0_wp .OR. .NOT. ANY ( pdb_2%refinable(1:11,iat) ) ) THEN 
!$OMP CRITICAL
                       IF ( run == 0 ) CALL allocate_atom_image(all_atom_images(iat), 0)
!$OMP END CRITICAL
                       CYCLE
                    ENDIF

!                   Save atom type:
                    ityp = pdb_2%atom_type(iat)

!                   Calculate number of gaussians in this entry of table:
                    ngauss = pdb_2%atomsf_table%record(ityp)%ngauss

                    IF ( ngauss /= SIZE ( pdb_2%atomsf_table%record(ityp) ) ) THEN
                        WRITE(*,*) ngauss, SIZE ( pdb_2%atomsf_table%record(ityp) ), ityp
                        CALL die('Programming error. Inconsistency in ATOMSF_TABLE records.', srname )
                    ENDIF

!                   Two possibilities here only:
                    IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                        WRITE(*,*) 'ngauss=', ngauss
                        CALL die('Programming error. Abnormal number of gaussians detected...',&
                                  srname )
                    ENDIF

                    a(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                    b(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%b(1:ngauss) 
                    Biso = pdb_2%biso(iat)
                    Us_present = .FALSE.

!                   Check whether Us are being used for this atom:
                    IF ( aniso_active ) THEN
                    
!                       Note that ">" operator was overloaded in a special way:
                        Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )

                        IF ( Us_present ) THEN

!$OMP CRITICAL
!                           Convert to matrix of aniso tensor elements:
                            UMAT = pdb_2%U(iat)

!                           Calculate eigenvalues:
                            ueigen = eigen_values ( UMAT )   ! excellent crash without CRITICAL
!$OMP END CRITICAL
                            Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )
                        ENDIF
                    ENDIF

!                   Finally we have sensible approach (using all gaussians to estimate radius):
!$OMP CRITICAL
                    r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + Biso + b_scale ) ** 2
                    r2_cut_off = MAX ( map_1%fft_radius ** 2 , r2_cut_off )

!                   Make sure sphere around atom fits into the cell:
                    r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                    r_cut_off  = SQRT( r2_cut_off )
!$OMP END CRITICAL
!                   Now we can continue transformations:
                    IF ( .NOT. Us_present ) THEN
                        a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                        be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )
                    ELSE
                        DO l = 1, ngauss

!$OMP CRITICAL
!                           Calculate array of anisotropic values for each gaussian:
                            U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )
!$OMP END CRITICAL
!                           Transform ellipsoid using Garib's method:
                            IF ( m /= 1 ) THEN
                                U(l) = RU(m) * U(l)
                            ENDIF

!                           Calculate determinant:
                            Udet(l) = udet_aniso( U(l) )

!                           Calculate inverted matrix:
                            v = minors_aniso(U(l)) / Udet(l)
                            VTS(l) = .CONVERT. v

                        ENDDO

!                       Scale constants in advance to reduce computation:
                        a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                    ENDIF

!                   Do NOT apply occupancy to modify atom image:
!                   a(1:ngauss) = pdb_2%occ(iat) * a(1:ngauss)

!                   Figure out box size around atom:
                    box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1%map_size ) )

!                   Convert to fractional:
!$OMP CRITICAL
                    IF ( m > 1 ) THEN
                        centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat) 
                    ELSE
                        centre_of_the_atom_in_fractions = map_1%DEORT        * pdb_2%xyz(iat)
                    ENDIF
!$OMP END CRITICAL
!                   This is our atomic centre on the grid:
                    closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1%map_size )
                    box_lower_limit    = closest_grid_point - box_size
                    box_upper_limit    = closest_grid_point + box_size

!                   Initialize atom grid points counter:
                    n = 0

!                   Create a box around the atom and fill it with atom e.d. values:
                    boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                        DO jv = box_lower_limit(2), box_upper_limit(2)
                            DO ju = box_lower_limit(1), box_upper_limit(1)
                                uvw  = (/ju, jv, jw/)
                                duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions   

!                               High point of electron density generation:
                                r2 = map_1%REAL_TENSOR * duvw
                                IF ( r2 <= r2_cut_off ) THEN

                                    n = n + 1

                                    IF ( run == 1 ) THEN

!                                       Fit uvw vector into map dimensions:
                                        jut = uvw .MMOD. map_1%map_size

!                                       Most CPU expensive part -- saving atomic image:
                                        IF ( Us_present ) THEN
                                            duvwort = map_1%ORT * duvw
                                            DO l = 1, ngauss
                                                qform(l) = duvwort .DOT. ( VTS(l) * duvwort )
                                            ENDDO
                                            all_atom_images(iat)%ro(n) = SUM (a(1:ngauss) * EXP ( -0.5_wp * qform(1:ngauss) ))
                                        ELSE
                                            all_atom_images(iat)%ro(n) = SUM ( a(1:ngauss) * EXP ( -be(1:ngauss) * r2 ) )
                                        ENDIF

!                                       1-D array of indices will start with 1 (conventional FORTRAN style),
!                                       note that we can't use map_1%nu variables since they reflect map dimensions
!                                       but map_1%array dimensions. See FFTW 3.0 manual.
                                        all_atom_images(iat)%uvw(n) = convert_to_1d ( jut, dims(1), dims(2), 1 )


                                        IF ( debug > 1019 ) THEN
                                            IF ( iat == 1 .AND. n <=10 ) THEN
                                                WRITE(*,"(' index(1d)=', I10, ' atomic density=', ES9.2, I4)") &
                                                          all_atom_images(iat)%uvw(n), all_atom_images(iat)%ro(n), map_1%nu
                                            ENDIF
                                        ENDIF

                                    ENDIF

                                ENDIF
                            ENDDO    
                        ENDDO
                    ENDDO boxz

!$OMP CRITICAL 
                    IF ( run == 0 ) THEN
                        CALL allocate_atom_image(all_atom_images(iat), n)
                    ENDIF
                    dgp = dgp + n
!$OMP END CRITICAL
                ENDDO coord

            ENDDO symmetry
        ENDDO runs

        WRITE(*,"(' PREPARE_ATOM_IMAGES> ', 'Number of grid points per atom: ', I12)") INT ( dgp / REAL(natom, KIND=wp) )
        WRITE(*,"(' PREPARE_ATOM_IMAGES> ', 'Number of grid points involved: ', I12)") dgp
        WRITE(*,"(' PREPARE_ATOM_IMAGES> ', 'Total memory needed (for speed up) :', I6, ' MB')") &
        NINT ( dgp * 8 * 1.5_wp / 2 ** 20 )

!       Free memory if necessary:
        IF ( ALLOCATED ( RU ) ) CALL deallocate_array ( RU )
    END SUBROUTINE prepare_atom_images 

    SUBROUTINE generate_density_using_images(map_1, pdb_2, all_atom_images)
!
!       Purpose:
!       =======
!       Merges atomic images forming a density map
!       This considerably may speed up calculations.
!       Simplicity of the code and almost total absence
!       of loops supports this idea.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2008       B.Strokopytov          Original code
!       Mar 2009       B.Strokopytov          This routine is bandwidth dependent
!                                             since there is almost nothing to compute.
!                                             OPENMP does not help to improve speed.
!
        TYPE(map),                      INTENT(INOUT) :: map_1
        TYPE(pdb),                      INTENT(IN)    :: pdb_2
        TYPE(atom_image), DIMENSION(:), INTENT(IN)    :: all_atom_images 
!       Local variables:
        INTEGER                                       :: natom
!       Counters:
        INTEGER                                       :: iat
        INTEGER                                       :: n 
        CHARACTER(LEN=32),                    SAVE    :: srname = 'generate_density_using_images'
        INTEGER, EXTERNAL                             :: OMP_GET_NUM_THREADS
        
!       Start with basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.', &
                     srname) 
        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' GENERATE_DENSITY> ', 'natom=', A)" ) TRIM ( int_to_c ( natom ) )
            CALL die('Zero number of atoms on input.', srname)
        ENDIF

!$OMP PARALLEL
!$OMP MASTER
      PRINT *,'----------------------------------------------'
!$    PRINT *,'Number of Threads = ',OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
!$OMP CRITICAL
      PRINT *,'----------------------------------------------'
!$OMP END CRITICAL
!$OMP PARALLEL
      PRINT *,'Printing one line per active thread....'
!$OMP END PARALLEL

!       Initialize map:
!        map_1%array = 0.0_wp
        CALL set_array_to_zero(map_1%array, SIZE ( map_1%array ) )

!OMP PARALLEL DO DEFAULT (PRIVATE) SHARED (pdb_2, all_atom_images, map_1, natom, srname)
        coord:DO iat = 1, natom

             IF ( pdb_2%occ(iat) == 0.0_wp ) CYCLE

             n = SIZE ( all_atom_images(iat)%ro )
             IF ( n > 0 ) THEN
                 CALL genden_no_1(n, all_atom_images(iat)%ro, all_atom_images(iat)%uvw,&
                                  pdb_2%occ(iat), map_1%array)
             ENDIF
        ENDDO coord
    END SUBROUTINE generate_density_using_images

    SUBROUTINE generate_density_for_map_array_using_images(map_1, pdb_2, q, all_atom_images, true_sp_group)
!
!       Purpose:
!       =======
!       Generates density for a few maps using 2- or 5-Gaussian approximation
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2007       B.Strokopytov          Original code based on GENERATE_DENSITY routine.
!                                             We are dealing with array of maps now.
!
        TYPE(map),     DIMENSION(:), ALLOCATABLE,     INTENT(INOUT)        :: map_1
        TYPE(pdb),                                    INTENT(IN)           :: pdb_2
        REAL(KIND=wp), DIMENSION(:),                  INTENT(IN)           :: q
        TYPE(atom_image), DIMENSION(:), ALLOCATABLE,  INTENT(IN)           :: all_atom_images
        TYPE(space_group),                            INTENT(IN), OPTIONAL :: true_sp_group
!       Local variables:
        INTEGER                                                            :: natom
        INTEGER                                                            :: nsym
!       Refinement:
        REAL(KIND=wp)                                                      :: occ
        INTEGER                                                            :: number_of_maps
        INTEGER                                                            :: my_first_allocated_map
!       Pointers:
        INTEGER                                                            :: indxyz
!       Counters:
        INTEGER                                                            :: iat
        INTEGER                                                            :: imap
        INTEGER                                                            :: l 
        INTEGER                                                            :: m
        INTEGER                                                            :: n 
!       OMP business:
        INTEGER                                                            :: max_threads

        CHARACTER(LEN=32)                                                  :: srname
        
        srname = 'generate_density_for_map_array_using_images'

!       Start with basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN

            natom = SIZE( pdb_2 )

        ELSE

            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname) 
        ENDIF

        IF ( natom < 1 ) THEN

            WRITE ( *, "(' GENERATE_DENSITY_FOR_MAP_ARRAY> ', 'natom=', A)" ) TRIM ( int_to_c ( natom ) )
            CALL die ( ' Zero or negative number of atoms.', srname )

        ENDIF

!       number_of_maps is either 5 or 11 ( when Us are being refined):
        number_of_maps = SIZE ( map_1 )
        IF ( number_of_maps /= 5 .AND. number_of_maps /= 11 ) THEN
            WRITE(*,*) ' number_of_maps=', number_of_maps
            CALL die('Programming error. Unusual number of maps...', srname)
        ENDIF

!       Long and hard check with respect to number of refinable parameters:
        IF ( COUNT ( pdb_2%refinable ) /= SIZE ( q ) ) THEN
            WRITE(*,*)  COUNT ( pdb_2%refinable ), SIZE ( q )
            CALL die('Programming error. Inconsistent size of Q (or REFINABLE) array...',&
                     srname)
        ENDIF

!       Figure out first allocated map in the array of MAPS:
        my_first_allocated_map = first_allocated_map ( map_1 )


!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( MY_FIRST_ALLOCATED_MAP, NUMBER_OF_MAPS, MAP_1)
        DO imap = my_first_allocated_map, number_of_maps 
            IF ( ALLOCATED ( map_1(imap) ) ) THEN
                map_1(imap)%array = 0.0_wp
            ENDIF
        ENDDO
!$OMP END PARALLEL DO

!       Check parameters:
        MAX_THREADS = 0
        DO imap = my_first_allocated_map, SIZE ( map_1 )
            IF ( ALLOCATED ( map_1(imap) ) ) THEN
                IF (  map_1(my_first_allocated_map)%b_scale /= map_1(imap)%b_scale )  THEN
                    WRITE(*,"(2(I2,F6.2))") my_first_allocated_map, map_1(my_first_allocated_map)%b_scale,&
                                            imap, map_1(imap)%b_scale
                   CALL die('Programming error. Unequal B_SCALE in different maps.',&
                            srname )
                ENDIF
!               Number of threads should not exceed number of allocated maps:
                MAX_THREADS = MAX_THREADS + 1

            ENDIF
        ENDDO

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1(my_first_allocated_map)%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1(my_first_allocated_map)%cell
            CALL die ( 'Programming error. Unit cell has not been initialized for current map.',&
                       srname )
        ENDIF

!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
        ELSE
            nsym = 1
        ENDIF

        symmetry:DO m = 1, nsym

            coord:DO iat = 1, natom

!               Do NOT touch this:
                l = COUNT ( pdb_2%refinable(1:11, iat) )

!               No refinable parameters -- skip this atom, e.g. hydrogen:
                IF ( l == 0 ) CYCLE

!               Save atom occupancy:
                occ = pdb_2%occ(iat)

!               Refinement of zero occncies is possible for occncy-only refinement:
                IF ( occ == 0.0_wp .AND. .NOT. pdb_2%refinable(5,iat) ) CYCLE

                indxyz = pdb_2%indxyz(iat)

!               Speed it up in case quasi-occupancies for this atom are ALL zero:
                IF ( ALL ( q(indxyz+1:indxyz+l) == 0.0_wp ) ) CYCLE ! really L not 11

!               Size of atom image:
                n = SIZE ( all_atom_images(iat)%ro )

!               atom_ed_value = SUM ( a(1:ngauss) * EXP ( -be(1:ngauss) * r2 ) )

!               That's perhaps the best point for parallelization since atoms will never intersect.
!               The drawback is that the threads will have to be initialized too often and number of
!               threads should be equal to number of active maps.

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (PDB_2, IAT, N, ALL_ATOM_IMAGES, Q, INDXYZ, MAP_1, OCC, NUMBER_OF_MAPS) 
                DO imap = 1, number_of_maps 
!                    WRITE(*,*) ' iat=', iat, n, indxyz
                    IF ( pdb_2%refinable(imap, iat) ) THEN

!                       This is good enough for parallel computing:
                        l = COUNT ( pdb_2%refinable(1:imap,iat) )
                        IF ( imap /= 5  ) THEN
                            CALL genden_no_2(n, all_atom_images(iat)%ro, all_atom_images(iat)%uvw, &
                                             q(indxyz+l) * occ, map_1(imap)%array)
                        ELSE
                            CALL genden_no_2(n, all_atom_images(iat)%ro, all_atom_images(iat)%uvw, &
                                             q(indxyz+l), map_1(imap)%array)
                        ENDIF

                    ENDIF
                ENDDO
!$OMP END PARALLEL DO

            ENDDO coord

        ENDDO symmetry

    END SUBROUTINE generate_density_for_map_array_using_images
 
END MODULE fast_genden 
