MODULE genden
USE aniso_manip
USE aniso_symmetry_manip
! New entry MAR 2010
USE basic_adjacency
USE basic_pdb_manip
USE map_fft
USE logical_map_manip
USE refinement_util
!USE vectors
IMPLICIT NONE
CONTAINS
    SUBROUTINE generate_density(map_1, pdb_2, true_sp_group)
!
!       Purpose:
!       =======
!       Generates density using 2- or 5-Gaussian approximation
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2003       B.Strokopytov          Original code
!
!       Apr 2005       B.Strokopytov          When calculating raidus cutoff occupancies shall be ignored
!                                             from now on to eliminate negative effects associated with
!                                             small occupancies. We need this for certain algorithms.
!
!       Nov 2005       B.Strokopytov          e_limit has been changes to 18.0 A in constants.f90
!                                             instead of original value of 7.0 which allows
!                                             in conjuction with appropriate anti-alias b_scale
!                                             extremely accurate calculation of structure factors
!                                             (error in phases of 0.001 degr or less ) at any resolution
!                                             using Shannon rate of just 1.5. See (Navaza, 2002, Acta Cryst.)
!                                             Maximum sphere radius is now checked because at low
!                                             resolution it may not even fit into cell
!
!       Oct 2007      B.Strokopytov           e_limit has been replaced with fft_radius now
!                                             present in both MTZ and MAP TYPE.
!
!       Oct 2007      B.Strokopytov           Major rearrangemnet of PDB type. Atomic scattering
!                                             factors are now in PDB_2%ATOMSF_TABLE
!
!       Oct 2007      B.Strokopytov           Improved calculations for max radii (see variable R2_CELL_MAX)
!       Mar 2009      B.Strokopytov           Trying to introduce ADPs
!       Aug 2009      B.Strokopytov           DEBUG was missing in the description of OMP loop. CORRECTED.
!       
        TYPE(map),                   INTENT(INOUT)        :: map_1
        TYPE(pdb),                   INTENT(IN)           :: pdb_2
        TYPE(space_group),           INTENT(IN), OPTIONAL :: true_sp_group
!       Local variables:
        INTEGER                                           :: natom
        INTEGER                                           :: nsym
!       e.d. stuff:
        INTEGER                                           :: ityp      ! atom type 
        INTEGER                                           :: ngauss   
        REAL(KIND=wp), DIMENSION(5)                       :: a
        REAL(KIND=wp), DIMENSION(5)                       :: b 
        REAL(KIND=wp), DIMENSION(5)                       :: be
        REAL(KIND=wp)                                     :: b_scale
        REAL(KIND=wp)                                     :: r_cut_off
        REAL(KIND=wp)                                     :: r2_cut_off
        REAL(KIND=wp)                                     :: r2
        TYPE(vector_int)                                  :: uvw
        TYPE(vector)                                      :: duvw
!       Atom/box parameters:
        TYPE(vector)                                      :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                                  :: closest_grid_point
        TYPE(vector_int)                                  :: box_size
        INTEGER,       DIMENSION(3)                       :: jut
        INTEGER,       DIMENSION(3)                       :: box_lower_limit
        INTEGER,       DIMENSION(3)                       :: box_upper_limit
        REAL(KIND=wp), DIMENSION(3)                       :: ort_diag
        TYPE(vector)                                      :: deort_diag
        REAL(KIND=wp)                                     :: atom_ed_value
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                     :: r2_cell_max
!       Aniso:
        LOGICAL                                           :: aniso_active
        LOGICAL                                           :: Us_present
!       Need to calculate eigenvalues:
        TYPE(matrix)                                      :: UMAT
!       Need U for each gaussian:
        TYPE(aniso),   DIMENSION(5)                       :: U
        TYPE(matrix),  DIMENSION(5)                       :: VTS
        REAL(KIND=wp), DIMENSION(5)                       :: udet
        REAL(KIND=wp), DIMENSION(5)                       :: qform
        REAL(KIND=wp)                                     :: Biso
        REAL(KIND=wp), DIMENSION(3)                       :: ueigen
        REAL(KIND=wp), DIMENSION(6)                       :: v
        TYPE(vector)                                      :: duvwort
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE      :: RU
!       Counters:
        INTEGER                                           :: iat
        INTEGER                                           :: l 
        INTEGER                                           :: m
!       Box counters:
        INTEGER                                           :: ju
        INTEGER                                           :: jv
        INTEGER                                           :: jw
        CHARACTER(LEN=32),                        SAVE    :: srname = 'generate_density'
!       Figure out how many density grid points is involved:
        INTEGER(KIND=eb)                                  :: dgp
!       CPU:
        REAL(KIND=sp)                                     :: t0
        REAL(KIND=sp)                                     :: time0
        REAL(KIND=sp)                                     :: time1
!        INTEGER, EXTERNAL                                 :: omp_get_num_threads
!        INTEGER, EXTERNAL                                 :: omp_get_thread_num
        INTEGER                                           :: num1
        INTEGER                                           :: num2
        INTEGER                                           :: nthreads
!        
        CALL CPU_TIME(time0)
        t0 = SECNDS(0.0_sp)

        dgp = 0

!       Start with basic checks:
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

!       Initialize map:
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

!       Check the number of available threads/processors:
        CALL get_num_threads(nthreads)

        symmetry:DO m = 1, nsym
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( M, NSYM, NATOM, MAP_1, PDB_2, B_SCALE, ANISO_ACTIVE, &
!OMP                                       ORT_DIAG, DEORT_DIAG, R2_CELL_MAX, DGP, RU, SRNAME, DEBUG)
            coord:DO iat = 1, natom

                IF ( pdb_2%occ(iat) == 0.0_wp ) CYCLE

!               Save atom type:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss = pdb_2%atomsf_table%record(ityp)%ngauss

                IF ( ngauss /= SIZE ( pdb_2%atomsf_table%record(ityp) ) ) THEN
                    WRITE(*,*) ngauss, SIZE ( pdb_2%atomsf_table%record(ityp) ), ityp
                    CALL die('Programming error. Inconsistency in ATOMSF_TABLE records.', srname )
                ENDIF

!               Two possibilities here only:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             srname )
                ENDIF

                a(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)
                Biso = pdb_2%biso(iat)
               
!               Check whether Us have been read (appropriate record present):
                Us_present = .FALSE.

                IF ( aniso_active ) THEN

!                   Note that ">" was overloaded in a special way:
                    Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )

                    IF ( Us_present ) THEN

!                       Convert to matrix of aniso elements:
                        UMAT = pdb_2%U(iat)

!                       Calculate eigenvalues:
!OMP CRITICAL
                        ueigen = eigen_values ( UMAT )
!OMP END CRITICAL
                        Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )
                    ENDIF
                ENDIF

                IF ( ANY ( b(1:ngauss) + Biso + b_scale > 6.0_wp * 10.0_wp ** 2 )  ) THEN
                    WRITE(*,*) ' Biso=', Biso, ' b_scale=', b_scale
                    WRITE(*,*) ' iat=', iat
                    IF ( Us_present ) THEN
                        WRITE(*,"(' eigenvalues=', 3F10.3)") ueigen
                        CALL warn('Unreasonable Uiso.', srname)
                    ENDIF
                ENDIF

!               Finally we have sensible approach (avoid problems with radius estimation during parallelization):
!OMP CRITICAL
                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + Biso + b_scale ) ** 2
                r2_cut_off = MAX ( map_1%fft_radius ** 2 , r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT( r2_cut_off )
!OMP END CRITICAL

!               Now we can continue transformations:
                IF ( .NOT. Us_present ) THEN

                    a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                    be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

                ELSE

                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian:
                        U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )

!                       Transform ellipsoid using Garib's method:
                        IF ( m /= 1 ) THEN
!                           This is not dangerous for parallel computations since U(l) is of type ANISO:
                            U(l) = RU(m) * U(l)
                        ENDIF

!                       Calculate determinant:
                        Udet(l)  = udet_aniso(U(l))

!                       Calculate inverted matrix:
                        v = minors_aniso(U(l)) / Udet(l)
                        VTS(l) = .CONVERT. v
                        
                    ENDDO

!                   Scale constants in advance to reduce computation:
                    a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF
             
!               Global debugging:
                IF ( debug >= 10015 ) THEN
                    WRITE(*,*) iat, ' r_cut_off=', r_cut_off, ' fft_radius=', map_1%fft_radius
                ENDIF

!               Apply occupancy AFTER radius estimation:
                a(1:ngauss) = pdb_2%occ(iat) * a(1:ngauss)

!               Figure out box size around atom:
                box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1%map_size ) )

!               Convert to fractional:
                IF ( m > 1 ) THEN
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat) 
                ELSE
                    centre_of_the_atom_in_fractions = map_1%DEORT        * pdb_2%xyz(iat)
                ENDIF
!
!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .x. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size

!               THis might be relevant only at low resolution and very large number of processors(48 or more):
                num1 = nthreads 
                num2 = 1

!               At 2.5 A the box size is normally should be about 24 grid points:
                IF ( nthreads >= 4 ) THEN
  9999              IF ( box_upper_limit(3) - box_lower_limit(3) + 1 < num1 .AND. MOD ( num1 , 2) == 0 &
                         .AND. num1 > num2 ) THEN

                        IF ( debug > 10000 ) THEN
                            CALL messag('Rearranging threads...', srname)
                            WRITE(*,*) box_upper_limit(3) - box_lower_limit(3) + 1, ' nthreads=', nthreads
                        ENDIF
!                       This must be divisible by 2:
                        num1 = num1 / 2
                        num2 = num2 * 2
                        IF ( debug > 10000 ) THEN
                            WRITE(*,*) num1, num2
                        ENDIF

!                       This IF OPERATOR maybe repeated twice or more times if necessary, putting more load on the second loop:
                        GOTO 9999
                    ENDIF
                ENDIF
!               
!               Create a box around the atom and fill it with atom e.d. values:
!
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( BOX_LOWER_LIMIT, BOX_UPPER_LIMIT, MAP_1, CENTRE_OF_THE_ATOM_IN_FRACTIONS, &
!$OMP                                       R2_CUT_OFF, US_PRESENT, VTS, A, BE, NGAUSS, DEBUG) NUM_THREADS(NUM1)
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( BOX_LOWER_LIMIT, BOX_UPPER_LIMIT, MAP_1, CENTRE_OF_THE_ATOM_IN_FRACTIONS, &
!$OMP                                       R2_CUT_OFF, US_PRESENT, VTS, A, BE, NGAUSS, DEBUG) NUM_THREADS(NUM2)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)
                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions   

!                           High point of electron density generation:
                            r2 = map_1%real_tensor * duvw
                            IF ( r2 <= r2_cut_off ) THEN
!OMP CRITICAL
!                                dgp = dgp + 1
!OMP END CRITICAL
!                               Fit uvw vector into map dimensions:
                                jut = uvw .MMOD. map_1%map_size

                                IF ( us_present ) THEN

                                    duvwort = map_1%ORT * duvw

                                    DO l = 1, ngauss
                                        qform(l) = duvwort .DOT. ( VTS(l) * duvwort )
                                    ENDDO

                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -0.5_wp * qform(1:ngauss) ) )

                                ELSE
!                                   Most CPU expensive part:
                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -be(1:ngauss) * r2 ) )
                                ENDIF

!                               Accumulate contibutions:
!OMP CRITICAL
                                map_1%array(jut(1), jut(2), jut(3)) = map_1%array(jut(1), jut(2), jut(3)) + atom_ed_value
!OMP END CRITICAL
                                IF ( debug > 10019 ) THEN
                                    WRITE (*, "(' jut=', 3I5, ' map_1%array=', ES12.5)")&
                                    jut, map_1%array( jut(1), jut(2), jut(3) )
                                ENDIF

                            ENDIF

                        ENDDO    
                    ENDDO
                ENDDO boxz
!$OMP END PARALLEL DO
            ENDDO coord
!OMP END PARALLEL DO
        ENDDO symmetry

        IF ( ALLOCATED ( RU ) ) CALL deallocate_array ( RU )

!        WRITE(*,"(' GENERATE_DENSITY> ', 'Number of grid points per atom: ', I12)") INT ( dgp / REAL(natom, KIND=wp) )
!        WRITE(*,"(' GENERATE_DENSITY> ', 'Number of grid points involved: ', I12)") dgp
!        WRITE(*,"(' GENERATE_DENSITY> ', 'Total memory needed (for speed up) :', I6, ' MB')") &
!        NINT ( dgp * 8 * 1.5 / 10.0 ** 6 )
        CALL CPU_TIME(time1)
        WRITE(*,"(' GENERATE_DENSITY> ', 'CPU time=', F7.3, ' s Elapsed time=', F7.3, ' s')") &
        time1 - time0, SECNDS(t0)
    END SUBROUTINE generate_density

    SUBROUTINE generate_density_nogauss(map_1, pdb_2, selected_atoms)
!
!       Purpose:
!       =======
!       Generates density using 2- or 5-Gaussian approximation
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2003       B.Strokopytov          Original code
!
!       Apr 2005       B.Strokopytov          When calculating raidus cutoff occupancies shall be ignored
!                                             from now on to eliminate negative effects associated with
!                                             small occupancies. We need this for certain algorithms.
!
!       Nov 2005       B.Strokopytov          e_limit has been changes to 18.0 A in constants.f90
!                                             instead of original value of 7.0 which allows
!                                             in conjuction with appropriate anti-alias b_scale
!                                             extremely accurate calculation of structure factors
!                                             (error in phases of 0.001 degr or less ) at any resolution
!                                             using Shannon rate of just 1.5. See (Navaza, 2002, Acta Cryst.)
!                                             Maximum sphere radius is now checked because at low
!                                             resolution it may not even fit into cell
!
!       Oct 2007      B.Strokopytov           e_limit has been replaced with fft_radius now
!                                             present in both MTZ and MAP TYPE.
!
!       Oct 2007      B.Strokopytov           Major rearrangemnet of PDB type. Atomic scattering
!                                             factors are now in PDB_2%ATOMSF_TABLE
!
!       Oct 2007      B.Strokopytov           Improved calculations for max radii (see variable R2_CELL_MAX)
!       Mar 2009      B.Strokopytov           Trying to introduce ADPs
!       
        TYPE(map),                   INTENT(INOUT)        :: map_1
        TYPE(pdb),                   INTENT(IN)           :: pdb_2
        INTEGER,      DIMENSION(:),  INTENT(IN)           :: selected_atoms 
!       Local variables:
        INTEGER                                           :: natom
!       e.d. stuff:
        INTEGER                                           :: ityp      ! atom type 
        INTEGER                                           :: ngauss   
        REAL(KIND=wp), DIMENSION(5)                       :: a
        REAL(KIND=wp), DIMENSION(5)                       :: b 
        REAL(KIND=wp), DIMENSION(5)                       :: be
        REAL(KIND=wp)                                     :: b_scale
        REAL(KIND=wp)                                     :: r_cut_off
        REAL(KIND=wp)                                     :: r2_cut_off
        REAL(KIND=wp)                                     :: r2
        TYPE(vector_int)                                  :: uvw
        TYPE(vector)                                      :: duvw
!       Atom/box parameters:
        TYPE(vector)                                      :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                                  :: closest_grid_point
        TYPE(vector_int)                                  :: box_size
        INTEGER,       DIMENSION(3)                       :: jut
        INTEGER,       DIMENSION(3)                       :: box_lower_limit
        INTEGER,       DIMENSION(3)                       :: box_upper_limit
        REAL(KIND=wp), DIMENSION(3)                       :: ort_diag
        TYPE(vector)                                      :: deort_diag
        REAL(KIND=wp)                                     :: atom_ed_value
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                     :: r2_cell_max
!       Aniso:
        LOGICAL                                           :: aniso_active
        LOGICAL                                           :: Us_present
!       Need to calculate eigenvalues:
        TYPE(matrix)                                      :: UMAT
!       Need U for each gaussian:
        TYPE(aniso),   DIMENSION(5)                       :: U
        TYPE(matrix),  DIMENSION(5)                       :: VTS
        REAL(KIND=wp), DIMENSION(5)                       :: udet
        REAL(KIND=wp), DIMENSION(5)                       :: qform
        REAL(KIND=wp)                                     :: Biso
        REAL(KIND=wp), DIMENSION(3)                       :: ueigen
        REAL(KIND=wp), DIMENSION(6)                       :: v
        TYPE(vector)                                      :: duvwort
!       Counters:
        INTEGER                                           :: iv
        INTEGER                                           :: iat
        INTEGER                                           :: l 
        INTEGER                                           :: m
!       Box counters:
        INTEGER                                           :: ju
        INTEGER                                           :: jv
        INTEGER                                           :: jw
        CHARACTER(LEN=32)                                 :: srname
!       Figure out how many density grid points is involved:
        INTEGER(KIND=eb)                                  :: dgp
        
        dgp = 0
        srname =  'generate_density_nogauss'

!       Start with basic checks:
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

!       Initialize map:
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
        
            coord:DO iv = 1, SIZE ( selected_atoms )
                iat = selected_atoms(iv)

                IF ( pdb_2%occ(iat) == 0.0_wp ) CYCLE

!               Save atom type:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss = pdb_2%atomsf_table%record(ityp)%ngauss

                IF ( ngauss /= SIZE ( pdb_2%atomsf_table%record(ityp) ) ) THEN
                    WRITE(*,*) ngauss, SIZE ( pdb_2%atomsf_table%record(ityp) ), ityp
                    CALL die('Programming error. Inconsistency in ATOMSF_TABLE records.', srname )
                ENDIF

!               Two possibilities here only:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             srname )
                ENDIF

                a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a
                b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b
                Biso = pdb_2%biso(iat)
               
!               Check whether Us have been read (appropriate record present):
                Us_present = .FALSE.

                IF ( aniso_active ) THEN

!                   Note that ">" was overloaded in a special way:
                    Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )

!                   Convert to matrix of aniso elements:
                    IF ( Us_present ) THEN
                        UMAT = pdb_2%U(iat)

!                       Calculate eigenvalues:
                        ueigen = eigen_values ( UMAT )
                        Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )
                        IF ( Biso == 0.0_wp ) THEN
                            WRITE(*,*) ' iat=', iat, ' Biso=', Biso
                            CALL die('Something is very wrong here: Biso == 0.', srname)
                        ENDIF
                    ENDIF

                ENDIF

!               Finally we have sensible approach:
                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + Biso + b_scale ) ** 2
                r2_cut_off = MAX ( map_1%fft_radius ** 2 , r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT( r2_cut_off )

!               Reset gaussians after the radius has been calculated:
                ngauss = 1
                a(1:ngauss)  = 1.0_wp
                b(1:ngauss)  = 0.0_wp

!               Now we can continue transformations:
                IF ( .NOT. Us_present ) THEN

                    a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                    be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

                ELSE

                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian:
                        U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )

!                       Calculate determinant:
                        Udet(l)  = udet_aniso(U(l))

!                       Calculate inverted matrix:
                        v = minors_aniso(U(l)) / Udet(l)
                        VTS(l) = .CONVERT. v
                        
                    ENDDO

!                   Scale constants in advance to reduce computation:
                    a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF
             
!               Global debugging:
                IF ( debug >= 10015 ) THEN
                    WRITE(*,*) iat, ' r_cut_off=', r_cut_off, ' fft_radius=', map_1%fft_radius
                ENDIF

!               Apply occupancy AFTER radius estimation:
                a(1:ngauss) = pdb_2%occ(iat) * a(1:ngauss)

!               Figure out box size around atom:
                box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1%map_size ) )

!               Convert to fractional:
                IF ( m > 1 ) THEN
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat) 
                ELSE
                    centre_of_the_atom_in_fractions = map_1%DEORT        * pdb_2%xyz(iat)
                ENDIF
!
!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .x. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size
!
!               Create a box around the atom and fill it with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)
                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions   

!                           High point of electron density generation:
                            r2 = map_1%real_tensor * duvw
                            IF ( r2 <= r2_cut_off ) THEN
                                dgp = dgp + 1

!                               Fit uvw vector into map dimensions:
                                jut = uvw .MMOD. map_1%map_size

                                IF ( us_present ) THEN

                                    duvwort = map_1%ORT * duvw

                                    DO l = 1, ngauss
                                        qform(l) = duvwort .DOT. ( VTS(l) * duvwort )
                                    ENDDO

                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -0.5_wp * qform(1:ngauss) ) )

                                ELSE
!                                   Most CPU expensive part:
                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -be(1:ngauss) * r2 ) )
                                ENDIF

!                               Accumulate contibutions:
                                map_1%array(jut(1), jut(2), jut(3)) = map_1%array(jut(1), jut(2), jut(3)) + atom_ed_value

                                IF ( debug > 10019 ) THEN
                                    WRITE (*, "(' jut=', 3I5, ' map_1%array=', ES12.5)")&
                                    jut, map_1%array( jut(1), jut(2), jut(3) )
                                ENDIF

                            ENDIF

                        ENDDO    
                    ENDDO
                ENDDO boxz

            ENDDO coord


        WRITE(*,"(' GENERATE_DENSITY> ', 'Number of grid points per atom: ', I12)") &
                  INT ( dgp / REAL ( SIZE ( selected_atoms ), KIND=wp) )
        WRITE(*,"(' GENERATE_DENSITY> ', 'Number of grid points involved: ', I12)") dgp
        WRITE(*,"(' GENERATE_DENSITY> ', 'Total memory needed (for speed up) :', I6, ' MB')") &
        NINT ( dgp * 8 * 1.5 / 10.0 ** 6 )
    END SUBROUTINE generate_density_nogauss

    SUBROUTINE generate_density_for_map_array(map_1, pdb_2, q, true_sp_group)
!
!       Purpose:
!       =======
!       Generates density for a few maps using 2- or 5-Gaussian approximation
!       This routine is needed for normal matrix-vector multiplication.
!       Q is array of quasi occupancies which may assume negative values.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2007       B.Strokopytov          Original code based on GENERATE_DENSITY routine.
!                                             We are dealing with array of maps now.
!
!       Mar 2009       B.Strokopytov          ADPs added.
!       Mar 2009          -"-                 Now estimates r_cut_off using all gaussians.
!
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: map_1
        TYPE(pdb),                                INTENT(IN)           :: pdb_2
        REAL(KIND=wp), DIMENSION(:),              INTENT(IN)           :: q
        TYPE(space_group),                        INTENT(IN), OPTIONAL :: true_sp_group
!       Local variables:
        INTEGER                                                        :: natom
        INTEGER                                                        :: nsym
!       e.d. stuff:
        INTEGER                                                        :: ityp      ! atom type 
        INTEGER                                                        :: ngauss   
        REAL(KIND=wp), DIMENSION(5)                                    :: a
        REAL(KIND=wp), DIMENSION(5)                                    :: b 
        REAL(KIND=wp), DIMENSION(5)                                    :: be
        REAL(KIND=wp)                                                  :: b_scale
        REAL(KIND=wp)                                                  :: r_cut_off
        REAL(KIND=wp)                                                  :: r2_cut_off
        REAL(KIND=wp)                                                  :: r2
        TYPE(vector_int)                                               :: uvw
        TYPE(vector)                                                   :: duvw
!       Atom/box parameters:
        TYPE(vector)                                                   :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                                               :: closest_grid_point
        TYPE(vector_int)                                               :: box_size
        INTEGER,       DIMENSION(3)                                    :: jut
        INTEGER,       DIMENSION(3)                                    :: box_lower_limit
        INTEGER,       DIMENSION(3)                                    :: box_upper_limit
        REAL(KIND=wp), DIMENSION(3)                                    :: ort_diag
        TYPE(vector)                                                   :: deort_diag
        REAL(KIND=wp)                                                  :: atom_ed_value
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                                  :: r2_cell_max
!       Aniso:
        LOGICAL                                                        :: aniso_active
        LOGICAL                                                        :: Us_present
!       Need to calculate eigenvalues:
        TYPE(matrix)                                                   :: UMAT
!       Need U for each gaussian:
        TYPE(aniso),   DIMENSION(5)                                    :: U
        TYPE(matrix),  DIMENSION(5)                                    :: VTS
        REAL(KIND=wp), DIMENSION(5)                                    :: udet
        REAL(KIND=wp), DIMENSION(5)                                    :: qform
        REAL(KIND=wp)                                                  :: Biso
        REAL(KIND=wp), DIMENSION(3)                                    :: ueigen
        REAL(KIND=wp), DIMENSION(6)                                    :: v
        TYPE(vector)                                                   :: duvwort
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE                   :: RU
!       Refinement:
        REAL(KIND=wp)                                                  :: occ
        LOGICAL, DIMENSION(3)                                          :: refine_atomic_xyz      
        LOGICAL                                                        :: refine_atomic_Biso        
        LOGICAL                                                        :: refine_atomic_occ
        LOGICAL, DIMENSION(6)                                          :: refine_atomic_Us
        INTEGER                                                        :: number_of_maps
        INTEGER                                                        :: my_first_allocated_map
!       Pointers:
        INTEGER                                                        :: indxyz
        INTEGER                                                        :: indbiso
        INTEGER                                                        :: indocc
        INTEGER                                                        :: indu
!       Counters:
        INTEGER                                                        :: iat
        INTEGER                                                        :: imap
        INTEGER                                                        :: l 
        INTEGER                                                        :: m
!       Box counters:
        INTEGER                                                        :: ju
        INTEGER                                                        :: jv
        INTEGER                                                        :: jw
        CHARACTER(LEN=32)                                              :: srname
        
        srname = 'generate_density_for_map_array'

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

        IF ( ALLOCATED ( pdb_2%atomsf_table ) ) THEN

            IF ( SIZE ( pdb_2%atomsf_table ) < 1 ) THEN

                WRITE(*,*) SIZE ( pdb_2%atomsf_table )
                CALL die('Programming error. PDB_2 has not been initialized properly.',&
                         srname )
            ENDIF

        ELSE

            CALL die('ATOMSF_TABLE has not been initialized.', srname )

        ENDIF

        IF ( .NOT. ALLOCATED ( pdb_2%atom_type ) ) THEN

            CALL die('Programming error. ATOM_TYPE has not been initialized properly.',&
                     srname )

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

!       Initialize maps:
        DO imap = my_first_allocated_map, SIZE ( map_1 )
            IF ( ALLOCATED ( map_1(imap) ) ) THEN
                map_1(imap)%array = 0.0_wp
            ENDIF
        ENDDO

!       Check parameters:
        DO imap = my_first_allocated_map, SIZE ( map_1 )
            IF ( ALLOCATED ( map_1(imap) ) ) THEN
                IF (  map_1(my_first_allocated_map)%b_scale /= map_1(imap)%b_scale )  THEN
                    WRITE(*,"(2(I2,F6.2))") my_first_allocated_map, map_1(my_first_allocated_map)%b_scale,&
                                            imap, map_1(imap)%b_scale
                   CALL die('Programming error. Unequal B_SCALE in different maps.',&
                            srname )
                ENDIF
            ENDIF
        ENDDO

!       Initialize:
        b_scale    = map_1(my_first_allocated_map)%b_scale

        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )
        ort_diag   = extract_diagonal ( map_1(my_first_allocated_map)%ORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1(my_first_allocated_map)%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1(my_first_allocated_map)%cell
            CALL die ( 'Programming error. Unit cell has not been initialized for current map.',&
                       srname )
        ENDIF

!       Calculate maximum sphere radius around any atom:
        r2_cell_max = 0.25_wp * MINVAL ( ort_diag(1:3) ** 2 )
        
!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
            CALL allocate_array ( RU, nsym )
!           Convert symmetry operators to 6x6 matrix:
            DO m = 1, nsym
                RU(m) = .SIXBYSIX. ( map_1(my_first_allocated_map)%ORT &
                      * (.SYMA. map_1(my_first_allocated_map)%sp_group%sym(m) &
                      * map_1(my_first_allocated_map)%DEORT ) )
            ENDDO

        ELSE
            nsym = 1
        ENDIF

        aniso_active = ALLOCATED ( pdb_2%U )

        symmetry:DO m = 1, nsym

!           Initialise/re-initialise pointer:
            indxyz = 0

            coord:DO iat = 1, natom

!               Number of refinable parameters for this atom:
                l = COUNT ( pdb_2%refinable(1:11, iat) )

!               No refinable parameters -- skip this atom:
                IF ( l == 0 ) CYCLE

!               Save atom occupancy:
                occ = pdb_2%occ(iat)

!               Refinement of zero occncies is possible for occncy only refinement:
                IF ( occ == 0.0_wp .AND. .NOT. pdb_2%refinable(5,iat) ) THEN

!                   If we skip this atom altogether we still need to increment the pointer:
                    indxyz = indxyz + l 
                    CYCLE
                ENDIF

!               Speed it up in case quasi-occupancies for this atom are ALL zero:
                IF ( ALL ( q(indxyz+1:indxyz+l) == 0.0_wp ) ) THEN

!                   Don't forget to increment the pointer though before cycling:
                    indxyz = indxyz + l
                    CYCLE 
                ENDIF

!               Save type of atom:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss = pdb_2%atomsf_table%record(ityp)%ngauss

                IF ( ngauss /= SIZE ( pdb_2%atomsf_table%record(ityp) ) ) THEN

                    WRITE(*,*) ngauss, SIZE ( pdb_2%atomsf_table%record(ityp) )
                    CALL die('Programming error. Inconsistency in ATOMSF_TABLE records.', srname )

                ENDIF

!               Check:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN

                    WRITE(*,*) 'ngauss=', ngauss
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             srname )
                ENDIF

                a(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss)  = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)

!               Always calculate Biso:
                Biso = pdb_2%biso(iat)

!               Check whether Us have been read (appropriate record present):
                Us_present = .FALSE.

                IF ( aniso_active ) THEN

!                   Note that ">" was overloaded in a special way:
                    Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )
                    IF ( Us_present ) THEN
!                       Convert to matrix of aniso elements:
                        UMAT = pdb_2%U(iat)

!                       Calculate eigenvalues:
                        ueigen = eigen_values ( UMAT )
                        Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )
                    ENDIF
                ENDIF


!               Finally we have sensible approach:
                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + Biso + b_scale ) ** 2 
                r2_cut_off = MAX ( map_1(my_first_allocated_map)%fft_radius ** 2 , r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT( r2_cut_off )

                IF ( .NOT. Us_present ) THEN

!                   Note that no occupancies here:
                    a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                    be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

                ELSE


                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian:
                        U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )

!                       Transform ellipsoid using Garib's method:
                        IF ( m /= 1 ) THEN
!                           This is not dangerous for parallel computation since U(l) is of type ANISO
!                           and no array allocation takes place:
                            U(l) = RU(m) * U(l)
                        ENDIF

!                       Calculate determinant:
                        Udet(l)  = udet_aniso(U(l))

!                       Calculate inverted matrix:
                        v = minors_aniso(U(l)) / Udet(l)
                        VTS(l) = .CONVERT. v

                    ENDDO

!                   Scale constants in advance to reduce computation:
                    a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF

!               Figure out box size around atom:
                box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!               Convert to fractional:
                IF ( m > 1 ) THEN
                    centre_of_the_atom_in_fractions = map_1(my_first_allocated_map)%SYM_DEORT(m) * pdb_2%xyz(iat) 
                ELSE
                    centre_of_the_atom_in_fractions = map_1(my_first_allocated_map)%DEORT        * pdb_2%xyz(iat)
                ENDIF

!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size

!               Initialise logicals:
                refine_atomic_xyz  = pdb_2%refinable(1:3, iat)
                refine_atomic_Biso = pdb_2%refinable(4, iat)
                refine_atomic_occ  = pdb_2%refinable(5,iat)
                refine_atomic_Us   = pdb_2%refinable(6:11,iat)

                IF ( refine_atomic_biso ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF
!               
                IF ( refine_atomic_occ ) THEN
                    indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

                IF ( ANY ( refine_atomic_Us ) ) THEN

                    indu = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )

!                   Check consistency of refinable parameters:
                    IF ( refine_atomic_biso ) THEN
                        WRITE(*,*) ' iat=', iat, ' refinable params=', pdb_2%refinable(1:11,iat)
                        CALL die('O-o-o-p-s. Programming error. Both Biso and U''s are refinable.', srname)
                    ENDIF

                ENDIF

!               Create a box around the atom and fill it with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw  = (/ju, jv, jw/)

                            duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_the_atom_in_fractions   

!                           High point of electron density generation:
                            r2 = map_1(my_first_allocated_map)%real_tensor * duvw
                            IF ( r2 <= r2_cut_off ) THEN

!                               Fit uvw vector into map dimensions:
                                jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                IF ( Us_present ) THEN

                                    duvwort = map_1(my_first_allocated_map)%ORT * duvw

                                    DO l = 1, ngauss
                                        qform(l) = duvwort .DOT. ( VTS(l) * duvwort )
                                    ENDDO

                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -0.5_wp * qform(1:ngauss) ) )

                                ELSE

!                                   Isotropic case:
                                    atom_ed_value = SUM ( a(1:ngauss) * EXP ( -be(1:ngauss) * r2 ) )
                                ENDIF

!                               Accumulate contibutions for xyz:
                                any_xyz:IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                   x,y,z maps 1 till 3:
                                    IF ( ALL ( refine_atomic_xyz ) ) THEN

                                        DO imap = 1, 3 
                                            map_1(imap)%array(jut(1), jut(2), jut(3)) = &
                                            map_1(imap)%array(jut(1), jut(2), jut(3)) + q(indxyz+imap) * occ * atom_ed_value
                                        ENDDO

                                    ELSE

!                                       If some coordinates unrefinable do it one by one:
                                        l = 0
                                        DO imap = 1, 3
                                            IF ( refine_atomic_xyz(imap) ) THEN
                                                l = l + 1

!                                               "X"-q  always goes to X-map (map #1), "Y"-q -- to Y-map (map #2) etc.:
                                                map_1(imap)%array(jut(1), jut(2), jut(3)) = &
                                                map_1(imap)%array(jut(1), jut(2), jut(3)) + q(indxyz+l) * occ * atom_ed_value
                                            ENDIF
                                        ENDDO
                                    
                                    ENDIF
                                ENDIF any_xyz

!                               4th map:
                                IF ( refine_atomic_biso ) THEN

                                    map_1(4)%array(jut(1), jut(2), jut(3)) = &
                                    map_1(4)%array(jut(1), jut(2), jut(3)) + q(indbiso) * occ * atom_ed_value

                                ENDIF

!                               5th map:
                                IF ( refine_atomic_occ ) THEN

                                    map_1(5)%array(jut(1), jut(2), jut(3)) = &
                                    map_1(5)%array(jut(1), jut(2), jut(3)) + q(indocc) * atom_ed_value

                                ENDIF

!                               6th - 11th map:
                                aniso_all:IF ( ALL ( refine_atomic_Us ) ) THEN

                                    DO  imap = 1, 6
                                        map_1(imap+5)%array(jut(1), jut(2), jut(3)) = &
                                        map_1(imap+5)%array(jut(1), jut(2), jut(3)) + q(indu+imap) * occ * atom_ed_value
                                    ENDDO

                                ELSE IF ( ANY ( refine_atomic_Us ) ) THEN

                                    l = 0
                                    DO  imap = 1, 6
                                        IF ( refine_atomic_Us(imap) ) THEN
                                            l = l + 1
                                            map_1(imap+5)%array(jut(1), jut(2), jut(3)) = &
                                            map_1(imap+5)%array(jut(1), jut(2), jut(3)) + q(indu+l) * occ * atom_ed_value
                                        ENDIF
                                    ENDDO

                                ENDIF aniso_all

                                IF ( debug > 10015 ) THEN
                                    WRITE (*, "(' jut=', 3I5, ' map_1%array=', ES12.5)")&
                                    jut, map_1(my_first_allocated_map)%array( jut(1), jut(2), jut(3) )
                                ENDIF

                            ENDIF

                        ENDDO    
                    ENDDO
                ENDDO boxz

!               Increment pointer:
                indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )

            ENDDO coord

        ENDDO symmetry

!       Free memory:
        IF ( ALLOCATED ( RU ) ) CALL deallocate_array ( RU )

    END SUBROUTINE generate_density_for_map_array

    SUBROUTINE generate_logical_density ( map_1, pdb_2, nsym_start, nsym_finish )
!
!       Purpose:
!       =======
!
!       Calculates a mask surrounding atoms within a certain distance
!       defined by map_1%radius (3.0 is close to optimal we think)
!       but side chain atoms beyond CB should be ignored (see below).
!
!       map_1 must be allocated before calling this subroutine
!       pdb_2 must be allocated and initialized in advance too

!       nsym_start & nsym_finish can be used as following:
!       a) For calculation in P1 space group           : nsym_start = 1; nsym_finish = 1
!       b) For calculation in TRUE space group         : nsym_start = 1; nsym_finish = map_1%sp_group%number_of_symops
!       c) For calculation in REDUCED TRUE space group : nsym_start = 2; nsym_finish = map_1%sp_group%number_of_symops
!          Last case is useful to estimate self-ovelap when comparing maps from case a) and c)                      
!
!       Date:           Programmer:             Description of changes:
!       ====            ==========              ======================
!       Dec 2003        B.Strokopytov           Original code
!       Nov 2005        B.Strokopytov           map_1%atom_name added and used now
!
        TYPE(map_logical), INTENT(INOUT)  :: map_1
        TYPE(pdb),         INTENT(IN)     :: pdb_2
        INTEGER,           INTENT(IN)     :: nsym_start
        INTEGER,           INTENT(IN)     :: nsym_finish
!       Local variables:
        INTEGER                           :: natom
        REAL(KIND=wp)                     :: r_cut_off
        REAL(KIND=wp)                     :: r2_cut_off
        REAL(KIND=wp)                     :: r2
        TYPE(vector)                      :: duvw
        TYPE(vector_int)                  :: uvw
!       Atom box params:
        TYPE(vector)                      :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                  :: closest_grid_point
        TYPE(vector_int)                  :: box_size
        INTEGER, DIMENSION(3)             :: jut
        INTEGER, DIMENSION(3)             :: box_lower_limit
        INTEGER, DIMENSION(3)             :: box_upper_limit
        TYPE(vector)                      :: deort_diag
        INTEGER, DIMENSION(3)             :: v
!       Counters:
        INTEGER                           :: m
        INTEGER                           :: iat
!       Box counters:
        INTEGER                           :: ju
        INTEGER                           :: jv
        INTEGER                           :: jw

!       Checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die ( 'Programming error. PDB_2 has not been allocated properly...',&
                       'generate_logical_density' ) 
        ENDIF
        IF ( natom < 1 ) THEN
            WRITE ( *, "(' GENERATE_LOGICAL_DENSITY> ', 'natom=',I6)") natom 
            CALL die ( 'Programming error. Zero or negative number of atoms.', &
                       'generate_logical_density' )
        ENDIF

!       Initialize map:
        map_1%ntrue = 0
        IF ( ALLOCATED( map_1%array ) ) THEN
            map_1%array = .FALSE.
        ELSE
            CALL die ( 'Programming error. Memory for logical map_1 has to be allocated in advance.',&
                       'generate_logical_density' )
        ENDIF

!       Initialize parameters:
        deort_diag = extract_diagonal ( map_1%DEORT )
        r_cut_off  = map_1%radius
        r2_cut_off = r_cut_off ** 2
        box_size   = NINT( r_cut_off * ( deort_diag .x. map_1%map_size ) )
        symm:DO m = nsym_start, nsym_finish
            atom:DO iat = 1, natom
                IF ( pdb_2%occ(iat) == 0.0_wp) CYCLE atom
                IF ( map_1%atom_name(1)(1:3) /= 'ALL' ) THEN
                    IF ( .NOT. ANY ( pdb_2%atom_name(iat) == map_1%atom_name ) ) CYCLE atom
                ENDIF

!               Convert to fractions:
                centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat)

!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .x. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size
                IF ( debug > 10006 ) THEN
                    WRITE(*,*) ' r_cut_off                           = ', r_cut_off
                    v = box_size
                    WRITE(*,*) ' box_size                            = ', v
                    v = closest_grid_point
                    WRITE(*,*) ' closest grid point to atomic centre = ', v
                    WRITE(*,*) ' box_lower_limit                     = ', box_lower_limit
                    WRITE(*,*) ' box_upper_limit                     = ', box_upper_limit
                ENDIF

!               Create a box around the atom and fill it with atom e.d. values:
                DO jw = box_lower_limit(3), box_upper_limit(3)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)

!                           Set uvw vector:
                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions

!                           Check distance:
                            r2 = map_1%REAL_TENSOR * duvw
                            IF ( r2 <= r2_cut_off ) THEN
                                jut  = uvw .MMOD. map_1%map_size

!                               We want to calculate true element only once:
                                IF ( .NOT. map_1%array ( jut(1), jut(2), jut(3) ) ) THEN
                                     map_1%ntrue =  map_1%ntrue +  1
                                     map_1%array( jut(1), jut(2), jut(3) ) = .TRUE.
                                ENDIF
                            ENDIF
                        ENDDO    
                    ENDDO
                ENDDO

            ENDDO atom
        ENDDO symm
    END SUBROUTINE generate_logical_density

    SUBROUTINE generate_constant_density ( map_1, pdb_2, sp_group, radius )
!
!       Purpose:
!       =======
!       Calculates a mask surrounding atoms within a certain distance
!       defined by map_1%radius. 2.5-2.6 A is close to optimal we think for complete models.
!
!       Note:
!       ====
!       map_1 must be allocated before calling this subroutine
!       pdb_2 must be allocated and initialized in advance too
!
!       Date:           Programmer:               Description of changes:
!       =====           ==========                ======================
!       Apr 2004        B.Strokopytov             Original code
!
        TYPE(map),         INTENT(INOUT)        :: map_1
        TYPE(pdb),         INTENT(IN)           :: pdb_2
        TYPE(space_group), INTENT(IN), OPTIONAL :: sp_group
        REAL(KIND=wp),     INTENT(IN), OPTIONAL :: radius
!       Local variables:
        INTEGER                                 :: natom
        INTEGER                                 :: nsym_finish
!       Electron density atom parameters:
        REAL(KIND=wp)                           :: r_cut_off
        REAL(KIND=wp)                           :: r2_cut_off
        REAL(KIND=wp)                           :: r2
        TYPE(vector_int)                        :: uvw
        TYPE(vector)                            :: duvw
        TYPE(vector)                            :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                        :: closest_grid_point
        TYPE(vector_int)                        :: box_size
        INTEGER, DIMENSION(3)                   :: jut
        INTEGER, DIMENSION(3)                   :: box_lower_limit
        INTEGER, DIMENSION(3)                   :: box_upper_limit
        TYPE(vector)                            :: deort_diag
        INTEGER, DIMENSION(3)                   :: v
!       Counters:
        INTEGER                                 :: m
        INTEGER                                 :: iat
!       Box counters
        INTEGER                                 :: ju
        INTEGER                                 :: jv
        INTEGER                                 :: jw
        CHARACTER(LEN=32),                 SAVE :: srname='generate_constant_density'
!       Experimeting with threads:
        INTEGER                                 :: nthreads
        INTEGER                                 :: num1
        INTEGER                                 :: num2
!       Value of constant density:
        REAL(KIND=wp)                           :: const
        REAL(KIND=sp)                           :: t0
        REAL(KIND=sp)                           :: cpu0
        REAL(KIND=sp)                           :: cpu1
        REAL(KIND=wp)                           :: r_cut

!       Timing:
        t0 = SECNDS(0.0)
        CALL CPU_TIME(cpu0)

!       Checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die ( 'PDB arrays have not been allocated properly...', srname ) 
        ENDIF

        IF ( natom < 1 ) THEN
            WRITE(*,"(' GENERATE_CONSTANT_DENSITY> ', 'natom= ',A)") TRIM ( int_to_c (natom ) )
            CALL die(' Zero or negative number of atoms.', srname)
        ENDIF

!       Initialize map:
        IF ( ALLOCATED( map_1%array ) ) THEN
            map_1%array = 0.0_wp
        ELSE
            CALL die('Memory for  map_1 has to be allocated in advance.', srname)
        ENDIF

        IF ( PRESENT ( sp_group ) ) THEN
            nsym_finish = sp_group%number_of_symops
        ELSE
            nsym_finish = 1
        ENDIF

!       Initialize parameters:
        IF ( PRESENT ( radius ) ) THEN
            r_cut_off = radius
        ELSE
            r_cut_off = 2.50_wp
        ENDIF

        deort_diag = extract_diagonal ( map_1%DEORT )
!        r2_cut_off = r_cut_off ** 2
!        box_size   = NINT( r_cut_off * ( deort_diag .x. map_1%map_size ) )

        CALL get_num_threads(nthreads)

        DO m = 1, nsym_finish
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( M, NSYM_FINISH, NATOM, MAP_1, PDB_2, DEORT_DIAG, R_CUT_OFF, DEBUG, SRNAME) 
            DO iat = 1, natom

!                const = pdb_2%occ(iat)

!               Prepare for variable radius:
                r_cut = r_cut_off * (pdb_2%occ(iat) ** (1.0_wp / 3.0_wp))
                box_size   = NINT( r_cut * ( deort_diag .x. map_1%map_size ) )
                r2_cut_off = r_cut ** 2

!               Need to divide to avoid scaling of other arrays, e.g. FPART:
                const = 1.0_wp / nsym_finish

!               Another possibility is:
!                const= pdb_2%occ(iat)/nsym_finish
!               But it is not clear what to do with overlapping atoms...

                IF ( const == 0.0_wp) CYCLE

!               Convert to fractions:
                centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat)

!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .x. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size

                IF ( debug > 10006 ) THEN
                    WRITE(*,*) ' r_cut_off                           = ', r_cut_off
                    v = box_size
                    WRITE(*,*) ' box_size                            = ', v
                    v = closest_grid_point
                    WRITE(*,*) ' closest grid point to atomic centre = ', v
                    WRITE(*,*) ' box_lower_limit                     = ', box_lower_limit
                    WRITE(*,*) ' box_upper_limit                     = ', box_upper_limit
                ENDIF

!               THis might be relevant only at low resolution and very large number of processors(48 or more):
                num1 = nthreads 
                num2 = 1

!               At 2.5 A the box size is normally should be about 24 grid points:
                IF ( nthreads >= 4 ) THEN

  9999              IF ( box_upper_limit(3) - box_lower_limit(3) + 1 < num1 .AND. MOD ( num1 , 2) == 0 &
                         .AND. num1 > num2 ) THEN

                        IF ( debug > 10000 ) THEN
                            CALL messag('Rearranging threads...', srname)
                            WRITE(*,*) box_upper_limit(3) - box_lower_limit(3) + 1, ' nthreads=', nthreads 
                        ENDIF

!                       This must be divisible by 2:
                        num1 = num1 / 2
                        num2 = num2 * 2
                        IF ( debug > 10000 ) THEN
                            WRITE(*,*) num1, num2
                        ENDIF
!                       This IF OPERATOR maybe repeated twice or more times if necessary, putting more load on the second loop:
                        GOTO 9999

                    ENDIF
                ENDIF
!               

!               Create a box around the atom and fill it with atom e.d. values:
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( BOX_LOWER_LIMIT, BOX_UPPER_LIMIT, MAP_1, CENTRE_OF_THE_ATOM_IN_FRACTIONS, &
!OMP                                       R2_CUT_OFF, DEBUG, SRNAME, CONST) NUM_THREADS(NUM1)
                DO jw = box_lower_limit(3), box_upper_limit(3)
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( BOX_LOWER_LIMIT, BOX_UPPER_LIMIT, MAP_1, CENTRE_OF_THE_ATOM_IN_FRACTIONS, &
!OMP                                       R2_CUT_OFF, DEBUG, SRNAME, CONST) NUM_THREADS(NUM2)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)
                            uvw = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions
                            r2 = map_1%REAL_TENSOR * duvw
                            IF ( r2 <= r2_cut_off ) THEN

                                jut = uvw .MMOD. map_1%map_size

                                IF ( map_1%array ( jut(1), jut(2), jut(3) ) == 0.0_wp )  THEN

!                                    Set density to occupancy of this particular atom:
!                                     map_1%array( jut(1), jut(2), jut(3) ) = const
                                     map_1%array( jut(1), jut(2), jut(3) ) = const
                                ENDIF

                            ENDIF
                        ENDDO    
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

!       Execution time:
        CALL CPU_TIME(cpu1)
        WRITE(*,"(' GENERATE_CONSTANT_DENSITY> ', 'CPU time=', F7.3, ' s Elapsed time=', F7.3, ' s')")&
        cpu1 - cpu0, SECNDS(t0)
    END SUBROUTINE generate_constant_density

    SUBROUTINE shave_constant_density(map_1, map_2, pdb_2, neighbors, radius, shrink_rad, sp_group)
!
!       Purpose:
!       =======
!       Calculates a mask surrounding atoms within a certain distance
!       defined by map_1%radius. 2.5-2.6 A is close to optimal we think for complete models.
!
!       Note:
!       ====
!       map_1 must be allocated before calling this subroutine
!       pdb_2 must be allocated and initialized in advance too
!
!       Date:           Programmer:               Description of changes:
!       =====           ==========                ======================
!       Apr 2004        B.Strokopytov             Original code
!
        TYPE(map),                          INTENT(INOUT)        :: map_1
        TYPE(map_logical),                  INTENT(INOUT)        :: map_2
        TYPE(pdb),                          INTENT(IN)           :: pdb_2
        TYPE(adjacency_list), DIMENSION(:), INTENT(IN)           :: neighbors
!       Need info on VDW individual atom radius:        
        REAL(KIND=wp),                      INTENT(IN)           :: radius
        REAL(KIND=wp),                      INTENT(IN)           :: shrink_rad
        TYPE(space_group),                  INTENT(IN), OPTIONAL :: sp_group
!       Local variables:
        INTEGER                                                  :: natom
        INTEGER                                                  :: nsym_finish
        INTEGER                                                  :: neigh
!       Electron density atom parameters:
        REAL(KIND=wp)                                            :: r_cut_off
        REAL(KIND=wp)                                            :: r_cut_off_min
        REAL(KIND=wp)                                            :: r2_cut_off
        REAL(KIND=wp)                                            :: r2_cut_off_min
        REAL(KIND=wp)                                            :: r2
        TYPE(vector_int)                                         :: uvw
        TYPE(vector)                                             :: duvw
        TYPE(vector)                                             :: centre_of_the_atom_in_fractions 
        TYPE(vector_int)                                         :: closest_grid_point
        TYPE(vector_int)                                         :: box_size
        INTEGER, DIMENSION(3)                                    :: jut
        INTEGER, DIMENSION(3)                                    :: box_lower_limit
        INTEGER, DIMENSION(3)                                    :: box_upper_limit
        TYPE(vector)                                             :: deort_diag
        INTEGER, DIMENSION(3)                                    :: v
!       Counters:
        INTEGER                                                  :: m
        INTEGER                                                  :: iat
        INTEGER                                                  :: j
        INTEGER                                                  :: jat
!       Box counters
        INTEGER                                                  :: ju
        INTEGER                                                  :: jv
        INTEGER                                                  :: jw
!       Test:
        INTEGER                                                  :: painted
        INTEGER                                                  :: repainted
        INTEGER                                                  :: shaved
        REAL(KIND=wp)                                            :: solvent
!       Checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die ( 'PDB arrays have not been allocated properly...', 'generate_constant_density' ) 
        ENDIF

        IF ( natom < 1 ) THEN
            WRITE(*,"(' GENERATE_CONSTANT_DENSITY> ', 'natom =',I6)") natom 
            CALL die(' Zero or negative number of atoms.', 'generate_constant_density')
        ENDIF

!       Initialize map:
        IF ( .NOT. ALLOCATED( map_1%array ) ) THEN
            CALL die('Memory for  map_1 has to be allocated in advance.', 'generate_constant_density')
        ENDIF

        IF ( PRESENT ( sp_group ) ) THEN
            nsym_finish = sp_group%number_of_symops
        ELSE
            nsym_finish = 1
        ENDIF

        v = map_1%map_size
        solvent = REAL ( v(1) * v(2) * v(3), KIND=wp)
        WRITE(*,*) ' grid volume=', solvent
        WRITE(*,*) ' model grid=', COUNT(map_1%array == 1.0_wp)
        solvent = 1.0_wp - COUNT(map_1%array == 1.0_wp) * map_1%sp_group%number_of_symops  / solvent
        
        WRITE(*,"(' Solvent content (%):', F8.1)") solvent * 100

        IF ( .NOT. ALLOCATED ( map_2%array ) ) THEN
            CALL die('Memory for  map_1 has to be allocated in advance.', 'generate_constant_density') 
        ELSE
            map_2%array = .FALSE.
        ENDIF

!       FIXME: This should depend upon atom type:
        r_cut_off     = radius
        r_cut_off_min = radius - shrink_rad

        deort_diag = extract_diagonal ( map_1%DEORT )

        r2_cut_off     = r_cut_off ** 2
        r2_cut_off_min = r_cut_off_min ** 2

        box_size   = NINT ( r_cut_off * ( deort_diag .x. map_1%map_size ) ) 

        DO iat = 1, natom

            painted = 0
            repainted = 0
            shaved = 0

            neigh =  SIZE ( neighbors(iat)%v )

            DO m = 1, nsym_finish

!               Loop over neighboring atoms:
                paint:DO j = 1, neigh

                    jat = neighbors(iat)%v(j)

                    IF ( pdb_2%occ(jat) == 0.0_wp) CYCLE

!                   Convert to fractions:
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(jat)

!                   This is our atomic centre on the grid:
                    closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1%map_size )

                    box_lower_limit = closest_grid_point - box_size
                    box_upper_limit = closest_grid_point + box_size

                    IF ( debug > 10006 ) THEN
                        WRITE(*,*) ' r_cut_off                           = ', r_cut_off
                        v = box_size
                        WRITE(*,*) ' box_size                            = ', v
                        v = closest_grid_point
                        WRITE(*,*) ' closest grid point to atomic centre = ', v
                        WRITE(*,*) ' box_lower_limit                     = ', box_lower_limit
                        WRITE(*,*) ' box_upper_limit                     = ', box_upper_limit
                    ENDIF

!                   Create a box around the atom and fill it with atom e.d. values:
                    DO jw = box_lower_limit(3), box_upper_limit(3)
                        DO jv = box_lower_limit(2), box_upper_limit(2)
                            DO ju = box_lower_limit(1), box_upper_limit(1)
                                uvw = (/ju, jv, jw/)
                                duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions
                                r2 = map_1%REAL_TENSOR * duvw
                                IF ( r2 <= r2_cut_off ) THEN
                                    jut = uvw .MMOD. map_1%map_size
                                    IF ( map_1%array ( jut(1), jut(2), jut(3) ) == 1.0_wp )  THEN
                                         map_1%array( jut(1), jut(2), jut(3) ) = 2.0_wp
                                         painted = painted + 1
                                    ENDIF
                                ENDIF
                            ENDDO    
                        ENDDO
                    ENDDO

                ENDDO paint

                IF ( pdb_2%occ(iat) == 0.0_wp) CYCLE

!               Convert to fractions:
                centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat)

!               This is our atomic centre on the grid:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .x. map_1%map_size )
                box_lower_limit = closest_grid_point - box_size
                box_upper_limit = closest_grid_point + box_size

                IF ( debug > 10006 ) THEN
                    WRITE(*,*) ' r_cut_off                           = ', r_cut_off
                    v = box_size
                    WRITE(*,*) ' box_size                            = ', v
                    v = closest_grid_point
                    WRITE(*,*) ' closest grid point to atomic centre = ', v
                    WRITE(*,*) ' box_lower_limit                     = ', box_lower_limit
                    WRITE(*,*) ' box_upper_limit                     = ', box_upper_limit
                ENDIF

!               Create a box around the atom and shave it if the density is not from neighboring atoms:
                DO jw = box_lower_limit(3), box_upper_limit(3)
                    DO jv = box_lower_limit(2), box_upper_limit(2)
                        DO ju = box_lower_limit(1), box_upper_limit(1)
                            uvw = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions
                            r2 = map_1%REAL_TENSOR * duvw
                            IF ( r2 > r2_cut_off_min .AND. r2 <= r2_cut_off ) THEN
                                jut = uvw .MMOD. map_1%map_size
                                IF ( map_1%array ( jut(1), jut(2), jut(3) ) == 1.0_wp )  THEN

!                                   Save position of grid points to be removed from the mask:
                                    map_2%array( jut(1), jut(2), jut(3) ) = .TRUE.
                                    shaved = shaved + 1
                                ENDIF
                             ENDIF
                         ENDDO    
                    ENDDO
                ENDDO

!               Repeat the procedure and repaint neighboring atoms:
!               Loop over neighboring atoms:
                repaint:DO j = 1, neigh
                    jat = neighbors(iat)%v(j)

                    IF ( pdb_2%occ(jat) == 0.0_wp) CYCLE

!                   Convert to fractions:
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(jat)

!                   This is our atomic centre on the grid:
                    closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1%map_size )
                    box_lower_limit = closest_grid_point - box_size
                    box_upper_limit = closest_grid_point + box_size

                    IF ( debug > 10006 ) THEN
                        WRITE(*,*) ' r_cut_off                           = ', r_cut_off
                        v = box_size
                        WRITE(*,*) ' box_size                            = ', v
                        v = closest_grid_point
                        WRITE(*,*) ' closest grid point to atomic centre = ', v
                        WRITE(*,*) ' box_lower_limit                     = ', box_lower_limit
                        WRITE(*,*) ' box_upper_limit                     = ', box_upper_limit
                    ENDIF

!                   Create a box around the atom and fill it with atom e.d. values:
                    DO jw = box_lower_limit(3), box_upper_limit(3)
                        DO jv = box_lower_limit(2), box_upper_limit(2)
                            DO ju = box_lower_limit(1), box_upper_limit(1)
                                uvw = (/ju, jv, jw/)
                                duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions
                                r2 = map_1%REAL_TENSOR * duvw
                                IF ( r2 <= r2_cut_off ) THEN
                                    jut = uvw .MMOD. map_1%map_size
                                    IF ( map_1%array ( jut(1), jut(2), jut(3) ) == 2.0_wp )  THEN
                                         map_1%array( jut(1), jut(2), jut(3) ) = 1.0_wp
                                         repainted = repainted + 1
                                    ENDIF
                                ENDIF
                            ENDDO    
                        ENDDO
                    ENDDO

                ENDDO repaint

!           Symmetry:
            ENDDO

!            IF ( debug > 10000 ) THEN
                WRITE(*,"(' atom no.=', I6, 10X, ' painted=', I6, ' repainted=', I6, ' shaved=', I6)")&
                           iat, painted, repainted, shaved
!            ENDIF
 
        ENDDO
        
        v = map_1%map_size
!       FIXME move to the start of the routine:
        IF ( map_1%map_size /= map_2%map_size ) THEN
            CALL die('Map sizes differ.', 'shave_const_den')
        ENDIF
        DO ju = 0, v(1) - 1
            DO jv = 0, v(2) - 1
                DO jw = 0, v(3) - 1
                    IF ( map_2%array(ju, jv, jw) ) THEN
                        map_1%array(ju, jv, jw) = 0.0_wp
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

!       FIXME do this in scaling routine:
        solvent = v(1) * v(2) * v(3)
        solvent = 1.0_wp - COUNT(map_1%array == 1.0_wp) * map_1%sp_group%number_of_symops  / solvent
        WRITE(*,"(' Solvent content (%) after shaving:', F8.1)") solvent * 100

        IF ( COUNT (map_1%array == 2.0_wp ) > 0 ) THEN
            STOP 'STOP BUG'
        ENDIF

    END SUBROUTINE shave_constant_density

END MODULE genden 
