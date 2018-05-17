MODULE convolution
USE aniso_manip
USE aniso_symmetry_manip
USE atom_images
USE convolution_util
USE map_operations
USE refinement_util
!USE sparse_basic, ONLY:coo_matrix
USE sparse_basic
USE OMP_LIB
IMPLICIT NONE
CONTAINS
    SUBROUTINE first_derivatives_wrt_xyzbqU ( g, map_1, pdb_2, mode, true_sp_group )
!
!       Purpose:
!       =======
!       Gradient calculation. This sub convolutes atomic density derivatives
!       with Fourier map using `ngauss' Gaussian approximation.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2003       B.Strokopytov          Original code
!
!       Apr 2005       B.Strokopytov          When calculating radius cut-off occupancies shall be ignored
!                                             from now on to eliminate negative effects associated with
!                                             small occupancies. We need this for certain algorithms.
!
!       Nov 2005       B.Strokopytov          e_limit has been changes to 18.0 A in constants.f90
!                                             instead of original value of 7.0 which allows
!                                             in conjuction with appropriate anti-alias b_scale
!                                             extremely accurate calculation of structure factors
!     
!       Oct 2007       B.Strokopytov          atomsf_table introduced ( memory economy )
!                                             plus flexible number of gaussians
!
!       Oct 2007       B.Strokopytov          J.Navaza's FFT radius finally added                   
!                                             e_limit tumbai.
!
!       Oct 2007       B.Strokopytov          Derivatives against Biso and occupancies added.                    
!
!       Oct 2007       B.Strokopytov          ort_diag added for calculation of maximum
!                                             inscribed sphere size
!
!       Oct 2007       B.Strokopytov          EXP_BE array added for speed up of Biso/occ refinement
!
!       Oct 2007       B.Strokopytov          More careful treatment of occupancies ( division by occ removed).
!
!       Oct 2007       B.Strokopytov          PDB_2%REFINABLE array added for generality and simplicity
!
!       Mar 2007       B.Strokopytov          ADPs have been added (anisotropic refinement)   
!
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT)        :: g
        TYPE(map),                   INTENT(IN)           :: map_1
        TYPE(pdb),                   INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),            INTENT(IN)           :: mode
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
        REAL(KIND=wp), DIMENSION(5)                       :: bb
        REAL(KIND=wp), DIMENSION(5)                       :: exp_be 
        REAL(KIND=wp)                                     :: b_scale
        REAL(KIND=wp)                                     :: r_cut_off
        REAL(KIND=wp)                                     :: r2_cut_off
        REAL(KIND=wp)                                     :: r2
        TYPE(vector_int)                                  :: uvw
        TYPE(vector)                                      :: duvw
        REAL(KIND=wp)                                     :: map_value
        REAL(KIND=wp)                                     :: occ
!       Atom/box parameters:
        TYPE(vector)                                      :: centre_of_the_atom_in_fractions
        TYPE(vector_int)                                  :: closest_grid_point
        TYPE(vector_int)                                  :: box_size
        INTEGER, DIMENSION(3)                             :: jut
        INTEGER, DIMENSION(3)                             :: box_lower_limit
        INTEGER, DIMENSION(3)                             :: box_upper_limit
        TYPE(vector)                                      :: deort_diag
        TYPE(vector)                                      :: ort_diag
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                     :: r2_cell_max
!       Aniso:
        LOGICAL                                           :: aniso_active
        LOGICAL                                           :: Us_present
        TYPE(matrix)                                      :: UMAT
        TYPE(aniso),  DIMENSION(5)                        :: U
        TYPE(matrix), DIMENSION(5)                        :: VTS
        REAL(KIND=wp), DIMENSION(5)                       :: udet
        REAL(KIND=wp), DIMENSION(5)                       :: aexp_qform
        REAL(KIND=wp)                                     :: Biso
        REAL(KIND=wp), DIMENSION(3)                       :: ueigen
!       Need this for every grid point:
        REAL(KIND=wp), DIMENSION(3)                       :: duvwort
        REAL(KIND=wp), DIMENSION(3,5)                     :: VTS_duvwort
        REAL(KIND=wp), DIMENSION(3)                       :: dRodX
!       Aniso ADPs derivatives section:
        REAL(KIND=wp), DIMENSION(6)                       :: dRodU
        REAL(KIND=wp), DIMENSION(6,5)                     :: v
        REAL(KIND=wp), DIMENSION(6,6,5)                   :: vtsd
        REAL(KIND=wp), DIMENSION(6)                       :: xxt
        REAL(KIND=wp), DIMENSION(6,5)                     :: utemp
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE      :: RU
!       Refinement variables:
        LOGICAL                                           :: refine_xyz
        LOGICAL                                           :: refine_biso
        LOGICAL                                           :: refine_occ
        LOGICAL                                           :: refine_Us
        LOGICAL                                           :: refine_atomic_occ
        LOGICAL                                           :: refine_atomic_Biso
        LOGICAL, DIMENSION(3)                             :: refine_atomic_xyz
        LOGICAL, DIMENSION(6)                             :: refine_atomic_Us
!       Pointers:
        INTEGER                                           :: indxyz
        INTEGER                                           :: indbiso
        INTEGER                                           :: indocc
        INTEGER                                           :: indu
!       Counters:
        INTEGER                                           :: iat
        INTEGER                                           :: k 
        INTEGER                                           :: l 
        INTEGER                                           :: m
!       Box counters:
        INTEGER                                           :: ju
        INTEGER                                           :: jv
        INTEGER                                           :: jw
        CHARACTER(LEN=32)                                 :: srname

!       Historic moment:
        srname = 'first_derivatives_wrt_xyzbqU'

        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

!       Start with basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname)
        ENDIF

        IF ( .NOT. ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL die('Programming error. PDB_2%REFINABLE has not been initialized properly.',&
                     srname)            
        ELSE
            IF ( SIZE ( g ) /= COUNT ( pdb_2%refinable ) ) THEN
                CALL die('Programming error. Inconsistent sizes for G and PDB_2%REFINABLE arrays.',&
                srname)
            ENDIF
        ENDIF

!       Check whether aniso records present:
        aniso_active = ALLOCATED ( pdb_2%U )

!       Initialize parameters:
        g          = 0.0_wp
        b_scale    = map_1%b_scale
        deort_diag = extract_diagonal ( map_1%DEORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1%cell
            CALL die('Programming error. Unit cell has not been initialized for current map.',&
                     srname)
        ENDIF

!       Calculate maximum sphere radius around any atom:
        ort_diag = extract_diagonal ( map_1%ORT )

!       Max radius should not exceed 1/2 * cell in any direction:
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

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

        symmetry:DO m = 1, nsym

! Last line in OMP instruction is not exactly necessary...
!$OMP PARALLEL DO DEFAULT ( PRIVATE ) SHARED ( M, NSYM, NATOM, MAP_1, PDB_2, B_SCALE, ANISO_ACTIVE, &
!$OMP                                          ORT_DIAG, DEORT_DIAG, R2_CELL_MAX, RU, G, DEBUG, SRNAME, &
!$OMP                                          REFINE_XYZ, REFINE_BISO, REFINE_OCC, REFINE_US )
            coord:DO iat = 1, natom

!               Need this because of occupancy refinement:
                occ = pdb_2%occ(iat)
                IF ( occ == 0.0_wp .AND. .NOT. refine_occ ) THEN
                    CYCLE
                ENDIF

!               Save type of atom:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss =  pdb_2%atomsf_table%record(ityp)%ngauss 

!               Check:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             'generate_density')
                ENDIF

!               Copy gaussian constants:
                a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)
                Biso = pdb_2%biso(iat)

                Us_present = .FALSE.

                IF ( aniso_active ) THEN
                    Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )
                    UMAT = pdb_2%U(iat)
!$OMP CRITICAL
                    ueigen = eigen_values ( UMAT )
!$OMP END CRITICAL
                    Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )
                ENDIF

!               Finally we have sensible approach:
!                l = MAXLOC ( ABS ( a(1:ngauss) ), DIM=1 )
                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + biso + b_scale ) ** 2

!               Make sure that minimal radius is no less than Navaza's:
                r2_cut_off = MAX ( map_1%fft_radius ** 2, r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )
                
!               This is for Biso refinement only:
                IF ( refine_biso ) THEN
                    bb(1:ngauss) = 1.0_wp / ( b(1:ngauss) + biso + b_scale )
                ENDIF

                IF ( .NOT. Us_present ) THEN

                    a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                    be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

                ELSE

                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian:
                        U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )

!                       Transform ellipsoid using Garib's method:
                        IF ( m /= 1 ) THEN
                            U(l) = RU(m) * U(l)
                        ENDIF
                        
!                       Calculate determinant:
                        Udet(l)  = udet_aniso ( U(l) )

!                       Calculate inverted matrix:
                        VTS(l) = .CONVERT. ( minors_aniso ( U(l) ) / Udet(l) )

!                       U derivative section:
                        DO k = 1, 6
                            VTSD(k,1:6,l) = VTS_first_deriv ( U(l), k )
                        ENDDO
                        V(1:6,l) = det_first_deriv_aniso ( U(l) ) / UDET(l)              
                    ENDDO

!                   Scale constants in advance to reduce computation:
                    a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF

!               Make sure that minimal radius is no less than Navaza's:
                r2_cut_off = MAX ( map_1%fft_radius ** 2, r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Figure out box size around atom:
                box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1%map_size ) )

!               Convert to fractional:
                IF ( m > 1 ) THEN
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat)
                ELSE
                    centre_of_the_atom_in_fractions = map_1%DEORT        * pdb_2%xyz(iat)
                ENDIF

!               This is our atomic centre on the grid: THIS CAN BE OUT OF MAP RANGE:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size

!               Calculate gradient pointers:
                indxyz = pdb_2%indxyz(iat)
                refine_atomic_xyz  = pdb_2%refinable(1:3, iat)

                refine_atomic_biso = pdb_2%refinable(4, iat)

                IF ( refine_atomic_biso ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF

                refine_atomic_occ = pdb_2%refinable(5,iat)
                refine_atomic_Us  = pdb_2%refinable(6:11,iat)
                
!               
                IF ( refine_atomic_occ ) THEN
                    indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

                IF ( ANY ( refine_atomic_Us ) ) THEN
                    indu = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

                IF ( pdb_2%refinable(4,iat) .AND. ANY ( pdb_2%refinable(6:11,iat) ) ) THEN
                    WRITE(*,*) pdb_2%refinable(4,iat)
                    WRITE(*,*)  pdb_2%refinable(6:11,iat)
                    CALL die('Cannot refine Biso and Us for the same atom simultaneously.', srname)
                ENDIF

                IF ( debug > 1000 ) THEN
                    WRITE(*,*) ' indxyz=', indxyz, ' iat=',iat 
                ENDIF

!               Create a box around the atom and convolute density with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    
                    DO jv = box_lower_limit(2), box_upper_limit(2)

                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw  = (/ju, jv, jw/)

!                           This is the correct interpretatin of convolution theorem (r(i)-r):
                            duvw = centre_of_the_atom_in_fractions - REAL ( uvw ) / map_1%map_size

!                           High point of electron density generation:
                            r2 = map_1%real_tensor * duvw

                            r_cut:IF ( r2 <= r2_cut_off ) THEN

                                jut       = uvw .MMOD. map_1%map_size
                                map_value = map_1%array(jut(1), jut(2), jut(3))

!                               Convert to orthogonal system of coodinates (easier to calculate derivatives):
                                IF ( m == 1 ) THEN
                                    duvwort = map_1%ORT * duvw
                                ELSE

!                                   This is inverted matrix of a sym operator multiplied by orthogonalization matrix:
                                    duvwort = map_1%SYM_ORT_INVERTED(m) * duvw
                                ENDIF



!                               Calculate appropriate exponenets:
                                IF ( Us_present ) THEN

!                                   Going to speed up aniso calculations:
                                    DO l = 1, ngauss
                                        VTS_duvwort(1:3,l) = VTS(l) * duvwort
                                        aexp_qform(l) = a(l) * EXP ( -0.5_wp * DOT_PRODUCT ( duvwort, VTS_duvwort(1:3,l) ) )
                                    ENDDO

                                ELSE

!                                   Going to speed up in case of Biso/occncy refinement:
                                    exp_be(1:ngauss) = EXP ( -be(1:ngauss) * r2 )

                                ENDIF

                                IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                   Accumulate gradient:
                                    IF ( .NOT. Us_present ) THEN

!                                       If some coordinates unrefinable do it one by one:
                                        l = 0
                                        DO k = 1, 3
                                            IF ( refine_atomic_xyz(k) ) THEN
                                                l = l + 1
                                                g(indxyz+l) = g(indxyz+l) - 2.0_wp * map_value * occ &
                                                            * SUM ( be(1:ngauss) * a(1:ngauss) * exp_be (1:ngauss) )       &
                                                            * duvwort(k)
                                            ENDIF
                                         ENDDO

                                    ELSE

!                                       Us_present -> do appropiate aniso calculations even assuming Us are not refined:
                                        dRodX = -MATMUL ( VTS_duvwort(1:3,1:ngauss), aexp_qform(1:ngauss) ) 

!                                       Accumulate refined components only:
                                        l = 0

!                                       Convert vector to array:
                                        DO k = 1, 3
                                           IF ( refine_atomic_xyz(k) ) THEN
                                               l = l + 1
                                               g(indxyz+l) = g(indxyz+l) + map_value * occ * dRodX(k)
                                           ENDIF
                                        ENDDO

                                    ENDIF

                                ENDIF
                                
                                aniso_or_biso:IF ( refine_atomic_biso ) THEN

!                                   Accumulate derivatives:
                                    g(indbiso) = g(indbiso) + occ * map_value * SUM ( bb(1:ngauss) * a(1:ngauss) &
                                               * ( -1.5_wp + be(1:ngauss) * r2 ) * exp_be(1:ngauss)  ) 

                                ELSE IF ( ANY ( refine_atomic_Us ) ) THEN

!                                   Prepare xxt vector:
                                    xxt(1) = duvwort(1) ** 2
                                    xxt(2) = duvwort(2) ** 2
                                    xxt(3) = duvwort(3) ** 2
                                    xxt(4) = 2.0_wp * duvwort(1) * duvwort(2)
                                    xxt(5) = 2.0_wp * duvwort(1) * duvwort(3)
                                    xxt(6) = 2.0_wp * duvwort(2) * duvwort(3)
                                   
!                                   Loop over all gaussians to get dRo/dU derivs:
                                    DO l = 1, ngauss
                                        utemp(1:6,l) = v(1:6,l) + MATMUL ( VTSD(1:6,1:6,l), xxt )
                                    ENDDO
                                    dRodU = -0.5_wp * MATMUL ( utemp(1:6,1:ngauss), aexp_qform(1:ngauss) )

                                    l = 0
                                    DO k = 1, 6
                                        IF ( refine_atomic_Us(k) ) THEN
                                            l = l + 1
                                            g(indu + l) = g(indu + l) + map_value * occ * dRodU(k)
                                        ENDIF
                                    ENDDO

                                ENDIF aniso_or_biso

!                               "Varying blocks" scheme:
                                IF ( refine_atomic_occ ) THEN

!                                   Accumulate derivatives for occupancy:
                                    IF ( Us_present ) THEN
                                        g(indocc) = g(indocc) + map_value * SUM ( aexp_qform(1:ngauss) )
                                    ELSE
                                        g(indocc) = g(indocc) + map_value * SUM (  a(1:ngauss) * exp_be(1:ngauss) )
                                    ENDIF
                                ENDIF

                            ENDIF r_cut

                        ENDDO
                    ENDDO
                ENDDO boxz

            ENDDO coord
        ENDDO symmetry

    END SUBROUTINE first_derivatives_wrt_xyzbqu

    SUBROUTINE matrix_vector_convolution(g, map_1, pdb_2, true_sp_group)
!
!       Purpose:
!       =======
!       Calculates matrix vector product using array of maps.
!
!       Just convolutes modified atomic density using Tronrud's method with Fourier map using `ngauss' Gaussian approximation
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2007       B.Strokopytov          Original code
!
!       Oct 2007       FIXME                  Needs some research to figure symmetry of resulting maps.
!                                             It will definitely work in P1 sp.group.
!                                           
!      VERY IMPORTANT NOTE: TO KEEP MATRIX SYMMETRICAL R2_CUT_OFF SHOULD CHANGE IN SYNC
!                           WITH RADIUS IN GENDEN
!
        REAL(KIND=wp), DIMENSION(:),              INTENT(INOUT)        :: g
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: map_1
        TYPE(pdb),                                INTENT(IN)           :: pdb_2
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
        REAL(KIND=wp), DIMENSION(5)                                    :: exp_be 
        REAL(KIND=wp), DIMENSION(5)                                    :: abexp 
        REAL(KIND=wp)                                                  :: b_scale
        REAL(KIND=wp)                                                  :: gaussian_sum 
        REAL(KIND=wp), DIMENSION(6)                                    :: gg
!       Effective atomic radii:
        REAL(KIND=wp)                                                  :: r_cut_off
        REAL(KIND=wp)                                                  :: r2_cut_off
        REAL(KIND=wp)                                                  :: r2
        TYPE(vector_int)                                               :: uvw
        TYPE(vector)                                                   :: duvw
!        REAL(KIND=wp), DIMENSION(3)                                    :: duvwort
!       Can handle 11 maps only:
        REAL(KIND=wp), DIMENSION(11)                                   :: map_value
        REAL(KIND=wp)                                                  :: occ
!       Atom/box parameters:
        TYPE(vector)                                                   :: centre_of_the_atom_in_fractions
        TYPE(vector_int)                                               :: closest_grid_point
        TYPE(vector_int)                                               :: box_size
        INTEGER, DIMENSION(3)                                          :: jut
        INTEGER, DIMENSION(3)                                          :: box_lower_limit
        INTEGER, DIMENSION(3)                                          :: box_upper_limit
        TYPE(vector)                                                   :: deort_diag
        TYPE(vector)                                                   :: ort_diag
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

!       Refinement variables:
        LOGICAL, DIMENSION(3)                                          :: refine_atomic_xyz
        LOGICAL                                                        :: refine_atomic_Biso
        LOGICAL                                                        :: refine_atomic_occ
        LOGICAL, DIMENSION(6)                                          :: refine_atomic_Us
        INTEGER                                                        :: indxyz
        INTEGER                                                        :: indbiso
        INTEGER                                                        :: indocc
        INTEGER                                                        :: indu
        INTEGER                                                        :: number_of_maps
        INTEGER                                                        :: my_first_allocated_map
!       Counters:
        INTEGER                                                        :: iat
        INTEGER                                                        :: imap
        INTEGER                                                        :: l
        INTEGER                                                        :: m
!       Box counters:
        INTEGER                                                        :: ju
        INTEGER                                                        :: jv
        INTEGER                                                        :: jw
        CHARACTER(LEN=32),                                        SAVE :: srname = 'matrix_vector_convolution'

!       Start with basic checks:
        number_of_maps = SIZE ( map_1 )

!       Basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname)
        ENDIF

!       Another important check:
        IF ( .NOT. ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL die('Programming error. PDB_2%REFINABLE has not been initialized properly.',&
                     srname)
        ELSE
            IF ( SIZE ( g ) /= COUNT ( pdb_2%refinable ) ) THEN
                CALL die('Programming error. Inconsistent sizes for G and PDB_2%REFINABLE arrays.',&
                srname)
            ENDIF
        ENDIF

        aniso_active = ALLOCATED ( pdb_2%U )

!       Initialize parameters:
        g = 0.0_wp
        my_first_allocated_map = first_allocated_map ( map_1 )

!       Copy b_scale from the first allocated map:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1(my_first_allocated_map)%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1(my_first_allocated_map)%cell
            CALL die('Programming error. Unit cell has not been initialized for current map.',&
                     srname)
        ENDIF

!       Calculate maximum sphere radius around any atom:
        ort_diag = extract_diagonal ( map_1(my_first_allocated_map)%ORT )

!       Make sure the whole sphere fits in unit cell in any case:
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
        ELSE
            nsym = 1
        ENDIF

        symmetry:DO m = 1, nsym

            indxyz = 0
            coord:DO iat = 1, natom

                occ = pdb_2%occ(iat)

!               Save type of atom:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss = pdb_2%atomsf_table%record(ityp)%ngauss 

!               Check:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             'generate_density')
                ENDIF

!               Copy gaussian constants:
                a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)
                Biso = pdb_2%biso(iat)

!               Check whether Us have been read (appropriate record present):
                Us_present = .FALSE.

                IF ( aniso_active ) THEN

!                   Note that ">" was overloaded in a special way:
                    Us_present = ANY ( pdb_2%U(iat) > 0.0_wp )

!                   Convert to matrix of aniso elements:
                    UMAT = pdb_2%U(iat)

!                   Calculate eigenvalues:
                    ueigen = eigen_values ( UMAT )
                    Biso = 8.0_wp * pi_squared * MAXVAL ( ueigen )

                ENDIF


!               Finally we have sensible approach:
!                l = MAXLOC ( ABS ( a(1:ngauss) ), DIM=1 )

!               Superior quality radius for atoms with high B-values, gives excellent result on RNAse:
                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + Biso + b_scale ) ** 2              

!               Make sure radius can only increase in case of high B-values:
                r2_cut_off = MAX ( map_1(my_first_allocated_map)%fft_radius ** 2 , r2_cut_off ) 
!               Make sure sphere around any atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

                IF ( .NOT. Us_present ) THEN
                    a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                    be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )
                ELSE

                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian:
                        U(l) = pdb_2%U(iat) + ( ( b(l) + b_scale ) / ( 8.0_wp * pi_squared ) )

!                       Transform ellipsoid using Garib's method:
!                        IF ( m /= 1 ) THEN
!                            U(l) = RU(m) * U(l)
!                        ENDIF

!                       Calculate determinant:
                        Udet(l)  = udet_aniso(U(l))

!                       Calculate inverted matrix:
                        v = minors_aniso(U(l)) / Udet(l)
                        VTS(l) = .CONVERT. v
                    ENDDO

!                   Scale constants in advance to reduce computation:
                    a(1:ngauss) = a(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF

                IF ( debug > 25 ) THEN
                    WRITE(*,*) ' potential matrix vector r_cut_off=', SQRT(r2_cut_off)
                ENDIF

                IF ( debug == 22 ) THEN
                    WRITE(*,*) ' matrix vector r_cut_off=', r_cut_off
                ENDIF

!               Figure out box size around atom:
                box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

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

                refine_atomic_xyz = pdb_2%refinable(1:3,iat)
                refine_atomic_biso = pdb_2%refinable(4,iat)
                refine_atomic_occ = pdb_2%refinable(5,iat)
                refine_atomic_Us  = pdb_2%refinable(6:11,iat)

!               Calculate matrix pointers:
                IF ( refine_atomic_biso ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF

                IF ( refine_atomic_occ ) THEN
                    indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

                IF ( ANY ( refine_atomic_Us ) ) THEN
                    indu = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF                

!               Special debugging:
                IF ( debug > 22 ) THEN
                    WRITE(*,*) iat, indxyz, refine_atomic_xyz
                ENDIF

!               Create a box around the atom and convolute density with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    
                    DO jv = box_lower_limit(2), box_upper_limit(2)

                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_the_atom_in_fractions

!                           High point of electron density generation:
                            r2 = map_1(my_first_allocated_map)%real_tensor * duvw

                            IF ( r2 <= r2_cut_off ) THEN

                                jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                IF ( Us_present ) THEN
                                    
                                    duvwort = map_1(my_first_allocated_map)%ORT * duvw
                                    DO l = 1, ngauss
                                        qform(l) = duvwort .DOT. ( VTS(l) * duvwort )
                                    ENDDO

                                    gaussian_sum = SUM ( a(1:ngauss) * EXP ( -0.5_wp * qform(1:ngauss) ) )

                                ELSE
!                                   Going to speed up in case of Biso/occncy refinement:
                                    exp_be(1:ngauss) = EXP ( -be(1:ngauss) * r2 )
                                    abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)
                                    gaussian_sum     = SUM ( abexp(1:ngauss) )
                                ENDIF

                                IF ( ANY ( refine_atomic_xyz ) ) THEN                                    
                                    DO l = 1, 3
                                        IF ( refine_atomic_xyz(l) ) THEN
                                            map_value(l) = map_1(l)%array(jut(1), jut(2), jut(3))
                                        ENDIF
                                    ENDDO
                                ENDIF
 
                                IF ( refine_atomic_biso ) THEN
                                    map_value(4) = map_1(4)%array(jut(1), jut(2), jut(3)) 
                                ENDIF

                                IF ( refine_atomic_occ ) THEN
                                    map_value(5) = map_1(5)%array(jut(1), jut(2), jut(3)) 
                                ENDIF

                                IF ( ANY ( refine_atomic_Us ) ) THEN
                                    DO l = 6, 11
                                        map_value(l) = map_1(l)%array(jut(1), jut(2), jut(3))
                                    ENDDO
                                ENDIF

!                               See if we want XYZ refinement:
                                IF ( ANY ( refine_atomic_xyz ) ) THEN

                                    gg(1:3) = 0.0_wp
                                    DO imap = 1, 3
                                        IF ( refine_atomic_xyz(imap) ) THEN
                                            gg(imap) =  occ * map_value(imap) * gaussian_sum
                                        ENDIF
                                    ENDDO

!                                   Apply translationless inverted symmetry operator:
                                    IF ( m > 1 ) THEN
!                                       duvwort = map_1(1)%DEORT * gg
                                        gg(1:3) = map_1(my_first_allocated_map)%SYM_ORT_INVERTED_DEORT(m) * gg(1:3)
                                    ENDIF

!                                   Initialise pointer:
                                    l = 0
                                    DO imap = 1, 3
                                        IF ( refine_atomic_xyz(imap) ) THEN
                                            l = l + 1
                                            g(indxyz+l) = g(indxyz+l) + gg(imap)
                                        ENDIF
                                    ENDDO
                                ENDIF

                                IF ( refine_atomic_biso ) THEN

!                                   Accumulate matrix-vector using map #4:
                                    g(indbiso) = g(indbiso) + occ * map_value(4) * gaussian_sum
                                    
                                ENDIF

                                IF ( refine_atomic_occ ) THEN

!                                   Accumulate matrix_vector product using map #5 without atomic occncy:
                                    g(indocc) = g(indocc) + map_value(5) * gaussian_sum
                                         
                                ENDIF

                                IF ( ANY ( refine_atomic_Us ) ) THEN

                                    gg(1:6) = 0.0_wp

                                    DO imap = 1, 6
                                        IF ( refine_atomic_Us(imap) ) THEN
                                            gg(imap) = occ * map_value(5+imap) * gaussian_sum
                                        ENDIF
                                    ENDDO

!                                   Hope this will never be used:
                                    IF ( m > 1 ) THEN
                                        gg = TRANSPOSE ( RU(m) ) * gg
                                    ENDIF

                                    l = 0
                                    DO imap = 1, 6
                                        IF ( refine_atomic_Us(imap) ) THEN
                                            l = l + 1
                                            g(indu+l) = g(indu+l) + gg(imap)
                                        ENDIF
                                    ENDDO
                                ENDIF
                            ENDIF

                        ENDDO
                    ENDDO
                ENDDO boxz

!               Increment pointer to atomic block:
                indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )

            ENDDO coord
        ENDDO symmetry

    END SUBROUTINE matrix_vector_convolution 

    SUBROUTINE fast_matrix_vector_convolution(g, map_1, pdb_2, all_atom_images, true_sp_group)
!
!       Purpose:
!       =======
!       This modification of MATRIX_VECTOR_CONVOLUTION routine calculates matrix vector product for array of maps using
!       precalculated atom density images which results in approximately 10x speed-up on our DDR2 architecture (Intel 
!       945GM chip).
!
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2008       B.Strokopytov          Original code
!
!
        REAL(KIND=wp),    DIMENSION(:),              INTENT(INOUT)        :: g
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: map_1
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        TYPE(atom_image), DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: all_atom_images
        TYPE(space_group),                           INTENT(IN), OPTIONAL :: true_sp_group
!       Local variables:
        INTEGER                                                           :: natom
        INTEGER                                                           :: nsym
        REAL(KIND=wp)                                                     :: occ
!       Refinement variables:
        INTEGER                                                           :: indxyz
        INTEGER                                                           :: number_of_maps
        INTEGER                                                           :: my_first_allocated_map
!       Counters:
        INTEGER                                                           :: iat
        INTEGER                                                           :: imap
        INTEGER                                                           :: l
        INTEGER                                                           :: m
        INTEGER                                                           :: n
        CHARACTER(LEN=32),                                       SAVE     :: srname = 'fast_matrix_vector_convolution'
        REAL(KIND=wp), EXTERNAL                                           :: convolution_no_1

!       Start with basic checks:
        number_of_maps = SIZE ( map_1 )

!       Basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname)
        ENDIF

!       Another important check:
        IF ( .NOT. ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL die('Programming error. PDB_2%REFINABLE has not been initialized properly.',&
                     srname)
        ELSE
            IF ( SIZE ( g ) /= COUNT ( pdb_2%refinable ) ) THEN
                CALL die('Programming error. Inconsistent sizes for G and PDB_2%REFINABLE arrays.',&
                srname)
            ENDIF
        ENDIF

!       Initialize parameters:
        g = 0.0_wp
        my_first_allocated_map = first_allocated_map ( map_1 )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1(my_first_allocated_map)%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1(my_first_allocated_map)%cell
            CALL die('Programming error. Unit cell has not been initialized for current map.',&
                     srname)
        ENDIF
!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
        ELSE
            nsym = 1
        ENDIF

!       This is obsolete:
        symmetry:DO m = 1, nsym

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED ( PDB_2, ALL_ATOM_IMAGES, MAP_1, G, &
!$OMP                                       NATOM, NUMBER_OF_MAPS, DEBUG, M, NSYM )
            coord:DO iat = 1, natom

                occ = pdb_2%occ(iat)

                indxyz = pdb_2%indxyz(iat)
!               Special debugging:
                IF ( debug > 22 ) THEN
                    WRITE(*,*) iat, indxyz
                ENDIF

                n = SIZE ( all_atom_images(iat)%ro  )

!               Initialise pointer:
                l = 0
                DO imap = 1, number_of_maps
                    IF ( pdb_2%refinable(imap,iat) ) THEN
                        l = l + 1
                        IF ( imap /= 5 ) THEN
                            g(indxyz+l) = convolution_no_1(n, all_atom_images(iat)%ro, all_atom_images(iat)%uvw, &
                                                           occ, map_1(imap)%array)
                        ELSE
                            g(indxyz+l) = convolution_no_1(n, all_atom_images(iat)%ro, all_atom_images(iat)%uvw, &
                                                           1.0_wp, map_1(imap)%array)
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO coord
        ENDDO symmetry

    END SUBROUTINE fast_matrix_vector_convolution

    SUBROUTINE convolution_via_second_derivatives(g, map_1, pdb_2, mode, true_sp_group)
!
!       Purpose:
!       =======
!       Calculates matrix vector product using array of maps.
!       Just convolutes modified atomic density using Tronrud's method with appropriare Fourier map
!       achieving summation.
!
!       Note:
!       ====
!       This method seems numerically less accurate than algorithm which uses MATRIX_VECTOR_CONVOLUTION
!       routine above. Logic of this routine is much more complex.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2007       B.Strokopytov          Original code
!
!       Dec 2007       B.Strokopytov          Needs some work to improve speed.
!                                             Too many similar sums being calculated.
!     
!       Dec 2007       B.Strokopytov          Heavily tested. Some loss of precision is still observed.
!                                             This results in symmetry factors from SYMMLQ about 10^-7 - 10^-9
!                                             which is worse compared to MATRIX_VECTOR_CONVOLUTION
!                                             which gives 10^-14 - 10^-16 symmetry quality factors.
!
!                                             This imprecision affects smallest eigenvalues.
!                                            
!       Jan 2008                              BUG associated with missing global refinement parameters CORRECTED.
!                                             Subroutine which params has been added.
!
!       Jan 2008       B.Strokopytov          Some work on speed up has been done. SUM_ABEXP_BE has been added.
!                                             Requires further work on speed and precision (unless precision cannot
!                                             be improved due to asymmetry in the calculations). Note
!                                             also that we no longer have REAL joint atom. Just image
!                                             of single atom is being used.
!
!
        REAL(KIND=wp), DIMENSION(:),              INTENT(INOUT)        :: g
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: map_1
        TYPE(pdb),                                INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                         INTENT(IN)           :: mode
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
        REAL(KIND=wp), DIMENSION(5)                                    :: bbb
        REAL(KIND=wp), DIMENSION(5)                                    :: ber2
        REAL(KIND=wp), DIMENSION(5)                                    :: exp_be 
        REAL(KIND=wp), DIMENSION(5)                                    :: abexp 
        REAL(KIND=wp), DIMENSION(5)                                    :: abexp_be 
        REAL(KIND=wp)                                                  :: b_scale
!       Effective atomic radii:
        REAL(KIND=wp)                                                  :: r_cut_off
        REAL(KIND=wp)                                                  :: r2_cut_off
        REAL(KIND=wp)                                                  :: r2
        TYPE(vector_int)                                               :: uvw
        TYPE(vector)                                                   :: duvw
        REAL(KIND=wp), DIMENSION(3)                                    :: duvwort
        REAL(KIND=wp), DIMENSION(3)                                    :: gg
        REAL(KIND=wp), DIMENSION(3,3)                                  :: BMAT
        REAL(KIND=wp), DIMENSION(3,3)                                  :: GMAT
        REAL(KIND=wp), DIMENSION(3,3)                                  :: XXT
!       Can handle 11 maps only:
        REAL(KIND=wp), DIMENSION(11)                                   :: map_value
        REAL(KIND=wp)                                                  :: occ
        REAL(KIND=wp)                                                  :: gaussian_sum
        REAL(KIND=wp)                                                  :: sum_abexp_be
!       Atom/box parameters:
        TYPE(vector)                                                   :: centre_of_the_atom_in_fractions
        TYPE(vector_int)                                               :: closest_grid_point
        TYPE(vector_int)                                               :: box_size
        INTEGER, DIMENSION(3)                                          :: jut
        INTEGER, DIMENSION(3)                                          :: box_lower_limit
        INTEGER, DIMENSION(3)                                          :: box_upper_limit
        TYPE(vector)                                                   :: deort_diag
        TYPE(vector)                                                   :: ort_diag
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                                  :: r2_cell_max
!       Refinement variables:
        LOGICAL, DIMENSION(3)                                          :: refine_atomic_xyz
        LOGICAL                                                        :: refine_atomic_Biso
        LOGICAL                                                        :: refine_atomic_occ
        INTEGER                                                        :: indxyz
        INTEGER                                                        :: indbiso
        INTEGER                                                        :: indocc
        INTEGER                                                        :: my_first_allocated_map
!       Global refinement variables:
        LOGICAL                                                        :: refine_xyz
        LOGICAL                                                        :: refine_Biso
        LOGICAL                                                        :: refine_occ
        LOGICAL                                                        :: refine_Us
!       Counters:
        INTEGER                                                        :: iat
        INTEGER                                                        :: l
        INTEGER                                                        :: k 
        INTEGER                                                        :: m
!       Box counters:
        INTEGER                                                        :: ju
        INTEGER                                                        :: jv
        INTEGER                                                        :: jw
        CHARACTER(LEN=32),                                        SAVE :: srname='convolution_via_second_derivatives'

!       Basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname)
        ENDIF

!       Another important check:
        IF ( .NOT. ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL die('Programming error. PDB_2%REFINABLE has not been initialized properly.',&
                     srname)
        ELSE
            IF ( SIZE ( g ) /= COUNT ( pdb_2%refinable ) ) THEN
                CALL die('Programming error. Inconsistent sizes for G and PDB_2%REFINABLE arrays.',&
                srname)
            ENDIF
        ENDIF

        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

!       Set this from the start:
        my_first_allocated_map = first_allocated_map ( map_1 )

!       Initialize parameters:
        g = 0.0_wp

!       GMAT:
        GMAT = 0.0_wp
        DO k = 1, 3
            GMAT(k,k) = 2.0_wp
        ENDDO

!       Copy b_scale from the first map:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1(my_first_allocated_map)%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1(my_first_allocated_map)%cell
            CALL die('Unit cell is too small for calculations. Use SHELXL.',&
                     srname)
        ENDIF

!       Calculate maximum sphere radius around any atom:
        ort_diag = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.25_wp * MINVAL ( ort_diag .X. ort_diag )

!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
        ELSE
            nsym = 1
        ENDIF

        symmetry:DO m = 1, nsym

!           Initialise pointer:
            indxyz = 0

!           Loop over all atoms:
            coord:DO iat = 1, natom

                occ = pdb_2%occ(iat)
                 
!               BUG CORRECTED this how it should be:
!                IF ( COUNT ( pdb_2%refinable(1:11, iat) ) == 0 ) THEN
!                    CYCLE
!                ENDIF

!               Save type of atom:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss =  pdb_2%atomsf_table%record(ityp)%ngauss 

!               Check:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             srname)
                ENDIF

!               Copy gaussian constants:
                a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)

!               Additional arrays:
                bbb(1:ngauss) = 1.0_wp / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

                r2_cut_off = integration_radius ( a(1:ngauss), b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ** 2

!               Make sure radius can only increase in case of high B-values:
                r2_cut_off = MAX ( map_1(my_first_allocated_map)%fft_radius ** 2 , r2_cut_off )

!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Debugging:
                IF ( debug > 40 ) THEN
                    WRITE(*,*) ' actual r_cut_off=',  SQRT ( r2_cut_off )
                ENDIF

!               Prepare atomic gaussians:
                a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

!               Figure out box size around atom:
                box_size  = NINT( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

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

!               Calculate matrix pointers:
                refine_atomic_xyz  = pdb_2%refinable(1:3,iat)
                refine_atomic_biso = pdb_2%refinable(4,iat)
                refine_atomic_occ =  pdb_2%refinable(5,iat)
            
                IF ( refine_atomic_biso ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF

                IF ( refine_atomic_occ ) THEN
                    indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

!               Create a box around the atom and convolute density with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    
                    DO jv = box_lower_limit(2), box_upper_limit(2)

                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_the_atom_in_fractions

!                           High point of electron density generation:
                            r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw

                            cut_off:IF ( r2 <= r2_cut_off ) THEN

                                jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

!                               Convert to orthogonal system of coodinates, need this for all derivatives:
                                IF ( refine_xyz ) THEN
                                    duvwort = 2.0_wp * map_1(my_first_allocated_map)%ORT * duvw
                                ENDIF

!                               Going to speed up in case of Biso/occncy refinement:
                                ber2(1:ngauss)   = be(1:ngauss) * r2
                                exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)

!                               Additional array for speed up:
                                abexp_be(1:ngauss) = abexp(1:ngauss) * be(1:ngauss)
                                sum_abexp_be = SUM ( abexp_be(1:ngauss) )

!                               Must use global definitions here: 
                                IF ( refine_xyz ) THEN

!                                   Even if XYZ(iat) is unrefinable we may need this for cross derivatives (e.g.in P1):
                                    DO l = 1, 3
                                        map_value(l) = map_1(l)%array(jut(1), jut(2), jut(3))
                                    ENDDO
                                ENDIF

                                IF ( refine_biso ) THEN
                                    map_value(4) = map_1(4)%array(jut(1), jut(2), jut(3))
                                ENDIF

                                IF ( refine_occ ) THEN
                                    map_value(5) = map_1(5)%array(jut(1), jut(2), jut(3))
                                ENDIF

                                anyxyz:IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                   Calculate XXT matrix:
                                    DO l = 1, 3
                                        DO k = 1, 3

!                                           Note that DUVWORT has been scaled by a factor of 2 already:
                                            XXT(l,k) = duvwort(k) * duvwort(l)
                                        ENDDO
                                    ENDDO

!                                   Accumulate matrix-vector product wrt XYZ(11):
                                    BMAT = sum_abexp_be * GMAT - SUM ( abexp_be(1:ngauss) * be(1:ngauss) ) * XXT

!                                   Accumulate product for coordinates via scaled BMAT:
                                    gg = occ * MATMUL ( BMAT, map_value(1:3) )
                                    CALL convert_gg(g, refine_atomic_xyz, gg, indxyz)

!                                   All contibutions from coordinates and their interactions has been added...

                                ENDIF anyxyz

!                               Atom may not have refinable biso but we allow calculation of cross terms:
                                my_refine_biso:IF ( refine_biso ) THEN

!                                   Need to use global refinement parameters here:
                                    IF ( refine_xyz ) THEN

                                        gaussian_sum = SUM ( abexp_be(1:ngauss) * bbb(1:ngauss) &
                                                     * (5.0_wp / 2.0_wp - ber2(1:ngauss)) ) 

!                                       Accumulate matrix-vector product wrt XB+YB+ZB(14) using map #4:
                                        IF ( ANY ( refine_atomic_xyz ) ) THEN
                                            gg = -occ * map_value(4) * gaussian_sum * duvwort
                                            CALL convert_gg(g, refine_atomic_xyz, gg, indxyz)
                                        ENDIF

!                                       --- Biso row start here:
                                        
!                                       Accumulate matrix-vector product wrt BX+BY+BZ (14) using map #1-3:
                                        IF ( refine_atomic_biso ) THEN
                                            g(indbiso) = g(indbiso) &
                                                       + occ * gaussian_sum * DOT_PRODUCT ( map_value(1:3), duvwort )
                                        ENDIF
                                    ENDIF

!                                   Continue with Biso row:  

!                                   Accumulate BB matrix-vector(12) if atom has refinable Biso:
                                    IF ( refine_atomic_biso ) THEN
                                        g(indbiso) = g(indbiso)                                                        &
                                                   + occ * map_value(4)                                                &
                                                   * SUM (  abexp(1:ngauss) * bbb(1:ngauss) ** 2                       &
                                                   * ( ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp) + 15.0_wp / 4.0_wp) )     

                                    ENDIF

                                ENDIF my_refine_biso

!                               Atom may not have refinable occupancy but we must allow calculation of cross terms
!                               therefore using global refinement definition here:
                                my_refine_occ:IF ( refine_occ ) THEN

!                                   Note: Maps #1-4 remove one occupancy (occ(jat)),
!                                   but map #5 does not remove any 

!                                   Finish accumulation for XYZ rows:

!                                   Accumulate derivatives XO, YO, ZO(15) for XYZ rows using map #5:
                                    IF ( ANY ( refine_atomic_xyz ) ) THEN
                                        gg = occ * map_value(5) * sum_abexp_be * duvwort
                                        CALL convert_gg(g, refine_atomic_xyz, gg, indxyz)
                                    ENDIF

!                                   Finish Biso row:

!                                   Accumulate derivatives BQ(16) for Biso row using map #5:
                                    IF ( refine_atomic_Biso ) THEN
                                        g(indbiso) = g(indbiso) + occ * map_value(5)        &     
                                                   * SUM ( abexp(1:ngauss) * bbb(1:ngauss)  &
                                                   *     ( -1.5_wp + ber2(1:ngauss) ) )
                                    ENDIF

!                                   --- Now occupancy row itself, note OCC is not used anymore ---

!                                   If we are here then occupancies are being refined but we need to check whether this
!                                   particular atom is being refined against occupancy:
                                    atomic_occ:IF ( refine_atomic_occ ) THEN

!                                       Accumulate OX+OY+OZ interaction (15) using global refinement parameter check:
                                        IF ( refine_xyz )  THEN
                                            g(indocc) = g(indocc) - sum_abexp_be * DOT_PRODUCT ( map_value(1:3), duvwort )
                                        ENDIF

!                                       Summation of matrix elements does not requir that this atom has refinable Biso,
!                                       using only global refinement definition:
                                        IF ( refine_Biso ) THEN

!                                           QB interaction(16), note plus sign here:    
                                            g(indocc) = g(indocc) + map_value(4)               &
                                                      * SUM ( abexp(1:ngauss) * bbb(1:ngauss)  &   
                                                      * (-1.5_wp + ber2(1:ngauss)) )
                                        ENDIF

!                                       Finally QQ interaction (13):
                                        g(indocc) = g(indocc) + map_value(5) * SUM ( abexp(1:ngauss) )

                                    ENDIF atomic_occ
     
                                ENDIF my_refine_occ

                            ENDIF cut_off

                        ENDDO
                    ENDDO
                ENDDO boxz

!               Increment pointer to atomic block:
                indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )

            ENDDO coord
        ENDDO symmetry

    END SUBROUTINE convolution_via_second_derivatives

    SUBROUTINE simple_convolution ( g, map_1, pdb_2, mode, true_sp_group )
!
!       Purpose:
!       =======
!       Gradient calculation. This sub convolutes atomic density derivatives
!       with Fourier map using `ngauss' Gaussian approximation.
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Oct 2003       B.Strokopytov          Original code
!
!       Apr 2005       B.Strokopytov          When calculating radius cut-off occupancies shall be ignored
!                                             from now on to eliminate negative effects associated with
!                                             small occupancies. We need this for certain algorithms.
!
!       Nov 2005       B.Strokopytov          e_limit has been changes to 18.0 A in constants.f90
!                                             instead of original value of 7.0 which allows
!                                             in conjuction with appropriate anti-alias b_scale
!                                             extremely accurate calculation of structure factors
!     
!       Oct 2007       B.Strokopytov          atomsf_table introduced ( memory economy )
!                                             plus flexible number of gaussians
!
!       Oct 2007       B.Strokopytov          J.Navaza's FFT radius finally added                   
!                                             e_limit tumbai.
!
!       Oct 2007       B.Strokopytov          Derivatives against Biso and occupancies added.                    
!
!       Oct 2007       B.Strokopytov          ort_diag added for calculation of maximum
!                                             inscribed sphere size
!
!       Oct 2007       B.Strokopytov          EXP_BE array added for speed up of Biso/occ refinement
!
!       Oct 2007       B.Strokopytov          More careful treatment of occupancies ( division by occ removed).
!
!       Oct 2007       B.Strokopytov          PDB_2%REFINABLE array added for generality and simplicity
!   
!
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT)        :: g
        TYPE(map),                   INTENT(IN)           :: map_1
        TYPE(pdb),                   INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),            INTENT(IN)           :: mode
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
        REAL(KIND=wp), DIMENSION(5)                       :: bb
        REAL(KIND=wp), DIMENSION(5)                       :: exp_be 
        REAL(KIND=wp)                                     :: b_scale
        REAL(KIND=wp)                                     :: r_cut_off
        REAL(KIND=wp)                                     :: r2_cut_off
        REAL(KIND=wp)                                     :: r2
        TYPE(vector_int)                                  :: uvw
        TYPE(vector)                                      :: duvw
        REAL(KIND=wp)                                     :: map_value
        REAL(KIND=wp)                                     :: occ
!       Atom/box parameters:
        TYPE(vector)                                      :: centre_of_the_atom_in_fractions
        TYPE(vector_int)                                  :: closest_grid_point
        TYPE(vector_int)                                  :: box_size
        INTEGER, DIMENSION(3)                             :: jut
        INTEGER, DIMENSION(3)                             :: box_lower_limit
        INTEGER, DIMENSION(3)                             :: box_upper_limit
        TYPE(vector)                                      :: deort_diag
        TYPE(vector)                                      :: ort_diag
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                     :: r2_cell_max
!       Refinement variables:
        LOGICAL                                           :: refine_xyz
        LOGICAL                                           :: refine_biso
        LOGICAL                                           :: refine_occ
        LOGICAL                                           :: refine_Us
        LOGICAL                                           :: refine_atomic_occ
        LOGICAL                                           :: refine_atomic_Biso
        LOGICAL, DIMENSION(3)                             :: refine_atomic_xyz
!       Pointers:
        INTEGER                                           :: indxyz
        INTEGER                                           :: indbiso
        INTEGER                                           :: indocc
!       Counters:
        INTEGER                                           :: iat
        INTEGER                                           :: k 
        INTEGER                                           :: l 
        INTEGER                                           :: m
!       Box counters:
        INTEGER                                           :: ju
        INTEGER                                           :: jv
        INTEGER                                           :: jw
        CHARACTER(LEN=32), SAVE                           :: srname = 'simple_convolution'

        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

!       Start with basic checks:
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE( pdb_2 )
        ELSE
            CALL die('Programming error. PDB_2 has not been initialized properly.',&
                     srname)
        ENDIF

        IF ( .NOT. ALLOCATED ( pdb_2%refinable ) ) THEN
            CALL die('Programming error. PDB_2%REFINABLE has not been initialized properly.',&
                     srname)            
        ELSE
            IF ( SIZE ( g ) /= COUNT ( pdb_2%refinable ) ) THEN
                CALL die('Programming error. Inconsistent sizes for G and PDB_2%REFINABLE arrays.',&
                srname)
            ENDIF
        ENDIF

!       Initialize parameters:
        g          = 0.0_wp
        b_scale    = map_1%b_scale
        deort_diag = extract_diagonal ( map_1%DEORT )

!       Check whether sensible cell parameters have been defined:
        IF ( ALL ( map_1%cell(1:3) < 2.0_wp ) ) THEN
            WRITE(*,*) map_1%cell
            CALL die('Programming error. Unit cell has not been initialized for current map.',&
                     srname)
        ENDIF

!       Calculate maximum sphere radius around any atom:
        ort_diag = extract_diagonal ( map_1%ORT )
        r2_cell_max = 0.25_wp * MINVAL ( ort_diag .X. ort_diag )

!       Decide about symmetry:
        IF ( PRESENT ( true_sp_group ) ) THEN
            nsym = true_sp_group%number_of_symops
        ELSE
            nsym = 1
        ENDIF

        symmetry:DO m = 1, nsym
            indxyz = 0
            coord:DO iat = 1, natom

!               Need this because of occupancy refinement:
                occ = pdb_2%occ(iat)
                IF ( occ == 0.0_wp .AND. .NOT. refine_occ ) THEN

!                   If we skip an atom we still need to increment the pointer:
                    indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )
                    CYCLE
                ENDIF

!               Save type of atom:
                ityp = pdb_2%atom_type(iat)

!               Calculate number of gaussians in this entry of table:
                ngauss =  pdb_2%atomsf_table%record(ityp)%ngauss 

!               Check:
                IF ( ngauss /= 2 .AND. ngauss /= 5 ) THEN
                    WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                             'generate_density')
                ENDIF

!               Copy gaussian constants:
                a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
                b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)

!               This is for Biso refinement only:
                IF ( refine_biso ) THEN
                    bb(1:ngauss) = 1.0_wp / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )
                ENDIF

                a(1:ngauss)  = a(1:ngauss) * ( SQRT ( pi4 / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale ) ) ) ** 3
                be(1:ngauss) = 4.0_wp * pi_squared        / ( b(1:ngauss) + pdb_2%biso(iat) + b_scale )

!               Calculate floating r2_cut_off:
                r2_cut_off = MAXVAL ( map_1%fft_radius ** 2 &
                           + LOG ( MAX ( ABS ( a(1:ngauss) ), small ) ) / be(1:ngauss) )
                r2_cut_off = MAX ( 0.0_wp, r2_cut_off )


!               Make sure sphere around atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Figure out box size around atom:
                box_size  = NINT ( r_cut_off * ( deort_diag .X. map_1%map_size ) )

!               Convert to fractional:
                IF ( m > 1 ) THEN
                    centre_of_the_atom_in_fractions = map_1%SYM_DEORT(m) * pdb_2%xyz(iat)
                ELSE
                    centre_of_the_atom_in_fractions = map_1%DEORT        * pdb_2%xyz(iat)
                ENDIF

!               This is our atomic centre on the grid: THIS CAN BE OUT OF MAP RANGE:
                closest_grid_point = NINT ( centre_of_the_atom_in_fractions .X. map_1%map_size )
                box_lower_limit    = closest_grid_point - box_size
                box_upper_limit    = closest_grid_point + box_size

!               Calculate gradient pointers:
                refine_atomic_xyz  = pdb_2%refinable(1:3, iat)

                refine_atomic_biso = refine_biso .AND. pdb_2%refinable(4, iat)

                IF ( refine_atomic_biso ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF

                refine_atomic_occ = refine_occ .AND. pdb_2%refinable(5,iat)

!               
                IF ( refine_atomic_occ ) THEN
                    indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

!               Create a box around the atom and convolute density with atom e.d. values:
                boxz:DO jw = box_lower_limit(3), box_upper_limit(3)
                    
                    DO jv = box_lower_limit(2), box_upper_limit(2)

                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw  = (/ju, jv, jw/)
                            duvw = REAL ( uvw ) / map_1%map_size - centre_of_the_atom_in_fractions

!                           High point of electron density generation:
                            r2 = map_1%real_tensor * duvw

                            IF ( r2 <= r2_cut_off ) THEN

                                jut       = uvw .MMOD. map_1%map_size
                                map_value = map_1%array(jut(1), jut(2), jut(3))

!                               Going to speed up in case of Biso/occncy refinement:
                                exp_be(1:ngauss) = EXP ( -be(1:ngauss) * r2 )

                                IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                   Accumulate gradient:
                                        l = 0
                                        DO k = 1, 3
! FIXME JUST FOR TEST: WANT ZEROS for g(2),g(3)
                                            IF ( k > 1 ) CYCLE
                                            IF ( refine_atomic_xyz(k) ) THEN
                                                l = l + 1
                                                g(indxyz+l) = g(indxyz+l) + 2.0_wp * map_value * occ                 &
                                                            * SUM (  a(1:ngauss) * exp_be (1:ngauss) )
                                            ENDIF
                                        ENDDO    
                                ENDIF

                                IF ( refine_atomic_biso ) THEN

!                                   Accumulate derivatives:
                                    g(indbiso) = g(indbiso) + occ * map_value * SUM ( bb(1:ngauss) * a(1:ngauss) &
                                               * ( -1.5_wp + be(1:ngauss) * r2 ) * exp_be(1:ngauss)  ) 
                                    
                                ENDIF

!                               "Varying blocks" scheme:
                                IF ( refine_atomic_occ ) THEN

!                                   Accumulate derivatives for occupancy:
                                    g(indocc) = g(indocc) + map_value * SUM (  a(1:ngauss) * exp_be(1:ngauss) )
                                ENDIF
                            ENDIF

                        ENDDO
                    ENDDO
                ENDDO boxz

!               Calculate pointers to gradient array:
                indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )

            ENDDO coord
        ENDDO symmetry

    END SUBROUTINE simple_convolution 

    SUBROUTINE normal_matrix_diagonal(g, map_1, pdb_2)
!
!       Purpose:
!       =======
!
!       Calculates diagonal elements G of normal matrix (H1 + H2 terms) with respect to
!       X,Y,Z,B and Q params.
!
!       Supply two maps a la Tronrud and you are done. Just one map will result in
!       calculation of H1 terms only.
!
!       Date:                 Programmer:             History of changes:
!       ====                  ==========              ==================
!       October  31, 2007     B.Strokopytov           Original code
!
!       November 14, 2007     B.Strokopytov           Accurate estimation of integration radius added.
!
!       November 27, 2007     B.Strokopytov           Occupancy only refinement speed up.
!
!       December 13, 2007     B.Strokopytov           Computation of H2 terms finally corrected. Code
!                                                     partially rewritten and reduced.
!
!       December 14, 2007     B.Strokopytov           Corrected chain rule for non-P1 sp.groups.
!
!       December 20, 2007     B.Strokopytov           Got rid of MTZ_1.
!     
!       December 20, 2007     B.Strokopytov           Added some notes on chain rule.
!
!       Theoretical note about chain rule:
!       =================================
!       If we have two block matrices (each block is 3x3):

!       |A11 A12 A13 A14|   |B11|
!       |A21 A22 A23 A24|   |B21|   
!       |A31 A32 A33 A34| & |B31|
!       |A41 A42 A43 A44|   |B41|
!
!       Then their product will be:
!
!       |A11 * B11 + A12 * B21 + A13 * B31 + A14 * B41|
!       !A21 * B11 + A22 * B21 + A23 * B31 + A24 * B41|
!       |A31 * B11 + A32 * B21 + A33 * B31 + A34 * B41|
!       |A41 * B11 + A42 * B21 + A43 * B31 + A44 * B41|
!
!       Chain rule note:
!       ===============
!       In orthogonal coordinate system we have:
!
!                     T      -1    -1 -1 -1       -1
!       [Asym] = [OAD] = [OAD] = [D  A  O  ] = [OA  D]
!
!       Obviously:
!             T
!       [Asym] = [OAD] 
!       
!       Transpose is due to chain rule, if we going to multiply from the right:
!         2    2      T
!       [d R/dx][Asym]
!
!       Remaining problems:
!       ==================
!       How to get rid of the symmetry loop? Looks like a horrible problem
!       in cubic space groups. 
!
        REAL(KIND=wp), DIMENSION(:),              INTENT(INOUT) :: g
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: map_1
        TYPE(pdb),                                INTENT(IN)    :: pdb_2
!       Local variables:
        LOGICAL                                                 :: h2
        TYPE(vector)                                            :: centre_of_atom_in_fractions
        TYPE(vector_int)                                        :: closest_grid_point
        TYPE(vector_int)                                        :: box_size
        INTEGER,       DIMENSION(3)                             :: jut
        INTEGER,       DIMENSION(3)                             :: box_lower_limit
        INTEGER,       DIMENSION(3)                             :: box_upper_limit
        TYPE(vector)                                            :: deort_diag
        TYPE(vector)                                            :: ort_diag
!       Normal atom params:
        INTEGER                                                 :: ityp
        INTEGER                                                 :: jtyp
        INTEGER                                                 :: ngauss
        INTEGER                                                 :: jgauss
        REAL(KIND=wp)                                           :: b_scale
!       Scattering:
!       For diagonal terms 15 gaussians would be the MAX number,
!       but for non-diagonal terms need 25 (looking into the future):
        INTEGER,       PARAMETER                                :: MAXGAU = 25
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: a
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: b
!       Joint atom:
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: ajo 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bjo
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: be
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bbb 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bbb2 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: ber2
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: exp_be
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: abexp
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: abexp_be
        REAL(KIND=wp)                                           :: gaussian_scale
!       Distance:
        REAL(KIND=wp)                                           :: r2
        TYPE(vector_int)                                        :: uvw
        TYPE(vector)                                            :: duvw
        REAL(KIND=wp), DIMENSION(3)                             :: duvwort
!       Cut-off:
        REAL(KIND=wp)                                           :: r2_cut_off
        REAL(KIND=wp)                                           :: r_cut_off
        REAL(KIND=wp)                                           :: max_r_cut_off
        REAL(KIND=wp)                                           :: min_r_cut_off
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                           :: r2_cell_max
!       Can handle 2 maps only:
        REAL(KIND=wp)                                           :: map_value
        REAL(KIND=wp)                                           :: occ
        REAL(KIND=wp)                                           :: occ_i
        REAL(KIND=wp)                                           :: occ_j
!       Maps variables:
        INTEGER                                                 :: number_of_maps
        INTEGER                                                 :: my_first_allocated_map
!       Refinement:        
        INTEGER                                                 :: np
        INTEGER                                                 :: natom
        LOGICAL,       DIMENSION(3)                             :: refine_atomic_xyz
        LOGICAL                                                 :: refine_atomic_biso
        LOGICAL                                                 :: refine_atomic_occ
!       Pointers:
        INTEGER                                                 :: indxyz
        INTEGER                                                 :: indbiso
        INTEGER                                                 :: indocc
!       Counters:
        INTEGER                                                 :: iat
        INTEGER                                                 :: jat
        INTEGER                                                 :: k 
        INTEGER                                                 :: l
        INTEGER                                                 :: run
!       Symmetry:
        TYPE(vector)                                            :: xyz_iat
        TYPE(vector)                                            :: xyz_jat
        TYPE(vector),  DIMENSION(2)                             :: xyz_frac
        REAL(KIND=wp), DIMENSION(3,3)                           :: GMAT
        REAL(KIND=wp), DIMENSION(3,3)                           :: XXT
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE              :: AMAT
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE              :: ASYM
        REAL(KIND=wp), DIMENSION(3,3)                           :: BMAT
!       Box counters:
        INTEGER                                                 :: i
        INTEGER                                                 :: j
!       Box counters:
        INTEGER                                                 :: ju
        INTEGER                                                 :: jv
        INTEGER                                                 :: jw
!       CPU monitoring:
        REAL(KIND=sp)                                           :: time0
        REAL(KIND=sp)                                           :: time1
!       Printing:
        INTEGER(KIND=eb)                                        :: modh
        INTEGER(KIND=eb)                                        :: pair
        CHARACTER(LEN=32),                               SAVE   :: srname = 'normal_matrix_diagonal'
!       Test:
        REAL(KIND=wp), DIMENSION(3)                             :: xyz

!       Checkz: 
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die('Programming error. pdb_2 has not been initialized.', &
                     srname)
        ENDIF

!       This must be set from the start:
        my_first_allocated_map = 1

!       Check map allocation:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Programming error. Map array has not been allocated properly.', &
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( map_1(my_first_allocated_map) ) ) THEN
            CALL die('Programming error. MAP_1(1) has not been allocated properly.',&
                     srname)
        ELSE

            number_of_maps = SIZE ( map_1 )
            IF ( number_of_maps < 2 ) THEN
                WRITE(*,*) number_of_maps
                CALL die('Programming error. Need 2 maps....', srname)
            ENDIF

        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' NORMAL_MATRIX_DIAGONAL> ', 'natom=', I6)") natom
            CALL die ( 'Programming error. Need at least one atom for density generation.', srname)
        ENDIF

!       Long and hard check: 
        np = COUNT ( pdb_2%refinable )
        IF ( np /=  SIZE ( g ) ) THEN
            WRITE(*,*) np, SIZE ( g )
            CALL die('Programming error. Incorrect size of diag G array.', srname)
        ENDIF
       
        CALL messag('Going to calculate '//TRIM ( int_to_c ( SIZE ( g ) ) )//' matrix elements... Wait...',&
                    srname)

!       Initialize diag array:
        g = 0.0_wp

!       Allocate additional matrices for chain rule:
        CALL allocate_array(ASYM, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops, 3)
        CALL allocate_array(AMAT, 3, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)

!       Build 3Nsym x 3 matrix applying chain rule (see notes above):
        DO j = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops

!           Calculate position of the block:  
            l = 3 * (j - 1)
!                                        T 
!           Calculating transpose [Asym ] :
            BMAT = map_1(my_first_allocated_map)%ORT             &
                 * .SYMA. map_1(my_first_allocated_map)%sp_group%SYM(j) &
                 *        map_1(my_first_allocated_map)%DEORT

!           To avoid Intel compiler warnings in debugging mode:
            ASYM(l+1:l+3,1:3) = BMAT
        ENDDO

!       Initialize 3x3 diagonal matrix:
        GMAT = 0.0_wp
        DO l = 1, 3
            GMAT(l,l) = 2.0_wp
        ENDDO

!       Other params:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Calculate maximum sphere radius around any atom:
        ort_diag    = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Set for statistics:
        max_r_cut_off = -10.0_wp ** 4 
        min_r_cut_off =  10.0_wp ** 4
        
!       Initialise pointers:
        indxyz = 0

!       Timing:
        CALL CPU_TIME(time0)
        modh = MAX ( SIZE ( g ) / 20, 1 )

!       Counter for printing:
        pair  = 0 

!       Loop through all atoms:
        coord:DO iat = 1, natom

!           Speed up for occupancy only refinement:
            IF ( COUNT ( pdb_2%refinable(1:11, iat ) ) == 0 ) THEN
!               Nothing to refine:
                CYCLE coord
            ENDIF

!           Just diagonal for now:
            jat  = iat

            occ_i = pdb_2%occ(iat)
            occ_j = pdb_2%occ(jat)

!           Save atom types for this pair:
            ityp     = pdb_2%atom_type(iat)
            ngauss   = pdb_2%atomsf_table%record(ityp)%ngauss
            jtyp     = pdb_2%atom_type(jat)
            jgauss   = pdb_2%atomsf_table%record(jtyp)%ngauss

!           FIXME: additional code (lower) is needed to treat this situation correctly:
            IF ( ngauss /=  jgauss ) THEN
                WRITE(*,*) iat, jat
                WRITE(*,*) ngauss,  pdb_2%atomsf_table%record(jtyp)%ngauss 
                CALL die('Oops... Cannot deal with differing NGAUSS for atom pair yet...', srname)
            ENDIF

!           Check:
            IF ( ngauss /=2 .AND. ngauss /= 5 ) THEN
                WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                CALL die('Programming error. Abnormal number of gaussians detected...',&
                          srname)
            ENDIF

!           Copy gaussian constants:
            a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
            b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss) + ( pdb_2%biso(iat) + b_scale )

            IF ( debug > 22 ) THEN
                WRITE(*,"(' a(1:ngauss)=', 15F7.2)") a(1:ngauss)
                WRITE(*,"(' b(1:ngauss)=', 15F7.2)") b(1:ngauss)
                WRITE(*,"(' pdb_2%biso(iat)=',F7.2)") pdb_2%biso(iat)
            ENDIF
  
!           Joint atom constants:
            occ = occ_i * occ_j

            gauss_outer:DO i = 1, ngauss
                gauss_inner:DO j = 1, i

                   IF ( i /= j ) THEN
                       gaussian_scale = 2.0_wp
                   ELSE
                       gaussian_scale = 1.0_wp
                   ENDIF

!                  Convert to one dimensional table (lower triange of 2-dim table):
                   l = (i * (i - 1)) / 2 + j

!                  Joint atom constants:
                   ajo(l) = gaussian_scale * a(i) * a(j) 
                   bjo(l) = b(i) + b(j) 

                ENDDO gauss_inner
            ENDDO gauss_outer

!           Redefine ngauss:
            ngauss = ngauss * ( ngauss + 1 ) / 2

!           Debug:
            IF ( debug > 22 ) THEN
                WRITE(*,"(' ajo=', 15F8.2)") ajo(1:ngauss)
                WRITE(*,"(' bjo=', 15F8.2)") bjo(1:ngauss)
            ENDIF

!           Final joint atom constants (switching  back to usual names):
            a(1:ngauss)  = ajo(1:ngauss) * ( SQRT ( pi4 / bjo(1:ngauss) ) ) ** 3
            be(1:ngauss) = 4.0 * pi_squared / bjo(1:ngauss)

            IF ( debug > 22 ) THEN
                WRITE(*,"(' a =', 15F8.4)") a(1:ngauss)
                WRITE(*,"(' be=', 15F8.4)") be(1:ngauss)
            ENDIF

!           Arrays of type 1.0/( b(igauss) + biso(iat) + b_scale ) to speed up:
            bbb(1:ngauss) = 1.0_wp / bjo(1:ngauss)

            IF ( debug > 22 ) THEN
                WRITE(*,"(' bbb=', 15F8.4)") bbb(1:ngauss)
            ENDIF 

!           Superior floating radius, gives excellent results on RNAse for high B-values:
            r2_cut_off = integration_radius ( ajo(1:ngauss), bjo(1:ngauss) ) ** 2

            IF ( debug > 22 ) THEN
                WRITE(*,*) ' optimal diag r_cut_off=', SQRT ( r2_cut_off )
            ENDIF

!           Make sure sphere around any atom fits into the cell:
            r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
            r_cut_off  = SQRT ( r2_cut_off )

!           Accumulate statics for R_CUT_OFF:
            max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
            min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )         

            IF ( debug > 22 ) THEN
                WRITE(*,*) ' my_first_allocated_map=', my_first_allocated_map
                WRITE(*,*) ' diag r_cut_off=', r_cut_off
            ENDIF

!           Figure out box size around joint atom:
            box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!           Figure out refinable parameters:
            refine_atomic_xyz  = pdb_2%refinable(1:3,iat)
            refine_atomic_biso = pdb_2%refinable(4,iat)
            refine_atomic_occ  = pdb_2%refinable(5,iat)

!           Calculate matrix pointers:
            IF ( refine_atomic_biso ) THEN
                indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
!               Speedup array:
                bbb2(1:ngauss) = bbb(1:ngauss) ** 2
            ENDIF

            IF ( refine_atomic_occ ) THEN
                indocc = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
            ENDIF

!           Special debugging:
            IF ( debug > 22 ) THEN
                WRITE(*,*) iat, indxyz, refine_atomic_xyz
            ENDIF

!           Convert to fractions:
            xyz_frac = (/ map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(iat), &
                          map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(jat) /)
         
            xyz_iat = xyz_frac(1)

!           Initialize AMAT for chain rule:
            AMAT = 0.0_wp

!           Loop over all symmetry operators:
            symmetry:DO j = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops
                xyz_jat = map_1(my_first_allocated_map)%sp_group%SYM(j) * xyz_frac(2)

!               Initialsing for each symmetry operator to build large matrix:
                BMAT = 0.0_wp

!               Calculate joint atom position for H1 term:
                Hterms:DO run = 1, 2

!                   Better not to run cycle against MY_FIRST_ALLOCATED_MAP:
                    my_first_allocated_map = run

!                   Second map is used to obtain H2 terms:
                    H2 = my_first_allocated_map == 2

!                   Figure out how to calculate centre of joint atom:
                    IF ( .NOT. H2 ) THEN
                        centre_of_atom_in_fractions = xyz_iat - xyz_jat
                    ELSE
                        centre_of_atom_in_fractions = xyz_iat + xyz_jat
                    ENDIF

                    IF ( debug > 22 ) THEN
                        xyz = centre_of_atom_in_fractions
                        IF ( .NOT. H2 ) THEN
                            WRITE(*,"( 'H1 centre= ', 3F8.5)") xyz
                        ELSE
                            WRITE(*,"( 'H2 centre= ', 3F8.5)") xyz
                        ENDIF
                    ENDIF

!                   Figure out closest grid point:
                    closest_grid_point = NINT ( centre_of_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )

!                   Box limits:
                    box_lower_limit = closest_grid_point - box_size
                    box_upper_limit = closest_grid_point + box_size

!                   Create a box around the joint atom and calculate its density and convolute with the map:
                    boxzz:DO jw = box_lower_limit(3), box_upper_limit(3)

                        DO jv = box_lower_limit(2), box_upper_limit(2)

                            DO ju = box_lower_limit(1), box_upper_limit(1)

                                uvw = (/ ju, jv, jw /)
                                duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_atom_in_fractions

!                               High point of electron density convolution:
                                r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw   ! corresponds to r - r(n) + r(m)

                                IF ( r2 <= r2_cut_off ) THEN

                                    jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                    map_value = map_1(my_first_allocated_map)%array(jut(1), jut(2), jut(3))

!                                   Going to speed up in case of Biso/occncy refinement:
                                    ber2(1:ngauss)   = be(1:ngauss) * r2
                                    exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                    abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)

!                                   Just diagonal elements:
                                    IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                       Convert to orthogonal system of coodinates (easier to calculate derivatives):
                                        duvwort = map_1(my_first_allocated_map)%ORT * duvw

!                                       Calculate Tronrud XXT matrix:
                                        DO l = 1, 3
                                            DO k = 1, 3
                                                XXT(l,k) = 4.0_wp * duvwort(k) * duvwort(l)
                                            ENDDO
                                        ENDDO


!                                       Additional array for speed up:
                                        abexp_be(1:ngauss) =  (occ * map_value) * abexp(1:ngauss) * be(1:ngauss)

!                                       Build 3x3 orthogonal matrix:
                                        DO k = 1, ngauss
                                            BMAT = BMAT + abexp_be(k) * (GMAT - be(k) * XXT) 
                                        ENDDO

                                    ENDIF

!                                   Switch sign in case of H2 terms for Biso and occncy:
                                    IF ( H2 ) map_value = -map_value

                                    IF ( refine_atomic_biso ) THEN
                                        g(indbiso) = g(indbiso) + occ * map_value                                  &
                                                                * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)          &
                                                                * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp)      &
                                                                + 15.0_wp / 4.0_wp ) ) 
                                                    
                                    ENDIF   
                           
                                    IF ( refine_atomic_occ ) THEN
                                        g(indocc) = g(indocc) + map_value * SUM ( abexp(1:ngauss) )
                                    ENDIF

                                ENDIF
                            ENDDO
                        ENDDO
                    ENDDO boxzz
                ENDDO Hterms

                IF ( ANY ( refine_atomic_xyz ) ) THEN

!                   Build 3x3Nsym matrix:
                    AMAT(1:3,3*(j-1)+1:3*(j-1)+3) = BMAT
                ENDIF

!               A small caveat here -> my_first_allocated_map will have value of "2" in outer loops.
            ENDDO symmetry

            chainrule:IF ( ANY ( refine_atomic_xyz ) ) THEN

!               Finish with XYZ chain rule for this atom, Asym is on the right side:
                BMAT = MATMUL ( AMAT, ASYM )

!               Extract and add diagonal terms only:
                l = 0
                DO k = 1, 3
                    IF ( refine_atomic_xyz(k) ) THEN
                        l = l + 1
                        g(indxyz+l) = g(indxyz+l) + BMAT(k,k)
                    ENDIF
                ENDDO

!               Debug if necessary:
                debugg:IF ( debug > 25 ) THEN

                    IF ( map_1(my_first_allocated_map)%sp_group%number_of_symops <= 24 ) THEN

                        DO k = 1, 3
                            WRITE(*,"(' AMAT=',72ES9.2)") &
                            (AMAT(k,l), l = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)
                        ENDDO

                    ENDIF
   
                    DO k = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops
                        WRITE(*,"(' ASYM=', 3F9.5)") (ASYM(k,l), l = 1, 3)
                    ENDDO

                    DO k = 1, 3
                        WRITE(*,"(' BMAT=', 3ES9.2)") (BMAT(k,l), l = 1, 3)
                    ENDDO

                ENDIF debugg

            ENDIF chainrule

!           Monitor CPU:
            l = 0
            DO k = 1, 11

                IF ( pdb_2%refinable(k,iat) ) THEN
                    l = l + 1
                    pair = pair + 1
                    IF ( MOD ( pair, modh ) == 0 ) THEN
                        CALL CPU_TIME ( time1 )
                        WRITE(*,"(' NORMAL_MATRIX_DIAGONAL> ',                                            &
                       &'matrix element[', A, ',', A, ']=', T62, ES12.5, 2X,'CPU time=',F8.1, ' s', I10)")&
                        TRIM ( int_to_c ( pair ) ), TRIM ( int_to_c ( pair ) ), g(pair), time1 - time0, pair
                    ENDIF

!                   Check matrix elements for non-positive values:
                    IF ( g(indxyz+l) <= 0.0_wp ) THEN
                        WRITE(*,"(' NORMAL_MATRIX_DIAGONAL> ', 'pointer=', I8, ' diag=', ES12.5, ' iat=', I6, ' ref. param=',I3)")&
                                  indxyz+l, g(indxyz+l), iat, l
                        WRITE(*,"(' NORMAL_MATRIX_DIAGONAL> ', 'Biso=', F8.2, ' occ=', F6.2, ' r_cut_off=', F8.2)")&
                                  pdb_2%biso(iat), pdb_2%occ(iat), r_cut_off

!                       Bailing out:
                        CALL die('Programming error. Zero or negative normal matrix diagonal element. Refinement halted.',&
                                 srname)

                    ENDIF

                ENDIF

            ENDDO

!           Increment pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11,iat) )

        ENDDO coord

!       Report stats:
        WRITE(*,"(' NORMAL_MATRIX_DIAGONAL> ', 'R_CUT_OFF: max=', F8.1, ' min=', F8.1)")&
        max_r_cut_off, min_r_cut_off

!       Free memory:
        CALL deallocate_array(AMAT)
        CALL deallocate_array(ASYM)

    END SUBROUTINE normal_matrix_diagonal

    SUBROUTINE dense_normal_matrix(HMAT, map_1, pdb_2)
!
!       Purpose:
!       =======
!
!       Calculates diagonal elements G of normal matrix (H1 + H2 terms) with respect to
!       X,Y,Z,B and Q params.
!
!       Supply two maps a la Tronrud and you are done. Just one map will result in
!       calculation of H1 terms only.
!
!       Date:                 Programmer:             History of changes:
!       ====                  ==========              ==================
!       January 10, 2008      B.Strokopytov           Original code
!       January 15, 2008      B.Strokopytov           Added BXYZ_DERIV and OXYZ_DERIV arrays for 
!                                                     chain rule speed-up.
!
!       Chain rule note:
!       ===============
!       In orthogonal coordinate system we have:
!
!                     T      -1    -1 -1 -1       -1
!       [Asym] = [OAD] = [OAD] = [D  A  O  ] = [OA  D]
!
!       Obviously:
!             T
!       [Asym] = [OAD] 
!       
!       Transpose is due to chain rule, if we going to multiply from the right:
!         2    2      T
!       [d R/dx][Asym]
!
!       Remaining problems:
!       ==================
!       How to get rid of the symmetry loop? Looks like a horrible CPU problem
!       in cubic space groups. 
!
        REAL(KIND=wp), DIMENSION(:,:),            INTENT(INOUT) :: HMAT
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: map_1
        TYPE(pdb),                                INTENT(IN)    :: pdb_2
!       Local variables:
        LOGICAL                                                 :: h2
        TYPE(vector)                                            :: centre_of_atom_in_fractions
        TYPE(vector_int)                                        :: closest_grid_point
        TYPE(vector_int)                                        :: box_size
        INTEGER,       DIMENSION(3)                             :: jut
        INTEGER,       DIMENSION(3)                             :: box_lower_limit
        INTEGER,       DIMENSION(3)                             :: box_upper_limit
        TYPE(vector)                                            :: deort_diag
        TYPE(vector)                                            :: ort_diag
!       Normal atom params:
        INTEGER                                                 :: ityp
        INTEGER                                                 :: jtyp
        INTEGER                                                 :: ngauss
        INTEGER                                                 :: igauss
        INTEGER                                                 :: jgauss
        REAL(KIND=wp)                                           :: b_scale
!       Scattering:
!       For diagonal terms 15 gaussians would be the MAX number,
!       but for non-diagonal terms need 25 (looking into the future):
        INTEGER,       PARAMETER                                :: MAXGAU = 25
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: a
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: b
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: ai
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bi
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: aj
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bj
!       Joint atom:
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: ajo 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bjo
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: be
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bbb 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: bbb2 
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: ber2
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: exp_be
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: abexp
        REAL(KIND=wp), DIMENSION(MAXGAU)                        :: abexp_be
        REAL(KIND=wp)                                           :: gaussian_scale
        REAL(KIND=wp)                                           :: sum_bbb_15
        REAL(KIND=wp)                                           :: sum_bbb_52
        REAL(KIND=wp)                                           :: sum_abexp_be
        REAL(KIND=wp), DIMENSION(3)                             :: gg
!       Distance:
        REAL(KIND=wp)                                           :: r2
        TYPE(vector_int)                                        :: uvw
        TYPE(vector)                                            :: duvw
        REAL(KIND=wp), DIMENSION(3)                             :: duvwort
!       Cut-off:
        REAL(KIND=wp)                                           :: r2_cut_off
        REAL(KIND=wp)                                           :: r_cut_off
        REAL(KIND=wp)                                           :: max_r_cut_off
        REAL(KIND=wp)                                           :: min_r_cut_off
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                           :: r2_cell_max
!       Can handle 2 maps only:
        REAL(KIND=wp)                                           :: map_value
        REAL(KIND=wp)                                           :: occ
        REAL(KIND=wp)                                           :: occ_i
        REAL(KIND=wp)                                           :: occ_j
!       Maps variables:
        INTEGER                                                 :: number_of_maps
        INTEGER                                                 :: my_first_allocated_map
!       Refinement:        
        INTEGER                                                 :: np
        INTEGER                                                 :: natom
        LOGICAL,       DIMENSION(3,2)                           :: refine_atomic_xyz
        LOGICAL,       DIMENSION(2)                             :: refine_atomic_biso
        LOGICAL,       DIMENSION(2)                             :: refine_atomic_occ
!       Pointers:
        INTEGER                                                 :: indxyz
        INTEGER                                                 :: indbiso
        INTEGER                                                 :: indocc
        INTEGER                                                 :: jndxyz
        INTEGER                                                 :: jndbiso
        INTEGER                                                 :: jndocc 
!       Counters:
        INTEGER                                                 :: iat
        INTEGER                                                 :: jat
        INTEGER                                                 :: k 
        INTEGER                                                 :: kk
        INTEGER                                                 :: l
        INTEGER                                                 :: ll
        INTEGER                                                 :: run
!       Symmetry:
        TYPE(vector)                                            :: xyz_iat
        TYPE(vector)                                            :: xyz_jat
        TYPE(vector),  DIMENSION(2)                             :: xyz_frac
        REAL(KIND=wp), DIMENSION(3,3)                           :: GMAT
        REAL(KIND=wp), DIMENSION(3,3)                           :: XXT
!       Additional chain rule arrays:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE              :: AMAT
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE              :: ASYM
        REAL(KIND=wp), DIMENSION(3,3)                           :: BMAT
        REAL(KIND=wp), DIMENSION(3,3)                           :: SMAT
        REAL(KIND=wp), DIMENSION(3)                             :: bxyz_deriv
        REAL(KIND=wp), DIMENSION(3)                             :: oxyz_deriv
!       General counters:
        INTEGER                                                 :: i
        INTEGER                                                 :: j
!       Symmetry counter: 
        INTEGER                                                 :: m 
!       Box counters:
        INTEGER                                                 :: ju
        INTEGER                                                 :: jv
        INTEGER                                                 :: jw
!       CPU monitoring:
        REAL(KIND=sp)                                           :: time0
        REAL(KIND=sp)                                           :: time1
!       Printing:
        INTEGER(KIND=eb)                                        :: modh
        INTEGER(KIND=eb)                                        :: pair
        CHARACTER(LEN=32),                               SAVE   :: srname = 'dense_normal_matrix'
!       Test:
        REAL(KIND=wp), DIMENSION(3)                             :: xyz

!       Checkz: 
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die('Programming error. pdb_2 has not been initialized.', &
                     srname)
        ENDIF

!       This must be set from the start:
        my_first_allocated_map = 1

!       Check map allocation:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Programming error. Map array has not been allocated properly.', &
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( map_1(my_first_allocated_map) ) ) THEN
            CALL die('Programming error. MAP_1(1) has not been allocated properly.',&
                     srname)
        ELSE

            number_of_maps = SIZE ( map_1 )
            IF ( number_of_maps < 2 ) THEN
                WRITE(*,*) number_of_maps
                CALL die('Programming error. Need 2 maps....', srname)
            ENDIF

        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' NORMAL_MATRIX_DIAGONAL> ', 'natom=', I6)") natom
            CALL die ( 'Programming error. Need at least one atom for density generation.', srname)
        ENDIF

!       Long and hard check: 
        np = COUNT ( pdb_2%refinable )
        IF ( np ** 2 /=  SIZE ( HMAT ) ) THEN
            WRITE(*,*) np, np ** 2, SIZE ( HMAT )
            CALL die('Programming error. Incorrect size of HMAT array.', srname)
        ENDIF
       
        CALL messag('Going to calculate '//TRIM ( int_to_c ( SIZE ( HMAT ) ) )//' matrix elements... Wait...',&
                    srname)

!       Initialize dense matrix array:
        HMAT = 0.0_wp

!       Allocate additional matrices for chain rule:
        CALL allocate_array(ASYM, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops, 3)
        CALL allocate_array(AMAT, 3, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)

!       Build 3Nsym x 3 matrix applying chain rule (see notes above):
        DO m = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops

!           Calculate position of the block:  
            l = 3 * (m - 1)
!                                        T 
!           Calculating transpose [Asym ] :
            BMAT = map_1(my_first_allocated_map)%ORT             &
                 * .SYMA. map_1(my_first_allocated_map)%sp_group%SYM(m) &
                 *        map_1(my_first_allocated_map)%DEORT

!           To avoid Intel compiler warnings (in debugging mode):
            ASYM(l+1:l+3,1:3) = BMAT
        ENDDO

!       Initialize 3x3 diagonal matrix:
        GMAT = 0.0_wp
        DO l = 1, 3
            GMAT(l,l) = 2.0_wp
        ENDDO

!       Other params:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Calculate maximum sphere radius around any atom:
        ort_diag    = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Set for statistics:
        max_r_cut_off = -10.0_wp ** 4 
        min_r_cut_off =  10.0_wp ** 4
        
!       Initialise matirx pointer along IAT atoms:
        indxyz = 0

!       Timing:
        CALL CPU_TIME(time0)
        modh = MAX ( SIZE ( HMAT ) / 100, 1 )

!       Counter for printing:
        pair  = 0 

!       Loop through all atoms:
        coord:DO iat = 1, natom

            occ_i = pdb_2%occ(iat)

!           Initialize matrix pointer for JAT atoms:
            jndxyz = 0
            ityp     = pdb_2%atom_type(iat)
            igauss   = pdb_2%atomsf_table%record(ityp)%ngauss

!           Copy gaussian constants for IAT atom:
            ai(1:igauss) = pdb_2%atomsf_table%record(ityp)%a(1:igauss)
            bi(1:igauss) = pdb_2%atomsf_table%record(ityp)%b(1:igauss) + ( pdb_2%biso(iat) + b_scale )
!           Check:

            IF ( igauss /=2 .AND. igauss /= 5 ) THEN
                    WRITE(*,*) 'igauss=', igauss, ' atom=', pdb_2%atom_name(iat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                              srname)
            ENDIF

            DO jat = 1, natom

                IF ( COUNT ( pdb_2%refinable(1:11, jat ) ) == 0 ) THEN
!                   Nothing to refine:
                    CYCLE
                ENDIF

                occ_j = pdb_2%occ(jat)

!               Save atom types for this pair:
                jtyp     = pdb_2%atom_type(jat)
                jgauss   = pdb_2%atomsf_table%record(jtyp)%ngauss

                IF ( jgauss /=2 .AND. jgauss /= 5 ) THEN
                    WRITE(*,*) 'jgauss=', jgauss, ' atom=', jat, pdb_2%atom_name(jat)
                    CALL die('Programming error. Abnormal number of gaussians detected...',&
                              srname)
                ENDIF


                IF ( debug > 22 ) THEN
                    WRITE(*,"(' a(1:ngauss)=', 15F7.2)") a(1:ngauss)
                    WRITE(*,"(' b(1:ngauss)=', 15F7.2)") b(1:ngauss)
                    WRITE(*,"(' pdb_2%biso(iat)=',F7.2)") pdb_2%biso(iat)
                ENDIF

!               Copy gaussian constants for JAT atom:
                aj(1:jgauss) = pdb_2%atomsf_table%record(jtyp)%a(1:jgauss)
                bj(1:jgauss) = pdb_2%atomsf_table%record(jtyp)%b(1:jgauss) + ( pdb_2%biso(jat) + b_scale )

  
!               Joint atom constants:
                occ = occ_i * occ_j

                IF ( jat == iat ) THEN

!                   Atom on the diagonal:
                    gauss_outer:DO i = 1, igauss
                        gauss_inner:DO j = 1, i

                            IF ( i /= j ) THEN
                                gaussian_scale = 2.0_wp
                            ELSE
                                gaussian_scale = 1.0_wp
                            ENDIF

!                           Convert to one dimensional table (lower triange of 2-dim table):
                            l = (i * (i - 1)) / 2 + j

!                           Joint atom constants:
                            ajo(l) = gaussian_scale * ai(i) * ai(j) 
                            bjo(l) = bi(i) + bi(j) 

                        ENDDO gauss_inner
                    ENDDO gauss_outer
            
!                   Calculate ngauss:
                    ngauss = (igauss * ( igauss + 1 )) / 2

                ELSE

!                   Different atoms:
                    l = 0
                    DO i = 1, igauss
                        DO j = 1, jgauss
                            l = l + 1

!                           Joint atom constants:
                            ajo(l) = ai(i) * aj(j)
                            bjo(l) = bi(i) + bj(j)
                        ENDDO
                    ENDDO

!                   Redefine ngauss:
                    ngauss = l

!                   Check number of gaussians:
                    IF ( ngauss /= 4 .AND. ngauss /= 10 .AND. ngauss /= 25 ) THEN
                        WRITE(*,*) ngauss
                        CALL die('Programming error. Number of gaussians is not equal to the expected value.',&
                                 srname)
                    ENDIF
                ENDIF

!               Debug:
                IF ( debug > 22 ) THEN
                    WRITE(*,"(' ajo=', 15F8.2)") ajo(1:ngauss)
                    WRITE(*,"(' bjo=', 15F8.2)") bjo(1:ngauss)
                ENDIF

!               Final joint atom constants (switching  back to usual names):
                a(1:ngauss)  = ajo(1:ngauss) * ( SQRT ( pi4 / bjo(1:ngauss) ) ) ** 3
                be(1:ngauss) = 4.0_wp * pi_squared / bjo(1:ngauss)

                IF ( debug > 22 ) THEN
                    WRITE(*,"(' a =', 15F8.4)") a(1:ngauss)
                    WRITE(*,"(' be=', 15F8.4)") be(1:ngauss)
                ENDIF

!               Array to speed up:
                bbb(1:ngauss) = 1.0_wp / bjo(1:ngauss)

                IF ( debug > 22 ) THEN
                    WRITE(*,"(' bbb=', 15F8.4)") bbb(1:ngauss)
                ENDIF 

!               Superior floating radius, gives excellent results for RNAse in case of high B-values:
                r2_cut_off = integration_radius ( ajo(1:ngauss), bjo(1:ngauss) ) ** 2

                IF ( debug > 22 ) THEN
                    WRITE(*,*) ' optimal diag r_cut_off=', SQRT ( r2_cut_off )
                ENDIF

!               Make sure sphere around any atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Accumulate statics for R_CUT_OFF:
                max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
                min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )         

                IF ( debug > 22 ) THEN
                    WRITE(*,*) ' my_first_allocated_map=', my_first_allocated_map
                    WRITE(*,*) ' diag r_cut_off=', r_cut_off
                ENDIF

!               Figure out box size around joint atom:
                box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!               Figure out refinable parameters:
                DO k = 1, 2
                    IF ( k == 1 ) THEN
                        l = iat
                    ELSE
                        l = jat
                    ENDIF
                    refine_atomic_xyz(1:3,k)  = pdb_2%refinable(1:3,l)
                    refine_atomic_biso(k)     = pdb_2%refinable(4,l)
                    refine_atomic_occ(k)      = pdb_2%refinable(5,l)
                ENDDO
                

!               Calculate matrix pointers:
                IF ( refine_atomic_biso(1) ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                ENDIF

                IF ( refine_atomic_occ(1) ) THEN
                    indocc = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                ENDIF

!               Using JAT now:
                IF ( refine_atomic_biso(2) ) THEN
                    jndbiso = jndxyz + COUNT ( pdb_2%refinable(1:4, jat) )
                ENDIF

                IF ( refine_atomic_occ(2) ) THEN
                    jndocc = jndxyz + COUNT ( pdb_2%refinable(1:5, jat) )
                ENDIF

!               Special arrays:
                IF ( ANY ( refine_atomic_biso ) ) THEN
!                   Speedup array:
                    bbb2(1:ngauss) = bbb(1:ngauss) ** 2
                ENDIF


!               Special debugging:
                IF ( debug > 22 ) THEN
                    WRITE(*,*) iat, indxyz, refine_atomic_xyz
                    WRITE(*,*) jat, jndxyz, pdb_2%refinable(1:3,jat)
                ENDIF

!               Convert to fractions:
                xyz_frac = (/ map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(iat), &
                              map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(jat) /)
         
                xyz_iat = xyz_frac(1)

!               Initialize AMAT for chain rule:
                AMAT = 0.0_wp

!               Loop over all symmetry operators:
                symmetry:DO m = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops

!                   This will cause some asymmetry in calculations (see below):
                    xyz_jat = map_1(my_first_allocated_map)%sp_group%SYM(m) * xyz_frac(2)

!                   Prepare symmetry matrix for chain rule for certain derivatives:
                    l = 3 * ( m - 1)
                    SMAT = TRANSPOSE ( ASYM(l+1:l+3,1:3) )

!                   Initialising coordinate matrix for each symmetry operator to build even larger matrix:
                    BMAT = 0.0_wp

!                   Small arrays to speed up chain rule application for BX,BY,BZ & OX,OY,OZ derivatives:
                    bxyz_deriv = 0.0_wp
                    oxyz_deriv = 0.0_wp 

!                   Calculate joint atom position for H1 term:
                    Hterms:DO run = 1, 2

!                       Better not to run cycle against MY_FIRST_ALLOCATED_MAP:
                        my_first_allocated_map = run

!                       Second map is used to obtain H2 terms:
                        H2 = my_first_allocated_map == 2

!                       Figure out how to calculate centre of joint atom:
                        IF ( .NOT. H2 ) THEN
                            centre_of_atom_in_fractions = xyz_iat - xyz_jat
                        ELSE
                            centre_of_atom_in_fractions = xyz_iat + xyz_jat
                        ENDIF

                        IF ( debug > 22 ) THEN

                            xyz = centre_of_atom_in_fractions
                            IF ( .NOT. H2 ) THEN
                                WRITE(*,"( 'H1 centre= ', 3F8.5)") xyz
                            ELSE
                                WRITE(*,"( 'H2 centre= ', 3F8.5)") xyz
                            ENDIF

                        ENDIF

!                       Figure out closest grid point:
                        closest_grid_point = NINT ( centre_of_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )

!                       Box limits:
                        box_lower_limit = closest_grid_point - box_size
                        box_upper_limit = closest_grid_point + box_size

!                       Create a box around the joint atom and calculate its density and convolute with the map:
                        boxzz:DO jw = box_lower_limit(3), box_upper_limit(3)

                            DO jv = box_lower_limit(2), box_upper_limit(2)

                                DO ju = box_lower_limit(1), box_upper_limit(1)

                                    uvw = (/ ju, jv, jw /)
                                    duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_atom_in_fractions

!                                   High point of electron density convolution:
                                    r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw   ! corresponds to r - r(n) + r(m)

                                    IF ( r2 <= r2_cut_off ) THEN

                                        jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                        map_value = map_1(my_first_allocated_map)%array(jut(1), jut(2), jut(3))

!                                       Going to speed up in case of Biso/occncy refinement:
                                        ber2(1:ngauss)   = be(1:ngauss) * r2
                                        exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                        abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)

!                                       Additional array for speed up:
                                        IF ( ANY ( refine_atomic_xyz ) ) THEN
                                            abexp_be(1:ngauss) = abexp(1:ngauss) * be(1:ngauss)
                                            sum_abexp_be = SUM ( abexp_be(1:ngauss) )

!                                           Convert to orthogonal system of coodinates (easier to calculate derivatives):
                                            duvwort = 2.0_wp * map_1(my_first_allocated_map)%ORT * duvw
                                        ENDIF
 
!                                       Special array for speed-up:
                                        IF ( ANY ( refine_atomic_biso ) ) THEN 
                                            sum_bbb_52 = SUM ( abexp_be(1:ngauss) * bbb(1:ngauss) &
                                                             * (5.0_wp / 2.0_wp - ber2(1:ngauss)) )
                                        ENDIF

                                        IF ( ( refine_atomic_biso(1) .AND. refine_atomic_occ(2) ) .OR. &
                                             ( refine_atomic_biso(2) .AND. refine_atomic_occ(1) ) ) THEN
                                             sum_bbb_15 = SUM ( abexp(1:ngauss) * bbb(1:ngauss) * (-1.5_wp + ber2(1:ngauss)) )
                                        ENDIF
                                    
                                        anyxyz:IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

!                                           Makes sense to calculate if JAT atom has refinable XYZ:
                                            IF (  ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                               BUG CORRECTED Jan 2008 BVS -> calcn of duvwort moved out from here:

!                                               Calculate Tronrud XXT matrix:
                                                DO l = 1, 3
                                                    DO k = 1, 3

!                                                       DUVWORT has been scaled by 2 already:
                                                        XXT(l,k) = duvwort(k) * duvwort(l)
                                                    ENDDO
                                                ENDDO


!                                               Build 3x3 matrix for orthogonal coordinates:
                                                BMAT = BMAT + map_value * occ * (sum_abexp_be * GMAT           &
                                                            - SUM ( abexp_be(1:ngauss) * be(1:ngauss) ) * XXT)

                                            ENDIF


!                                           Check whether Biso refinable for JAT atom:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Accumulate matrix values wrt XB+YB+ZB(14):

                                                gg = -(occ * map_value * sum_bbb_52) * duvwort

!                                               Switch sign for H2 terms:
                                                IF ( H2 ) gg = -gg

!                                               General case -> not every XYZ(iat) is refinable:
                                                IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN
                                                    l = 0
                                                    DO k = 1, 3
                                                        IF ( refine_atomic_xyz(k,1) ) THEN
                                                            l = l + 1
                                                            HMAT(indxyz+l,jndbiso) = HMAT(indxyz+l,jndbiso) + gg(k)
                                                        ENDIF
                                                    ENDDO
                                                ENDIF

                                            ENDIF

!                                           Accumulate derivatives XO, YO, ZO(15) for XYZ rows:
                                            IF ( refine_atomic_occ(2) ) THEN

!                                               JAT atom seems to have refinable occupancy:
                                                gg = occ_i * map_value * sum_abexp_be * duvwort

!                                               Switch sign for H2 terms:
                                                IF ( H2 ) gg = -gg

!                                               Considering general case -> not every XYZ(iat) is refinable:
                                                IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN                                                

                                                    l = 0
                                                    DO k = 1, 3
                                                        IF ( refine_atomic_xyz(k,1) ) THEN
                                                            l = l + 1
                                                            HMAT(indxyz+l,jndocc) = HMAT(indxyz+l,jndocc) + gg(k)
                                                        ENDIF
                                                    ENDDO
                                                ENDIF

                                            ENDIF

                                        ENDIF anyxyz

!                                       -------------- Biso row ---------------------------

                                        IF ( refine_atomic_biso(1) ) THEN

!                                           Accumulate derivatives for BX, BY, BZ(14):
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

                                                bxyz_deriv = bxyz_deriv + (occ * map_value * sum_bbb_52) * duvwort

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Accumalate BB matrix elements (12):
                                                IF ( .NOT. H2 ) THEN
                                                    HMAT(indbiso,jndbiso) = HMAT(indbiso,jndbiso) + occ * map_value      &
                                                                          * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)      &
                                                                          * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp)  &
                                                                          + 15.0_wp / 4.0_wp ) )
                                                ELSE
                                                    HMAT(indbiso,jndbiso) = HMAT(indbiso,jndbiso) - occ * map_value      &
                                                                          * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)      &
                                                                          * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp)  &
                                                                          + 15.0_wp / 4.0_wp ) )

                                                ENDIF

                                            ENDIF                                                    

!                                           Check whether JAT atom has refinable occ:
                                            IF ( refine_atomic_occ(2) ) THEN

!                                               Accumulate BQ matrix elements (16):
                                                IF ( .NOT. H2 ) THEN
                                                    HMAT(indbiso,jndocc) = HMAT(indbiso,jndocc) + occ_i * map_value * sum_bbb_15 
                                                ELSE
                                                    HMAT(indbiso,jndocc) = HMAT(indbiso,jndocc) - occ_i * map_value * sum_bbb_15
                                                ENDIF

                                            ENDIF

                                        ENDIF   
                           
!                                       ------------- Now occupancy row itself, note OCC_J is used --------------:
                                        occncy:IF ( refine_atomic_occ(1) ) THEN

!                                           Check whether JAT atom has refinable XYZ:
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                               OX+OY+OZ interaction (15):
                                                oxyz_deriv = oxyz_deriv - (occ_j * map_value * sum_abexp_be) * duvwort

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               QB interaction(16), note plus sign here:
                                                IF ( .NOT. H2 ) THEN
                                                    HMAT(indocc, jndbiso) = HMAT(indocc, jndbiso) + occ_j * map_value * sum_bbb_15
                                                ELSE
!                                                   Switch sign for H2 terms:
                                                    HMAT(indocc, jndbiso) = HMAT(indocc, jndbiso) - occ_j * map_value * sum_bbb_15
                                                ENDIF

                                            ENDIF

!                                           Check whether JAT atom has refinable occncy and apply(13):
                                            IF ( refine_atomic_occ(2) ) THEN

!                                               Switch sign for H2 terms:
                                                IF ( .NOT. H2 ) THEN
                                                    HMAT(indocc,jndocc) = HMAT(indocc,jndocc) + map_value * SUM ( abexp(1:ngauss) )
                                                ELSE
                                                    HMAT(indocc,jndocc) = HMAT(indocc,jndocc) - map_value * SUM ( abexp(1:ngauss) )
                                                ENDIF

                                            ENDIF

                                        ENDIF occncy


                                    ENDIF
                                ENDDO
                            ENDDO
                        ENDDO boxzz
                    ENDDO Hterms

                    IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) )  THEN

!                       Build 3x3Nsym matrix:
                        AMAT(1:3,3*(m-1)+1:3*(m-1)+3) = BMAT
                    ENDIF

!                   BUG CORRECTED JAN 2008 BVS (CANNOT MERGE WITH PREVIOUS IF):
                    IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                       Take care of chain rule for BXYZ_DERIV & OXYZ_DERIV:
                        IF ( refine_atomic_biso(1) ) THEN
                            bxyz_deriv = MATMUL ( SMAT, bxyz_deriv )
                        ENDIF

                        IF ( refine_atomic_occ(1) ) THEN
                             oxyz_deriv = MATMUL ( SMAT, oxyz_deriv )
                        ENDIF                      

!                       Finish with chain rule for two classes of second_derivatives:
                        IF ( refine_atomic_biso(1) ) THEN

!                           Consider general case -> not every XYZ(jat) is refinable:
                            l = 0
                            DO k = 1, 3
                                IF ( refine_atomic_xyz(k,2) ) THEN
                                    l = l + 1
                                    HMAT(indbiso,jndxyz+l) = HMAT(indbiso,jndxyz+l) + bxyz_deriv(k)
                                ENDIF
                            ENDDO
                        ENDIF

                        IF ( refine_atomic_occ(1) ) THEN

!                           Consider general case -> not every XYZ(jat) is refinable:
                            l = 0
                            DO k = 1, 3
                                IF ( refine_atomic_xyz(k,2) ) THEN
                                    l = l + 1
                                    HMAT(indocc, jndxyz+l) = HMAT(indocc, jndxyz+l) + oxyz_deriv(k)
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
!                   A small caveat here -> my_first_allocated_map will have value of "2" in outer loops.

                ENDDO symmetry

                chainrule:IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                   Finish with XYZ chain rule for this atom, ASYM is on the right side:
                    BMAT = MATMUL ( AMAT, ASYM )

!                   Extract and add refinable items only:
                    l = 0
                    DO k = 1, 3
                        IF ( refine_atomic_xyz(k,1) ) THEN
                            l = l + 1
                            ll = 0
                            DO kk = 1, 3
                                IF ( refine_atomic_xyz(kk,2) ) THEN
                                    ll = ll + 1
                                    HMAT(indxyz+l,jndxyz+ll) = BMAT(k,kk)
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDDO

!                   Debug if necessary:
                    debugg:IF ( debug > 25 ) THEN

                        IF ( map_1(my_first_allocated_map)%sp_group%number_of_symops <= 24 ) THEN

                            DO k = 1, 3
                                WRITE(*,"(' AMAT=',72ES9.2)") &
                                (AMAT(k,l), l = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)
                            ENDDO

                        ENDIF
   
                        DO k = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops
                            WRITE(*,"(' ASYM=', 3F9.5)") (ASYM(k,l), l = 1, 3)
                        ENDDO

                        DO k = 1, 3
                            WRITE(*,"(' BMAT=', 3ES9.2)") (BMAT(k,l), l = 1, 3)
                        ENDDO

                    ENDIF debugg

                ENDIF chainrule

!               Monitor CPU:
                l = 0
                monitor:DO k = 1, 11

                    outer_if:IF ( pdb_2%refinable(k,iat) ) THEN
                        l = l + 1
                        ll = 0
                        DO kk = 1, 11
                            IF ( pdb_2%refinable(kk,jat) ) THEN
                                ll = ll + 1
                                pair = pair + 1

                                IF ( MOD ( pair, modh ) == 0 ) THEN
                                    CALL CPU_TIME ( time1 )
                                    WRITE(*,"(' DENSE_NORMAL_MATRIX> ',                                                &
                                   &'matrix element[', A, ',', A, ']=', T62, ES12.5, 2X,'CPU time=',F8.1, ' s', I10)") &
                                    TRIM ( int_to_c ( indxyz+l ) ), TRIM ( int_to_c ( jndxyz+ll ) ),&
                                    HMAT(indxyz+l,jndxyz+ll), time1 - time0, pair
                                ENDIF

!                               Check diagonal matrix elements for non-positive values:
                                IF ( indxyz+l == jndxyz+ll .AND. HMAT(indxyz+l,jndxyz+ll) <= 0.0_wp ) THEN
                                    WRITE(*,"(' DENSE_NORMAL_MATRIX> ', 'matirx element=', 2I8, ' HMAT=', ES12.5, &
                                   &' iat=', I6, ' ref. param=', I3)") &
                                    indxyz+l, jndxyz+ll, HMAT(indxyz+l,jndxyz+ll), iat, l

                                    WRITE(*,"(' DENSE_NORMAL_MATRIX> ', 'Biso=', ES9.2, ' occ=', ES9.2, ' r_cut_off=', F8.2)") &
                                    pdb_2%biso(iat), pdb_2%occ(iat), r_cut_off

!                                   Bailing out:
                                    CALL die('Programming error. Zero or negative normal matrix diagonal element. &
                                             & Refinement halted. Sorry.', srname)

                                ENDIF


                            ENDIF
                        ENDDO

                    ENDIF outer_if

                ENDDO monitor

!               Increment pointer for JAT atom loop:
                jndxyz = jndxyz + COUNT ( pdb_2%refinable(1:11,jat) )

            ENDDO

!           Increment pointer for IAT atom loop:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11,iat) )

        ENDDO coord

!       Report stats:
        WRITE(*,"(' DENSE_NORMAL_MATRIX> ', 'R_CUT_OFF: max=', F8.1, ' min=', F8.1)")&
        max_r_cut_off, min_r_cut_off

!       Free memory:
        CALL deallocate_array(AMAT)
        CALL deallocate_array(ASYM)

    END SUBROUTINE dense_normal_matrix

    SUBROUTINE sparse_normal_matrix(H, map_1, pdb_2, pair_list, nnz, diag)
!
!       Purpose:
!       =======
!
!       Calculates sparse matrix H (H1 + H2 terms) with respect to
!       X,Y,Z,B and Q params.
!
!       Supply two maps a la Tronrud and you are done. Just one map will result in
!       calculation of H1 terms only.
!
!       Date:                 Programmer:             History of changes:
!       ====                  ==========              ==================
!       January 10, 2008      B.Strokopytov           Original code
!       January 15, 2008      B.Strokopytov           Added BXYZ_DERIV and OXYZ_DERIV arrays for 
!                                                     chain rule speed-up.
!
!       Chain rule note:
!       ===============
!       In orthogonal coordinate system we have:
!
!                     T      -1    -1 -1 -1       -1
!       [Asym] = [OAD] = [OAD] = [D  A  O  ] = [OA  D]
!
!       Obviously:
!             T
!       [Asym] = [OAD] 
!       
!       Transpose is due to chain rule, if we going to multiply from the right:
!         2    2      T
!       [d R/dx][Asym]
!
!       Remaining problems:
!       ==================
!       How to get rid of the symmetry loop? Looks like a horrible CPU problem
!       in cubic space groups. 
!
        TYPE(coo_matrix),                                 INTENT(INOUT) :: H
        TYPE(map),           DIMENSION(:),   ALLOCATABLE, INTENT(IN)    :: map_1
        TYPE(pdb),                                        INTENT(IN)    :: pdb_2
        INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(IN)    :: pair_list
        INTEGER(KIND=eb),                                 INTENT(OUT)   :: nnz
        REAL(KIND=wp),       DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: diag
!       Local variables:
        LOGICAL                                                         :: h2
        TYPE(vector)                                                    :: centre_of_atom_in_fractions
        TYPE(vector_int)                                                :: closest_grid_point
        TYPE(vector_int)                                                :: box_size
        INTEGER,             DIMENSION(3)                               :: jut
        INTEGER,             DIMENSION(3)                               :: box_lower_limit
        INTEGER,             DIMENSION(3)                               :: box_upper_limit
        TYPE(vector)                                                    :: deort_diag
        TYPE(vector)                                                    :: ort_diag
!       Normal atom params:
        INTEGER                                                         :: ngauss
        REAL(KIND=wp)                                                   :: b_scale
!       Scattering:
!       For diagonal terms 15 gaussians would be the MAX number,
!       but for non-diagonal terms need 25 (looking into the future):
        INTEGER,             PARAMETER                                  :: MAXGAU = 25
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: a
!       Joint atom:
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: ajo 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bjo
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: be
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bbb 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bbb2 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: ber2
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: exp_be
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: abexp
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: abexp_be
        REAL(KIND=wp)                                                   :: sum_bbb_15
        REAL(KIND=wp)                                                   :: sum_bbb_52
        REAL(KIND=wp)                                                   :: sum_abexp_be
        REAL(KIND=wp),       DIMENSION(3)                               :: gg
!       Distance:
        REAL(KIND=wp)                                                   :: r2
        TYPE(vector_int)                                                :: uvw
        TYPE(vector)                                                    :: duvw
        REAL(KIND=wp),       DIMENSION(3)                               :: duvwort
!       Cut-off:
        REAL(KIND=wp)                                                   :: r2_cut_off
        REAL(KIND=wp)                                                   :: r_cut_off
        REAL(KIND=wp)                                                   :: max_r_cut_off
        REAL(KIND=wp)                                                   :: min_r_cut_off
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                                   :: r2_cell_max
!       Aniso section:
        LOGICAL                                                         :: aniso_active
        LOGICAL,       DIMENSION(2)                                     :: Us_present
        TYPE(aniso)                                                     :: Ujoint
        TYPE(matrix)                                                    :: UMAT
        TYPE(aniso),   DIMENSION(MAXGAU)                                :: U
        TYPE(matrix),  DIMENSION(MAXGAU)                                :: VTS
        REAL(KIND=wp), DIMENSION(MAXGAU)                                :: udet
        REAL(KIND=wp), DIMENSION(MAXGAU)                                :: aexp_qform
        REAL(KIND=wp)                                                   :: Biso
        REAL(KIND=wp), DIMENSION(3)                                     :: ueigen
        REAL(KIND=wp), DIMENSION(3,MAXGAU)                              :: VTS_duvwort
!       Aniso ADPs derivatives section:
!        REAL(KIND=wp), DIMENSION(6)                                     :: dRodU
!       Cofactors:
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: minors
        REAL(KIND=wp), DIMENSION(6)                                     :: minors_divided_by_det
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: det_scaled_deriv
!       Deriv. wrt aniso U
        REAL(KIND=wp), DIMENSION(6,6,MAXGAU)                            :: VTSD
        REAL(KIND=wp), DIMENSION(3,3)                                   :: VTSD_ROW
        REAL(KIND=wp), DIMENSION(3,1)                                   :: vtemp
        REAL(KIND=wp), DIMENSION(6,1)                                   :: utemp
        REAL(KIND=wp), DIMENSION(1,6)                                   :: wtemp
!       Mixed xxt vector with last three components doubled:
        REAL(KIND=wp), DIMENSION(6)                                     :: wxxt
!       Aniso matrices:
        TYPE(matrix)                                                    :: VXXT        
        REAL(KIND=wp), DIMENSION(6,6,6,MAXGAU)                          :: CTS
        REAL(KIND=wp), DIMENSION(6,6,MAXGAU)                            :: DZERO
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DSECOND
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DSECOND_DERIV
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DFOURTH
        REAL(KIND=wp), DIMENSION(6,6)                                   :: MINORS_SEC_DER
!       Additional matrices:
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DDT
        REAL(KIND=wp), DIMENSION(6,6)                                   :: VVT
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: VTSD_wxxt
!       Areas of normal matrix:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: XB
        REAL(KIND=wp), DIMENSION(6)                                     :: UO
        REAL(KIND=wp), DIMENSION(6)                                     :: UO_ALL_GAUSS
!       Symmetry:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: XU_SYM
        REAL(KIND=wp),       DIMENSION(3)                               :: BX_SYM
        REAL(KIND=wp),       DIMENSION(3)                               :: OX_SYM
        REAL(KIND=wp), DIMENSION(3,6)                                   :: UX_SYM
        REAL(KIND=wp), DIMENSION(6,6)                                   :: UU_SYM
        REAL(KIND=wp), DIMENSION(6,6)                                   :: BU_SYM
        REAL(KIND=wp), DIMENSION(6)                                     :: OU_SYM ! just 36 multiplications
!       Temp:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: TEMP36
        REAL(KIND=wp), DIMENSION(3,6)                                   :: TEMP36_ALL_GAUSS
        REAL(KIND=wp), DIMENSION(6,6)                                   :: UUMAT
        INTEGER                                                         :: ts
        INTEGER                                                         :: ij
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE                    :: RU 
!       Can handle 2 maps only:
        REAL(KIND=wp)                                                   :: map_value
        REAL(KIND=wp)                                                   :: occ
        REAL(KIND=wp)                                                   :: occ_i
        REAL(KIND=wp)                                                   :: occ_j
!       Maps variables:
        INTEGER                                                         :: number_of_maps
        INTEGER                                                         :: my_first_allocated_map
!       Refinement:        
        INTEGER                                                         :: natom
        INTEGER                                                         :: np
        LOGICAL,             DIMENSION(3,2)                             :: refine_atomic_xyz
        LOGICAL,             DIMENSION(2)                               :: refine_atomic_biso
        LOGICAL,             DIMENSION(2)                               :: refine_atomic_occ
        LOGICAL,             DIMENSION(6,2)                             :: refine_atomic_Us
        LOGICAL                                                         :: at_least_one_aniso
        LOGICAL                                                         :: ub_ref
        INTEGER                                                         :: indxyz
        INTEGER                                                         :: jndxyz
!       Counters:
        INTEGER                                                         :: iat
        INTEGER                                                         :: jat
        INTEGER                                                         :: k 
        INTEGER                                                         :: kk
        INTEGER                                                         :: l
        INTEGER                                                         :: run
!       Symmetry:
        TYPE(vector)                                                    :: xyz_iat
        TYPE(vector)                                                    :: xyz_jat
        TYPE(vector),        DIMENSION(2)                               :: xyz_frac
        REAL(KIND=wp),       DIMENSION(3,3)                             :: GMAT
        REAL(KIND=wp),       DIMENSION(3,3)                             :: XXT
!       Additional chain rule arrays:
        REAL(KIND=wp),       DIMENSION(:,:), ALLOCATABLE                :: AMAT
        REAL(KIND=wp),       DIMENSION(:,:), ALLOCATABLE                :: ASYM
        REAL(KIND=wp),       DIMENSION(3,3)                             :: BMAT
        REAL(KIND=wp),       DIMENSION(3,3)                             :: SMAT
!       General counters:
        INTEGER                                                         :: i
        INTEGER                                                         :: j
!       Symmetry counter: 
        INTEGER                                                         :: m
!       Non-zero elements counter:
        INTEGER                                                         :: atom_pair
        INTEGER                                                         :: elem
!       Additional squared pointer matrix to sparse matrix array:
        INTEGER,             DIMENSION(max_block_size,max_block_size)   :: table
        INTEGER                                                         :: block_size
!       Box counters:
        INTEGER                                                         :: ju
        INTEGER                                                         :: jv
        INTEGER                                                         :: jw
!       CPU monitoring:
        REAL(KIND=sp)                                                   :: time0
        REAL(KIND=sp)                                                   :: time1
!       Printing:
        INTEGER(KIND=eb)                                                :: modh
        INTEGER(KIND=eb)                                                :: pair
        CHARACTER(LEN=32),                                       SAVE   :: srname = 'sparse_normal_matrix'
!       Test:
        REAL(KIND=wp), DIMENSION(3)                                     :: xyz

!       Checkz: 
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die('Programming error. pdb_2 has not been initialized.', &
                     srname)
        ENDIF

!       This must be set from the start:
        my_first_allocated_map = 1

!       Check map allocation:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Programming error. Map array has not been allocated properly.', &
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( map_1(my_first_allocated_map) ) ) THEN
            CALL die('Programming error. MAP_1(1) has not been allocated properly.',&
                     srname)
        ELSE

            number_of_maps = SIZE ( map_1 )
            IF ( number_of_maps < 2 ) THEN
                WRITE(*,*) number_of_maps
                CALL die('Programming error. Need 2 maps....', srname)
            ENDIF

        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' NORMAL_MATRIX_DIAGONAL> ', 'natom=', I6)") natom
            CALL die ( 'Programming error. Need at least one atom for density generation.', srname)
        ENDIF


        IF ( H%nnz /= SIZE ( H%val ) ) THEN
            WRITE(*,*)  H%nnz, SIZE ( H%val )
            CALL die('Programming error. Matrix H has been incorrectly allocated.', srname)
        ENDIF

        CALL messag('Going to calculate '//TRIM ( int_to_c ( H%nnz ) )//' matrix elements... Wait...',&
                    srname)

!       Check whether aniso records present:
        aniso_active = ALLOCATED ( pdb_2%U )

!       Initialize sparse matrix array:
        H%val = 0.0_wp
        H%ia  = 0
        H%ja  = 0

!       Initialize pointer to atomic pair block:
        nnz   = 0

!       Allocate additional matrices for chain rule:
        CALL allocate_array(ASYM, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops, 3)
        CALL allocate_array(AMAT, 3, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)
        IF ( aniso_active ) THEN
            CALL allocate_array ( RU,  map_1(my_first_allocated_map)%sp_group%number_of_symops)
        ENDIF

!       Build 3Nsym x 3 matrix applying chain rule (see notes above):
        DO m = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops

!           Calculate position of the block:  
            l = 3 * (m - 1)
!                                        T 
!           Calculating transpose [Asym ] :
            BMAT = map_1(my_first_allocated_map)%ORT                    &
                 * .SYMA. map_1(my_first_allocated_map)%sp_group%SYM(m) &
                 *        map_1(my_first_allocated_map)%DEORT

!           To avoid Intel compiler warnings (in debugging mode):
            ASYM(l+1:l+3,1:3) = BMAT
            IF ( aniso_active ) THEN
                RU(m) =  .SIXBYSIX. BMAT
            ENDIF
        ENDDO

!       Initialize 3x3 diagonal matrix:
        GMAT = 0.0_wp
        DO l = 1, 3
            GMAT(l,l) = 2.0_wp
        ENDDO

!       Other params:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Calculate maximum sphere radius around any atom:
        ort_diag    = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Set for statistics:
        max_r_cut_off = -10.0_wp ** 4 
        min_r_cut_off =  10.0_wp ** 4
        
!       Timing:
        CALL CPU_TIME(time0)
        modh = MAX ( SIZE ( H%val ) / 100, 1 )

!       Counter for matrix printing:
        pair  = 0
 
!       Prepare matrix pointers:
!        CALL allocate_array(ind, SIZE ( pdb_2 ) )

!       Need to check this:
!        ind(1) = 0
!        DO iat = 2, natom 
!            ind(iat) = ind(iat-1) + COUNT ( pdb_2%refinable(1:max_block_size, iat-1) )
!        ENDDO

!        IF ( debug > 40 ) THEN
!            WRITE(*,*) ' ind=', ind
!        ENDIF

!       Loop through all atoms:
        coord:DO atom_pair = 1, SIZE ( pair_list, DIM=2 )

            iat = pair_list(1, atom_pair)
            jat = pair_list(2, atom_pair)

!           Check jat, iat order:
            IF ( jat < iat ) THEN
                WRITE(*,*) iat, jat
                CALL warn('Incorrect indices. JAT < IAT. Swapping...', srname)
                CALL swap (IAT, JAT)
                WRITE(*,*) iat, jat
                CALL sleep(1)
            ENDIF

!           Check whether it is worthwhile to calculate derivatives for this block:
            IF ( COUNT ( pdb_2%refinable(1:max_block_size, jat ) ) == 0 .OR. jat < iat ) THEN
                IF ( jat < iat ) THEN
                    CALL die('O-o-o-o-p-s-s... Swap does not work.', srname)
                ENDIF  
                CYCLE
            ENDIF

            occ_i = pdb_2%occ(iat)
            occ_j = pdb_2%occ(jat)

!           Joint atom constants:
            occ = occ_i * occ_j
!           Joint biso (always calculated):
            Biso = pdb_2%biso(iat) + pdb_2%biso(jat)


!           Prepare additional indexing table:
            CALL compose_atom_block_table (table, block_size, pdb_2, iat, jat)

            Us_present = .FALSE.
            at_least_one_aniso = .FALSE.

            IF ( aniso_active ) THEN
                Us_present(1) = ANY ( pdb_2%U(iat) > 0.0_wp ) 
                Us_present(2) = ANY ( pdb_2%U(jat) > 0.0_wp )
                at_least_one_aniso = Us_present(1) .OR. Us_present(2)
            ENDIF

            CALL compose_joint_atom (ajo, bjo, ngauss, pdb_2, b_scale, iat, jat)

            IF ( .NOT. at_least_one_aniso ) THEN

!                WRITE(*,"( ' bjo=', 25F6.1)") bjo
!               Increment joint gaussian B constants by joint Biso (see above). Had to interchange 2 lines below.
!               BUG corrected MAR 2009 BVS:
                bjo(1:ngauss) = bjo(1:ngauss) + Biso 
                a(1:ngauss)  = ajo(1:ngauss) * ( SQRT ( pi4 / bjo(1:ngauss) ) ) ** 3

!               Final joint atom constants (switching  back to usual names):
                be(1:ngauss) = 4.0_wp * pi_squared / bjo(1:ngauss)


                IF ( debug > 42 ) THEN
                    WRITE(*,*) ' ngauss=', ngauss
                    WRITE(*,"(' a =', 25F8.4)") a(1:ngauss)
                    WRITE(*,"(' be=', 25F8.4)") be(1:ngauss)
                ENDIF

!               Array to speed up:
                bbb(1:ngauss) = 1.0_wp / bjo(1:ngauss)

                IF ( debug > 42 ) THEN
                    WRITE(*,"(' bbb=', 15F8.4)") bbb(1:ngauss)
                ENDIF 

!               Superior floating radius, gives excellent results for RNAse in case of high B-values:
                r2_cut_off = integration_radius ( ajo(1:ngauss), bjo(1:ngauss) ) ** 2

                IF ( debug > 42 ) THEN
                    WRITE(*,*) ' optimal diag r_cut_off=', SQRT ( r2_cut_off )
                ENDIF

!               Make sure sphere around any atom fits into the cell:
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Accumulate statistics for R_CUT_OFF:
                max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
                min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )         

                IF ( debug > 42 ) THEN
                    WRITE(*,*) ' my_first_allocated_map=', my_first_allocated_map
                    WRITE(*,*) ' diag r_cut_off=', r_cut_off
                ENDIF

!               Figure out box size around joint atom:
                box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

            ENDIF

!           Figure out refinable parameters for atom pair:
            DO k = 1, 2
                l = pair_list(k,atom_pair)
                refine_atomic_xyz(1:3,k)  = pdb_2%refinable(1:3,l)
                refine_atomic_biso(k)     = pdb_2%refinable(4,l)
                refine_atomic_occ(k)      = pdb_2%refinable(5,l)
                refine_atomic_Us(1:6,k)   = pdb_2%refinable(6:11,l)
            ENDDO

!           Check whether we have pure aniso or mixed aniso/b refinement:
            IF ( at_least_one_aniso ) THEN
                ub_ref = ( ANY ( refine_atomic_Us(1:6,1 ) .AND. ANY ( refine_atomic_Us(1:6,2) ) ) ) .OR. &
                         ( ANY ( refine_atomic_Us(1:6,1 ) .AND. refine_atomic_biso(2) ) ) .OR.           &
                         ( ANY ( refine_atomic_Us(1:6,2 ) .AND. refine_atomic_biso(1) ) ) 
            ENDIF
          
!           Consult IND arrat to get matrix pointers:
            indxyz = pdb_2%indxyz(iat)
            jndxyz = pdb_2%indxyz(jat)

! -----------> STANDING HERE
                
!           Calculate H%IA, H%JA indices in advance (coordinate format):
            CALL assign_indices (H, table, pdb_2, nnz, indxyz, jndxyz, iat, jat)

!           Special arrays:
            IF ( ANY ( refine_atomic_biso ) ) THEN
!               Speedup array:
                bbb2(1:ngauss) = bbb(1:ngauss) ** 2
            ENDIF

!           Special debugging:
            IF ( debug > 22 ) THEN
                WRITE(*,*) iat, indxyz, refine_atomic_xyz
                WRITE(*,*) jat, jndxyz, pdb_2%refinable(1:3,jat)
            ENDIF

!           Convert to fractions:
            xyz_frac = (/ map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(iat), &
                          map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(jat) /)
         
            xyz_iat = xyz_frac(1)

!           Initialize AMAT for chain rule:
            AMAT = 0.0_wp

!           Loop over all symmetry operators:
            symmetry:DO m = 1, map_1(my_first_allocated_map)%sp_group%number_of_symops

!               This will cause some asymmetry in calculations (see below):
                xyz_jat = map_1(my_first_allocated_map)%sp_group%SYM(m) * xyz_frac(2)

                aniso_symmetry:IF ( at_least_one_aniso ) THEN

!                   COnsider a few cases how to deal with mixed refinement:
                    IF ( Us_present(1) .AND. Us_present(2) ) THEN
!                       Need to rotate:
                        Ujoint = pdb_2%U(iat) + RU(m) * pdb_2%U(jat)
                    ELSE IF ( Us_present(1) .AND. .NOT. Us_present(2) ) THEN
!                       No need to rotate:
                        Ujoint = pdb_2%U(iat) + pdb_2%biso(jat) / ( 8.0_wp * pi_squared )
                    ELSE IF ( .NOT. Us_present(1) .AND. Us_present(2) ) THEN
!                       Need to rotate:
                        Ujoint = RU(m) * pdb_2%U(jat) + pdb_2%biso(iat) / ( 8.0_wp * pi_squared )
                    ENDIF

!                   Figure out effective Biso of joint atom:
                    UMAT = Ujoint
                    ueigen = eigen_values ( UMAT )
                    Biso =  8.0_wp * pi_squared * MAXVAL ( ueigen )

!                   Superior floating radius just in right place:
                    r2_cut_off = estimate_integration_radius ( MAXVAL ( bjo(1:ngauss) ) + Biso, xaccuracy = 0.1_wp ** 7 ) ** 2
                    r2_cut_off = MIN  ( r2_cut_off, r2_cell_max )
                    r_cut_off  = SQRT ( r2_cut_off )

!                   Accumulate statistics for R_CUT_OFF:
                    max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
                    min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )

!                   Figure out box size around joint atom:
                    box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!                   Prepare all aniso crap that does not depend on coordinates (i.e. DUVWORT vector):
                    CTS = 0.0_wp
                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian (b_scale has been included in bjo:
                        U(l) = Ujoint + ( bjo(l) / ( 8.0_wp * pi_squared ) )

!                       Calculate minors for each aniso gauss:
                        minors(1:6,l) = minors_aniso ( U(l) )

!                       Calculate determinant:
                        Udet(l) = udet_aniso ( U(l) )

!                       Almost VTS matrix:
                        minors_divided_by_det = minors(1:6,l) / udet(l)

!                       Convert to VTS matrix (U**-1):
                        VTS(l) = .CONVERT. minors_divided_by_det
!                        CALL print_mat(VTS(l)) 

!                       U derivative section:
                        DO k = 1, 6
                            VTSD(k,1:6,l) = VTS_first_deriv ( U(l), k )
                        ENDDO

                        det_scaled_deriv(1:6,l) = det_first_deriv_aniso ( U(l) ) / UDET(l)
                        utemp(1:6,1) = det_scaled_deriv(1:6,l)
                        wtemp(1,1:6) = det_scaled_deriv(1:6,l)
                        VVT = MATMUL ( utemp, TRANSPOSE ( utemp ) )                  ! VVT
                        DSECOND_DERIV = det_second_deriv_mat (U(l) ) / udet(l)       ! D2U
                        DZERO(1:6,1:6,l) = 1.5_wp * VVT - DSECOND_DERIV

                        formula_24:DO ts = 1, 6

!                           Derivative of all U(ij) against single TS minor:
                            utemp(1:6,1) = minors_first_deriv_test(U(l), ts)
                            DDT          = MATMUL ( utemp, wtemp )
                            MINORS_SEC_DER = minors_second_deriv ( ts )        ! M2D

!                           Accumulate CTS for a given gaussian:
                            CTS(1:6,1:6,ts,l) =  CTS(1:6,1:6,ts,l) + ( minors_divided_by_det(ts)         &
                                              * (-3.0_wp * VVT + DSECOND_DERIV)                          &
                                              + (1.0_wp / udet(l) )                                      &
                                              * ( 1.5_wp * ( DDT + TRANSPOSE (DDT) ) - MINORS_SEC_DER ) )

                        ENDDO formula_24

                   ENDDO

!                  Scale constants in advance to reduce computation:
                   a(1:ngauss) = ajo(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF aniso_symmetry
               
!               Prepare symmetry matrix for chain rule for certain derivatives:
                l    = 3 * ( m - 1)
                SMAT = TRANSPOSE ( ASYM(l+1:l+3,1:3) )

!               Initialising coordinate matrix for each symmetry operator to build even larger matrix:
                BMAT = 0.0_wp

!               Small arrays to speed up chain rule application for BX,BY,BZ etc. derivatives:
                BX_SYM = 0.0_wp
                OX_SYM = 0.0_wp 
                XU_SYM = 0.0_wp
                UX_SYM = 0.0_wp
                UU_SYM = 0.0_wp
                BU_SYM = 0.0_wp
                OU_SYM = 0.0_wp

!               Calculate joint atom position for H1 term:
                Hterms:DO run = 1, 2

!                   Better not to run cycle against MY_FIRST_ALLOCATED_MAP:
                    my_first_allocated_map = run

!                   Second map is used to obtain H2 terms:
                    H2 = my_first_allocated_map == 2

!                   Figure out how to calculate centre of joint atom:
                    IF ( .NOT. H2 ) THEN
                        centre_of_atom_in_fractions = xyz_iat - xyz_jat
                    ELSE
                        centre_of_atom_in_fractions = xyz_iat + xyz_jat
                    ENDIF

                    IF ( debug > 22 ) THEN

                        xyz = centre_of_atom_in_fractions
                        IF ( .NOT. H2 ) THEN
                            WRITE(*,"( 'H1 centre= ', 3F8.5)") xyz
                        ELSE
                            WRITE(*,"( 'H2 centre= ', 3F8.5)") xyz
                        ENDIF

                    ENDIF

!                   Figure out closest grid point:
                    closest_grid_point = NINT ( centre_of_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )

!                   Box limits:
                    box_lower_limit = closest_grid_point - box_size
                    box_upper_limit = closest_grid_point + box_size

!                   Create a box around the joint atom and calculate its density and convolute with the map:
                    boxzz:DO jw = box_lower_limit(3), box_upper_limit(3)

                        DO jv = box_lower_limit(2), box_upper_limit(2)

                            DO ju = box_lower_limit(1), box_upper_limit(1)

                                uvw = (/ ju, jv, jw /)
                                duvw = REAL ( uvw ) / REAL( map_1(my_first_allocated_map)%map_size) &
                                                         - centre_of_atom_in_fractions

!                               High point of electron density convolution:
                                r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw   ! corresponds to r - r(n) + r(m)

                                rcut2:IF ( r2 <= r2_cut_off ) THEN

                                    jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                    map_value = map_1(my_first_allocated_map)%array(jut(1), jut(2), jut(3))

                                    aniso_business:IF ( at_least_one_aniso ) THEN

!                                       Note change of sign corresponing to (r(i) - r) as convolution theorem suggests:
                                        duvwort = -map_1(my_first_allocated_map)%ORT * duvw

!                                       Prepare major arrays depending on DUVWORT:
                                        IF ( ub_ref ) THEN 

                                            wxxt(1) = duvwort(1) ** 2
                                            wxxt(2) = duvwort(2) ** 2
                                            wxxt(3) = duvwort(3) ** 2
                                            wxxt(4) = 2.0_wp * duvwort(1) * duvwort(2)
                                            wxxt(5) = 2.0_wp * duvwort(1) * duvwort(3)
                                            wxxt(6) = 2.0_wp * duvwort(2) * duvwort(3)

                                        ENDIF

                                        DO l = 1, ngauss
                                            VTS_duvwort(1:3,l) = VTS(l) * duvwort
                                            aexp_qform(l) = a(l) * EXP ( -0.5_wp * DOT_PRODUCT ( duvwort, VTS_duvwort(1:3,l) ) )
                                            IF ( ub_ref ) VTSD_wxxt(1:6,l) = MATMUL ( VTSD(1:6,1:6,l), wxxt )
                                        ENDDO


!                                       Seems the right place to calculate this since all constants are in place:
                                        TEMP36_ALL_GAUSS = 0
                                        DO l = 1, ngauss
                                            utemp(1:6,1) = det_scaled_deriv(1:6,l) + VTSD_wxxt(1:6,l) 
                                            vtemp(1:3,1) = VTS_duvwort(1:3,l)
!                                           Initialise:
                                            TEMP36 = 0.5_wp * MATMUL ( vtemp, TRANSPOSE ( utemp ) )
!                                           Accumulate:
                                            DO ij = 1, 6
!                                               Convert to squared matrix:
                                                VTSD_ROW = .CONVERT. VTSD(ij,1:6,l)
                                                TEMP36(1:3,ij) = TEMP36(1:3,ij) - MATMUL ( VTSD_ROW, duvwort )

                                            ENDDO
                                            TEMP36_ALL_GAUSS = TEMP36_ALL_GAUSS + TEMP36 * aexp_qform(l)
                                        ENDDO

                                        IF ( ub_ref ) THEN
                                            UUMAT = 0

                                            DO l = 1, ngauss
                                                DSECOND = 0.0_wp
                                                DO ts = 1, 6
!                                                   THis is 256 multiplications * ngauss
                                                    DSECOND = DSECOND +  CTS(1:6,1:6,ts,l) * wxxt(ts)
                                                ENDDO
                                                utemp(1:6,1) = VTSD_wxxt(1:6,l)

!                                               Fourth order coefs:
                                                DFOURTH = 0.5_wp * MATMUL ( utemp, TRANSPOSE ( utemp ) )

                                                UUMAT = UUMAT + (DZERO(1:6,1:6,l) + DSECOND + DFOURTH) * aexp_qform(l)
                                            ENDDO

                                            UUMAT = ( 0.5_wp * map_value * occ ) * UUMAT ! correct common scale aniso

                                        ENDIF

                                        UO_ALL_GAUSS = 0
                                        DO l = 1, ngauss
                                            utemp(1:6,1) = det_scaled_deriv(1:6,l) + VTSD_wxxt(1:6,l) 
                                            UO_ALL_GAUSS = UO_ALL_GAUSS - utemp(1:6,1) * aexp_qform(l)
                                        ENDDO

                                        UO_ALL_GAUSS = 0.5_wp * map_value * UO_ALL_GAUSS


                                        any_xyz1_aniso:IF ( ANY ( refine_atomic_xyz(1:3,1 ) ) ) THEN

                                            XX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
!                                               ---------------------------------------------------------
!                                               XX AREA
                                                DO l = 1, ngauss
                                                    vtemp(1:3,1) = VTS_duvwort(1:3,l)
                                                    VXXT = MATMUL ( vtemp, TRANSPOSE ( vtemp ) )

!                                                   Note change of sign (tupi * i * h) ** 2 results in minus sign...
                                                    BMAT = BMAT - occ * map_value * (VXXT - VTS(l)) * aexp_qform(l) 
                                                ENDDO

                                            ENDIF XX_area

!                                           Check whether Biso refinable for JAT atom:
                                            XB_area:IF ( refine_atomic_biso(2) ) THEN

!                                               Accumulate matrix values wrt XB+YB+ZB:
                                                XB = TEMP36_ALL_GAUSS

!                                               Average along rows (just first 3 elements for each row):
                                                DO k = 1, 3
                                                    gg(k) = SUM ( XB(k,1:3) ) / ( 8.0_wp * pi_squared )
                                                ENDDO

!                                               Switch sign for H2 terms:
                                                IF ( H2 ) gg = -gg
                                                
!                                               Select:
                                                DO k = 1, 3
                                                    IF ( refine_atomic_xyz(k,1) ) THEN
                                                        elem = table(k,4)
                                                        IF ( elem > 0 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + occ * map_value * gg(k)
                                                        ENDIF
                                                    ENDIF
                                                ENDDO

                                            ENDIF XB_area 

!                                           Accumulate derivatives in XO area fro XYZ rows:
                                            XO_area:IF ( refine_atomic_occ(2) ) THEN
!                                                JAT atom seems to have refinable occupancy:
                                                 gg = -occ_i * map_value &
                                                    * MATMUL ( VTS_duvwort(1:3,1:ngauss), aexp_qform(1:ngauss) ) !BUG

!                                                Switch sign for H2 terms:
                                                 IF ( H2 ) gg = -gg
                                                 DO k = 1, 3
                                                     IF ( refine_atomic_xyz(k,1) ) THEN
                                                         elem = table(k,5)
                                                         IF ( elem > 0 ) THEN
                                                             H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                         ENDIF
                                                     ENDIF
                                                 ENDDO
                                            ENDIF XO_area

!                                           Entering XU area:
                                                      
                                            XU_area:IF ( ANY ( refine_atomic_us(1:6,2) ) )  THEN

                                                IF ( .NOT. H2 ) THEN
                                                    XU_SYM = XU_SYM + (occ * map_value ) * TEMP36_ALL_GAUSS
                                                ELSE
                                                    XU_SYM = XU_SYM - (occ * map_value ) * TEMP36_ALL_GAUSS
                                                ENDIF

                                            ENDIF XU_area

                                        ENDIF any_xyz1_aniso

!                                       ======= Biso row =======
                                        biso1:IF ( refine_atomic_biso(1) ) THEN

                                            BX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
                                                
!                                               Average along rows (just 3 elements for each row):
                                                DO k = 1, 3
                                                    BX_SYM(k) = BX_SYM(k) - occ * map_value &
                                                              * SUM ( TEMP36_ALL_GAUSS(k,1:3) ) / (8.0_wp * pi_squared)
                                                ENDDO
!                                               This will be transformed by symmetry later...

                                            ENDIF BX_area

!                                           Entering dummy BB area:
!                                           Check whether JAT atom has refinable Biso:
                                            BB_area:IF ( refine_atomic_biso(2) ) THEN
!                                                See isotropic area...
                                                 WRITE(*,"(' iat=',I6,' jat=',I6, ' table=',L1, ' at_least_one_aniso=',L1)")&
                                                             iat, jat, table(4,4), at_least_one_aniso
                                                 WRITE(*,"(' U(iat)=', 6F10.6)") pdb_2%U(iat)%u
                                                 WRITE(*,"(' U(jat)=', 6F10.6)") pdb_2%U(jat)%u
                                                 CALL die('O-o-o-p-s-s. Cannot get in here... Logic is wrong.', srname)
                                            ENDIF BB_area
                              
!                                           BO area: B isotropic, Oj atom has anizo values...

!                                           Need averaged gradients against U
!                                           Let's pretend we have UO refinement:
                                            BO_area:IF ( refine_atomic_occ(2) ) THEN

                                                elem = table(4,5)
!                                                UO = 0                                         

                                                IF ( elem > 0 ) THEN

!                                                    UO = occ_i * UO_ALL_GAUSS
!                                                    UO = occ_j *  MATMUL (  RU(m)%a, UO_ALL_GAUSS )
                                                    UO = occ_j * ( RU(m) * UO_ALL_GAUSS )

                                                    IF ( .NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + SUM ( UO(1:3) ) / ( 8 * pi ** 2 )
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - SUM ( UO(1:3) ) / ( 8 * pi ** 2 )
                                                    ENDIF
                                                ENDIF

                                            ENDIF BO_area


!                                           Check JAT atom for anisotropy refinement:
                                            BU_area:IF ( ANY ( refine_atomic_Us(2,1:6) ) ) THEN

!                                               Pretend we have UU refinement, then average. We need (1,1:6) vector.
!                                               Infamous UU 6x6 area:
                                               
                                                IF ( .NOT. H2 ) THEN
                                                    BU_SYM = BU_SYM + UUMAT
                                                ELSE
                                                    BU_SYM = BU_SYM - UUMAT
                                                ENDIF

                                            ENDIF BU_area

                                        ENDIF biso1

!                                       Finished with Biso row.

!                                       ===== Occupancy row start ======

                                        occ1:IF ( refine_atomic_occ(1 )  ) THEN

!                                           OX 1x3 area. Check JAT atom:
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                                Note "+" change of sign here:
                                                 OX_SYM = OX_SYM + occ_j * map_value &
                                                        * MATMUL ( VTS_duvwort(1:3,1:ngauss), aexp_qform(1:ngauss) )

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            OB_area:IF ( refine_atomic_biso(2) ) THEN

                                                elem = table(5,4)

!                                               Pretend it's OU area and average:
                                                IF ( elem > 0 ) THEN

                                                    UO = occ_j *  UO_ALL_GAUSS

                                                    IF (.NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + SUM( UO(1:3) ) / (8.0_wp * pi_squared)
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - SUM( UO(1:3) ) / (8.0_wp * pi_squared)
                                                    ENDIF    

                                                ENDIF

                                            ENDIF OB_area

!                                           Entering OO area:
                                            OO_area:IF ( refine_atomic_occ(2) ) THEN

                                                elem = table(5,5)
                                                IF ( elem > 0 ) THEN
                                                    IF (.NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + map_value * SUM ( aexp_qform(1:ngauss) )
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - map_value * SUM ( aexp_qform(1:ngauss) )
                                                    ENDIF 
                                                ENDIF
                                                
                                            ENDIF OO_area

!                                           Entering OU area:
                                            OU_area: iF ( any ( refine_atomic_us(1:6,2) ) )THEN

                                                IF ( .NOT. H2 ) THEN
                                                    OU_SYM = OU_SYM + occ_j * UO_ALL_GAUSS
                                                ELSE
                                                    OU_SYM = OU_SYM - occ_j * UO_ALL_GAUSS
                                                ENDIF

                                            ENDIF OU_area

                                        ENDIF occ1


!                                       ======= U row ======

                                        any_us1:IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN 

                                            UX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN 

!                                               Accumulate for all grid points for a given symmetry operator:
                                                UX_SYM = UX_SYM - (occ * map_value) * TEMP36_ALL_GAUSS

                                            ENDIF UX_area

!                                           Entering UB 6x1 area:
                                            UB_area:IF ( refine_atomic_biso(2) ) THEN

!                                               Pretend we have UU area, then average:

!                                               DO some averaging:
                                                DO k = 1, 6

!                                                   6:11,4 area:
                                                    elem = table(5+k,4)

                                                    IF ( elem > 0 ) THEN

                                                        IF ( .NOT. H2 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + SUM ( UUMAT(k,1:3) ) &
                                                                            / (8.0_wp * pi_squared)
                                                        ELSE
                                                            H%val(nnz+elem) = H%val(nnz+elem) - SUM ( UUMAT(k,1:3) ) &
                                                                            / (8.0_wp * pi_squared)
                                                        ENDIF

                                                    ENDIF

                                                ENDDO


                                            ENDIF UB_area

!                                           Entering 6x1 UO area:
                                            UO_area:IF ( refine_atomic_occ(2) ) THEN

!                                               It's important to have occ_i here, occ_j vanishes:
                                                UO =  occ_i * UO_ALL_GAUSS

                                                DO k = 1,6

                                                    elem = table(5+k, 5)

                                                    IF ( elem > 0 ) THEN

                                                        IF ( .NOT. H2 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + UO(k)
                                                        ELSE
                                                            H%val(nnz+elem) = H%val(nnz+elem) - UO(k)     
                                                        ENDIF

                                                     ENDIF
                                                ENDDO

                                                
                                            ENDIF UO_area

!                                           Entering infamous UU area:
                                            UU_area: IF ( ANY ( refine_atomic_Us(1:6,2)  ) ) THEN

                                                IF ( .NOT. H2 ) THEN
                                                    UU_SYM = UU_SYM + UUMAT
                                                ELSE
                                                    UU_SYM = UU_SYM - UUMAT
                                                ENDIF

                                            ENDIF UU_area

                                        ENDIF any_us1

!                                       ====== aniso business is over ======

                                    ELSE aniso_business

!                                       ======= Isotrpopic business =======
                                    
!                                       Going to speed up in case of Biso/occncy refinement:
                                        ber2(1:ngauss)   = be(1:ngauss) * r2
                                        exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                        abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)
                                        abexp_be(1:ngauss) = abexp(1:ngauss) * be(1:ngauss)


!                                       Special array for speed-up:
                                        IF ( ANY ( refine_atomic_biso ) ) THEN
                                                sum_bbb_52 = SUM ( abexp_be(1:ngauss) * bbb(1:ngauss) &
                                                           * (5.0_wp / 2.0_wp - ber2(1:ngauss)) )

                                            IF ( ( refine_atomic_biso(1) .AND. refine_atomic_occ(2) ) .OR. &
                                                ( refine_atomic_biso(2) .AND. refine_atomic_occ(1) ) ) THEN
                                                 sum_bbb_15 = SUM ( abexp(1:ngauss) * bbb(1:ngauss) * (-1.5_wp + ber2(1:ngauss)) )
                                            ENDIF

                                        ENDIF

!                                       Additional array for speed up:
                                        IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                           This constant is associated with XYZ refinement only:
                                            sum_abexp_be = SUM ( abexp_be(1:ngauss) )

!                                           Convert to orthogonal system of coodinates (easier to calculate derivatives):
                                            duvwort = 2.0_wp * map_1(my_first_allocated_map)%ORT * duvw
 
                                            anyxyz:IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

!                                               Makes sense to calculate if JAT atom has refinable XYZ:
                                                IF (  ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                                   Calculate Tronrud XXT matrix:
                                                    DO l = 1, 3
                                                        DO k = 1, 3
!                                                           DUVWORT has been scaled by a factor of 2 already:
                                                            XXT(l,k) = duvwort(k) * duvwort(l)
                                                        ENDDO
                                                    ENDDO

!                                                   Build 3x3 matrix for orthogonal coordinates:
                                                    BMAT = BMAT + map_value * occ * (sum_abexp_be * GMAT           &
                                                                - SUM ( abexp_be(1:ngauss) * be(1:ngauss) ) * XXT)

                                                ENDIF

!                                               Check whether Biso refinable for JAT atom:
                                                IF ( refine_atomic_biso(2) ) THEN

!                                                   Accumulate matrix values wrt XB+YB+ZB(14):
                                                    gg = -(occ * map_value * sum_bbb_52) * duvwort

!                                                   Switch sign for H2 terms:
                                                    IF ( H2 ) gg = -gg

!                                                   General case -> not every XYZ(iat) is refinable:                                            
                                                    IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

                                                        DO k = 1, 3

                                                            IF ( refine_atomic_xyz(k,1) ) THEN

!                                                               HMAT(indxyz+l,jndbiso):
                                                                elem = table(k,4)

                                                                IF ( elem > 0 ) THEN
                                                                    H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                                ENDIF        

                                                            ENDIF

                                                        ENDDO

                                                    ENDIF

                                                ENDIF

!                                               Accumulate derivatives XO, YO, ZO(15) for XYZ rows:
                                                IF ( refine_atomic_occ(2) ) THEN

!                                                   JAT atom seems to have refinable occupancy:
                                                    gg = occ_i * map_value * sum_abexp_be * duvwort

!                                                   Switch sign for H2 terms:
                                                    IF ( H2 ) gg = -gg

!                                                   Considering general case -> not every XYZ(iat) is refinable:
                                                    IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN                                                

                                                        DO k = 1, 3
                                                            IF ( refine_atomic_xyz(k,1) ) THEN

!                                                               HMAT(indxyz+l,jndocc) = HMAT(indxyz+l,jndocc) + gg(k)
                                                                elem = table(k,5)

                                                                IF ( elem > 0 ) THEN
                                                                    H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                                ENDIF

                                                            ENDIF
                                                        ENDDO

                                                    ENDIF

                                                ENDIF

                                            ENDIF anyxyz
                                        ENDIF


!                                   -------------- Biso row ---------------------------

                                        iso1:IF ( refine_atomic_biso(1) ) THEN

!                                           FIXME:
!                                           These derivatives seem to be in lower triangular in some circumstances:

!                                           Accumulate derivatives for BX, BY, BZ(14):
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) .AND. ANY ( table(4,1:3) > 0 ) ) THEN
                                                BX_SYM = BX_SYM + (occ * map_value * sum_bbb_52) * duvwort
                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Calculate pointer to position in sparse matrix:
                                                elem = table(4,4)
                                                
!                                               Accumulate BB matrix elements (12):
                                                IF ( .NOT. H2 ) THEN

!                                                   HMAT(indbiso,jndbiso):
                                                    H%val(nnz+elem) = H%val(nnz+elem) + occ * map_value           &
                                                                    * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)     &
                                                                    * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp) &
                                                                    + 15.0_wp / 4.0_wp ) )

                                                ELSE

                                                    H%val(nnz+elem) = H%val(nnz+elem) - occ * map_value           &
                                                                    * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)     &
                                                                    * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp) &
                                                                    + 15.0_wp / 4.0_wp ) )
                                               
                                                ENDIF

                                            ENDIF                                                    

!                                           Check whether JAT atom has refinable occ:
                                            IF ( refine_atomic_occ(2) ) THEN

!                                              Calculate pointer:
                                               elem = table(4,5)

!                                              Check that this is upper triangular part of the matrix:
                                               IF ( elem > 0 ) THEN

!                                                   Accumulate BQ matrix elements (16):
                                                    IF ( .NOT. H2 ) THEN
!                                                       HMAT(indbiso,jndocc):
                                                        H%val(nnz+elem) = H%val(nnz+elem) + occ_i * map_value * sum_bbb_15 
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - occ_i * map_value * sum_bbb_15
                                                    ENDIF

                                                ENDIF

                                            ENDIF

                                        ENDIF iso1  
                           
!                                       ------------- Now occupancy row itself, note OCC_J is used --------------:
                                        occncy:IF ( refine_atomic_occ(1) ) THEN

!                                           FIXME these derivatives maybe in lower triangular part in some circumstances:
!                                           The overhead calculation is just for atoms on the diagonal though:

!                                           Check whether JAT atom has refinable XYZ:
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
!                                                WRITE(*,*) ' Str.pl.at_least_one_aniso=', at_least_one_aniso,' iat=',iat,' jat=',jat
!                                               OX+OY+OZ interaction (15):
                                                OX_SYM = OX_SYM - (occ_j * map_value * sum_abexp_be) * duvwort

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Calculate pointer:
                                                elem = table(5,4)

!                                               This will implicitly indicate that iat /= jat:
                                                IF ( elem > 0 ) THEN

                                                
!                                                   QB interaction(16), note plus sign here:
                                                    IF ( .NOT. H2 ) THEN
!                                                       HMAT(indocc, jndbiso):
                                                        H%val(nnz+elem) = H%val(nnz+elem) + occ_j * map_value * sum_bbb_15
                                                    ELSE
!                                                       Switch sign for H2 terms:
                                                        H%val(nnz+elem) = H%val(nnz+elem) - occ_j * map_value * sum_bbb_15
                                                    ENDIF

                                                ENDIF

                                            ENDIF

!                                           Check whether JAT atom has refinable occncy and apply(13):
                                            IF ( refine_atomic_occ(2) ) THEN
!                                               Calculate pointer:
                                                elem = table(5,5)

!                                               Switch sign for H2 terms:
                                                IF ( .NOT. H2 ) THEN
!                                                   HMAT(indocc,jndocc):
                                                    H%val(nnz+elem) = H%val(nnz+elem) + map_value * SUM ( abexp(1:ngauss) )
                                                ELSE
                                                    H%val(nnz+elem) = H%val(nnz+elem) - map_value * SUM ( abexp(1:ngauss) )
                                                ENDIF

                                            ENDIF

                                        

                                        ENDIF occncy


                                    ENDIF aniso_business


                                ENDIF rcut2

                            ENDDO
                        ENDDO
                    ENDDO boxzz

                ENDDO Hterms

                IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) )  THEN

!                   Build 3x3Nsym matrix:
                    AMAT(1:3,3*(m-1)+1:3*(m-1)+3) = BMAT
                ENDIF

!               BUG CORRECTED JAN 2008 BVS (CANNOT MERGE WITH PREVIOUS IF):
                IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                   Take care of chain rule for BXYZ_DERIV & OXYZ_DERIV:
                    IF ( refine_atomic_biso(1) ) THEN
                        BX_SYM = MATMUL ( SMAT, BX_SYM )
                    ENDIF

                    IF ( refine_atomic_occ(1) ) THEN
                        OX_SYM = MATMUL ( SMAT, OX_SYM )
                    ENDIF                      

!                   Finish with chain rule for two classes of second_derivatives:
                    IF ( refine_atomic_biso(1) ) THEN

!                       Consider general case -> not every XYZ(jat) is refinable:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,2) ) THEN

!                               Calculate pointer:
                                elem = table(4,k)

!                               Upper triangular only:
                                IF ( elem > 0 ) THEN

!                                   HMAT(indbiso,jndxyz+l):
                                    H%val(nnz+elem) = H%val(nnz+elem) + BX_SYM(k)
                                ENDIF

                            ENDIF

                        ENDDO

                    ENDIF

                    IF ( refine_atomic_occ(1) ) THEN

!                       Consider general case -> not every XYZ(jat) is refinable:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,2) ) THEN

!                               Calculate pointer:
                                elem = table(5,k)

!                               Upper triangular:
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = H%val(nnz+elem) + OX_SYM(k)
                                ENDIF

                            ENDIF

                        ENDDO
                    ENDIF
                ENDIF

!               ==========
!               Aniso area:
!               ==========
                XU_sym_area:IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN
            
!                       XU_SYM contains accumulated XU matrix from all grid points in atom sphere:
                        CALL aniso_symop_matmul(XU_SYM, XU_SYM, RU(m))

!                       Select:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,1) ) THEN
                                DO ij = 1, 6
                                    elem = table(k,5+ij)
                                    IF ( elem > 0 ) THEN                                                            
                                        H%val(nnz+elem) = H%val(nnz+elem) + XU_SYM(k,ij)
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO

                    ENDIF

                ENDIF XU_sym_area

                Ux_sym_area:IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN
                    IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
                        UX_SYM = MATMUL ( SMAT, UX_SYM )
                   
        
                        DO ij = 1, 6
                            IF ( refine_atomic_Us(ij,1) ) THEN
                                DO k = 1, 3
                                    IF ( refine_atomic_xyz(k,2) ) THEN
                                        elem = table(5+ij,k)
                                        IF ( elem > 0 ) THEN
                                            H%val(nnz+elem) = H%val(nnz+elem) + UX_SYM(k,ij)
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF Ux_sym_area


                UU_sym_area: IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN

                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN
                        CALL aniso_symop_matmul(UU_SYM, UU_SYM, RU(m))
                        DO k = 1, 6
                            DO ij = 1, 6
                                elem = table(k+5,ij+5)
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = H%val(nnz+elem) + UU_SYM(k,ij)
                                ENDIF
                            ENDDO
                        ENDDO
                    ENDIF

                ENDIF UU_sym_area

                OU_sym_area: IF ( refine_atomic_occ(1) ) THEN
                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN

!                        OU_SYM = MATMUL ( OU_SYM, RU(m)%a )
                        OU_SYM = OU_SYM * RU(m)
                        DO k = 1,6
                            elem = table(5,5+k)
                            IF ( elem > 0 ) THEN
                                H%val(nnz+elem) = H%val(nnz+elem) + OU_SYM(k)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF OU_sym_area


                BU_sym_area: IF ( refine_atomic_biso(1) ) THEN
                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN

                        CALL aniso_symop_matmul(BU_SYM, RU(m), BU_SYM)
!                       Do some averaging:
                        DO k = 1, 6
                            elem = table(4,5+k)
                            IF ( elem > 0 ) THEN

!                               JUL 2009 BVS BUG corrected (changed order of indices in BU_SYM during summation: 
                                H%val(nnz+elem) = H%val(nnz+elem) + SUM ( BU_SYM(1:3,k) ) &
                                                                  / (8.0_wp * pi_squared)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF BU_sym_area

!               A small caveat here -> my_first_allocated_map will have value of "2" in outer loops.

            ENDDO symmetry

            chainrule:IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!               Finish with XYZ chain rule for this atom, ASYM is on the right side:
                BMAT = MATMUL ( AMAT, ASYM )

!               Extract and add refinable items only:
                DO k = 1, 3
                    IF ( refine_atomic_xyz(k,1) ) THEN
                        DO kk = 1, 3
                            IF ( refine_atomic_xyz(kk,2) ) THEN
!                               Calculate pointer to HMAT(indxyz+l,jndxyz+ll):
                                elem = table(k,kk)
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = BMAT(k,kk)
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO

!               Debug if necessary:
                debugg:IF ( debug > 25 ) THEN

                    IF ( map_1(my_first_allocated_map)%sp_group%number_of_symops <= 24 ) THEN

                        DO k = 1, 3
                            WRITE(*,"(' AMAT=',72ES9.2)") &
                            (AMAT(k,l), l = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops)
                        ENDDO

                    ENDIF
   
                    DO k = 1, 3 * map_1(my_first_allocated_map)%sp_group%number_of_symops
                        WRITE(*,"(' ASYM=', 3F9.5)") (ASYM(k,l), l = 1, 3)
                    ENDDO

                        DO k = 1, 3
                            WRITE(*,"(' BMAT=', 3ES9.2)") (BMAT(k,l), l = 1, 3)
                        ENDDO

                ENDIF debugg

            ENDIF chainrule

!           Monitor CPU:
            monitor:DO k = 1, max_block_size 

                DO kk = 1, max_block_size

!                   Calculate pointer:
                    elem = table(k,kk)

                    refinables:IF ( elem > 0 ) THEN
                        pair = pair + 1
                        IF ( MOD ( pair, modh ) == 0 ) THEN
                            CALL CPU_TIME ( time1 )
                            WRITE(*,"(' SPARSE_NORMAL_MATRIX> ',                                               &
                           &'matrix element[', A, ',', A, ']=', T62, ES12.5, 2X,'CPU time=',F8.1, ' s', I10)") &
                             TRIM ( int_to_c ( H%ia(nnz+elem ) ) ), TRIM ( int_to_c ( H%ja(nnz+elem ) ) ),     &
                             H%val(nnz+elem), time1 - time0, pair
                        ENDIF

                        IF ( H%ia(nnz+elem) == 0 .OR. H%ja(nnz+elem) == 0 ) THEN
                            WRITE(*,*) ' k=  ', k,   ' kk=', kk
                            WRITE(*,*) ' nnz=', nnz, ' elem=',elem
                            WRITE(*,*) ' H%ia(nnz+elem)=',  H%ia(nnz+elem)
                            WRITE(*,*) ' H%ja(nnz+elem)=',  H%ja(nnz+elem)
                        ENDIF

!                       Check diagonal matrix elements for non-positive values:
                        IF ( H%ia(nnz+elem) ==  H%ja(nnz+elem) .AND. H%val(nnz+elem) <= 0.0_wp ) THEN
                            WRITE(*,*) ' k=  ', k,   ' kk=', kk
                            WRITE(*,*) ' nnz=', nnz, ' elem=',elem
                            WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'matrix element=', 2I8, ' HMAT=', ES12.5, &
                           &' iat/jat=', 2I6, ' ref. params=', 2I3, ' nnz+elem=',I6)") &
                            H%ia(nnz+elem), H%ja(nnz+elem),  H%val(nnz+elem) ,iat, jat, k, kk, nnz+elem

                            WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Biso=', ES9.2, ' occ=', ES9.2, ' r_cut_off=', F8.2)") &
                            pdb_2%biso(iat), pdb_2%occ(iat), r_cut_off

!                           Bailing out:
                            CALL warn('Programming error. Zero or negative normal matrix diagonal element. &
                                     &Refinement halted. Sorry.', srname)

                        ENDIF


                    ENDIF refinables

                ENDDO

            ENDDO monitor

!           Increment non-zero element pointer:
            nnz = nnz + block_size

        ENDDO coord

!       Report integration radii stats:
        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'R_CUT_OFF: min=', F8.1, ' max=', F8.1, ' max possible=', F8.1)")&
        min_r_cut_off, max_r_cut_off, SQRT ( r2_cell_max )

!       Free memory:
        CALL deallocate_array(AMAT)
        CALL deallocate_array(ASYM)

!       Expand in coordinate matrix fmt:
        CALL messag('Doing in-place matrix expansion...', srname)

!       Count refinable parameters:

!       Number of diagonal elements before expansion:
        np = COUNT ( pdb_2%refinable )

        elem = 0
        DO l = 1, nnz
            i = H%ia(l)
            j = H%ja(l)            
            IF ( i == j ) THEN
                elem = elem + 1
                IF ( elem <= np ) THEN
                    diag(elem) = H%val(l)
                ELSE
                    CALL die('Programming error. Too many diagonal elements have been calculated.',&
                             srname)
                ENDIF              
            ENDIF 
        ENDDO

        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Current number of calculated nnz elements=           ', A)" ) &
        TRIM ( int_to_c ( nnz ) )
        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Calculated number of diagonal elements=              ', A)" ) & 
        TRIM ( int_to_c ( elem ) )
        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Exact wanted number of diagonal elements=            ', A)" ) &
        TRIM ( int_to_c ( np ) )
        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Expected number of additional non-diagonal elements= ', A)" ) &
        TRIM ( int_to_c ( nnz - np ) )
        WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 'Total expected number of elements=                   ', A)" ) &
        TRIM ( int_to_c ( 2 * nnz - np ) )

!       Additional checkz before expansion:
        IF ( np /= elem ) THEN
            WRITE(*,*) np, elem
            CALL warn('Programming error. Wrong number of diagonal elements has been calculated.', srname)
        ENDIF

        IF ( 2 * nnz - np /= H%nnz ) THEN
            WRITE(*,*) 2 * nnz - np, H%nnz
            CALL warn('Programming error. Inconsitent number of elements has been calculated.',&
                      srname)
        ENDIF

        IF ( debug > 40 ) THEN
            WRITE(*,*) ' testing nnz=', nnz
        ENDIF

    END SUBROUTINE sparse_normal_matrix

    SUBROUTINE fast_sparse_normal_matrix(H, map_1, pdb_2, pair_list, nnz, ind)
!
!       Purpose:
!       =======
!
!       Calculates sparse matrix H (H1 + H2 terms) with respect to
!       X,Y,Z,B and Q params.
!
!       Supply two maps a la Tronrud and you are done. Just one map will result in
!       calculation of H1 terms only.
!
!       Date:                 Programmer:             History of changes:
!       ====                  ==========              ==================
!       January 10, 2008      B.Strokopytov           Original code
!       January 15, 2008      B.Strokopytov           Added BXYZ_DERIV and OXYZ_DERIV arrays for 
!                                                     chain rule speed-up.
!
!       MAR 25, 2009             -"-                  nnz has been changed to INOUT
!                                                     H matrix will zeroized from the outside
!                                                     nnz will be zeroized from the outside
!                
!       MAY 2010                                      Components of RU(m) are no longer
!                                                     exposed due to calls to ANISO_SYMOP_MATMUL
!                                                     BLAS calls maybe used inside module ANISO_SYMOP_MANIP. 
!
!       Chain rule note:
!       ===============
!       In orthogonal coordinate system we have:
!
!                     T      -1    -1 -1 -1       -1
!       [Asym] = [OAD] = [OAD] = [D  A  O  ] = [OA  D]
!
!       Obviously:
!             T
!       [Asym] = [OAD] 
!       
!       Transpose is due to chain rule, if we going to multiply from the right:
!         2    2      T
!       [d R/dx][Asym]
!
!       Remaining problems:
!       ==================
!       How to get rid of the symmetry loop? Looks like a horrible CPU problem
!       in cubic space groups. 
!
        TYPE(coo_matrix),                                 INTENT(INOUT) :: H
        TYPE(map),           DIMENSION(:),   ALLOCATABLE, INTENT(IN)    :: map_1
        TYPE(pdb),                                        INTENT(IN)    :: pdb_2
        INTEGER,             DIMENSION(:,:),              INTENT(IN)    :: pair_list
        INTEGER(KIND=eb),                                 INTENT(INOUT) :: nnz
        INTEGER(KIND=eb),    DIMENSION(:),                INTENT(IN)    :: ind
!       Local variables:
        LOGICAL                                                         :: h2
        TYPE(vector)                                                    :: centre_of_atom_in_fractions
        TYPE(vector_int)                                                :: closest_grid_point
        TYPE(vector_int)                                                :: box_size
        INTEGER,             DIMENSION(3)                               :: jut
        INTEGER,             DIMENSION(3)                               :: box_lower_limit
        INTEGER,             DIMENSION(3)                               :: box_upper_limit
        TYPE(vector)                                                    :: deort_diag
        TYPE(vector)                                                    :: ort_diag
!       Normal atom params:
        INTEGER                                                         :: ngauss
        REAL(KIND=wp)                                                   :: b_scale
!       Scattering:
!       For diagonal terms 15 gaussians would be the MAX number,
!       but for non-diagonal terms need 25 (looking into the future):
        INTEGER,             PARAMETER                                  :: MAXGAU = 25
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: a
!       Joint atom:
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: ajo 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bjo
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: be
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bbb 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: bbb2 
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: ber2
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: exp_be
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: abexp
        REAL(KIND=wp),       DIMENSION(MAXGAU)                          :: abexp_be
        REAL(KIND=wp)                                                   :: sum_bbb_15
        REAL(KIND=wp)                                                   :: sum_bbb_52
        REAL(KIND=wp)                                                   :: sum_abexp_be
        REAL(KIND=wp),       DIMENSION(3)                               :: gg
!       Distance:
        REAL(KIND=wp)                                                   :: r2
        TYPE(vector_int)                                                :: uvw
        TYPE(vector)                                                    :: duvw
        REAL(KIND=wp),       DIMENSION(3)                               :: duvwort
!       Cut-off:
        REAL(KIND=wp)                                                   :: r2_cut_off
        REAL(KIND=wp)                                                   :: r_cut_off
        REAL(KIND=wp),                                          SAVE    :: max_r_cut_off  ! Don't touch this for Christ sake.
        REAL(KIND=wp),                                          SAVE    :: min_r_cut_off
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                                   :: r2_cell_max
!       Aniso section:
        LOGICAL                                                         :: aniso_active
        LOGICAL,       DIMENSION(2)                                     :: Us_present
        TYPE(aniso)                                                     :: Ujoint
        TYPE(matrix)                                                    :: UMAT
        TYPE(aniso),   DIMENSION(MAXGAU)                                :: U
        TYPE(matrix),  DIMENSION(MAXGAU)                                :: VTS
        REAL(KIND=wp), DIMENSION(MAXGAU)                                :: udet
        REAL(KIND=wp), DIMENSION(MAXGAU)                                :: aexp_qform
        REAL(KIND=wp)                                                   :: Biso
        REAL(KIND=wp), DIMENSION(MAXGAU)                                :: Uiso
        TYPE(aniso)                                                     :: U_iat
        TYPE(aniso)                                                     :: U_jat
        REAL(KIND=wp)                                                   :: biso_iat
        REAL(KIND=wp)                                                   :: biso_jat
        REAL(KIND=wp), DIMENSION(3)                                     :: ueigen
        REAL(KIND=wp), DIMENSION(3,MAXGAU)                              :: VTS_duvwort
!       Aniso ADPs derivatives section:
!        REAL(KIND=wp), DIMENSION(6)                                     :: dRodU
!       Cofactors:
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: minors
        REAL(KIND=wp), DIMENSION(6)                                     :: minors_divided_by_det
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: det_scaled_deriv
!       Deriv. wrt aniso U
        REAL(KIND=wp), DIMENSION(6,6,MAXGAU)                            :: VTSD
        REAL(KIND=wp), DIMENSION(3,3)                                   :: VTSD_ROW
        REAL(KIND=wp), DIMENSION(3,1)                                   :: vtemp
        REAL(KIND=wp), DIMENSION(6,1)                                   :: utemp
        REAL(KIND=wp), DIMENSION(1,6)                                   :: wtemp
!       Mixed xxt vector with last three components doubled:
        REAL(KIND=wp), DIMENSION(6)                                     :: wxxt
!       Aniso matrices:
        TYPE(matrix)                                                    :: VXXT        
        REAL(KIND=wp), DIMENSION(6,6,6,MAXGAU)                          :: CTS
        REAL(KIND=wp), DIMENSION(6,6,MAXGAU)                            :: DZERO
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DSECOND
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DSECOND_DERIV
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DFOURTH
        REAL(KIND=wp), DIMENSION(6,6)                                   :: MINORS_SEC_DER
!       Additional matrices:
        REAL(KIND=wp), DIMENSION(6,6)                                   :: DDT
        REAL(KIND=wp), DIMENSION(6,6)                                   :: VVT
        REAL(KIND=wp), DIMENSION(6,MAXGAU)                              :: VTSD_wxxt
!       Areas of normal matrix:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: XB
        REAL(KIND=wp), DIMENSION(6)                                     :: UO
        REAL(KIND=wp), DIMENSION(6)                                     :: UO_ALL_GAUSS
!       Symmetry:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: XU_SYM
        REAL(KIND=wp),       DIMENSION(3)                               :: BX_SYM
        REAL(KIND=wp),       DIMENSION(3)                               :: OX_SYM
        REAL(KIND=wp), DIMENSION(3,6)                                   :: UX_SYM
        REAL(KIND=wp), DIMENSION(6,6)                                   :: UU_SYM
        REAL(KIND=wp), DIMENSION(6,6)                                   :: BU_SYM
        REAL(KIND=wp), DIMENSION(6)                                     :: OU_SYM ! just 36 multiplications
!       Temp:
        REAL(KIND=wp), DIMENSION(3,6)                                   :: TEMP36
        REAL(KIND=wp), DIMENSION(3,6)                                   :: TEMP36_ALL_GAUSS
        REAL(KIND=wp), DIMENSION(6,6)                                   :: UUMAT
        INTEGER                                                         :: ts
        INTEGER                                                         :: ij
!       Aniso symmetry:
        TYPE(aniso_symop), DIMENSION(:), ALLOCATABLE                    :: RU 
!       Can handle 2 maps only:
        REAL(KIND=wp)                                                   :: map_value
        REAL(KIND=wp)                                                   :: occ
        REAL(KIND=wp)                                                   :: occ_i
        REAL(KIND=wp)                                                   :: occ_j
!       Maps variables:
        INTEGER                                                         :: number_of_maps
        INTEGER                                                         :: my_first_allocated_map
!       Refinement:        
!        INTEGER                                                         :: natom
        INTEGER                                                         :: np
        LOGICAL,             DIMENSION(3,2)                             :: refine_atomic_xyz
        LOGICAL,             DIMENSION(2)                               :: refine_atomic_biso
        LOGICAL,             DIMENSION(2)                               :: refine_atomic_occ
        LOGICAL,             DIMENSION(6,2)                             :: refine_atomic_Us
        LOGICAL                                                         :: at_least_one_aniso
        LOGICAL                                                         :: ub_ref
!       Pointers:
        INTEGER                                                         :: indxyz
        INTEGER                                                         :: jndxyz
!       Counters:
        INTEGER                                                         :: iat
        INTEGER                                                         :: jat
        INTEGER                                                         :: k 
        INTEGER                                                         :: kk
        INTEGER                                                         :: l
        INTEGER                                                         :: run
!       Symmetry:
        TYPE(vector)                                                    :: xyz_iat
        TYPE(vector)                                                    :: xyz_jat
        TYPE(vector),        DIMENSION(2)                               :: xyz_frac
        REAL(KIND=wp),       DIMENSION(3,3)                             :: GMAT
        REAL(KIND=wp),       DIMENSION(3,3)                             :: XXT
!       Additional chain rule arrays:
        REAL(KIND=wp),       DIMENSION(:,:), ALLOCATABLE                :: AMAT
        REAL(KIND=wp),       DIMENSION(:,:), ALLOCATABLE                :: ASYM
        REAL(KIND=wp),       DIMENSION(3,3)                             :: BMAT
        REAL(KIND=wp),       DIMENSION(3,3)                             :: SMAT
!       Symmetry counter: 
        INTEGER                                                         :: m
        INTEGER                                                         :: nsym
!       Non-zero elements counter:
        INTEGER                                                         :: atom_pair
        INTEGER                                                         :: elem
        INTEGER                                                         :: number_of_atom_pairs
!       Additional squared pointer matrix to sparse matrix array:
        INTEGER,             DIMENSION(max_block_size,max_block_size)   :: table
        INTEGER                                                         :: block_size
!       Box counters:
        INTEGER                                                         :: ju
        INTEGER                                                         :: jv
        INTEGER                                                         :: jw
!       CPU monitoring:
        REAL(KIND=sp),                                          SAVE    :: time0
        REAL(KIND=sp)                                                   :: time1
!       Printing:
        INTEGER(KIND=eb)                                                :: modh
        INTEGER(KIND=eb),                                        SAVE   :: pair
        INTEGER                                                         :: OMP_GET_THREAD_NUM
        CHARACTER(LEN=32),                                       SAVE   :: srname = 'fast_sparse_normal_matrix'
        LOGICAL                                                         :: output_stats
!       Test:
        REAL(KIND=wp), DIMENSION(3)                                     :: xyz

        number_of_atom_pairs =  SIZE ( pair_list, DIM=2 )
!       Return in case of an empty list:
        IF ( number_of_atom_pairs == 0 ) THEN
            RETURN
        ENDIF

!       Checkz: 
        IF ( .NOT. ALLOCATED ( pdb_2 ) ) THEN
            CALL die('Programming error. pdb_2 has not been initialized.', &
                     srname)
        ENDIF

!       This must be set from the start:
        my_first_allocated_map = 1

!       Check map allocation:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Programming error. Map array has not been allocated properly.', &
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( map_1(my_first_allocated_map) ) ) THEN
            CALL die('Programming error. MAP_1(1) has not been allocated properly.',&
                     srname)
        ELSE

            number_of_maps = SIZE ( map_1 )
            IF ( number_of_maps < 2 ) THEN
                WRITE(*,*) number_of_maps
                CALL die('Programming error. Need 2 maps....', srname)
            ENDIF

        ENDIF

        IF ( H%nnz /= SIZE ( H%val ) ) THEN
            WRITE(*,*)  H%nnz, SIZE ( H%val )
            CALL die('Programming error. Matrix H has been incorrectly allocated.', srname)
        ENDIF


        IF ( nnz == 0 ) THEN
            np = COUNT ( pdb_2%refinable )
            CALL messag('Going to calculate '//TRIM ( int_to_c ( H%nnz ) ) // '/' //&
                         &TRIM ( int_to_c ( (H%nnz + np) / 2 ) ) // ' matrix elements... Wait...',srname)
        ENDIF

!       Check whether aniso records present:
        aniso_active = ALLOCATED ( pdb_2%U )

!       Allocate additional matrices for chain rule:
        nsym = map_1(my_first_allocated_map)%sp_group%number_of_symops
        IF ( nsym < 1 .OR. nsym > 192 ) THEN
            WRITE(*,*) nsym
            CALL die('Programming error. Unreasonable value of NSYM.', srname)
        ENDIF
        CALL allocate_array(ASYM, 3 * nsym, 3)
        CALL allocate_array(AMAT, 3, 3 * nsym)
        IF ( aniso_active ) THEN
            CALL allocate_array ( RU,  map_1(my_first_allocated_map)%sp_group%number_of_symops)
        ENDIF

!       Build 3Nsym x 3 matrix applying chain rule (see notes above):
        DO m = 1, nsym 

!           Calculate position of the block:  
            l = 3 * (m - 1)
!                                        T 
!           Calculating transpose [Asym ] :
            BMAT = map_1(my_first_allocated_map)%ORT                    &
                 * .SYMA. map_1(my_first_allocated_map)%sp_group%SYM(m) &
                 *        map_1(my_first_allocated_map)%DEORT

!           To avoid Intel compiler warnings (in debugging mode):
            ASYM(l+1:l+3,1:3) = BMAT

            IF ( aniso_active ) THEN
!               To avoid Intel compiler warnings (in debugging mode) do this in two steps:
                DDT = .SIXBYSIX. BMAT
                RU(m) = DDT
            ENDIF

        ENDDO

!       Initialize 3x3 diagonal matrix:
        GMAT = 0.0_wp
        DO l = 1, 3
            GMAT(l,l) = 2.0_wp
        ENDDO

!       Other params:
        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Calculate maximum sphere radius around any atom:
        ort_diag    = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Set for statistics:
        IF ( nnz == 0 ) THEN
            max_r_cut_off = -10.0_wp ** 4 
            min_r_cut_off =  10.0_wp ** 4
            pair  = 0
            CALL CPU_TIME(time0)
        ENDIF
        output_stats = nnz < 0

!       Timing:
        modh = MAX ( SIZE ( H%val ) / 40, 1 )

!       Loop through all atoms:
!$OMP PARALLEL DO DEFAULT (PRIVATE) SHARED ( H, MAP_1, PDB_2, NUMBER_OF_MAPS, MODH, PAIR_LIST, PAIR, NSYM, SRNAME, &
!$OMP                                        IND, R2_CELL_MAX, B_SCALE, DEORT_DIAG, ORT_DIAG,  ASYM, GMAT, DEBUG,  &
!$OMP                                        MAX_R_CUT_OFF, MIN_R_CUT_OFF, TIME0, RU, ANISO_ACTIVE, NUMBER_OF_ATOM_PAIRS )
        coord:DO atom_pair = 1, number_of_atom_pairs

            iat = pair_list(1, atom_pair)
            jat = pair_list(2, atom_pair)

!           Much better to have it private:
            my_first_allocated_map = 1
!            WRITE(*,*) ' iat,jat=', iat,jat, my_first_allocated_map

!           Check jat, iat order:
            IF ( jat < iat ) THEN
                WRITE(*,*) iat, jat
                CALL warn('Incorrect indices. JAT < IAT. Swapping...', srname)
                CALL swap (IAT, JAT)
                WRITE(*,*) iat, jat
                CALL sleep(1)
            ENDIF

!           Check whether it is worthwhile to calculate derivatives for this block:
            IF ( COUNT ( pdb_2%refinable(1:11, jat ) ) == 0 .OR. jat < iat ) THEN
                IF ( jat < iat ) THEN
                    CALL die('O-o-o-o-p-s-s... Swap does not work.', srname)
                ENDIF  
                CYCLE
            ENDIF

            occ_i = pdb_2%occ(iat)
            occ_j = pdb_2%occ(jat)

!           Joint atom constants:
            occ = occ_i * occ_j

            biso_iat = pdb_2%biso(iat)
            biso_jat = pdb_2%biso(jat)


!           Joint biso (always calculated):
            Biso = biso_iat + biso_jat


!           Prepare additional indexing table:
            CALL compose_atom_block_table (table, block_size, pdb_2, iat, jat)

            Us_present = .FALSE.
            at_least_one_aniso = .FALSE.

            IF ( aniso_active ) THEN
                Us_present(1) = ANY ( pdb_2%U(iat) > 0.0_wp )
                Us_present(2) = ANY ( pdb_2%U(jat) > 0.0_wp )
                at_least_one_aniso = Us_present(1) .OR. Us_present(2)
                iF ( Us_present(1) ) U_iat = pdb_2%U(iat)
                IF ( Us_present(2) ) U_jat = pdb_2%U(jat)
            ENDIF

!            Reset number of gaussians due to new algorithm:
             ngauss = 1
             bjo    = b_scale
             ajo    = 0
             ajo(1) = 1.0_wp

            IF ( .NOT. at_least_one_aniso ) THEN

!               Increment joint gaussian B constants by joint Biso (see above). Had to interchange 2 lines below.
!               BUG corrected MAR 2009 BVS:
                bjo(1:ngauss) = bjo(1:ngauss) + Biso 
                a(1:ngauss)  = ajo(1:ngauss) * ( SQRT ( pi4 / bjo(1:ngauss) ) ) ** 3

!               Final joint atom constants (switching  back to usual names):
                be(1:ngauss) = 4.0_wp * pi_squared / bjo(1:ngauss)


                IF ( debug > 42 ) THEN
                    WRITE(*,*) ' ngauss=', ngauss
                    WRITE(*,"(' a =', 25F8.4)") a(1:ngauss)
                    WRITE(*,"(' be=', 25F8.4)") be(1:ngauss)
                ENDIF

!               Array to speed up:
                bbb(1:ngauss) = 1.0_wp / bjo(1:ngauss)

                IF ( debug > 42 ) THEN
                    WRITE(*,"(' bbb=', 15F8.4)") bbb(1:ngauss)
                ENDIF 

!               Superior floating radius, gives excellent results for RNAse in case of high B-values:
                r2_cut_off = estimate_integration_radius ( MAXVAL ( bjo(1:ngauss) ), xaccuracy = 0.1_wp ** 8 ) ** 2

                IF ( r2_cut_off == 0.0_wp ) THEN
                    WRITE(*,*) ' r2_cut_off=', r2_cut_off
                    CALL die('Programming error. R2_CUT_OFF has unreasonable value.', srname)
                ENDIF

                IF ( debug > 142 ) THEN
                    WRITE(*,*) ' optimal diag r_cut_off=', SQRT ( r2_cut_off )
                ENDIF

!               Make sure sphere around any atom fits into the cell:
!$OMP CRITICAL
                r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
                r_cut_off  = SQRT ( r2_cut_off )

!               Accumulate statistics for R_CUT_OFF:
                max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
                min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )         
!$OMP END CRITICAL
                IF ( debug > 42 ) THEN
                    WRITE(*,*) ' my_first_allocated_map=', my_first_allocated_map, OMP_GET_THREAD_NUM()
                    WRITE(*,*) ' diag r_cut_off=', r_cut_off
                ENDIF

!               Figure out box size around joint atom:
                box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

            ENDIF

!           Figure out refinable parameters for atom pair:
            DO k = 1, 2
                l = pair_list(k,atom_pair)
                refine_atomic_xyz(1:3,k)  = pdb_2%refinable(1:3,l)
                refine_atomic_biso(k)     = pdb_2%refinable(4,l)
                refine_atomic_occ(k)      = pdb_2%refinable(5,l)
                refine_atomic_Us(1:6,k)   = pdb_2%refinable(6:11,l)
            ENDDO

!           Check whether we have pure aniso or mixed aniso/b refinement:
            IF ( at_least_one_aniso ) THEN
!               Maybe this is not enterily correct since ref of occupancies may also require similar arrays:
                ub_ref = ( ANY ( refine_atomic_Us(1:6,1 ) .AND. ANY ( refine_atomic_Us(1:6,2) ) ) ) .OR. &
                         ( ANY ( refine_atomic_Us(1:6,1 ) .AND. refine_atomic_biso(2) ) ) .OR.           &
                         ( ANY ( refine_atomic_Us(1:6,2 ) .AND. refine_atomic_biso(1) ) ) 
            ENDIF
          
!           Consult INDXYZ array to get matrix pointers:
            indxyz = pdb_2%indxyz(iat)
            jndxyz = pdb_2%indxyz(jat)

!           Calculate H%IA, H%JA indices in advance (coordinate format):
            nnz = ind(atom_pair)
            CALL assign_indices (H, table, pdb_2, nnz, indxyz, jndxyz, iat, jat)

!           Special arrays:
            IF ( ANY ( refine_atomic_biso ) ) THEN
!               Speedup array:
                bbb2(1:ngauss) = bbb(1:ngauss) ** 2
            ENDIF

!           Special debugging:
            IF ( debug > 22 ) THEN
                WRITE(*,*) iat, jat,  OMP_GET_THREAD_NUM()
                WRITE(*,*) ' r_cut=', r_cut_off,  OMP_GET_THREAD_NUM()
!                WRITE(*,*) jat, jndxyz, pdb_2%refinable(1:3,jat)
            ENDIF

!           Convert to fractions (array allocation takes place) and be careful:
!$OMP CRITICAL
            xyz_frac(1) = map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(iat)
            xyz_frac(2) = map_1(my_first_allocated_map)%DEORT * pdb_2%xyz(jat)
!$OMP END CRITICAL
            xyz_iat = xyz_frac(1)

!           Initialize AMAT for chain rule (please don't touch $OMP CRITICAL):
!$OMP CRITICAL
            IF ( .NOT. ALLOCATED ( AMAT ) ) THEN
                CALL allocate_array(AMAT, 3, 3 * nsym)
            ENDIF

            AMAT = 0.0_wp
!$OMP END CRITICAL

!           Loop over all symmetry operators:
            symmetry:DO m = 1, nsym

!               This will cause some asymmetry in calculations (see below):
                xyz_jat = map_1(my_first_allocated_map)%sp_group%SYM(m) * xyz_frac(2)

                aniso_symmetry:IF ( at_least_one_aniso ) THEN

!                   Let's consider a few cases how to deal with mixed refinement:
! Especially nasty section with respect to parallelization:
!$OMP CRITICAL
                    IF ( Us_present(1) .AND. Us_present(2) ) THEN
!                       Need to rotate:
                        Ujoint = U_iat + RU(m) * U_jat
                    ELSE IF ( Us_present(1) .AND. .NOT. Us_present(2) ) THEN
!                       No need to rotate (add_u_to_b overloaded):
                        Ujoint = U_iat + biso_jat / ( 8.0_wp * pi_squared )
                    ELSE IF ( .NOT. Us_present(1) .AND. Us_present(2) ) THEN
!                       Need to rotate second ellipsoid:
                        Ujoint = biso_iat / ( 8.0_wp * pi_squared ) + RU(m) * U_jat
                    ENDIF

!                   Figure out effective Biso of joint atom:
                    UMAT = Ujoint
                    ueigen = eigen_values ( UMAT )   ! temporary created
                    Biso =  8.0_wp * pi_squared * MAXVAL ( ueigen )
                    Uiso(1:ngauss) =  bjo(1:ngauss) + Biso ! temporary created
!$OMP END CRITICAL
!                   Superior floating radius just in right place:
                    r2_cut_off = estimate_integration_radius ( MAXVAL ( Uiso(1:ngauss) ), xaccuracy = 0.1_wp ** 8 ) ** 2

!                   since ajo=1, ngauss=1 which result in signigicant accuracy loss...
!                   This is entirely incorrect BUG CORRECTED May 2009 BVS:
!                    r2_cut_off = integration_radius ( ajo(1:ngauss), Uiso(1:ngauss) ) ** 2
!$OMP CRITICAL
                    r2_cut_off = MIN  ( r2_cut_off, r2_cell_max )
!$OMP END CRITICAL
                    r_cut_off  = SQRT ( r2_cut_off )

!                   Accumulate statistics for R_CUT_OFF:
!$OMP CRITICAL
                    max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
                    min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )
!$OMP END CRITICAL
!                   Figure out box size around joint atom:
                    box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!                   Prepare all aniso crap that does not depend on coordinates (i.e. DUVWORT vector):
                    CTS = 0.0_wp
                    DO l = 1, ngauss

!                       Calculate array of anisotropic values for each gaussian (B_SCALE has been included in BJO):
                        U(l) = Ujoint + ( bjo(l) / ( 8.0_wp * pi_squared ) )

!                       Calculate minors for each aniso gauss:
                        minors(1:6,l) = minors_aniso ( U(l) )

!                       Calculate determinant:
                        Udet(l) = udet_aniso ( U(l) )

!                       Almost VTS matrix:
                        minors_divided_by_det = minors(1:6,l) / udet(l)

!                       Convert to VTS matrix (U**-1):
                        VTS(l) = .CONVERT. minors_divided_by_det
!                        CALL print_mat(VTS(l)) 

!                       U derivative section:
                        DO k = 1, 6
                            VTSD(k,1:6,l) = VTS_first_deriv ( U(l), k )
                        ENDDO

                        det_scaled_deriv(1:6,l) = det_first_deriv_aniso ( U(l) ) / UDET(l)
                        utemp(1:6,1) = det_scaled_deriv(1:6,l)
                        wtemp(1,1:6) = det_scaled_deriv(1:6,l)
                        VVT = MATMUL ( utemp, TRANSPOSE ( utemp ) )                  ! VVT
                        DSECOND_DERIV = det_second_deriv_mat (U(l) ) / udet(l)       ! D2U
                        DZERO(1:6,1:6,l) = 1.5_wp * VVT - DSECOND_DERIV

                        formula_24:DO ts = 1, 6

!                           Derivative of all U(ij) against single TS minor:
                            utemp(1:6,1) = minors_first_deriv_test(U(l), ts)
                            DDT          = MATMUL ( utemp, wtemp )
                            MINORS_SEC_DER = minors_second_deriv ( ts )        ! M2D

!                           Accumulate CTS for a given gaussian:
                            CTS(1:6,1:6,ts,l) =  CTS(1:6,1:6,ts,l) + ( minors_divided_by_det(ts)         &
                                              * (-3.0_wp * VVT + DSECOND_DERIV)                          &
                                              + (1.0_wp / udet(l) )                                      &
                                              * ( 1.5_wp * ( DDT + TRANSPOSE (DDT) ) - MINORS_SEC_DER ) )

                        ENDDO formula_24

                   ENDDO

!                  Scale constants in advance to reduce computation:
                   a(1:ngauss) = ajo(1:ngauss) / ( twopi ** 1.5_wp * SQRT ( udet(1:ngauss ) ) )

                ENDIF aniso_symmetry
               
!               Prepare symmetry matrix for chain rule for certain derivatives:
                l    = 3 * ( m - 1)
                SMAT = TRANSPOSE ( ASYM(l+1:l+3,1:3) )

!               Initialising coordinate matrix for each symmetry operator to build even larger matrix:
                BMAT = 0.0_wp

!               Small arrays to speed up chain rule application for BX,BY,BZ etc. derivatives:
                BX_SYM = 0.0_wp
                OX_SYM = 0.0_wp 
                XU_SYM = 0.0_wp
                UX_SYM = 0.0_wp
                UU_SYM = 0.0_wp
                BU_SYM = 0.0_wp
                OU_SYM = 0.0_wp

!               Calculate joint atom position for H1 term:
                Hterms:DO run = 1, 2

!                   Better not to run cycle against MY_FIRST_ALLOCATED_MAP:
                    my_first_allocated_map = run

!                   Second map is used to obtain H2 terms:
                    H2 = my_first_allocated_map == 2

!                   Figure out how to calculate centre of joint atom:
                    IF ( .NOT. H2 ) THEN
                        centre_of_atom_in_fractions = xyz_iat - xyz_jat
                    ELSE
                        centre_of_atom_in_fractions = xyz_iat + xyz_jat
                    ENDIF
                    IF ( debug > 22 ) THEN

                        xyz = centre_of_atom_in_fractions
                        IF ( .NOT. H2 ) THEN
                            WRITE(*,"( 'H1 centre= ', 3F8.5)") xyz
                        ELSE
                            WRITE(*,"( 'H2 centre= ', 3F8.5)") xyz
                        ENDIF

                    ENDIF

!                   Figure out closest grid point:
                    closest_grid_point = NINT ( centre_of_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )

!                   Box limits:
                    box_lower_limit = closest_grid_point - box_size
                    box_upper_limit = closest_grid_point + box_size

!                   Create a box around the joint atom and calculate its density and convolute with the map:
                    boxzz:DO jw = box_lower_limit(3), box_upper_limit(3)

                        DO jv = box_lower_limit(2), box_upper_limit(2)

                            DO ju = box_lower_limit(1), box_upper_limit(1)

                                uvw = (/ ju, jv, jw /)
                                duvw = REAL ( uvw ) / REAL( map_1(my_first_allocated_map)%map_size) &
                                                         - centre_of_atom_in_fractions

!                               High point of electron density convolution:
                                r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw   ! corresponds to r - r(n) + r(m)

                                rcut2:IF ( r2 <= r2_cut_off ) THEN

                                    jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                    map_value = map_1(my_first_allocated_map)%array(jut(1), jut(2), jut(3))

                                    aniso_business:IF ( at_least_one_aniso ) THEN

!                                       Note change of sign corresponing to (r(i) - r) as convolution theorem suggests:
                                        duvwort = -map_1(my_first_allocated_map)%ORT * duvw

!                                       Prepare major arrays depending on DUVWORT:
!                                        IF ( ub_ref ) THEN 

                                            wxxt(1) = duvwort(1) ** 2
                                            wxxt(2) = duvwort(2) ** 2
                                            wxxt(3) = duvwort(3) ** 2
                                            wxxt(4) = 2.0_wp * duvwort(1) * duvwort(2)
                                            wxxt(5) = 2.0_wp * duvwort(1) * duvwort(3)
                                            wxxt(6) = 2.0_wp * duvwort(2) * duvwort(3)

!                                        ENDIF

                                        DO l = 1, ngauss
                                            VTS_duvwort(1:3,l) = VTS(l) * duvwort
                                            aexp_qform(l) = a(l) * EXP ( -0.5_wp * DOT_PRODUCT ( duvwort, VTS_duvwort(1:3,l) ) )
!                                            IF ( ub_ref ) VTSD_wxxt(1:6,l) = MATMUL ( VTSD(1:6,1:6,l), wxxt )
                                            VTSD_wxxt(1:6,l) = MATMUL ( VTSD(1:6,1:6,l), wxxt )
                                        ENDDO


!                                       Seems the right place to calculate this since all constants are in place:
                                        TEMP36_ALL_GAUSS = 0
                                        DO l = 1, ngauss
                                            utemp(1:6,1) = det_scaled_deriv(1:6,l) + VTSD_wxxt(1:6,l) 
                                            vtemp(1:3,1) = VTS_duvwort(1:3,l)

!                                           Initialise:
                                            TEMP36 = 0.5_wp * MATMUL ( vtemp, TRANSPOSE ( utemp ) )

!                                           Accumulate:
                                            DO ij = 1, 6

!                                               Convert to squared matrix:
                                                VTSD_ROW = .CONVERT. VTSD(ij,1:6,l)
                                                TEMP36(1:3,ij) = TEMP36(1:3,ij) - MATMUL ( VTSD_ROW, duvwort )

                                            ENDDO
                                            TEMP36_ALL_GAUSS = TEMP36_ALL_GAUSS + TEMP36 * aexp_qform(l)
                                        ENDDO

                                        IF ( ub_ref ) THEN
                                            UUMAT = 0

                                            DO l = 1, ngauss
                                                DSECOND = 0.0_wp
                                                DO ts = 1, 6
!                                                   THis is 256 multiplications * ngauss
                                                    DSECOND = DSECOND +  CTS(1:6,1:6,ts,l) * wxxt(ts)
                                                ENDDO
                                                utemp(1:6,1) = VTSD_wxxt(1:6,l)

!                                               Fourth order coefs:
                                                DFOURTH = 0.5_wp * MATMUL ( utemp, TRANSPOSE ( utemp ) )

                                                UUMAT = UUMAT + (DZERO(1:6,1:6,l) + DSECOND + DFOURTH) * aexp_qform(l)
                                            ENDDO

                                            UUMAT = ( 0.5_wp * map_value * occ ) * UUMAT ! correct common scale aniso

                                        ENDIF

                                        UO_ALL_GAUSS = 0.0_wp
                                        DO l = 1, ngauss
                                            utemp(1:6,1) = det_scaled_deriv(1:6,l) + VTSD_wxxt(1:6,l) 
                                            UO_ALL_GAUSS = UO_ALL_GAUSS - utemp(1:6,1) * aexp_qform(l)
                                        ENDDO

                                        UO_ALL_GAUSS = 0.5_wp * map_value * UO_ALL_GAUSS


                                        any_xyz1_aniso:IF ( ANY ( refine_atomic_xyz(1:3,1 ) ) ) THEN

                                            XX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
!                                               ---------------------------------------------------------
!                                               XX AREA
                                                DO l = 1, ngauss
                                                    vtemp(1:3,1) = VTS_duvwort(1:3,l)
                                                    VXXT = MATMUL ( vtemp, TRANSPOSE ( vtemp ) )

!                                                   Note change of sign (tupi * i * h) ** 2 results in minus sign...
                                                    BMAT = BMAT - occ * map_value * (VXXT - VTS(l)) * aexp_qform(l) 
                                                ENDDO

                                            ENDIF XX_area

!                                           Check whether Biso refinable for JAT atom:
                                            XB_area:IF ( refine_atomic_biso(2) ) THEN

!                                               Accumulate matrix values wrt XB+YB+ZB:
                                                XB = TEMP36_ALL_GAUSS

!                                               Average along rows (just first 3 elements for each row):
                                                DO k = 1, 3
                                                    gg(k) = SUM ( XB(k,1:3) ) / ( 8.0_wp * pi_squared )
                                                ENDDO

!                                               Switch sign for H2 terms:
                                                IF ( H2 ) gg = -gg
                                                
!                                               Select:
                                                DO k = 1, 3
                                                    IF ( refine_atomic_xyz(k,1) ) THEN
                                                        elem = table(k,4)
                                                        IF ( elem > 0 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + occ * map_value * gg(k)
                                                        ENDIF
                                                    ENDIF
                                                ENDDO

                                            ENDIF XB_area 

!                                           Accumulate derivatives in XO area fro XYZ rows:
                                            XO_area:IF ( refine_atomic_occ(2) ) THEN
!                                                JAT atom seems to have refinable occupancy: BUG CORRECTED aexp_qform and
!                                                                                            VTS_duwort  used default.dims:
                                                 gg = -occ_i * map_value &
                                                    * MATMUL ( VTS_duvwort(1:3,1:ngauss), aexp_qform(1:ngauss) )

!                                                Switch sign for H2 terms:
                                                 IF ( H2 ) gg = -gg
                                                 DO k = 1, 3
                                                     IF ( refine_atomic_xyz(k,1) ) THEN
                                                         elem = table(k,5)
                                                         IF ( elem > 0 ) THEN
                                                             H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                         ENDIF
                                                     ENDIF
                                                 ENDDO
                                            ENDIF XO_area

!                                           Entering XU area:
                                                      
                                            XU_area:IF ( ANY ( refine_atomic_us(1:6,2) ) )  THEN

                                                IF ( .NOT. H2 ) THEN
                                                    XU_SYM = XU_SYM + (occ * map_value ) * TEMP36_ALL_GAUSS
                                                ELSE
                                                    XU_SYM = XU_SYM - (occ * map_value ) * TEMP36_ALL_GAUSS
                                                ENDIF

                                            ENDIF XU_area

                                        ENDIF any_xyz1_aniso

!                                       ======= Biso row =======
                                        biso1:IF ( refine_atomic_biso(1) ) THEN

                                            BX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
                                                
!                                               Average along rows (just 3 elements for each row):
                                                DO k = 1, 3
                                                    BX_SYM(k) = BX_SYM(k) - occ * map_value &
                                                              * SUM ( TEMP36_ALL_GAUSS(k,1:3) ) / (8.0_wp * pi_squared)
                                                ENDDO
!                                               This will be transformed by symmetry later...

                                            ENDIF BX_area

!                                           Entering dummy BB area:
!                                           Check whether JAT atom has refinable Biso:
                                            BB_area:IF ( refine_atomic_biso(2) ) THEN
!                                                See isotropic area...
                                                 WRITE(*,"(' iat=',I6,' jat=',I6, ' table=',L1, ' at_least_one_aniso=',L1)")&
                                                             iat, jat, table(4,4), at_least_one_aniso
                                                 WRITE(*,"(' U(iat)=', 6F10.6)") U_iat%u ! FIXME needs better printing for U
                                                 WRITE(*,"(' U(jat)=', 6F10.6)") U_jat%u
                                                 CALL die('O-o-o-p-s-s. Cannot get in here... Logic is wrong.', srname)
                                            ENDIF BB_area
                              
!                                           BO area: B isotropic, Oj atom has anizo values...

!                                           Need averaged gradients against U
!                                           Let's pretend we have UO refinement:
                                            BO_area:IF ( refine_atomic_occ(2) ) THEN

                                                elem = table(4,5)
!                                                UO = 0                                         

                                                IF ( elem > 0 ) THEN

                                                    UO = occ_i * UO_ALL_GAUSS

                                                    IF ( .NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + SUM ( UO(1:3) ) / ( 8.0_wp * pi ** 2 )
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - SUM ( UO(1:3) ) / ( 8.0_wp * pi ** 2 )
                                                    ENDIF
                                                ENDIF

                                            ENDIF BO_area


!                                           Check JAT atom for anisotropy refinement:
!                                           Severe boug discovered by Intel debugger JULY 2009 BVS:
!                                            BU_area:IF ( ANY ( refine_atomic_Us(2,1:6) ) ) THEN
                                            BU_area:IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN

!                                               Pretend we have UU refinement, then average. We need (1,1:6) vector.
!                                               Infamous UU 6x6 area:
                                               
                                                IF ( .NOT. H2 ) THEN
                                                    BU_SYM = BU_SYM + UUMAT
                                                ELSE
                                                    BU_SYM = BU_SYM - UUMAT
                                                ENDIF

                                            ENDIF BU_area

                                        ENDIF biso1

!                                       Finished with Biso row.

!                                       ===== Occupancy row start ======

                                        occ1:IF ( refine_atomic_occ(1 )  ) THEN

!                                           OX 1x3 area. Check JAT atom:
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                                Note "+" change of sign here:
                                                 OX_SYM = OX_SYM + occ_j * map_value &
                                                        * MATMUL ( VTS_duvwort(1:3,1:ngauss), aexp_qform(1:ngauss) )

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            OB_area:IF ( refine_atomic_biso(2) ) THEN

                                                elem = table(5,4)

!                                               Pretend it's OU area and average:
                                                IF ( elem > 0 ) THEN

!                                                   HEHE  We are in fast_sparse_normal_matrix routine: 
                                                    UO = occ_j *  UO_ALL_GAUSS

                                                    IF (.NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + SUM( UO(1:3) ) / (8.0_wp * pi_squared)
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - SUM( UO(1:3) ) / (8.0_wp * pi_squared)
                                                    ENDIF    

                                                ENDIF

                                            ENDIF OB_area

!                                           Entering OO area:
                                            OO_area:IF ( refine_atomic_occ(2) ) THEN

                                                elem = table(5,5)
                                                IF ( elem > 0 ) THEN
                                                    IF (.NOT. H2 ) THEN
                                                        H%val(nnz+elem) = H%val(nnz+elem) + map_value * SUM ( aexp_qform(1:ngauss) )
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - map_value * SUM ( aexp_qform(1:ngauss) )
                                                    ENDIF 
                                                ENDIF
                                                
                                            ENDIF OO_area

!                                           Entering OU area:
                                            OU_area: iF ( any ( refine_atomic_us(1:6,2) ) )THEN

                                                IF ( .NOT. H2 ) THEN
                                                    OU_SYM = OU_SYM + occ_j * UO_ALL_GAUSS
                                                ELSE
                                                    OU_SYM = OU_SYM - occ_j * UO_ALL_GAUSS
                                                ENDIF

                                            ENDIF OU_area

                                        ENDIF occ1


!                                       ======= U row ======

                                        any_us1:IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN 

                                            UX_area:IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN 

!                                               Accumulate for all grid points for a given symmetry operator:
                                                UX_SYM = UX_SYM - (occ * map_value) * TEMP36_ALL_GAUSS

                                            ENDIF UX_area

!                                           Entering UB 6x1 area:
                                            UB_area:IF ( refine_atomic_biso(2) ) THEN

!                                               Pretend we have UU area, then average:

!                                               DO some averaging:
                                                DO k = 1, 6

!                                                   6:11,4 area:
                                                    elem = table(5+k,4)

                                                    IF ( elem > 0 ) THEN

                                                        IF ( .NOT. H2 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + SUM ( UUMAT(k,1:3) ) &
                                                                            / (8.0_wp * pi_squared)
                                                        ELSE
                                                            H%val(nnz+elem) = H%val(nnz+elem) - SUM ( UUMAT(k,1:3) ) &
                                                                            / (8.0_wp * pi_squared)
                                                        ENDIF

                                                    ENDIF

                                                ENDDO


                                            ENDIF UB_area

!                                           Entering 6x1 UO area:
                                            UO_area:IF ( refine_atomic_occ(2) ) THEN

!                                               It's important to have occ_i here, occ_j vanishes:
                                                UO = occ_i * UO_ALL_GAUSS

                                                DO k = 1,6

                                                    elem = table(5+k, 5)

                                                    IF ( elem > 0 ) THEN

                                                        IF ( .NOT. H2 ) THEN
                                                            H%val(nnz+elem) = H%val(nnz+elem) + UO(k)
                                                        ELSE
                                                            H%val(nnz+elem) = H%val(nnz+elem) - UO(k)     
                                                        ENDIF

                                                     ENDIF
                                                ENDDO

                                                
                                            ENDIF UO_area

!                                           Entering infamous UU area:
                                            UU_area: IF ( ANY ( refine_atomic_Us(1:6,2)  ) ) THEN

                                                IF ( .NOT. H2 ) THEN
                                                    UU_SYM = UU_SYM + UUMAT
                                                ELSE
                                                    UU_SYM = UU_SYM - UUMAT
                                                ENDIF

                                            ENDIF UU_area

                                        ENDIF any_us1

!                                       ====== aniso business is over ======
                                    ELSE aniso_business

!                                       ======= Isotropic business =======

!                                       Arrays to speed up:
                                        ber2(1:ngauss)   = be(1:ngauss) * r2
                                        exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                        abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)
                                        abexp_be(1:ngauss) = abexp(1:ngauss) * be(1:ngauss)

!                                       Some constants for speed-up:
                                        IF ( ANY ( refine_atomic_biso ) ) THEN
                                                sum_bbb_52 = SUM ( abexp_be(1:ngauss) * bbb(1:ngauss) &
                                                           * (5.0_wp / 2.0_wp - ber2(1:ngauss)) )

                                            IF ( ( refine_atomic_biso(1) .AND. refine_atomic_occ(2) ) .OR. &
                                                ( refine_atomic_biso(2) .AND. refine_atomic_occ(1) ) ) THEN
                                                 sum_bbb_15 = SUM ( abexp(1:ngauss) * bbb(1:ngauss) * (-1.5_wp + ber2(1:ngauss)) )
                                            ENDIF
                                        ENDIF

!                                       ======= Isotropic business ======= BUG CORRECTED JULY 2009 ARRAYS MOVED ABOVE.

!                                       Additional array for speed up:
                                        IF ( ANY ( refine_atomic_xyz ) ) THEN

!                                           This constant is only associated with xyz refinement:
                                            sum_abexp_be = SUM ( abexp_be(1:ngauss) )

!                                           Convert to orthogonal system of coodinates (easier to calculate derivatives):
                                            duvwort = 2.0_wp * map_1(my_first_allocated_map)%ORT * duvw
 
                                            anyxyz:IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

!                                               Makes sense to calculate if JAT atom has refinable XYZ:
                                                IF (  ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                                                   Calculate Tronrud XXT matrix:
                                                    DO l = 1, 3
                                                        DO k = 1, 3
!                                                           DUVWORT has been scaled by a factor of 2 already:
                                                            XXT(l,k) = duvwort(k) * duvwort(l)
                                                        ENDDO
                                                    ENDDO

!                                                   Build 3x3 matrix for orthogonal coordinates:
                                                    BMAT = BMAT + map_value * occ * (sum_abexp_be * GMAT           &
                                                                - SUM ( abexp_be(1:ngauss) * be(1:ngauss) ) * XXT)

                                                ENDIF

!                                               Check whether Biso refinable for JAT atom:
                                                IF ( refine_atomic_biso(2) ) THEN

!                                                   Accumulate matrix values wrt XB+YB+ZB(14):
                                                    gg = -(occ * map_value * sum_bbb_52) * duvwort

!                                                   Switch sign for H2 terms:
                                                    IF ( H2 ) gg = -gg

!                                              General case -> not every XYZ(iat) is refinable:                                            
                                                    IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

                                                        DO k = 1, 3

                                                            IF ( refine_atomic_xyz(k,1) ) THEN

!                                                               HMAT(indxyz+l,jndbiso):
                                                                elem = table(k,4)

                                                                IF ( elem > 0 ) THEN
                                                                    H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                                ENDIF        

                                                            ENDIF

                                                        ENDDO

                                                    ENDIF

                                                ENDIF

!                                               Accumulate derivatives XO, YO, ZO(15) for XYZ rows:
                                                IF ( refine_atomic_occ(2) ) THEN

!                                                   JAT atom seems to have refinable occupancy:
                                                    gg = occ_i * map_value * sum_abexp_be * duvwort

!                                                   Switch sign for H2 terms:
                                                    IF ( H2 ) gg = -gg

!                                                   Considering general case -> not every XYZ(iat) is refinable:
                                                    IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN                                                

                                                        DO k = 1, 3
                                                            IF ( refine_atomic_xyz(k,1) ) THEN

!                                                               HMAT(indxyz+l,jndocc) = HMAT(indxyz+l,jndocc) + gg(k)
                                                                elem = table(k,5)

                                                                IF ( elem > 0 ) THEN
                                                                    H%val(nnz+elem) = H%val(nnz+elem) + gg(k)
                                                                ENDIF

                                                            ENDIF
                                                        ENDDO

                                                    ENDIF

                                                ENDIF

                                            ENDIF anyxyz
                                        ENDIF


!                                   -------------- Biso row ---------------------------

                                        iso1:IF ( refine_atomic_biso(1) ) THEN

!                                           FIXME:
!                                           These derivatives seem to be in lower triangular in some circumstances:

!                                           Accumulate derivatives for BX, BY, BZ(14):
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) .AND. ANY ( table(4,1:3) > 0 ) ) THEN
                                                BX_SYM = BX_SYM + (occ * map_value * sum_bbb_52) * duvwort
                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Calculate pointer to position in sparse matrix:
                                                elem = table(4,4)
                                                IF ( elem <= 0 ) STOP 'STOP ERROR'
                                                
!                                               Accumulate BB matrix elements (12):
                                                IF ( .NOT. H2 ) THEN

!                                                   HMAT(indbiso,jndbiso):
                                                    H%val(nnz+elem) = H%val(nnz+elem) + occ * map_value           &
                                                                    * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)     &
                                                                    * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp) &
                                                                    + 15.0_wp / 4.0_wp ) )

                                                ELSE

                                                    H%val(nnz+elem) = H%val(nnz+elem) - occ * map_value           &
                                                                    * SUM  ( abexp(1:ngauss) * bbb2(1:ngauss)     &
                                                                    * (ber2(1:ngauss) * (ber2(1:ngauss) - 5.0_wp) &
                                                                    + 15.0_wp / 4.0_wp ) )
                                               
                                                ENDIF

                                            ENDIF                                                    

!                                           Check whether JAT atom has refinable occ:
                                            IF ( refine_atomic_occ(2) ) THEN

!                                              Calculate pointer:
                                               elem = table(4,5)

!                                              Check that this is upper triangular part of the matrix:
                                               IF ( elem > 0 ) THEN

!                                                   Accumulate BQ matrix elements (16):
                                                    IF ( .NOT. H2 ) THEN
!                                                       HMAT(indbiso,jndocc):
                                                        H%val(nnz+elem) = H%val(nnz+elem) + occ_i * map_value * sum_bbb_15 
                                                    ELSE
                                                        H%val(nnz+elem) = H%val(nnz+elem) - occ_i * map_value * sum_bbb_15
                                                    ENDIF

                                                ENDIF

                                            ENDIF

                                        ENDIF iso1  
                           
!                                       ------------- Now occupancy row itself, note OCC_J is used --------------:
                                        occncy:IF ( refine_atomic_occ(1) ) THEN

!                                           FIXME these derivatives maybe in lower triangular part in some circumstances:
!                                           The overhead calculation is just for atoms on the diagonal though:

!                                           Check whether JAT atom has refinable XYZ:
                                            IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
!                                                WRITE(*,*) ' Str.pl.at_least_one_aniso=', at_least_one_aniso,' iat=',iat,' jat=',jat
!                                               OX+OY+OZ interaction (15):
                                                OX_SYM = OX_SYM - (occ_j * map_value * sum_abexp_be) * duvwort

                                            ENDIF

!                                           Check whether JAT atom has refinable Biso:
                                            IF ( refine_atomic_biso(2) ) THEN

!                                               Calculate pointer:
                                                elem = table(5,4)

!                                               This will implicitly indicate that iat /= jat:
                                                IF ( elem > 0 ) THEN

                                                
!                                                   QB interaction(16), note plus sign here:
                                                    IF ( .NOT. H2 ) THEN
!                                                       HMAT(indocc, jndbiso):
                                                        H%val(nnz+elem) = H%val(nnz+elem) + occ_j * map_value * sum_bbb_15
                                                    ELSE
!                                                       Switch sign for H2 terms:
                                                        H%val(nnz+elem) = H%val(nnz+elem) - occ_j * map_value * sum_bbb_15
                                                    ENDIF

                                                ENDIF

                                            ENDIF

!                                           Check whether JAT atom has refinable occncy and apply(13):
                                            IF ( refine_atomic_occ(2) ) THEN
!                                               Calculate pointer:
                                                elem = table(5,5)
                                                IF ( elem <= 0 ) STOP 'STOP ERROR OCC'
!                                               Switch sign for H2 terms:
                                                IF ( .NOT. H2 ) THEN
!                                                   HMAT(indocc,jndocc):
                                                    H%val(nnz+elem) = H%val(nnz+elem) + map_value * SUM ( abexp(1:ngauss) )
                                                ELSE
                                                    H%val(nnz+elem) = H%val(nnz+elem) - map_value * SUM ( abexp(1:ngauss) )
                                                ENDIF

                                            ENDIF

                                        

                                        ENDIF occncy


                                    ENDIF aniso_business


                                ENDIF rcut2

                            ENDDO
                        ENDDO
                    ENDDO boxzz
                ENDDO Hterms
                IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) )  THEN

!                   Build 3x3Nsym matrix:
                    AMAT(1:3,3*(m-1)+1:3*(m-1)+3) = BMAT
                ENDIF

!               BUG CORRECTED JAN 2008 BVS (CANNOT MERGE WITH PREVIOUS IF):
                IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!                   Take care of chain rule for BXYZ_DERIV & OXYZ_DERIV:
                    IF ( refine_atomic_biso(1) ) THEN
                        BX_SYM = MATMUL ( SMAT, BX_SYM )
                    ENDIF

                    IF ( refine_atomic_occ(1) ) THEN
                        OX_SYM = MATMUL ( SMAT, OX_SYM )
                    ENDIF                      

!                   Finish with chain rule for two classes of second_derivatives:
                    IF ( refine_atomic_biso(1) ) THEN

!                       Consider general case -> not every XYZ(jat) is refinable:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,2) ) THEN

!                               Calculate pointer:
                                elem = table(4,k)

!                               Upper triangular only:
                                IF ( elem > 0 ) THEN

!                                   HMAT(indbiso,jndxyz+l):
                                    H%val(nnz+elem) = H%val(nnz+elem) + BX_SYM(k)
                                ENDIF

                            ENDIF

                        ENDDO

                    ENDIF

                    IF ( refine_atomic_occ(1) ) THEN

!                       Consider general case -> not every XYZ(jat) is refinable:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,2) ) THEN

!                               Calculate pointer:
                                elem = table(5,k)

!                               Upper triangular:
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = H%val(nnz+elem) + OX_SYM(k)
                                ENDIF

                            ENDIF

                        ENDDO
                    ENDIF
                ENDIF

!               ==========
!               Aniso area:
!               ==========
                XU_sym_area:IF ( ANY ( refine_atomic_xyz(1:3,1) ) ) THEN

                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN
            
!                       XU_SYM contains accumulated XU matrix from all grid points in atom sphere:
                        CALL aniso_symop_matmul(XU_SYM, XU_SYM, RU(m))

!                       Select:
                        DO k = 1, 3
                            IF ( refine_atomic_xyz(k,1) ) THEN
                                DO ij = 1, 6
                                    elem = table(k,5+ij)
                                    IF ( elem > 0 ) THEN                                                            
                                        H%val(nnz+elem) = H%val(nnz+elem) + XU_SYM(k,ij)
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO

                    ENDIF

                ENDIF XU_sym_area

                Ux_sym_area:IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN
                    IF ( ANY ( refine_atomic_xyz(1:3,2) ) ) THEN
                        UX_SYM = MATMUL ( SMAT, UX_SYM )
                   
        
                        DO ij = 1, 6
                            IF ( refine_atomic_Us(ij,1) ) THEN
                                DO k = 1, 3
                                    IF ( refine_atomic_xyz(k,2) ) THEN
                                        elem = table(5+ij,k)
                                        IF ( elem > 0 ) THEN
                                            H%val(nnz+elem) = H%val(nnz+elem) + UX_SYM(k,ij)
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF Ux_sym_area


                UU_sym_area: IF ( ANY ( refine_atomic_Us(1:6,1) ) ) THEN

                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN
                        CALL aniso_symop_matmul(UU_SYM, UU_SYM, RU(m))
                        DO k = 1, 6
                            DO ij = 1, 6
                                elem = table(k+5,ij+5)
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = H%val(nnz+elem) + UU_SYM(k,ij)
                                ENDIF
                            ENDDO
                        ENDDO
                    ENDIF

                ENDIF UU_sym_area

                OU_sym_area: IF ( refine_atomic_occ(1) ) THEN
                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN

                        OU_SYM = OU_SYM * RU(m)
!                        OU_SYM = MATMUL ( OU_SYM, RU(m)%a )

                        DO k = 1,6
                            elem = table(5,5+k)
                            IF ( elem > 0 ) THEN
                                H%val(nnz+elem) = H%val(nnz+elem) + OU_SYM(k)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF OU_sym_area


                BU_sym_area: IF ( refine_atomic_biso(1) ) THEN
                    IF ( ANY ( refine_atomic_Us(1:6,2) ) ) THEN

                        CALL aniso_symop_matmul(BU_SYM, BU_SYM, RU(m))

!                       Do some averaging:
                        DO k = 1, 6
                            elem = table(4,5+k)
                            IF ( elem > 0 ) THEN

!                               HEHE2 changed order in SUM BVS JUL 2009:
                                H%val(nnz+elem) = H%val(nnz+elem) + SUM ( BU_SYM(1:3,k) ) &
                                                                  / (8.0_wp * pi_squared)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF BU_sym_area

!               A small caveat here -> my_first_allocated_map will have value of "2" in outer loops.

            ENDDO symmetry

!           Since this was private variable in parallel loop we need reasonable value for this:
            my_first_allocated_map = 1

            chainrule:IF ( ANY ( refine_atomic_xyz(1:3,1) ) .AND. ANY ( refine_atomic_xyz(1:3,2) ) ) THEN

!               Finish with XYZ chain rule for this atom, ASYM is on the right side:
                BMAT = MATMUL ( AMAT, ASYM ) 

!               Extract and add refinable items only:
                DO k = 1, 3
                    IF ( refine_atomic_xyz(k,1) ) THEN
                        DO kk = 1, 3
                            IF ( refine_atomic_xyz(kk,2) ) THEN
!                               Calculate pointer to HMAT(indxyz+l,jndxyz+ll):
                                elem = table(k,kk)
                                IF ( elem > 0 ) THEN
                                    H%val(nnz+elem) = BMAT(k,kk)
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO

!               Debug if necessary:
                debugg:IF ( debug > 25 ) THEN

                    IF ( nsym <= 24 ) THEN

                        DO k = 1, 3
                            WRITE(*,"(' AMAT=',72ES9.2)") &
                            (AMAT(k,l), l = 1, 3 * nsym)
                        ENDDO

                    ENDIF
   
                    DO k = 1, 3 * nsym
                        WRITE(*,"(' ASYM=', 3F9.5)") (ASYM(k,l), l = 1, 3)
                    ENDDO

                        DO k = 1, 3
                            WRITE(*,"(' BMAT=', 3ES9.2)") (BMAT(k,l), l = 1, 3)
                        ENDDO

                ENDIF debugg

            ENDIF chainrule


!           Monitor CPU:
            monitor:DO k = 1, 11 
                DO kk = 1, 11

!                   Calculate pointer:
                    elem = table(k,kk)

                    refinables:IF ( elem > 0 ) THEN

!$OMP CRITICAL
                        pair = pair + 1
                        IF ( MOD ( nnz+elem, modh ) == 0 ) THEN
                            CALL CPU_TIME ( time1 )
                            WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ',                                               &
                           &'matrix element[', A, ',', A, ']=', T62, ES12.5, 2X,'CPU time=',F8.1, ' s', I10, ' thread=',I3)") &
                             TRIM ( int_to_c ( H%ia(nnz+elem ) ) ), TRIM ( int_to_c ( H%ja(nnz+elem ) ) ),     &
                             H%val(nnz+elem), time1 - time0, pair, OMP_GET_THREAD_NUM()
                        ENDIF

                        IF ( H%ia(nnz+elem) == 0 .OR. H%ja(nnz+elem) == 0 ) THEN
                            WRITE(*,*) ' k=  ', k,   ' kk=', kk
                            WRITE(*,*) ' nnz=', nnz, ' elem=',elem
                            WRITE(*,*) ' H%ia(nnz+elem)=',  H%ia(nnz+elem)
                            WRITE(*,*) ' H%ja(nnz+elem)=',  H%ja(nnz+elem)
                            CALL die('O-o-o-o-p-s-s... Programming error. Incorrect coordinate matrix indices.',&
                                      srname)
                        ENDIF

!                       Check diagonal matrix elements for non-positive values:
                        IF ( H%ia(nnz+elem) ==  H%ja(nnz+elem) .AND. H%val(nnz+elem) <= 0.0_wp ) THEN
                            WRITE(*,*) ' k=  ', k,   ' kk=', kk
                            WRITE(*,*) ' nnz=', nnz, ' elem=',elem
                            WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ', 'matrix element=', 2I8, ' HMAT=', ES12.5, &
                           &' iat/jat=', 2I6, ' ref. params=', 2I3, ' nnz+elem=',I6)") &
                            H%ia(nnz+elem), H%ja(nnz+elem),  H%val(nnz+elem) ,iat, jat, k, kk, nnz+elem

                            WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ', 'Biso=', ES9.2, ' occ=', ES9.2, ' r_cut_off=', F8.2)") &
                            biso_iat + biso_jat, occ, r_cut_off

!                           Bailing out:
                            CALL die('Programming error. Zero or negative normal matrix diagonal element. &
                                     &Refinement halted. Sorry.', srname)

                        ENDIF

!                        Check diagonal matrix elements for abnormal values:
                        IF ( isNaN ( H%val(nnz+elem) ) ) THEN
                            WRITE(*,*) ' k=  ', k,   ' kk=', kk
                            WRITE(*,*) ' nnz=', nnz, ' elem=',elem
                            WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ', 'matrix element=', 2I8, ' HMAT=', ES12.5, &
                           &' iat/jat=', 2I6, ' ref. params=', 2I3, ' nnz+elem=',I6)") &
                            H%ia(nnz+elem), H%ja(nnz+elem),  H%val(nnz+elem) ,iat, jat, k, kk, nnz+elem

                            WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ', 'Biso=', ES9.2, ' occ=', ES9.2, ' r_cut_off=', F8.2)") &
                            biso_iat + biso_jat, occ, r_cut_off

!                           Bailing out:
                            CALL die('Programming error. Abnormal matrix element. &
                                     &Refinement halted. Sorry.', srname)

                        ENDIF
!$OMP END CRITICAL

                    ENDIF refinables

                ENDDO

            ENDDO monitor

!           This is not suitable for parallel computing:
!           nnz = nnz + block_size

! Seems like a SUPERBUG APR 2010 (please don't touch this):
!$OMP CRITICAL
            IF ( ALLOCATED ( AMAT ) ) CALL deallocate_array(AMAT)
!$OMP END CRITICAL

        ENDDO coord
!$OMP END PARALLEL DO 


!   FIXME should be moved to the calling routine:     
        IF ( output_stats ) THEN


!            Report integration radii stats:
             WRITE(*,"(' FAST_SPARSE_NORMAL_MATRIX> ', 'R_CUT_OFF: min=', F8.1, ' max=', F8.1, ' max possible=', F8.1)")&
                         min_r_cut_off, max_r_cut_off, SQRT ( r2_cell_max )
        ENDIF

!       Since nnz is not defined due to parallelization we must keep it positive:
        nnz = 1

!       Free memory:
        IF ( ALLOCATED ( ASYM ) ) CALL deallocate_array(ASYM)
        IF ( ALLOCATED ( RU ) )   CALL deallocate_array(RU)

    END SUBROUTINE fast_sparse_normal_matrix

    SUBROUTINE my_normal_matrix_diagonal(g, mtz_1, map_1, pdb_2, h2)
!
!       Purpose:
!       =======
!
!       Calculates diagonal elements G of normal matrix (H1 + H2 terms) with respect to
!       X,Y,Z,B and Q params.
!       Need to calc this routine twice since it requires a lot of memory.
!       In the first run only H1 terms will be calculated.
!       In the second run H2 terms will be added to H1 terms.

!       Supply 5 maps and you are done. Later we can increase this to 10 maps and calculate in one run.
!
!       Date:                 Programmer:             History of changes:
!       ====                  ==========              ==================
!       October  31, 2007     B.Strokopytov           Original code
!       November 14, 2007     B.Strokopytov           Accurate estimation of integration radius added.
!       November 27, 2007     B.Strokopytov           Occupancy only refinement speed up.
!       December 17, 2007     B.Strokopytov           Presently this algorithm will work correctly  in
!                                                     orthogonal space groups only or in P1
!                                                     Need to add cross terms HK, HL, KL in coefs
!                                                     But this will require 8(!!) maps and two runs
!                                                     for H1 and H2. This is obviously an overkill.
!
!       Note:
!       ====
!       This routine in experimental state now. Please AVOID using it for general (non-orthogonal) space groups. 
!
        REAL(KIND=wp), DIMENSION(:),              INTENT(INOUT) :: g
        TYPE(mtz),                                INTENT(IN)    :: mtz_1
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: map_1
        TYPE(pdb),                                INTENT(IN)    :: pdb_2
        LOGICAL,                                  INTENT(IN)    :: h2
!       Local variables:
        TYPE(vector)                                            :: centre_of_atom_in_fractions
        TYPE(vector_int)                                        :: closest_grid_point
        TYPE(vector_int)                                        :: box_size
        INTEGER,       DIMENSION(3)                             :: jut
        INTEGER,       DIMENSION(3)                             :: box_lower_limit
        INTEGER,       DIMENSION(3)                             :: box_upper_limit
        TYPE(vector)                                            :: deort_diag
        TYPE(vector)                                            :: ort_diag
!       Normal atom params:
        INTEGER                                                 :: ityp
        INTEGER                                                 :: jtyp
        INTEGER                                                 :: ngauss
        REAL(KIND=wp)                                           :: b_scale
        REAL(KIND=wp), DIMENSION(15)                            :: a
        REAL(KIND=wp), DIMENSION(15)                            :: b
!       Joint atom:
        REAL(KIND=wp), DIMENSION(15)                            :: ajo 
        REAL(KIND=wp), DIMENSION(15)                            :: bjo
        REAL(KIND=wp), DIMENSION(15)                            :: be
        REAL(KIND=wp), DIMENSION(15)                            :: bbb 
        REAL(KIND=wp), DIMENSION(15)                            :: ber2
        REAL(KIND=wp), DIMENSION(15)                            :: exp_be
        REAL(KIND=wp), DIMENSION(15)                            :: abexp
        REAL(KIND=wp)                                           :: gaussian_scale
!       Distance:
        REAL(KIND=wp)                                           :: r2
        TYPE(vector_int)                                        :: uvw
        TYPE(vector)                                            :: duvw
        REAL(KIND=wp), DIMENSION(3)                             :: gg 
!       Cut-off:
        REAL(KIND=wp)                                           :: r2_cut_off
        REAL(KIND=wp)                                           :: r_cut_off
        REAL(KIND=wp)                                           :: max_r_cut_off
        REAL(KIND=wp)                                           :: min_r_cut_off
!       Maximum allowed sphere radius:
        REAL(KIND=wp)                                           :: r2_cell_max
!       Can handle 11 maps only:
        REAL(KIND=wp), DIMENSION(11)                            :: map_value
        REAL(KIND=wp)                                           :: occ
        REAL(KIND=wp)                                           :: occ_i
        REAL(KIND=wp)                                           :: occ_j
!       Maps variables:
        INTEGER                                                 :: number_of_maps
        INTEGER                                                 :: my_first_allocated_map
!       Refinement:        
        INTEGER                                                 :: np
        INTEGER                                                 :: natom
        LOGICAL,       DIMENSION(3)                             :: refine_atomic_xyz
        LOGICAL                                                 :: refine_atomic_biso
        LOGICAL                                                 :: refine_atomic_occ
!       Pointers:
        INTEGER                                                 :: indxyz
        INTEGER                                                 :: indbiso
        INTEGER                                                 :: indocc
!       Counters:
        INTEGER                                                 :: iat
        INTEGER                                                 :: jat
        INTEGER                                                 :: k 
        INTEGER                                                 :: l
!       Symmetry:
        TYPE(vector)                                            :: xyz_iat
        TYPE(vector)                                            :: xyz_jat
        TYPE(vector), DIMENSION(2)                              :: xyz_frac
!       Box counters:
        INTEGER                                                 :: i
        INTEGER                                                 :: j
!       Box counters:
        INTEGER                                                 :: ju
        INTEGER                                                 :: jv
        INTEGER                                                 :: jw
!       CPU monitoring:
        REAL(KIND=sp)                                           :: time0
        REAL(KIND=sp)                                           :: time1
!       Printing:
        INTEGER(KIND=eb)                                        :: modh
        INTEGER(KIND=eb)                                        :: pair
        CHARACTER(LEN=32),                               SAVE   :: srname = 'normal_matrix_diagonal'
!       Test:
        REAL(KIND=wp), DIMENSION(3)                             :: xyz

!       Checkz: 
        IF ( ALLOCATED ( pdb_2 ) ) THEN
            natom = SIZE ( pdb_2 )
        ELSE
            CALL die('Programming error. pdb_2 has not been initialized.', &
                     srname)
        ENDIF

!       This must be set from the start:
        my_first_allocated_map = first_allocated_map ( map_1 )

!       Check map allocation:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Programming error. Map array has not been allocated properly.', &
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( map_1(my_first_allocated_map) ) ) THEN
            CALL die('Programming error. MAP_1(1) has not been allocated properly.',&
                     srname)
        ELSE
            number_of_maps = SIZE ( map_1 )
            IF ( number_of_maps < 1 ) THEN
                WRITE(*,*) number_of_maps
                CALL die('Programming error. Number of maps is too small.', srname)
            ENDIF

        ENDIF

        IF ( natom < 1 ) THEN
            WRITE ( *, "(' NORMAL_MATRIX_DIAGONAL> ', 'natom=', I6)") natom
            CALL die ( 'Programming error. Need at least one atom for density generation.', srname)
        ENDIF

!       Long and hard check: 
        np = COUNT ( pdb_2%refinable )
        IF ( np /=  SIZE ( g ) ) THEN
            WRITE(*,*) np, SIZE ( g )
            CALL die('Programming error. Incorrect size of diag G array.', srname)
        ENDIF
       
        CALL messag('Going to calculate '//TRIM ( int_to_c ( SIZE ( g ) ) )//' matrix elements... Wait...',&
                    srname)

!       Initialize diag array:
        IF ( .NOT. H2 ) THEN
!           Need to accumulate G array in case of H2:
            g = 0.0_wp
        ENDIF

        b_scale    = map_1(my_first_allocated_map)%b_scale
        deort_diag = extract_diagonal ( map_1(my_first_allocated_map)%DEORT )

!       Calculate maximum sphere radius around any atom:
        ort_diag    = extract_diagonal ( map_1(my_first_allocated_map)%ORT )
        r2_cell_max = 0.5_wp ** 2 * MINVAL ( ort_diag .X. ort_diag )

!       Set for statistics:
        max_r_cut_off = -10.0_wp ** 4 
        min_r_cut_off =  10.0_wp ** 4
        
!       Initialise pointers:
        indxyz = 0

!       Timing:
        CALL CPU_TIME(time0)
        modh = MAX ( SIZE ( g ) / 20, 1 )

!       Counter for printing:
        pair  = 0 

!       Loop through all atoms:
        coord:DO iat = 1, natom

!           Speed for occupancy only refinement:
            IF ( COUNT ( pdb_2%refinable(1:11, iat ) ) == 0 ) THEN
!               Nothing to refine:
                CYCLE coord
            ENDIF

            jat  = iat

            occ_i = pdb_2%occ(iat)
            occ_j = pdb_2%occ(jat)

!           Save type of atom:
            ityp = pdb_2%atom_type(iat)
            ngauss = pdb_2%atomsf_table%record(ityp)%ngauss
            jtyp = pdb_2%atom_type(jat)

            IF ( ngauss /=  pdb_2%atomsf_table%record(jtyp)%ngauss ) THEN
                WRITE(*,*) iat, jat
                WRITE(*,*) ngauss,  pdb_2%atomsf_table%record(jtyp)%ngauss 
                CALL die('Cannot deal with differing NGAUSS for atom pair.', srname)
            ENDIF

!           Check:
            IF ( ngauss /=2 .AND. ngauss /= 5 ) THEN
                WRITE(*,*) 'ngauss=', ngauss, ' atom=', pdb_2%atom_name(iat)
                CALL die('Programming error. Abnormal number of gaussians detected...',&
                          srname)
            ENDIF

!           Copy gaussian constants:
            a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
            b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss) + ( pdb_2%biso(iat) + b_scale )

            IF ( debug > 22 ) THEN
                WRITE(*,"(' a(1:ngauss)=', 15F7.2)") a(1:ngauss)
                WRITE(*,"(' b(1:ngauss)=', 15F7.2)") b(1:ngauss)
                WRITE(*,"(' pdb_2%biso(iat)=',F7.2)") pdb_2%biso(iat)
            ENDIF
  
!           Joint atom constants:
            occ = occ_i * occ_j

!           This will work only if two atoms have equal number of gaussians: FIXME:
            gauss_outer:DO i = 1, ngauss
                gauss_inner:DO j = 1, i

                   IF ( i /= j ) THEN
                       gaussian_scale = 2.0_wp
                   ELSE
                       gaussian_scale = 1.0_wp
                   ENDIF

!                  Convert to one dimensional table (lower triange of 2-dim table):
                   l = (i * (i - 1)) / 2 + j

!                  Joint atom constants:
                   ajo(l) = gaussian_scale * a(i) * a(j) 
                   bjo(l) = b(i) + b(j) 

                ENDDO gauss_inner
            ENDDO gauss_outer

!           Redefine ngauss:
            ngauss = ngauss * ( ngauss + 1 ) / 2

! Need to put some code in here for differing ngauss...FIXME

            IF ( debug > 22 ) THEN
                WRITE(*,"(' ajo=', 15F8.2)") ajo(1:ngauss)
                WRITE(*,"(' bjo=', 15F8.2)") bjo(1:ngauss)
            ENDIF

!           Final joint atom constants (switching  back to usual names):
            a(1:ngauss)  = ajo(1:ngauss) * ( SQRT ( pi4 / bjo(1:ngauss) ) ) ** 3
            be(1:ngauss) = 4.0 * pi_squared / bjo(1:ngauss)

            IF ( debug > 22 ) THEN
                WRITE(*,"(' a =', 15F8.4)") a(1:ngauss)
                WRITE(*,"(' be=', 15F8.4)") be(1:ngauss)
            ENDIF

!           Arrays of type 1.0/( b(igauss) + biso(iat) + b_scale ) to speed up:
            bbb(1:ngauss) = 1.0_wp / bjo(1:ngauss)

            IF ( debug > 22 ) THEN
                WRITE(*,"(' bbb=', 15F8.4)") bbb(1:ngauss)
            ENDIF 

!           Finally we have sensible approach:
            l = MAXLOC ( ABS ( a(1:ngauss) ), DIM=1 )

!           Superior floating radius, gives excellent results on RNAse for high B-values:
            r2_cut_off = integration_radius ( ajo(1:ngauss), bjo(1:ngauss) ) ** 2
            IF ( debug > 22 ) THEN
                WRITE(*,*) ' ajo(l), bjo(l)=',  ajo(l), bjo(l)
                WRITE(*,*) ' optimal diag r_cut_off=', SQRT ( r2_cut_off )
            ENDIF

!           Make sure sphere around any atom fits into the cell:
            r2_cut_off = MIN ( r2_cut_off, r2_cell_max )
            r_cut_off  = SQRT ( r2_cut_off )

!           Accumulate statics for R_CUT_OFF:
            max_r_cut_off = MAX ( max_r_cut_off, r_cut_off )
            min_r_cut_off = MIN ( min_r_cut_off, r_cut_off )         

            IF ( debug > 22 ) THEN
                WRITE(*,*) ' my_first_allocated_map=', my_first_allocated_map
                WRITE(*,*) ' diag r_cut_off=', r_cut_off
            ENDIF

!           Figure out box size around joint atom:
            box_size = NINT ( r_cut_off * ( deort_diag .X. map_1(my_first_allocated_map)%map_size ) )

!           Figure out refinable parameters:
            refine_atomic_xyz  = pdb_2%refinable(1:3,iat)
            refine_atomic_biso = pdb_2%refinable(4,iat)
            refine_atomic_occ  = pdb_2%refinable(5,iat)

!           Calculate matrix pointers:
            IF ( refine_atomic_biso ) THEN
                indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
            ENDIF

            IF ( refine_atomic_occ ) THEN
                indocc  = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
            ENDIF

!           Special debugging:
            IF ( debug > 22 ) THEN
                WRITE(*,*) iat, indxyz, refine_atomic_xyz
            ENDIF

!           Convert to fractions:
            xyz_frac = (/ mtz_1%DEORT * pdb_2%xyz(iat), mtz_1%DEORT * pdb_2%xyz(jat) /)
         
            xyz_iat  = xyz_frac(1)
            symmetry:DO j = 1, mtz_1%sp_group%number_of_symops
!            symmetry:DO j = 1, 1
                xyz_jat = mtz_1%sp_group%SYM(j) * xyz_frac(2)

!               Calculate joint atom position:
                IF ( .NOT. H2 ) THEN
                    centre_of_atom_in_fractions = xyz_iat - xyz_jat
                ELSE
                    centre_of_atom_in_fractions = xyz_iat + xyz_jat
                ENDIF

                IF ( debug > 22 ) THEN
                    xyz = centre_of_atom_in_fractions
                    IF ( .NOT. H2 ) THEN
                        WRITE(*,"('MY_NORMAL_MATRIX_DIAGONAL> ', 'H1 centre= ', 3F10.5)") xyz
                    ELSE
                        WRITE(*,"('MY_NORMAL_MATRIX_DIAGONAL> ', 'H2 centre= ', 3F10.5)") xyz
                    ENDIF
                ENDIF

                closest_grid_point = NINT ( centre_of_atom_in_fractions .X. map_1(my_first_allocated_map)%map_size )

!               Box limits:
                box_lower_limit = closest_grid_point - box_size
                box_upper_limit = closest_grid_point + box_size

!               Create a box around the joint atom and calculate its density and convolute with the map:
                boxzz:DO jw = box_lower_limit(3), box_upper_limit(3)

                    DO jv = box_lower_limit(2), box_upper_limit(2)

                        DO ju = box_lower_limit(1), box_upper_limit(1)

                            uvw = (/ ju, jv, jw /)
                            duvw = REAL ( uvw ) / map_1(my_first_allocated_map)%map_size - centre_of_atom_in_fractions

!                           High point of electron density convolution:
                            r2 = map_1(my_first_allocated_map)%REAL_TENSOR * duvw   ! corresponds to r - r(n) + r(m)

                            IF ( r2 <= r2_cut_off ) THEN

                                jut = uvw .MMOD. map_1(my_first_allocated_map)%map_size

                                IF ( refine_atomic_biso ) THEN
                                    map_value(4) = map_1(4)%array(jut(1), jut(2), jut(3))
                                ENDIF

                                IF ( refine_atomic_occ ) THEN
                                    map_value(5) = map_1(5)%array(jut(1), jut(2), jut(3))
                                ENDIF


!                               Going to speed up in case of Biso/occncy refinement:
                                ber2(1:ngauss)   = be(1:ngauss) * r2
                                exp_be(1:ngauss) = EXP ( -ber2(1:ngauss) )
                                abexp(1:ngauss)  = a(1:ngauss) * exp_be (1:ngauss)

!                               Just diagonal elements:
                                IF ( ANY ( refine_atomic_xyz ) ) THEN

                                    DO l = 1, 3
                                        IF ( refine_atomic_xyz(l) ) THEN
                                            map_value(l) = map_1(l)%array(jut(1), jut(2), jut(3))
                                        ENDIF
                                    ENDDO

                                    DO k = 1, 3
                                        IF ( refine_atomic_xyz(k) ) THEN
                                            gg(k) = occ * map_value(k) * SUM( abexp(1:ngauss) )
                                        ELSE
                                            gg(k) = 0.0_wp
                                        ENDIF
                                    ENDDO
                                                  
!                                   Apply translationless inverted symmetry operator:
                                    IF ( j > 1 ) THEN
!                                       duvwort = map_1(1)%DEORT * gg
                                        gg = map_1(my_first_allocated_map)%SYM_ORT_INVERTED_DEORT(j) * gg
                                    ENDIF              

!                                   Initialise pointer:
                                    l = 0
                                    DO k = 1, 3
                                        IF ( refine_atomic_xyz(k) ) THEN
                                            l = l + 1
                                            g(indxyz+l) = g(indxyz+l) + gg(k)
                                        ENDIF                               
                                    ENDDO
                                ENDIF
                                 
                                IF ( refine_atomic_biso ) THEN

!                                   Sign of map is different for H2, see coefs module sub matrix_diag_coefs_in_P1: 
                                    g(indbiso) = g(indbiso) + occ * map_value(4) * SUM ( abexp(1:ngauss) )
                                ENDIF   
                           
                                IF ( refine_atomic_occ ) THEN
                                    g(indocc) = g(indocc) + map_value(5) * SUM ( abexp(1:ngauss) )
                                    WRITE(*,"(' agarwal ',I4,' g=',ES9.2, ' map_value= ', ES9.2, ' gauss SUM=', ES9.2)")&
                                    iat, g(indocc), map_value(5), SUM ( abexp(1:ngauss) )

                                ENDIF

                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO boxzz

            ENDDO symmetry

!           Monitor CPU:
            l = 0
            DO k = 1, 11

                IF ( pdb_2%refinable(k,iat) ) THEN
                    l = l + 1
                    pair = pair + 1
                    IF ( MOD ( pair, modh ) == 0 ) THEN
                        CALL CPU_TIME ( time1 )
                        WRITE(*,"(' NORMAL_MATRIX_DIAGONAL> ',                                            &
                       &'matrix element[', A, ',', A, ']=', T62, ES12.5, 2X,'CPU time=',F8.1, ' s', I10)")&
                        TRIM ( int_to_c ( pair ) ), TRIM ( int_to_c ( pair ) ), g(pair), time1 - time0, pair
                    ENDIF

!                   Check matrix elements for non-positive values:
                    IF ( g(indxyz+l) <= 0.0_wp ) THEN
                        WRITE(*,"(' MY_NORMAL_MATRIX_DIAGONAL> ', 'pointer=', I8, ' diag=', ES12.5, ' iat=', I6, ' ref. param=',I3)")&
                                  indxyz+l, g(indxyz+l), iat, l
                        WRITE(*,"(' MY_NORMAL_MATRIX_DIAGONAL> ', 'Biso=', F8.2, ' occ=', F6.2, ' r_cut_off=', F8.2)")&
                                  pdb_2%biso(iat), pdb_2%occ(iat), r_cut_off

!                       Bailing out:
                        CALL die('Programming error. Zero or negative normal matrix diagonal element. Refinement halted.',&
                                 srname)

                    ENDIF

                ENDIF

            ENDDO

!           Increment pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11,iat) )

        ENDDO coord

!       Report stats:
        WRITE(*,"(' MY_NORMAL_MATRIX_DIAGONAL> ', 'R_CUT_OFF: max=', F8.1, ' min=', F8.1)")&
        max_r_cut_off, min_r_cut_off

    END SUBROUTINE my_normal_matrix_diagonal

    SUBROUTINE assign_indices (H, table, pdb_2, nnz, indxyz, jndxyz, iat, jat)
        TYPE(coo_matrix),          INTENT(INOUT) :: H
        INTEGER,   DIMENSION(:,:), INTENT(IN)    :: table
        TYPE(pdb),                 INTENT(IN)    :: pdb_2
        INTEGER(KIND=eb),          INTENT(IN)    :: nnz
        INTEGER,                   INTENT(IN)    :: indxyz
        INTEGER,                   INTENT(IN)    :: jndxyz
        INTEGER,                   INTENT(IN)    :: iat
        INTEGER,                   INTENT(IN)    :: jat
!       Local variables:
        INTEGER                                  :: elem
!       Counters:
        INTEGER                                  :: i
        INTEGER                                  :: j
        INTEGER                                  :: l
        INTEGER                                  :: ll

        IF ( SIZE ( table, DIM=1 ) /= MAX_BLOCK_SIZE ) THEN
            CALL die('Programming error. Incorrect table.', 'assign_indices')
        ENDIF

!       Calculate H%IA, H%JA indices in advance (coordinate format):
        l = 0
        DO i = 1, SIZE ( table, DIM=1 )

!           We must consider situation when some parameters may be UNREFINABLE:
            IF ( pdb_2%refinable(i,iat) ) THEN

                l  = l + 1
                ll = 0

                DO j = 1, SIZE ( table, DIM=2 )

                    IF ( pdb_2%refinable(j,jat) ) THEN
                        ll = ll + 1
                        elem = table(i,j)
                        IF ( elem > 0 ) THEN

                            H%ia(nnz+elem) = indxyz + l
                            H%ja(nnz+elem) = jndxyz + ll

                            IF ( debug > 40 ) THEN
                                 WRITE(*,*) ' nnz+elem= ', nnz+elem, H%ia(nnz+elem), H%ja(nnz+elem)
                            ENDIF
                            IF ( H%ia(nnz+elem) <=0 .OR. nnz+elem > h%nnz ) THEN
                                WRITE(*,*)  H%ia(nnz+elem), nnz+elem
                                CALL die('Unreasonable index IA.', 'assign_indices')
                            ENDIF

                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDDO

    END SUBROUTINE assign_indices


END MODULE convolution
