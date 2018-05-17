MODULE map_fft
USE basic_symmetry_operations
USE mtz_io
USE fft_util
USE refinement_util
IMPLICIT NONE
TYPE :: map
    INTEGER                                     :: nu = 0
    INTEGER                                     :: nv = 0
    INTEGER                                     :: nw = 0
    TYPE(vector_int)                            :: map_size
    REAL(KIND=wp)                               :: grid_factor = 4.0
    INTEGER                                     :: hmin
    INTEGER                                     :: hmax
    INTEGER                                     :: kmin
    INTEGER                                     :: kmax
    INTEGER                                     :: lmin
    INTEGER                                     :: lmax
    REAL(KIND=wp), DIMENSION(6)                 :: cell
    TYPE(space_group)                           :: sp_group
    TYPE(matrix)                                :: ORT
    TYPE(matrix)                                :: DEORT
    TYPE(symop),  DIMENSION(:), ALLOCATABLE     :: SYM_DEORT 
    TYPE(matrix), DIMENSION(:), ALLOCATABLE     :: SYM_ORT_INVERTED 
    TYPE(matrix), DIMENSION(:), ALLOCATABLE     :: SYM_ORT_INVERTED_DEORT 
    TYPE(tensor)                                :: real_tensor
    TYPE(tensor)                                :: recip_tensor
    REAL(KIND=wp)                               :: b_scale             ! antialias factor
    REAL(KIND=wp)                               :: badd                ! antialias factor
    REAL(KIND=wp)                               :: fft_radius          ! Navaza's radius
    REAL(KIND=wp)                               :: volume
!   fftw 3.0 plans, INTEGER*8 is used:
    INTEGER(KIND=eb)                            :: real_to_complex_in_place
    INTEGER(KIND=eb)                            :: real_to_complex_out_of_place
    INTEGER(KIND=eb)                            :: complex_to_real_in_place
    INTEGER(KIND=eb)                            :: complex_to_real_out_of_place
    INTEGER(KIND=eb)                            :: sine_to_real_in_place
!   peak picking part:
    TYPE(vector),  DIMENSION(:),    ALLOCATABLE :: peak_locations
    REAL(KIND=wp), DIMENSION(:),    ALLOCATABLE :: corr_coef
!
!   map 3D array:
    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: array
END TYPE

! Global FFTW constants:
INTEGER, PARAMETER :: FFTW_R2HC=0
INTEGER, PARAMETER :: FFTW_HC2R=1
INTEGER, PARAMETER :: FFTW_DHT=2
INTEGER, PARAMETER :: FFTW_REDFT00=3
INTEGER, PARAMETER :: FFTW_REDFT01=4
INTEGER, PARAMETER :: FFTW_REDFT10=5
INTEGER, PARAMETER :: FFTW_REDFT11=6
INTEGER, PARAMETER :: FFTW_RODFT00=7
INTEGER, PARAMETER :: FFTW_RODFT01=8
INTEGER, PARAMETER :: FFTW_RODFT10=9
INTEGER, PARAMETER :: FFTW_RODFT11=10
INTEGER, PARAMETER :: FFTW_FORWARD=-1
INTEGER, PARAMETER :: FFTW_BACKWARD=+1
INTEGER, PARAMETER :: FFTW_MEASURE=0
INTEGER, PARAMETER :: FFTW_DESTROY_INPUT=1
INTEGER, PARAMETER :: FFTW_UNALIGNED=2
INTEGER, PARAMETER :: FFTW_CONSERVE_MEMORY=4
INTEGER, PARAMETER :: FFTW_EXHAUSTIVE=8
INTEGER, PARAMETER :: FFTW_PRESERVE_INPUT=16
INTEGER, PARAMETER :: FFTW_PATIENT=32
INTEGER, PARAMETER :: FFTW_ESTIMATE=64
INTEGER, PARAMETER :: FFTW_ESTIMATE_PATIENT=128
INTEGER, PARAMETER :: FFTW_BELIEVE_PCOST=256
INTEGER, PARAMETER :: FFTW_DFT_R2HC_ICKY=512
INTEGER, PARAMETER :: FFTW_NONTHREADED_ICKY=1024
INTEGER, PARAMETER :: FFTW_NO_BUFFERING=2048
INTEGER, PARAMETER :: FFTW_NO_INDIRECT_OP=4096
INTEGER, PARAMETER :: FFTW_ALLOW_LARGE_GENERIC=8192
INTEGER, PARAMETER :: FFTW_NO_RANK_SPLITS=16384
INTEGER, PARAMETER :: FFTW_NO_VRANK_SPLITS=32768
INTEGER, PARAMETER :: FFTW_NO_VRECURSE=65536
INTEGER, PARAMETER :: FFTW_NO_SIMD=131072
!!INTEGER, PARAMETER :: NUMBER_OF_THREADS=2

INTERFACE allocated
    MODULE PROCEDURE allocated_map
END INTERFACE

INTERFACE size
    MODULE PROCEDURE map_size
END INTERFACE

INTERFACE allocate_map
    MODULE PROCEDURE allocate_real_map
END INTERFACE

INTERFACE deallocate_map
    MODULE PROCEDURE deallocate_real_map
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_of_real_maps
    MODULE PROCEDURE allocate_array_of_real_maps_2
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_array_of_real_maps
END INTERFACE

CONTAINS
    FUNCTION allocated_map ( map_1 )
        LOGICAL               :: allocated_map
        TYPE(map), INTENT(IN) :: map_1
        allocated_map = ALLOCATED(map_1%array)
    END FUNCTION allocated_map

    FUNCTION map_size ( map_1 )
        TYPE(vector_int)      :: map_size
        TYPE(map), INTENT(IN) :: map_1
        map_size = UBOUND(map_1%array) - LBOUND(map_1%array) + 1
    END FUNCTION map_size

    SUBROUTINE real_to_complex_plan_in_place ( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1
!       Local variables:
        INTEGER, DIMENSION(3)    :: temp
        INTEGER, DIMENSION(3)    :: map_size
        INTEGER, EXTERNAL        :: omp_get_num_threads

        temp = (/map_1%nu, map_1%nv, map_1%nw/)
        map_size = map_1%map_size
!
!       Check map structure:
        IF ( ANY( temp /= map_size ) ) THEN
            WRITE(*,"(' REAL_TO_COMPLEX_PLAN_IN_PLACE> ','nu, nv, nw= ', 3I5)")&
            map_1%nu, map_1%nv, map_1%nw 
            WRITE(*,"(' REAL_TO_COMPLEX_PLAN_IN_PLACE> ','map_size=   ', 3I5)")&
            map_size
            CALL die ('Internal map structure error.',&
                      'real_to_complex_plan_in_place')
        ENDIF

!       Check array bounds:
        temp  = UBOUND(map_1%array) - LBOUND(map_1%array) + 1
        WRITE(*,"(' REAL_TO_COMPLEX_PLAN_IN_PLACE> ', 'map array bounds: ', 3I5)")  temp 
        WRITE(*,"(' REAL_TO_COMPLEX_PLAN_IN_PLACE> ', 'map_size        : ', 3I5)")  map_size
!
!       We have 2 extra sections along x-axis  with double values:
        temp(1) = temp(1) - 2
        IF ( ANY ( temp /= map_size ) ) THEN
            CALL die ('Incorrect map dimensions.', 'real_to_complex_plan_in_place')
        ENDIF

!      fftw 3.0 version real to complex plan:
       CALL dfftw_plan_dft_r2c (map_1%real_to_complex_in_place, 3, map_size, &
                                map_1%array, map_1%array, FFTW_PATIENT+FFTW_UNALIGNED)
!                                map_1%array, map_1%array, FFTW_ESTIMATE)

    END SUBROUTINE real_to_complex_plan_in_place

    SUBROUTINE real_fft_in_place ( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1

        CALL dfftw_execute_dft_r2c ( map_1%real_to_complex_in_place, map_1%array, map_1%array )
    END SUBROUTINE real_fft_in_place

    SUBROUTINE destroy_real_fft_plan ( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1

        CALL dfftw_destroy_plan ( map_1%real_to_complex_in_place )
    END SUBROUTINE destroy_real_fft_plan

    SUBROUTINE complex_to_real_plan_in_place ( map_1 )
!
!       Purpose:
!       =======
!       Prepares complex-to-real plan for FFTW3
!
!       Date:          Programmer:               History of changes:
!       ====           ==========                ==================
!       Oct 2003       B.Strokopytov             Original code.
!      
!       Oct 2007       B.Strokopytov             DIMS array added (to avoid Intel compiler warnings).
!       Oct 2007       B.Strokopytov             Cosmetic changes.
!
        TYPE(map),             INTENT(INOUT) :: map_1
        INTEGER,  DIMENSION(3)               :: map_size
        INTEGER,  DIMENSION(3)               :: dims
        CHARACTER(LEN=32)                    :: srname
        INTEGER, EXTERNAL                    :: omp_get_num_threads

        srname = 'complex_to_real_plan_in_place'

!       Basic checks:
        IF ( .NOT. ALLOCATED ( map_1%sp_group%non_polar_axes ) ) THEN
            CALL die('Programming error: non_polar_axes must be initialized first.',&
                     srname)
        ELSE IF ( SIZE ( map_1%sp_group%non_polar_axes ) == 0 ) THEN
            CALL die('No plan can be created: all axes are polar.', &
                      srname)
        ENDIF

        map_size = map_1%map_size
        dims     = map_size(map_1%sp_group%non_polar_axes)
!        NUMBER_OF_THREADS = OMP_GET_NUM_THREADS()
!        IF ( NUMBER_OF_THREADS > 1 ) THEN
!            CALL dfftw_plan_with_nthreads(NUMBER_OF_THREADS)
!        ENDIF
          
        CALL dfftw_plan_dft_c2r (map_1%complex_to_real_in_place,  SIZE ( map_1%sp_group%non_polar_axes ), &
                                 dims, map_1%array, map_1%array, FFTW_PATIENT+FFTW_UNALIGNED)
!                                 dims, map_1%array, map_1%array, FFTW_ESTIMATE)

        WRITE(*,"(' COMPLEX_TO_REAL_PLAN> ', ' plan=', I15, ' size=', I2, ' dims=', 3I4)") &
                  map_1%complex_to_real_in_place,  SIZE ( map_1%sp_group%non_polar_axes ), dims
    END SUBROUTINE complex_to_real_plan_in_place

    SUBROUTINE complex_fft_in_place(map_1)
!
!
!       Purpose:
!       =======
!       Hermitian-to-real FFT transform
!       
!       Also fills in h=0 plane for FFTW
!
        TYPE (map), INTENT(INOUT) :: map_1
!       Local variables:
        TYPE(vector_int)                :: hkl
        INTEGER                         :: jk
        INTEGER                         :: jk1
        INTEGER                         :: jl
        INTEGER                         :: jl1
!       Counters
        INTEGER                         :: ik
        INTEGER                         :: il

!       Fill 0kl plane:
        DO il = map_1%lmin, map_1%lmax
            jl = lmod( il, map_1%nw)
            jl1= lmod(-il, map_1%nw)
            DO ik = map_1%kmin, map_1%kmax
                jk  = lmod( ik, map_1%nv)
                jk1 = lmod(-ik, map_1%nv)
                hkl = (/0,ik,il/)
                IF ( hemi_h ( hkl ) ) THEN
                    map_1%array(0,jk1,jl1) =  map_1%array(0,jk,jl)
                    map_1%array(1,jk1,jl1) = -map_1%array(1,jk,jl)
                ENDIF
            ENDDO
        ENDDO
!       Do complex Fourier transform:
        CALL dfftw_execute_dft_c2r ( map_1%complex_to_real_in_place, map_1%array, map_1%array )

    END SUBROUTINE complex_fft_in_place

    SUBROUTINE destroy_complex_fft_plan( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1

        CALL dfftw_destroy_plan(map_1%complex_to_real_in_place)
    END SUBROUTINE destroy_complex_fft_plan

    SUBROUTINE sine_to_real_plan_in_place ( map_1 )
!
!       Purpose:
!       =======
!       Prepares sine-to-real plan for FFTW3
!
!       Date:          Programmer:               History of changes:
!       ====           ==========                ==================
!       Oct 2003       B.Strokopytov             Original code.
!      
!       Oct 2007       B.Strokopytov             DIMS array added (to avoid Intel compiler warnings).
!       Oct 2007       B.Strokopytov             Cosmetic changes.
!
        TYPE(map),             INTENT(INOUT) :: map_1
        INTEGER,  DIMENSION(3)               :: map_size
        INTEGER,  DIMENSION(3)               :: dims
        CHARACTER(LEN=32)                    :: srname
        INTEGER,  DIMENSION(3)               :: fftw_r2r_kind

!       Set up sine fourier kind:
        fftw_r2r_kind = FFTW_RODFT00

        srname = 'sine_to_real_plan_in_place'

!       Basic checks:
        IF ( .NOT. ALLOCATED ( map_1%sp_group%non_polar_axes ) ) THEN
            CALL die('Programming error: non_polar_axes must be initialized first.',&
                     srname)
        ELSE IF ( SIZE ( map_1%sp_group%non_polar_axes ) == 0 ) THEN
            CALL die('No plan can be created: all axes are polar.', &
                      srname)
        ENDIF

        map_size = map_1%map_size
!        dims     = 0
        dims = map_size(map_1%sp_group%non_polar_axes)

!       Similar pointers => in_place transform:
        IF ( .NOT. ALLOCATED ( map_1 ) ) THEN
            CALL die('Progragramming error. Map_1 has not been allocated...', 'sine_to_real_plan_in_place')
        ENDIF
!        CALL dfftw_plan_dft (map_1%sine_to_real_in_place,  SIZE ( map_1%sp_group%non_polar_axes ), &
!                             dims, map_1%array, map_1%array, FFTW_RODFT00, FFTW_ESTIMATE)
        CALL dfftw_plan_r2r (map_1%sine_to_real_in_place,  SIZE ( map_1%sp_group%non_polar_axes ), &
                             dims, map_1%array, map_1%array, fftw_r2r_kind, FFTW_ESTIMATE)
 
        WRITE(*,*) ' map_1%sine_to_real_in_place=', map_1%sine_to_real_in_place, FFTW_RODFT00
    END SUBROUTINE sine_to_real_plan_in_place

    SUBROUTINE sine_to_real_fft_in_place(map_1)
!
!
!       Purpose:
!       =======
!       Sine-to-real FFTW transform
!       
!       Also fills in h=0 plane for FFTW
!
        TYPE (map), INTENT(INOUT) :: map_1
!       Local variables:
        TYPE(vector_int)                :: hkl
        INTEGER                         :: jk
        INTEGER                         :: jk1
        INTEGER                         :: jl
        INTEGER                         :: jl1
!       Counters
        INTEGER                         :: ik
        INTEGER                         :: il

!       Fill 0kl plane:
        DO il = map_1%lmin, map_1%lmax
            jl = lmod( il, map_1%nw)
            jl1= lmod(-il, map_1%nw)
            DO ik = map_1%kmin, map_1%kmax
                jk  = lmod( ik, map_1%nv)
                jk1 = lmod(-ik, map_1%nv)
                hkl = (/0,ik,il/)
                IF ( hemi_h ( hkl ) ) THEN
                    map_1%array(0,jk1,jl1) =  map_1%array(0,jk,jl)
                    map_1%array(1,jk1,jl1) = -map_1%array(1,jk,jl)
                ENDIF
            ENDDO
        ENDDO

!       Do complex Fourier transform:
        WRITE(*,*) ' before dfftw execute in sine_to_real_fft_in_place'
        CALL dfftw_execute ( map_1%sine_to_real_in_place )
        WRITE(*,*) ' after dfftw execute in sine_to_real_fft_in_place'

    END SUBROUTINE sine_to_real_fft_in_place

    SUBROUTINE destroy_sine_fft_plan( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1

        CALL dfftw_destroy_plan(map_1%sine_to_real_in_place)
    END SUBROUTINE destroy_sine_fft_plan
                                            
    SUBROUTINE allocate_real_map ( map_1, mtz_2, mpeaks, lpaxis, bmin )
!
!       Purpose:
!       =======
!       Allocates map for FFTW calculations
!
!       Note:
!       ====
!       lpaxis array serves as a signal to redefine map dimensions
!       We need this for calculation of T3 translation function
!       for first molecule in polar space groups only.
!       
!       However, in many cases this is undesirable behaviour
!       Subroutine reset_polar_axes should be used to avoid this.
!
!       Date:             Programmer:        History of changes:
!       ====              ==========         ==================
!       Nov 2003          B.Strokopytov      Original FORTRAN 90 code.
!       Oct 2006          B.Strokopytov      Bug corrected (see below).
!       Nov 2007          B.Strokopytov      New ORT_SYM_INVERTED_DEORT array added.
!       Nov 2007          B.Strokopytov      Small cosmetic changes.
!
        TYPE(map),             INTENT(INOUT)         :: map_1
        TYPE(mtz),             INTENT(IN)            :: mtz_2
        INTEGER,               INTENT(IN)            :: mpeaks
        LOGICAL, DIMENSION(3), INTENT(IN), OPTIONAL  :: lpaxis  ! signal to redefine map dimensions
        REAL(KIND=wp),         INTENT(IN), OPTIONAL  :: bmin    ! need to adjust b_scale
!       Local variables:
        INTEGER, DIMENSION(3)                        :: map_size                    
        INTEGER                                      :: m
!
        IF ( mpeaks >= 0 ) THEN
            CALL allocate_array ( map_1%peak_locations, mpeaks )
            CALL allocate_array ( map_1%corr_coef,      mpeaks )
        ENDIF
!
        map_1%map_size = gridres ( mtz_2, map_1%grid_factor )
        map_size       = map_1%map_size          
        map_1%nu       = map_size(1)
        map_1%nv       = map_size(2)
        map_1%nw       = map_size(3)

!       Need this for complex FFT transform:
        map_1%hmin         = mtz_2%hmin
        map_1%hmax         = mtz_2%hmax
        map_1%kmin         = mtz_2%kmin
        map_1%kmax         = mtz_2%kmax
        map_1%lmin         = mtz_2%lmin
        map_1%lmax         = mtz_2%lmax

        map_1%ORT          = mtz_2%ORT
        map_1%DEORT        = mtz_2%DEORT
        map_1%real_tensor  = mtz_2%real_tensor
        map_1%recip_tensor = mtz_2%recip_tensor
        map_1%b_scale      = mtz_2%b_scale
        map_1%badd         = mtz_2%badd
        IF ( PRESENT ( bmin ) ) THEN
            map_1%b_scale = mtz_2%badd - bmin
        ENDIF

!       BUG CORRECTED OCT 2007 BVS:
        map_1%fft_radius   = mtz_2%fft_radius
        map_1%volume       = mtz_2%volume
        map_1%sp_group     = mtz_2%sp_group
        map_1%cell         = mtz_2%cell

        CALL allocate_array ( map_1%SYM_DEORT, map_1%sp_group%number_of_symops )
        CALL allocate_array ( map_1%SYM_ORT_INVERTED, map_1%sp_group%number_of_symops )

!       Added new matrix array NOV 2007 BVS:
        CALL allocate_array ( map_1%SYM_ORT_INVERTED_DEORT, map_1%sp_group%number_of_symops )

!       BUG CORRECTED OCT 2007 BVS (these records were missing):
        DO m = 1, map_1%sp_group%number_of_symops

            map_1%SYM_DEORT(m)        = map_1%sp_group%SYM(m) * map_1%DEORT
            map_1%SYM_ORT_INVERTED(m) = map_1%ORT  * .INV. ( .SYMA. map_1%sp_group%SYM(m) )
            map_1%SYM_ORT_INVERTED_DEORT(m) = map_1%SYM_ORT_INVERTED(m) * map_1%DEORT

            IF ( debug > 30 ) THEN
                CALL print_symop(map_1%sp_group%SYM(m) * .INV. map_1%sp_group%SYM(m))
            ENDIF

        ENDDO
 
!       Check Shannon factor:
        IF ( mtz_2%hmax >= map_1%nu / 2 .OR. mtz_2%kmax >= map_1%nv / 2 .OR. mtz_2%lmax >= map_1%nw / 2 ) THEN

            WRITE(*, "(' ALLOCATE_MAP> ', 'Map size= ', 3I4)") map_size
            WRITE(*, "(' ALLOCATE_MAP> ', '2*hmax, 2*kmax, 2*lmax :',3I4)") &
            2*mtz_2%hmax, 2*mtz_2%kmax, 2*mtz_2%lmax

            CALL die ( 'Grid is too small for HKL data', 'allocate_map' )
        ENDIF

!       Change map dimension according to polarity of corresponding axis:
        IF ( PRESENT ( lpaxis ) ) THEN

            CALL messag('Changing map size according to polar axes definition...', 'allocate_map')

            IF ( lpaxis(1) ) map_1%nu = 1          ! results in 0:1 that's what we need for fftw ( one additional section, nu is odd)
            IF ( lpaxis(2) ) map_1%nv = 1          ! results in 0:0
            IF ( lpaxis(3) ) map_1%nw = 1          ! results in 0:0

!           Redefine map size:
            map_1%map_size = (/map_1%nu, map_1%nv, map_1%nw/)

            IF ( lpaxis(1) ) THEN
                CALL allocate_array ( map_1%array,            2, map_1%nv, map_1%nw, 'map_1%array' )
            ELSE
                CALL allocate_array ( map_1%array, map_1%nu + 2, map_1%nv, map_1%nw, 'map_1%array' )
            ENDIF

        ELSE

!           Allocate 2 additional sections along x-axis as required by FFTw.
!           We use inplace transforms only (memory economy):
            CALL allocate_array ( map_1%array, map_1%nu + 2, map_1%nv, map_1%nw, 'map_1%array' )
        ENDIF

!       Done:
        WRITE(*,"(' ALLOCATE_MAP> ', 'Map allocation completed. Map dimensions= ', 3I4)") map_size 
    END SUBROUTINE allocate_real_map

    SUBROUTINE deallocate_real_map( map_1 )
        TYPE(map), INTENT(INOUT) :: map_1
!
!       Deallocate huge 3D array first:
        IF ( ALLOCATED ( map_1%array ) ) CALL deallocate_array ( map_1%array )

!       Deallocate peak search arrays:
        IF ( ALLOCATED ( map_1%peak_locations ) ) CALL deallocate_array ( map_1%peak_locations )
        IF ( ALLOCATED ( map_1%corr_coef      ) ) CALL deallocate_array ( map_1%corr_coef )
        IF ( ALLOCATED ( map_1%SYM_DEORT      ) ) CALL deallocate_array ( map_1%SYM_DEORT )

        IF ( ALLOCATED ( map_1%SYM_ORT_INVERTED ) ) THEN
            CALL deallocate_array ( map_1%SYM_ORT_INVERTED )
        ENDIF

        IF ( ALLOCATED ( map_1%SYM_ORT_INVERTED_DEORT ) ) THEN
            CALL deallocate_array ( map_1%SYM_ORT_INVERTED_DEORT )
        ENDIF

!       Free space allocated for sp. group information:
        CALL deallocate_space_group ( map_1%sp_group )

    END SUBROUTINE deallocate_real_map

    SUBROUTINE allocate_array_of_real_maps ( maps, mtz_1, mode )
!
!       Purpose:
!       =======
!       Allocates array of maps depending on input mode.
!
!       Date:             Programmer:            History of changes:
!       ====              ==========             ==================
!       Oct 2007          B.Strokopytov          Original code
!
!
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        CHARACTER(LEN=*),                     INTENT(IN)             :: mode
!       Local variables:
        INTEGER                                                      :: number_of_maps
        INTEGER                                                      :: istat 
!       Refinement:
        LOGICAL                                                      :: refine_xyz
        LOGICAL                                                      :: refine_biso
        LOGICAL                                                      :: refine_occ
        LOGICAL                                                      :: refine_Us
!       Counters:
        INTEGER                                                      :: i
        CHARACTER(LEN=32), SAVE                                      :: srname =  'allocate_array_of_real_maps'
        LOGICAL                                                      :: plans_exist
        INTEGER                                                      :: my_first_allocated_map
  
        plans_exist = .FALSE.
        CALL which_params(mode, refine_xyz, refine_biso, refine_occ, refine_Us)

        IF ( .NOT. refine_Us ) THEN
            number_of_maps = 5
        ELSE
            number_of_maps = 11
        ENDIF

        ALLOCATE ( maps(1:number_of_maps), STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat, ' number of maps=', number_of_maps
            CALL die('Failed to allocate array MAPS array. Buy more memory...', srname)
        ENDIF

!       All maps are equal:
        maps(1:number_of_maps)%grid_factor = mtz_1%grid_factor

!       For XYZ refinement allocate map 1-3:
        xyz:IF ( refine_xyz ) THEN
            CALL messag('Allocating maps No. 1,2,3 for XYZ refinement...', srname)
            DO i = 1, 3
                CALL allocate_map(maps(i), mtz_1, 0)
            ENDDO

!           XYZ Plans:
            CALL messag('Preparing FFTW plans for XYZ refinement...', srname)

!            CALL dfftw_init_threads(i)
!             call dfftw_plan_with_nthreads(8)
!
!            IF ( i > 0 ) THEN
!                CALL messag('Threads initialization successfull return value= '//TRIM(int_to_c(i)), srname)
!            ELSE
!                CALL die('O-o-o-o-ps... Failed to initialize threads...', srname)
!            ENDIF

            CALL real_to_complex_plan_in_place(maps(1))
            CALL complex_to_real_plan_in_place(maps(1))

!           Trying to avoid unneccessary estimation of FFTW plans since all maps are equal:
            DO i = 2, 3
                maps(i)%complex_to_real_in_place = maps(1)%complex_to_real_in_place
                maps(i)%real_to_complex_in_place = maps(1)%real_to_complex_in_place
            ENDDO
            plans_exist = .TRUE.
        ENDIF xyz

!       Continue with plans preparation:
 
        IF ( refine_biso ) THEN
            
            CALL messag('Allocating map No. 4 for Biso refinement...', srname)
            CALL allocate_map(maps(4), mtz_1, 0)
            CALL messag ('Preparing FFTW plans for Biso refinement...', srname)
            IF ( plans_exist ) THEN
                maps(4)%real_to_complex_in_place  = maps(1)%real_to_complex_in_place
                maps(4)%complex_to_real_in_place  = maps(1)%complex_to_real_in_place
            ELSE
                CALL real_to_complex_plan_in_place(maps(4))
                CALL complex_to_real_plan_in_place(maps(4))
            ENDIF
        ENDIF

        IF ( refine_occ ) THEN

            CALL messag('Allocating map No. 5 for OCC refinement...', srname)
            CALL allocate_map(maps(5), mtz_1, 0)
            CALL messag ('Preparing FFTW plans for OCC refinement...', srname)
            IF ( plans_exist .AND. ALLOCATED ( maps(1) )  ) THEN
                maps(5)%real_to_complex_in_place  = maps(1)%real_to_complex_in_place
                maps(5)%complex_to_real_in_place  = maps(1)%complex_to_real_in_place
            ELSE IF (  plans_exist .AND. ALLOCATED ( maps(4) ) ) THEN
                maps(5)%real_to_complex_in_place  = maps(4)%real_to_complex_in_place
                maps(5)%complex_to_real_in_place  = maps(4)%complex_to_real_in_place
            ELSE
                CALL real_to_complex_plan_in_place(maps(5))
                CALL complex_to_real_plan_in_place(maps(5))
            ENDIF
        ENDIF

        IF ( refine_Us ) THEN
            CALL messag('Allocating maps No. 6 - 11 for ADP refinement...', srname)
            DO i = 6, 11
                CALL allocate_map(maps(i), mtz_1, 0)
            ENDDO

            CALL messag ('Preparing FFTW plans for Us refinement...', srname)

            IF ( plans_exist ) THEN
                my_first_allocated_map = first_allocated_map ( maps )
                DO i = 6,11
                    maps(i)%real_to_complex_in_place = maps(my_first_allocated_map)%real_to_complex_in_place
                    maps(i)%complex_to_real_in_place = maps(my_first_allocated_map)%complex_to_real_in_place
                ENDDO
            ELSE
                CALL real_to_complex_plan_in_place(maps(6))
                CALL complex_to_real_plan_in_place(maps(6))
                DO i = 7, 11
                    maps(i)%real_to_complex_in_place = maps(6)%real_to_complex_in_place
                    maps(i)%complex_to_real_in_place = maps(6)%complex_to_real_in_place
                ENDDO
            ENDIF
        ENDIF

        CALL messag ('ARRAY OF '//TRIM ( int_to_c ( number_of_maps ) )//' MAPS has been successfully allocated',&
                      srname)

    END SUBROUTINE allocate_array_of_real_maps

    SUBROUTINE allocate_array_of_real_maps_2 ( maps, mtz_1, number_of_maps )
!
!       Purpose:
!       =======
!       Allocates array of maps depending on input mode.
!
!       Date:             Programmer:            History of changes:
!       ====              ==========             ==================
!       Oct 2007          B.Strokopytov          Original code
!
!
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        INTEGER,                              INTENT(IN)             :: number_of_maps
!       Local variables:
        INTEGER                                                      :: istat 
!       Counters:
        INTEGER                                                      :: i
        CHARACTER(LEN=32), SAVE                                      :: srname =  'allocate_array_of_real_maps_2'

        ALLOCATE ( maps(1:number_of_maps), STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) ' istat=', istat, ' number of maps=', number_of_maps
            CALL die('Failed to allocate array MAPS array. Buy more memory...', srname)
        ENDIF

!       All maps are equal:
        maps(1:number_of_maps)%grid_factor = mtz_1%grid_factor

        CALL messag('Allocating maps...', srname)
        DO i = 1, number_of_maps 
            CALL allocate_map(maps(i), mtz_1, 0)
        ENDDO

!       XYZ Plans:
        CALL messag('Preparing FFTW plans...', srname)

        CALL real_to_complex_plan_in_place(maps(1))
        CALL complex_to_real_plan_in_place(maps(1))
        DO i = 2, number_of_maps
            maps(i)%real_to_complex_in_place = maps(1)%real_to_complex_in_place
            maps(i)%complex_to_real_in_place = maps(1)%complex_to_real_in_place
        ENDDO

        CALL messag ('Array of '//TRIM ( int_to_c ( number_of_maps ) )//' maps has been successfully allocated',&
                      srname)

    END SUBROUTINE allocate_array_of_real_maps_2

    SUBROUTINE deallocate_array_of_real_maps ( maps )
!
!       Purpose:
!       =======
!       Allocates array of maps depending on input mode.
!
!       Date:             Programmer:            History of changes:
!       ====              ==========             ==================
!       Oct 2007          B.Strokopytov          Original code
!
!
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: maps
!       Local variables:
        INTEGER                                             :: istat
!       Counters:
        INTEGER                                             :: i
        CHARACTER(LEN=32)                                   :: srname
        INTEGER                                             :: my_first_allocated_map

        srname = 'deallocate_array_of_real_maps'

        IF ( ALLOCATED ( maps ) ) THEN
            my_first_allocated_map = first_allocated_map ( maps )
            DO i = 1, SIZE ( maps )

                IF ( ALLOCATED ( maps(i) ) ) THEN
                    WRITE(*,"(' DEALLOCATE_ARRAY_OF_REAL_MAPS> ', 'Deallocating map #', A)")&
                    TRIM ( int_to_c ( i ) )
                    CALL deallocate_real_map(maps(i))
                ENDIF

            ENDDO

!           No need to destroy plan for each map anymore:            
            CALL destroy_real_fft_plan(maps(my_first_allocated_map))
            CALL destroy_complex_fft_plan(maps(my_first_allocated_map))

            DEALLOCATE ( maps, STAT=istat )
            IF ( istat /= 0 ) THEN
                CALL die('Failed to deallocate MAPS array', srname)
            ENDIF

        ELSE

            CALL warn('Nothing to deallocate.', srname)

        ENDIF
        
    END SUBROUTINE deallocate_array_of_real_maps

    FUNCTION first_allocated_map ( maps )
!
!       Purpose:
!       =======
!       Finds first allocated map in array of maps.
!
        INTEGER                                          :: first_allocated_map
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: maps
!       COunters:
        INTEGER                                          :: i

        DO i = 1, SIZE ( maps )
            IF ( ALLOCATED ( maps(i) ) ) THEN
                first_allocated_map = i
                RETURN
            ENDIF
        ENDDO
        CALL die ('Programming error. No maps have been properly allocated.',&
                  'first_allocated_map')

    END FUNCTION first_allocated_map

    SUBROUTINE apply_bmin_to_array_of_maps(maps, bmin)
        REAL(KIND=wp),                        INTENT(IN)    :: bmin
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: maps
!       Counters:
        INTEGER                                             :: i
        CHARACTER(LEN=32), SAVE                             :: srname = 'apply_bmin_to_array_of_maps'

        WRITE(*,"(' APPLY_BMIN_TO_ARRAY_OF_MAPS> ', ' Bmin=', F6.2, ' subtracting...')")&
        bmin

        IF ( .NOT. ALLOCATED ( maps ) ) THEN
            CALL die('Programming error. MAPS array has not been allocated properly.',&
                     srname)
        ENDIF

        DO i = 1, SIZE ( maps ) 
            IF ( ALLOCATED ( maps(i) ) ) THEN
                WRITE(*,"(' APPLY_BMIN_TO_ARRAY_OF_MAPS> ', ' map #', A,&
               &' Badd(old)=', F8.2,' Badd(new)=', F8.2)")&
                TRIM ( int_to_c ( i ) ), maps(i)%b_scale, maps(i)%badd - bmin
                maps(i)%b_scale = maps(i)%badd - bmin
            ENDIF
        ENDDO

    END SUBROUTINE apply_bmin_to_array_of_maps

END MODULE map_fft
