MODULE hermitian_map_manip
USE constants
USE mtz_io
USE fft_util
USE util
IMPLICIT NONE
!
!    Purpose:
!    =======
!    MAP_HERMITIAN type has been introduced in Oct 2003
!    for expanding set of reflection to hemisphere.
!    Rarely used in FFT calculations since most of
!    the time we deal with real maps and it would rather
!    inconvenient to do this COMPLEX to REAL conversion every
!    time.
!
TYPE :: map_hermitian
    INTEGER                                         :: nu
    INTEGER                                         :: nv
    INTEGER                                         :: nw
    TYPE(vector_int)                                :: map_size
    REAL(KIND=wp)                                   :: grid_factor
    INTEGER                                         :: hmin
    INTEGER                                         :: hmax
    INTEGER                                         :: kmin
    INTEGER                                         :: kmax
    INTEGER                                         :: lmin
    INTEGER                                         :: lmax
    TYPE(space_group)                               :: sp_group
    TYPE(matrix)                                    :: ORT
    TYPE(matrix)                                    :: DEORT
    TYPE(tensor)                                    :: real_tensor
    TYPE(tensor)                                    :: recip_tensor
    REAL(KIND=wp)                                   :: b_scale
    REAL(KIND=wp)                                   :: volume
!   fftw 3.0 plans, INTEGER*8 is used:
    INTEGER(KIND=eb)                                :: real_to_complex_in_place
    INTEGER(KIND=eb)                                :: real_to_complex_out_of_place
    INTEGER(KIND=eb)                                :: complex_to_real_in_place
    INTEGER(KIND=eb)                                :: complex_to_real_out_of_place
!   Map 3D array:
    COMPLEX(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: array 
END TYPE

CONTAINS
    SUBROUTINE reduce_fobs_to_hemisphere ( mtz_1 )
!
!       Purpose:
!       =======
!       Expands complex fobs array to hemisphere
!       Inteded use is phase translation function
!
!       Expands hkl to hkl_in_P1 
!       Recalcalates s_in_P1 for hemisphere
!       Allocates fc_in_P1 array
!
        TYPE(mtz), INTENT(INOUT)   :: mtz_1
!       Local variables:
        TYPE(map_hermitian)        :: map_2
        TYPE(vector_int)           :: hklm
        COMPLEX(KIND=wp)           :: fobs
        REAL(KIND=wp)              :: dphim
        INTEGER,      DIMENSION(3) :: uvw
        REAL(KIND=wp),DIMENSION(3) :: v
!       Counters:
        INTEGER                    :: i
        INTEGER                    :: m
        INTEGER                    :: iu
        INTEGER                    :: iv
        INTEGER                    :: iw
        INTEGER                    :: total_number_in_P1
        INTEGER                    :: total_old

        CALL messag ( ' ', 'reduce_fobs_to_hemisphere' )

        IF ( .NOT. ALLOCATED( mtz_1%phio ) ) THEN
            CALL die ( 'Expansion not possible: phio has not been initialized properly.',&
                       'reduce_fobs_to_hemisphere' )
        ENDIF

!       Allocate complex hermitian map:
        map_2%grid_factor = 3.0_wp   
        CALL allocate_hermitian_map ( map_2, mtz_1  )
 
!       Initialise map:
        map_2%array = (0.0_wp, 0.0_wp)
        total_number_in_P1 = 0
        DO m = 1, mtz_1%sp_group%number_of_primitive_symops
            DO i = 1, mtz_1%current_number_of_reflections
                hklm  =  mtz_1%sp_group%SYM_HKL(m) * mtz_1%hkl(i)
                dphim =  twopi * ( hklm .DOT. .SYMV. mtz_1%sp_group%SYM(m) )
!                dphim =  twopi * ( mtz_1%hkl(i) .DOT. .SYMV. mtz_1%sp_group%SYM(m) )
                fobs  =  mtz_1%fo(i) * phase ( mtz_1%phio(i) - dphim )

!               Multiply by figure-of-merit if appropriate:
                IF ( ALLOCATED ( mtz_1%fomo ) ) THEN
                    fobs = mtz_1%fomo(i) * fobs
                ENDIF

!               Using .HMOD. instead of .HTUMOD. because we are dealing with complex map here:
                IF ( hemi_h ( hklm ) ) THEN
                    uvw  =  hklm .HMOD. map_2%map_size
                ELSE 

!                   Move reflection index to unique hemisphere:  
                    uvw  = -hklm .HMOD. map_2%map_size
                    fobs =  CONJG(fobs)
                ENDIF

                IF ( map_2%array ( uvw(1), uvw(2), uvw(3) ) == ( 0.0_wp, 0.0_wp ) ) THEN
                    map_2%array ( uvw(1), uvw(2), uvw(3) ) = fobs                    
                    total_number_in_P1 = total_number_in_P1 + 1
                ELSE
                    IF ( ABS ( fobs - map_2%array ( uvw(1), uvw(2), uvw(3) ) ) > eps ) THEN
                        WRITE ( *, "(' REDUCE_FOBS_TO_HEMISPHERE> ', 'fobs=(', 2F12.6,&
                               &') map_2( ',I3,',',I3,',',I3,' ) =', 2F12.6)") &
                               fobs, uvw, map_2%array( uvw(1), uvw(2), uvw(3) )
                        WRITE(*,*) ' phio=', mtz_1%phio(i), ' dphim=', dphim, ' delta=', 180.0/pi* (mtz_1%phio(i)-dphim)
                        CALL warn('It seems we are producing rubbish and/or sysabs.', 'reduce_fobs_to_hemisphere')
                    ENDIF
                ENDIF  
            ENDDO
        ENDDO

        WRITE(*,"(' REDUCE_FOBS_TO_HEMISPHERE> ', A,' unique reflections in ', A,' space group',&
                 &' have been expanded to hemisphere containing ', A, ' reflections.')") &
                  TRIM(int_to_c(mtz_1%current_number_of_reflections)),&
                  TRIM(mtz_1%sp_group%space_group_name), TRIM(int_to_c(total_number_in_P1))

        CALL messag ( ' ', 'reduce_fobs_to_hemisphere' )

!       Allocate fobs:
        CALL allocate_array ( mtz_1%fobs, total_number_in_P1 )
!       Allocate fcalc:
        CALL allocate_array ( mtz_1%fc_in_P1, total_number_in_P1 )
 
!       Allocate hkl_in_P1:
        CALL allocate_array ( mtz_1%hkl_in_P1, total_number_in_P1 )

!       Allocate s_squared:
        CALL allocate_array ( mtz_1%s_in_P1, total_number_in_P1 )

        total_old = total_number_in_P1
!       Sorting in the following way : l slowest, k medium, h - fastest

!       Initialize counter again:
        total_number_in_P1 = 0
        DO iw = mtz_1%lmin, mtz_1%lmax
            DO iv = mtz_1%kmin, mtz_1%kmax
                DO iu = mtz_1%hmin, mtz_1%hmax

!                   Choosing only unique hemisphere:
                    hklm = (/ iu, iv, iw /)
                    IF ( hemi_h ( hklm ) ) THEN
                        uvw = hklm .HMOD. map_2%map_size
                    ELSE
                        CYCLE 
                    ENDIF

!                   We are interested only in non-zero entries in the map:
                    IF ( map_2%array ( uvw(1), uvw(2), uvw(3) ) /= (0.0_wp, 0.0_wp) ) THEN
                        total_number_in_P1 = total_number_in_P1 + 1

!                       Check array boundaries:
                        IF ( total_number_in_P1 > SIZE ( mtz_1%hkl_in_P1 ) ) THEN
                            CALL die ( 'Programming error. Total_number_in_P1 exceeds size of hkl_in_P1 array.',&
                                       'reduce_fobs_to_hemisphere' )
                        ENDIF

!                       Form final arrays:
                        mtz_1%hkl_in_P1(total_number_in_P1)  = hklm
                        mtz_1%fobs(total_number_in_P1)       = map_2%array(uvw(1), uvw(2), uvw(3))                        
                        mtz_1%s_in_P1(total_number_in_P1)    = mtz_1%RECIP_TENSOR * hklm

!                       Monitor results:
                        IF ( MOD ( total_number_in_P1, 1000 ) == 1 ) THEN
                            uvw = mtz_1%hkl_in_P1(total_number_in_P1)
                            WRITE ( *, "(' REDUCE_FOBS_TO_HEMISPERE_RESULT> ', ' index=', 3I5, ' s=', F9.6,&
                                     &' m*fobs=',2F9.2)") uvw, &
                                                        mtz_1%s_in_P1(total_number_in_P1), &
                                                        mtz_1%fobs(total_number_in_P1)
                        ENDIF

                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        CALL deallocate_hermitian_map ( map_2 )

        IF ( total_old /= total_number_in_P1 ) THEN
            WRITE(*,*)  ' total_old=', total_old, ' total_number_in_P1=', total_number_in_P1
            CALL die ( 'Programming error. Inconsistent total numbers.', 'reduce_fobs_to_hemisphere' )
        ENDIF

        CALL messag ( 'Done...', 'reduce_fobs_to_hemisphere' )

    END SUBROUTINE reduce_fobs_to_hemisphere

    SUBROUTINE reduce_fo_to_hemisphere ( mtz_1 )
!
!       Purpose:
!       =======
!       Expands arrays to hemispere:
!
!       a) Expands fo into fo_in_P1
!       b) Expands sigfo into sigfo_in_P1
!       c) Expands test into test_in_P1 ( for R-free )
!       d) Recalculates mult to mult_in_P1 from scratch (removed)
!       e) Allocates mtz_1%fc_in_P1
!       f) Allocates and recalculates mtz_1%s_in_P1 from scratch
!       g) Allocates and recalculates mtz_1%hkl_in_P1
!
!       Note:
!       ====
!       Not all arrays above are really necessary in all cases
!       Most important are fo_in_P1, fc_in_P1, s_in_P1, hkl_in_P1
!       Usage of sigo_in_P1, Rfree_flag_in_P1 depends upon particular algorithm
!
!       OCT 2007 this subroutine needs to be made more general for other types of
!                paramaters like weight, HL coefs, etc.
!
!
!       Date:         Programmer:               Description of changes:
!       ====          ==========                ======================
!       Jan 2004      B.Strokopytov             Original code
!       Nov 2005          -"-                   Removed fobs2 array altogether
!       Nov 2005          -"-                   allocate_array is used now where possible 
!       Dec 2005          -"-                   BUG CORRECTED iu,iv,iw -> uvw(1), uvw(2), uvw(3)
!       Sep 2007          -"-                   mult_in_P1 has been removed
!       Oct 2007      B.Strokopytov             logical*1 test replaced with Rfree_flag (integer)
!

        TYPE(mtz), INTENT(INOUT)    :: mtz_1
!       Local variables:
        TYPE(map_hermitian)         :: map_2
        TYPE(vector_int)            :: hklm
        COMPLEX(KIND=wp)            :: fobs
        INTEGER, DIMENSION(3)       :: uvw
!       Eps, centric, multiplicity, etc.:
        LOGICAL                     :: centre
        LOGICAL                     :: sysabs
        REAL(KIND=wp)               :: epsil
        REAL(KIND=wp)               :: phi_centric
        REAL(KIND=wp), DIMENSION(3) :: HO
!       Printing:
        INTEGER                     :: nmod
!       Counters:
        INTEGER                     :: i
        INTEGER                     :: m
        INTEGER                     :: iu
        INTEGER                     :: iv
        INTEGER                     :: iw
        INTEGER                     :: total_number_in_P1
        INTEGER                     :: total_old
        CHARACTER(LEN=32), SAVE     :: srname = 'reduce_fo_to_hemisphere'
        INTEGER, DIMENSION(3)       :: ihkl
        CALL messag(' ', srname)

!       Allocate complex hermitian map (this grid factor is almost arbitrary):
        map_2%grid_factor = 3.0_wp   
        CALL allocate_hermitian_map ( map_2, mtz_1  )

!       Initialize map:
!        map_2%array = CMPLX ( 0.0_wp, 0.0_wp, KIND=wp)
        map_2%array = (0.0_wp, 0.0_wp)

        total_number_in_P1 = 0

        DO m = 1, mtz_1%sp_group%number_of_primitive_symops
            DO i = 1, mtz_1%current_number_of_reflections
                hklm  =  mtz_1%sp_group%SYM_HKL(m) * mtz_1%hkl(i)
                fobs = CMPLX ( mtz_1%fo(i), mtz_1%sigfo(i), KIND=wp )

!               Set fobs negative if Rfree_flag array allocated:                  
!                IF ( ALLOCATED ( mtz_1%Rfree_flag ) ) THEN
!                    IF ( mtz_1%Rfree_flag(i) == 0 ) THEN
!                        fobs = CMPLX ( -ABS(mtz_1%fo(i)), mtz_1%sigfo(i), KIND=wp )
!                    ENDIF
!                ENDIF

!               Using .HMOD. instead of .HTUMOD. because we are dealing with complex map here:
                ihkl = mtz_1%hkl(i)
                uvw  = hklm
                IF ( ALL(ABS(ihkl) == (/0,27,25/) ) )THEN
                   WRITE(*,"(' before hemi orig hkl=',3I4,' m=',I2, ' uvw=',3I4  ' fobs=',2F10.2)") ihkl, m, uvw, fobs
                ENDIF

                IF ( hemi_h ( hklm ) ) THEN
                    uvw  =  hklm .HMOD. map_2%map_size
                ELSE   
                    uvw  = -hklm .HMOD. map_2%map_size
!                    fobs =  CONJG(fobs)
                ENDIF
                IF ( ALL(ABS(ihkl) == (/0,27,25/) ) )THEN
                   WRITE(*,"(' after hemi hkl=',3I4,' m=',I2, ' uvw=', 3I4, ' fobs=',2F10.2)") ihkl, m,  uvw, fobs
                ENDIF
                IF ( map_2%array ( uvw(1), uvw(2), uvw(3) ) == ( 0.0_wp, 0.0_wp ) ) THEN
                    map_2%array(uvw(1), uvw(2), uvw(3) ) = fobs               
                    total_number_in_P1                   = total_number_in_P1 + 1
                ELSE
                    IF ( ABS ( ABS (fobs) - ABS (map_2%array ( uvw(1), uvw(2), uvw(3) )) ) > eps ) THEN
                        WRITE(*,"(' trying to set=',2ES9.2, ' already in=',2ES9.2)") fobs, map_2%array ( uvw(1), uvw(2), uvw(3) ) 
!                        WRITE(*,*) ABS ( ABS (fobs) - ABS (map_2%array ( uvw(1), uvw(2), uvw(3) )) ), eps
                        WRITE (*, "(' REDUCE_FO_TO_HEMISPHERE> ', 'fobs=', 2F12.4,&
                               &' map_2[ ',I3,',',I3,',',I3,'] =', 2F12.4)")&
                               fobs, uvw, map_2%array( uvw(1), uvw(2), uvw(3) )
                        uvw = hklm
                        WRITE(*,*) ' hkl=', uvw, ' m=', m
                        CALL warn('Programming error. It seems we are producing rubbish and/or sysabs.',&
                                 srname)
                    ENDIF
                ENDIF  
            ENDDO
        ENDDO

        WRITE ( *, "(' REDUCE_FO_TO_HEMISPHERE> ', A,' unique reflections in ', A,' space group', &
       &             ' have been expanded to hemisphere containing ', A, ' reflections.')") &
              TRIM(int_to_c(mtz_1%current_number_of_reflections)), &
              TRIM(mtz_1%sp_group%space_group_name),               &
              TRIM(int_to_c(total_number_in_P1))

!       Printing parameter:
        nmod = MAX ( total_number_in_P1 / 20, 1 )
        
        CALL messag('Allocating corresponding arrays... ', &
                    srname)

!       Need to update and clean T3_translation_function routines (BVS SEP 2004):
        CALL allocate_array ( mtz_1%fo_in_P1,         total_number_in_P1, 'fo_in_P1'  )          ! real wp
        CALL allocate_array ( mtz_1%sigfo_in_P1,      total_number_in_P1, 'sigfo_in_P1' )        ! real_wp
!        CALL allocate_array ( mtz_1%mult_in_P1,       total_number_in_P1, 'mult_in_P1' )         ! real_wp
        CALL allocate_array ( mtz_1%fc_in_P1,         total_number_in_P1, 'fc_in_P1' )           ! complex wp
        CALL allocate_array ( mtz_1%hkl_in_P1,        total_number_in_P1, 'hkl_in_P1' )          ! integer vector
        CALL allocate_array ( mtz_1%s_in_P1,          total_number_in_P1, 's_in_P1' )            ! real wp
!        CALL allocate_array ( mtz_1%Rfree_flag_in_P1, total_number_in_P1, 'Rfree_flag_in_P1' )   ! now integer

!       Added Oct 2007 BVS:
        CALL allocate_array ( mtz_1%HO,          total_number_in_P1, 'HO' )            ! real wp
        CALL allocate_array ( mtz_1%KO,          total_number_in_P1, 'KO' )            ! real wp
        CALL allocate_array ( mtz_1%LO,          total_number_in_P1, 'LO' )            ! real wp

        total_old = total_number_in_P1

!       Sorting in the following way - l slowest, k medium, h fastest:
        total_number_in_P1 = 0
        DO iw = mtz_1%lmin, mtz_1%lmax
            DO iv = mtz_1%kmin, mtz_1%kmax
                DO iu = mtz_1%hmin, mtz_1%hmax

                    hklm = (/iu, iv, iw/)

!                   Hermitian map - using .HMOD. instead .HTUMOD. :
                    IF ( hemi_h ( hklm ) ) THEN
                        uvw = hklm .HMOD. map_2%map_size
                    ELSE
                        CYCLE
                    ENDIF

                    nonzero:IF ( map_2%array ( uvw(1), uvw(2), uvw(3) ) /= (0.0_wp, 0.0_wp) ) THEN

                        total_number_in_P1 = total_number_in_P1 + 1

                        IF ( total_number_in_P1 > SIZE ( mtz_1%fo_in_P1 ) ) THEN

                            CALL die('Total_number_in_P1 exceeds size of fobs array.',&
                                     srname)

                        ELSE IF ( total_number_in_P1 > SIZE(mtz_1%hkl_in_P1) ) THEN

                            CALL die('Total_number_in_P1 exceeds size of hkl_in_P1 array.',&
                                     srname)

                        ENDIF

                        mtz_1%hkl_in_P1(total_number_in_P1) = hklm

!                       Added OCT 2007 BVS, HKLs in orthogonal system of coordinates:
                        HO = TRANSPOSE ( mtz_1%DEORT ) * hklm
                        mtz_1%Ho(total_number_in_P1) = Ho(1)
                        mtz_1%Ko(total_number_in_P1) = Ho(2)
                        mtz_1%Lo(total_number_in_P1) = Ho(3)

!                       Need to use ABS to make fo positive, negative was used for test array expansion:
                        mtz_1%fo_in_P1(total_number_in_P1)   = ABS ( REAL ( map_2%array(uvw(1),uvw(2),uvw(3)), KIND=wp ) ) 

!                       BUG corrected Dec 15 2005 BVS:
                        mtz_1%sigfo_in_P1(total_number_in_P1)= AIMAG ( map_2%array(uvw(1),uvw(2),uvw(3)) )          

!                       Calculate multiplicity:
!                        CALL centric_and_epsilon_tests ( hklm, centre, sysabs, epsil, phi_centric, mtz_1%sp_group, &
!                                                         multiplicity=mtz_1%mult_in_P1(total_number_in_P1) ) 

                        mtz_1%s_in_P1(total_number_in_P1) = mtz_1%RECIP_TENSOR * hklm

!                       Rfree flag:
!                        mtz_1%Rfree_flag_in_P1(total_number_in_P1) = 1

!                        IF (  REAL ( map_2%array(uvw(1), uvw(2), uvw(3)) ) < 0.0_wp ) THEN
!                            mtz_1%Rfree_flag_in_P1(total_number_in_P1) = 0
!                        ENDIF

                        IF ( MOD ( total_number_in_P1, nmod ) == 0 ) THEN

!                           *** Need this conversion for printing:
                            uvw = mtz_1%hkl_in_P1(total_number_in_P1)
                            HO(1) = mtz_1%Ho(total_number_in_P1)
                            HO(2) = mtz_1%Ko(total_number_in_P1)
                            HO(3) = mtz_1%Lo(total_number_in_P1)

!                           Monitor results:
                            WRITE (*,"(' REDUCE_FO_TO_HEMISPERE_RESULT> ', 'Ho index=', 3F6.2, ' s=', F8.4, &
                            &' fo=', ES9.2, ' mult(P1)=', ES7.1)" ) &
                            HO,                              &
                            mtz_1%s_in_P1(total_number_in_P1),     &
                            mtz_1%fo_in_P1(total_number_in_P1),    &
                            1.0_wp
!                            mtz_1%Rfree_flag_in_P1(total_number_in_P1)

                        ENDIF

                    ENDIF nonzero

                ENDDO
            ENDDO
        ENDDO

        CALL deallocate_hermitian_map ( map_2 )

        IF ( total_old /= total_number_in_P1 ) THEN
            WRITE(*,*) total_old, total_number_in_P1
            CALL die ( 'Programming error. total_old inconsistent with total_number_in_P1.', 'reduce_to_hemisphere' )
        ENDIF

        CALL messag ( 'Done...', 'reduce_fo_to_hemisphere' )
        CALL messag ( ' ', 'reduce_fo_to_hemisphere' )

    END SUBROUTINE reduce_fo_to_hemisphere

    SUBROUTINE move_fo_to_new_index_list ( mtz_1 )
!
!       Purpose:
!       =======
!       Expands fo associated arrays to hemispere and chooses new ASU or whatever using ideal_hkl list.
!
!       a) Expands fo into fo_in_P1
!       b) Expands sigfo into sigfo_in_P1
!       c) Selects fo according to IDEAL_HKL list
!
!       Note:
!       ====
!       Use of hemi_l function (L>=0) instead of hemi_h
!
!       Date:         Programmer:               Description of changes:
!       ====          ==========                ======================
!       Apr 2010      B.Strokopytov             Original code
!
        TYPE(mtz), INTENT(INOUT)    :: mtz_1
!       Local variables:
        TYPE(map_hermitian)         :: map_2
        TYPE(vector_int)            :: hklm
        COMPLEX(KIND=wp)            :: fobs
        INTEGER, DIMENSION(3)       :: uvw
!       Eps, centric, multiplicity, etc.:
        LOGICAL                     :: centre
        LOGICAL                     :: sysabs
        REAL(KIND=wp)               :: epsil
        REAL(KIND=wp)               :: phi_centric
        REAL(KIND=wp), DIMENSION(3) :: HO
!       Printing:
        INTEGER                     :: nmod
!       Counters:
        INTEGER                     :: i
        INTEGER                     :: m
        INTEGER                     :: iu
        INTEGER                     :: iv
        INTEGER                     :: iw
        INTEGER                     :: total_number_in_P1
        INTEGER                     :: nr
        INTEGER                     :: nr_old
        INTEGER, DIMENSION(3)       :: ihkl
        CHARACTER(LEN=32), SAVE     :: srname = 'move_fo_to_new_index_list'
        CALL messag(' ', srname)

!       Allocate complex hermitian map (this grid factor is almost arbitrary):
        map_2%grid_factor = 3.0_wp   
        CALL allocate_hermitian_map ( map_2, mtz_1  )

!       Initialize map:
        map_2%array = (0.0_wp, 0.0_wp)

        total_number_in_P1 = 0

        DO m = 1, mtz_1%sp_group%number_of_primitive_symops
            DO i = 1, mtz_1%current_number_of_reflections
                hklm  =  mtz_1%sp_group%SYM_HKL(m) * mtz_1%hkl(i)
                fobs = CMPLX ( mtz_1%fo(i), mtz_1%sigfo(i), KIND=wp )

                ihkl = mtz_1%hkl(i)

!               Using .HMOD. instead of .HTUMOD. because we are dealing with complex map here:
                IF ( hemi_h ( hklm ) ) THEN
                    uvw  =  hklm .HMOD. map_2%map_size
                ELSE   
                    uvw  = -hklm .HMOD. map_2%map_size
                ENDIF

                IF ( map_2%array ( uvw(1), uvw(2), uvw(3) ) == ( 0.0_wp, 0.0_wp ) ) THEN
                    map_2%array(uvw(1), uvw(2), uvw(3) ) = fobs               
                    total_number_in_P1                   = total_number_in_P1 + 1
                ELSE
                    IF ( ABS ( ABS (fobs) - ABS (map_2%array ( uvw(1), uvw(2), uvw(3) )) ) > eps ) THEN
                        WRITE(*,"(' trying to set=',2ES9.2, ' already in=',2ES9.2)") fobs, map_2%array ( uvw(1), uvw(2), uvw(3) ) 
                        WRITE (*, "(' MOVE_FO_TO_NEW_INDEX_LIST> ', 'fobs=', 2F12.4,&
                               &' map_2[ ',I3,',',I3,',',I3,'] =', 2F12.4)")&
                               fobs, uvw, map_2%array( uvw(1), uvw(2), uvw(3) )
                        uvw = hklm
                        WRITE(*,*) ' hkl=', uvw, ' m=', m
                        CALL die('Programming error. It seems we are producing rubbish and/or sysabs.',&
                                 srname)
                    ENDIF
                ENDIF  
            ENDDO
        ENDDO

        WRITE ( *, "(' MOVE_FO_TO_NEW_INDEX_LIST> ', A,' unique reflections in ', A,' space group', &
       &             ' have been expanded to hemisphere containing ', A, ' reflections.')") &
              TRIM(int_to_c(mtz_1%current_number_of_reflections)), &
              TRIM(mtz_1%sp_group%space_group_name),               &
              TRIM(int_to_c(total_number_in_P1))

        nr_old = SIZE ( mtz_1%hkl )

!       Printing parameter:
        IF ( ALLOCATED ( mtz_1%ideal_hkl ) ) THEN
            nr = SIZE ( mtz_1%ideal_hkl )
            nmod = MAX ( nr / 20, 1 )
        ELSE
            CALL die('Programming error. Array IDEAL_HKL has not been allocated.', srname)
        ENDIF
!        WRITE(*,*) ' Completeness=', (100.0_wp * nr_old) / nr, ' %'

!       Massive deallocation:
        CALL deallocate_array(mtz_1%hkl)
        CALL deallocate_array(mtz_1%mult)
        CALL deallocate_array(mtz_1%fo)
        CALL deallocate_array(mtz_1%sigfo)
        CALL deallocate_array(mtz_1%s)

!       Reallocate appropriate arrays:
        CALL allocate_array(mtz_1%hkl, nr)
!       Move to new hkl array:
        mtz_1%hkl = mtz_1%ideal_hkl
!       Better to recalculate this bloody multiplicity array:
        CALL allocate_array(mtz_1%mult, nr)
        CALL allocate_array(mtz_1%fo, nr)
        CALL allocate_array(mtz_1%sigfo, nr)        
        CALL allocate_array(mtz_1%s, nr) 

        DO i = 1, nr

            hklm = mtz_1%ideal_hkl(i)
 

!           Hermitian map - using .HMOD. instead .HTUMOD. :
            IF ( hemi_h ( hklm ) ) THEN
                uvw =  hklm .HMOD. map_2%map_size
            ELSE
                uvw = -hklm .HMOD. map_2%map_size
            ENDIF

            mtz_1%fo(i)    = ABS ( REAL ( map_2%array(uvw(1), uvw(2), uvw(3)), KIND=wp ) ) 
            mtz_1%sigfo(i) = AIMAG ( map_2%array(uvw(1),uvw(2),uvw(3)) )
            mtz_1%s(i)     = mtz_1%RECIP_TENSOR * hklm
            CALL centric_and_epsilon_tests ( hklm, centre, sysabs, epsil, phi_centric, mtz_1%sp_group, &
                                             multiplicity=mtz_1%mult(i) )

!           Flag missing reflection:          
            IF ( mtz_1%sigfo(i) == 0.0_wp .AND. mtz_1%fo(i) == 0.0_wp ) THEN
                mtz_1%sigfo(i) = -1.0_wp
            ENDIF


!           Monitor results:
            IF ( MOD (i, nmod) == 0 ) THEN
               HO = hklm
               WRITE (*,"(' MOVE_FO_TO_NEW_INDEX_LIST> ', 'HKL index=', 3F6.0, &
             &            ' fo=', ES9.2, ' sifg(fo)=',  ES9.2, ' mult=', F3.0)" ) &
                           HO, mtz_1%fo(i), mtz_1%sigfo(i), mtz_1%mult(i)
            ENDIF

        ENDDO

!       Reset number of current reflections:
        mtz_1%current_number_of_reflections = SIZE ( mtz_1%hkl )

        CALL deallocate_hermitian_map ( map_2 )
        CALL deallocate_array(mtz_1%ideal_hkl)
        CALL messag('Done...', srname)
        CALL messag(' ', srname )

    END SUBROUTINE move_fo_to_new_index_list

    SUBROUTINE allocate_hermitian_map ( map_1, mtz_2 )
!     
!
!       Purpose:
!       =======
!       Allocates hermitian map 
!
!
!       Note:
!       ===== 
!       1st index is made twice as smaller
!       This is fftw requirement when calling from fortran.
!       This subroutine version is intended for array  expansion to hemisphere
!
!       ************************************************       
!       To use for complex-to-complex Fourier transforms
!       we need complex_map_manip module. Hermitian maps 
!       assume symmetry for complex array which is not 
!       appopriate in general case
!       ************************************************       
!
        TYPE(map_hermitian), INTENT(INOUT) :: map_1
        TYPE(mtz),           INTENT(IN)    :: mtz_2
!       Local variables:
        INTEGER, DIMENSION(3)              :: ext(3)
        INTEGER, DIMENSION(3)              :: map_size
        INTEGER                            :: istat

        map_1%map_size = gridres ( mtz_2, map_1%grid_factor )
        map_size = map_1%map_size           
        map_1%nu = map_size(1)
        map_1%nv = map_size(2)
        map_1%nw = map_size(3)

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
        map_1%sp_group     = mtz_2%sp_group

        WRITE ( *, "(' ALLOCATE_MAP> ', 'Map size : ', 3I4)") map_size

!       Need to modify hmax test ( nu instead of nu/2 ):
        IF ( mtz_2%hmax >= map_1%nu .OR. mtz_2%kmax >= map_1%nv/2 .OR. mtz_2%lmax >= map_1%nw/2 ) THEN
            WRITE ( *, "(' ALLOCATE_MAP> ', 'Map size : ', 3I4)") map_size
            WRITE ( *, "(' ALLOCATE_MAP> ', 'hmax, 2*kmax, 2*lmax :', 3I4)" ) &
            mtz_2%hmax, 2*mtz_2%kmax, 2*mtz_2%lmax
            CALL die ( 'Grid is too small for HKL data', 'allocate_map' )
        ENDIF

        IF ( .NOT. ALLOCATED ( map_1%array ) ) THEN

            ALLOCATE ( map_1%array(0:map_1%nu/2+1, 0:map_1%nv-1, 0:map_1%nw-1 ), STAT = istat)

            IF ( istat /= 0 ) THEN
                WRITE(*,'('' ALLOCATE_HERMITIAN_MAP> '', ''Map dims = '', 3I4)') map_1%nu, map_1%nv, map_1%nw/2 
                CALL die ( ' Failed to allocate this bloody map.', 'allocate_hermitian_map' )
            ELSE
                WRITE(*,'('' ALLOCATE_HERMITIAN_MAP> '', ''Map allocation completed. Map dimensions= '', 3I4)') &
                map_1%nu/2 + 2, map_1%nv, map_1%nw
            ENDIF

        ELSE
            ext = UBOUND(map_1%array) - LBOUND(map_1%array) + 1

            IF ( ext(1) /= map_1%nu/2 + 1 .OR. ext(2) /= map_1%nv .OR. ext(3) /= map_1%nw ) THEN

!               map_size is too small - we need to reallocate the map:
                DEALLOCATE ( map_1%array, STAT = istat )
                IF ( istat /=0 ) CALL die (' Failed to deallocate the map','allocate_hermitian_map')
                ALLOCATE ( map_1%array(0:map_1%nu/2, 0:map_1%nv-1, 0:map_1%nw-1 ), STAT=istat)

                IF ( istat /=0 ) THEN
                    CALL die (' Failed to reallocate the map', 'allocate_hermitian_map')
                ELSE
                    CALL messag( 'Map was successfully reallocated.', 'allocate_hermitian_map')
                ENDIF

            ELSE
                CALL messag( 'No need to allocate a map.', 'allocate_hermitian_map')
                RETURN
            ENDIF

        ENDIF
    END SUBROUTINE allocate_hermitian_map

    SUBROUTINE deallocate_hermitian_map ( map_1 )
!
!       Date:           Programmer:            Description of changes:
!       ====            ==========             ======================
!       Oct 2003        B.Strokopytov          Original code
!       Nov 2005        B.Strokopytov          After revision a long standing bug corrected:
!                                              deallocation of sp_group added
!
        TYPE(map_hermitian), INTENT(INOUT) :: map_1
!       Local variables:
        INTEGER                          :: istat

        IF ( ALLOCATED ( map_1%array ) ) THEN
            DEALLOCATE ( map_1%array, STAT=istat )
            IF ( istat /= 0 ) THEN
                CALL die ( 'Failed to deallocate MAP_1...', 'deallocate_hermitian_map' )
            ENDIF
        ENDIF

!       Bug corrected:
        CALL deallocate_space_group ( map_1%sp_group )

    END SUBROUTINE deallocate_hermitian_map
END MODULE hermitian_map_manip
