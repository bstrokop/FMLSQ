MODULE first_derivatives
USE coefs
USE corr_coef_manip
USE convolution, ONLY:first_derivatives_wrt_xyzbqu, simple_convolution
USE extract_fo_fc
USE fast_genden, ONLY: prepare_atom_images
USE fft_util
USE genden 
USE scaling
IMPLICIT NONE
CONTAINS
    SUBROUTINE new_lsq_first_deriv_wrt_all ( g, mtz_1, map_1, binary_map, logical_map, pdb_2, mode, fun, only_fun )
!
!       Purpose:
!       =======
!       Calculates gradients for usual least squares target:
!
!       f(xyz) = SUM w(h)(|Fo|-k*|Fc|)**2
!                 h
!
!       and generates images 
!
!       Date:           Programmer:             History of changes:
!       ====            ==========              ==================
!       OCT 2007        B.Strokopytov           Original code
!
        REAL(KIND=wp),    DIMENSION(:),              INTENT(INOUT)           :: g
        TYPE(mtz),                                   INTENT(INOUT)           :: mtz_1
        TYPE(map),                                   INTENT(INOUT)           :: map_1
        TYPE(map),                                   INTENT(INOUT)           :: binary_map
        TYPE(map_logical),                           INTENT(INOUT)           :: logical_map
        TYPE(pdb),                                   INTENT(INOUT)           :: pdb_2
        CHARACTER(LEN=*),                            INTENT(IN)              :: mode
        REAL(KIND=wp),                               INTENT(INOUT), OPTIONAL :: fun
        LOGICAL,                                     INTENT(IN),    OPTIONAL :: only_fun
!       Local variables:
        REAL(KIND=wp)                                                        :: f
        REAL(KIND=wp)                                                        :: R
        INTEGER,          DIMENSION(3)                                       :: uvw
!       Parameters for anisotropic scaling (including bulk solvent) moved to MTZ structure:
!        REAL(KIND=wp),    DIMENSION(9)                                       :: p
        INTEGER                                                              :: np
!       Counter:
        INTEGER                                                              :: i
        CHARACTER(LEN=32),                                          SAVE     :: srname = 'new_lsq_first_deriv_wrt_all'
        REAL(KIND=sp)                                                        :: t0, t1
        REAL(KIND=wp)                                                        :: gnorm

!       Despite having images at hand it is better to generate density using standard routines (memory economy,
!       see comments in PREPARE_ATOM_IMAGES (fast genden module). Using this approach we will NOT store images of
!       unrefinable atoms (or atoms with zero occupancies). BOOK_LINE_SEARCH will have no difficulties in calling these
!       routines:
        CALL generate_density(map_1, pdb_2)

!! IMPORTANT CORRECTION:
!!        CALL generate_density_using_images(map_1, pdb_2, pdb_2%all_atom_images)

        CALL real_fft_in_place(map_1)
        CALL fcalc_in_P1_space_group(mtz_1%fc_in_P1, mtz_1, map_1)

!       Just for tests:   FIXME REMOVE 2 line below 
        IF ( .NOT. ALLOCATED ( mtz_1%fobs ) ) THEN 
            CALL allocate_array(mtz_1%fobs, SIZE ( mtz_1%hkl ) )
        ENDIF

!       Need this for stats:
        IF ( .NOT. ALLOCATED ( mtz_1%fc ) ) THEN
            CALL allocate_array(mtz_1%fc, SIZE ( mtz_1%hkl ) )
        ENDIF
        CALL fcalc_in_true_space_group(mtz_1%fc, mtz_1, map_1)

!       Save phases for second derivatives calculations:
        IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL allocate_array(mtz_1%fc_in_P1_init, SIZE ( mtz_1%fc_in_P1 ) )
        ENDIF



!       Just copy:
        mtz_1%fc_in_P1_init = mtz_1%fc_in_P1

!       Need to make sure some approximation of scale factor exists:
        mtz_1%fofc_scale = SUM ( mtz_1%fo_in_P1 * ABS ( mtz_1%fc_in_P1 ) ) &
                         / SUM ( REAL ( mtz_1%fc_in_P1 * CONJG ( mtz_1%fc_in_P1 ), KIND=wp) )

!       Initialize bulk solvent/aniso params:
        IF ( .NOT. ALLOCATED ( mtz_1%pbulk ) ) THEN
            CALL allocate_array(mtz_1%pbulk, 9)

!           2.0-2.1 A should be around optimum:
            mtz_1%radius = 2.1_wp

            mtz_1%pbulk(1:3) = (/mtz_1%fofc_scale,0.35_wp,70.0_wp/)
            mtz_1%pbulk(4:9) = 0.0_wp
        ENDIF

!       Determine bulk solvent params, overall ADP and calculate partional contrib. from solvent model:
        CALL scale_aniso(mtz_1%pbulk, binary_map, logical_map, mtz_1, pdb_2)

!       Debugging:
!       CALL scale_aniso_test(mtz_1%pbulk, mtz_1, pdb_2)

!       Add contirbution from bulk solvent partial model:
        mtz_1%fc_in_P1 = mtz_1%fc_in_P1 + mtz_1%fpart_in_P1
        mtz_1%fc_in_P1_init = mtz_1%fc_in_P1


!       Add solvent contribution in true space group:
        mtz_1%fc = mtz_1%fc + mtz_1%fpart

!       Prepare for statistics calcns:
        IF ( .NOT. ALLOCATED ( mtz_1%rres ) ) THEN
            CALL resbins(mtz_1, 20, 50)
        ENDIF

        np = COUNT ( pdb_2%refinable ) + 9
        IF ( .NOT. PRESENT ( only_fun ) ) THEN
            CALL shell_r_factor(mtz_1%fo, mtz_1%fc, mtz_1%sigfo, mtz_1%sc, mtz_1%weight, mtz_1%mult, &
                                mtz_1%bin, mtz_1%rres, np)
        ENDIF

!            test_hkl = (/12,2,0/)
!            DO i = 1, SIZE ( mtz_1%hkl_in_P1 )
!                IF ( mtz_1%hkl(i) == test_hkl ) THEN
!                   WRITE(*,*) ' SUPER TEST:'
!                   WRITE(*,"(3i5,2F10.2) 12,2,0,mtz_1%fo_in_P1(i), mtz_1%sc_in_P1(i) * ABS ( mtz_1%fc_in_P1(i) )
!                ENDIF
!            ENDDO


        IF ( .NOT. ALLOCATED ( mtz_1%wphases ) ) THEN
            CALL allocate_array ( mtz_1%wphases, SIZE ( mtz_1%fc_in_P1_init ) )
        ENDIF

!       Using elemental functions here:
        mtz_1%wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )


        IF ( debug > 22 ) THEN
            DO i = 1, 10
                uvw = mtz_1%hkl_in_P1(i)
                WRITE(*,"(' LSQ_FIRST_DERIV_WRT_ALL> ', 'hkl', 3I4, ' phase=', F8.2, ' fo=',F9.1)")&
                uvw, 180.0 / pi * phase_of ( mtz_1%fc_in_P1(i) ), mtz_1%fo_in_P1(i)
                IF  ( isNaN ( phase_of ( mtz_1%fc_in_P1(i) ) ) ) THEN
                    WRITE(*,*) ' biso=', pdb_2%biso
                    CALL die('Programming error. Undefined phases...', srname)
                ENDIF
            ENDDO
        ENDIF

!       A factor of two is necessary to account for the whole sphere of observations:
        f = 2.0_wp * SUM ( mtz_1%weight_in_P1 &
          * ( mtz_1%fo_in_P1 - mtz_1%sc_in_P1 * ABS ( mtz_1%fc_in_P1 ) ) ** 2 )

        IF ( PRESENT ( fun ) ) THEN
            fun = f
        ENDIF

        IF ( PRESENT ( only_fun ) ) THEN
            IF ( only_fun ) RETURN
        ENDIF

        CALL lsq_first_deriv_wrt_fc(mtz_1, map_1)
        CALL simple_map_coefs_in_P1(mtz_1, map_1)
        CALL complex_fft_in_place(map_1)
        WRITE(*,"(' Map maximum=   ', ES9.2)") MAXVAL ( map_1 ) / mtz_1%volume 

!       Final convolution:
        CALL first_derivatives_wrt_xyzbqu(g, map_1, pdb_2, mode)

!       Scale up with number_of_symops constant:
        g = mtz_1%sp_group%number_of_symops * g
        gnorm = SQRT ( SUM ( g * g ) )
        WRITE(*,"(' NEW_LSQ_FIRST_DERIV_WRT_ALL> ', 'Gradient norm=',ES9.2)") gnorm
        IF ( gnorm == 0.0_wp ) THEN
            CALL die('Programming error. Gradient norm is equal exactly to zero.', srname)
        ENDIF

!       More detailed gradient check:
        DO i = 1, SIZE ( g )
            IF ( g(i) == 0.0_wp ) THEN
                WRITE(*,*) i, g(i)
                CALL warn('Possible programming error. Certain components of gradient vector are equal to zero.',&
                          srname)
            ENDIF
            IF ( isNaN ( g(i)) ) THEN
                WRITE(*,*) i, g(i)
            ENDIF
        ENDDO
 
    END SUBROUTINE new_lsq_first_deriv_wrt_all

    SUBROUTINE generate_ideal_fcalc_in_true_space_group(mtz_1, map_1, binary_map, logical_map, pdb_2)
!
!       Purpose:
!       =======
!       mtz_1%hkl must be 100% complete before calling this routine.
!       Then:
!       The routine can fill missing mtz_1%fo values with mtz_1%fc
!
!       Firstly, it regenerates FC values from current model.
!       Regenerates FPART (both arrays are in true space group ASU).
!
!       Secondly, it produces final FC values and fills in missing FO values 
!       thus artificially increasing correlation coefficients.
!
!       Function
!       f(xyz) = SUM w(h)(|Fo|-k*|Fc|)**2
!                 h
!
!       remains the same but number of observations are now slighly higher,
!       therefore correlation coef. must increase. If this is not the
!       case something is wrong.
!
!
!       Date:           Programmer:             History of changes:
!       ====            ==========              ==================
!       APR 2010        B.Strokopytov           Original code
!       MAY 2010        B.Strokopytov           True space group introduced when generating
!                                               mask in real space. BUG CORRECTED.
!       MAY 2010        B.Strokopytov           Check for artiticiall CC increase added.
!
        TYPE(mtz),                                   INTENT(INOUT)           :: mtz_1
        TYPE(map),                                   INTENT(INOUT)           :: map_1
        TYPE(map),                                   INTENT(INOUT)           :: binary_map
        TYPE(map_logical),                           INTENT(INOUT)           :: logical_map
        TYPE(pdb),                                   INTENT(INOUT)           :: pdb_2
!       Local variables:
        REAL(KIND=wp),    DIMENSION(9)                                       :: p
        REAL(KIND=wp)                                                        :: exp_qform
        REAL(KIND=wp)                                                        :: b_scale 
        REAL(KIND=wp)                                                        :: R
        INTEGER,          DIMENSION(3)                                       :: uvw
        REAL(KIND=wp),    DIMENSION(3)                                       :: v
        REAL(KIND=wp),    DIMENSION(3)                                       :: temp
        REAL(KIND=wp),    DIMENSION(3,3)                                     :: BCART
!       There are some plans to use more sophisticated approach than simple radius:
        REAL(KIND=wp)                                                        :: radius
        REAL(KIND=wp)                                                        :: sc
!       Parameters for anisotropic scaling (including bulk solvent):
        INTEGER                                                              :: np
        INTEGER                                                              :: nr
!       Counters:
        INTEGER                                                              :: i
        INTEGER                                                              :: j 
        CHARACTER(LEN=32),                                          SAVE     :: srname='generate_ideal_fcalc_in_true_space_group'
!        TYPE(vector_int)                                                     :: test_hkl=(/5,12,19/)

!       Set number of reflections for output:
        nr = SIZE ( mtz_1%hkl )

!       Copy parameters:
        p(1:9) = mtz_1%pbulk(1:9)

        CALL generate_density(map_1, pdb_2)
        CALL real_fft_in_place(map_1)

!       MTZ_1%FC may have different size (unless we have 100% completness in the data set):
        IF ( ALLOCATED ( mtz_1%fc ) ) THEN
            CALL deallocate_array ( mtz_1%fc )
        ENDIF

        CALL allocate_array(mtz_1%fc, nr )
        CALL fcalc_in_true_space_group(mtz_1%fc, mtz_1, map_1)

!       Generate constant model density (FIXME radius):
        radius = mtz_1%radius 

!       BUG CORRECTED MAY 2010 BVS: sp_group added to subroutine call,
!       This must be in sync with ANISO_SCALE routine otherwise incorrect results maybe produced:
        CALL generate_constant_density(binary_map, pdb_2, sp_group=mtz_1%sp_group, radius=radius)
        WRITE(*,"(' SCALE_ANISO> ', 'Mask generation with radius=', F5.1,' has been generated')") &
                  radius
        WRITE(*,"(' SCALE_ANISO> ', 'Solvent content=', F5.1, '%')") &
        100.0_wp * COUNT ( binary_map%array == 0.0_wp ) / REAL ( binary_map%nu * binary_map%nv * binary_map%nw )

!       Shrink constant density:
!       CALL shave_constant_density(binary_map, logical_map, pdb_2, neighbours, radius=radius, shrink_rad=shrink_rad)

!       Calculate partial structure factors in P1:
        CALL real_fft_in_place(binary_map)

!       Reallocate fpart (it has changed its dimension):
        IF ( ALLOCATED ( mtz_1%fpart ) ) THEN
            CALL deallocate_array ( mtz_1%fpart )
        ENDIF
        CALL allocate_array(mtz_1%fpart, nr)

!       Extract fpart using current hkl list:
        b_scale = mtz_1%b_scale
        mtz_1%b_scale = 0.0_wp
        CALL fcalc_in_true_space_group(mtz_1%fpart, mtz_1, binary_map)
        mtz_1%b_scale = b_scale

!       Swap phases:
        mtz_1%fpart = -mtz_1%fpart

!       Determine bulk solvent params, overall ADP and calculate partional contrib. from solvent model:
!        CALL scale_aniso(p, binary_map, logical_map, mtz_1, pdb_2)

!       Report bulk/aniso parameters:
        WRITE(*,"(' GENERATE_IDEAL_FCALC> ', 'Overall scale :', F10.4, ' Bulk solvent scale:', F8.4, &
                 &' Bulk solvent temp. factor:', F8.2)") p(1:3)
        WRITE(*,"(' GENERATE_IDEAL_FCALC> ', 'Overall temp. factor:', F8.2)") SUM ( p(4:6) ) / 3.0_wp
        WRITE(*,"(' GENERATE_IDEAL_FCALC> ', 'Overall ADP:', 6F10.4)") p(4:6) - SUM ( p(4:6) ) / 3.0_wp, p(7:9)


        BCART = .CONVERT. p(4:9)
!        IF( debug > 10000 ) THEN
            DO i = 1, 3
                WRITE(*,"(3F6.2)") (BCART(i,j),j=1,3)
            ENDDO
!        ENDIF

!       Reallocate mtz_1%sc:
        IF ( ALLOCATED ( mtz_1%sc ) ) THEN
            CALL deallocate_array(mtz_1%sc)
        ENDIF
        CALL allocate_array(mtz_1%sc, nr)

!       Loop over all reflections:
        DO i = 1, nr
            mtz_1%fpart(i) = p(2) * mtz_1%fpart(i) * EXP ( -0.25_wp * p(3) * mtz_1%s(i) )
            v = TRANSPOSE ( mtz_1%DEORT ) * mtz_1%hkl(i)
            temp = MATMUL ( BCART, v )
            exp_qform = EXP ( -0.25_wp * DOT_PRODUCT ( v, temp ) )

!           Save overall scale multiplied by aniso contribution:
            mtz_1%sc(i) = p(1) * exp_qform
        ENDDO
         
!       Add bulk solvent contribution and apply scaling:
        mtz_1%fc =  mtz_1%sc * (mtz_1%fc + mtz_1%fpart)

        DO i = 1, nr 
            IF ( mtz_1%sigfo(i) < 0.0_wp ) THEN
                mtz_1%fo(i) = ABS ( mtz_1%fc(i) )
            ENDIF
        ENDDO

!       Global check:
        IF ( debug > 10000 ) THEN
        DO i = 1, SIZE(mtz_1%fo)
            uvw = mtz_1%hkl(i)
            WRITE(*,"(3I5,2ES11.4,4X,2ES12.4)") uvw, mtz_1%fo(i), ABS(mtz_1%fc(i)), mtz_1%sc(i), ABS(mtz_1%fpart(i))
            DO j = 1, SIZE(mtz_1%hkl_in_P1)
                IF ( mtz_1%hkl_in_P1(j) == mtz_1%hkl(i)) THEN 
                    WRITE(*,"(15X,2ES11.4,4x,2ES12.4)") mtz_1%fo_in_P1(j), mtz_1%sc_in_P1(j)*ABS(mtz_1%fc_in_P1_init(j)),&
                                                            mtz_1%sc_in_P1(j),  ABS(mtz_1%fpart_in_P1(j))
                    IF ( ABS(mtz_1%sc(i)*ABS(mtz_1%fc_in_P1_init(j)) - ABS(mtz_1%fc(i))) > eps ) THEN
                        CALL die('Problems. problems.',srname)
                    ENDIF
                    EXIT
                ENDIF
            ENDDO
            IF ( j == 1 + SIZE(mtz_1%Hkl_in_p1) )THEN
                WRITE(*,*) ' MIssing STUFF'
            ENDIF
        ENDDO
        ENDIF

!       Calculate sums just for test:
        sc =  SUM(mtz_1%fo); R= SUM(ABS(mtz_1%fc))
        WRITE(*,"(' Sums (fo/fc)=',3ES12.4)") sc, R, sc / R
        mtz_1%fc_in_P1 = mtz_1%sc_in_P1 * mtz_1%fc_in_P1_init

        sc = fofc_corr_coef(mtz_1%fo_in_P1, mtz_1%fc_in_P1)
        WRITE(*,"(' Check final CORR. COEF.(P1)=',F10.5)") sc

!       Do not need scale factors anymore:        
        mtz_1%sc = 1.0_wp
        R = r_factor(mtz_1%fo, mtz_1%fc, mtz_1%sc, mtz_1%mult)
        WRITE(*,"(' Check final R-value(mult)=',F10.5)") R

        R = r_factor(mtz_1%fo, mtz_1%fc, mtz_1%sc)
        WRITE(*,"(' Check final R-value(no mult)=',F10.5)") R
        R = fofc_corr_coef(mtz_1%fo, mtz_1%fc, mtz_1%mult)
        WRITE(*,"(' Check final CORR. COEF.(mult)=',F10.5)") R

!       R should always be larger then sc otherwise something is incorrect:
        IF ( R < sc ) THEN
            CALL warn('O-o-o-o-o-p-s-s. Possible programming error.', srname)
        ENDIF

    END SUBROUTINE generate_ideal_fcalc_in_true_space_group

    SUBROUTINE lsq_first_deriv_wrt_fc(mtz_1, map_1)
!
!       Purpose:
!       =======
!       Calculates matrix of second derivatives for SUM w(h)(Fo-Fc)**2  function wrt Fcalc
!       This matrix IS DIAGONAL due to simplicity of the expression for this function.
!
!       Date:           Programmer:                 History of changes:
!       ====            ==========                  ==================
!       Oct 2007        B.Strokopytov               Original code
!
        TYPE(mtz),        INTENT(INOUT) :: mtz_1
        TYPE(map),        INTENT(INOUT) :: map_1
!       Local
        REAL(KIND=wp)                   :: fft_scale
        REAL(KIND=wp)                   :: fo2_aver
!       Test:
        INTEGER                         :: i
        INTEGER, DIMENSION(3)           :: uvw
        CHARACTER(LEN=32), SAVE         :: srname = 'lsq_first_deriv_wrt_fc'

!       Checkz:
        IF ( SIZE ( mtz_1%fc_in_P1 ) /= SIZE ( mtz_1%fo_in_P1 ) ) THEN
            WRITE(*,*)  SIZE ( mtz_1%fc_in_P1 ), SIZE ( mtz_1%fo_in_P1 )
            CALL messag('Programming error. Inconsistent sizes for FC_IN_P1 and FO_IN_P1',&
                        srname)
        ENDIF

!       Need this to scale coefficients properly otherwise the result will depend on map grid:
        fft_scale = map_1%volume / ( map_1%nu * map_1%nv * map_1%nw )

        IF ( fft_scale <= 0.0_wp ) THEN
            WRITE(*,*) ' fft_scale=', fft_scale
            CALL die('Programming error. FFT scale factor is non-positive.', srname)

!       FIXME. Better to check mtz_1%sc_in_P1:
        ELSE IF ( mtz_1%fofc_scale <= 0.0_wp ) THEN
            WRITE(*,*) ' fofc_scale=',  mtz_1%fofc_scale
            CALL die('Programming error. Scale factor is undefined.',&
                     srname)
        ENDIF

!       Multiply (|Fo| - |Fc|) by phase of Fc, apply k and fft_scale:
        mtz_1%fc_in_P1 = (2.0_wp * fft_scale) * mtz_1%sc_in_P1                              &
                       * ( mtz_1%fo_in_P1 - mtz_1%sc_in_P1 * ABS (  mtz_1%fc_in_P1_init ) ) &
                       * mtz_1%fc_in_P1_init / ABS (mtz_1%fc_in_P1_init )

!       Apply weighting scheme:
        IF ( ALLOCATED ( mtz_1%weight_in_P1 ) ) THEN
            mtz_1%fc_in_P1 = mtz_1%weight_in_P1 * mtz_1%fc_in_P1
        ENDIF

!       Apply MTZ_1%B_SCALE (Badd) for derivative (antialias):
        mtz_1%fc_in_P1 = mtz_1%fc_in_P1 * EXP ( 0.25_wp * mtz_1%b_scale * mtz_1%s_in_P1 )

!       Volume scale up:
        mtz_1%fc_in_P1 = map_1%volume * mtz_1%fc_in_P1

        IF ( debug > 15 ) THEN
            WRITE(*,*) ' Test for phases survival...'
            DO i = 1, 10
                uvw = mtz_1%hkl_in_P1(i)
                WRITE(*,"('LSQ_FIRST_DERIVATIVES_WRT_FC> ', 'hkl_in_P1: ', 3I4, ' phase=', F8.2)") &
                uvw, 180.0 / pi * phase_of ( mtz_1%fc_in_P1_init(i) )
            ENDDO
        ENDIF                      
    END SUBROUTINE lsq_first_deriv_wrt_fc

    SUBROUTINE test_grad_wrt_all ( p, my_mtz, my_map, random_pdb, f0 )
!
!       THIS WILL WORK ONLY IN MODE='XYZBQ' WITH 5*natom variables
!            
!

        REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: p
        TYPE(mtz),                   INTENT(INOUT) :: my_mtz
        TYPE(map),                   INTENT(INOUT) :: my_map
        TYPE(pdb),                   INTENT(INOUT) :: random_pdb
        REAL(KIND=wp),               INTENT(IN)    :: f0
!       Local variables:
        REAL(KIND=wp)                              :: delta_x
        REAL(KIND=wp)                              :: fdx
        TYPE(vector)                               :: dxyz
        REAL(KIND=wp), DIMENSION(3)                :: delta_xyz
        INTEGER                                    :: natom
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: k 
        INTEGER                                    :: indxyz

!       Makes no sense to test smth less than 10**-6 due to deteriorating accuracy of calculations:
        delta_x = 0.1_wp
        DO WHILE ( delta_x >= 0.1 ** 8 )

            delta_x = delta_x * 0.1_wp
            WRITE(*,"(' TEST_GRAD_WRT_x> ', 'delta_x=', ES12.5)") delta_x

            natom = SIZE ( random_pdb )

!           Test occupancy derivatives for all atoms:
            DO i = 1, natom
                DO k = 1, 5

!               Change occupancy by a small amount:
                delta_xyz = 0.0_wp
                IF ( k <=3 ) THEN
                   delta_xyz(k) = delta_x
                   dxyz = delta_xyz
                   random_pdb%xyz(i) = random_pdb%xyz(i) + dxyz
                ELSE IF ( k == 4 ) THEN
                    random_pdb%biso(i) = random_pdb%biso(i) + delta_x
                ELSE IF ( k == 5 ) THEN
                    random_pdb%occ(i) = random_pdb%occ(i) + delta_x
                ENDIF

                CALL generate_density ( my_map, random_pdb )
                CALL real_fft_in_place ( my_map )
                CALL fcalc_in_P1_space_group ( my_mtz%fc_in_P1, my_mtz, my_map )
                indxyz = 5 * (i-1) + k

!               Gradient is calculated using the whole sphere of reflection
!               Function is calculated using the hemi-sphere hence factor of 2 to sum the whole sphere up:
                fdx = 2.0_wp * SUM ( ( my_mtz%fo_in_P1 - my_mtz%fofc_scale * ABS ( my_mtz%fc_in_P1 ) ) ** 2 )
                WRITE ( *, "(' TEST_GRAD_WRT_OCC> ', I4, ' grad num=', ES12.5, ' grad=', ES12.5, ' ratio=', ES12.5)" ) &
                indxyz, (fdx - f0) / delta_x, p(indxyz), (fdx - f0) / delta_x / ( p(indxyz))

!                   Set params back to initial values:
                    IF ( k <= 3 ) THEN
                        random_pdb%xyz(i) = random_pdb%xyz(i) - dxyz
                    ELSE IF ( k == 4 ) THEN
                        random_pdb%biso(i) = random_pdb%biso(i) - delta_x
                    ELSE IF ( k == 5 ) THEN
                        random_pdb%occ(i) = random_pdb%occ(i) - delta_x
                    ENDIF
                ENDDO

            ENDDO
        ENDDO
    END SUBROUTINE test_grad_wrt_all

END MODULE first_derivatives
