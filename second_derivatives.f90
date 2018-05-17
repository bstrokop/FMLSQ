MODULE second_derivatives
USE corr_coef_manip
USE coefs
USE convolution
USE distance_manip
USE extract_fo_fc
USE fast_genden
USE GALAHAD_QPB_double
USE genden
USE sparse_basic
IMPLICIT NONE
CONTAINS
    SUBROUTINE cond_Ax(np, y, x, mtz_1, maps, pdb_2, mode, sk, Q)
        INTEGER,                                     INTENT(IN)           :: np
        REAL(KIND=wp),    DIMENSION(np),             INTENT(OUT)          :: y
        REAL(KIND=wp),    DIMENSION(np),             INTENT(IN)           :: x
        TYPE(mtz),                                   INTENT(INOUT)        :: mtz_1
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                            INTENT(IN)           :: mode
        REAL(KIND=wp),    DIMENSION(:),              INTENT(IN), OPTIONAL :: sk
        TYPE(csr_matrix),                            INTENT(IN), OPTIONAL :: Q
        CHARACTER(LEN=32),                           SAVE                 :: srname='cond_Ax'
!       Local array:
        REAL(KIND=wp), DIMENSION(np)                                      :: r1
!       Calculate product of matrices Q * H * Q:
        IF ( PRESENT ( Q ) ) THEN
            CALL matvec(np, y, x, Q)
            CALL Ax(np, r1, y, mtz_1, maps, pdb_2, mode)
            CALL matvec(np, y, r1, Q)
        ELSE
            CALL die('Programming error. Q matrix must be present.', srname)
        ENDIF 
    END SUBROUTINE cond_Ax

    SUBROUTINE Ax(np, y, x, mtz_1, maps, pdb_2, mode, sk, Q)
!
!       Purpose:
!       =======
!       Calculates matrix/vector product y = Ax using convolution:
!
!       Date:         Programmer:            History of changes:
!       ====          ==========             ==================
!       Oct 2007      B.Strokopytov          Original code
!       Oct 2007      B.Strokopytov          Known bugs corrected 
!       Oct 2007      B.Strokopytov          H2 LOGICAL introduced basically for tests
!       Oct 2007      B.Strokopytov          FC_IN_P1 now used as temporary array
!       Oct 2007      B.Strokopytov          WPHASES and FCQ arrays are candidates for MTZ structure.
!       
        INTEGER,                                     INTENT(IN)           :: np
        REAL(KIND=wp),    DIMENSION(np),             INTENT(OUT)          :: y
        REAL(KIND=wp),    DIMENSION(np),             INTENT(IN)           :: x
        TYPE(mtz),                                   INTENT(INOUT)        :: mtz_1
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
!        TYPE(atom_image), DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: all_atom_images
        CHARACTER(LEN=*),                            INTENT(IN)           :: mode
        REAL(KIND=wp),    DIMENSION(:),              INTENT(IN), OPTIONAL :: sk
        TYPE(csr_matrix),                            INTENT(IN), OPTIONAL :: Q
!       Local variables:
        INTEGER                                                           :: n
        INTEGER                                                           :: number_of_maps
        INTEGER                                                           :: nrefl
        LOGICAL                                                           :: H2
!       Counters:
        INTEGER                                                           :: imap
!       Test:
        INTEGER,          DIMENSION(3)                                    :: uvw
        INTEGER                                                           :: i
!       CPU:
        REAL(KIND=sp)                                                     :: time0
        REAL(KIND=sp)                                                     :: total_time0
        REAL(KIND=sp)                                                     :: fft_time
        REAL(KIND=sp)                                                     :: coef_time1
        REAL(KIND=sp)                                                     :: coef_time2
        REAL(KIND=sp)                                                     :: real_space_time
!       Name:
        CHARACTER(LEN=32),                                           SAVE :: srname = 'Ax'
        INTEGER                                                           :: nout

        total_time0 = SECNDS (0.0)
        fft_time = 0
        real_space_time = 0

!       PDB params:
        n  = SIZE ( pdb_2 )

!       Checkz:
        IF ( np /= COUNT ( pdb_2%refinable ) ) THEN

            WRITE(*,*) np, COUNT ( pdb_2%refinable )
            CALL die('Programming error. Inconsistent array sizes.', srname)

        ENDIF

        IF ( mtz_1%fofc_scale <= 0.0_wp ) THEN

            WRITE(*,*)  mtz_1%fofc_scale
            CALL die('Programming error. Scale has not been initalized...', srname)

        ENDIF

!       X-ray data allocation check:
        IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1 ) ) THEN
            CALL die('Programming error. Array FC_IN_P1 should have been allocated...',&
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL die('Programming error. Array FC_IN_P1_INIT must have been allocated...',&
                     srname)
        ELSE IF ( SIZE ( mtz_1%fc_in_P1 ) /= SIZE ( mtz_1%fc_in_P1_init ) ) THEN
            WRITE(*,*)  SIZE ( mtz_1%fc_in_P1 ), SIZE ( mtz_1%fc_in_P1_init )
            CALL die('Programming error. Inconsistent sizes for FC_IN_P1 and FC_IN_P1_INIT arrays.',&
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fc_temp) ) THEN
            CALL die('O-o-o-o-p-s.Programming error. MTZ_1%FC_TEMP array is not allocated.', srname)
        ELSE IF ( .NOT. ALLOCATED ( pdb_2%all_atom_images ) ) THEN
            CALL die('O-o-o-o-p-s.Programming error. ALL_ATOM_IMAGES array has not been allocated.', srname)
        ENDIF

!       Figure out number of maps:
        number_of_maps = SIZE ( maps )
        nrefl = SIZE ( mtz_1%fc_in_P1 )

        time0 = SECNDS(0.0)

        IF ( PRESENT ( sk ) ) THEN
            CALL generate_density_for_map_array_using_images(maps, pdb_2, sk * x, pdb_2%all_atom_images)
        ELSE
            CALL generate_density_for_map_array_using_images(maps, pdb_2, x, pdb_2%all_atom_images)
        ENDIF

        real_space_time = real_space_time + secnds(time0)
        IF ( debug > 15 ) THEN
            DO i = 1, 10
                uvw = mtz_1%hkl_in_P1(i)
                WRITE(*,"('Ax> ', 'hkl_in_P1: ', 3I4, ' phase=', F8.2)") &
                uvw, 180.0 / pi * phase_of ( mtz_1%fc_in_P1(i) )
            ENDDO
        ENDIF

!       CPU for fft:
        time0 = SECNDS(0.0)

!
!       Real FFT for each allocated map:
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(NUMBER_OF_MAPS, MAPS)
        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN
                CALL real_fft_in_place(maps(imap))
            ENDIF
        ENDDO

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (NUMBER_OF_MAPS, MAPS, MTZ_1, NREFL, DEBUG)
        DO imap = 1, number_of_maps

            IF ( ALLOCATED ( maps(imap) ) ) THEN
!INTEL                CALL real_fft_in_place(maps(imap))
!
!               Cannot calculate coefs before all fcq_in_P1, are at hand:
                CALL fcalc_in_P1_space_group(mtz_1%fcq_in_P1(1:nrefl, imap), mtz_1, maps(imap), xnsym=1)

!               Temporary debugging lines:
                IF ( debug > 15 ) THEN
                    IF ( imap == 5 ) THEN
                        DO i = 1, 10
                            uvw = mtz_1%hkl_in_P1(i)
                            WRITE(*,"('Ax> ', 'hkl_in_P1: ', 3I4, ' FCQ(5)=', 2ES9.2)") &
                            uvw,  mtz_1%fcq_in_P1(i,imap)
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
        fft_time = fft_time + SECNDS(time0)
        H2 = .TRUE.
!       FIXME:
!        H2 = .FALSE.

        time0 = SECNDS(0.0)

! FIXME
! Need to modify Ax_matrix_vector_coefs introducin

! FIXME
! This parallelization looks like a big mess. Desperately need something less ugly. Think it over, Boris (MAR 2010):

!$OMP PARALLEL SECTIONS DEFAULT(PRIVATE) SHARED(mtz_1, nrefl, mode, H2, maps, number_of_maps)
!$OMP SECTION
         imap = 1
         IF ( ALLOCATED ( maps(imap) ) ) THEN
             CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
         ENDIF
!$OMP SECTION
         imap = 2
         IF ( ALLOCATED ( maps(imap) ) ) THEN
             CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
         ENDIF
!$OMP SECTION
         imap = 3
         IF ( ALLOCATED ( maps(imap) ) ) THEN
             CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
         ENDIF
!$OMP SECTION
         imap = 4
         IF ( ALLOCATED ( maps(imap) ) ) THEN
             CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
         ENDIF
!$OMP SECTION
         imap = 5
         IF ( ALLOCATED ( maps(imap) ) ) THEN
             CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
         ENDIF

!        Aniso maps:
!$OMP SECTION
         IF ( number_of_maps > 5 ) THEN
             imap = 6
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                             imap, mode, H2=H2)
             ENDIF
         ENDIF
!$OMP SECTION
         IF ( number_of_maps > 5 ) THEN
             imap = 7
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                             imap, mode, H2=H2)
             ENDIF
         ENDIF
!$OMP SECTION
         IF ( number_of_maps > 5 ) THEN
             imap = 8
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                         imap, mode, H2=H2)
             ENDIF
         ENDIF
!$OMP SECTION
         IF ( number_of_maps == 11 ) THEN
             imap = 9
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                             imap, mode, H2=H2)
             ENDIF
         ENDIF
!$OMP SECTION
         IF ( number_of_maps == 11 ) THEN
             imap = 10
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                             imap, mode, H2=H2)
             ENDIF
         ENDIF
!$OMP SECTION
         IF ( number_of_maps == 11 ) THEN
             imap = 11
             IF ( ALLOCATED ( maps(imap) ) ) THEN
                 CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_temp(1:nrefl,imap),&
                                             imap, mode, H2=H2)
             ENDIF
         ENDIF

!$OMP END PARALLEL SECTIONS

!       Need ALL FSQ_SYM arrays before calling this:
        coef_time1 = SECNDS(time0)
        time0 = SECNDS(0.0)

!       Apply symmetry and determine final coefs:
        CALL matrix_magic_coefs(mtz_1, maps, mode)

        coef_time2 = SECNDS(time0)

        time0 = SECNDS(0.0)
!INTEL
!OMP PARALLEL DO DEFAULT ( PRIVATE ) SHARED ( MAPS, NUMBER_OF_MAPS, DEBUG )
        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN
                CALL complex_fft_in_place ( maps(imap) )
                IF ( debug > 20 ) THEN
                    WRITE(*,"(' Ax> ', 'map type:',I2, ' max=', ES9.2, ' min=', ES9.2)") &
                    imap, MAXVAL( maps(imap) ), MINVAL ( maps(imap) )
                ENDIF
            ENDIF
        ENDDO
!INTEL
!OMP END PARALLEL DO

!       Accumulate time for complex fft:
        fft_time = fft_time + SECNDS(time0)

        time0 = SECNDS(0.0)

        CALL fast_matrix_vector_convolution(y, maps, pdb_2, pdb_2%all_atom_images)

        real_space_time = real_space_time + SECNDS ( time0 )

!   FIXME: Need more sophistication with scaling:
!       Applying final scaling to the matrix-vector product:
!        y = (mtz_1%fofc_scale ** 2 *  mtz_1%sp_group%number_of_symops) * y
        y =  mtz_1%sp_group%number_of_symops * y

!       Scale the product if SK is present:
        IF ( PRESENT ( sk ) ) THEN
            y = sk * y
        ENDIF

!       CPU stats:
        nout = 0
        IF ( nout == 6 ) THEN
            WRITE(*,"(' Ax> Coef time= ', 2F8.2, ' s')") coef_time1, coef_time2
            WRITE(*,"(' Ax> FFT time= ', F8.2, ' s')") fft_time
            WRITE(*,"(' Ax> Real Space Operations time= ', F8.2, ' s')") real_space_time
            WRITE(*,"(1X,40('-'))")
            WRITE(*,"(' Ax> Total time= ', F8.2, ' s')") SECNDS(total_time0)
            WRITE(*,"(1X,40('-'),/)")
        ENDIF
    END SUBROUTINE Ax

    SUBROUTINE Ax_slow(np, y, x, mtz_1, maps, pdb_2, mode, sk)
!
!       Purpose:
!       =======
!       Calculates matrix/vector product y = Ax using convolution:
!
!       Date:         Programmer:            History of changes:
!       ====          ==========             ==================
!       Oct 2007      B.Strokopytov          Original code
!       Oct 2007      B.Strokopytov          Known bugs corrected 
!       Oct 2007      B.Strokopytov          H2 LOGICAL introduced basically for tests
!       Oct 2007      B.Strokopytov          FC_IN_P1 now used as temporary array
!       Oct 2007      B.Strokopytov          WPHASES and FCQ arrays are candidates for MTZ structure.
!       
        INTEGER,                                     INTENT(IN)           :: np
        REAL(KIND=wp),    DIMENSION(np),             INTENT(OUT)          :: y
        REAL(KIND=wp),    DIMENSION(np),             INTENT(IN)           :: x
        TYPE(mtz),                                   INTENT(INOUT)        :: mtz_1
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                            INTENT(IN)           :: mode
        REAL(KIND=wp),    DIMENSION(:),              INTENT(IN), OPTIONAL :: sk
!       Local variables:
        INTEGER                                                           :: n
        INTEGER                                                           :: number_of_maps
        INTEGER                                                           :: nrefl
        LOGICAL                                                           :: H2
!       Counters:
        INTEGER                                                           :: imap
!       Test:
        INTEGER,          DIMENSION(3)                                    :: uvw
        INTEGER                                                           :: i 
!       Name:
        CHARACTER(LEN=32),                                           SAVE :: srname = 'Ax'

!       PDB params:
        n  = SIZE ( pdb_2 )

!       Checkz:
        IF ( np /= COUNT ( pdb_2%refinable ) ) THEN

            WRITE(*,*) np, COUNT ( pdb_2%refinable )
            CALL die('Programming error. Inconsistent array sizes.', srname)

        ENDIF

        IF ( mtz_1%fofc_scale <= 0.0_wp ) THEN

            WRITE(*,*)  mtz_1%fofc_scale
            CALL die('Programming error. Scale has not been initalized...', srname)

        ENDIF

!       X-ray data allocation check:
        IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1 ) ) THEN
            CALL die('Programming error. Array FC_IN_P1 should have been allocated...',&
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL die('Programming error. Array FC_IN_P1_INIT must have been allocated...',&
                     srname)
        ELSE IF ( SIZE ( mtz_1%fc_in_P1 ) /= SIZE ( mtz_1%fc_in_P1_init ) ) THEN
            WRITE(*,*)  SIZE ( mtz_1%fc_in_P1 ), SIZE ( mtz_1%fc_in_P1_init )
            CALL die('Programming error. Inconsistent sizes for FC_IN_P1 and FC_IN_P1_INIT arrays.',&
                     srname)
        ENDIF

!       Figure out number of maps:
        number_of_maps = SIZE ( maps )
        nrefl = SIZE ( mtz_1%fc_in_P1 )

        IF ( PRESENT ( sk ) ) THEN
            CALL generate_density_for_map_array(maps, pdb_2, sk * x)
        ELSE
            CALL generate_density_for_map_array(maps, pdb_2, x)
        ENDIF

!       Save phases before it is too late:
        IF ( .NOT. ALLOCATED ( mtz_1%wphases ) ) THEN
            CALL allocate_array ( mtz_1%wphases, SIZE ( mtz_1%fc_in_P1_init ) )
        ENDIF

!       Using elemental functions here:
        mtz_1%wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

        IF ( debug > 15 ) THEN
            DO i = 1, 10
                uvw = mtz_1%hkl_in_P1(i)
                WRITE(*,"('Ax> ', 'hkl_in_P1: ', 3I4, ' phase init=', F8.2)") &
                uvw, 180.0 / pi * phase_of ( mtz_1%fc_in_P1_init(i) )
            ENDDO
        ENDIF

!       Real FFT for each allocated map:
        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN

                CALL real_fft_in_place(maps(imap))

!               Cannot calculate coefs before all fcq_in_P1, are at hand:
                CALL fcalc_in_P1_space_group(mtz_1%fcq_in_P1(1:nrefl, imap), mtz_1, maps(imap), xnsym=1)

!               Temporary debugging lines:
                IF ( debug > 15 ) THEN
                    IF ( imap == 5 ) THEN
                        DO i = 1, 10
                            uvw = mtz_1%hkl_in_P1(i)
                            WRITE(*,"('Ax> ', 'hkl_in_P1: ', 3I4, ' FCQ(5)=', 2ES9.2)") &
                            uvw,  mtz_1%fcq_in_P1(i,imap) 
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
        ENDDO

        H2 = .TRUE.
!       FIXME:
!        H2 = .FALSE.

!       Here fc_in_P1 will be destroyed hence loop division:
        DO imap = 1, number_of_maps

            IF ( ALLOCATED ( maps(imap) ) ) THEN

!               Using FC_IN_P1 as temporary array to avoid unnecessary allocations/deallocations:
                CALL Ax_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), mtz_1%fc_in_P1, &
                                            imap, mode, H2=H2)


            ENDIF

        ENDDO

!       Need all FSQ_SYM before calling this:
        CALL matrix_magic_coefs(mtz_1, maps, mode)

        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN
                CALL complex_fft_in_place ( maps(imap) )
                IF ( debug > 20 ) THEN
                    WRITE(*,"(' Ax> ', 'map type:',I2, ' max=', ES9.2, ' min=', ES9.2)") &
                    imap, MAXVAL( maps(imap) ), MINVAL ( maps(imap) )
                ENDIF
            ENDIF
        ENDDO

!       No longer need that:
!        CALL deallocate_array ( mtz_1%wphases )

        CALL matrix_vector_convolution(y, maps, pdb_2)
       
!       Scale has been applied to matrix-vector product:
        y = (mtz_1%fofc_scale ** 2 *  mtz_1%sp_group%number_of_symops) * y

!       Scale the product if SK is present:
        IF ( PRESENT ( sk ) ) THEN
            y = sk * y
        ENDIF

    END SUBROUTINE Ax_slow

    SUBROUTINE Bx(np, y, x, mtz_1, maps, pdb_2, mode, sk)
!
!       Purpose:
!       =======
!       Calculates matrix/vector product y = Ax using convolution:
!
!       Date:         Programmer:            History of changes:
!       ====          ==========             ==================
!       Oct 2007      B.Strokopytov          Original code
!       Oct 2007      B.Strokopytov          Known bugs corrected 
!       Oct 2007      B.Strokopytov          H2 LOGICAL introduced basically for tests
!       Oct 2007      B.Strokopytov          FC_IN_P1 now used as temporary array
!       Oct 2007      B.Strokopytov          WPHASES and FCQ arrays are candidates for MTZ structure.
!       
        INTEGER,                                     INTENT(IN)           :: np
        REAL(KIND=wp),    DIMENSION(np),             INTENT(OUT)          :: y
        REAL(KIND=wp),    DIMENSION(np),             INTENT(IN)           :: x
        TYPE(mtz),                                   INTENT(INOUT)        :: mtz_1
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                   INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                            INTENT(IN)           :: mode
        REAL(KIND=wp),    DIMENSION(:),              INTENT(IN), OPTIONAL :: sk
!       Local variables:
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE                       :: wphases 
        INTEGER                                                           :: n
        INTEGER                                                           :: number_of_maps
        INTEGER                                                           :: nrefl
        LOGICAL                                                           :: H2
!       Counters:
        INTEGER                                                           :: imap
!       Test:
        INTEGER,          DIMENSION(3)                                    :: uvw
        INTEGER                                                           :: i 
!       Name:
        CHARACTER(LEN=32),                                           SAVE :: srname = 'Bx'

!       PDB params:
        n  = SIZE ( pdb_2 )
!        np = SIZE ( x )
!       Checkz:
        IF ( np /= COUNT ( pdb_2%refinable ) ) THEN

            WRITE(*,*) np, COUNT ( pdb_2%refinable )
            CALL die('Programming error. Inconsistent array sizes.', srname)

        ENDIF

        IF ( mtz_1%fofc_scale <= 0.0_wp ) THEN

            WRITE(*,*)  mtz_1%fofc_scale
            CALL die('Programming error. Scale has not been initalized...', srname)

        ENDIF

!       X-ray data allocation check:
        IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1 ) ) THEN
            CALL die('Programming error. Array FC_IN_P1 should have been allocated...',&
                     srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL die('Programming error. Array FC_IN_P1_INIT must have been allocated...',&
                     srname)
        ELSE IF ( SIZE ( mtz_1%fc_in_P1 ) /= SIZE ( mtz_1%fc_in_P1_init ) ) THEN
            WRITE(*,*)  SIZE ( mtz_1%fc_in_P1 ), SIZE ( mtz_1%fc_in_P1_init )
            CALL die('Programming error. Inconsistent sizes for FC_IN_P1 and FC_IN_P1_INIT arrays.',&
                     srname)
        ENDIF

!       Figure out number of maps:
        number_of_maps = SIZE ( maps )
        nrefl = SIZE ( mtz_1%fc_in_P1 )

        IF ( PRESENT ( sk ) ) THEN
            CALL generate_density_for_map_array(maps, pdb_2, sk * x)
        ELSE
            CALL generate_density_for_map_array(maps, pdb_2, x)
        ENDIF

!       Save phases before it is too late:
        IF ( .NOT. ALLOCATED ( wphases ) ) THEN
            CALL allocate_array ( wphases, SIZE ( mtz_1%fc_in_P1_init ) )
        ENDIF

!       Using elemental functions here:
        wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

        IF ( debug > 15 ) THEN
            DO i = 1, 10
                uvw = mtz_1%hkl_in_P1(i)
                WRITE(*,"('Bx> ', 'hkl_in_P1: ', 3I4, ' phase=', F8.2)") &
                uvw, 180.0 / pi * phase_of ( mtz_1%fc_in_P1(i) )
            ENDDO
        ENDIF

!       Real FFT for each allocated map:
        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN

                CALL real_fft_in_place(maps(imap))

!               Cannot calculate coefs before all fcq_in_P1, are at hand:
                CALL fcalc_in_P1_space_group(mtz_1%fcq_in_P1(1:nrefl, imap), mtz_1, maps(imap), xnsym=1)

!               Temporary debugging lines:
                IF ( debug > 15 ) THEN
                    IF ( imap == 5 ) THEN
                        DO i = 1, 10
                            uvw = mtz_1%hkl_in_P1(i)
                            WRITE(*,"(' Bx> ', 'hkl_in_P1: ', 3I4, ' FCQ(5)=', 2ES9.2)") &
                            uvw,  mtz_1%fcq_in_P1(i,imap) 
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
        ENDDO

        H2 = .TRUE.

!       Here fc_in_P1 will be destroyed hence loop division:
        DO imap = 1, number_of_maps

            IF ( ALLOCATED ( maps(imap) ) ) THEN
                CALL Bx_matrix_vector_coefs(mtz_1%fcq_sym(1:nrefl,imap), mtz_1, maps(imap), &
                                             wphases, imap)
            ENDIF

        ENDDO

!       Need all FCQ_SYM arrays before calling this:
        CALL matrix_magic_coefs(mtz_1, maps, mode)

        DO imap = 1, number_of_maps
            IF ( ALLOCATED ( maps(imap) ) ) THEN

                CALL complex_fft_in_place ( maps(imap) )
                IF ( debug > 20 ) THEN
                    WRITE(*,"(' Bx> ', 'map type:',I2, ' max=', ES9.2, ' min=', ES9.2)") &
                    imap, MAXVAL ( maps(imap) ), MINVAL ( maps(imap) )
                ENDIF
            ENDIF
        ENDDO

!       Free memory:
        CALL deallocate_array ( wphases )

        CALL convolution_via_second_derivatives(y, maps, pdb_2, mode)
       
!       Apply scale to matrix-vector product:
        y = (mtz_1%fofc_scale ** 2 *  mtz_1%sp_group%number_of_symops) * y

!       Scale the product if SK is present:
        IF ( PRESENT ( sk ) ) THEN
            y = sk * y
        ENDIF

    END SUBROUTINE Bx

    SUBROUTINE xyzbq_wrt_fc ( mtz_1, map_1, imap )
!
!       Purpose:
!       =======
!       Calculates coefficients for calculation of convolution of sum of second derivatives using
!       TRONRUD(1999) approach.
!
!       Input:
!       =====
!
!       Output:
!       ======
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Oct 2007           B.Strokopytov           Highly original code.
!       Oct 2007           B.Strokopytov           SGN array added to account for differing sings of
!                                                  mixed second derivatives
!
!       Oct 2007           B.Strokopytov           DEVELOPMENT OF THIS APPROACH TEMPORARILY FROZEN.
!
!       Note:
!       ====
!
!       a) Contrary to Agarwal and Tronrud we do not care which symmetry we obtain - we are in P1 always. 
!          Hence no need to calculate the second map in a separate step ( see Tronrud(1999) ).
!          We just sum up coefs for H1 and H2 terms here:
!
!       c) No need to have several maps in this subroutine since we only use common parameters for each map: 
!
!       d) Normal matrix only. No second derivatives yet.
!
!       Coefficient addition/subtraction of H1 and H2 terms definetely will work for XYZ business
!       Not so sure about other derivatives. NEEDS further work/research.
!
!
        TYPE(mtz),    INTENT(INOUT)  :: mtz_1
!       Just one map is sufficient:
        TYPE(map),    INTENT(IN)     :: map_1
        INTEGER,      INTENT(IN)     :: imap
!       Local vars:
        REAL(KIND=wp)                :: fft_scale
        INTEGER                      :: nrefl
        REAL(KIND=wp), DIMENSION(10) :: sgn
        CHARACTER(LEN=32)            :: srname
        
        sgn(1:10) =  1.0_wp

!       Looks right but needs a lot of checking:
        sgn(4:5)  = -1.0_wp

        srname = 'xyzbq_wrt_fc'

        nrefl = SIZE ( mtz_1%fc_in_P1 )

!       Checkz:
        IF ( SIZE ( mtz_1%fo_in_P1 ) /= nrefl ) THEN
            CALL messag('Programming error. Inconsistent sizes for fcalc and fo_in_P1',&
                        srname)
        ENDIF

!       Need this to scale fcalc properly otherwise the result will depend on
!       map grid:

        fft_scale = map_1%volume / ( map_1%nu * map_1%nv * map_1%nw )

        IF ( fft_scale <= 0.0_wp ) THEN
            CALL die ( 'Programming error. FFT scale factor is non-positive.',&
                       'patt_overlap_2nd_deriv_wrt_fc' )
        ENDIF

!       Famous Agarwal(1978) - Strokopytov(2008) terms:
        mtz_1%fcq_in_P1(1:nrefl,imap) =                     mtz_1%fcq_in_P1(1:nrefl,imap)   &  ! H1 term
                                      - sgn(imap) * CONJG ( mtz_1%fcq_in_P1(1:nrefl,imap) ) &  ! H2 term
                                      * phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1 ) )         ! H2 term
                       
!       Prepare for anti-alias:
        mtz_1%fcq_in_P1(1:nrefl, imap) = mtz_1%fcq_in_P1(1:nrefl, imap) * EXP ( 0.25_wp * mtz_1%b_scale * mtz_1%s_in_P1 ) 

!       Scale up coefs:
        mtz_1%fcq_in_P1(1:nrefl, imap) = fft_scale * mtz_1%weight_in_P1 * map_1%volume * mtz_1%fcq_in_P1(1:nrefl, imap)

!       And fc_in_P1 survives again. We used only its phase here.

    END SUBROUTINE xyzbq_wrt_fc

    SUBROUTINE tronrud_diagonal( diag, mtz_1, pdb_2)
!
!       Purpose:
!       =======
!       Calculates diagonal of normal matrix using Tronrud method:
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Oct 2007           B.Strokopytov           original code based on D.Tronrud paper.
!
        REAL(KIND=wp),    DIMENSION(:), INTENT(INOUT) :: diag
        TYPE(mtz),                      INTENT(INOUT) :: mtz_1
        TYPE(pdb),                      INTENT(IN)    :: pdb_2
!       Local vars:
        INTEGER                                       :: number_of_maps
        TYPE(map),        DIMENSION(:), ALLOCATABLE   :: maps
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE   :: wphases
        REAL(KIND=wp)                                 :: bmin
!       Counters:
        CHARACTER(LEN=32), SAVE                       :: srname = 'tronrud_diagonal'

        IF ( .NOT. ALLOCATED ( pdb_2 ) ) THEN
            CALL die('Programming error. PDB_2 has not been allocated.', srname)
        ENDIF

!       Allocate just 2 maps:
        CALL allocate_array(maps, mtz_1, 2)
        number_of_maps = SIZE ( maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' TRONRUD_DIAGONAL>',  ' B_SCALE for the 1st map=', ES9.2)") maps(1)%b_scale
            WRITE(*,"(' TRONRUD_DIAGONAL>',  ' B_SCALE for the 2nd map=', ES9.2)") maps(2)%b_scale
        ENDIF

!       Calculate and apply Bmin to all allocated maps:
        bmin = get_bmin ( pdb_2 )
        CALL apply_bmin_to_array_of_maps(maps, bmin)

        IF ( ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL allocate_array(wphases, SIZE (  mtz_1%fc_in_P1_init ) )
        ELSE
            CALL die('Programming error. FC_IN_P1_INIT array has not been allocated.',&
                     srname)
        ENDIF

        wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

!       GOING TO DESTROY mtz_1%fc_in_P1:
        CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(1), wphases, H2=.FALSE.)
        CALL simple_map_coefs_in_P1(mtz_1, maps(1))
        CALL complex_fft_in_place(maps(1))

        WRITE(*,"(' TRONRUD_DIAGONAL> ', 'map 1:', ' min=', ES9.2, ' max=', ES9.2)")&
                   MINVAL( maps(1) ), MAXVAL ( maps(1) )

!       Calculate H2 coefs:
        CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(2), wphases, H2=.TRUE.)
        CALL simple_map_coefs_in_P1(mtz_1, maps(2))
        CALL complex_fft_in_place(maps(2))

        WRITE(*,"(' TRONRUD_DIAGONAL> ', 'map 2:', ' min=', ES9.2, ' max=', ES9.2)")&
        MINVAL( maps(2) ), MAXVAL ( maps(2) )

!       Free memory:
        CALL deallocate_array(wphases)

!       Calculate diagonal via convolution:
        CALL normal_matrix_diagonal(diag, maps, pdb_2)
 
!       Scale normal matrix diagonal:
        diag = (mtz_1%fofc_scale ** 2 * mtz_1%sp_group%number_of_symops) * diag

!       Free map memory:
        CALL deallocate_array ( maps )

    END SUBROUTINE tronrud_diagonal

    SUBROUTINE calculate_dense_matrix(H, mtz_1, pdb_2)
!
!       Purpose:
!       =======
!       Calculates diagonal of normal matrix using Tronrud method:
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Jan 2008           B.Strokopytov           original code based on D.Tronrud paper.
!
        REAL(KIND=wp),    DIMENSION(:,:), INTENT(INOUT) :: H
        TYPE(mtz),                        INTENT(INOUT) :: mtz_1
        TYPE(pdb),                        INTENT(IN)    :: pdb_2
!       Local vars:
        INTEGER                                         :: number_of_maps
        TYPE(map),        DIMENSION(:),   ALLOCATABLE   :: maps
        COMPLEX(KIND=wp), DIMENSION(:),   ALLOCATABLE   :: wphases
        REAL(KIND=wp)                                   :: bmin
!       Counters:
        CHARACTER(LEN=32), SAVE                         :: srname = 'calculate_dense_matrix'

        IF ( .NOT. ALLOCATED ( pdb_2 ) ) THEN
            CALL die('Programming error. PDB_2 has not been allocated.', srname)
        ENDIF

!       Allocate just 2 maps:
        CALL allocate_array(maps, mtz_1, 2)
        number_of_maps = SIZE ( maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' CALCULATE_DENSE_MATRIX>',  ' B_SCALE for the 1st map=', ES9.2)") maps(1)%b_scale
            WRITE(*,"(' CALCULATE_DENSE_MATRIX>',  ' B_SCALE for the 2nd map=', ES9.2)") maps(2)%b_scale
        ENDIF

!       Calculate and apply Bmin to all allocated maps:
        bmin = get_bmin ( pdb_2 )
        CALL apply_bmin_to_array_of_maps(maps, bmin)

        IF ( ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL allocate_array(wphases, SIZE (  mtz_1%fc_in_P1_init ) )
        ELSE
            CALL die('Programming error. FC_IN_P1_INIT array has not been allocated.',&
                     srname)
        ENDIF

        wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

!       GOING TO DESTROY mtz_1%fc_in_P1:
        CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(1), wphases, H2=.FALSE.)
        CALL simple_map_coefs_in_P1(mtz_1, maps(1))
        CALL complex_fft_in_place(maps(1))
        WRITE(*,"(' CALCULATE_DENSE_MATRIX> ', 'map 1:', ' max=', ES9.2, ' min=', ES9.2)")&
        MAXVAL( maps(1) ), MINVAL ( maps(1) )

!       Calculate H2 coefs:
        CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(2), wphases, H2=.TRUE.)
        CALL simple_map_coefs_in_P1(mtz_1, maps(2))
        CALL complex_fft_in_place(maps(2))
        WRITE(*,"(' CALCULATE_DENSE_MATRIX> ', 'map 2:', ' max=', ES9.2, ' min=', ES9.2)")&
        MAXVAL( maps(2) ), MINVAL ( maps(2) )

!       Free memory:
        CALL deallocate_array(wphases)

!       Calculate dense normal matrix via convolution:
        CALL dense_normal_matrix(H, maps, pdb_2)
 
!       Scale normal matrix:
        H = (mtz_1%fofc_scale ** 2 * mtz_1%sp_group%number_of_symops) * H

!       Free map memory:
        CALL deallocate_array ( maps )

    END SUBROUTINE calculate_dense_matrix

    SUBROUTINE calculate_sparse_matrix(H, mtz_1, maps, pdb_2, atom_pairs, diag)
!
!       Purpose:
!       =======
!       Calculates diagonal of normal matrix using Tronrud method:
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Jan 2008           B.Strokopytov           Original code based on CALCULATE_DENSE_MATRIX code
!
        TYPE(coo_matrix),                                 INTENT(INOUT) :: H
        TYPE(mtz),                                        INTENT(INOUT) :: mtz_1
        TYPE(map),        DIMENSION(:),   ALLOCATABLE                   :: maps
        TYPE(pdb),                                        INTENT(INOUT) :: pdb_2
        INTEGER,          DIMENSION(:,:), ALLOCATABLE,    INTENT(IN)    :: atom_pairs
        REAL(KIND=wp),    DIMENSION(:),                   INTENT(INOUT) :: diag
!       Local vars:
        TYPE(map),        DIMENSION(:),   ALLOCATABLE                   :: tronrud_maps
        
        IF ( ALLOCATED ( maps(1) ) .AND. ALLOCATED ( maps(2) ) ) THEN
            CALL sparse_matrix(H, mtz_1, maps, pdb_2, atom_pairs, diag)
        ELSE
            CALL allocate_array(tronrud_maps, mtz_1, 2)
            CALL sparse_matrix(H, mtz_1, tronrud_maps, pdb_2, atom_pairs, diag)
            CALL deallocate_array(tronrud_maps)
        ENDIF

    END SUBROUTINE calculate_sparse_matrix

    SUBROUTINE sparse_matrix(H, mtz_1, maps, pdb_2, atom_pairs, diag)
!
!       Purpose:
!       =======
!       Calculates diagonal of normal matrix using Tronrud method:
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Jan 2008           B.Strokopytov           Original code based on CALCULATE_DENSE_MATRIX code
!
        TYPE(coo_matrix),                                 INTENT(INOUT) :: H
        TYPE(mtz),                                        INTENT(INOUT) :: mtz_1
        TYPE(map),        DIMENSION(:),   ALLOCATABLE                   :: maps
        TYPE(pdb),                                        INTENT(INOUT) :: pdb_2
        INTEGER,          DIMENSION(:,:), ALLOCATABLE,    INTENT(IN)    :: atom_pairs
        REAL(KIND=wp),    DIMENSION(:),                   INTENT(INOUT) :: diag
!       Local vars:
        INTEGER                                                         :: number_of_maps
        REAL(KIND=wp)                                                   :: bmin
!       Current number of non-zero elements in the upper triangular part of sparse matrix:
        INTEGER(KIND=eb)                                                :: nnz
        INTEGER(KIND=eb)                                                :: number_of_atom_pairs
!       Counters:
        INTEGER                                                         :: i
        INTEGER                                                         :: j
        INTEGER                                                         :: imap
        INTEGER(KIND=eb)                                                :: atom_pair
        INTEGER(KIND=eb)                                                :: elem
        INTEGER(KIND=eb)                                                :: l
!       Need this for pointers:
        INTEGER                                                         :: iat
        INTEGER                                                         :: jat
        INTEGER                                                         :: block_size
        INTEGER(KIND=eb), DIMENSION(:), ALLOCATABLE                     :: ind
        INTEGER                                                         :: np
!       Mudified maps:
        INTEGER                                                         :: ityp
        INTEGER                                                         :: jtyp
        CHARACTER(LEN=32), SAVE                                         :: srname = 'calculate_sparse_matrix'
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE                     :: ai 
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE                     :: bi
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE                     :: aj
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE                     :: bj
        INTEGER                                                         :: ngauss
        INTEGER,          DIMENSION(:), ALLOCATABLE                     :: list_of_atom_types
        INTEGER(KIND=eb), DIMENSION(:), ALLOCATABLE                     :: vector
!       CPU:
        REAL(KIND=sp)                                                   :: t0
        REAL(KIND=sp)                                                   :: time0
        REAL(KIND=sp)                                                   :: time1

        t0 = SECNDS(0.0)
        CALL CPU_TIME(time0)

        IF ( .NOT. ALLOCATED ( pdb_2 ) ) THEN
            CALL die('Programming error. PDB_2 has not been allocated.', srname)
        ENDIF

!       Allocate just 2 maps:
!        CALL allocate_array(maps, mtz_1, 2)
        number_of_maps = SIZE ( maps )

        IF ( debug > 3 ) THEN
            WRITE(*,"(' CALCULATE_SPARSE_MATRIX>',  ' B_SCALE for the 1st map=', ES9.2)") maps(1)%b_scale
            WRITE(*,"(' CALCULATE_SPARSE_MATRIX>',  ' B_SCALE for the 2nd map=', ES9.2)") maps(2)%b_scale
        ENDIF

!       Calculate and apply Bmin to all allocated maps:
        bmin = get_bmin ( pdb_2 )
        CALL apply_bmin_to_array_of_maps(maps, bmin)

        IF ( ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            IF ( .NOT. ALLOCATED ( mtz_1%wphases ) ) THEN
                CALL allocate_array(mtz_1%wphases, SIZE (  mtz_1%fc_in_P1_init ) )
            ENDIF
        ELSE
            CALL die('Programming error. FC_IN_P1_INIT array has not been allocated.',&
                     srname)
        ENDIF

! FIXME wphases are supposed to be in mtz_1:
        mtz_1%wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

!       Statistics on atom types
        CALL count_atom_types ( pdb_2, list_of_atom_types )
        number_of_atom_pairs = SIZE ( atom_pairs, DIM=2) 
!        diag = 0.0_wp
        nnz = 0
        H%val = 0
        H%ja = 0
        H%ia = 0
!       Need to predict correct pointers for atom blocks to successfully parallelize:
        CALL allocate_array(ind, number_of_atom_pairs )
        ind(1) = 0
        DO atom_pair = 2, number_of_atom_pairs
            iat = atom_pairs(1, atom_pair-1)
            jat = atom_pairs(2, atom_pair-1)
            IF ( iat /= jat ) THEN
                ind(atom_pair) = ind(atom_pair-1) + COUNT ( pdb_2%refinable(1:11,iat) ) &
                                                  * COUNT ( pdb_2%refinable(1:11,jat) )
            ELSE
                block_size = COUNT ( pdb_2%refinable(1:11,iat) )
                block_size = (block_size * ( block_size + 1 )) / 2
                ind(atom_pair) = ind(atom_pair-1) + block_size
            ENDIF

            IF ( debug > 1000 ) THEN
                WRITE(*,"(' atom pair=',I8, ' iat=', I6, ' jat=', I6, ' ind=', I8)") &
                atom_pair, atom_pairs(1:2, atom_pair),  ind(atom_pair)
            ENDIF

        ENDDO

        DO i = 1, SIZE ( list_of_atom_types )
            ityp = list_of_atom_types(i)
            ngauss =  pdb_2%atomsf_table%record(ityp)%ngauss
            IF ( ngauss == 5 ) THEN
                WRITE(*,"(1X,A4, 5F8.4, 2X, 5F8.4, 4X, I4)")  pdb_2%atomsf_table%record(ityp)%id, &
                      pdb_2%atomsf_table%record(ityp)%a, pdb_2%atomsf_table%record(ityp)%b, ngauss
            ELSE IF ( ngauss == 2 ) THEN
                WRITE(*,"(1X,A4, 2F8.4, 2X, 2F8.4, 4X, I4)")  pdb_2%atomsf_table%record(ityp)%id, &
                      pdb_2%atomsf_table%record(ityp)%a, pdb_2%atomsf_table%record(ityp)%b, ngauss

            ELSE
                WRITE(*,*) ' ngauss=', ngauss
                CALL warn('O-o-o-o-p-s-s. Possible programming error. Incorrect number of gaussians.', srname)
            ENDIF
        ENDDO
 
        outer:DO i = 1, SIZE ( list_of_atom_types )
             ityp = list_of_atom_types(i)
             ngauss =  pdb_2%atomsf_table%record(ityp)%ngauss
             CALL allocate_array(ai, ngauss)
             CALL allocate_array(bi, ngauss)
             ai = pdb_2%atomsf_table%record(ityp)%a
             bi = pdb_2%atomsf_table%record(ityp)%b
             inner:DO j = i, SIZE ( list_of_atom_types )
                jtyp = list_of_atom_types(j)
                ngauss =  pdb_2%atomsf_table%record(jtyp)%ngauss
                CALL allocate_array(aj, ngauss)
                CALL allocate_array(bj, ngauss)
                aj = pdb_2%atomsf_table%record(jtyp)%a
                bj = pdb_2%atomsf_table%record(jtyp)%b

!               Debugging:
                IF ( debug > 1000 ) THEN
                    WRITE(*,"(1X,A4, 5F8.4, 2X, 5F8.4, 4X, I4)")  pdb_2%atomsf_table%record(ityp)%id, &
                              pdb_2%atomsf_table%record(ityp)%a, pdb_2%atomsf_table%record(ityp)%b, ngauss
                    WRITE(*,"(1X,A4, 5F8.4, 2X, 5F8.4, 4X, I4)")  pdb_2%atomsf_table%record(jtyp)%id, &
                              pdb_2%atomsf_table%record(jtyp)%a, pdb_2%atomsf_table%record(jtyp)%b, ngauss
!                    WRITE(*,*) ' Calculating coefs for', ityp, jtyp, ' pairs'
                ENDIF

                CALL prepare_vector_for_pairs ( vector, pdb_2, atom_pairs, ityp, jtyp )

!               This may happen for single atom blocks and in some other cases (one atom type present):
                IF ( SIZE ( vector ) == 0 ) THEN
!                   No need to call anything:
                    CYCLE inner
                ENDIF

                IF ( debug > 1000 ) THEN
                    WRITE(*,"(1X,2I6)" ) atom_pairs(1:2,vector)
                ENDIF

!               GOING TO DESTROY mtz_1%fc_in_P1:
                CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(1), mtz_1%wphases, H2=.FALSE., &
                                                ai=ai, bi=bi, aj=aj, bj=bj)
                CALL simple_map_coefs_in_P1(mtz_1, maps(1))
                IF ( debug > 100 ) THEN
                    WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'map 1:', ' min=', ES9.2, ' max=', ES9.2)")&
                                MINVAL( maps(1) ), MAXVAL ( maps(1) )
                ENDIF

!               Calculate H2 coefs:
                CALL tronrud_matrix_coefs_in_P1(mtz_1%fc_in_P1, mtz_1, maps(2), mtz_1%wphases, H2=.TRUE., &
                                                ai=ai, bi=bi, aj=aj, bj=bj)

                CALL simple_map_coefs_in_P1(mtz_1, maps(2))

! Parallelization will be done by Intel MKL library:
!OMP PARALLEL DO
                DO imap = 1, 2
                    CALL complex_fft_in_place(maps(imap))
                ENDDO
!OMP END PARALLEL DO
                IF ( debug > 100 ) THEN
                    WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'map 2:', ' min=', ES9.2, ' max=', ES9.2)")&
                    MINVAL( maps(2) ), MAXVAL ( maps(2) )
                ENDIF

!               Set special flag for output using nnz:
                IF ( i == SIZE ( list_of_atom_types .AND. i == j ) ) nnz = -1 

!               Calculate sparse normal matrix via convolution:
                CALL fast_sparse_normal_matrix(H, maps, pdb_2, atom_pairs(1:2,vector), nnz, ind(vector))
            ENDDO inner
            CALL deallocate_array(aj)
            CALL deallocate_array(bj)
            CALL deallocate_array(ai)
            CALL deallocate_array(bi)
        ENDDO outer


!       Get last pair:
        iat = atom_pairs(1, number_of_atom_pairs )
        jat = atom_pairs(2, number_of_atom_pairs )
        IF ( iat /= jat ) THEN
            block_size = COUNT ( pdb_2%refinable(1:11,iat) ) * COUNT ( pdb_2%refinable(1:11,jat) )
        ELSE
            block_size = COUNT ( pdb_2%refinable(1:11,iat) )
            block_size = ( block_size * ( block_size + 1 ) ) / 2
        ENDIF

!       Add last block size to last pointer to get proper NNZ:
        nnz = ind ( number_of_atom_pairs ) + block_size

!       Check number of matrix elements again:
        IF ( 2 * nnz - COUNT ( pdb_2%refinable ) /= H%nnz ) THEN
            WRITE(*,*) ' nnz=', nnz
            WRITE(*,*) 2 * nnz - COUNT ( pdb_2%refinable ), H%nnz
            CALL die ('Inconsistent number of matrix elements calculated.', srname)
        ENDIF

!       Check number of refinable parameters:
        IF ( COUNT ( pdb_2%refinable ) /= H%nrow ) THEN
            WRITE(*,*) COUNT ( pdb_2%refinable ),  H%nrow
            CALL die('Inconsistent matrix parameters.', srname)
        ENDIF

!        H%val(1:nnz) = (mtz_1%fofc_scale ** 2 * mtz_1%sp_group%number_of_symops) * H%val(1:nnz)
        H%val(1:nnz) =  mtz_1%sp_group%number_of_symops * H%val(1:nnz)

        np = COUNT ( pdb_2%refinable )
        elem = 0

        DO l = 1, nnz
            i = H%ia(l)
            j = H%ja(l)
            IF ( i == j ) THEN
                elem = elem + 1
!                diag(i) = H%val(l)
            ENDIF
        ENDDO

        WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'Current number of calculated nnz elements=           ', A)" ) &
                    TRIM ( int_to_c ( nnz ) )
        WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'Calculated number of diagonal elements=              ', A)" ) &
                    TRIM ( int_to_c ( elem ) )
        WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'Exact wanted number of diagonal elements=            ', A)" ) &
                    TRIM ( int_to_c ( np ) )
        WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'Expected number of additional non-diagonal elements= ', A)" ) &
                    TRIM ( int_to_c ( nnz - np ) )
        WRITE(*,"(' CALCULATE_SPARSE_MATRIX> ', 'Total expected number of elements=                   ', A)" ) &
                    TRIM ( int_to_c ( 2 * nnz - np ) )



!       diag is no longer needed: FIXME
!        diag = (mtz_1%fofc_scale ** 2 * mtz_1%sp_group%number_of_symops) * diag

!       Debugging:
        IF ( debug > 100 ) THEN
            WRITE(*,*) ' H%val=', H%val(1:10)
            WRITE(*,*) ' diag=', diag
        ENDIF

!       Free memory:
        CALL deallocate_array(vector)
        CALL deallocate_array(list_of_atom_types)
        IF ( ALLOCATED ( ind ) ) CALL deallocate_array(ind)

!        CALL dfftw_cleanup()
        CALL CPU_TIME(time1)
        WRITE(*,"(' CALCULATE SPARSE MATRIX> ', 'Done. Total CPU time=', F10.1, 's Elapsed time=', F9.1, ' s')") &
                  time1 - time0, SECNDS(t0)

    END SUBROUTINE sparse_matrix

    SUBROUTINE agarwal_diagonal( diag, maps, mtz_1, pdb_2)
!
!       Purpose:
!       =======
!       Calculates diagonal of normal matrix using Agarwal method:
!
!       Date:              Programmer:             History of changes:
!       ====               ==========              ==================
!       Oct 2007           B.Strokopytov           original code based on Agarwal's paper...
!
        REAL(KIND=wp),    DIMENSION(:),              INTENT(INOUT) :: diag
        TYPE(map),        DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: maps
        TYPE(mtz),                                   INTENT(INOUT) :: mtz_1
        TYPE(pdb),                                   INTENT(IN)    :: pdb_2
!       Local vars:
        LOGICAL                                                    :: h2
        INTEGER                                                    :: number_of_maps
        INTEGER                                                    :: nrefl
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE                :: wphases
!       Counters:
        INTEGER                                                    :: imap
        INTEGER                                                    :: run 
        CHARACTER(LEN=32), SAVE                                    :: srname = 'agarwal_diagonal'

        IF ( .NOT. ALLOCATED ( pdb_2 ) ) THEN
            CALL die('Programming error. PDB_2 has not been allocated.', srname)
        ENDIF

        nrefl = SIZE ( mtz_1%fc_in_P1 )
        number_of_maps = SIZE ( maps )
        
        IF ( debug > 3 ) THEN
            WRITE(*,"(' AGARWAL_DIAGONAL>',  ' B_SCALE for the 1st map=', ES9.2)") maps(1)%b_scale
            WRITE(*,"(' AGARWAL_DIAGONAL>',  ' B_SCALE for the 2nd map=', ES9.2)") maps(2)%b_scale
        ENDIF

        IF ( ALLOCATED ( mtz_1%fc_in_P1_init ) ) THEN
            CALL allocate_array(wphases, SIZE (  mtz_1%fc_in_P1_init ) )
        ELSE
            CALL die('Programming error. FC_IN_P1_INIT array has not been allocated.',&
                     srname)
        ENDIF

!       Using elemental function here:
        wphases = phase ( 2.0_wp * phase_of ( mtz_1%fc_in_P1_init ) )

        DO run = 0, 1
!        DO run = 0, 0

!           Important to start with H2 = .false. for correct H1+H2 summation:
            h2 = ( run == 1 )

!           This will produce fcq_in_P1 coefs:
            CALL matrix_diagonal_coefs_in_P1(mtz_1, maps, wphases, h2) 
 
            DO imap = 1, number_of_maps
 
                IF ( ALLOCATED ( maps(imap) ) ) THEN 
                    mtz_1%fc_in_P1 = mtz_1%fcq_in_P1(1:nrefl, imap)
                    CALL simple_map_coefs_in_P1(mtz_1, maps(imap))
                    CALL complex_fft_in_place(maps(imap))
                    WRITE(*,"(' AGARWAL_DIAGONAL> ', 'map ', A,':', ' max=', ES9.2, ' min=', ES9.2)")&
                    TRIM ( int_to_c ( imap ) ), MAXVAL ( maps(imap) ), MINVAL ( maps(imap) )
                ENDIF 

            ENDDO

            CALL my_normal_matrix_diagonal(diag, mtz_1, maps, pdb_2, h2)
        ENDDO

!       Free memory:
        CALL deallocate_array(wphases)

!       Scale normal matrix diagonal:
        diag = (mtz_1%fofc_scale ** 2 * mtz_1%sp_group%number_of_symops) * diag

    END SUBROUTINE agarwal_diagonal

    FUNCTION joint_atom_density_value ( r2, aa, bb, joint_occ, joint_b, joint_b_scale )
        REAL(KIND=wp)                             :: joint_atom_density_value
        REAL(KIND=wp),                 INTENT(IN) :: r2
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: aa
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: bb
        REAL(KIND=wp),                 INTENT(IN) :: joint_occ                
        REAL(KIND=wp),                 INTENT(IN) :: joint_b
        REAL(KIND=wp),                 INTENT(IN) :: joint_b_scale          
!       Local variables:
        REAL(KIND=wp)                             :: sum
        REAL(KIND=wp)                             :: ae
        REAL(KIND=wp)                             :: be
!       Counters
        INTEGER                                   :: iscatt
        INTEGER                                   :: jscatt
!
        sum = 0.0_wp
        DO iscatt = 1, 5
            DO jscatt = 1, 5
                ae  = aa(iscatt, jscatt) *  ( SQRT( pi4/ ( bb(iscatt, jscatt) + joint_b + joint_b_scale ) ) ) ** 3
                be  = 4.0_wp * pi_squared / ( bb(iscatt, jscatt)              + joint_b + joint_b_scale )
                sum = sum + joint_occ * ae * EXP(-be * r2)
            ENDDO
        ENDDO
        joint_atom_density_value = sum
    END FUNCTION joint_atom_density_value

    FUNCTION atom_scatt ( pdb_2, iat, s )
        REAL(KIND=wp)             :: atom_scatt
        TYPE(pdb),     INTENT(IN) :: pdb_2
        INTEGER,       INTENT(IN) :: iat
        REAL(KIND=wp), INTENT(IN) :: s
!       Local variables:
        INTEGER                     :: ityp
        REAL(KIND=wp), DIMENSION(5) :: a
        REAL(KIND=wp), DIMENSION(5) :: b
        INTEGER                     :: ngauss
        REAL(KIND=wp)               :: sthol

        sthol = s/4.0_wp
        ityp  = pdb_2%atom_type(iat)
        ngauss = pdb_2%atomsf_table%record(ityp)%ngauss
        a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a
        b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b
        
!       Note the absence of occupancy (it's unity): this is the correct result
        atom_scatt  = SUM ( a(1:ngauss) * EXP ( -b(1:ngauss) * sthol ) * EXP ( -pdb_2%biso(iat)  * sthol ) )
    END FUNCTION atom_scatt   

    SUBROUTINE test_matrix_wrt_occ ( p, my_mtz, my_map, one_atom_pdb, random_pdb, n_test )
        TYPE(QPT_problem_type),      INTENT(IN)    :: p
        TYPE(mtz),                   INTENT(IN)    :: my_mtz
        TYPE(map),                   INTENT(INOUT) :: my_map                
        TYPE(pdb),                   INTENT(IN)    :: one_atom_pdb
        TYPE(pdb),                   INTENT(IN)    :: random_pdb
        INTEGER,                     INTENT(IN)    :: n_test
!       Local variables:
        COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: fcalc_iat
        COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: fcalc_jat
        REAL(KIND=wp),   DIMENSION(3)              :: xyz_iat
        REAL(KIND=wp),   DIMENSION(3)              :: xyz_jat
        INTEGER                                    :: iat
        INTEGER                                    :: jat
        REAL(KIND=wp)                              :: Hij
        REAL(KIND=wp)                              :: ratio
!       Counters
        INTEGER                                    :: pair

!       Checkz:
        IF ( .NOT. ALLOCATED ( one_atom_pdb ) ) THEN
            CALL die ( 'Programming error. one_atom_pdb has not been allocated properly.', 'test_matrix_wrt_occ' )
        ELSE IF ( .NOT. ALLOCATED ( random_pdb ) ) THEN
            CALL die ( 'Programming error. random_pdb has not been allocated properly.', 'test_matrix_wrt_occ' )
        ENDIF

        IF ( SIZE ( one_atom_pdb ) /= 1 ) THEN
            CALL die ( 'Programming error. one_atom_pdb size is not 1.', 'test_matrix_wrt_occ' )
        ELSE IF ( SIZE(one_atom_pdb) > SIZE(random_pdb) ) THEN
            CALL die ( 'Programming error. Size of random pdb is too small.', 'test_matrix_wrt_occ' )
        ENDIF

        IF ( one_atom_pdb%occ(1) /= 1.00_wp ) THEN
            CALL die('Programming error. Incorrect occupancy in one_atom_pdb.', 'test_matrix_wrt_occ')
        ENDIF

        IF ( .NOT. ALLOCATED ( my_map ) ) THEN
            CALL die ( 'Programming error. My_map has not been initialized properly.', 'test_matrix_wrt_occ' )
        ENDIF
!       ==============================================Test values of matrix elements=================================================

!       Allocate vectors:
        CALL allocate_array ( fcalc_iat, SIZE ( my_mtz%fo ) )
        CALL allocate_array ( fcalc_jat, SIZE ( my_mtz%fo ) )

        CALL generate_density ( my_map, one_atom_pdb )
        CALL real_fft_in_place ( my_map )
        DO pair = 1, MIN ( n_test, SIZE(p%H%val) )    ! to be safe
            iat = p%H%col(pair)
            jat = p%H%row(pair)

!           Calculate fcalc's from two single atoms:
            CALL fcalc_in_true_space_group ( fcalc_iat, my_mtz, my_map, my_mtz%DEORT * random_pdb%xyz(iat) )
            CALL fcalc_in_true_space_group ( fcalc_jat, my_mtz, my_map, my_mtz%DEORT * random_pdb%xyz(jat) )

!           Make convolution now:
            Hij = ioic_overlap ( my_mtz%fo, fcalc_iat, fcalc_jat, my_mtz%mult )
            WRITE(*,"(' TEST_MATRIX> ', 'iat=', I8, ' jat=', I8,' ioic_overlap=', ES12.5, ' convo=', ES12.5, ' rat=',ES12.5)")&
            iat, jat, Hij, p%H%val(pair), Hij/p%H%val(pair) ! mult in true sp.group BVS

!           Print elements deviating most from hij which is supposed to be much more accurate:
            ratio = Hij / p%H%val(pair)
            IF ( ABS(ratio) > 2.0_wp .OR. ABS(ratio) < 0.5_wp .OR. ratio < 0.0_wp ) THEN
                xyz_iat = my_mtz%DEORT * random_pdb%xyz(iat)
                xyz_jat = my_mtz%DEORT * random_pdb%xyz(jat)
                WRITE(*,"(' xyz_iat=', 3F7.3, ' xyz_jat=', 3F7.3, ' min_dist=', F9.4)")                    &
                xyz_iat, xyz_jat,                       &
                SQRT(min_dist(my_mtz, my_mtz%DEORT*random_pdb%xyz(iat), my_mtz%DEORT*random_pdb%xyz(jat)))
            ENDIF
        ENDDO                        

        CALL deallocate_array ( fcalc_iat )
        CALL deallocate_array ( fcalc_jat )
        
    END SUBROUTINE test_matrix_wrt_occ

END MODULE second_derivatives
