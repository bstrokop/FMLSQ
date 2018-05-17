MODULE fast_distance
USE basic_adjacency
USE basic_symmetry_operations
USE basic_pdb_manip
USE fail
USE fft_util
USE distance_manip
USE sorting_facilities
USE util
USE vectors
IMPLICIT NONE
TYPE :: cube
    TYPE(vector),     DIMENSION(:), ALLOCATABLE :: xyz
    INTEGER,          DIMENSION(:), ALLOCATABLE :: atom_number
!   Introducing symmetry operator:
    TYPE(symop),      DIMENSION(:), ALLOCATABLE :: atom_symop 
END TYPE
INTERFACE allocate_array
    MODULE PROCEDURE allocate_3d_cube_array
END INTERFACE
INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_3d_cube_array
END INTERFACE
CONTAINS
    SUBROUTINE allocate_3d_cube_array ( cubes, box_grid, cube_size )
        TYPE(cube), DIMENSION(:,:,:), ALLOCATABLE :: cubes
        INTEGER,    DIMENSION(:)                  :: box_grid
        INTEGER,    DIMENSION(:,:,:), ALLOCATABLE :: cube_size
!       Local variables:
        INTEGER                                   :: istat
!       Counters:
        INTEGER                                   :: iu
        INTEGER                                   :: iv
        INTEGER                                   :: iw

!       Checkz:
        IF ( SIZE ( box_grid ) /= 3 ) THEN
            CALL die('Programming error. Size of box_grid is supposed to be 3', 'allocate_3d_cube_array')
        ENDIF

        IF ( ALLOCATED ( cubes ) ) THEN
            CALL deallocate_array ( cubes )
        ENDIF
        ALLOCATE ( cubes(0:box_grid(1)-1,0:box_grid(2)-1,0:box_grid(3)-1), STAT=istat )

        IF ( istat /= 0 ) THEN
            CALL die ( 'Failed to allocate cubes array.', 'allocate_3d_cube_array' )
        ENDIF

        DO iw = 0, box_grid(3) - 1
            DO iv = 0, box_grid(2) - 1
                DO iu = 0, box_grid(1) - 1
                    CALL allocate_array ( cubes(iu, iv, iw)%atom_number, cube_size(iu, iv, iw) )
                    CALL allocate_array ( cubes(iu, iv, iw)%xyz,         cube_size(iu, iv, iw) )
                    CALL allocate_array ( cubes(iu, iv, iw)%atom_symop,  cube_size(iu, iv, iw) )
                ENDDO
            ENDDO
        ENDDO

        CALL messag ( TRIM ( int_to_c (box_grid(1)*box_grid(2)*box_grid(3)))//' cubes have been allocated.',&
                     'allocate_3d_cube_array' )

    END SUBROUTINE allocate_3d_cube_array

    SUBROUTINE deallocate_3d_cube_array ( cubes )
        TYPE(cube), DIMENSION(:,:,:), ALLOCATABLE :: cubes
!       Local variables:
        INTEGER                                   :: istat
        INTEGER                                   :: number_of_cubes
!       Counters:
        INTEGER                                   :: iu
        INTEGER                                   :: iv
        INTEGER                                   :: iw

        number_of_cubes = SIZE(cubes, DIM=3) * SIZE(cubes, DIM=2) * SIZE(cubes, DIM=1)

        CALL messag ( 'Deallocating array of '//TRIM(int_to_c(number_of_cubes))//' cubes', 'deallocate_3d_cube_array' )

        DO iw = 0, SIZE(cubes, DIM=3) - 1
            DO iv = 0, SIZE(cubes, DIM=2) - 1
                DO iu = 0, SIZE(cubes, DIM=1) - 1
                    CALL deallocate_array ( cubes(iu, iv, iw)%atom_number )
                    CALL deallocate_array ( cubes(iu, iv, iw)%xyz         )
                    CALL deallocate_array ( cubes(iu, iv, iw)%atom_symop  )
                ENDDO
            ENDDO
        ENDDO

        DEALLOCATE ( cubes, STAT=istat )

        IF ( istat /= 0 ) THEN
            CALL die ( 'Failed to deallocate cubes array.', 'deallocate_3d_cube_array' )
        ENDIF
        CALL messag('Done. ','deallocate_3d_cube_array')
    END SUBROUTINE deallocate_3d_cube_array

    SUBROUTINE make_box ( bl, bu, mtz_1, pdb_2, r_crit, nout )
!
!       Purpose: 
!       =======
!       makes artifical rectangular cell whose dimension is by r_crit bigger then native one.
!
        REAL(KIND=wp), DIMENSION(:), INTENT(OUT) :: bl
        REAL(KIND=wp), DIMENSION(:), INTENT(OUT) :: bu
        TYPE(mtz),                   INTENT(IN)  :: mtz_1
        TYPE(pdb),                   INTENT(IN)  :: pdb_2
        REAL(KIND=wp),               INTENT(IN)  :: r_crit
        INTEGER,                     INTENT(IN)  :: nout
!       Local variables:
        REAL(KIND=wp), DIMENSION(3)              :: mod_bu_bl
        INTEGER,       DIMENSION(3)              :: box_grid
        REAL(KIND=wp), DIMENSION(3)              :: xyz_ort
!       Counters
        INTEGER                                  :: iat
        INTEGER                                  :: k
!
        bl = 1.E20_wp
        bu = 0.0_wp
        DO iat = 1, SIZE ( pdb_2 )

!           PDB_2 must contain fractional coordinates:
            xyz_ort = mtz_1%ORT * pdb_2%xyz(iat)
            DO k = 1, 3
                bl(k) = MIN ( bl(k), xyz_ort(k) )
                bu(k) = MAX ( bu(k), xyz_ort(k) )
            ENDDO
        ENDDO

        IF ( nout > 0 ) THEN
            WRITE(nout, "(' MAKE_BOX> ', 'r_crit=', 3F5.2, ' A')") r_crit
            WRITE(nout, "(' MAKE_BOX> ', 'lower bound=', 3F8.2)")  bl
            WRITE(nout, "(' MAKE_BOX> ', 'upper bound=', 3F8.2)")  bu
            WRITE(nout, "(' MAKE_BOX> ', 'difference  =', 3F8.2)") bu - bl
        ENDIF

        mod_bu_bl = MOD ( bu - bl, r_crit )
        bu = bu - 0.5_wp * mod_bu_bl
        bl = bl + 0.5_wp * mod_bu_bl
        IF ( nout > 0 ) THEN
            WRITE(nout, "(' MAKE_BOX> ', 'modified lower bound=', 3F8.2)")  bl
            WRITE(nout, "(' MAKE_BOX> ', 'modified upper bound=', 3F8.2)")  bu
            WRITE(nout, "(' MAKE_BOX> ', 'modified difference  =', 3F8.2)") bu - bl
            WRITE(nout, "(' MAKE_BOX> ', 'mod difference now =', 3F8.2)")   MOD ( bu - bl + 0.000001_wp, r_crit )
        ENDIF

        bu = bu + 2.0_wp * r_crit
        bl = bl - 2.0_wp * r_crit
        IF ( nout > 0 ) THEN
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'final lower bound=', 3F8.2)")  bl
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'final upper bound=', 3F8.2)")  bu
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'final difference= ', 3F8.2)")  bu - bl
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'mod difference now (should be zero vector)=', 3F5.2)") &
            MOD(bu - bl + 0.000001_wp, r_crit)
        ENDIF

        box_grid = NINT ( ( bu - bl ) / r_crit )
        IF ( nout > 0 ) THEN 
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'xyz box grid=', 3I5)") box_grid 
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'total number of ',F5.2,' Angstrom cubes= ',A)") r_crit, &
            TRIM ( int_to_c ( box_grid(1) * box_grid(2) * box_grid(3) ) )
            WRITE(nout, "(' MAKE_BOX_RESULT> ', 'likely speed up=', F7.1, ' times')") &
            box_grid(1) * box_grid(2) * box_grid(3) / 27.0_wp
        ENDIF

        IF ( debug > 100 ) THEN
            WRITE(*,"(' MAKE_BOX_RESULT> ', 'likely speed up=', F7.1, ' times')") &
            (mtz_1%volume / mtz_1%sp_group%number_of_symops) / (3.0_wp * r_crit) ** 3
        ENDIF

    END SUBROUTINE make_box
  
    SUBROUTINE generate_box_bound_sym_atoms ( my_cube, bl, bu, mtz_1, pdb_1, r_crit, max_cube_size, nout )
!       
!       Purpose:
!       =======
!       Generate symmetry related atoms falling inside rectangular box
!       limited by lower bound bl(1:3) and upper bound bu(1:3)
!
        TYPE(cube),    DIMENSION(:,:,:), ALLOCATABLE :: my_cube
        REAL(KIND=wp), DIMENSION(:),     INTENT(IN)  :: bl
        REAL(KIND=wp), DIMENSION(:),     INTENT(IN)  :: bu
        TYPE(mtz),                       INTENT(IN)  :: mtz_1
        TYPE(pdb),                       INTENT(IN)  :: pdb_1
        REAL(KIND=wp),                   INTENT(IN)  :: r_crit
        INTEGER,                         INTENT(OUT) :: max_cube_size
        INTEGER,                         INTENT(IN)  :: nout
!       Local variables:
        INTEGER,       DIMENSION(:,:,:), ALLOCATABLE :: cube_size                
        INTEGER,       DIMENSION(3)                  :: ind
        INTEGER,       PARAMETER                     :: trmax = 2 
        TYPE(vector)                                 :: translation
        TYPE(vector)                                 :: xyz
        TYPE(vector)                                 :: xyz_sym
        REAL(KIND=wp), DIMENSION(:,:),   ALLOCATABLE :: xyz_temp
        REAL(KIND=wp), DIMENSION(3)                  :: xyz_ort
        REAL(KIND=wp), DIMENSION(3)                  :: v
        REAL(KIND=wp)                                :: box_volume
        INTEGER,       DIMENSION(3)                  :: box_grid
        INTEGER                                      :: istat
        INTEGER                                      :: n
!       Counters:
        INTEGER                                      :: atom_counter
        INTEGER                                      :: i 
        INTEGER                                      :: iu
        INTEGER                                      :: iv
        INTEGER                                      :: iw
        INTEGER                                      :: k
        INTEGER                                      :: m
        INTEGER                                      :: run

!       Checkz:
        IF ( .NOT. ALLOCATED ( pdb_1 ) ) THEN
            CALL die ( 'Programming error. pdb_1 has not been initialized.', 'generate_box_bound_sym_atoms' )
        ELSE IF ( SIZE ( bl ) /= SIZE ( bu ) ) THEN
            CALL die ( 'Programming error.  bl and bu arrays have different size.', 'generate_box_bound_sym_atoms' )
        ENDIF

        n = SIZE ( pdb_1 )

        IF ( nout > 0 ) THEN
            WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'lower bound=', 3F8.2)")  bl
            WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'upper bound=', 3F8.2)")  bu
        ENDIF

        box_volume = 1.0_wp
        DO k = 1, 3
            box_volume = box_volume * ( bu(k) - bl(k) )
        ENDDO

!       Estimate number of atoms which would fit in the rectangular box:
        IF ( nout > 0 ) THEN
            WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'rectangular box volume=',F12.1)") box_volume
            WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'expected number of generated atoms= ',A)") &
            TRIM ( int_to_c ( NINT ( ( n * mtz_1%sp_group%number_of_symops ) * box_volume / mtz_1%volume ) ) )
        ENDIF

        box_grid = NINT ( (bu - bl) / r_crit )

        IF ( nout > 0 ) THEN
            WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'box grid=', 3I5)" ) box_grid
        ENDIF

        CALL allocate_array ( cube_size, box_grid, 'cube_size' )

        outer_loop:DO run = 0, 1

        cube_size    = 0
        atom_counter = 0

            atoms:DO i = 1, n

                IF ( nout > 0 ) THEN
                    IF ( MOD ( i, 1000 ) == 1) THEN
                        WRITE(nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ','atom number= ', A)") &
                        TRIM ( int_to_c ( i ) )
                    ENDIF
                ENDIF

                xyz = pdb_1%xyz(i)  ! assumes that pdb_1 contains fractional coordinates as needed
                symmetry:DO m = 1, mtz_1%sp_group%number_of_symops
                    xyz_sym = mtz_1%sp_group%SYM(m) * xyz
                    along_z:DO iw = -trmax, trmax
                        v(3) = REAL(iw, KIND=wp)
                        along_y:DO iv = -trmax, trmax
                            v(2) = REAL(iv, KIND=wp)
                            along_x:DO iu = -trmax, trmax
                                v(1) = REAL(iu, KIND=wp)
                                translation = v
                                xyz_ort = mtz_1%ORT * ( xyz_sym + translation )
                                DO k = 1, 3
                                    IF ( xyz_ort(k) > bu(k) ) CYCLE along_x
                                    IF ( xyz_ort(k) < bl(k) ) CYCLE along_x
                                ENDDO

!                               If we are here atom successfully passed bound test above:
                                atom_counter = atom_counter + 1

!                               Cube index by any means must not exceed box_grid in any direction:
                                ind = INT ( (xyz_ort - bl) / r_crit )
                                DO k = 1,3
                                    IF (ind(k) == box_grid(k) ) THEN
                                        ind(k) = box_grid(k) - 1
                                    ENDIF
                                ENDDO

!                               Here we can also estimate number of atoms in each cube:
                                cube_size(ind(1), ind(2), ind(3)) = cube_size(ind(1), ind(2), ind(3)) + 1

!                               If arrays already allocated arrays then make assignments:
                                IF ( run == 1) THEN
                                    my_cube(ind(1), ind(2), ind(3))%atom_number(cube_size(ind(1), ind(2), ind(3))) = i
                                    my_cube(ind(1), ind(2), ind(3))%xyz        (cube_size(ind(1), ind(2), ind(3))) = xyz_ort
!                                   Probably need to orthogonalize this:
                                    my_cube(ind(1), ind(2), ind(3))%atom_symop (cube_size(ind(1), ind(2), ind(3))) = &
                                    mtz_1%sp_group%SYM(m) + translation
                                ENDIF
                            ENDDO along_x
                        ENDDO along_y
                    ENDDO along_z
                ENDDO symmetry
            
            ENDDO atoms

!           Ready to allocate cube array:
            myrun:IF ( run == 0 ) THEN
                CALL allocate_array ( my_cube, box_grid, cube_size )
                IF ( nout > 0 ) THEN
                    WRITE(nout,"(' GENERATE_BOX_BOUND_SYM_ATOMS> ',&
                                &' *** Number of generated box symmetry related atoms= ', A)") &
                                TRIM ( int_to_c ( atom_counter ) )
                ENDIF

!               Check whether we need extensive printing:
                IF ( nout > 0 ) THEN

                    DO iw = 0, box_grid(3) - 1
                        DO iv = 0, box_grid(2) - 1
                            IF ( box_grid(1) - 1 > 2*5 ) THEN
                                WRITE(*,"(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 2I5,' cube sizes=', 5I6, '   ... ', 5I6)") &
                                iw, iv, (cube_size(iu, iv, iw), iu=              0, MIN(4, box_grid(1)-1)), &
                                        (cube_size(iu, iv, iw), iu=box_grid(1)-1-4,        box_grid(1)-1 )
                            ELSE
                                WRITE(*,"(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 2I5,' cube sizes=', 5I6, '   ... ', 5I6)") &
                                iw, iv, (cube_size(iu, iv, iw), iu=              0, MIN(4, box_grid(1)-1))
                            ENDIF
                        ENDDO
                    ENDDO

                ENDIF
            ELSE

                printing:IF ( nout > 0 ) THEN
                    CALL messag('Printing 0,0,0 cube...', 'generate_box_bound_sym_atoms')
                    WRITE(nout,"(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'size(0,0,0)=',I6,                     &
                               & ' atom numbers=', 3(12X,I10),  ' ...')")                                  &
                    SIZE ( my_cube(0,0,0)%atom_number ),                                                   &
                    (my_cube(0,0,0)%atom_number(i), i = 1, MIN ( 3, SIZE ( my_cube(0,0,0)%atom_number ) ))
!                   XYZ_TEMP is needed for safe printing:
                    IF ( SIZE ( my_cube(0,0,0)%xyz ) > 0 ) THEN
                        CALL allocate_array ( xyz_temp, 3, MIN( 3, SIZE( my_cube(0,0,0)%xyz ) ) )
                        DO i = 1, SIZE ( xyz_temp, DIM=2 )
                            xyz_temp(1:3,i) = my_cube(0,0,0)%xyz(i)
                        ENDDO
                        WRITE (nout, "(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'size(0,0,0)=',I6,&
                                 & ' xyz=         ', 3(1X,3F7.2), ' ...')")&
                        SIZE(my_cube(0,0,0)%xyz), (xyz_temp(1:3,i), i=1,SIZE ( xyz_temp, DIM=2 ))
                    
                        CALL deallocate_array ( xyz_temp )
                    ENDIF
                ENDIF printing

            ENDIF myrun

        ENDDO outer_loop

!       Calculate max_cube_size before deallocation:
        max_cube_size = MAXVAL ( cube_size )
        WRITE(*,"(' GENERATE_BOX_BOUND_SYM_ATOMS> ', 'Max. number of atoms in cube=', I8)") max_cube_size

        CALL deallocate_array ( cube_size, 'cube_size' )
        CALL messag ( 'Done.', 'generate_box_bound_sym_atoms' )

    END SUBROUTINE generate_box_bound_sym_atoms

    SUBROUTINE fast_cube_dist ( adj_list, my_cube, bl, bu, mtz_1, pdb_2, r_crit, max_cube_size )
!
!       Purpose:
!       =======
!       a) Prepare matrix in sparse co-ordinate format. r_short must be absent
!          when calling routine for this purpose.
!       b) If there is a need for list of non-bonded contacts then we must produce a list of pairs of
!          atoms together with their symmetry operators. Perhaps we need to modify this routine for
!          this particular purpose. This has not been done yet.
!
!       Note:
!       ====
!       On input PDB_2 MUST contain fractional coordinates.
!
!       Date:             Programmer:          Description of changes:
!       ====              ==========           ======================
!       Nov 2004          B.Strokopytov        Original code
!       Dec 2005          B.Strokopytov        ALLOCATABLE arrays introduced.
!       Dec 2005          B.Strokopytov        Swapped to lower triangle of H_val matrix
!                                              SEE "BUG CORRECTED" below.
!
        TYPE(adjacency_list), DIMENSION(:),     ALLOCATABLE, INTENT(INOUT)         :: adj_list
        TYPE(cube),           DIMENSION(:,:,:), ALLOCATABLE                        :: my_cube
        REAL(KIND=wp),        DIMENSION(:),                  INTENT(IN)            :: bl
        REAL(KIND=wp),        DIMENSION(:),                  INTENT(IN)            :: bu
        TYPE(mtz),                                           INTENT(IN)            :: mtz_1
        TYPE(pdb),                                           INTENT(INOUT)         :: pdb_2
        REAL(KIND=wp),                                       INTENT(IN)            :: r_crit
        INTEGER,                                             INTENT(IN)            :: max_cube_size
!       Local variables:
        REAL(KIND=wp),        DIMENSION(3)                                         :: vec_diff
        REAL(KIND=wp),        DIMENSION(3)                                         :: xyz_ort
        REAL(KIND=wp),        DIMENSION(3)                                         :: xyz
        INTEGER,              DIMENSION(3)                                         :: ind
        INTEGER,              DIMENSION(3)                                         :: box_grid
        REAL(KIND=wp)                                                              :: r_crit2
        REAL(KIND=wp)                                                              :: r2
        INTEGER                                                                    :: modcon
!       Duplicates elimination:
        INTEGER                                                                    :: maxjat
        INTEGER                                                                    :: duplicates
        INTEGER,              DIMENSION(9*max_cube_size)                           :: jat_number
        REAL(KIND=wp),        DIMENSION(9*max_cube_size)                           :: dist
!       Test:
        LOGICAL                                                                    :: printing
!       Counters:
        INTEGER                                                                    :: contact
        INTEGER                                                                    :: iat
        INTEGER                                                                    :: iu
        INTEGER                                                                    :: iv
        INTEGER                                                                    :: iw
        INTEGER                                                                    :: j
        INTEGER                                                                    :: jat
        INTEGER                                                                    :: jat_atom_contacts
        INTEGER                                                                    :: total_duplicates

!       Checkz:
        IF ( r_crit <= 0.0_wp ) THEN
            CALL die ('Programming error. r_crit should be positive number.', 'fast_cube_dist' )
        ENDIF

!       For diagonal elements only we may use r_crit = 0.1 A but then this routine is not needed: 
        CALL messag ( ' ', 'fast_cube_dist' )
        r_crit2 = r_crit ** 2

!       Recalculate box grid first:
        box_grid = NINT ( (bu - bl) / r_crit )
        WRITE(*, "(' FAST_CUBE_DIST> ', 'box_grid=', 3I4)") box_grid

        WRITE(*, "(' FAST_CUBE_DIST> ', 'r_crit= ', F5.2, ' r_crit**2=', F7.2)")&
        r_crit, r_crit2

        IF ( ANY ( pdb_2%occ < 0.0_wp) ) THEN
            CALL die ( 'Programming error. pdb_2 contains negative occupancies.', 'fast_cube_dist' )
        ENDIF

        modcon = 1000
        maxjat = 9 * max_cube_size
        CALL allocate_array ( adj_list, SIZE ( pdb_2 ) )

!       Initialize counters:
        contact          = 0
        total_duplicates = 0

!       Loop over all atoms:
        atoms_loop:DO iat = 1, SIZE ( pdb_2 )

            IF ( pdb_2%occ(iat) <= 0.0_wp ) CYCLE atoms_loop

!           pdb_2 co-ordinates MUST be in fractions:
            xyz_ort  = mtz_1%ORT * pdb_2%xyz(iat)       
            ind      = INT ( (xyz_ort - bl) / r_crit )
            jat_atom_contacts = 0 ! need to know number of contacts for iat atom number

            DO iw = ind(3) - 1, ind(3) + 1
                IF ( iw > box_grid(3) - 1 ) THEN
                    WRITE(*,*) iw, box_grid(3)
                    WRITE(*,*) ' iat=',iat, ' xyz_ort=', xyz_ort
                    CALL die('Programming error. Box grid bound violated along_z.', 'fast_cube_dist')
                ENDIF

                DO iv = ind(2) - 1, ind(2) + 1

                    IF ( iv > box_grid(2) - 1 ) THEN
                        WRITE(*,*) iv, box_grid(2)
                        WRITE(*,*) ' iat=',iat, ' xyz_ort=', xyz_ort
                        WRITE(*,*) ' ind =', ind
                        CALL die('Programming error. Box grid bound violated along_y.', 'fast_cube_dist')
                    ENDIF

                    along_x:DO iu = ind(1) - 1, ind(1) + 1
                        IF ( iu > box_grid(1) - 1 ) THEN
                            WRITE(*,*) iu, box_grid(1)
                            CALL die('Programming error. Box grid bound violated along_x.', 'fast_cube_dist')
                        ENDIF

!                       Loop along all atoms in this cube:
                        cube_loop:DO j = 1, SIZE ( my_cube(iu,iv,iw)%atom_number )
                            jat = my_cube(iu,iv,iw)%atom_number(j)

!                           Atom numbers are always sorted in increasing order inside cubes
!                           There is no need to look for atoms with numbers less than iat:

!                           Probably removed because we want an adjacency list:
!BVS                        IF ( jat > iat ) CYCLE along_x


!                           No need to form an atom pair (zero or negative occupnacy):
                            IF ( pdb_2%occ(jat) <= 0.0_wp ) CYCLE cube_loop

!                           Coordinates in Angstrom units in my_cube:
                            xyz = my_cube(iu,iv,iw)%xyz(j)
                            vec_diff = xyz - xyz_ort

!                           Some speed up:
                            IF ( ANY ( ABS ( vec_diff ) > r_crit ) ) CYCLE cube_loop

!                           Calculate squared distance between atoms:
                            r2 = DOT_PRODUCT ( vec_diff, vec_diff )

                            IF ( r2 <= r_crit2 ) THEN
!
                                    
                                jat_atom_contacts = jat_atom_contacts + 1   ! counting contacts for iat atom
                                IF ( jat_atom_contacts > maxjat ) THEN
                                    WRITE(*,*) ' iat=', iat, ' jat_atom_contacts=', jat_atom_contacts, ' > maxjat=', maxjat
                                    CALL die ( 'Recompile. Too many contacts for current atom.',&
                                               'fast_cube_dist' )
                                ENDIF

!                               Prepare to check for duplicates: 
                                jat_number(jat_atom_contacts) = jat
                                dist(jat_atom_contacts) = SQRT ( r2 )

!                               Generally we need to figure out how many short contacts we have:
                                contact  = contact  + 1
                                IF ( MOD ( contact, modcon ) == 1 ) THEN
                                    WRITE(*,"(' FAST_CUBE_DIST> ', 'distance ', A,  ' ... ', A,' =', F5.2, &
                                   &T60,'short distances found=', I12)") &
                                    TRIM(int_to_c(iat)), TRIM(int_to_c(jat)), SQRT(r2), contact
                                ENDIF
                            ENDIF
                        ENDDO cube_loop
                    ENDDO along_x
                ENDDO
            ENDDO

!           Eliminate duplicates (they may arise due to, e.g., proper 2-fold symmetry):
            printing =.FALSE.
!            IF ( iat == 1041 ) printing = .TRUE.
            CALL eliminate_duplicates ( maxjat, jat_number, dist, jat_atom_contacts, duplicates, printing )
            total_duplicates = total_duplicates + duplicates
            contact = contact - duplicates  ! necessary for proper matrix allocation and usage

!           Save numbers as adjacency list for each atom:
            IF ( jat_atom_contacts > 0 ) THEN

!               Allocation of ADJ_LIST component will be done during copy operation:
!               BUG CORRECTED FEB 2008 BVS: "-duplicates" was missing:
!               We need to keep only unique elements for matrix but for NB contacts we need
!               quite different approach (interaction of atom with all symmetry mates):
                adj_list(iat) = jat_number(1:jat_atom_contacts - duplicates)

                IF ( debug > 40 ) THEN
                    WRITE(*,"(' atom ', I5, ' degree ', I2, ' list', 20I4)") &
                    iat, SIZE(adj_list(iat)), adj_list(iat)%v
                ENDIF
            ELSE
                IF ( debug > 40 ) THEN
                    WRITE(*,"(' atom ', I5, ' degree ', I2, ' list', 20I4)") &
                    iat, SIZE( adj_list(iat) )
                ENDIF
            ENDIF

        ENDDO atoms_loop

!   Check number of atoms with positive occupancy:
    jat = COUNT ( pdb_2%occ > 0.0_wp )
    WRITE ( *, "(' FAST_CUBE_DIST> ', 'total number of atomic pairs         =   ', A)") &
    TRIM ( int_to_c ( contact ) )
    WRITE ( *, "(' FAST_CUBE_DIST> ', 'total number of duplicates eliminated=   ', A)") &
    TRIM ( int_to_c ( total_duplicates ) )
    WRITE ( *, "(' FAST_CUBE_DIST> ', 'number of atoms with positive occupancy= ', A)") &
    TRIM ( int_to_c ( jat ) )

    END SUBROUTINE fast_cube_dist

    SUBROUTINE eliminate_duplicates ( maxjat, arr, dist, n, duplicates, printing )
!
!       Purpose:
!       =======
!       Eliminates duplicates in an integer array
!       Length of array arr "n" is modified on output
!       as well as its contents if duplicates are found
        INTEGER,                          INTENT(IN)    :: maxjat
        INTEGER,       DIMENSION(maxjat), INTENT(INOUT) :: arr
        REAL(KIND=wp), DIMENSION(maxjat), INTENT(INOUT) :: dist
        INTEGER,                          INTENT(INOUT) :: n
        INTEGER,                          INTENT(OUT)   :: duplicates
        LOGICAL,                          INTENT(IN)    :: printing
!       Automatic arrays:
        INTEGER,       DIMENSION(n)                     :: temp
        INTEGER,       DIMENSION(n)                     :: tag
        REAL(KIND=wp), DIMENSION(n)                     :: criteria
!       Counters:
        INTEGER                                         :: i
        INTEGER                                         :: positive

!       Initialize duplicates counter:
        duplicates = 0

!       Just one element in array - no need to bother:
        IF ( n <= 1 ) RETURN

        IF ( maxjat < n ) THEN
            WRITE(*,*) maxjat, n
            CALL die ( 'Programming error. Oops... Array ARR has insufficient size.', 'eliminate_duplicates' )
        ENDIF

        criteria = arr(1:n)
        tag = (/ (i,i=1,n) /)
        CALL sortag ( criteria, n, tag )

        IF ( printing ) THEN
            WRITE(*,*) ' criteria=', criteria
            WRITE(*,*) ' tag=' ,tag
        ENDIF

        DO i = 1, n
            temp(i)     = arr(tag(i))
            criteria(i) = dist(tag(i))
        ENDDO

        IF ( printing ) THEN
            WRITE(*,*) ' temp=', temp
            WRITE(*,*) ' criteria=', criteria
              
        ENDIF

!       At this point we have sorted arr in temp and sorted distances in criteria:
        DO i = 2, n
            IF ( temp(i-1) == temp(i) ) THEN

!               Important to change previous, i-1 th element:
                temp(i-1) = -temp(i-1) 

!               Save minimal distance when eliminating duplicates (not absolutely necessary):
                criteria(i) = MIN ( criteria(i-1), criteria(i) ) 
                duplicates = duplicates + 1
            ENDIF
        ENDDO

        positive = COUNT ( temp > 0 )
        IF ( printing ) THEN
            WRITE(*,*) ' temp=', temp
            WRITE(*,*) ' positive=', positive
        ENDIF

!       Use FORTRAN 90 PACK function:    
        dist(1:positive) = PACK ( criteria(1:n), temp > 0 )
        arr(1:positive)  = PACK ( temp,          temp > 0 )
        IF ( printing ) THEN
            WRITE(*,*) ' dist=', dist(1:positive)
            WRITE(*,*) ' arr=',  arr(1:positive)
            STOP 'STOP in eliminate duplicates'
        ENDIF
!       Array temp, criteria & tag should be deallocated automatically in f95:
    END SUBROUTINE eliminate_duplicates

    FUNCTION number_of_bad_contacts(mtz_1, convergence_radius, r_short)
        INTEGER                   :: number_of_bad_contacts
        TYPE(mtz),     INTENT(IN) :: mtz_1
        REAL(KIND=wp), INTENT(IN) :: convergence_radius
        REAL(KIND=wp), INTENT(IN) :: r_short              
!       Local variables:
        INTEGER                   :: expected_number_of_atoms_in_ASU
        REAL(KIND=wp)             :: ASU_volume
        REAL(KIND=wp)             :: vol_ratio

!       Estimate total number of contacts:
        ASU_volume = mtz_1%volume / mtz_1%sp_group%number_of_symops
        expected_number_of_atoms_in_ASU = ASU_volume/sphere_volume(convergence_radius)
        vol_ratio = sphere_volume(r_short) / ASU_volume
        vol_ratio = vol_ratio + vol_ratio ** 2 + vol_ratio ** 3 ! enough is enough 
        number_of_bad_contacts = 0.5 * vol_ratio * REAL(expected_number_of_atoms_in_ASU, KIND=wp) ** 2   
    END FUNCTION number_of_bad_contacts

    SUBROUTINE check_fast_distance ( adj_list, mtz_1, pdb_2, r_crit )
        TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: adj_list
        TYPE(mtz),                                       INTENT(IN)    :: mtz_1
        TYPE(pdb),                                       INTENT(IN)    :: pdb_2
        REAL(KIND=wp),                                   INTENT(IN)    :: r_crit
!       Local variables:
        REAL(KIND=wp)                                                  :: r_crit2
        REAL(KIND=wp)                                                  :: r2
        INTEGER                                                        :: maxjat = 10000
        INTEGER,              DIMENSION(:), ALLOCATABLE                :: jat_number   ! for duplicate elimination
        INTEGER                                                        :: iat
        INTEGER                                                        :: jat
        INTEGER                                                        :: jat_atom_contacts

        CALL allocate_array ( jat_number, maxjat )
        CALL allocate_array ( adj_list, SIZE ( pdb_2 ) )

        r_crit2 = r_crit ** 2

        DO iat = 1, SIZE ( pdb_2 )
            jat_atom_contacts = 0
            DO jat = 1, SIZE ( pdb_2 )
                IF ( iat /= jat ) THEN 
                    r2 = min_dist ( mtz_1, pdb_2%xyz(iat), pdb_2%xyz(jat) )
                    IF ( r2  <= r_crit2 ) THEN
                        jat_atom_contacts = jat_atom_contacts + 1
                        jat_number(jat_atom_contacts) = jat 
                    ENDIF
                ENDIF
            ENDDO

!           Here we may eliminate duplicates:
            IF ( jat_atom_contacts > 0 ) adj_list(iat) = jat_number(1:jat_atom_contacts)
        ENDDO

    END SUBROUTINE check_fast_distance

    SUBROUTINE pairs(pair, mtz_1, pdb_2, r_crit, nout)
        INTEGER,              DIMENSION(:,:),   ALLOCATABLE, INTENT(INOUT) :: pair
        TYPE(mtz),                                           INTENT(IN)    :: mtz_1
        TYPE(pdb),                                           INTENT(INOUT) :: pdb_2
        REAL(KIND=wp),                                       INTENT(IN)    :: r_crit
        INTEGER,                                             INTENT(IN)    :: nout
!       Local variables:
        REAL(KIND=wp),        DIMENSION(3)                                 :: bl
        REAL(KIND=wp),        DIMENSION(3)                                 :: bu
        INTEGER                                                            :: max_cube_size
!       Cubing algorithm params:
        TYPE(cube),           DIMENSION(:,:,:), ALLOCATABLE                :: my_cube
!       Array of adjacency lists:
        TYPE(adjacency_list), DIMENSION(:),     ALLOCATABLE                :: adj_list
        REAL(KIND=sp)                                                      :: average_number_of_contacts
!       Counters:
        INTEGER                                                            :: i
        INTEGER                                                            :: j
        INTEGER                                                            :: modcon
        INTEGER                                                            :: number_of_entries
        CHARACTER(LEN=32)                                                  :: srname='pairs'

!       Convert PDB_2 to fractions:
        pdb_2 = pdb_2%DEORT * pdb_2
        CALL make_box ( bl, bu, mtz_1, pdb_2, r_crit, nout )
        CALL generate_box_bound_sym_atoms ( my_cube, bl, bu, mtz_1, pdb_2, r_crit, max_cube_size, nout )
        CALL fast_cube_dist ( adj_list, my_cube, bl, bu, mtz_1, pdb_2, r_crit, max_cube_size )

!       Free memory:
        CALL deallocate_array ( my_cube )

!       Convert back to ANGSTROMS:
        pdb_2 = pdb_2%ORT * pdb_2

!       Statistics:
        number_of_entries = 0
        DO i = 1, SIZE ( adj_list )
            number_of_entries = number_of_entries + SIZE ( adj_list(i) )
        ENDDO

        IF ( nout > 0 ) THEN
            WRITE(nout, "(' PAIRS> ', 'number of entries in adjacency list=', I10)") &
            number_of_entries
            WRITE(nout, "(' PAIRS> ', 'average number of contacts per atom=', F7.2)" ) &
            REAL ( number_of_entries ) / SIZE ( adj_list )
        ENDIF

!       Allocate PAIR array:
        CALL allocate_array ( pair, 2, number_of_entries )

        modcon = MAX ( number_of_entries / 10, 1)
        number_of_entries = 0

!       Debugging:
        IF ( debug > 100 ) THEN
            WRITE(*,*) ' Printing adjacency list:'
        ENDIF

        DO i = 1, SIZE ( adj_list )

!           Debugging:
            IF ( debug > 100 ) THEN
                WRITE(*,"(I6, 10X, 100I6)") i, adj_list(i)%v
            ENDIF

            DO j = 1, SIZE ( adj_list(i) )

                number_of_entries = number_of_entries + 1
                pair(1,number_of_entries) = i
                pair(2,number_of_entries) = adj_list(i)%v(j)

!               Monitor atom pairs:
                IF ( nout > 0 ) THEN
                    IF ( MOD ( number_of_entries, modcon ) == 0 )  THEN
                        WRITE(*,"(' PAIRS> ', 'entry number=', I9, ' atomic pair=', 2I6)") &
                        number_of_entries, pair(1:2,number_of_entries)
                    ENDIF
                ENDIF

            ENDDO

        ENDDO

        CALL messag('Finished adjacency list processing...', srname)

!       No need to keep it anymore:
        CALL deallocate_array ( adj_list )

    END SUBROUTINE pairs

    SUBROUTINE fast_neighbours_list(adj_list, mtz_1, pdb_2, r_crit, nout)
        TYPE(adjacency_list), DIMENSION(:),     ALLOCATABLE, INTENT(INOUT) :: adj_list
        TYPE(mtz),                                           INTENT(IN)    :: mtz_1
        TYPE(pdb),                                           INTENT(INOUT) :: pdb_2
        REAL(KIND=wp),                                       INTENT(IN)    :: r_crit
        INTEGER,                                             INTENT(IN)    :: nout
!       Local variables:
        REAL(KIND=wp),        DIMENSION(3)                                 :: bl
        REAL(KIND=wp),        DIMENSION(3)                                 :: bu
        INTEGER                                                            :: max_cube_size
!       Cubing algorithm params:
        TYPE(cube),           DIMENSION(:,:,:), ALLOCATABLE                :: my_cube
!       Array of adjacency lists:
        REAL(KIND=sp)                                                      :: average_number_of_contacts
        INTEGER,              DIMENSION(:), ALLOCATABLE                    :: temp
!       Counters:
        INTEGER                                                            :: i
        INTEGER                                                            :: j
!        INTEGER                                                            :: modcon
        INTEGER                                                            :: number_of_entries
        CHARACTER(LEN=32)                                                  :: srname='pairs'

!       Convert PDB_2 to fractions:
        pdb_2 = pdb_2%DEORT * pdb_2
        CALL make_box ( bl, bu, mtz_1, pdb_2, r_crit, nout )
        CALL generate_box_bound_sym_atoms ( my_cube, bl, bu, mtz_1, pdb_2, r_crit, max_cube_size, nout )
        CALL fast_cube_dist ( adj_list, my_cube, bl, bu, mtz_1, pdb_2, r_crit, max_cube_size )

!       Free memory:
        CALL deallocate_array ( my_cube )

!       Convert back to ANGSTROMS:
        pdb_2 = pdb_2%ORT * pdb_2

!       Statistics:
        number_of_entries = 0
        DO i = 1, SIZE ( adj_list )
            number_of_entries = number_of_entries + SIZE ( adj_list(i) )
        ENDDO

        IF ( nout > 0 ) THEN
            WRITE(nout, "(' PAIRS> ', 'number of entries in adjacency list=', I10)") &
            number_of_entries
            WRITE(nout, "(' PAIRS> ', 'average number of contacts per atom=', F7.2)" ) &
            REAL ( number_of_entries ) / SIZE ( adj_list )
        ENDIF

!       Debugging:
!        IF ( debug > 100 ) THEN
            WRITE(*,*) ' Printing adjacency list:'
!        ENDIF

        DO i = 1, SIZE ( adj_list )

!           Remove atom with index I from the list (self-interaction):
            CALL remove_node(adj_list(i), i)

!           Debugging:
!           IF ( debug > 100 ) THEN
                WRITE(*,"(I6, 10X, 100I6)") i, adj_list(i)%v
!           ENDIF

        ENDDO

        CALL messag('Finished adjacency list processing...', srname)

    END SUBROUTINE fast_neighbours_list 

END MODULE fast_distance
