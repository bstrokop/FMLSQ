MODULE logical_map_manip
USE constants
USE fail
USE fft_util
USE mtz_io
USE util
USE vectors
IMPLICIT NONE

TYPE :: map_logical
    INTEGER                                     :: nu
    INTEGER                                     :: nv
    INTEGER                                     :: nw
    TYPE(vector_int)                            :: map_size
    REAL(KIND=wp)                               :: grid_factor
    TYPE(space_group)                           :: sp_group
    TYPE(matrix)                                :: ORT
    TYPE(matrix)                                :: DEORT
    TYPE(symop),      DIMENSION(:), ALLOCATABLE :: SYM_DEORT 
    TYPE(tensor)                                :: real_tensor
    TYPE(tensor)                                :: recip_tensor
    REAL(KIND=wp)                               :: radius
    CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: atom_name

!   map 3D array:
    LOGICAL(KIND=1), DIMENSION(:,:,:), ALLOCATABLE :: array 
    INTEGER                                        :: ntrue
END TYPE
PRIVATE
PUBLIC :: ALLOCATED
PUBLIC :: SIZE
PUBLIC :: allocate_map
PUBLIC :: deallocate_map
PUBLIC :: map_logical
INTERFACE ALLOCATED
    MODULE PROCEDURE allocated_logical_map
END INTERFACE

INTERFACE SIZE
    MODULE PROCEDURE logical_map_size
END INTERFACE

INTERFACE allocate_map
    MODULE PROCEDURE allocate_logical_map
END INTERFACE

INTERFACE deallocate_map
    MODULE PROCEDURE deallocate_logical_map
END INTERFACE

CONTAINS
    FUNCTION allocated_logical_map ( map_1 )
        LOGICAL                       :: allocated_logical_map
        TYPE(map_logical), INTENT(IN) :: map_1

        allocated_logical_map = ALLOCATED ( map_1%array )
    END FUNCTION allocated_logical_map

    FUNCTION logical_map_size ( map_1 )
        TYPE(vector_int)              :: logical_map_size
        TYPE(map_logical), INTENT(IN) :: map_1

        logical_map_size = UBOUND ( map_1%array ) - LBOUND ( map_1%array ) + 1
    END FUNCTION logical_map_size

    SUBROUTINE allocate_logical_map(map_1, mtz_2)
!
!       Purpose:
!       =======
!       Allocates logical byte map for overlap calculation
!
!       Date:        Programmer:         Description of changes:
!       ====         ==========          ======================
!       Dec 2003     B. Strokopytov      Original code
!
!
!
        TYPE(map_logical), INTENT(INOUT)        :: map_1
        TYPE(mtz),         INTENT(INOUT)        :: mtz_2
!        TYPE(cards),       INTENT(IN)           :: cards_1
!       Local variables:
        INTEGER,  DIMENSION(3)                  :: map_size
        INTEGER                                 :: m
        INTEGER                                 :: istat
        REAL(KIND=wp)                           :: resol_max
!        
!       Need to reset mtz_2%resol_max to use artificial resolution. Save resol_max:
        resol_max = mtz_2%resol_max
!
!       Assigned arbitrary grid_factor 2.5 seems as good compromise
!       between speed and accuracy:
!        map_1%grid_factor = 2.5_wp
!
!       Introducing artifical 3.0 A resolution:
!        mtz_2%resol_max = 3.0_wp
        map_1%map_size  = gridres ( mtz_2, map_1%grid_factor )
!
!       Restoring old value of resol_max:
        mtz_2%resol_max = resol_max
!
!       Just duplicate map_size:
        map_size  = map_1%map_size
        map_1%nu  = map_size(1)
        map_1%nv  = map_size(2)
        map_1%nw  = map_size(3)
!
        map_1%ORT          = mtz_2%ORT
        map_1%DEORT        = mtz_2%DEORT
        map_1%REAL_TENSOR  = mtz_2%REAL_TENSOR
        map_1%RECIP_TENSOR = mtz_2%RECIP_TENSOR
        map_1%sp_group     = mtz_2%sp_group

!       Setup modified symmetry operators:
        CALL allocate_array ( map_1%SYM_DEORT, map_1%sp_group%number_of_symops )
        DO m = 1, map_1%sp_group%number_of_symops 
            map_1%SYM_DEORT(m) = map_1%sp_group%SYM(m) * map_1%DEORT
        ENDDO

!       Setup uniform atomic radius.
!       2.5-2.6 ANGSTROMS seems to be close to optimal for complete models:
!        map_1%radius = cards_1%atom_radius
!        CALL allocate_array ( map_1%atom_name, SIZE ( cards_1%atom_name ) )
!        map_1%atom_name = cards_1%atom_name

!        WRITE(*, "(' ALLOCATE_LOGICAL_MAP> ', 'atom radius=', F6.1,' ANGSTROMS')" ) map_1%radius
!        WRITE(*, "(' ALLOCATE_LOGICAL_MAP> ', 'atom names= ', 100(1X,A))" ) map_1%atom_name
        map_1%ntrue = 0

!       Memory allocation, no additional sections needed:
        IF ( .NOT. ALLOCATED ( map_1%array ) ) THEN        
            ALLOCATE ( map_1%array( 0:map_1%nu-1, 0:map_1%nv-1, 0:map_1%nw-1 ), STAT=istat )
            IF ( istat == 0 ) THEN
                CALL messag( 'Logical map has been successfully allocated.', 'allocate_logical_map' )
            ELSE
                CALL die( 'Failed to allocate logical map.', 'allocate_logical_map' )
            ENDIF
        ELSE
            CALL deallocate_logical_map ( map_1 )
            ALLOCATE ( map_1%array( 0:map_1%nu-1, 0:map_1%nv-1, 0:map_1%nw-1 ), STAT=istat )
            IF ( istat == 0 ) THEN
                CALL messag( 'Logical map has been successfully allocated.', 'allocate_logical_map' )
            ELSE
                WRITE(*, '('' ALLOCATE_LOGICAL_MAP> istat= '', A)' ) TRIM ( int_to_c ( istat ) )
                CALL die( 'Failed to allocate logical map.', 'allocate_logical_map' )
            ENDIF
        ENDIF

    END SUBROUTINE allocate_logical_map

    SUBROUTINE deallocate_logical_map ( map_1 )
!
!       Purpose:
!       =======
!       Deallocates logical byte map
!
!       Date:          Programmer:       Description of changes:
!       ====           ==========        ======================
!       Dec 2003       Strokopytov B.    Original code
!       Nov 2005       Strokopytov B.    allocate_array added
!       Nov 2005       Strokopytov B.    ALLOCATED(sp_group) added
!
        TYPE(map_logical), INTENT(INOUT) :: map_1
!       Local variables:
        INTEGER                          :: istat

        IF ( ALLOCATED ( map_1%array ) ) THEN
            DEALLOCATE ( map_1%array, STAT = istat )

            IF ( istat == 0 ) THEN
                CALL messag ( 'Logical map was successfully deallocated.', 'deallocate_logical_map')
            ELSE
                WRITE(*,"(' DEALLOCATE_LOGICAL_MAP> ', 'istat= ', A)") TRIM( int_to_c ( istat ) )
                CALL die('Failed to deallocate locial map.', 'deallocate_logical_map')
            ENDIF

!           Free memory for other arrays:
            IF ( ALLOCATED ( map_1%SYM_DEORT ) ) CALL deallocate_array ( map_1%SYM_DEORT )
            IF ( ALLOCATED ( map_1%sp_group ) )  CALL deallocate_space_group ( map_1%sp_group )
            IF ( ALLOCATED ( map_1%atom_name ) ) CALL deallocate_array ( map_1%atom_name )
       ELSE
           CALL warn('Logical map has already been allocated.', 'deallocate_logical_map')
       ENDIF

    END SUBROUTINE deallocate_logical_map

END MODULE logical_map_manip
