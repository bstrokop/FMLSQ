MODULE atom_images
USE constants
USE fail
USE select_kinds
USE util
USE vectors
IMPLICIT NONE

INTERFACE ALLOCATED
    MODULE PROCEDURE allocated_atom_image
END INTERFACE
INTERFACE SIZE
    MODULE PROCEDURE size_of_atom_image
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_of_atom_images
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_array_of_atom_images
END INTERFACE

TYPE :: atom_image
    REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE :: ro
!   This assumes that map array does not contain more than 2**31 - 1 grid points:
    INTEGER(KIND=fb), DIMENSION(:), ALLOCATABLE :: uvw
END TYPE

PRIVATE
PUBLIC :: atom_image
PUBLIC :: ALLOCATED
PUBLIC :: SIZE
PUBLIC :: allocate_array
PUBLIC :: allocate_atom_image
PUBLIC :: deallocate_array
PUBLIC :: convert_to_1d

CONTAINS
    FUNCTION allocated_atom_image ( atom_image_1 )
        LOGICAL                      :: allocated_atom_image
        TYPE(atom_image), INTENT(IN) :: atom_image_1
        allocated_atom_image = ALLOCATED ( atom_image_1%ro )
    END FUNCTION allocated_atom_image

    FUNCTION size_of_atom_image ( atom_image_1 )
        INTEGER                      :: size_of_atom_image
        TYPE(atom_image), INTENT(IN) :: atom_image_1
        size_of_atom_image = SIZE ( atom_image_1%ro )
    END FUNCTION size_of_atom_image

    SUBROUTINE allocate_atom_image(atom_image_1, n)
        TYPE(atom_image), INTENT(INOUT) :: atom_image_1
        INTEGER,          INTENT(IN)    :: n
        IF ( ALLOCATED ( atom_image_1 ) ) THEN
            CALL warn('Atom image is already allocated. Inaccurate programming.','allocate_atom_image')
            CALL deallocate_atom_image(atom_image_1)
        ENDIF
        CALL allocate_array(atom_image_1%ro, n)
        CALL allocate_array(atom_image_1%uvw, n)
    END SUBROUTINE allocate_atom_image

    SUBROUTINE deallocate_atom_image(atom_image_1 )
        TYPE(atom_image), INTENT(INOUT)        :: atom_image_1
        IF ( ALLOCATED ( atom_image_1%ro ) ) THEN
            CALL deallocate_array(atom_image_1%ro)
        ENDIF
        IF ( ALLOCATED ( atom_image_1%uvw ) ) THEN
            CALL deallocate_array(atom_image_1%uvw) 
        ENDIF 
    END SUBROUTINE deallocate_atom_image

    SUBROUTINE allocate_array_of_atom_images(all_atom_images, natom)
        TYPE(atom_image), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: all_atom_images
        INTEGER,                                     INTENT(IN)    :: natom
!       Local variables:
        INTEGER                                                    :: istat
        CHARACTER(LEN=32), SAVE                                    :: srname = 'allocate_array_of_atom_images'
        IF ( ALLOCATED ( all_atom_images ) ) THEN
            CALL deallocate_array_of_atom_images(all_atom_images)
        ENDIF        
        ALLOCATE ( all_atom_images(natom), STAT=istat )
        IF ( istat /= 0 ) THEN
            WRITE(*,*) natom
            CALL die('OOPS. Failed to allocate array of atomic images...', srname)
        ENDIF
    END SUBROUTINE allocate_array_of_atom_images

    SUBROUTINE deallocate_array_of_atom_images(all_atom_images)
        TYPE(atom_image), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: all_atom_images
!       Local variables:
        INTEGER                                                    :: istat
        CHARACTER(LEN=32), SAVE                                    :: srname='deallocate_array_of_atom_images'
!       Counters:
        INTEGER                                                    :: i
        INTEGER(KIND=eb)                                           :: n
        n = 0
        IF ( ALLOCATED ( all_atom_images ) ) THEN
            IF ( debug > 1000 ) THEN
                CALL messag('Deallocating ALL_ATOM_IMAGES.', srname)
            ENDIF
            DO i = 1, SIZE ( all_atom_images )
                IF ( ALLOCATED ( all_atom_images(i) ) ) THEN
                    n = n + SIZE ( all_atom_images(i) ) 
                    CALL deallocate_atom_image ( all_atom_images(i) )
                ENDIF
            ENDDO
            DEALLOCATE ( all_atom_images, STAT=istat)
            IF ( istat /= 0 ) THEN
                CALL die('OOPS. Failed to allocate array of atomic images...', srname)
            ENDIF
        ELSE
            CALL die('Programming error. Array ALL_ATOM_IMAGES has been deallocated already...', &
                     srname)
        ENDIF

        IF ( debug > 1000 ) THEN
            WRITE(*,"(' DEALLOCATE_ARRAY_OF_ATOM_IMAGES> ', I8, ' MB of memory has been deallocated.')")&
                      NINT (n * 12.0_wp / 2 ** 20)
        ENDIF
    END SUBROUTINE deallocate_array_of_atom_images

    SUBROUTINE load_atom_image(atom_image_1, ro, uvw, dims)
!
!       Purpose:
!       =======
!       Allocates and fills arrays for ATOM_IMAGE type variable.
!       
!       Note that density array shoud have been allocated as:
!       den(0:nx-1,0:ny-1,0:nz-1)
!
!
        TYPE(atom_image),               INTENT(INOUT)  :: atom_image_1
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)     :: ro 
        TYPE(vector_int), DIMENSION(:), INTENT(IN)     :: uvw
        INTEGER,          DIMENSION(:), INTENT(IN)     :: dims
        IF ( .NOT. ALLOCATED ( atom_image_1 ) ) THEN
            CALL allocate_atom_image(atom_image_1, SIZE ( ro ) )
        ENDIF
        atom_image_1%ro  = ro

!       Fortran style of converting to 1D arrays:
        atom_image_1%uvw = convert_to_1d ( uvw, dims(1), dims(2), 1)
!       For C:
!       nuvw = uvw; temp = nuvw(1); nuvw(1) = nuvw(3); nuvw(3) = temp;
!       atom_image_1%uvw = convert_to_1d ( uvw, dims(3), dims(2), 1)
!       

    END SUBROUTINE load_atom_image

    SUBROUTINE dump_atom_image_content(atom_image_1)
        TYPE(atom_image), INTENT(IN) :: atom_image_1
!       Counter: 
        INTEGER                      :: i

        IF ( ALLOCATED ( atom_image_1 ) ) THEN
            DO i = 1, SIZE ( atom_image_1 ) 
                WRITE(*,"(' PRINT_ATOM_IMAGE_CONTENT> ', 'index(1d)=', I10, ' ro= ', ES9.2)") &
                          atom_image_1%uvw(i), atom_image_1%ro(i)
            ENDDO
        ENDIF
    END SUBROUTINE dump_atom_image_content

    ELEMENTAL FUNCTION convert_to_1d ( uvw, nx, ny, base )
!
!       Purpose:
!       =======
!       Packs three indices of map grid point to one-dimensional
!       index to be used as one-dimensional array.
!       nx * ny * iz + nx * iy + ix + base =
!       nx * ( ny * iz  + iy ) + ix + base
!
!       N.B. base should be zero for array starting from zero
!       and 1 for standard fortran array.
!
!       For C routines both UVW and DIMS must be inverted/swapped.
!
!       Date:          Programmer:        History of changes:
!       ====           ==========         ==================
!       Oct 2008       B.Strokopytov      Original code
!
!
        INTEGER(KIND=fb)                  :: convert_to_1d
        TYPE(vector_int),      INTENT(IN) :: uvw
        INTEGER,               INTENT(IN) :: nx
        INTEGER,               INTENT(IN) :: ny
        INTEGER,               INTENT(IN) :: base
!       Local variables:
        INTEGER, DIMENSION(3)             :: nuvw
        INTEGER                           :: temp

        nuvw = uvw
        convert_to_1d = nx * ny * nuvw(3) + nx * nuvw(2) + nuvw(1) + base
    END FUNCTION convert_to_1d

    SUBROUTINE check_convert_to_1d
        INTEGER,          PARAMETER                             :: nx=4
        INTEGER,          PARAMETER                             :: ny=6
        INTEGER,          PARAMETER                             :: nz=4
        REAL(KIND=wp),    DIMENSION(0:nx-1,0:ny-1,0:nz-1)       :: ro
        INTEGER,          DIMENSION(3)                          :: dims
        INTEGER,          DIMENSION(:), ALLOCATABLE             :: ind
        TYPE(vector_int), DIMENSION(:), ALLOCATABLE             :: uvw
        INTEGER,          DIMENSION(3)                          :: nuvw
        REAL(KIND=wp)                                           :: my_ro
!       Counters:
        INTEGER                                                 :: ix
        INTEGER                                                 :: iy
        INTEGER                                                 :: iz
        INTEGER                                                 :: k 
!       External:
        REAL(KIND=wp), EXTERNAL                                 :: get_ro
        REAL(KIND=wp), EXTERNAL                                 :: convolution_no_1
        REAL(KIND=sp)                                           :: t1,t0        

!       Fill array RO with random numbers:      
        CALL RANDOM_NUMBER(ro)
        CALL allocate_array ( uvw, nx*ny*nz )
        CALL allocate_array ( ind, nx*ny*nz )
        k = 0
!       Intentionally use "incorrect" looping to inroduce some disorder in 1-D indices:
        DO ix=0,nx-1
            DO iy=0,ny-1
                DO iz=0,nz-1
                    k = k + 1
                    uvw(k) = (/ix,iy,iz/)
                ENDDO
            ENDDO
        ENDDO

        ind = convert_to_1d ( uvw, nx, ny, 1)
        WRITE(*,*) ind

        DO k = 1, nx*ny*nz
            my_ro = get_ro(ro,ind(k))
            nuvw = uvw(k)
            WRITE(*,"(2ES9.2, I10, 5X, 3I4)") ro(nuvw(1),nuvw(2),nuvw(3)), my_ro, ind(k),  nuvw
            IF ( ro(nuvw(1),nuvw(2),nuvw(3)) /= my_ro ) STOP 'Error'
        ENDDO
!       Check convolution in libro.a:
        ind = (/ (k,k=1,nx*ny*nz)  /)

        CALL CPU_TIME ( t0 )
!        DO k = 1, 400000 * 50
        DO k = 1, 18500 * 110
            my_ro = convolution_no_1(nx * ny * nz, ro, ind, ro )
        ENDDO
        CALL CPU_TIME ( t1 )

        WRITE(*,*) ' convolution_no_1:', my_ro
        WRITE(*,*) t1 - t0, ' s'

        CALL CPU_TIME ( t0 )
        DO k = 1, 400000 * 50
            my_ro = SUM ( ro * ro )
        ENDDO
        CALL CPU_TIME ( t1 )
        WRITE(*,*) ' F90 dot product: ', SUM ( ro * ro )
        WRITE(*,*) t1 - t0, ' s '
    END SUBROUTINE check_convert_to_1d

END MODULE atom_images 
