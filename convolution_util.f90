MODULE convolution_util
USE basic_pdb_manip
USE constants
USE fail
USE select_kinds
IMPLICIT NONE
CONTAINS
    SUBROUTINE compose_joint_atom (ajo, bjo, ngauss, pdb_2, b_scale, iat, jat)
        INTEGER, PARAMETER                         :: MAXGAU=25
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: ajo
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: bjo
        INTEGER,                     INTENT(OUT)   :: ngauss
        TYPE(pdb),                   INTENT(IN)    :: pdb_2
        REAL(KIND=wp),               INTENT(IN)    :: b_scale
        INTEGER,                     INTENT(IN)    :: iat
        INTEGER,                     INTENT(IN)    :: jat
!       Local variables:
        INTEGER                                    :: ityp
        INTEGER                                    :: jtyp
        INTEGER                                    :: igauss
        INTEGER                                    :: jgauss
        REAL(KIND=wp)                              :: gaussian_scale
!       Aniso:
        LOGICAL                                    :: Us_present
        LOGICAL                                    :: aniso_active
!       Local arrays:
        REAL(KIND=wp), DIMENSION(5)                :: ai
        REAL(KIND=wp), DIMENSION(5)                :: bi
        REAL(KIND=wp), DIMENSION(5)                :: aj
        REAL(KIND=wp), DIMENSION(5)                :: bj
!       Counters:
        INTEGER                                    :: i
        INTEGER                                    :: j
        INTEGER                                    :: l
        CHARACTER(LEN=32),                    SAVE :: srname = 'compose_joint_atom'

        ityp   = pdb_2%atom_type(iat)
        igauss = pdb_2%atomsf_table%record(ityp)%ngauss

!       Copy gaussian constants for IAT atom:
        ai(1:igauss) = pdb_2%atomsf_table%record(ityp)%a(1:igauss)
        bi(1:igauss) = pdb_2%atomsf_table%record(ityp)%b(1:igauss) + b_scale 

!       Check:
        IF ( igauss /=2 .AND. igauss /= 5 ) THEN
            WRITE(*,*) 'igauss=', igauss, ' atom=', pdb_2%atom_name(iat)
            CALL die('Programming error. Abnormal number of gaussians detected...',&
                     srname)
        ENDIF

!       Save atom types for this pair:
        jtyp   = pdb_2%atom_type(jat)
        jgauss = pdb_2%atomsf_table%record(jtyp)%ngauss

        IF ( jgauss /=2 .AND. jgauss /= 5 ) THEN
            WRITE(*,*) 'jgauss=', jgauss, ' atom=', jat, pdb_2%atom_name(jat)
            CALL die('Programming error. Abnormal number of gaussians detected...',&
                      srname)
        ENDIF

!       Copy gaussian constants for JAT atom:
        aj(1:jgauss) = pdb_2%atomsf_table%record(jtyp)%a(1:jgauss)
        bj(1:jgauss) = pdb_2%atomsf_table%record(jtyp)%b(1:jgauss) + b_scale 

!       More safe to set ajo, bjo to zero though not absolutely necessary:
!       if we forget to set ngauss somewhere... that it's much better...
        ajo = 0.0_wp
        bjo = 0.0_wp        
        IF ( jat == iat ) THEN

!               Atom on the diagonal:
                gauss_outer:DO i = 1, igauss
                    gauss_inner:DO j = 1, i

                        IF ( i /= j ) THEN
                            gaussian_scale = 2.0_wp
                        ELSE
                            gaussian_scale = 1.0_wp
                        ENDIF

!                       Convert to one dimensional table (lower triange of 2-dim table):
                        l = (i * (i - 1)) / 2 + j

!                       Joint atom constants:
                        ajo(l) = gaussian_scale * ai(i) * ai(j)
                        bjo(l) = bi(i) + bi(j)

                    ENDDO gauss_inner
                ENDDO gauss_outer

!               Calculate ngauss:
                ngauss = (igauss * ( igauss + 1 )) / 2

!               Check number of gaussians:
                IF ( ngauss /= 3 .AND. ngauss /= 15 ) THEN
                    WRITE(*,*) ngauss
                    WRITE(*,*) iat, jat
                    CALL die('Programming error. Number of gaussians is not equal to the expected value.',&
                             srname)
                ENDIF

            ELSE

!               Different atoms:
                l = 0
                DO i = 1, igauss
                    DO j = 1, jgauss
                        l = l + 1

!                       Joint atom constants:
                        ajo(l) = ai(i) * aj(j)
                        bjo(l) = bi(i) + bj(j)
                    ENDDO
                ENDDO
!               Redefine ngauss:
                ngauss = l

!               Check number of gaussians:
                IF ( ngauss /= 4 .AND. ngauss /= 10 .AND. ngauss /= MAXGAU ) THEN
                    WRITE(*,*) ngauss
                    CALL die('Programming error. Number of gaussians is not equal to the expected value.',&
                             srname)
                ENDIF

            ENDIF

!           Debug:
            IF ( debug > 42 ) THEN
                WRITE(*,"(' ajo=', 15F8.2)") ajo(1:ngauss)
                WRITE(*,"(' bjo=', 15F8.2)") bjo(1:ngauss)
            ENDIF

    END SUBROUTINE compose_joint_atom

    SUBROUTINE compose_atom_block_table (table, block_size, pdb_2, iat, jat)
        INTEGER, DIMENSION(:,:), INTENT(INOUT) :: table
        INTEGER,                 INTENT(OUT)   :: block_size
        TYPE(pdb),               INTENT(IN)    :: pdb_2
        INTEGER,                 INTENT(IN)    :: iat
        INTEGER,                 INTENT(IN)    :: jat 
!       Local variables:
!       Counters:
        INTEGER                                :: elem
        INTEGER                                :: i
        INTEGER                                :: j

!       Prepare additional indexing table:
        table = 0
        elem  = 0

        IF ( debug > 40 ) THEN
            WRITE(*,*) ' indexing table:'
        ENDIF

        outer:DO i = 1, SIZE ( table, DIM = 1 ) 

            IF ( pdb_2%refinable(i,iat) ) THEN

                inner:DO j = 1, SIZE ( table, DIM = 1 )

!                   In case of atom block on main diagonal calculate only n*(n+1)/2 elements:
                    IF ( iat == jat .AND. i > j ) CYCLE inner

                    IF ( pdb_2%refinable(j,jat) ) THEN
                        elem = elem + 1
                        IF ( elem > 11 ** 2 ) THEN
                            WRITE(*,*) ' elem=', elem
                            CALL die('Unreasonable elem value.','compose_atom_block_table')
                        ENDIF
                        table(i,j) = elem
                    ENDIF

                ENDDO inner

            ENDIF

!           Debugging inside outer loop:
            IF ( debug > 40 ) THEN
                WRITE(*,"(' SPARSE_NORMAL_MATRIX> ', 11I4)") &
                (TABLE(i,j),j=1,SIZE (table,DIM=2) )
            ENDIF
        ENDDO outer

        block_size = elem
        
    END SUBROUTINE compose_atom_block_table

END MODULE convolution_util

