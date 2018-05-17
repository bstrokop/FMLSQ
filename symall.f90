MODULE symall 
USE aniso_manip
USE mtz_io
USE basic_pdb_manip
USE vectors
IMPLICIT NONE
CONTAINS
    SUBROUTINE sfcalc(hkl, sinthl2, sp_group, pdb_2, apart, bpart, dadp, dbdp)
!
!       Purpose:
!       =======
!       Calculates structure factors and vector of DADP, DBDP derivatives
!       w.r.t. single HKL reflection.
!
!       Date:              Programmer:               History of changes:
!       ====               ==========                ==================
!       Oct 2008           B.Strokopytov             Original code
!       Mar 2009           B.Strokopytov             Anisotropic Displacement factors added.
!
        TYPE(vector_int),            INTENT(IN)    :: hkl
        REAL(KIND=wp),               INTENT(IN)    :: sinthl2
        TYPE(space_group),           INTENT(IN)    :: sp_group
        TYPE(pdb),                   INTENT(IN)    :: pdb_2
        REAL(KIND=wp),               INTENT(INOUT) :: apart
        REAL(KIND=wp),               INTENT(INOUT) :: bpart
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: dadp
        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: dbdp
!        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: aaa
!        REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: bbb

!       Local variables:
        TYPE(vector_int)                           :: hklm
!       Reflection index in orthogonal system:
        REAL(KIND=wp), DIMENSION(3)                :: Ho
        TYPE(vector)                               :: H
!       Mixed indices:
        REAL(KIND=wp), DIMENSION(6)                :: HH
!       Gaussian constants:
        REAL(KIND=wp), DIMENSION(5)                :: a
        REAL(KIND=wp), DIMENSION(5)                :: b
        INTEGER                                    :: ngauss
!       Phases:
        REAL(KIND=wp)                              :: aa
        REAL(KIND=wp)                              :: bb
        REAL(KIND=wp)                              :: trans
        TYPE(vector)                               :: xyz
        REAL(KIND=wp)                              :: phase 
        REAL(KIND=wp)                              :: scattering_power
        REAL(KIND=wp)                              :: sq_factor
!       Atom types:
        INTEGER                                    :: ityp
!       Pointers:
        INTEGER                                    :: indxyz
        INTEGER                                    :: indbiso
        INTEGER                                    :: indocc
        INTEGER                                    :: indu
        LOGICAL                                    :: aniso
        LOGICAL                                    :: refine_aniso_only
        TYPE(matrix)                               :: U
        REAL(KIND=wp)                              :: qform
!       Counters:
        INTEGER                                    :: iat
        INTEGER                                    :: i
        INTEGER                                    :: k
        INTEGER                                    :: l 
        INTEGER                                    :: m

        aniso = ALLOCATED ( pdb_2%U )

        apart = 0.0_wp
        bpart = 0.0_wp

        dadp = 0.0_wp
        dbdp = 0.0_wp
        indxyz = 0

!       Initialize derivative arrays:
        DO iat = 1, SIZE ( pdb_2 )

!           If aniso array is allocated and non-zero convert to 3x3 matrix:            
            IF ( aniso ) THEN
                IF ( ANY ( pdb_2%U(iat)%U > 0.0_wp ) ) THEN
                    U = pdb_2%U(iat)
                ENDIF
            ENDIF

!           De-orthogonalise:
            xyz = pdb_2%DEORT * pdb_2%xyz(iat)

!           Atom type:
            ityp   = pdb_2%atom_type(iat)
            ngauss = pdb_2%atomsf_table%record(ityp)%ngauss

            a(1:ngauss) = pdb_2%atomsf_table%record(ityp)%a(1:ngauss)
            b(1:ngauss) = pdb_2%atomsf_table%record(ityp)%b(1:ngauss)
  
            scattering_power = SUM ( a(1:ngauss) * EXP ( -b(1:ngauss) * sinthl2 ) )

!           We must have this all the time for atoms with no aniso values (mixed refinement):
            sq_factor = pdb_2%occ(iat) * scattering_power * EXP( -pdb_2%biso(iat) * sinthl2 )

!           Apply symmetry operations (computing in true space group):                 
            symmetry:DO m = 1, sp_group%number_of_symops

                hklm  = sp_group%SYM_HKL(m) * hkl

!               Calculate orthogonal derivatives:
                H = TRANSPOSE ( pdb_2%DEORT ) * hklm

!               Move to array:
                Ho  = H

!               Make sure aniso arrays are allocated before further checking/computing:
                IF ( aniso ) THEN
!                   If aniso thermal factors differ from zeroes use them:
                    IF ( ANY ( pdb_2%U(iat)%u /= 0.0_wp ) )  THEN
                        qform = H .DOT. ( U * H )
                        sq_factor = pdb_2%occ(iat) * scattering_power * EXP ( -2.0_wp * pi_squared * qform )
                    ENDIF
                ENDIF

                trans = hkl .DOT. .SYMV. sp_group%sym(m)
                phase  = twopi * ( (hklm .DOT. xyz) + trans )

                aa = sq_factor * COS ( phase )
                bb = sq_factor * SIN ( phase )

!                aaa(iat) = aaa(iat) + aa
!                bbb(iat) = bbb(iat) + bb

                apart = apart + aa
                bpart = bpart + bb


!               Derivatives dadp (current reflection):
                IF ( ANY ( pdb_2%refinable(1:3,iat) ) ) THEN

                    l = 0
                    DO k = 1, 3
                        IF ( pdb_2%refinable(k,iat) ) THEN
                            l = l + 1
                            dadp(indxyz+l) = dadp(indxyz+l) - twopi * Ho(k) * bb
                            dbdp(indxyz+l) = dbdp(indxyz+l) + twopi * Ho(k) * aa
                        ENDIF
                    ENDDO

                ENDIF

                IF ( pdb_2%refinable(4, iat ) ) THEN
                    indbiso = indxyz + COUNT ( pdb_2%refinable(1:4, iat) )
                    dadp(indbiso) = dadp(indbiso) - sinthl2 * aa
                    dbdp(indbiso) = dbdp(indbiso) - sinthl2 * bb
                ENDIF

                IF ( pdb_2%refinable(5,iat) ) THEN
                    indocc = indxyz + COUNT ( pdb_2%refinable(1:5, iat) )
                    dadp(indocc) = dadp(indocc) + aa / pdb_2%occ(iat)
                    dbdp(indocc) = dbdp(indocc) + bb / pdb_2%occ(iat)
                ENDIF

!               Aniso arrays are not always allocated:
                IF ( aniso ) THEN

                    indu = indxyz + COUNT ( pdb_2%refinable(1:5,iat) )

                    IF ( ANY ( pdb_2%refinable(6:11,iat) ) ) THEN

                        HH(1) = -2.0_wp * pi_squared * Ho(1) * Ho(1)
                        HH(2) = -2.0_wp * pi_squared * Ho(2) * Ho(2)
                        HH(3) = -2.0_wp * pi_squared * Ho(3) * Ho(3)

                        HH(4) = -4.0_wp * pi_squared * Ho(1) * Ho(2)
                        HH(5) = -4.0_wp * pi_squared * Ho(1) * Ho(3)
                        HH(6) = -4.0_wp * pi_squared * Ho(2) * Ho(3)

                        l = 0 
                        DO k = 1, 6
                            IF ( pdb_2%refinable(5+k,iat) ) THEN 
                                l = l + 1
                                dadp(indu+l) = dadp(indu+l) + HH(k) * aa
                                dbdp(indu+l) = dbdp(indu+l) + HH(k) * bb
                            ENDIF
                        ENDDO

                    ENDIF
                ENDIF
            ENDDO symmetry

!           Increment pointer:
            indxyz = indxyz + COUNT ( pdb_2%refinable(1:11, iat) )
        ENDDO

    END SUBROUTINE sfcalc
END MODULE symall
