MODULE cell_manip
USE constants
USE fail
USE symmetry_manip
USE select_kinds
USE vectors
IMPLICIT NONE
CONTAINS    
    SUBROUTINE calcmetric ( cell, ORT, DEORT, TREAL, TRECIP, volume )
!
!       Purpose:
!       =======
!       Calculate the orthonormal - crystallographic transform
!
!       Date:          Programmer:          Description of changes:
!       ====           ==========           ======================
!       Sep 2003       B. Strokopytov       Based on the code by K. Cowtan
!       Oct 2003                            Adapted to F95
!       Sep 2005                            removed mtz type
!       Nov 2005                            Long standing bug with TRECIP tensor corrected
!                                           TRANSPOSE must be used
!
!       Note: cell angles assumed to be in degrees
!
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: cell
        TYPE(matrix),                INTENT(OUT)           :: ORT
        TYPE(matrix),                INTENT(OUT)           :: DEORT
        TYPE(tensor),                INTENT(OUT), OPTIONAL :: TREAL
        TYPE(tensor),                INTENT(OUT), OPTIONAL :: TRECIP
        REAL(KIND=wp),               INTENT(OUT), OPTIONAL :: volume
!
!       Local variables:
        REAL(KIND=wp)                                      :: a
        REAL(KIND=wp)                                      :: b
        REAL(KIND=wp)                                      :: c
        REAL(KIND=wp)                                      :: alph
        REAL(KIND=wp)                                      :: beta
        REAL(KIND=wp)                                      :: gamm
        REAL(KIND=wp)                                      :: sum
        REAL(KIND=wp)                                      :: cosas
        REAL(KIND=wp)                                      :: q123
        REAL(KIND=wp), DIMENSION(3,3)                      :: q
        REAL(KIND=wp), DIMENSION(3,3)                      :: qt
!
        q  = 0.0_wp
        qt = 0.0_wp
!
        a = cell(1)
        b = cell(2)
        c = cell(3)

        alph   = pi / 180.0_wp * cell(4)
        beta   = pi / 180.0_wp * cell(5)
        gamm   = pi / 180.0_wp * cell(6)

        sum = 0.5_wp * (alph + beta + gamm)
        IF ( PRESENT ( volume ) ) THEN
            volume = 2.0_wp * a * b * c &
                   * SQRT ( SIN ( sum ) * SIN ( sum - alph ) * SIN ( sum - beta ) * SIN ( sum - gamm ) )
        ENDIF
        cosas  = (COS ( beta ) * COS ( gamm ) - COS ( alph )) / (SIN ( beta ) * SIN ( gamm ))

!       orthogonalization matrix
        q(1,1) =  a
        q(1,2) =  b * COS ( gamm )
        q(1,3) =  c * COS ( beta )
        q(2,2) =  b * SIN ( gamm )
        q(2,3) = -c * SIN ( beta ) * cosas
        q(3,3) =  c * SIN ( beta ) * SQRT ( 1.0 - cosas**2 )

        ORT = q

!       and its inverse
        q123    = q(1,1) * q(2,2) * q(3,3)
        qt(1,1) = 1.0 / q(1,1)
        qt(1,2) = -q(1,2) * q(3,3) / q123
        qt(1,3) = (q(1,2) * q(2,3) - q(2,2) * q(1,3) ) / q123
        qt(2,2) = 1.0 / q(2,2)
        qt(2,3) = -q(1,1) * q(2,3) / q123
        qt(3,3) = 1.0 / q(3,3)

        DEORT = qt
!
        IF ( PRESENT ( volume ) ) THEN
            WRITE (*, 4) volume
        ENDIF
        CALL print_mat ( ORT, DEORT, ORT * DEORT )

!       Calculate real space metric tensor:
        IF ( PRESENT ( TREAL ) ) TREAL  = q

!       Calculate reciprocal space metric tensor:
        IF ( PRESENT ( TRECIP ) ) TRECIP = TRANSPOSE ( qt )

!       Print both tensors:
        IF ( PRESENT ( TREAL ) .AND. PRESENT ( TRECIP ) ) THEN
            CALL print_tensor( TREAL, TRECIP )
        ENDIF
 4      FORMAT ( ' CALCMETRIC> Unit cell volume:', F10.1, ' A**3' )
    END SUBROUTINE calcmetric

    SUBROUTINE check_sym_cell_consistency ( REAL_TENSOR, sp_group )
!
!       Purpose:
!       =======
!       Tests whether cell parameters consistent with space group symmetry operators
!
!       Date:              Programmer:       Description of changes:
!       ====               ==========        ======================
!       Oct 2004           B. Strokopytov    Lifted from X-plor and rewritten for F95
!       Oct 2005           B. Strokopytov    .SYMA. operator added to extract matrix part of symop
!
        TYPE(tensor),      INTENT(IN)           :: REAL_TENSOR
        TYPE(space_group), INTENT(IN)           :: sp_group
!       Local variables:
        REAL(KIND=wp), PARAMETER                :: relative_error = 0.1_wp ** 6
        TYPE(vector)                            :: test_vector
        REAL(KIND=wp)                           :: test_distance
        REAL(KIND=wp)                           :: delta
!       Counters:
        INTEGER                                 :: i
!
        IF ( SIZE( sp_group%sym ) == 1 ) RETURN
        CALL messag( ' ', 'check_sym_cell_consistency' )
        test_vector   = (/ 0.13_wp, 0.17_wp, 0.19_wp /)
        test_distance = REAL_TENSOR * test_vector

!       Loop over all symmetry operators:
        DO i = 2, SIZE(sp_group%sym)
!           Test whether test_distance has changed in any significant way when rotational operators are applied:
            delta  = ABS ( REAL_TENSOR * ( ( .SYMA. sp_group%sym(i) ) * test_vector ) - test_distance )
            IF ( delta > relative_error * test_distance ) THEN
                CALL die( 'Spacegroup is incompatible with cell dims.', 'check_sym_cell_consistency' )
            ENDIF
        ENDDO
        CALL messag( 'Looks good...', 'check_sym_cell_consistency' )
    END SUBROUTINE check_sym_cell_consistency

END MODULE cell_manip
