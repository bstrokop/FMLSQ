MODULE testing 
USE mtz_io
USE basic_pdb_manip
USE corr_coef_manip
USE extract_fo_fc
USE genden
USE map_fft
USE map_operations
!USE sf_table_manip
USE symall
USE vectors
!USE mr_list_manip
IMPLICIT NONE
CONTAINS
    SUBROUTINE dv(am, g, mtz_1, pdb_2)        
        REAL(KIND=wp), DIMENSION(:,:),          INTENT(INOUT) :: am
        REAL(KIND=wp), DIMENSION(:),            INTENT(INOUT) :: g
        TYPE(mtz),                              INTENT(INOUT) :: mtz_1
        TYPE(pdb),                              INTENT(IN)    :: pdb_2
!       Local variables:
        TYPE(vector_int)                                      :: hkl
        REAL(KIND=wp)                                         :: a
        REAL(KIND=wp)                                         :: b
        REAL(KIND=wp)                                         :: ca
        REAL(KIND=wp)                                         :: cb
        REAL(KIND=wp)                                         :: dy
        REAL(KIND=wp)                                         :: rho
        REAL(KIND=wp)                                         :: sc
        REAL(KIND=wp)                                         :: yo
        REAL(KIND=wp)                                         :: yc
        REAL(KIND=wp)                                         :: w
        REAL(KIND=wp)                                         :: wdy
!       Refinement:
        INTEGER                                               :: np
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE              :: dadp
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE              :: dbdp
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE              :: dfdp 
!       Pointers:
        INTEGER                                               :: indxyz                
        INTEGER                                               :: indbiso
        INTEGER                                               :: indocc
!       Counters:
        INTEGER                                               :: i
        INTEGER                                               :: iat
        INTEGER                                               :: k
        INTEGER                                               :: l
        INTEGER, DIMENSION(3)                                 :: uvw 
!       Test:
        REAL(KIND=wp)                                         :: r_factor

        np = COUNT ( pdb_2%refinable )

!       Allocate derivative arrays: 
        CALL allocate_array(dadp, np)
        CALL allocate_array(dbdp, np)
        CALL allocate_array(dfdp, np)
!
        CALL allocate_array(mtz_1%fc, SIZE ( mtz_1%hkl ) )

!       Set scaling:
        sc = mtz_1%fofc_scale
!        WRITE(*,*) ' sc=', sc
        CALL print_mat(pdb_2%DEORT)

!       Initialize:
        g  = 0.0_wp
        am = 0.0_wp

!       Loop over all reflections:
        DO i = 1, SIZE ( mtz_1%hkl )
            hkl = mtz_1%hkl(i)
            rho = mtz_1%s(i) / 4.0_wp
            yo = mtz_1%fo(i)

            CALL sfcalc(hkl, rho, mtz_1%sp_group, pdb_2, a, b, dadp, dbdp)

            yc = SQRT ( a ** 2 + b ** 2 )
            yc = sc * yc

!           Calculate ABS(fo) - ABS(fc) difference: 
            dy = 2.0_wp * ( yo - yc )
            
!           Weighting scheme (no such array yet):
            w   = mtz_1%weight(i)
!            w   = 1.0_wp
            wdy = w * dy

            a = sc * a
            b = sc * b

            ca = a / yc
            cb = b / yc

!           A and B has already been scaled:
            mtz_1%fc(i) = CMPLX(a,b,KIND=wp) / sc

!           dF/dp = (dF/dA) * (dA/dp) + (dF/dB) * (dB/dp):
!           Positional derivatives. Factor of 2 to account for the whole sphere and multiplicity of reflection:
            dfdp = ( 2.0_wp * SQRT ( mtz_1%mult(i) ) * sc ) * ( dadp * ca + dbdp * cb ) 

!           Accumulate gradient:
            g = g + wdy * SQRT ( mtz_1%mult(i) ) * dfdp

!           Accumulate full matrix:
            DO k = 1, np
                DO l = 1, np
                    am(k,l) = am(k,l) + dfdp(k) * dfdp(l)
                ENDDO
            ENDDO
            

        ENDDO
!        DO k = 1, np
!            WRITE(*,"(I4,' row SUM=',ES12.5)") k, SUM( am(k,1:np) )
!        ENDDO
        IF ( ALLOCATED ( mtz_1%fobs ) .AND. ALLOCATED ( mtz_1%fc ) ) THEN

!           mtz_1%fobs in this context is actually FFT Fcalcs:
            WRITE(*,*) ' fofc_corr_coef=', fofc_corr_coef(mtz_1%fobs, mtz_1%fc, mtz_1%mult)
            r_factor = simple_Rvalue(ABS(mtz_1%fobs), mtz_1%fc)
            WRITE(*,*) ' R-factor=', r_factor
        ENDIF

        CALL deallocate_array(dfdp)
        CALL deallocate_array(dadp)
        CALL deallocate_array(dbdp)

        IF ( ALLOCATED ( mtz_1%fc ) ) THEN
            CALL deallocate_array (mtz_1%fc)
        ENDIF

    END SUBROUTINE dv

END MODULE testing 
