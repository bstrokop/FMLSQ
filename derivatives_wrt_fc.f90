MODULE derivatives_wrt_fc
USE fail 
USE mtz_io
USE select_kinds
IMPLICIT NONE
CONTAINS
    SUBROUTINE fofc_cc_deriv_wrt_fc ( mtz_1, x_cc )
!
!       Purpose:
!       =======
!       Calculate correlation coefficient derivatives with respect to Fc, 
!       constructing coefficients for map to be used for convolution with derivative of atomic density
!
!       Other targets can be introduced later.
!
!       Date:
!       ====
!       September 2004
!
        TYPE(mtz),               INTENT(INOUT)         :: mtz_1
        REAL(KIND=wp),           INTENT(OUT), OPTIONAL :: x_cc
!       Local variables:
        INTEGER                                        :: n
        INTEGER                                        :: n1
        REAL(KIND=wp)                                  :: yo
        REAL(KIND=wp)                                  :: yc
!       Counters:
        INTEGER                                        :: i
        INTEGER                                        :: itest
!       Sums:
        REAL(KIND=wp)                                  :: sumoc
        REAL(KIND=wp)                                  :: sumo2
        REAL(KIND=wp)                                  :: sumc2
        REAL(KIND=wp)                                  :: sumyo
        REAL(KIND=wp)                                  :: sumyc
!       Local CC:
        REAL(KIND=wp)                                  :: a
        REAL(KIND=wp)                                  :: c
        REAL(KIND=wp)                                  :: cc
        REAL(KIND=wp)                                  :: d
        REAL(KIND=wp)                                  :: dCCdfc
!       REAL(KIND=wp)                                  :: p_prime
!       REAL(KIND=wp)                                  :: q_prime
!
        IF (.NOT. ALLOCATED(mtz_1%fo_in_P1)) THEN
            CALL die('mtz_1%fo_in_P1 array has not been initialized properly.',  'fofc_corr_coef_fc_deriv_wp')
        ELSE IF (.NOT. ALLOCATED(mtz_1%fc_in_P1)) THEN
            CALL die('fcalc array has not been initialized properly.', 'fofc_corr_coef_fc_deriv_wp')
        ELSE
            n = SIZE(mtz_1%fo_in_P1)
            IF ( SIZE (mtz_1%fc_in_P1) /= n ) THEN
                WRITE(*,'('' FOFC_CORR_COEF_FC_DERIV_SP> '','' fobs size ='', I6, '' fcalc size ='', I6)')&
                n, SIZE(mtz_1%fc_in_P1)
                CALL die ( 'Size of fobs and fcalc arrays differ.', 'fofc_corr_coef_fc_deriv_wp' )
            ENDIF
        ENDIF

        sumoc = 0.0_wp
        sumo2 = 0.0_wp
        sumc2 = 0.0_wp
        sumyo = 0.0_wp
        sumyc = 0.0_wp
!       sumw  = 0.0_wp

        IF( .NOT.ALLOCATED (mtz_1%test_in_P1) ) THEN
            CALL die ( 'test_in_P1 array has not been allocated properly.', 'fofc_corr_coef_fc_deriv_wp' )
        ENDIF

!       Summation over hemisphere - no multiplicity factors needed:
        itest = 0
        DO i = 1, n

!           No need to bother about derivative if it is test reflection:
            IF (mtz_1%test_in_P1(i)) THEN
                itest = itest + 1
                CYCLE
            ENDIF

!           Accumalate unweighted sums for correlation:
            yo    = mtz_1%fo_in_P1(i)
            yc    = ABS(mtz_1%fc_in_P1(i))
!           sumw  = sumw  + mtz_1%mult_in_P1(i)
            sumyo = sumyo + yo
            sumyc = sumyc + yc
            sumc2 = sumc2 + yc * yc
            sumo2 = sumo2 + yo * yo
            sumoc = sumoc + yo * yc
        ENDDO
        n1 = n - itest
        a = n1 * sumoc - sumyo * sumyc
        c = n1 * sumo2 - sumyo * sumyo
        d = n1 * sumc2 - sumyc * sumyc

!       This won't affect CC but will affect derivatives

        cc = a / SQRT(c * d)
        IF ( PRESENT(x_cc) ) x_cc = cc 
        DO i = 1, n

!           No need to bother about derivative if it is test reflection:
            IF (mtz_1%test_in_P1(i)) THEN
                mtz_1%fc_in_P1(i) = 0.0_wp        ! BUG FIXED OCT 2004
                CYCLE
            ENDIF
            yo = mtz_1%fo_in_P1(i)
            yc = ABS(mtz_1%fc_in_P1(i))

!           Beautiful symmetric form for the first derivative wrt Fc(h_prime) 
!           Minus sign arises since we want to minimize 1 - cc, i.e maximize cc
            dCCdfc  = -cc * ( (n1 * yo - sumyo)/a - (n1 * yc  - sumyc)/d ) * mtz_1%mult_in_P1(i)

!           Multiply with the phase of fcalc and store it in mtz_1%fc_in_P1 array
!           dccdfc is a real function. In case of complex function derivative we would need to take complex conjugate. 
!           See "A memory-efficient FFT algorithm" by A.T.Brunger (1989). Acta Cryst. A45

!           This is by far more efficient than mtz_1%fc_in_P1(i) = dCCdfc * phase( phase_of( mtz_1%fc_in_P1(i) ) ):
            IF ( yc /= 0.0_wp ) THEN
                mtz_1%fc_in_P1(i) = (dCCdfc / yc ) * mtz_1%fc_in_P1(i)
            ELSE
                mtz_1%fc_in_P1(i) = 0.0_wp   ! phase undefined for zero amplitude (can be anything)
            ENDIF
!           WRITE(*,*) ' i=', i, ' yo =', yo, ' yc=', yc, ' mtz_1%fc_in_P1(i)=', mtz_1%fc_in_P1(i)
        ENDDO
        
    END SUBROUTINE fofc_cc_deriv_wrt_fc

END MODULE derivatives_wrt_fc
