MODULE aniso_derivatives
USE aniso_manip
USE mtz_io
IMPLICIT NONE
CONTAINS
    SUBROUTINE bulk_aniso_deriv(p, np, mtz_1, f0, only_fun, A, dv)
!
!       Purpose:
!       =======
!       Calculate first derivatives, linearized second derivatives and target function value.
!
!       Date:              Programmer:               History of changes:
!       ====               ==========                ==================
!       MAR 2010           B.Strokopytov             Original code
!       APR 2010           B.Strokopytov             Added calculation of SC_IN_P1, FPART_IN_P1 ETC.
!   
!       P(1)  is an overall scaling factor
!       P(2)  ksol (solvent scale factor around 0.35 e/A**3)
!       P(3)  Bsol (solvent thermal factor around 46 A**2)
!       P(4) - P(9) overall anisothermal parameters arranged as the following matrix:
!
!       |P(4)P(7)P(8)|
!       |..  P(5)P(9)|
!       |        P(6)|
!
!
!       DV is a vector of first derivatives.
!       A  normal matrix
!       f0 target function value
!     
        REAL(KIND=wp), DIMENSION(np),    INTENT(INOUT) :: p
        INTEGER,                         INTENT(IN)    :: np
        TYPE(mtz),                       INTENT(INOUT) :: mtz_1
        REAL(KIND=wp),                   INTENT(out)   :: f0
        LOGICAL,                         INTENT(IN)    :: only_fun
        REAL(KIND=wp), DIMENSION(np),    INTENT(INOUT) :: dv
        REAL(KIND=wp), DIMENSION(np,np), INTENT(INOUT) :: A
!       Local arrays and variables:
        REAL(KIND=wp), DIMENSION(3,3)                  :: BCART
        REAL(KIND=wp), DIMENSION(3)                    :: v
        REAL(KIND=wp), DIMENSION(3)                    :: temp
        REAL(KIND=wp)                                  :: qform
        REAL(KIND=wp)                                  :: exp_qform
        REAL(KIND=wp), DIMENSION(np,1)                 :: dfdp
        REAL(KIND=wp)                                  :: dy
        REAL(KIND=wp)                                  :: yo
        REAL(KIND=wp)                                  :: rf
!
        COMPLEX(KIND=wp)                               :: yc
        COMPLEX(KIND=wp)                               :: ybulk
!
        REAL(KIND=wp)                                  :: sqrtw
        INTEGER                                        :: nr
!       Counters:
        INTEGER                                        :: i
        INTEGER                                        :: j
!       
        CHARACTER(LEN=32),                     SAVE    :: srname='bulk_aniso_deriv'
               

!       Prepare BCART (transform 1d array in 2d):
        BCART = .CONVERT. p(4:9)
        IF( debug > 10000 ) THEN
            DO i = 1, 3
                WRITE(*,"(3F6.2)") (BCART(i,j),j=1,3)
            ENDDO
        ENDIF

!       Initialize function and derivative arrays
        f0 = 0.0_wp
        IF ( .NOT. only_fun ) THEN
            dv = 0.0_wp
            A = 0.0_wp
        ENDIF

        nr = SIZE ( mtz_1%fo_in_P1 )

!       FIXME: Maybe need parallelization?:
        rf = 0.0
        DO i = 1, nr

!           Calculate bulk solvent contribution applying ksol and Bsol:
            ybulk = p(2) * mtz_1%fpart_in_P1(i) * EXP ( -0.25_wp * p(3) * mtz_1%s_in_P1(i) )

!           Save bulk contribution:
!!            IF ( only_fun ) THEN
!!                mtz_1%fpart_in_P1(i) = ybulk
!!            ENDIF

!           Aniso exponent:
            v         = (/ mtz_1%HO(i), mtz_1%KO(i), mtz_1%LO(i) /)
            temp      = MATMUL ( BCART, v )
            exp_qform = EXP ( -0.25_wp * DOT_PRODUCT ( v, temp ) )

!           Save overall scale multiplited by aniso contirbution:
            IF ( only_fun ) THEN
                mtz_1%sc_in_P1(i) = p(1) * exp_qform
            ENDIF

            yo = ABS ( mtz_1%fo_in_P1(i) )

!           Apply overall scale and aniso scaling adding contribution from solvent mask:
            yc = p(1) * (mtz_1%fc_in_P1(i) + ybulk) * exp_qform

            dy = yo - ABS ( yc )

            sqrtw = SQRT (  2.0_wp * mtz_1%weight_in_P1(i) )

!           Accumulate function (half of a sphere):
            f0 = f0 + mtz_1%weight_in_P1(i) * dy ** 2
            rf = rf + ABS ( dy )

            IF ( .NOT. only_fun ) THEN

!               Overall scale factor and bulk solvent derivatives:
                dfdp(1,1) = ABS ( yc ) / p(1)
                dfdp(2,1) = 0.5_wp * p(1) * exp_qform &
                          * REAL ( yc * CONJG ( ybulk ) + CONJG ( yc ) * ybulk, KIND=wp ) &
                          / ( ABS ( yc ) * p(2) )
                dfdp(3,1) =  -0.25_wp * dfdp(2,1) * p(2) * mtz_1%s_in_P1(i)

                dfdp(4,1) = -0.25_wp * ABS ( yc ) * v(1) ** 2
                dfdp(5,1) = -0.25_wp * ABS ( yc ) * v(2) ** 2
                dfdp(6,1) = -0.25_wp * ABS ( yc ) * v(3) ** 2

                dfdp(7,1) = -0.5_wp * ABS ( yc ) * v(1) * v(2)
                dfdp(8,1) = -0.5_wp * ABS ( yc ) * v(1) * v(3)
                dfdp(9,1) = -0.5_wp * ABS ( yc ) * v(2) * v(3)

                dfdp = sqrtw * dfdp

!               Accumulate LSQ matrix (possible problems with parallelization):
                A = A + MATMUL ( dfdp, TRANSPOSE ( dfdp ) )

!               Accumulate derivatives:
                dv(1:9) = dv(1:9) + sqrtw * dy * dfdp(1:9,1)

            ENDIF
        ENDDO
        IF ( debug > 10000 ) THEN
            WRITE(*,"( ' R-factor=',F10.4)") rf / SUM ( mtz_1%fo_in_P1) 
        ENDIF
    END SUBROUTINE bulk_aniso_deriv

    SUBROUTINE calculate_fpart(p, mtz_1)
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)        :: p
        TYPE(mtz),                   INTENT(INOUT)     :: mtz_1
!       Local variables:
        REAL(KIND=wp), DIMENSION(3,3)                  :: BCART
        REAL(KIND=wp), DIMENSION(3)                    :: v
        REAL(KIND=wp), DIMENSION(3)                    :: temp
        REAL(KIND=wp)                                  :: exp_qform
        INTEGER                                        :: nr
!       Counters:
        INTEGER                                        :: i
        INTEGER                                        :: j 
        CHARACTER(LEN=32),                      SAVE   :: srname='calculate_fpart'
                  
!       Checkz:
        IF ( .NOT. ALLOCATED ( mtz_1%fpart ) ) THEN
            CALL die('Programming error. Array MTZ_1%FPART has not been allocated.', srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%fpart_in_P1 ) ) THEN
            CALL die('Programming error. Array MTZ_1%FPART_IN_P1 has not been allocated.', srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%sc ) ) THEN
            CALL die('Programming error. Array MTZ_1%SC has not been allocated.', srname)
        ELSE IF ( .NOT. ALLOCATED ( mtz_1%sc_in_P1 ) ) THEN
            CALL die('Programming error. Array MTZ_1%SC_IN_P1 has not been allocated.', srname)
        ENDIF

!       Prepare BCART (transform 1d array in 2d):
        BCART = .CONVERT. p(4:9)
        IF( debug > 10000 ) THEN
            DO i = 1, 3
                WRITE(*,"(3F6.2)") (BCART(i,j),j=1,3)
            ENDDO
        ENDIF

!       Initialize arrays mtz_1%fpart and mtz_1%sc (for statistics):
        nr = SIZE ( mtz_1%hkl )
        DO i = 1, nr
            mtz_1%fpart(i) = p(2) * mtz_1%fpart(i) * EXP ( -0.25_wp * p(3) * mtz_1%s(i) )
            v = TRANSPOSE ( mtz_1%DEORT ) * mtz_1%hkl(i)
            temp = MATMUL ( BCART, v )
            exp_qform = EXP ( -0.25_wp * DOT_PRODUCT ( v, temp ) )

!           Save overall scale multiplied by aniso contribution:
            mtz_1%sc(i) = p(1) * exp_qform

        ENDDO

!       DO the same in P1:
        nr = SIZE ( mtz_1%hkl_in_P1 )
        DO i = 1, nr
            mtz_1%fpart_in_P1(i) = p(2) * mtz_1%fpart_in_P1(i) * EXP ( -0.25_wp * p(3) * mtz_1%s_in_P1(i) )
            v = TRANSPOSE ( mtz_1%DEORT ) * mtz_1%hkl_in_P1(i)
            temp = MATMUL ( BCART, v )
            exp_qform = EXP ( -0.25_wp * DOT_PRODUCT ( v, temp ) )

!           Save overall scale multiplied by aniso contribution:
            mtz_1%sc_in_P1(i) = p(1) * exp_qform

        ENDDO

!       FIXME parallelization is on the way...

    END SUBROUTINE calculate_fpart
 
END MODULE aniso_derivatives
