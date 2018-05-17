MODULE corr_coef_manip
!
!  Purpose:
!  =======
!  Calculates various correlation coefficients
!  plus some overlap functions
!
!  Date:         Programmer:           Description of changes:
!  ====          ==========            ======================
!  Nov 2003      B.Strokopytov         Original code
!  Nov 2005      B.Strokopytov         Almost complete rewrite
!                                      `weight' array are now optional in almost every routine
!                                      patt_overlap_function_in_P1_mtz removed                                      
!  Sep 2007      B.Strokopytov         Uncentered_ioic_corr_coef function added
!
USE fail
USE select_kinds
USE util
IMPLICIT NONE

INTERFACE fofc_corr_coef
    MODULE PROCEDURE abs_complex_fofc_corr_coef
    MODULE PROCEDURE abs_abs_fofc_corr_coef
    MODULE PROCEDURE int_weighted_fofc_corr_coef
    MODULE PROCEDURE test_abs_complex_fofc_corr_coef
    MODULE PROCEDURE complex_complex_fofc_corr_coef
END INTERFACE

INTERFACE ioic_corr_coef
    MODULE PROCEDURE abs_complex_ioic_corr_coef
END INTERFACE

INTERFACE ioic_overlap
    MODULE PROCEDURE abs_complex_ioic_overlap
END INTERFACE

INTERFACE fofc_overlap
    MODULE PROCEDURE complex_complex_fofc_overlap
END INTERFACE

INTERFACE r_factor
    MODULE PROCEDURE abs_complex_r_factor
END INTERFACE

CONTAINS
    SUBROUTINE shell_r_factor(fobs, fcalc, sig, sc, weight, mult, bin, rres, np)
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: sc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: sig
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: weight 
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: mult
        INTEGER,          DIMENSION(:), INTENT(IN)           :: bin
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: rres 
        INTEGER,                        INTENT(IN)           :: np
        INTEGER,                        SAVE                 :: icyc=0 
!       Local arrays:
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE          :: znum
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE          :: zden
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE          :: zsig
        REAL(KIND=wp),    DIMENSION(:), ALLOCATABLE          :: zw
        INTEGER,          DIMENSION(:), ALLOCATABLE          :: zn
        INTEGER,          DIMENSION(:), ALLOCATABLE          :: zmult
        COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE          :: scaled_fcalc
!       Local variables:
        INTEGER                                              :: nr
        INTEGER                                              :: n
        INTEGER                                              :: nshell
        REAL(KIND=wp)                                        :: yo
        COMPLEX(KIND=wp)                                     :: yc
        REAL(KIND=wp)                                        :: dy
        REAL(KIND=wp)                                        :: wdy
        REAL(KIND=wp)                                        :: sqrtw
        REAL(KIND=wp)                                        :: znm
        REAL(KIND=wp)                                        :: zdn
        REAL(KIND=wp)                                        :: wrznum
        REAL(KIND=wp)                                        :: wrzden
!       Counters:
        INTEGER                                              :: i
        CHARACTER(LEN=32),                          SAVE     :: srname='shell_r_factor'

        icyc = icyc + 1
        nr = SIZE ( fobs )

        IF ( SIZE ( fcalc ) /= SIZE ( fobs ) ) THEN
            WRITE(*,*) SIZE ( fcalc ), SIZE ( fobs )
            CALL die('Programming error. Sizes of arrays FCALC and FOBS differ.', srname)
        ENDIF

        IF ( SIZE ( fobs ) /= SIZE ( sc ) ) THEN
            WRITE(*,*) SIZE ( fcalc ), SIZE ( fobs )
            CALL die('Programming error. Sizes of arrays FCALC and FOBS differ.', srname)
        ENDIF

        IF ( ANY ( sc < 0.0_wp ) ) THEN
            CALL die('Programming error. SC array has negative components.', srname)
        ENDIF

        IF ( ANY ( weight <= 0.0_wp ) ) THEN
            WRITE(*,*) MINVAL ( weight )
            CALL die('Programming error. Cannot deal with negative or zero weights.', srname)
        ENDIF

        nshell = SIZE ( rres )

!       Allocate arrays for resolution breakdown:
        CALL allocate_array(zden, nshell)
        CALL allocate_array(znum, nshell)
        CALL allocate_array(zsig, nshell)
        CALL allocate_array(zn,   nshell)
        CALL allocate_array(zmult, nshell)
        CALL allocate_array(zw,    nshell)
!
        CALL allocate_array(scaled_fcalc, nr)
        scaled_fcalc = sc * fcalc

!       Initialized arrays:
        znum = 0.0_wp
        zden = 0.0_wp
        zsig = 0.0_wp
        zw   = 0.0_wp
        wrznum = 0.0_wp
        wrzden = 0.0_wp
        zn = 0
        zmult = 0
        
        DO i =1, nr
            yo = ABS ( fobs(i) )
            yc = sc(i) * ABS ( fcalc(i) )
            dy = ABS ( yo - yc )
            sqrtw = SQRT ( weight(i) )
            wdy = sqrtw * dy
            wrznum = wrznum + mult(i) * wdy ** 2
            wrzden = wrzden + mult(i) * (sqrtw * yo ) ** 2
!           Binning:
            n = bin(i)
            znum(n) = znum(n) + mult(i) * dy
            zden(n) = zden(n) + mult(i) * yo
            zsig(n) = zsig(n) + mult(i) * ABS ( sig(i) )
            zw(n) = zw(n) + mult(i) / sqrtw
            zmult(n) = zmult(n) + mult(i)

!           Number of reflections in the shell:
            zn(n) = zn(n) + 1
        ENDDO

        CALL messag(' ', srname)        

!       BUG CORRECTED APR 2008: non-linear scaling has to be applied:
        WRITE(*,2020) SUM ( mult * ABS ( fcalc ) ), fofc_corr_coef(fobs, scaled_fcalc, mult)
        WRITE(*,2030)
        WRITE(*,2040) SUM ( znum ), SUM ( zden ), SUM ( znum ) / SUM( zden )
        WRITE(*,2050) SQRT ( wrznum), SQRT ( wrzden ), SQRT ( wrznum / wrzden )
        WRITE(*,2062) TRIM ( int_to_c ( SIZE ( fobs ) ) )

!       WRZNUM is a half of target function:
        WRITE(*,2070) icyc, wrznum * 2, SQRT ( wrznum * 2 / ( SIZE ( fobs ) - np ) )

!       NO = SIZE ( fobs ) + geom restraints:FIXME
        WRITE(*,2080) SIZE ( fobs ), np, np, SUM ( znum ) / SUM ( zmult )
        WRITE(*,2100)

!       Initialize for sphere R-factors:
        znm = 0.0_wp
        zdn = 0.0_wp

        DO n = 1, nshell
            znm = znm + znum(n)
            zdn = zdn + zden(n)
            WRITE(*,"(1X, F6.2, I8, 3F10.2, 2F10.4)") &
            rres(n), zn(n), zsig(n)/zmult(n), zw(n)/zmult(n), znum(n)/zmult(n), znum(n) / zden(n), znm / zdn
        ENDDO
        WRITE(*,2110)

!       Deallocation:
        CALL deallocate_array(znum)
        CALL deallocate_array(zden)
        CALL deallocate_array(zsig)
        CALL deallocate_array(zw)
        CALL deallocate_array(zmult)
        CALL deallocate_array(zn)

        CALL deallocate_array(scaled_fcalc)
        CALL messag(' ', srname)
 
 2020 FORMAT(1X,'SUMYC ---------> ',T57,F17.4,/,&
     &       1X,'CORR COEFF ----> ',T66, F8.4)
 2030 FORMAT(1X,41X,'NUMERATOR   DENOMINATOR     R')
 2040 FORMAT(1X,'R FACTOR INCLUDING ZEROS',7X,F19.0,F13.0,1X,F8.4)
 2050 FORMAT(1X,'WEIGHTED R FACTOR INCLUDING ZEROS', F17.0,F13.0,1X,F8.4)
 2062 FORMAT(1X,'NUMBER OF STRUCTURE FACTORS READ FROM MTZ/CARDS = ',A)
 2070 FORMAT(1X,'DISCREPANCY FACTORS BASED ON PARAMETERS BEFORE CYCLE',I3, &
     &  /,T42,' SUM(W*(O-C)**2) IS ',T62,'=>',1PE10.2,' <=',&
     &  /,T34,' SQRT(SUM(W*(O-C)**2/NO-NP)   ',1PE10.2,&
     &    2X,'= ERROR IN AN OBSERVATION OF UNIT WEIGHT')
 2080 FORMAT(1X,'NO = ',I6,3X,'NV = ',I7,3X,'NP = ',I7,/, &
     &       ' THE AVERAGE |FO-FC| DISCREPANCY = ',F10.2)
 2100 FORMAT('   ',10(' -'),' RESOLUTION BREAKDOWN',10(' -')/&
     &'   DMIN    REFS     <SIGF>    <SIGF>   <FO-FC>    SHELL',&
     &'    SPHERE',/,                                           &
     &                 21X,'DATA      USED', 17X,          'R         R')
 2110 FORMAT(3X,30(' -'))
    END SUBROUTINE shell_r_factor
 
    FUNCTION abs_complex_r_factor( fobs, fcalc, sc, weight )
        REAL(KIND=wp)                                        :: abs_complex_r_factor
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: sc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
        CHARACTER(LEN=32),                          SAVE     :: srname='abs_complex_r_factor'

        IF ( SIZE ( fcalc ) /= SIZE ( fobs ) ) THEN
            WRITE(*,*) SIZE ( fcalc ), SIZE ( fobs )
            CALL die('Programming error. Sizes of arrays FCALC and FOBS differ.', srname)
        ENDIF

        IF ( SIZE ( fobs ) /= SIZE ( sc ) ) THEN
            WRITE(*,*) SIZE ( fcalc ), SIZE ( fobs )
            CALL die('Programming error. Sizes of arrays FCALC and FOBS differ.', srname)
        ENDIF

        IF ( ANY ( sc < 0.0_wp ) ) THEN
            CALL die('Programming error. SC array has negative components.', srname)
        ENDIF

        IF ( .NOT. PRESENT ( weight ) ) THEN
            abs_complex_r_factor = SUM ( ABS ( fobs - sc * ABS ( fcalc ) ) ) / SUM ( ABS ( fobs ) )
        ELSE
            abs_complex_r_factor = SUM ( weight * ABS ( fobs - sc * ABS ( fcalc ) ) ) / SUM ( weight * ABS ( fobs ) )
        ENDIF

    END FUNCTION abs_complex_r_factor

    FUNCTION abs_complex_fofc_corr_coef ( fobs, fcalc, weight )
        REAL(KIND=wp)                                        :: abs_complex_fofc_corr_coef
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        INTEGER                                              :: n
!       Sums
        REAL(KIND=wp)                                        :: sumoc
        REAL(KIND=wp)                                        :: sumo2
        REAL(KIND=wp)                                        :: sumc2
        REAL(KIND=wp)                                        :: sumyo
        REAL(KIND=wp)                                        :: sumyc
        REAL(KIND=wp)                                        :: sumw
!
        n = SIZE ( fobs )
        IF ( SIZE ( fcalc ) /= n ) THEN
            WRITE ( *, "(' ABS_COMPLEX_FOFC_CORR_COEF> ', ' fobs size=', I6, ' fcalc size=', I6)" )&
            n, SIZE(fcalc)
            CALL die ( 'Size of fobs and fcalc arrays differ.', 'abs_complex_fofc_corr_coef' )
        ENDIF
 
        IF ( PRESENT ( weight ) ) THEN

            IF ( SIZE ( weight ) /= SIZE ( fcalc ) ) THEN
                WRITE(*,*) SIZE ( weight ), SIZE ( fcalc )
                CALL die ( 'Size of fobs and fcalc arrays differ.', 'abs_complex_fofc_corr_coef' )
            ENDIF

            sumoc = SUM ( weight * fobs * ABS ( fcalc ) )
            sumo2 = SUM ( weight * fobs * fobs )
            sumc2 = SUM ( weight * REAL ( fcalc * CONJG ( fcalc ) ) )
            sumyo = SUM ( weight * ABS ( fobs  ) )
            sumyc = SUM ( weight * ABS ( fcalc ) )
            sumw  = SUM ( weight )

            abs_complex_fofc_corr_coef = ( sumw * sumoc - sumyo * sumyc ) &
                                       / SQRT ( ( sumw * sumo2 - sumyo * sumyo  ) * ( sumw * sumc2 - sumyc * sumyc ) )
        ELSE

            sumoc = SUM ( fobs * ABS ( fcalc ) )
            sumo2 = SUM ( fobs ** 2 )
            sumc2 = SUM ( fcalc * CONJG ( fcalc ) )
            sumyo = SUM ( fobs )
            sumyc = SUM ( ABS ( fcalc ) )

            abs_complex_fofc_corr_coef = ( n * sumoc - sumyo * sumyc ) &
                                       / SQRT ( ( n * sumo2 - sumyo * sumyo ) * ( n * sumc2 - sumyc * sumyc ) )
        ENDIF

    END FUNCTION abs_complex_fofc_corr_coef

    FUNCTION abs_abs_fofc_corr_coef ( fobs, fcalc, weight )
        REAL(KIND=wp)                                     :: abs_abs_fofc_corr_coef
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)           :: fobs
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp), DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        INTEGER                                           :: n
!       Sums:
        REAL(KIND=wp)                                     :: sumoc
        REAL(KIND=wp)                                     :: sumo2
        REAL(KIND=wp)                                     :: sumc2
        REAL(KIND=wp)                                     :: sumyo
        REAL(KIND=wp)                                     :: sumyc
        REAL(KIND=wp)                                     :: sumw

        n = SIZE ( fobs )
        IF ( SIZE ( fcalc ) /= SIZE( fobs ) ) THEN
            WRITE ( *, "(' ABS_ABS_FOFC_CORR_COEF> ', ' fobs size=', I6, ' fcalc size=',I6)" )&
            n, SIZE ( fcalc )
            CALL die ( 'Programming error. Size of fobs and fcalc arrays differ.', 'abs_abs_fofc_corr_coef' )
        ENDIF

        IF ( PRESENT ( weight ) ) THEN
            IF ( n /= SIZE ( weight ) ) THEN
                WRITE(*,*) n, SIZE ( weight )
                CALL die ( 'Programming error. Size of fobs and weight arrays differ.', 'abs_abs_fofc_corr_coef' )
            ENDIF
        ENDIF

        IF ( PRESENT ( weight ) ) THEN

            sumoc = SUM ( weight * fobs *  fcalc )
            sumo2 = SUM ( weight * fobs * fobs )
            sumc2 = SUM ( weight * fcalc * fcalc )
            sumyo = SUM ( weight * fobs  )
            sumyc = SUM ( weight * fcalc )
            sumw  = SUM ( weight )

            abs_abs_fofc_corr_coef = ( sumw * sumoc - sumyo * sumyc ) &
                                   / SQRT ( ( sumw * sumo2 - sumyo * sumyo  ) * ( sumw * sumc2 - sumyc * sumyc ) )

        ELSE
         
!           Note absence of ABS:
            sumoc = SUM ( fobs * fcalc ) 
            sumo2 = SUM ( fobs * fobs )
            sumc2 = SUM ( fcalc * fcalc )
            sumyo = SUM ( fobs )
            sumyc = SUM ( fcalc )

            abs_abs_fofc_corr_coef = ( n * sumoc - sumyo * sumyc ) &
                                   / SQRT ( ( n * sumo2 - sumyo * sumyo ) * ( n * sumc2 - sumyc * sumyc ) )
        ENDIF

    END FUNCTION abs_abs_fofc_corr_coef

    FUNCTION test_abs_complex_fofc_corr_coef ( fobs, fcalc, test, weight )
        COMPLEX(KIND=wp)                                     :: test_abs_complex_fofc_corr_coef
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        LOGICAL(KIND=1),  DIMENSION(:), INTENT(IN)           :: test
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        INTEGER,          DIMENSION(:), ALLOCATABLE          :: vector
        REAL(KIND=wp)                                        :: cc
        REAL(KIND=wp)                                        :: cc_free
!       Counters:
        INTEGER                                              :: i        

        IF ( SIZE ( fcalc ) /= SIZE ( fobs ) ) THEN
            WRITE ( *, "(' TEST_ABS_COMPLEX_FOFC_CORR_COEF> ', ' fobs size =',I6, ' fcalc size =', I6)" )&
            SIZE(fobs), SIZE(fcalc)
            CALL die ( 'Programming error. Size of fobs and fcalc arrays differ.', 'test_abs_complex_fofc_corr_coef' )
        ENDIF

        IF ( SIZE(test) /= SIZE ( fobs ) ) THEN
            WRITE(*,*) ' test size=', SIZE(test), ' fobs size=', SIZE(fobs)
            CALL die ( 'Programming error. Size of fobs and test arrays differ.', 'test_abs_complex_fofc_corr_coef' )
        ENDIF

        IF ( COUNT ( test ) == SIZE ( test ) .OR. COUNT ( test ) == 0 ) THEN
            WRITE(*,*) COUNT ( test ), SIZE ( test )
            CALL die ( 'Programming error. Cannot calculate CC-free.', 'test_abs_complex_fofc_corr_coef' )
        ENDIF

        IF ( PRESENT ( weight ) ) THEN       
            IF ( SIZE ( fobs ) /= SIZE ( weight ) ) THEN
                 WRITE(*,*) SIZE ( fobs ), SIZE ( weight )
                 CALL die ( 'Programming error. Size of fobs and weight arrays differ.', 'test_abs_complex_fofc_corr_coef' )
            ENDIF
        ENDIF

        CALL allocate_array ( vector, COUNT ( test ) )
        vector = PACK ( (/(i,i=1,SIZE(test))/), test )

        IF ( PRESENT ( weight ) ) THEN
            cc_free = fofc_corr_coef ( fobs(vector), fcalc(vector), weight(vector) )
        ELSE
            cc_free = fofc_corr_coef ( fobs(vector), fcalc(vector) )
        ENDIF
       
        CALL allocate_array ( vector, COUNT ( .NOT. test ) )
        vector = PACK ( (/(i,i=1,SIZE(test))/), .NOT. test )

        IF ( PRESENT ( weight ) ) THEN
            cc =  fofc_corr_coef ( fobs(vector), fcalc(vector), weight(vector) )
        ELSE
            cc = fofc_corr_coef ( fobs(vector), fcalc(vector) )
        ENDIF

!       Bug corrected -> deallocation added Nov 2005:

!       Free memory:
        CALL deallocate_array ( vector )

        test_abs_complex_fofc_corr_coef = CMPLX ( cc, cc_free, KIND=wp ) 

    END FUNCTION test_abs_complex_fofc_corr_coef

    FUNCTION complex_complex_fofc_corr_coef ( fobs, fcalc, weight )
!
!       Date:         Programmer:            Description of changes:
!       ====          ==========             ======================
!       Dec 2003      B.Strokopytov          Original code
!       Nov 2005      B.Strokopytov          Weights added
!
!
        REAL(KIND=wp)                                        :: complex_complex_fofc_corr_coef
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        INTEGER                                              :: n
!       Sums:
        REAL(KIND=wp)                                        :: sumoc
        REAL(KIND=wp)                                        :: sumo2
        REAL(KIND=wp)                                        :: sumc2

        n = SIZE ( fobs )
        IF ( SIZE ( fcalc ) /= n ) THEN
            WRITE ( *, "(' COMPLEX_COMPLEX_FOFC_CORR_COEF> ', ' fobs vector size=', I6, ' fcalc vector size=', I6)")&
            n, SIZE ( fcalc )
            CALL die ( 'Programming error: size of fobs and fcalc arrays differ.', 'complex_complex_fofc_corr_coef' )
        ENDIF

!       Added check for completeness FEB 2007 BVS:
        IF ( PRESENT ( weight ) ) THEN
            IF ( SIZE ( weight ) /= n )  THEN
                WRITE(*,*) SIZE ( weight ), n
                CALL die ( 'Programming error: size of fobs and weight arrays differ.', 'complex_complex_fofc_corr_coef' )
            ENDIF
        ENDIF

        IF ( PRESENT ( weight ) ) THEN
            sumoc = SUM ( REAL ( weight * fobs * CONJG(fcalc) ) )
            sumo2 = SUM ( REAL ( weight * fobs * CONJG(fobs) ) )
            sumc2 = SUM ( REAL ( weight * fcalc * CONJG(fcalc) ) )
        ELSE
            sumoc = SUM ( REAL ( fobs * CONJG(fcalc) ) )
            sumo2 = SUM ( REAL ( fobs * CONJG(fobs) ) )
            sumc2 = SUM ( REAL ( fcalc * CONJG(fcalc) ) )
        ENDIF

        complex_complex_fofc_corr_coef = sumoc / SQRT ( sumo2 * sumc2 )

    END FUNCTION complex_complex_fofc_corr_coef

    FUNCTION int_weighted_fofc_corr_coef ( fobs, fcalc, weight )
        REAL(KIND=wp)                              :: int_weighted_fofc_corr_coef
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN) :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN) :: fcalc
        INTEGER,          DIMENSION(:), INTENT(IN) :: weight
!       Local variables:
        INTEGER                                    :: n
!       Sums:
        REAL(KIND=wp)                              :: sumoc
        REAL(KIND=wp)                              :: sumo2
        REAL(KIND=wp)                              :: sumc2
        REAL(KIND=wp)                              :: sumyo
        REAL(KIND=wp)                              :: sumyc
        REAL(KIND=wp)                              :: sumw

!       Checkz:
        n = SIZE(fobs)
        IF ( SIZE ( fcalc ) /= n ) THEN
            WRITE ( *, "(' INT_WEIGHTED_FOFC_CORR_COEF> ', ' fobs size=', I6,' fcalc size=', I6)" )&
            n, SIZE ( fcalc )
            CALL die ( 'Size of fobs n fcalc arrays differ.', 'int_weighted_fofc_corr_coef' )
        ELSE IF ( SIZE ( weight ) /= n ) THEN
            WRITE ( *, "(' INT_WEIGHTED_FOFC_CORR_COEF> ', ' fobs size=', I6, ' weight size=', I6)" )&
            n, SIZE ( weight )
            CALL die ( 'Size of fobs n fcalc arrays differ.', 'int_weighted_fofc_corr_coef' )
        ENDIF

        sumoc = SUM ( weight * ABS ( fobs ) * ABS ( fcalc ) )
        sumo2 = SUM ( weight * fobs ** 2 )
        sumc2 = SUM ( weight * REAL ( fcalc * CONJG ( fcalc ), KIND=wp) )
        sumyo = SUM ( weight * ABS ( fobs ) )
        sumyc = SUM ( weight * ABS ( fcalc ) )

!       To avoid overflow use REAL:
        sumw  = SUM ( REAL ( weight, KIND=wp) )

        int_weighted_fofc_corr_coef = ( sumw * sumoc - sumyo * sumyc ) &
                                    / SQRT ( ( sumw * sumo2 - sumyo * sumyo  ) * ( sumw * sumc2 - sumyc * sumyc ) )

    END FUNCTION int_weighted_fofc_corr_coef

    FUNCTION simple_Rvalue (fo, fc)
        REAL(KIND=wp)                              :: simple_Rvalue
        REAL(KIND=wp), DIMENSION(:),    INTENT(IN) :: fo
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN) :: fc
!       Local variables:
        REAL(KIND=wp)                              :: fofc_scale
        fofc_scale    = SUM ( ABS ( fo ) * ABS ( fc ) ) / SUM ( ABS ( fc ) ** 2 )
        WRITE(*,*) ' simple_Rvalue fofc_scale=', fofc_scale
        simple_rvalue = SUM ( ABS ( fo - fofc_scale * ABS ( fc ) ) ) / SUM ( ABS ( fo ) )
    END FUNCTION simple_Rvalue

    FUNCTION abs_complex_ioic_corr_coef ( fobs, fcalc, weight )
        REAL(KIND=wp)                                        :: abs_complex_ioic_corr_coef
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        INTEGER                                              :: n
!       Sums:
        REAL(KIND=wp)                                        :: sumoc
        REAL(KIND=wp)                                        :: sumo2
        REAL(KIND=wp)                                        :: sumc2
        REAL(KIND=wp)                                        :: sumyo
        REAL(KIND=wp)                                        :: sumyc
        REAL(KIND=wp)                                        :: sumw

        n = SIZE ( fobs )
        IF ( SIZE ( fcalc ) /= n ) THEN
            WRITE ( *, "(' ABS_COMPLEX_IOIC_CORR_COEF> ', ' fobs size=', I6, ' fcalc size=', I6 )" ) &
            n, SIZE ( fcalc )
            CALL die ( 'Size of fobs n fcalc arrays differ.', 'abs_complex_ioic_corr_coef' )
        ENDIF

!       Bug corrected FEB 2007 BVS (must check presence of weight parameter):
        IF ( PRESENT ( weight ) ) THEN
            IF ( SIZE ( weight ) /= n ) THEN
                WRITE ( *, "(' ABS_COMPLEX_IOIC_CORR_COEF> ', ' fobs size=', I6, ' weight size=', I6 )") &
                n, SIZE ( weight )
                CALL die ( 'Size of fobs & weight arrays differ.', 'abs_complex_ioic_corr_coef' )
            ENDIF
        ENDIF
 
        IF ( PRESENT ( weight ) ) THEN
!           BUG corrected weighting was missing in 5 lines below FEB 2007 BVS:
            sumoc = SUM ( weight * fobs * fobs * REAL ( fcalc * CONJG(fcalc), KIND=wp ) )
            sumo2 = SUM ( weight * fobs ** 4 )
            sumc2 = SUM ( weight * REAL ( fcalc * CONJG(fcalc), KIND=wp ) ** 2 )
            sumyo = SUM ( weight * fobs * fobs )
            sumyc = SUM ( weight * REAL ( fcalc * CONJG(fcalc), KIND=wp ) )
            sumw  = SUM ( weight )
      
            abs_complex_ioic_corr_coef = ( sumw * sumoc - sumyo * sumyc ) &
                                       / SQRT ( ( sumw * sumo2 - sumyo * sumyo  ) * ( sumw * sumc2 - sumyc * sumyc ) )

        ELSE

            sumoc = SUM ( fobs * fobs * REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) )
            sumo2 = SUM ( fobs ** 4 )
            sumc2 = SUM ( REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) ** 2 )
            sumyo = SUM ( fobs * fobs )
            sumyc = SUM ( fcalc * CONJG ( fcalc ) )

            abs_complex_ioic_corr_coef = ( n * sumoc - sumyo * sumyc ) &
                                       / SQRT ( ( n * sumo2 - sumyo * sumyo ) * ( n * sumc2 - sumyc * sumyc ) )
           
        ENDIF

    END FUNCTION abs_complex_ioic_corr_coef

    FUNCTION uncentered_ioic_corr_coef ( fobs, fcalc, weight )
        REAL(KIND=wp)                                        :: uncentered_ioic_corr_coef
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight
!       Local variables:
        REAL(KIND=wp)                                        :: sumoc
        REAL(KIND=wp)                                        :: sumo2
        REAL(KIND=wp)                                        :: sumc2

!       CHeckz:
        IF ( SIZE ( fobs ) /= SIZE ( fcalc ) ) THEN
            WRITE(*,*) SIZE ( fobs), SIZE ( fcalc )
            CALL die ( 'Programming error. Size of FOBS and FCALC differ.', 'complex_complex_fofc_overlap' )
        ENDIF

        IF ( PRESENT ( weight ) ) THEN
!           Start with checks:
            IF ( SIZE ( fobs ) /= SIZE ( weight ) ) THEN
                WRITE(*,*) SIZE ( fobs), SIZE ( weight )
                CALL die ( 'Programming error. Size of FOBS and WEIGHT differ.', 'complex_complex_fofc_overlap' )
            ENDIF
            sumoc = SUM ( weight * fobs ** 2 * REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) )
            sumo2 = SUM ( weight * fobs ** 4 )
            sumc2 = SUM ( weight * REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) ** 2 )
        ELSE
            sumoc = SUM ( fobs ** 2 * REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) )
            sumo2 = SUM ( fobs ** 4 )
            sumc2 = SUM ( REAL ( fcalc * CONJG ( fcalc ), KIND=wp ) ** 2 )

        ENDIF

        uncentered_ioic_corr_coef = sumoc / SQRT ( sumo2 * sumc2 )

    END FUNCTION uncentered_ioic_corr_coef

    FUNCTION complex_complex_fofc_overlap ( fobs, fcalc, weight ) RESULT ( f )
!
!       Purpose:
!       =======
!       Calculates Fourier overlap function
!
!       Date:          Programmer:            Description of changes:
!       ====           ==========             ======================
!       Nov 2005       B.Strokopytov          Original code
!
!
        REAL(KIND=wp)                                        :: f
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight

!       Checkz:
        IF ( SIZE ( fobs ) /= SIZE ( fcalc ) ) THEN
            WRITE(*,*) SIZE ( fobs), SIZE ( fcalc )
            CALL die ( 'Programming error. Size of FOBS and FCALC differ.', 'complex_complex_fofc_overlap' )        
        ENDIF

        IF ( PRESENT ( weight ) ) THEN

!           Start with checks:
            IF ( SIZE ( fobs ) /= SIZE ( weight ) ) THEN
                WRITE(*,*) SIZE ( fobs), SIZE ( weight )
                CALL die ( 'Programming error. Size of FOBS and WEIGHT differ.', 'complex_complex_fofc_overlap' )
            ENDIF

            f = SUM ( weight * REAL ( fobs * CONJG ( fcalc ) ) )

        ELSE

            f = SUM ( REAL ( fobs * CONJG ( fcalc ) ) )

        ENDIF
    END FUNCTION complex_complex_fofc_overlap

    FUNCTION abs_complex_ioic_overlap ( fobs, fcalc_1, fcalc_2, weight ) RESULT ( f )
!
!       Purpose:
!       =======
!       Calculates patterson overlap function in reciprocal space (half sphere due to multiplicity)
!
!       f= SUM ( w * Fo**2 * Fci * CONJG ( Fcj ) )
!           h
!
!       The whole sphere will require additional factor of 2
!
        REAL(KIND=wp)                                        :: f
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN)           :: fobs
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc_1
        COMPLEX(KIND=wp), DIMENSION(:), INTENT(IN)           :: fcalc_2
        REAL(KIND=wp),    DIMENSION(:), INTENT(IN), OPTIONAL :: weight

        IF ( PRESENT ( weight ) ) THEN
            f = 2.0_wp * SUM ( ( weight * fobs ** 2 ) * REAL ( fcalc_1 * CONJG ( fcalc_2 ) ) )
        ELSE
            f = 2.0_wp * SUM ( fobs ** 2 * REAL ( fcalc_1 * CONJG ( fcalc_2 ) ) )
        ENDIF

    END FUNCTION abs_complex_ioic_overlap

    FUNCTION cc_from_matrix ( A, LT, node )
        REAL(KIND=wp)                                          :: cc_from_matrix
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: A
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: LT
        INTEGER,       DIMENSION(:),                INTENT(IN) :: node
!       Local variables:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE             :: q

!       Checkz:
        IF ( MAXVAL ( node ) > SIZE ( A, DIM=1 ) ) THEN
            WRITE(*,*) MAXVAL ( node ), SIZE ( A, DIM=1 )
            CALL die ( 'Programming error. Index of node exceeds size of matrix A.', 'cc_from_matrix' )
        ENDIF
        CALL allocate_array ( q, SIZE ( A, DIM=1 ) )
        q = 0.0_wp
        q(node) = 1.0_wp
        q = MATMUL ( LT, q )
        q = .NORM. q
        cc_from_matrix = SQRT ( 1.0_wp - DOT_PRODUCT ( q, MATMUL ( A, q ) ) )
        CALL deallocate_array ( q )

    END FUNCTION cc_from_matrix

END MODULE corr_coef_manip
