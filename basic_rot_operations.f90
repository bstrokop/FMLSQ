MODULE basic_rot_operations
USE constants
USE fail
USE vectors
USE select_kinds
USE string_manip
IMPLICIT NONE
INTEGER, PARAMETER :: accuracy = 6
INTERFACE OPERATOR ( .CCP. )
    MODULE PROCEDURE ccp4_matrix_from_angle_vector
    MODULE PROCEDURE ccp4_angle_vector_from_matrix
END INTERFACE
INTERFACE OPERATOR ( .XPLOR. )
    MODULE PROCEDURE xplor_matrix_from_angle_vector
    MODULE PROCEDURE xplor_angle_vector_from_matrix
END INTERFACE
INTERFACE OPERATOR ( .DIRCOS. )
    MODULE PROCEDURE find_dircos
    MODULE PROCEDURE find_dircos_2
END INTERFACE
PRIVATE
PUBLIC :: dircos_kappa
PUBLIC :: OPERATOR ( .CCP. )
PUBLIC :: OPERATOR ( .XPLOR. )
PUBLIC :: OPERATOR ( .DIRCOS. )
CONTAINS
      FUNCTION ccp4_matrix_from_angle_vector ( mode, vec_1 ) RESULT ( new )
!
!
!
!          Note:
!          ====
!          OMEGA/PSI
!                      gives inclination of rotation axis to K axis;
!
!          PHI
!                      gives anticlockwise rotation from I to projection
!                      of rotation axis onto I-J plane;
!
!          KAPPA
!                      is the rotation about the rotation axis.
!
!          see $CCP4/doc/rotmat.doc
!
!          Note that the rotation matrix generated from
!          (OMEGA,PHI,KAPPA) is identical to that generated from
!          (PI-OMEGA,PI+PHI,-KAPPA) so it is conventional to
!          restrict KAPPA to range: 0 to PI
!
!
          TYPE(matrix)                  :: new
          CHARACTER(len=*),  INTENT(IN) :: mode
          TYPE(vector),      INTENT(IN) :: vec_1
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: A
          REAL(KIND=wp), DIMENSION(3)   :: v
          CHARACTER(LEN=LEN(mode))              :: string
          
          string = ADJUSTL(mode)
          CALL ucase(string)
          v = vec_1        
          IF ( string == 'EULE' .OR. string == 'SPHE' ) THEN
              CALL ccprot ( A, v(1), v(2), v(3), string )
          ELSE IF ( string == 'LATT' ) THEN
!             Converting Lattman to Euler angles during ccprot call:
              CALL ccprot ( A, 0.5_wp * (v(1) + v(3)), v(2), 0.5_wp * (v(1) - v(3)), 'EULE' )
          ELSE
              WRITE(*,*) ' mode= ', mode
              CALL die('Programming error. Unknown mode.', 'ccp4_matrix_from_angle_vector')
          ENDIF
          new = A 
      END FUNCTION ccp4_matrix_from_angle_vector

      FUNCTION ccp4_angle_vector_from_matrix ( mode, mat_1 ) RESULT ( new )
          TYPE(vector)                  :: new
          CHARACTER(LEN=*),  INTENT(IN) :: mode
          TYPE(matrix),      INTENT(IN) :: mat_1
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: A
          REAL(KIND=wp), DIMENSION(3)   :: v
          CHARACTER(LEN=LEN(mode))      :: string

          string = ADJUSTL(mode)
          CALL ucase(string)
          A = mat_1
          IF ( string(1:4) == 'EULE' ) THEN
              CALL eulera(A, v(1), v(2), v(3))
          ELSE IF ( string(1:4) == 'SPHE' ) THEN
              CALL matpol(A, v(1), v(2), v(3))
          ELSE
              WRITE(*,*) ' mode=', mode
              CALL die('Programming error. Unknown mode.', 'ccp4_angle_vector_from_matrix')
          ENDIF
          new = v
      END FUNCTION ccp4_angle_vector_from_matrix

      FUNCTION xplor_matrix_from_angle_vector ( mode, vec_1 ) RESULT ( new )
          TYPE(matrix)                  :: new
          CHARACTER(LEN=*),  INTENT(IN) :: mode
          TYPE(vector),      INTENT(IN) :: vec_1
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: A
          REAL(KIND=wp), DIMENSION(3)   :: v
          REAL(KIND=wp), DIMENSION(3)   :: axis
          CHARACTER(LEN=LEN(mode))      :: string

          string = ADJUSTL(mode)
          CALL ucase(string)
          v = vec_1
          IF ( string(1:4) == 'EULE' .OR. string(1:4) == 'SPHE' .OR. string(1:4) == 'LATT' ) THEN
              CALL rotmat ( A, v(1), v(2), v(3), axis, string )
          ELSE
              WRITE(*,*) ' mode= ', mode
              CALL die('Programming error. Unknown mode.', 'ccp4_matrix_from_angle_vector')
          ENDIF
          new = A
      END FUNCTION xplor_matrix_from_angle_vector

      FUNCTION xplor_angle_vector_from_matrix ( mode, mat_1 ) RESULT ( new )
          TYPE(vector)                  :: new
          CHARACTER(LEN=*),  INTENT(IN) :: mode
          TYPE(matrix),      INTENT(IN) :: mat_1
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: A
          REAL(KIND=wp), DIMENSION(3)   :: v
          REAL(KIND=wp), DIMENSION(3)   :: axis
          CHARACTER(LEN=LEN(mode))      :: string

          string = ADJUSTL(mode)
          CALL ucase(string)
          A = mat_1
          IF ( string(1:4) == 'EULE' .OR. string(1:4) == 'SPHE' .OR. string(1:4) == 'LATT' ) THEN
              CALL matrot ( A, v(1), v(2), v(3), axis, string )
          ELSE
              CALL die ( 'Programming error. Unknown mode.', 'ccp4_angles_from_matrix' )
          ENDIF
          new = v
      END FUNCTION xplor_angle_vector_from_matrix

      FUNCTION find_dircos( mat_1 ) RESULT ( new )
          TYPE(vector)                  :: new
          TYPE(matrix), INTENT(IN)      :: mat_1
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: ROT
          REAL(KIND=wp), DIMENSION(3)   :: dircos
          REAL(KIND=wp)                 :: kappa

          ROT = mat_1
          CALL dcosfd(ROT, dircos, kappa)
          new = dircos
      END FUNCTION find_dircos

      SUBROUTINE dircos_kappa ( mat_1, vec_2, kappa)
          TYPE(matrix),  INTENT(IN)  :: mat_1
          TYPE(vector),  INTENT(OUT) :: vec_2
          REAL(KIND=wp), INTENT(OUT) :: kappa
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: ROT
          REAL(KIND=wp), DIMENSION(3)   :: dircos

          ROT = mat_1
          CALL dcosfd(ROT, dircos, kappa)
          vec_2 = dircos
      END SUBROUTINE dircos_kappa

      FUNCTION find_dircos_2 ( mat_1, mat_2 ) RESULT (angle)
          REAL(KIND=wp)                 :: angle
          TYPE(matrix), INTENT(IN)      :: mat_1
          TYPE(matrix), INTENT(IN)      :: mat_2
!         Local variables:
          REAL(KIND=wp), DIMENSION(3,3) :: ROT
          REAL(KIND=wp), DIMENSION(3)   :: dircos_1
          REAL(KIND=wp), DIMENSION(3)   :: dircos_2
          REAL(KIND=wp)                 :: kappa
          REAL(KIND=wp)                 :: cosa
          
          ROT = mat_1
          CALL dcosfd(ROT, dircos_1, kappa)
          ROT = mat_2
          CALL dcosfd(ROT, dircos_2, kappa)
          cosa  = DOT_PRODUCT( dircos_1, dircos_2 ) /&
                  SQRT ( DOT_PRODUCT(dircos_1, dircos_1 ) * DOT_PRODUCT( dircos_2, dircos_2 ) )

!         Some important precautions:
          cosa  = MAX(-1.0_wp, cosa)
          cosa  = MIN( 1.0_wp, cosa)
          angle = 180.0_wp / pi * ACOS ( cosa ) 
      END FUNCTION find_dircos_2

!     ==============================================
      SUBROUTINE CCPROT(TRMAT,ALPHA,BETA,GAMMA,MODE)
!     ==============================================
!     .. Array Arguments ..
      REAL(KIND=wp), DIMENSION(3,3), INTENT(INOUT) :: TRMAT
!
!     .. Scalar Arguments ..
      REAL(KIND=wp)                                :: ALPHA
      REAL(KIND=wp)                                :: BETA
      REAL(KIND=wp)                                :: GAMMA
      CHARACTER(LEN=4)                             :: MODE
!     ..
!     .. Local Scalars ..
      REAL(KIND=wp)                                :: AKAPPA
      REAL(KIND=wp)                                :: COSA
      REAL(KIND=wp)                                :: COSB
      REAL(KIND=wp)                                :: COSG
      REAL(KIND=wp)                                :: COSK
      REAL(KIND=wp)                                :: DEGTOR
      REAL(KIND=wp)                                :: PHI
      REAL(KIND=wp)                                :: PSI
      REAL(KIND=wp)                                :: SINA
      REAL(KIND=wp)                                :: SINB
      REAL(KIND=wp)                                :: SING
      REAL(KIND=wp)                                :: SINK
!     ..
!     .. Local Arrays ..
      REAL(KIND=wp), DIMENSION(3)                  :: DC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC                                    :: ATAN,COS,SIN
!     ..
!
      DEGTOR = pi/180.0_wp
!
      IF ( MODE=='EULE' ) THEN
        COSB = COS(BETA*DEGTOR)
        SINB = SIN(BETA*DEGTOR)
        COSA = COS(ALPHA*DEGTOR)
        SINA = SIN(ALPHA*DEGTOR)
        COSG = COS(GAMMA*DEGTOR)
        SING = SIN(GAMMA*DEGTOR)
!
        TRMAT(1,1) = COSG*COSB*COSA - SING*SINA
        TRMAT(2,1) = COSG*COSB*SINA + SING*COSA
        TRMAT(3,1) = -COSG*SINB
        TRMAT(1,2) = -SING*COSB*COSA - COSG*SINA
        TRMAT(2,2) = -SING*COSB*SINA + COSG*COSA
        TRMAT(3,2) = SING*SINB
        TRMAT(1,3) = SINB*COSA
        TRMAT(2,3) = SINB*SINA
        TRMAT(3,3) = COSB
      END IF
!
      IF ( MODE== 'SPHE' ) THEN
        PSI = ALPHA*DEGTOR
        PHI = BETA*DEGTOR
        AKAPPA = GAMMA*DEGTOR
        DC(1) = COS(PHI)*SIN(PSI)
        DC(2) = SIN(PHI)*SIN(PSI)
        DC(3) = COS(PSI)
        COSK = COS(AKAPPA)
        SINK = SIN(AKAPPA)
!
        TRMAT(1,1) = DC(1)*DC(1)* (1.0-COSK) + COSK
        TRMAT(1,2) = -DC(3)*SINK + DC(1)*DC(2)* (1.0-COSK)
        TRMAT(1,3) = +DC(2)*SINK + DC(3)*DC(1)* (1.0-COSK)
        TRMAT(2,1) = +DC(3)*SINK + DC(1)*DC(2)* (1.0-COSK)
        TRMAT(2,2) = DC(2)*DC(2)* (1.0-COSK) + COSK
        TRMAT(2,3) = -DC(1)*SINK + DC(2)*DC(3)* (1.0-COSK)
        TRMAT(3,1) = -DC(2)*SINK + DC(3)*DC(1)* (1.0-COSK)
        TRMAT(3,2) = +DC(1)*SINK + DC(2)*DC(3)* (1.0-COSK)
        TRMAT(3,3) = DC(3)*DC(3)* (1.0-COSK) + COSK
      END IF
      RETURN
!
      END SUBROUTINE ccprot
!
!     =========================================
      SUBROUTINE MATROT(ROT,T1,T2,T3,AXIS,MODE)
!     =========================================
!
! Routine computes Eulerian angles (MODE="EULE"), Lattman angles
! (MODE="LATT"), spherical polar angles (MODE="SPHE"), or rotation
! axis and angle (MODE="AXIS") from unitary matrix ROT.  The angular
! definitions are described in routine ROTMAT.
!
! Restrictions:
!    Matrix ROT has to be unitary, i.e., det(ROT)=+1 otherwise a fatal
!    warning will be issued.
!
! Conventions:
!    In Eulerian angle mode T2 can be forced without restriction of generality
!   to be located in the interval 0<= T2 <= pi (this is a consequence of
!   the identity operation 
!     t1,t2,t3 -> pi+t1,-t2,pi+t3 in Eulerian angle space).
!
!   T1 is forced into 0 <= t1 < 2*pi and T3 is forced into
!   0<=T3< 2*pi.  If T2 is 0 or PI we set T3 to 0.0 without restriction
!   of generality.
!
!   Lattman angles are forced into the intervals (0 <= theta- <= 2*pi,
!   0<= theta2 <= pi, 0 <= theta+ < 4*pi).
!
!   For spherical polar angles we force psi (=t1) into 0<= t1 < pi and
!   phi (=t2) into 0<= t2 < pi, kappa (=t3) into 0 <= t3 < 2*pi
!   without restriction of generality.  If kappa is equal to 0.0,
!    all angles will be 0.0.  If psi is equal to 0.0, phi will be set to 0.0.
!
!  Spherical polar angles are used to compute the rotation axis in AXIS mode.
!  Without restriction of generality, kappa is forced into 0 <= kappa <= pi.
!
! Input:
!    MODE specifies angle mode
!    ROT(3,3) contains the rotation matrix.  The matrix is defines
!    the rotation according to r'(i)=sum_j ROT(i,j)*r(j)
!
! Ouput:
!  T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
!  T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
!  T1,T2,T3 are psi (incl. vs. y), phi (azimuthal), kappa for MODE="SPHE"
!             and MODE="AXIS"
!  T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
! Note: all rotations are counter-clockwise
!
! Author: Axel T. Brunger
! =======================
!
!
!     .. Scalar Arguments ..
      REAL(KIND=wp)                 :: T1
      REAL(KIND=wp)                 :: T2
      REAL(KIND=wp)                 :: T3
      CHARACTER(LEN=4)              :: MODE
!     ..
!     .. Array Arguments ..
      REAL(KIND=wp), DIMENSION(3)   :: AXIS
      REAL(KIND=wp), DIMENSION(3,3) :: ROT
!     ..
!     .. Local Scalars ..
      REAL(KIND=wp)                 :: C1
      REAL(KIND=wp)                 :: C12
      REAL(KIND=wp)                 :: C2
      REAL(KIND=wp)                 :: C22
      REAL(KIND=wp)                 :: C3
      REAL(KIND=wp)                 :: DET
      REAL(KIND=wp)                 :: PII
      REAL(KIND=wp)                 :: RAD
      REAL(KIND=wp)                 :: S1
      REAL(KIND=wp)                 :: S2
      REAL(KIND=wp)                 :: S3
      REAL(KIND=wp)                 :: TM
      REAL(KIND=wp)                 :: TP
      INTEGER                       :: I
      INTEGER                       :: J
      LOGICAL                       :: COND
!     ..
!     .. Local Arrays ..
      REAL(KIND=wp), DIMENSION(3,3) :: ROT2
!     ..
!     .. External Subroutines ..
      EXTERNAL                      :: CCPERR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC                     :: ABS,ACOS,ATAN,COS,MAX,MIN,SIGN,SQRT
!     ..
!
      PII = ATAN(1.0_wp)*4.0_wp
      RAD = ATAN(1.0_wp)/45.0_wp
!
!---- test to make sure that the matrix is unitary
!
      DET = (ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3) + &
            (ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3) + &
            (ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)

      IF (ABS(DET-1.0_wp) > 0.0001_wp) THEN
!
!            *************************************
        CALL CCPERR(1,'matrix ROT is not unitary')
!            *************************************
!
      ELSE
!
        IF (MODE=='EULE' .OR. MODE=='LATT') THEN
!
!---- first, let's get t2:
!
          C2 = MAX(-1.0_wp,MIN(1.0_wp,ROT(3,3)))
          T2 = ACOS(C2)/RAD
!
!---- second, let's compute SIN(T2) (OE positive)
!
          S2 = SQRT(MAX(0.0_wp,1.0_wp-ROT(3,3)**2))
!
!---- if COS(T2) not equal +-1 then T1 and T3 are uniquely determined
!
          IF (ABS(ABS(C2)-1.0_wp) > 0.1_wp ** accuracy ) THEN  ! BVS reduced constant by a factor of 10
            S3 = ROT(1,3)/S2
            C3 = ROT(2,3)/S2
            C3 = MAX(-1.0_wp,MIN(1.0_wp,C3))
            T3 = ACOS(C3)
            IF (S3 < 0.0_wp) T3 = 2.0_wp*PII - T3
            T3 = T3/RAD
            S1 = ROT(3,1)/S2
            C1 = -ROT(3,2)/S2
            C1 = MAX(-1.0_wp,MIN(1.0_wp,C1))
            T1 = ACOS(C1)
            IF (S1 < 0.0_wp) T1 = 2.0_wp*PII - T1
            T1 = T1/RAD
          ELSE
!
!---- without restriction T3 can be set to 0.0
!
            T3 = 0.0
            C1 = ROT(1,1)
            S1 = ROT(1,2)
            C1 = MAX(-1.0_wp,MIN(1.0_wp,C1))
            T1 = ACOS(C1)
            IF (S1<0.0_wp) T1 = 2.0_wp*PII - T1
            T1 = T1/RAD
          END IF
!
          IF (MODE=='LATT') THEN
!
!---- in "Lattman" mode compute theta+, theta- from theta1, theta3
!
            TP = T1 + T3
            TM = T1 - T3
            T1 = TP
            T3 = TM
!
!---- we know that 0 <= theta1 < 2*pi and 0 <= theta3 < 2*pi.
!     We have to project this into the Lattman space
!     asymmetric unit (0<= theta+ < 4*pi,
!      0 <= theta- < 2*pi).
!
            IF (T3 < 0.0_wp) THEN
              IF (T1 < 360.0_wp) THEN
                T1 = T1 + 360.0_wp
                T3 = T3 + 360.0_wp
              ELSE
                T1 = T1 - 360.0_wp
                T3 = T3 + 360.0_wp
              END IF
            END IF
!
          END IF
!
        ELSE IF (MODE=='SPHE' .OR. MODE=='AXIS') THEN
!
!---- first lets get COS (kappa)
!
          C3 = (ROT(1,1)+ROT(2,2)+ROT(3,3)-1.0_wp)/2.0_wp
          C3 = MAX(-1.0_wp,MIN(1.0_wp,C3))
!
!---- special case COS(kappa)=1
!
          IF (ABS(C3-1.0_wp) < 0.1_wp ** accuracy ) THEN
            T1 = 0.0_wp
            T2 = 0.0_wp
            T3 = 0.0_wp
          ELSE
!
!---- determine COS^2(psi)
!
            C12 = MAX(0.0_wp, (ROT(2,2)-C3)/ (1.0_wp-C3))
!
!---- special case COS(psi)=+-1
!
            IF (ABS(C12-1.0_wp) < 0.1_wp ** accuracy ) THEN
              T1 = 0.0_wp
              T2 = 0.0_wp
              T3 = ACOS(C3)
              IF (ROT(3,1)<0.0_wp) T3 = 2.0_wp*PII - T3
              T3 = T3/RAD
            ELSE
!
!---- determine COS^2(phi)
!
              C22 = MAX(0.0_wp, (ROT(1,1)-C3)/ ((1.0_wp-C3)* (1.0_wp-C12)))
!
!---- special case COS(phi)=+-1
!
              IF (ABS(C22-1.0_wp) < 0.1_wp ** accuracy ) THEN
                C1 = SIGN(1.0_wp, (ROT(1,2)+ROT(2,1))/ (1.0_wp-C3))* &
                     SQRT(C12)
                C1 = MAX(-1.0_wp,MIN(1.0_wp,C1))
                T1 = ACOS(C1)/RAD
                T2 = 0.0_wp
                T3 = ACOS(C3)
                IF (ROT(2,3)-ROT(3,2) < 0.0_wp) T3 = 2.0_wp*PII - T3
                T3 = T3/RAD
              ELSE
!
!---- determine sign of COS(psi) and then determine psi
!
                C1 = SIGN(1.0_wp, (-ROT(2,3)-ROT(3,2))/ (1.0_wp-C3))* &
                     SQRT(C12)
                C1 = MAX(-1.0_wp,MIN(1.0_wp,C1))
                T1 = ACOS(C1)/RAD
!
!---- determine sign of COS(phi) and then determine phi
!
                C2 = SIGN(1.0_wp, (-ROT(1,3)-ROT(3,1))/ (1.0_wp-C3))* &
                     SQRT(C22)
                C2 = MAX(-1.0_wp,MIN(1.0_wp,C2))
                T2 = ACOS(C2)/RAD
!
!---- determine kappa
!
                T3 = ACOS(C3)
                IF (ABS(C2-1.0_wp) < 0.1_wp ** accuracy ) THEN
                  IF (ROT(2,3)-ROT(3,2) < 0.0_wp) T3 = 2.0_wp*PII - T3
                ELSE
                  IF (ROT(1,2)-ROT(2,1) > 0.0_wp) T3 = 2.0_wp*PII - T3
                END IF
                T3 = T3/RAD
              END IF
            END IF
          END IF
!
          IF (MODE=='AXIS') THEN
            AXIS(2) = COS(T1*RAD)
            AXIS(1) = SQRT(MAX(0.0_wp,1.0_wp-AXIS(2)**2))*COS(T2*RAD)
            AXIS(3) = -SQRT(MAX(0.0_wp,1.0_wp-AXIS(1)**2-AXIS(2)**2))
!
!---- without restriction of generality we force T3 in 0 <= T3 <= pi
!
            IF (T3 > 180.0_wp) THEN
              T3 = 2.0_wp*180.0_wp - T3
              AXIS(1) = -AXIS(1)
              AXIS(2) = -AXIS(2)
              AXIS(3) = -AXIS(3)
            END IF
          END IF
        END IF
!
!---- now back-compute the matrix as an internal consistency check
!
!            ******************************
        CALL ROTMAT(ROT2,T1,T2,T3,AXIS,MODE)
!            ******************************
!
        COND = ANY( ABS(ROT-ROT2) > 0.0001_wp )
!
        IF (COND) THEN
          WRITE (6,FMT='(A,3G12.4)') &
            ' %ROTMAT-ERR: inconsistent T1,T2,T3=',T1,T2,T3
          WRITE (6,FMT='(/A,3(/3F12.6))') ' ROT =',&
            ((ROT(I,J),J=1,3),I=1,3)
          WRITE (6,FMT='(/A,3(/3F12.6))') ' ROT2 =',&
            ((ROT2(I,J),J=1,3),I=1,3)
!              ***********************************************
          CALL CCPERR(1,'Error in internal consistency check')
!              ***********************************************
        END IF
      END IF
      END SUBROUTINE matrot
!
!     ====================================
      SUBROUTINE RORDER(SYMMAT,NAX,NORDER)
!     ====================================
!
!---- Determine rotational order of matrix SYMMAT about
!     the NAX'th axis
!
! On entry:
!  SYMMAT(3,3) symmetry matrix acting on orthogonal coordinates
!  NAX         axis to test
!
! On exit:
!  NORDER      order of axis if valid (2,3,4,6 only),
!              else = 1
!
!
!     .. Scalar Arguments ..
      INTEGER              :: NAX
      INTEGER              :: NORDER
!     ..
!     .. Array Arguments ..
      REAL(KIND=wp), DIMENSION(3,3) :: SYMMAT
!     ..
!     .. Local Scalars ..
      REAL(KIND=wp)        :: COSKAP,DET,KAPPA,ORDER,RTODEG,TRACE
      INTEGER              :: I
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC            :: ABS,ACOS,ATAN,MAX,MIN,NINT
!     ..
!
      RTODEG = 45.0_wp/ATAN(1.0_wp)
!
!---- First check that determinant = +1
!
!
      DET = (SYMMAT(1,1)*SYMMAT(2,2)-SYMMAT(1,2)*SYMMAT(2,1))*SYMMAT(3,3) + &
            (SYMMAT(2,1)*SYMMAT(3,2)-SYMMAT(2,2)*SYMMAT(3,1))*SYMMAT(1,3) + &
            (SYMMAT(3,1)*SYMMAT(1,2)-SYMMAT(3,2)*SYMMAT(1,1))*SYMMAT(2,3)

!
      IF (DET >= 0.999_wp) THEN
!
!---- Check for pure rotation (ie rotation leaves value of NAX'th
!     axis unchanged). Check NAX'th row of matrix
!
        TRACE = 0.0_wp
!
        DO 10 I = 1,3
          IF (I == NAX) THEN
            IF (ABS(SYMMAT(NAX,I)-1.0) > 0.1_wp ** accuracy ) GO TO 30
          ELSE IF (ABS(SYMMAT(NAX,I)) > 0.1_wp ** accuracy ) THEN
            GO TO 20
          END IF
          TRACE = SYMMAT(I,I) + TRACE
   10   CONTINUE
!
!---- Success, this is a rotation about NAX
!     Get rotation angle Kappa:  Trace = 1 + 2 cos Kappa
!
        COSKAP = MIN(MAX((TRACE-1.0_wp)*0.5_wp,-1.0_wp),+1.0_wp)
        KAPPA = ACOS(COSKAP)*RTODEG
!
        IF (ABS(KAPPA) < 0.01_wp) THEN
          GO TO 40
        ELSE
          ORDER = 360.0_wp/KAPPA
          NORDER = NINT(ORDER)
!
          IF (ABS(ORDER-NORDER) > 0.001_wp) THEN
!
!---- Not close enough to integral
!
            GO TO 40
!
          ELSE IF (NORDER.NE.2 .AND. NORDER.NE.3 .AND. NORDER.NE.4 .AND. &
                   NORDER.NE.6) THEN
!
!---- Not acceptable order
!
            GO TO 40
          ELSE
!
!---- Order .gt. 0 found
!
            RETURN
          END IF
        END IF
!
!---- Failed
!
   20   CONTINUE
        GO TO 40
!
!---- Failed
!
   30   CONTINUE
      END IF
!
!---- Failed to find acceptable order
!
   40 NORDER = 1
!
      END SUBROUTINE rorder
!
!     =========================================
      SUBROUTINE ROTMAT(ROT,T1,T2,T3,AXIS,MODE)
!     =========================================
!
! Routine computes unitary rotation matrix ROT using Eulerian angles
! (MODE="EULE"), Lattman angles (MODE="LATT"), spherical polar angles
! (MODE="SPHE") or a rotation about the specified axis (MODE="AXIS").
!
! Input:
!    MODE specifies angle mode
!    T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
!    T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
!    T1,T2,T3 are psi (incl. vs. y), phi (azimuthal angle, that is,
!             the angle between the x-axis and the projection of the
!             axis into the x,z plane ), kappa for MODE="SPHE"
! T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
!  Note: all rotations are counter-clockwise
!
! Output:
!    ROT(3,3) contains the rotation matrix.  Should be applied as
!    r'(i)=sum_j ROT(i,j)*r(j)
!
!    AXIS contains the normalized input AXIS vector and T1,T2,T3 will
!    contain the spherical polar angles for MODE="SPHE"
!
! Author: Axel T. Brunger
! =======================
!
!     .. Scalar Arguments ..
      REAL(KIND=wp)                 :: T1
      REAL(KIND=wp)                 :: T2
      REAL(KIND=wp)                 :: T3
      CHARACTER(LEN=4)              :: MODE
!     ..
!     .. Array Arguments ..
      REAL(KIND=wp), DIMENSION(3)   :: AXIS
      REAL(KIND=wp), DIMENSION(3,3) :: ROT
!     ..
!     .. Local Scalars ..
      REAL(KIND=wp)                 :: C1
      REAL(KIND=wp)                 :: C2
      REAL(KIND=wp)                 :: C3
      REAL(KIND=wp)                 :: C3C
      REAL(KIND=wp)                 :: DTEMP1
      REAL(KIND=wp)                 :: DTEMP2
      REAL(KIND=wp)                 :: DTEMP3
      REAL(KIND=wp)                 :: DTEMP4
      REAL(KIND=wp)                 :: NN
      REAL(KIND=wp)                 :: RAD
      REAL(KIND=wp)                 :: S1
      REAL(KIND=wp)                 :: S1SQ
      REAL(KIND=wp)                 :: S2
      REAL(KIND=wp)                 :: S3
      REAL(KIND=wp)                 :: TM
      REAL(KIND=wp)                 :: TP
      REAL(KIND=wp)                 :: TT1
      REAL(KIND=wp)                 :: TT3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC                     :: ACOS,ATAN,COS,MAX,MIN,SIN,SQRT
!     ..
!
      RAD = ATAN(1.0_wp)/45.0_wp
!
      IF (MODE=='EULE' .OR. MODE=='LATT') THEN
!
! Computes rotation matrix corresponding to the three
! Eulerian angles theta1=T1, theta2=T2, theta3=T3.  This uses
! the Rossmann and Blow convention, i.e., a theta1 rotation
! around the Z axis, followed by
! a theta2 rotation around the moved X axis, and followed by a
! theta3 rotation around the moved Z axis.  These angles are positive
! if they are anti-clockwise when looking down the relevant axis.
!
        IF (MODE=='LATT') THEN
!
!---- In "Lattman" mode, theta+=theta1+theta3=T1, theta2=T2,
!     theta-=theta1-theta3=T3
!
          TP = T1
          TM = T3
          TT1 = (TP+TM)/2.0_wp
          TT3 = (TP-TM)/2.0_wp
        ELSE
          TT1 = T1
          TT3 = T3
        END IF
!
        S1 = SIN(RAD*TT1)
        S2 = SIN(RAD*T2)
        S3 = SIN(RAD*TT3)
        C1 = COS(RAD*TT1)
        C2 = COS(RAD*T2)
        C3 = COS(RAD*TT3)
        ROT(1,1) = -S1*C2*S3 + C1*C3
        ROT(1,2) = C1*C2*S3 + S1*C3
        ROT(1,3) = S2*S3
        ROT(2,1) = -S1*C2*C3 - C1*S3
        ROT(2,2) = C1*C2*C3 - S1*S3
        ROT(2,3) = S2*C3
        ROT(3,1) = S1*S2
        ROT(3,2) = -C1*S2
        ROT(3,3) = C2
!
      ELSE IF (MODE=='SPHE' .OR. MODE=='AXIS') THEN
!
        IF (MODE=='AXIS') THEN
!
! In axis mode we obtain the psi, phi spherical polar angles from
! the AXIS vector.
! The rotation kappa corresonds to an anticlockwise rotation (t3) about
! the specified axis. The axis is desribed as a vector (AXIS).
!
          NN = SQRT(AXIS(1)**2+AXIS(2)**2+AXIS(3)**2)
!
          IF (NN < 0.0001_wp) THEN
            WRITE (6,FMT=*) 'rotation vector has 0.0 length'
            T1 = 0.0_wp
            T2 = 0.0_wp
            T3 = 0.0_wp
          ELSE
            AXIS(1) = AXIS(1)/NN
            AXIS(2) = AXIS(2)/NN
            AXIS(3) = AXIS(3)/NN
            T1 = ACOS(AXIS(2))/RAD
!
            IF (AXIS(1)**2+AXIS(3)**2 < 0.1_wp ** accuracy ) THEN
              T2 = 0.0_wp
            ELSE
!
! -- alliant change intrinsics need to match args
!
!      T2=ACOS(MAX(-1.0,MIN(1.0,AXIS(1)/SQRT(AXIS(1)**2+AXIS(3)**2))))
!     &    /RAD
!
              DTEMP1 = SQRT(AXIS(1)**2+AXIS(3)**2)
              DTEMP2 = AXIS(1)/DTEMP1
              DTEMP3 = MIN(1.0_wp,DTEMP2)
              DTEMP4 = MAX(-1.0_wp,DTEMP3)
              T2 = ACOS(DTEMP4)/RAD
!
!
              IF (AXIS(3) > 0.0_wp) T2 = -T2
            END IF
          END IF
        END IF
!
! compute rotation matrix corresponding to the three
! spherical polar angles psi=t1, phi=t2 and kappa=t3.  The rotation is
! described by specification of the direction of an axis through
! phi (azimutal angle between the x axis and the projection of the
! rotation axis on the x-z plane) and psi (inclination versus y axis).
! The angle kappa specifies the rotation around the specified axis.
!  The kappa angle is anti-clockwise when looking along the rotation axis.
! The phi angle is anti-clockwise when looking along y.
!
        S1 = SIN(RAD*T1)
        S2 = SIN(RAD*T2)
        S3 = SIN(RAD*T3)
        C1 = COS(RAD*T1)
        C2 = COS(RAD*T2)
        C3 = COS(RAD*T3)
        S1SQ = S1*S1
        C3C = 1.0_wp - C3
        ROT(1,1) = S1SQ*C2*C2*C3C + C3
        ROT(1,2) = S1*C1*C2*C3C - S1*S2*S3
        ROT(1,3) = -S1SQ*C2*S2*C3C - C1*S3
        ROT(2,1) = S1*C1*C2*C3C + S1*S2*S3
        ROT(2,2) = C1*C1*C3C + C3
        ROT(2,3) = -S1*C1*S2*C3C + S1*C2*S3
        ROT(3,1) = -S1SQ*S2*C2*C3C + C1*S3
        ROT(3,2) = -S1*C1*S2*C3C - S1*C2*S3
        ROT(3,3) = S1SQ*S2*S2*C3C + C3
!
      END IF
!
      RETURN
      END SUBROUTINE rotmat

      SUBROUTINE EULERA(TRSYM,ALPHA,BETA,GAMMA)
!     =========================================
!
! Get Euler angles ALPHA,BETA,GAMMA from matrix TRSYM
!
      REAL(KIND=wp), DIMENSION(3,3) :: TRSYM
      REAL(KIND=wp)                 :: ALPHA
      REAL(KIND=wp)                 :: BETA
      REAL(KIND=wp)                 :: GAMMA
!
      REAL(KIND=wp)                 :: RTD
      REAL(KIND=wp)                 :: R13
      REAL(KIND=wp)                 :: R32
      REAL(KIND=wp)                 :: R23
      REAL(KIND=wp)                 :: R31
      RTD = 180.0_wp/pi
!
      IF (ABS(TRSYM(3,3)) >= 1.0 - 0.1_wp ** accuracy ) GOTO 13
!
      BETA=RTD*ACOS(TRSYM(3,3))
!
      R13=TRSYM(1,3)
      R23=TRSYM(2,3)
      R31=TRSYM(3,1)
      R32=TRSYM(3,2)
      GAMMA=RTD*ATAN2(R32,-R31)
      ALPHA=RTD*ATAN2(R23,R13)
      RETURN
!
!
!  BETA = 0 OR 180  CAN ONLY FIND ALPHA + OR - GAMMA.
13    BETA=0.
      IF (TRSYM(3,3).LT.0.) BETA=180.
      GAMMA=RTD*ATAN2(TRSYM(2,1),TRSYM(2,2))
      ALPHA=0.
      END SUBROUTINE eulera

      SUBROUTINE matpol(ROT, phi, psi, kappa, dircos)
!         
!         Purpose:
!         =======
!         Given matrix ROT determines phi psi, kappa and direction cosines
!
!         Date:         Programmer:               History of changes:
!         ====          ==========                ==================
!         Oct 2005      Strokopytov B.            Original code
!         Nov 2005      Strokopytov B.            Double precision introduced
!         Nov 2005      Strokopytov B.            Internal consistency check added
!
          REAL(KIND=wp), DIMENSION(3,3), INTENT(IN)            :: ROT
          REAL(KIND=wp),                 INTENT(OUT)           :: phi
          REAL(KIND=wp),                 INTENT(OUT)           :: psi
          REAL(KIND=wp),                 INTENT(OUT)           :: kappa
          REAL(KIND=wp), DIMENSION(3),   INTENT(OUT), OPTIONAL :: dircos
!
          REAL(KIND=wp), DIMENSION(3,3)                        :: ROT2
          REAL(KIND=wp), DIMENSION(3)                          :: dircos2                                    
!         Counters:
          INTEGER                                              :: i
          INTEGER                                              :: j

          CALL dcosfd(ROT, dircos2, kappa)
          IF ( PRESENT(dircos) ) dircos = dircos2
          CALL polang(dircos2, phi, psi)

!         Internal consistency checks:
          dircos2 = (/phi, psi, kappa/)
          CALL polmat(ROT2, dircos2)
          IF ( ANY ( ABS(ROT - ROT2) > 0.0001_wp ) ) THEN
              WRITE(*,"(1X,'ROT',31X,'ROT2')")
              DO i = 1, 3
                  WRITE(*,"(1X,3F10.5, 4X, 3F10.5)") ( rot(i,j), j=1,3), ( rot2(i,j), j=1,3)
              ENDDO
              CALL ccperr(1, '*** Internal consistency check failed for matpol ***')
          ENDIF
      END SUBROUTINE matpol
!
!
!
      SUBROUTINE DCOSFD(ROT,DIRCOS,KAPPA)
!     ===================================
!******************************************************
!     Given a rotation matrix ROTN, expressing a rotation
!     through an angle KAPPA right-handedly about an axis
!     with direction cosines DIRCOS(I),
!     DCOSFD determines DIRCOS() and KAPPA
!     KAPPA is returned in degrees in the range 0 to 180
!
!     Expression for rotation matrix is (see J & J, "Mathl. Phys.",
!     p. 122):-
!
!  cw + n(1)n(1)(1-cw)    n(1)n(2)(1-cw)-n(3)sw  n(1)n(3)(1-cw)+n(3)sw
!  n(1)n(2)(1-cw)+n(3)sw  cw + n(2)n(2)(1-cw)    n(2)n(3)(1-cw)-n(1)sw
!  n(3)n(1)(1-cw)-n(2)sw  n(3)n(2)(1-cw)+n(1)sw  cw + n(3)n(3)(1-cw)
!
!     where cw = cos(KAPPA), sw = sin(KAPPA),
!     & n(1), n(2), n(3) are the direction cosines.
!
!*******************************************************}
!
      REAL(KIND=wp), DIMENSION(3,3), INTENT(IN)  :: ROT
      REAL(KIND=wp), DIMENSION(3),   INTENT(OUT) :: DIRCOS
      REAL(KIND=wp),                 INTENT(OUT) :: KAPPA
      REAL(KIND=wp), DIMENSION(3)                :: S
      REAL(KIND=wp), DIMENSION(3)                :: P
      REAL(KIND=wp), DIMENSION(3)                :: D
      REAL(KIND=wp)                              :: DIFF
      REAL(KIND=wp)                              :: RTD
      REAL(KIND=wp)                              :: TRACE
      REAL(KIND=wp)                              :: R2
      REAL(KIND=wp)                              :: R3
      REAL(KIND=wp)                              :: SINKAP
      REAL(KIND=wp)                              :: COSKAP
      REAL(KIND=wp)                              :: TERM
      INTEGER,       DIMENSION(3)                :: NS
      INTEGER,       DIMENSION(3)                :: NP
      INTEGER,       DIMENSION(3)                :: NCS
      INTEGER,       DIMENSION(3)                :: ND
      INTEGER                                    :: KNTRL
      INTEGER                                    :: J
      INTEGER                                    :: JMX

      RTD   = 180.0_wp/pi
      DIFF  = 0.1_wp ** 4
      KNTRL = 0
!     Trace = 1 + 2.cos w
      TRACE = ROT(1,1)+ROT(2,2)+ROT(3,3)
!     S(1)=2.n(1).sin w
!     S(2)=2.n(2).sin w
!     S(3)=2.n(3).sin w
      S(1)=ROT(3,2)-ROT(2,3)
      S(2)=ROT(1,3)-ROT(3,1)
      S(3)=ROT(2,1)-ROT(1,2)
      CALL ORDR3(KNTRL,S,NS)
!--Is biggest S() zero (i.e. is sin w = 0)?
      IF (ABS(S(NS(1))).GT.DIFF) THEN
        R2=S(NS(2))/S(NS(1))
        R3=S(NS(3))/S(NS(1))
        DIRCOS(NS(1))=1./SQRT(1. + R2**2 + R3**2)
        DIRCOS(NS(2))=R2*DIRCOS(NS(1))
        DIRCOS(NS(3))=R3*DIRCOS(NS(1))
        GOTO 123
      ENDIF
!--Calculation when sin w = 0 (esp. when w = 180').
!     P(1)=n(2).n(3).{0(W=0)or2(w=180)}
!     P(2)=n(3).n(1).{0(W=0)or2(w=180)}
!     P(3)=n(1).n(2).{0(W=0)or2(w=180)}
      P(1)=ROT(3,2)
      P(2)=ROT(1,3)
      P(3)=ROT(2,1)
      CALL ORDR3(KNTRL,P,NP)
!--   Is biggest P() zero (i.e. are all off-diag. terms zero)?
      IF (ABS(P(NP(1))).LT.DIFF) THEN
        IF (TRACE.GT.0.) THEN
!--Matrix is unit matrix.
          KAPPA=0.0
          DIRCOS(1)=0.
          DIRCOS(2)=0.
          DIRCOS(3)=1.
          RETURN
        ENDIF
!--Trace -ve, so dyad about x,y or z.
        KAPPA=PI
        DO J=1,3
          D(J)=ROT(J,J)+1.
        ENDDO
        CALL ORDR3(KNTRL,D,ND)
        DIRCOS(ND(1))=1.
        DIRCOS(ND(2))=0.
        DIRCOS(ND(3))=0.
        RETURN
      ENDIF
      IF (ABS(P(NP(2))).LT.DIFF) THEN
!--One d.c. is zero.
        DIRCOS(NP(1))=0.
        DIRCOS(NP(2))=SQRT(MAX(.5*(1.+ROT(NP(2),NP(2))),0.))
        DIRCOS(NP(3))=SQRT(MAX(.5*(1.+ROT(NP(3),NP(3))),0.))
        IF (P(NP(1)).LT.0.) DIRCOS(NP(2))=-DIRCOS(NP(2))
        GOTO 123
      ENDIF
!--All d.c's are nonzero.
!     R2=n(NP(1))/n(NP(2))
!     R3=n(NP(1))/n(NP(3))
      R2=P(NP(2))/P(NP(1))
      R3=P(NP(3))/P(NP(1))
      DIRCOS(NP(2))=1./SQRT(1. + R2**2 + (R2/R3)**2)
      DIRCOS(NP(3))=1./SQRT(1. + R3**2 + (R3/R2)**2)
!--Check: don't take sqrt of negative values
      TERM=1.-DIRCOS(NP(2))**2-DIRCOS(NP(3))**2
      IF (TERM.LT.0.) THEN
        IF (ABS(TERM).GT.DIFF)&
        CALL CCPERR (1,'Negative value in DCOSFD')
        TERM=0.
      ENDIF
      DIRCOS(NP(1))=SQRT(TERM)
!--Adjust signs of DIRCOS().
      JMX=0
      DO J=1,3
        IF (P(J).GT.0.) THEN
          IF (JMX.GT.0) JMX=-1
          IF (JMX.EQ.0) JMX=J
        ENDIF
      ENDDO
      IF (JMX.GT.0) DIRCOS(JMX)=-DIRCOS(JMX)
      IF (JMX.EQ.0) WRITE(6,1358)
1358  FORMAT(' All P()s. are negative or zero')
!--Given d.c's, calculate KAPPA.
123   KNTRL=1
      CALL ORDR3(KNTRL,DIRCOS,NCS)
!     Find KAPPA.
      SINKAP=ROT(NCS(3),NCS(2))-ROT(NCS(2),NCS(3))
      SINKAP=.5*SINKAP/DIRCOS(NCS(1))
      COSKAP=.5*(TRACE-1.)
      KAPPA=RTD*ATAN2(SINKAP,COSKAP)
!  If kappa negative, negate kappa and invert DIRCOS
      IF (KAPPA.LT.0.) THEN
        KAPPA = -KAPPA
        DIRCOS(1) = -DIRCOS(1)
        DIRCOS(2) = -DIRCOS(2)
        DIRCOS(3) = -DIRCOS(3)
      ENDIF
      END SUBROUTINE dcosfd
!
!
!
!
!
      SUBROUTINE ORDR3(KNTRL,Q,NQ)
!     ============================
!
!*************************************************************
!     Find NQ(1) so that |Q(NQ(1)| is biggest Q.
!     If KNTRL=1, make NQ(1),NQ(2),NQ(3) a positive permutn.
!     If KNTRL= anything else,
!     find nos. NQ() so that|Q(NQ(1)|>|Q(NQ(2))|>|Q(NQ(3))|
!
!*************************************************************
     
      REAL(KIND=wp), DIMENSION(3) :: Q
      INTEGER,       DIMENSION(3) :: NQ
      INTEGER                     :: KNTRL
      INTEGER                     :: J

!     Find no. NQ(1) of Q() with largest |Q()|.
      NQ(1)=1
      IF (ABS(Q(2)).GT.ABS(Q(1))) NQ(1)=2
      IF (ABS(Q(3)).GT.ABS(Q(NQ(1)))) NQ(1)=3
!     Find nos. for other two Q()'s & order so |Q(NQ(2))|>|Q(NQ(3))|
      NQ(2)=MOD(NQ(1),3)+1
      NQ(3)=MOD(NQ(2),3)+1
      IF (KNTRL.NE.1 .AND. ABS(Q(NQ(3))).GT.ABS(Q(NQ(2)))) THEN
        J=NQ(3)
        NQ(3)=NQ(2)
        NQ(2)=J
      ENDIF
      END SUBROUTINE ordr3
!
!
!
      SUBROUTINE POLANG(DC,OMEGA,PHI)
!     ===============================
!
! Find polar angles OMEGA & PHI (degrees) from direction cosines DC
!  OMEGA in range 0 to 180, PHI in range 0 to 360
!
      REAL(KIND=wp), DIMENSION(3) :: DC
      REAL(KIND=wp)               :: OMEGA, PHI
!
      REAL(KIND=wp)               :: RTD
!
      RTD = 180.0_wp/pi
      OMEGA=RTD*ACOS(DC(3))
      IF (DC(1).EQ.0. .AND. DC(2).EQ.0.) THEN
        PHI=0.
      ELSE
        PHI=RTD*ATAN2(DC(2),DC(1))
      ENDIF
      END SUBROUTINE polang
!
!
!
      SUBROUTINE RNGANG(A)
!     ====================
!
! Put angle A into range 0-360 degrees
!
      REAL(KIND=wp), INTENT(INOUT) :: A
!
1     IF (A.LT.0.) THEN
        A = A + 360.
        GOTO 1
      ENDIF
2     IF (A.GE.360.) THEN
        A = A - 360.
        GOTO 2
      ENDIF
      END SUBROUTINE rngang
!
!
!
      SUBROUTINE POLMAT(ROTN,ANGLES)
!     ==============================
!
!**************************************************************
!       Sets rotation matrix ROTN expressing a rotation through an
!       angle AKAPPA right-handedly about an axis
!       with polar angles OMEGA, PHI (OMEGA with Z-axis, projection
!       giving angle PHI with X-axis).
!       These angles give direction cosines DIRCOS(I).
!***************************************************************
!
      REAL(KIND=wp), DIMENSION(3,3), INTENT(OUT) :: ROTN
      REAL(KIND=wp), DIMENSION(3),   INTENT(IN)  :: ANGLES
!     Local variables:
      REAL(KIND=wp), DIMENSION(3)                :: DIRCOS
      REAL(KIND=wp)                              :: dtr
      REAL(KIND=wp)                              :: snom
      REAL(KIND=wp)                              :: csom
      REAL(KIND=wp)                              :: snph
      REAL(KIND=wp)                              :: csph
      REAL(KIND=wp)                              :: snka
      REAL(KIND=wp)                              :: cska
!     Counters:
      INTEGER                                    :: i
      INTEGER                                    :: j
      INTEGER                                    :: k
      INTEGER                                    :: k1
      REAL(KIND=wp)                              :: epsijk

      dtr = pi/180.0_wp
!
      SNOM=SIN(DTR*ANGLES(1))
      CSOM=COS(DTR*ANGLES(1))
      SNPH=SIN(DTR*ANGLES(2))
      CSPH=COS(DTR*ANGLES(2))
      SNKA=SIN(DTR*ANGLES(3))
      CSKA=COS(DTR*ANGLES(3))
      DIRCOS(1)=SNOM*CSPH
      DIRCOS(2)=SNOM*SNPH
      DIRCOS(3)=CSOM
      DO I=1,3
        DO J=1,3
          K1=6-I-J
          K=K1
          IF (K1.LT.1 .OR. K1.GT.3) K=3
          EPSIJK=((I-J)*(J-K)*(K-I))/2
          ROTN(I,J)=DIRCOS(I)*DIRCOS(J)*(1.-CSKA)-EPSIJK*DIRCOS(K)*SNKA
          IF (I.EQ.J) ROTN(I,J)=ROTN(I,J)+CSKA
        ENDDO
      ENDDO
      END SUBROUTINE polmat

END MODULE basic_rot_operations
