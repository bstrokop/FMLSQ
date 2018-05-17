MODULE estimators
USE arpack 
USE mkl95_lapack, ONLY : syevd, gttrf, gttrs, gtcon, gtrfs
USE second_derivatives
!USE implicit_matrix_inversion
!USE sparse_basic
IMPLICIT NONE
CONTAINS
    SUBROUTINE dense_1_norm ( a, est, diag )
!
!   Purpose:
!   =======
!   Estimates 1-norm of a symmetric matrix
!
!
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)            :: A
        REAL(KIND=wp),                 INTENT(OUT)           :: est
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN), OPTIONAL  :: diag
!       Allocatable arrays:
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE   :: x
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE   :: v
        INTEGER,       DIMENSION(:), ALLOCATABLE   :: isgn
!       Local variables:
        INTEGER                                    :: kase
        INTEGER                                    :: n
        CHARACTER(LEN=32)                          :: srname = 'dense_1_norm'
!       Counters:
        INTEGER                                    :: itns

        n = SIZE ( A, DIM=1 )
        IF ( n /= SIZE ( A, DIM=2 ) ) THEN
            WRITE(*,*) n, SIZE ( A, DIM=2 )
            CALL die ('Programming error. Matrix A is not a square matrix.', srname)
        ENDIF

!       Allocate temporary arrays:
        CALL allocate_array(x, n)
        CALL allocate_array(v, n)
        CALL allocate_array(isgn, n)

!       Initialize:
        kase = 0

!       Use reverse communication interface:
        itns = 0
        est  = 1.0_wp
        DO WHILE ( .TRUE. )

            CALL sonest ( n, v, x, isgn, est, kase )

!           For symmetric matrices only:
            IF ( kase == 1 .OR. kase == 2 ) THEN
                WRITE(*,"(' DENSE_1-NORM> ', 'iteration=', I2, ' current estimate=', ES9.2)") itns, est

!               Overwrite X with Ax matrix-vector product:
                IF ( .NOT. PRESENT ( diag ) ) THEN

!                   Matrix is assumed to be symmetric here:
                    x = MATMUL ( A, x )
                ELSE
                    x = MATMUL ( A, x * diag )
                    x = diag * x
                ENDIF

                itns = itns + 1
            
            ELSE IF ( kase == 0 ) THEN

!               Algorithm has converged:
                EXIT
            ENDIF
        ENDDO

!       Free memory:
        CALL deallocate_array (x)
        CALL deallocate_array (v)
        CALL deallocate_array (isgn)  

        WRITE(*,"(' DENSE_1_NORM> ', 'Estimated 1-Norm=', ES9.2)") est
    END SUBROUTINE dense_1_norm

    FUNCTION norm_1(np, mtz_1, maps, pdb_2, mode, sk)
!
!       Purpose:
!       =======
!       Calculates 1-norm of normal matrix
!       Could be used as a very rough largest lambda approximation.
!       
! 
!       Note:
!       ====
!       We deal with symmetrical positive definite matrix here only.
!
!       Date:              Programmer:              History of changes:
!       ====               ==========               ==================
!       Mar 2008           B.Strokopytov            Original code
!
        REAL(KIND=wp)                                                :: norm_1
        INTEGER,                              INTENT(IN)             :: np
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(pdb),                            INTENT(IN)             :: pdb_2
        CHARACTER(LEN=*),                     INTENT(IN)             :: mode
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN),   OPTIONAL :: sk
!       Local automatic arrays:
        REAL(KIND=wp), DIMENSION(np)                                 :: x
        REAL(KIND=wp), DIMENSION(np)                                 :: y 
        REAL(KIND=wp), DIMENSION(np)                                 :: v
        INTEGER,       DIMENSION(np)                                 :: isgn
!       Local variables:
        INTEGER                                                      :: kase
        REAL(KIND=wp)                                                :: est
!       Counters:
        INTEGER                                                      :: itns
        CHARACTER(LEN=32)                                            :: srname = 'norm_1'

!       Initialize:
        kase = 0

!       Use reverse communication interface:
        itns = 0
        est  = 1.0_wp
        DO WHILE ( .TRUE. )

            CALL sonest(np, v, x, isgn, est, kase)

!           For symmetric matrices only:
            IF ( kase == 1 .OR. kase == 2 ) THEN
                WRITE(*,"(' NORM_1> ', 'iteration=', I2, ' current estimate=', ES9.2)") itns, est

                IF ( .NOT. PRESENT ( sk ) ) THEN
                    CALL Ax(np, y, x, mtz_1, maps, pdb_2, mode )
                ELSE
                    CALL Ax(np, y, x, mtz_1, maps, pdb_2, mode, sk )
                ENDIF
!               Overwrite x:
                x = y
 
                itns = itns + 1

            ELSE IF ( kase == 0 ) THEN

!               Algorithm has converged:
                EXIT

            ENDIF

        ENDDO

!       Final result:
        norm_1 = est

    END FUNCTION norm_1

    FUNCTION inv_norm_1(np, mtz_1, maps, pdb_2, mode, sk)
!
!       Purpose:
!       =======
!       Calculates 1-norm of normal matrix
!       Could be used as a very rough largest lambda approximation.
!       
! 
!       Note:
!       ====
!       We deal with symmetrical positive definite matrix here only.
!
!       Date:              Programmer:              History of changes:
!       ====               ==========               ==================
!       Mar 2008           B.Strokopytov            Original code
!
        REAL(KIND=wp)                                                :: inv_norm_1
        INTEGER,                              INTENT(IN)             :: np
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(pdb),                            INTENT(IN)             :: pdb_2
        CHARACTER(LEN=*),                     INTENT(IN)             :: mode
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN),   OPTIONAL :: sk
!       Local automatic arrays:
        REAL(KIND=wp), DIMENSION(np)                                 :: x
        REAL(KIND=wp), DIMENSION(np)                                 :: y 
        REAL(KIND=wp), DIMENSION(np)                                 :: v
        INTEGER,       DIMENSION(np)                                 :: isgn
!       SYMMLQ:
        REAL(KIND=wp)                                                :: shift
        LOGICAL                                                      :: checkA
        LOGICAL                                                      :: precon
        INTEGER                                                      :: itnlim
        INTEGER                                                      :: nout
        REAL(KIND=wp)                                                :: rtol
        INTEGER                                                      :: istop
        INTEGER                                                      :: itn
        REAL(KIND=wp)                                                :: Anorm
        REAL(KIND=wp)                                                :: Acond
        REAL(KIND=wp)                                                :: rnorm 
        REAL(KIND=wp)                                                :: ynorm
!       Local variables:
        INTEGER                                                      :: kase
        REAL(KIND=wp)                                                :: est
!       Counters:
        INTEGER                                                      :: itns
        CHARACTER(LEN=32)                                            :: srname = 'inv_norm_1'

!       Initialise SYMMLQ:
        checkA = .TRUE.

!       Preconditioning will used but will be applied to matrix directly (see subroutine SYMMLQ):
        precon = .FALSE.
 
        itnlim = 2 * np
        nout   = 0 
        rtol   = 0.1_wp ** 3

!       Initialize:
        kase = 0

!       Use reverse communication interface:
        itns = 0
        est  = 1.0_wp

        DO WHILE ( .TRUE. )

            CALL sonest(np, v, x, isgn, est, kase)

!           For symmetric matrices only:
            IF ( kase == 1 .OR. kase == 2 ) THEN
                WRITE(*,"('INV_NORM_1> ', 'iteration=', I2, ' current estimate=', ES9.2)") itns, est

                IF ( PRESENT ( sk ) ) THEN
!                    x = x * sk
                    CALL symmlq ( np, Ax, Diagonal_Msolve, x, shift, checkA, precon,              &
                                  y,  itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm, &
                                  mtz_1, maps, pdb_2, mode, sk)
!                    y = y / sk
                ELSE
                    CALL warn('Inaccurate programming...', srname)
                    CALL symmlq ( np, Ax, Diagonal_Msolve, x, shift, checkA, precon,              &
                                  y,  itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm, &
                                  mtz_1, maps, pdb_2, mode, sk)

                ENDIF

!               Overwrite x:
                x = y
 
                itns = itns + 1

            ELSE IF ( kase == 0 ) THEN

!               Algorithm has converged:
                EXIT

            ENDIF

        ENDDO

!       Final result:
        inv_norm_1 = 1.0_wp / est

    END FUNCTION inv_norm_1

    SUBROUTINE sonest (n, v, x, isgn, est, kase)
 
!     Code converted using TO_F90 by Alan Miller
!     Date: 2008-02-26  Time: 14:04:15
!
!     Note:
!     ====
!     Need to get rid of those GOTOs.

        INTEGER,       INTENT(IN)    :: n
        REAL(KIND=wp), INTENT(INOUT) :: v(n)
        REAL(KIND=wp), INTENT(INOUT) :: x(n)
        INTEGER,       INTENT(INOUT) :: isgn(n)
        REAL(KIND=wp), INTENT(INOUT) :: est
        INTEGER,       INTENT(INOUT) :: kase
!
!     SONEST ESTIMATES THE 1-NORM OF A SQUARE, REAL MATRIX  A.
!     REVERSE COMMUNICATION IS USED FOR EVALUATING
!     MATRIX-VECTOR PRODUCTS.

!     ON ENTRY

!        N       INTEGER
!                THE ORDER OF THE MATRIX.  N .GE. 1.

!        ISGN    INTEGER(N)
!                USED AS WORKSPACE.

!        KASE    INTEGER
!                = 0.

!     ON INTERMEDIATE RETURNS

!        KASE    = 1 OR 2.

!        X       REAL(N)
!                MUST BE OVERWRITTEN BY

!                     A*X,             IF KASE=1,
!                     TRANSPOSE(A)*X,  IF KASE=2,

!                AND SONEST MUST BE RE-CALLED, WITH ALL THE OTHER
!                PARAMETERS UNCHANGED.

!     ON FINAL RETURN

!        KASE    = 0.

!        EST     REAL
!                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).

!        V       REAL(N)
!                = A*W,   WHERE  EST = NORM(V)/NORM(W)
!                         (W  IS NOT RETURNED).

!     THIS VERSION DATED MARCH 16, 1988.
!     NICK HIGHAM, UNIVERSITY OF MANCHESTER.

!     REFERENCE
!     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
!     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
!     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
!     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.

!     SUBROUTINES AND FUNCTIONS
!     BLAS     ISAMAX, SASUM, SCOPY
!     GENERIC  ABS, NINT, REAL, SIGN

        INTEGER,       PARAMETER :: itmax = 5
        REAL(KIND=wp), PARAMETER :: zero = 0.0_wp
        REAL(KIND=wp), PARAMETER :: one = 1.0_wp
        REAL(KIND=wp), PARAMETER :: two = 2.0_wp

!       INTERNAL VARIABLES:
        INTEGER       :: i, iter, j, jlast, jump
        REAL(KIND=wp) :: altsgn, estold, temp
!       Save all between calls:
        SAVE

        IF ( kase == 0 ) THEN
            x = one / REAL ( n, KIND=wp )
            kase = 1
            jump = 1
            RETURN
        END IF

        SELECT CASE ( jump )
            CASE (1)
            GO TO 100
            CASE (2)
            GO TO  200
            CASE (3)
            GO TO  300
            CASE (4)
            GO TO  400
            CASE(5)
            GO TO  500
        END SELECT

!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

    100 CONTINUE
        IF ( n == 1 ) THEN
            v(1) = x(1)
            est = ABS ( v(1) )
!           ... QUIT
            GO TO 510
        END IF

!       est = sasum(n,x,1)
        est = SUM ( ABS ( x ) )

        DO i = 1,n
            x(i)    = SIGN (one, x(i) )
            isgn(i) = NINT ( x(i) )
        ENDDO
        kase = 2
        jump = 2
        RETURN

!       ................ ENTRY   (JUMP = 2)
!       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

    200 CONTINUE

!       j = isamax(n,x,1)
        j = MAXLOC(x, DIM=1)
        iter = 2

!       MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

    220 CONTINUE
        DO i = 1,n
            x(i) = zero
        ENDDO
        x(j) = one
        kase = 1
        jump = 3
        RETURN

!       ................ ENTRY   (JUMP = 3)
!       X HAS BEEN OVERWRITTEN BY A*X.

    300 CONTINUE
!       CALL scopy(n,x,1,v,1)
        v = x
        estold = est
!       est = sasum(n,v,1)
        est = SUM ( ABS ( v ) )
        DO i = 1,n
            IF ( NINT ( SIGN(one,x(i)) ) /= isgn(i) ) GO TO 320
        ENDDO

!       REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED:
        GO TO 410

    320 CONTINUE

!       TEST FOR CYCLING.
        IF (est <= estold) GO TO 410

        DO i = 1,n
            x(i) = SIGN ( one, x(i) )
            isgn(i) = NINT ( x(i) )
        ENDDO
        kase = 2
        jump = 4
        RETURN

!       ................ ENTRY   (JUMP = 4)
!       X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

    400 CONTINUE
        jlast = j
!       j = isamax(n,x,1)
        j = MAXLOC ( x, DIM=1 )
        IF ( (  x(jlast) /= ABS ( x(j) )  ) .AND. (iter < itmax)   ) THEN
            iter = iter + 1
            GO TO 220
        ENDIF

!       ITERATION COMPLETE.  FINAL STAGE.

    410 CONTINUE
        altsgn = one
        DO i = 1,n
            x(i) = altsgn * (one + REAL ( i-1, KIND=wp ) / REAL ( n-1, KIND=wp ))
            altsgn = -altsgn
        ENDDO
        kase = 1
        jump = 5
        RETURN

!       ................ ENTRY   (JUMP = 5)
!       X HAS BEEN OVERWRITTEN BY A*X.

    500 CONTINUE
!       temp = two*sasum(n,x,1)/REAL(3*n)
        temp = two * SUM ( ABS ( x ) ) / REAL ( 3 * n, KIND=wp )
        IF ( temp > est ) THEN
!           CALL scopy(n,x,1,v,1)
            v = x
            est = temp
        ENDIF

    510 kase = 0

    END SUBROUTINE sonest

    SUBROUTINE diag_dense_matrix_estimate ( ntrials, A, diag_estimate )
!
!       Purpose:
!       =======
!       Estimate diagonal elements (and trace) of squared matrix.
!       
!       Discussion:
!       ==========
!       Addition of linear solver instead of MATMUL will allow estimation
!       of Tr(f(A)), e.g. Tr(A^-1) !!
!
!       This is major achievement of Bekas et al. (not mine).
!
        INTEGER, INTENT(IN)                          :: ntrials
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)    :: A
        REAL(KIND=wp), DIMENSION(:),   INTENT(INOUT) :: diag_estimate
!       Local arrays:
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE   :: av
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE   :: v
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE   :: tk
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE   :: qk
!       Local variables:
        INTEGER                                      :: n ! matrix size
!       Counters:
        INTEGER                                      :: k
        CHARACTER(LEN=32)                            :: srname = 'diag_dense_matrix_estimate'

        n = SIZE ( A, DIM=1 )
        CALL init_random_number()
         
        CALL allocate_array(v, n)
        CALL allocate_array(av,n)
        CALL allocate_array(qk,n)
        CALL allocate_array(tk,n)

        diag_estimate = 0.0_wp
        tk = 0.0_wp
        qk = 0.0_wp
        DO k = 1, ntrials
            CALL RANDOM_NUMBER(v)
!           Create random number distribution between -1 and +1:
            v = 2.0_wp * v - 1.0_wp

            av = MATMUL(A,v)
            tk = tk + av * v 
            qk = qk + v * v
        ENDDO
        IF ( ALL ( qk > 0.0_wp ) ) THEN
            diag_estimate = tk / qk
            WRITE(*,"(' bekas diag=', 15ES9.2)") diag_estimate
        ELSE
            CALL die('Programming error. Zero components are not allowed in QK vector.',&
                      srname)
        ENDIF

        CALL deallocate_array(av)
        CALL deallocate_array(v)
        CALL deallocate_array(tk)
        CALL deallocate_array(qk)

    END SUBROUTINE diag_dense_matrix_estimate

    SUBROUTINE hadamard_diag_estimate(np, mtz_1, maps, pdb_2, mode, sk, diag_estimate)
!
!       Purpose:
!       =======
!       Estimate diagonal elements (and trace) of squared matrix.
!       
!       Discussion:
!       ==========
!       Addition of linear solver instead of MATMUL will allow estimation
!       of Tr(f(A)), e.g. Tr(A^-1) !!
!
!       This is major achievement of Bekas et al. (not mine).
!
!       Note:
!       ====
!       Before calling this routine need to calculate gradients to prepare phases, etc.
!
        INTEGER,                                   INTENT(IN)            :: np
        TYPE(mtz),                                 INTENT(INOUT)         :: mtz_1
        TYPE(map),     DIMENSION(:),  ALLOCATABLE, INTENT(INOUT)         :: maps
        TYPE(pdb),                                 INTENT(INOUT)         :: pdb_2
        CHARACTER(LEN=*),                          INTENT(IN)            :: mode
        REAL(KIND=wp), DIMENSION(:),               INTENT(IN),  OPTIONAL :: sk
        REAL(KIND=wp), DIMENSION(:),               INTENT(INOUT)         :: diag_estimate
!       Local arrays:
        REAL(KIND=wp), DIMENSION(np)                                     :: av
        REAL(KIND=wp), DIMENSION(np)                                     :: v
        REAL(KIND=wp), DIMENSION(np)                                     :: tk
        REAL(KIND=wp), DIMENSION(np)                                     :: qk
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE                       :: H
!       Local variables:
        INTEGER                                                          :: m 
        INTEGER                                                          :: n
        INTEGER                                                          :: ntrials
        INTEGER                                                          :: non_positive
!       Counters:
        INTEGER                                                          :: k
        CHARACTER(LEN=32)                                                :: srname = 'hadamard_diag_estimate'

!       Initialize arrays:         
        diag_estimate = 0.0_wp
        tk = 0.0_wp
        qk = 0.0_wp

        m = 128 
        WRITE(*,"(' DIAG_IMPLICIT_MATRIX_ESTIMATE> ', 'Just ', A, &
        &' matrix-vector multiplications needed to esimate main diagonal.')") &
        TRIM ( int_to_c ( m ) )
!       Figure out n:
        n = 1
        DO k = 1, 30
            n = 2 * n

!           Found suitable n:
            IF ( n > np ) EXIT 
        ENDDO

        WRITE(*,*) k, n, ' for Hadamard matrix'
        CALL allocate_array ( H, m, n)
        CALL hadamard(m, n, H)
        WRITE(*,*) ' Hadamard matrix has been generated. Occupying ', m*n*8, ' bytes.'

        DO ntrials = 1, m

            v(1:np) =  H(ntrials,1:np)

!           Scaling is crucial to get there faster and get more accurate result:
            IF ( PRESENT ( sk ) ) THEN
                v(1:np) =  sk(1:np) * v(1:np)
            ENDIF
            WRITE(*,*) ' positive entries=', COUNT(v > 0), ' negative=', COUNT (v < 0)

!           Cycle if too many unbalances signs:
            IF ( ABS ( COUNT ( v >= 0 ) -  COUNT ( v < 0 ) ) > 4 ) CYCLE
!           Implicit matrix multiplication:
            CALL Ax(np, av, v, mtz_1, maps, pdb_2, mode )
            tk = tk + av * v 
            qk = qk + v * v
            diag_estimate = tk / qk
            non_positive = COUNT ( diag_estimate < 0.0_wp )
            WRITE(*,"('Number of matrix-vector multiplications=', I6, ' trace=', ES11.4, ' non-positive=',I4)") &
            ntrials, SUM ( diag_estimate ), non_positive

        ENDDO

        IF ( ALLOCATED ( H ) ) CALL deallocate_array(H)

    END SUBROUTINE hadamard_diag_estimate

    SUBROUTINE random_diag_estimate(np, mtz_1, maps, pdb_2, mode, sk, diag_estimate)
!
!       Purpose:
!       =======
!       Estimate diagonal elements (and trace) of squared matrix.
!       
!       Discussion:
!       ==========
!       Addition of linear solver instead of MATMUL will allow estimation
!       of Tr(f(A)), e.g. Tr(A^-1) !!
!
!       This is major achievement of Bekas et al. (not mine).
!
!       Note:
!       ====
!       Before calling this routine need to calculate gradients to prepare phases, etc.
!
        INTEGER,                                   INTENT(IN)            :: np
        TYPE(mtz),                                 INTENT(INOUT)         :: mtz_1
        TYPE(map),     DIMENSION(:),  ALLOCATABLE, INTENT(INOUT)         :: maps
        TYPE(pdb),                                 INTENT(INOUT)         :: pdb_2
        CHARACTER(LEN=*),                          INTENT(IN)            :: mode
        REAL(KIND=wp), DIMENSION(:),               INTENT(IN),  OPTIONAL :: sk
        REAL(KIND=wp), DIMENSION(:),               INTENT(INOUT)         :: diag_estimate
!       Local automatic arrays:
        REAL(KIND=wp), DIMENSION(np)                                     :: av
        REAL(KIND=wp), DIMENSION(np)                                     :: v
        REAL(KIND=wp), DIMENSION(np)                                     :: tk
        REAL(KIND=wp), DIMENSION(np)                                     :: qk
!       Local variables:
        INTEGER                                                          :: m 
        INTEGER                                                          :: ntrials
        INTEGER                                                          :: non_positive
!       Counters:
        INTEGER                                                          :: k
        CHARACTER(LEN=32)                                                :: srname = 'random_diag_estimate'

        CALL init_random_number()
!
!       Initialize arrays:         
        diag_estimate = 0.0_wp
        tk = 0.0_wp
        qk = 0.0_wp

!       Going to carry out M matirx-vector multiplications:
        m = 64 
        CALL messag('Need '// TRIM ( int_to_c ( m ) ) //' matrix-vector multiplications to estimate main diagonal.&
                   & Wait...', srname)
        DO ntrials = 1, m

            CALL RANDOM_NUMBER(v)

!           Create random number distribution between -1 and +1:
            v = v + v - 1.0_wp

!           Scaling is VERY important:
            IF ( PRESENT ( sk ) ) THEN
                v(1:np) =  sk(1:np) * v(1:np)
            ENDIF

!           Implicit matrix multiplication:
            CALL Ax(np, av, v, mtz_1, maps, pdb_2, mode )
            tk = tk + av * v 
            qk = qk + v * v
            diag_estimate = tk / qk

!           Check how many non-positive elements we have now:
            non_positive = COUNT ( diag_estimate < 0.0_wp )
            WRITE(*,"(' RANDOM_DIAG_ESTIMATE> ', 'Number of matrix-vector multiplications=', I6, &
           &          ' trace=', ES11.4, ' non-positive=',I4)") &
            ntrials, SUM ( diag_estimate ), non_positive

        ENDDO

    END SUBROUTINE random_diag_estimate

    SUBROUTINE hadamard ( m, n, a )

!*****************************************************************************80
!
!! HADAMARD returns the HADAMARD matrix.
!
!  Definition:
!
!    A Hadamard matrix is a square matrix A of order N, whose entries are
!    only +1's or -1's, with the property that:
!
!      A * A' = N * I.
!
!  Notes:
!
!    A Hadamard matrix must be of order 1, 2, or else a multiple of 4.
!    It is not known whether a Hadamard matrix exists for every multiple
!    of 4.
!
!    The method used here allows the user to request a Hadamard matrix
!    of any rectangular order, M by N.  The algorithm then essentially
!    finds the largest powers of 2 that are less than or equal to M and
!    N, and produces a Hadamard-like matrix in that space, setting the
!    rest of the matrix to 0.  Thus, the matrix returned by this routine
!    is only a Hadamard matrix if M = N = a power of 2.
!
!  Formula:
!
!    The following recursive formula is used to produce a series of
!    Hadamard matrices of increasing size.
!
!    H(0) = [1]
!
!    H(1) = [ H(0)  H(0) ] = [ 1  1]
!           [ H(0) -H(0) ]   [ 1 -1]
!
!    H(2) = [ H(1)  H(1) ] = [ 1  1  1  1]
!           [ H(1) -H(1) ]   [ 1 -1  1 -1]
!                            [ 1  1 -1 -1]
!                            [ 1 -1 -1  1]
!
!    and so on.
!
!  Properties:
!
!    All entries of a Hadamard matrix are either +1 or -1.  Matrices
!    produced by this routine will be +1 or -1 up to a certain row
!    and column, beyond which the entries will be zero.
!
!    The Hadamard matrices produced by this routine have the property
!    that the first row and column are entirely 1's, although this
!    is not a requirement for a Hadamard matrix.
!
!    The matrices produced by this algorithm are (loosely) symmetric,
!    although that is not required for a Hadamard matrix.
!
!    Hadamard matrices exhibit the maximum possible relative growth of pivot
!    elements during Gaussian elimination with complete pivoting.
!
!    The inverse of a Hadamard matrix of order N is A itself,
!    scaled by 1.0/N.  Thus 1.0/sqrt(N) times a Hadamard matrix
!    yields a symmetric matrix which is its own inverse, or
!    "involutional".
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!  Modified:
!
!    28 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory and David Karney,
!    Example 3.14,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969, page 42, QA263 G862.
!
!    William Pratt,
!    Digital Image Processing,
!    John Wiley and Sons, 1978.
!
!    Herbert Ryser,
!    Combinatorial Mathematics,
!    John Wiley and Sons, 1963.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        real    ( kind = wp ) a(m,n)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) nn

        if ( m <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'HADAMARD - Fatal error!'
            write ( *, '(a,i8)' ) '  Input value of M = ', m
            write ( *, '(a)' ) '  but M must be positive.'
            stop
        end if

        if ( n <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'HADAMARD - Fatal error!'
            write ( *, '(a,i8)' ) '  Input value of N = ', n
            write ( *, '(a)' ) '  but N must be positive.'
            stop
        end if

        a(1,1) = 1.0D+00

        nn = 1

        do while ( nn < n .or. nn < m )

             do i = 1, nn
                do j = 1, nn

                    if ( i <= m .and. j + nn <= n ) then
                        if ( 2 * nn <= n ) then
                            a(i,j+nn) = a(i,j)
                        else
                            a(i,j+nn) = 0.0D+00
                        end if
                    end if

                    if ( i + nn <= m .and. j <= n ) then
                        if ( 2 * nn <= m ) then
                            a(i+nn,j) = a(i,j)
                        else
                            a(i+nn,j) = 0.0D+00
                        end if
                    end if

                    if ( i + nn <= m .and. j + nn <= n ) then
                        if ( 2 * nn <= m .and. 2 * nn <= n ) then
                            a(i+nn,j+nn) = - a(i,j)
                        else
                            a(i+nn,j+nn) = 0.0D+00
                        end if
                    end if

                end do
            end do

            nn = 2 * nn

        end do

    END SUBROUTINE hadamard

    SUBROUTINE golub_dense(np, A, shift, u, bounds, nout, sk, xpower)
        INTEGER,                       INTENT(IN)             :: np
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)             :: A
!       Lowest and largest lambda:
        REAL(KIND=wp), DIMENSION(2),   INTENT(IN)             :: shift
!       Vector to esimate u'f(A)u product:
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN)             :: u
        REAL(KIND=wp), DIMENSION(2),   INTENT(OUT)            :: bounds
        INTEGER,                       INTENT(IN)             :: nout
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN),   OPTIONAL :: sk                  
!       Normally -1 to estimate inverse matrix element or A^-1 product:
        INTEGER,                       INTENT(IN),   OPTIONAL :: xpower
!       Local arrays:
        REAL,    PARAMETER                                    :: eps = SQRT ( EPSILON ( 1.0_wp ) )
        INTEGER, PARAMETER                                    :: maxit = 50
        REAL(KIND=wp), DIMENSION(maxit+1,maxit+1)             :: Tj
        REAL(KIND=wp), DIMENSION(np)                          :: x 
        REAL(KIND=wp), DIMENSION(np)                          :: xold
        REAL(KIND=wp), DIMENSION(np)                          :: rj
        REAL(KIND=wp), DIMENSION(maxit+1)                     :: ej
        REAL(KIND=wp), DIMENSION(maxit+1)                     :: theta
        REAL(KIND=wp), DIMENSION(maxit+1)                     :: xtemp 
!       TT matrix to avoid Intel compiler warnings:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE            :: TT
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE            :: first_ev_components
!       Local variables:
        REAL(KIND=wp)                                         :: alpha
        REAL(KIND=wp)                                         :: gamma
        REAL(KIND=wp), DIMENSION(2)                           :: psi
        REAL(KIND=wp), DIMENSION(2)                           :: Ij
        REAL(KIND=wp), DIMENSION(2)                           :: bounds_old
        REAL(KIND=wp)                                         :: unorm2
        INTEGER                                               :: power
        REAL(KIND=wp)                                         :: error
        REAL(KIND=wp)                                         :: p_error
!       Counters:
        INTEGER                                               :: i 
        INTEGER                                               :: j
        INTEGER                                               :: j1
        INTEGER                                               :: k
        INTEGER                                               :: m 
        CHARACTER(LEN=32)                                     :: srname = 'golub_dense'

        unorm2 = DOT_PRODUCT ( u, u )

        IF ( debug > 40 ) THEN
            WRITE(*,*) ' unorm2=', unorm2
        ENDIF

!       This corresponds to element estimation of inverse matrix (default):
        power = -1  

!       What we gonna do:
        IF ( PRESENT ( xpower ) ) THEN        
            power = xpower
        ENDIF

!       Initialise arrays:
        x = u / SQRT ( unorm2 ) 
        xold  = 0.0_wp
        gamma = 0.0_wp

!       Initialise Tj matrix:
        Tj = 0.0_wp

!       Some pretty printing:
        IF ( nout > 0 ) THEN
            CALL messag(' ', srname)
            CALL messag('   itns    Lower Bound   Upper Bound  Difference       error', srname)
        ENDIF

        outer_loop:DO j = 1, maxit
            j1    = j + 1

!           sk is needed to simulate diagonal conditioning:            
            IF ( PRESENT ( sk )  ) THEN
                rj    = MATMUL ( A, sk * x )
                rj    = sk * rj
            ELSE
                rj =  MATMUL ( A, x )
            ENDIF

            alpha = DOT_PRODUCT ( x, rj )
            rj = rj - alpha * x - gamma * xold
            gamma = SQRT ( DOT_PRODUCT ( rj, rj ) )
            Tj(j,j)    = alpha

!           Prepare ej, only j-th component will be non-zero:
            ej      = 0.0_wp
            ej(j:j) = gamma ** 2

!           Lower and upper bound PSI:
            gauss_radau_rule:DO k = 1, 2
                IF ( j == 1 ) THEN
                    xtemp(j:j) = gamma ** 2 / ( Tj(j,j) - shift(k) )
                ELSE
                    error = shift(k)

!                   Need some copying to avoid Inter compiler warnings:
                    CALL allocate_array ( TT, j, j )
                    TT = Tj(1:j,1:j)

                    CALL dense_linear_solver(j, xtemp(1:j), ej(1:j), error, TT)

!                   Free temporary array:
                    CALL deallocate_array ( TT )
                ENDIF

                psi(k:k) = shift(k:k) + xtemp(j:j)

!               Construct Tj_prime:
                Tj(j,j1)  = gamma
                Tj(j1,j)  = gamma
                Tj(j1,j1) = psi(k)

!               Debugging:
                IF ( debug > 40 ) THEN
                    DO m = 1, j1
                        WRITE(*,"(' Tj: ', 100ES9.2)") Tj(m,1:j1)
                    ENDDO
                ENDIF

!               Calculate eigenvalues theta plus first components of eigenvectors (Tj'): 
                CALL allocate_array ( TT, j1, j1 )
                TT = Tj(1:j1,1:j1)
                CALL ev_comp(j1, TT, theta(1:j1), first_ev_components)
                CALL deallocate_array(TT)

!               Debugging:
                IF ( debug > 40 ) THEN
                    WRITE(*,"(' eigenvalues=                     ', 100ES9.2)") theta(1:j1)
                    WRITE(*,"(' first components of eigenvectors=', 100ES9.2)") first_ev_components(1:j1)
                ENDIF

!               f(theta_k) = theta_k^power:
                Ij(k) = SUM ( first_ev_components(1:j1) ** 2 * theta(1:j1) ** power )

            ENDDO gauss_radau_rule

!           We want bounds(2) contain upper bound:
            bounds(2) = Ij(1)

!           BOUNDS(1) will contain lower bound u'f(A)u:
            bounds(1) = Ij(2)

!           Calculate renormalized final bounds for u'f(A)u product:
            bounds(1:2) = unorm2 * bounds(1:2) !* DOT_PRODUCT ( sk * u, sk * u ) 

            IF ( nout > 0 ) THEN
                error   =  ABS ( bounds(1) - bounds(2) )
                p_error = (100.0_wp * error) /  ABS ( bounds(1) )
                IF ( p_error < 100.0_wp ) THEN
                    WRITE(*,"(' GOLUB_DENSE> ', 3X, I4, 2X, 3(2X,ES11.4), F12.2,'%')") &
                    j, bounds(1:2), error, p_error
                ELSE
                    WRITE(*,"(' GOLUB_DENSE> ', 3X, I4, 2X, 3(2X,ES11.4), 5X,' >100.0%')") &
                    j, bounds(1:2), error
                ENDIF
            ENDIF

            IF ( j > 1 ) THEN
                IF ( ALL ( ABS ( bounds - bounds_old ) < eps * bounds ) ) THEN
                    EXIT
                ENDIF
            ENDIF

            xold = x
            x = rj / gamma

!           Copy to old bounds:
            bounds_old = bounds

        ENDDO outer_loop

!       In rare cases we will have this problem because of precision loss:
        IF ( bounds(1) > bounds(2) ) THEN
            error = bounds(1)
            bounds(1) = bounds(2)
            bounds(2) = error
            IF ( debug > 100 ) THEN
                CALL messag('Ooops... Swap has been made...', srname)
            ENDIF
        ENDIF

        IF ( nout > 0 ) THEN
            CALL messag(' ', srname)
            WRITE(*,"(' GOLUB_DENSE> ', 'Final lower bound=', ES13.6, ' Final Upper bound=', ES13.6)")&
            bounds(1:2) !* DOT_PRODUCT ( sk, sk )  
            WRITE(*,*) ' unorm2=', unorm2
            CALL messag('Done.', srname)
        ENDIF


!       Free memory:
        IF ( ALLOCATED ( TT ) ) CALL deallocate_array (TT)
        IF ( ALLOCATED ( first_ev_components ) ) CALL deallocate_array (first_ev_components)
            
    END SUBROUTINE golub_dense

    SUBROUTINE dense_linear_solver ( n, x, b, shift, A )
        INTEGER,                       INTENT(IN)    :: n
        REAL(KIND=wp), DIMENSION(n),   INTENT(INOUT) :: x
        REAL(KIND=wp), DIMENSION(n),   INTENT(IN)    :: b
        REAL(KIND=wp),                 INTENT(IN)    :: shift
        REAL(KIND=wp), DIMENSION(n,n), INTENT(IN)    :: A
!       Local variables:
        LOGICAL                                      :: precon
        LOGICAL                                      :: checkA
        INTEGER                                      :: itnlim
        INTEGER                                      :: nout
        REAL(KIND=wp)                                :: rtol
        INTEGER                                      :: istop
        INTEGER                                      :: itn
        REAL(KIND=wp)                                :: Anorm
        REAL(KIND=wp)                                :: Acond
        REAL(KIND=wp)                                :: rnorm
        REAL(KIND=wp)                                :: ynorm
        REAL(KIND=wp), DIMENSION(n)                  :: y 
        REAL(KIND=wp), DIMENSION(n)                  :: sk ! dummy SK array 
!       Test:
        REAL(KIND=wp), DIMENSION(n)                  :: r1
!       LAPACK:
        REAL(KIND=wp), EXTERNAL                      :: dnrm2
        REAL(KIND=wp)                                :: r1norm
!       Counters:
        INTEGER                                      :: i
        INTEGER                                      :: j 
        CHARACTER(LEN=32)                            :: srname = 'dense_linear_solver'

!       Use silent mode (no printing):
        nout = 0 
        checkA = .TRUE.
        precon = .TRUE.
        itnlim = 5 * n
        rtol = EPSILON ( 1.0_wp )
!
!        WRITE(*,*) '  Parameter     Diagonal element     Shift     Diagonal - shift   RHS b' 
        DO i = 1, n
            sk(i) = ABS ( A(i,i) - shift )
            WRITE(*,"(' Tj:', I4, 500F10.5)") i, (A(i,j),j=1,n)
!           SK(i) MUST BE POSITIVE:
!!            IF ( A(i,i) - shift <= 0.0_wp ) THEN
!                WRITE(*,"(1X, I8, 4(5X, ES9.2))") i,  A(i,i), shift, A(i,i) - shift, b(i)
!!!                CALL die('Possible programming error. Negative diagonal element.', srname)
!!            ENDIF
        ENDDO

!       Solve Mx = b:
        CALL symmlq_dense(n, Aprod_dense, Diagonal_Msolve, b, shift, checkA, precon,     &
                          x, itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm, &
                          sk,  A)

        IF ( itn > 3 * n .OR. istop /= 1) THEN
            WRITE(*,*) ' istop=', istop, ' itns made in linear solver=', itn
            nout = 6
        ENDIF

        IF ( nout > 0 ) THEN
            IF ( shift /= 0.0_wp ) THEN
                WRITE(*,*) ' shift=', shift
            ENDIF

!           Check solution:
            y(1:n) = MATMUL ( A(1:n,1:n), x(1:n) ) - shift * x(1:n)

            r1  = b - y
            r1norm = SQRT ( DOT_PRODUCT ( r1, r1 ) )

!           Report quality for first 10 eqns:
            DO i = 1, MIN ( 10, n )
                WRITE(*,"(' b=', ES12.5, ' y=', ES12.5, ' diff= ', ES9.2)") &
                            b(i), y(i), r1(i)
            ENDDO

!           Print residual:
            WRITE(nout, 2000) r1norm
    2000    FORMAT(/ ' Final residual =', 1p, e8.1)

!           Print solution:
            WRITE(nout, 2100) (i, x(i), i = 1, n)
    2100    FORMAT(/ ' Solution  x' / 1p, 4(i6, e14.6))
        ENDIF

    END SUBROUTINE dense_linear_solver

    SUBROUTINE ev_comp(n, A, lambdas, first_ev_components )
        INTEGER,                                  INTENT(IN)    :: n
        REAL(KIND=wp), DIMENSION(:,:),            INTENT(IN)    :: A
        REAL(KIND=wp), DIMENSION(:),              INTENT(INOUT) :: lambdas
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: first_ev_components
!       Local arrays:
        REAL(KIND=wp), DIMENSION(n,n)                           :: eigenvectors

        IF ( .NOT. ALLOCATED ( first_ev_components ) ) THEN
            CALL allocate_array ( first_ev_components, n)
        ELSE
            CALL deallocate_array ( first_ev_components )
            CALL allocate_array ( first_ev_components, n)
        ENDIF

        CALL dense_eigen(n, A, lambdas, eigenvectors)
        first_ev_components(1:n) = eigenvectors(1,1:n)        

    END SUBROUTINE ev_comp

    SUBROUTINE dense_eigen ( n, A, lambdas, eigen_vectors )
        INTEGER,                       INTENT(IN)    :: n
        REAL(KIND=wp), DIMENSION(n,n), INTENT(IN)    :: A
        REAL(KIND=wp), DIMENSION(n),   INTENT(INOUT) :: lambdas
        REAL(KIND=wp), DIMENSION(n,n), INTENT(INOUT) :: eigen_vectors
!       Local arrays:
        INTEGER                                      :: info
        CHARACTER(LEN=32)                            :: srname='dense_eigen'

!       Preserve matrix A:        
        eigen_vectors = A

!       Figure out eigenvectors:
        CALL syevd(eigen_vectors, lambdas, jobz='V', UPLO='U', info=info)

        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('Possible programming error. INFO is not zero on return form SYEVD',&
                      srname)
        ENDIF

    END SUBROUTINE dense_eigen

    SUBROUTINE dense_frobenius_norm ( n, A, maxit, trace, bounds, hutchinson, sk, xp )
        INTEGER,                       INTENT(IN)             :: n
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)             :: A
        INTEGER,                       INTENT(IN)             :: maxit
        REAL(KIND=wp),                 INTENT(OUT)            :: trace
        REAL(KIND=wp), DIMENSION(2),   INTENT(INOUT)          :: bounds
        LOGICAL,                       INTENT(IN)             :: hutchinson
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN)             :: sk
        REAL(KIND=wp),                 INTENT(IN),   OPTIONAL :: xp

        IF ( PRESENT ( xp ) ) THEN
            CALL dense_trace_estimation ( n, A, maxit, trace, bounds, hutchinson, sk, xp, xpower=2 )
        ELSE
            CALL dense_trace_estimation ( n, A, maxit, trace, bounds, hutchinson, sk, xpower=2 )
        ENDIF

    END SUBROUTINE dense_frobenius_norm

    SUBROUTINE dense_trace_estimation ( n, A, maxit, trace, bounds, hutchinson, sk, xp, xpower )
        INTEGER,                       INTENT(IN)             :: n
        REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)             :: A
        INTEGER,                       INTENT(IN)             :: maxit
        REAL(KIND=wp),                 INTENT(OUT)            :: trace
        REAL(KIND=wp), DIMENSION(2),   INTENT(INOUT)          :: bounds
        LOGICAL,                       INTENT(IN)             :: hutchinson
        REAL(KIND=wp), DIMENSION(:),   INTENT(IN),   OPTIONAL :: sk
        REAL(KIND=wp),                 INTENT(IN),   OPTIONAL :: xp
        INTEGER,                       INTENT(IN),   OPTIONAL :: xpower
!       Local arrays:
        REAL(KIND=wp), DIMENSION(n)                           :: z
!       This depends on MAXIT:
        REAL(KIND=wp), DIMENSION(maxit)                       :: L
        REAL(KIND=wp), DIMENSION(maxit)                       :: U
        REAL(KIND=wp), DIMENSION(2)                           :: extreme_lambdas
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE            :: H
!       Local variables:
        INTEGER                                               :: power
        REAL(KIND=wp)                                         :: p
        REAL(KIND=wp)                                         :: Lmin
        REAL(KIND=wp)                                         :: Umax
        REAL(KIND=wp)                                         :: Lp 
        REAL(KIND=wp)                                         :: Up
        REAL(KIND=wp)                                         :: nu
        REAL(KIND=wp)                                         :: Ij
!       Counters:
        INTEGER                                               :: i
        INTEGER                                               :: j
        INTEGER                                               :: positive
        INTEGER                                               :: total_sign
        CHARACTER(LEN=32)                                     :: srname='dense_trace_estimation'

        CALL messag(' ', srname)

!       Set default probability level:
        p = 0.95_wp
        IF ( PRESENT ( xp ) ) THEN
            p = xp
            IF ( p < 0.75_wp ) THEN
                CALL warn('Too low probablity level.', srname)
            ENDIF
        ENDIF
        
        power = -1
        IF ( PRESENT ( xpower ) ) THEN
            power = xpower
        ENDIF

!       Cannot spend to much time on estimating lambdas accurately:
        CALL dense_1_norm ( A, extreme_lambdas(2) )
        extreme_lambdas(1) = 0.1_wp ** 10 * extreme_lambdas(2)

        IF ( Hutchinson ) THEN
            CALL init_random_number()
        ELSE
            i = 1
            DO j = 1, 30
                i = 2 * i
                IF ( i > n ) EXIT
            ENDDO
            CALL allocate_array ( H, maxit, i)
            CALL hadamard(maxit, i, H)
        ENDIF

        CALL messag(' ', srname)
        CALL messag('      itns      trace(est.)', srname)


!       Random number quality check:
        total_sign = 0
        positive   = 0

        DO j = 1, maxit

            IF ( HUTCHINSON ) THEN

!               Using Hutchinson method:
                CALL RANDOM_NUMBER(z)
                DO i = 1, n
                    IF ( z(i) < 0.5_wp ) THEN
                        z(i) = -1.0_wp 
                    ELSE
                        z(i) = 1.0_wp
                    ENDIF
                ENDDO
            ELSE

!               Use Hadamard row vectors:
                z(1:n) = H(j,1:n)

                IF ( debug > 40 ) THEN
                    IF ( n < 20 ) THEN
                        WRITE(*,"(' z=', 100F4.1)") z(1:n)
                    ENDIF
                ENDIF
            ENDIF

            total_sign = total_sign + n
            positive = positive + COUNT ( z(1:n) > 0 )

!           Estimate z'Tr(A**-1)z:
            IF ( PRESENT ( sk ) ) THEN
                z = sk * z
                CALL golub_dense(n, A, extreme_lambdas, z, bounds, 0, sk, xpower=power)
            ELSE
                CALL golub_dense(n, A, extreme_lambdas, z, bounds, 0, xpower=power)
            ENDIF

            L(j) = bounds(1)
            U(j) = bounds(2)

            IF ( L(j) <= 0.0_wp .OR. U(j) <= 0.0_wp ) THEN
                WRITE(*,*) L(j), U(j)
                CALL die('Matrix is not positive definite.', srname)
            ENDIF

            Ij = ( 0.5_wp / j ) * SUM ( L(1:j) + U(1:j) )
            WRITE(*,"(' DENSE_TRACE_ESTIMATION> ', 5X, I5, 5X, ES11.4)")  j, Ij

            Lmin = MIN ( Lmin, L(j) )
            Umax = MAX ( Umax, U(j) )
                
            nu = SQRT ( -0.5_wp * j * ( Umax - Lmin ) ** 2 * LOG ( ( 1.0_wp - p ) / 2.0_wp ) )

            Lp = ( 1.0_wp / j ) * SUM ( L(1:j) ) - nu / j
            Up = ( 1.0_wp / j ) * SUM ( U(1:j) ) + nu / j

        ENDDO
 
        trace = Ij

        bounds(1) = Lp
        bounds(2) = Up

!       Report quality of random number generation:
        CALL messag(' ',srname)

        WRITE(*,"(' DENSE_TRACE_ESTIMATION> ', 'estimated trace=', ES11.4, /,                                    &
       &          ' DENSE_TRACE_ESTIMATION> ' , I2, '%-probability lower and upper bounds for trace=', 2ES9.2)") &
        trace, NINT ( 100 * p ), bounds
        WRITE(*,"(' DENSE_TRACE_ESTIMATION> ', 'RG quality: +/- sign=', 2I8)") &
        positive, total_sign - positive

        IF ( ALLOCATED ( H ) ) CALL deallocate_array(H)
        
        CALL messag(' ',srname)
    END SUBROUTINE dense_trace_estimation

    SUBROUTINE golub_implicit(np, mtz_1, maps, pdb_2, mode, shift, u, bounds, nout, sk, QINV, O, itn)
        INTEGER,                              INTENT(IN)             :: np
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(pdb),                            INTENT(IN)             :: pdb_2
        CHARACTER(LEN=*),                     INTENT(IN)             :: mode
        REAL(KIND=wp), DIMENSION(2),          INTENT(IN)             :: shift
!       Vector to esimate u'f(A)u product:
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN)             :: u
        REAL(KIND=wp), DIMENSION(2),          INTENT(OUT)            :: bounds
        INTEGER,                              INTENT(IN)             :: nout
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN),   OPTIONAL :: sk
        TYPE(csr_matrix),                     INTENT(IN)             :: QINV
        REAL(KIND=wp), DIMENSION(:,:),        INTENT(INOUT),OPTIONAL :: O
!       Normally -1 to estimate inverse matrix element or A^-1 product:
        INTEGER,                              INTENT(OUT),  OPTIONAL :: itn
!       Local arrays:
        REAL,    PARAMETER                                           :: eps = 0.1_wp ** 5 
        INTEGER, PARAMETER                                           :: maxit = 500
        INTEGER, PARAMETER                                           :: enough_itns = 2
!        REAL(KIND=wp), DIMENSION(maxit+1,maxit+1)                    :: Tj
!       Tridiagonal matrix:
        REAL(KIND=wp), DIMENSION(maxit)                              :: d
        REAL(KIND=wp), DIMENSION(maxit)                              :: dl
        REAL(KIND=wp), DIMENSION(maxit)                              :: du
!       Additional arrays for Lanczos:
        REAL(KIND=wp), DIMENSION(np)                                 :: x 
        REAL(KIND=wp), DIMENSION(np)                                 :: g 
        REAL(KIND=wp), DIMENSION(np)                                 :: rj
        REAL(KIND=wp), DIMENSION(maxit+1)                            :: ej
        REAL(KIND=wp)                                                :: sigma_j
!       Local variables:
        REAL(KIND=wp)                                                :: alpha
        REAL(KIND=wp)                                                :: gamma
        REAL(KIND=wp), DIMENSION(2)                                  :: psi
        REAL(KIND=wp), DIMENSION(2)                                  :: Ij
        REAL(KIND=wp), DIMENSION(2)                                  :: bounds_old
        REAL(KIND=wp)                                                :: unorm2
        INTEGER                                                      :: power
        REAL(KIND=wp)                                                :: error
        REAL(KIND=wp)                                                :: p_error
!       Counters:
        INTEGER                                                      :: enough
        INTEGER                                                      :: i 
        INTEGER                                                      :: j
        INTEGER                                                      :: j1
        INTEGER                                                      :: k
        INTEGER                                                      :: m 
        CHARACTER(LEN=32)                                            :: srname = 'golub_implicit'

        unorm2 = DOT_PRODUCT ( u, u )
        IF ( debug > 1040 ) THEN
            WRITE(*,*) ' unorm2=', unorm2
        ENDIF

!       This corresponds to element estimation of inverse matrix (default):
        power = -1  

!       Initialise arrays:
        x = u / SQRT ( unorm2 )
        IF ( ALL ( x == 0.0_wp ) ) THEN
            CALL die('Programming error. Zero norm vector on output', srname)
        ENDIF

        gamma = 0.0_wp
        rj = 0.0_wp

!       Initialise Tj matrix:
!        Tj = 0.0_wp
        dl = 0.0_wp
        du = 0.0_wp
        d = 0.0_wp

!       Some pretty printing:
        IF ( nout > 0 ) THEN
            WRITE(*,*)
            WRITE(*,*) '  itns     Lower Bound  Upper Bound  Difference        error'
        ENDIF

        enough = 0
        outer_loop:DO j = 1, maxit
            j1    = j + 1
            IF ( PRESENT ( O ) ) THEN
                O(1:np,j) = x
            ENDIF

!           Implicit matrix-vector multiplication:
!            IF ( .NOT. PRESENT ( sk ) ) THEN
                CALL cond_Ax(np, g, x, mtz_1, maps, pdb_2, mode, sk=sk, Q=QINV)
!            ELSE
!                CALL Ax(np, rj, x, mtz_1, maps, pdb_2, mode, sk)
!            ENDIF


            alpha = DOT_PRODUCT ( x, g )

!           Get Ax on the first iteration or  Ax - gamma * xold otherwise:
            rj = rj + g 
!            alpha = DOT_PRODUCT ( x, rj )

!           Get rj = Ax - gamma * x - alpha x:
            rj = rj - alpha * x

!           Full orthogonalization stage:
            IF ( PRESENT ( O ) ) THEN
                g = rj
                IF ( j > 1 .AND. j <= SIZE (O, DIM=2) ) THEN
                    DO m = 1, j 
                        rj = rj - DOT_PRODUCT ( g, O(1:np,m) ) * O(1:np,m)
                    ENDDO
                ENDIF
            ENDIF

            gamma = SQRT ( DOT_PRODUCT ( rj, rj ) )

!           Build main diagonal:
            d(j)    = alpha

!           Have to do clean up all the time because of tridiagonal solver which uses ej:
            ej(1:j) = 0.0_wp            
            ej(j:j) = gamma ** 2

!           Lower and upper bound PSI:
            gauss_radau_rule:DO k = 1, 2
                IF ( j == 1 ) THEN
                    sigma_j = gamma ** 2 / ( d(j) - shift(k) )
                ELSE
                    error = shift(k)

!                   Need some copying to avoid Inter compiler warnings:
!                    CALL allocate_array ( TT, j, j )
!                    TT = Tj(1:j,1:j)
                    CALL general_tridiagonal_solver(j, dl, d, du, ej, sigma_j, shift(k) )
!                    CALL dense_linear_solver(j, xtemp(1:j), ej(1:j), error, TT)

!                   Free temporary array:
!                    CALL deallocate_array ( TT )
                ENDIF

!               Shift(k) = a or b, xtemp(j) = sigma_j:
                psi(k) = shift(k) + sigma_j

!               Construct Tj_prime:
!                Tj(j,j1)  = gamma
!                Tj(j1,j)  = gamma
!                Tj(j1,j1) = psi(k)
                dl(j) = gamma
                du(j) = gamma
                d(j1) = psi(k)

!               Calculate eigenvalues theta plus first components of eigenvectors (Tj'):
                ej(1)    = 1.0_wp
                ej(2:j1) = 0.0_wp
                CALL general_tridiagonal_solver(j1, dl, d, du, ej, sigma_j)

!               Calculated TT(1,1)**-1 as suggested by Golub:
                Ij(k) = sigma_j

            ENDDO gauss_radau_rule

!           We want bounds(2) contain upper bound:
            bounds(2) = Ij(1)

!           BOUNDS(1) will contain lower bound u'f(A)u:
            bounds(1) = Ij(2)

!           Calculate renormalized final bounds for u'f(A)u product:
            bounds(1:2) = unorm2 * bounds(1:2)

            IF ( nout > 0 ) THEN
                error   =  ABS ( bounds(1) - bounds(2) )
                p_error = (100.0_wp * error) /  ABS ( bounds(1) )
                IF ( j < 10 .OR. MOD ( j, 10 ) == 0 .OR. &
                     ALL ( ABS ( bounds - bounds_old ) < 10 * eps * bounds ) ) THEN
                    IF ( p_error < 100.0_wp ) THEN
                        WRITE(*,"(3X, I4, 2X, 3(2X,ES11.4), F12.2,'%')") &
                          j, bounds(1:2), error, p_error
                    ELSE
                        WRITE(*,"(3X, I4, 2X, 3(2X,ES11.4), 5X,' >100.0%')") &
                          j, bounds(1:2), error
                    ENDIF
                    IF ( MOD (j,10) == 0 ) WRITE(*,*)
                ENDIF
            ENDIF

            IF ( j > 1 ) THEN
                IF ( ALL ( ABS ( bounds - bounds_old ) < eps * bounds .AND. bounds(1) < bounds(2) ) ) THEN
!                    enough = enough + 1
!                    IF ( enough == enough_itns ) EXIT
                    EXIT
                ENDIF
            ENDIF

!           Save x:
            g = x

            IF ( gamma == 0.0_wp ) EXIT

!           Generate new x:
            x = rj / gamma

!           Use old x:
            rj = -gamma * g  

!           Copy to old bounds:
            bounds_old = bounds

        ENDDO outer_loop

!       In rare cases we will have this problem because of precision loss:
        IF ( bounds(1) > bounds(2) ) THEN
            error = bounds(1)
            bounds(1) = bounds(2)
            bounds(2) = error
            IF ( debug > 100 ) THEN
                CALL messag('Ooops... Swap has been made...', srname)
            ENDIF
        ENDIF

        IF ( nout > 0 ) THEN
            WRITE(*,"(/,1X,84('+'))")
            WRITE(*,"(53X,'Final             Final')")
            WRITE(*,"(53X,'Lower bound       Upper bound')")
            WRITE(*,"(48X, 2(4X, ES14.7))") bounds(1:2)
            WRITE(*,"(1X,84('+'),/)")
            CALL messag('Done.', srname)
        ENDIF

        IF ( PRESENT ( itn ) ) itn = j

    END SUBROUTINE golub_implicit

    SUBROUTINE golub_implicit_anypower(np, mtz_1, maps, pdb_2, mode, shift, u, bounds, nout, sk, QINV, O, xpower, itn)
        INTEGER,                              INTENT(IN)             :: np
        TYPE(mtz),                            INTENT(INOUT)          :: mtz_1
        TYPE(map), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: maps
        TYPE(pdb),                            INTENT(IN)             :: pdb_2
        CHARACTER(LEN=*),                     INTENT(IN)             :: mode
        REAL(KIND=wp), DIMENSION(2),          INTENT(IN)             :: shift
!       Vector to esimate u'f(A)u product:
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN)             :: u
        REAL(KIND=wp), DIMENSION(2),          INTENT(OUT)            :: bounds
        INTEGER,                              INTENT(IN)             :: nout
        REAL(KIND=wp), DIMENSION(:),          INTENT(IN),   OPTIONAL :: sk
        TYPE(csr_matrix),                     INTENT(IN)             :: QINV
        REAL(KIND=wp), DIMENSION(:,:),        INTENT(INOUT),OPTIONAL :: O
!       Normally -1 to estimate inverse matrix element or A^-1 product:
        INTEGER,                              INTENT(IN),   OPTIONAL :: xpower
        INTEGER,                              INTENT(OUT),  OPTIONAL :: itn
!       Local arrays:
        REAL,    PARAMETER                                           :: eps = 0.1_wp ** 5 
        INTEGER, PARAMETER                                           :: maxit = 500
        INTEGER, PARAMETER                                           :: enough_itns = 2
        REAL(KIND=wp), DIMENSION(maxit+1,maxit+1)                    :: Tj
        REAL(KIND=wp), DIMENSION(np)                                 :: x 
        REAL(KIND=wp), DIMENSION(np)                                 :: g 
        REAL(KIND=wp), DIMENSION(np)                                 :: rj
        REAL(KIND=wp), DIMENSION(maxit+1)                            :: ej
        REAL(KIND=wp), DIMENSION(maxit+1)                            :: theta
        REAL(KIND=wp), DIMENSION(maxit+1)                            :: xtemp 
!       TT matrix to avoid Intel compiler warnings:
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE                   :: TT
        REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE                   :: first_ev_components
!       Local variables:
        REAL(KIND=wp)                                                :: alpha
        REAL(KIND=wp)                                                :: gamma
        REAL(KIND=wp), DIMENSION(2)                                  :: psi
        REAL(KIND=wp), DIMENSION(2)                                  :: Ij
        REAL(KIND=wp), DIMENSION(2)                                  :: bounds_old
        REAL(KIND=wp)                                                :: unorm2
        INTEGER                                                      :: power
        REAL(KIND=wp)                                                :: error
        REAL(KIND=wp)                                                :: p_error
!       Counters:
        INTEGER                                                      :: enough
        INTEGER                                                      :: i 
        INTEGER                                                      :: j
        INTEGER                                                      :: j1
        INTEGER                                                      :: k
        INTEGER                                                      :: m 
        CHARACTER(LEN=32)                                            :: srname = 'golub_implicit_anypower'

        unorm2 = DOT_PRODUCT ( u, u )
        IF ( debug > 40 ) THEN
            WRITE(*,*) ' unorm2=', unorm2
        ENDIF

!       This corresponds to element estimation of inverse matrix (default):
        power = -1  

!       What we gonna do:
        IF ( PRESENT ( xpower ) ) THEN        
            power = xpower
        ENDIF

!       Initialise arrays:
        x = u / SQRT ( unorm2 )
        IF ( ALL ( x == 0.0_wp ) ) THEN
            CALL die('Programming error. Zero norm vector on output', srname)
        ENDIF
        gamma = 0.0_wp
        rj = 0.0_wp

!       Initialise Tj matrix:
        Tj = 0.0_wp

!       Some pretty printing:
        IF ( nout > 0 ) THEN
            CALL messag(' ', srname)
            WRITE(*,*)
            WRITE(*,*) '   itns    Lower Bound   Upper Bound  Difference       error'
        ENDIF

        enough = 0
        outer_loop:DO j = 1, maxit
            j1    = j + 1
            IF ( PRESENT ( O ) ) THEN
                O(1:np,j) = x
            ENDIF

!           Implicit matrix-vector multiplication:
!            IF ( .NOT. PRESENT ( sk ) ) THEN
                CALL cond_Ax(np, g, x, mtz_1, maps, pdb_2, mode, sk=sk, Q=QINV)
!            ELSE
!                CALL Ax(np, rj, x, mtz_1, maps, pdb_2, mode, sk)
!            ENDIF


            alpha = DOT_PRODUCT ( x, g )

!           Get Ax on the first iteration or  Ax - gamma * xold otherwise:
            rj = rj + g 
!            alpha = DOT_PRODUCT ( x, rj )

!           Get rj = Ax - gamma * x - alpha x:
            rj = rj - alpha * x

!           Full orthogonalization stage:
            IF ( PRESENT ( O ) ) THEN
                g = rj
                IF ( j > 1 .AND. j <= SIZE (O, DIM=2) ) THEN
                    DO m = 1, j 
                        rj = rj - DOT_PRODUCT ( g, O(1:np,m) ) * O(1:np,m)
                    ENDDO
                ENDIF
            ENDIF

            gamma = SQRT ( DOT_PRODUCT ( rj, rj ) )
            Tj(j,j)    = alpha

!           Prepare ej, only j-th component will be non-zero:
            ej      = 0.0_wp
            ej(j:j) = gamma ** 2

!           Lower and upper bound PSI:
            gauss_radau_rule:DO k = 1, 2
                IF ( j == 1 ) THEN
                    xtemp(1) = gamma ** 2 / ( Tj(1,1) - shift(k) )
                ELSE
                    error = shift(k)

!                   Need some copying to avoid Inter compiler warnings:
                    CALL allocate_array ( TT, j, j )
                    TT = Tj(1:j,1:j)

                    CALL dense_linear_solver(j, xtemp(1:j), ej(1:j), error, TT)

!                   Free temporary array:
                    CALL deallocate_array ( TT )
                ENDIF

!               Shift(k) = a or b, xtemp(j) = sigma_j:
                psi(k) = shift(k) + xtemp(j)

!               Construct Tj_prime:
                Tj(j,j1)  = gamma
                Tj(j1,j)  = gamma
                Tj(j1,j1) = psi(k)

!               Debugging:
                IF ( debug > 40 ) THEN
                    DO m = 1, j1
                        WRITE(*,"(' Tj: ', 100ES9.2)") Tj(m,1:j1)
                    ENDDO
                ENDIF

!               Calculate eigenvalues theta plus first components of eigenvectors (Tj'): 
                CALL allocate_array ( TT, j1, j1 )
                TT = Tj(1:j1,1:j1)
                CALL ev_comp(j1, TT, theta(1:j1), first_ev_components)
                IF ( ANY ( theta(1:j1) < 0.0_wp ) ) THEN
                    WRITE(*,"(' theta(1:j1)=', 500ES11.3)") theta(1:j1)
                    WRITE(*,*) ' xtemp(j)=sigma_j=', xtemp(j)
                    DO m = 1, j1
                        WRITE(*,"(' Tj:', I4, 500F10.5)") i, TT(m,1:j1)
                    ENDDO
                    CALL die('O-o-o-o-p-s-s. Negative lambdas detected...', srname)
                ENDIF
                CALL deallocate_array(TT)

!               Debugging:
                IF ( debug > 40 ) THEN
                    WRITE(*,"(' eigenvalues=                     ', 100ES9.2)") theta(1:j1)
                    WRITE(*,"(' first components of eigenvectors=', 100ES9.2)") first_ev_components(1:j1)
                ENDIF

!               f(theta_k) = theta_k^power:
                Ij(k) = SUM ( first_ev_components(1:j1) ** 2 * theta(1:j1) ** power )

            ENDDO gauss_radau_rule

!           We want bounds(2) contain upper bound:
            bounds(2) = Ij(1)

!           BOUNDS(1) will contain lower bound u'f(A)u:
            bounds(1) = Ij(2)

!           Calculate renormalized final bounds for u'f(A)u product:
            bounds(1:2) = unorm2 * bounds(1:2)

            IF ( nout > 0 ) THEN
                error   =  ABS ( bounds(1) - bounds(2) )
                p_error = (100.0_wp * error) /  ABS ( bounds(1) )
                IF ( j < 10 .OR. MOD ( j, 10 ) == 0 .OR. &
                     ALL ( ABS ( bounds - bounds_old ) < 10 * eps * bounds ) ) THEN
                    IF ( p_error < 100.0_wp ) THEN
                        WRITE(*,"(3X, I4, 2X, 3(2X,ES11.4), F12.2,'%')") &
                          j, bounds(1:2), error, p_error
                    ELSE
                        WRITE(*,"(3X, I4, 2X, 3(2X,ES11.4), 5X,' >100.0%')") &
                          j, bounds(1:2), error
                    ENDIF
                    IF ( MOD (j,10) == 0 ) WRITE(*,*)
                ENDIF
            ENDIF

            IF ( j > 1 ) THEN
                IF ( ALL ( ABS ( bounds - bounds_old ) < eps * bounds .AND. bounds(1) < bounds(2) ) ) THEN
!                    enough = enough + 1
!                    IF ( enough == enough_itns ) EXIT
                    EXIT
                ENDIF
            ENDIF

!           Save x:
            g = x

            IF ( gamma == 0.0_wp ) EXIT

!           Generate new x:
            x = rj / gamma

!           Use old x:
            rj = -gamma * g  

!           Copy to old bounds:
            bounds_old = bounds

        ENDDO outer_loop

!       In rare cases we will have this problem because of precision loss:
        IF ( bounds(1) > bounds(2) ) THEN
            error = bounds(1)
            bounds(1) = bounds(2)
            bounds(2) = error
            IF ( debug > 100 ) THEN
                CALL messag('Ooops... Swap has been made...', srname)
            ENDIF
        ENDIF

        IF ( nout > 0 ) THEN
            CALL messag(' ', srname)
            WRITE(*,"(' GOLUB_IMPLICIT> ', 'Final lower bound=', ES13.6, ' Final Upper bound=', ES13.6)")&
            bounds(1:2)
            CALL messag('Done.', srname)
        ENDIF

        IF ( PRESENT ( itn ) ) itn = j

!       Free memory:
        IF ( ALLOCATED ( TT ) )                  CALL deallocate_array (TT)
        IF ( ALLOCATED ( first_ev_components ) ) CALL deallocate_array (first_ev_components)
            
    END SUBROUTINE golub_implicit_anypower
   
    SUBROUTINE implicit_trace_estimation ( n, mtz_1, maps, pdb_2, mode,  maxit, trace, bounds, &
                                           hutchinson,  sk, QINV, xp, xpower )
        INTEGER,                                     INTENT(IN)             :: n
        TYPE(mtz),                                   INTENT(INOUT)          :: mtz_1
        TYPE(map),     DIMENSION(:),   ALLOCATABLE,  INTENT(INOUT)          :: maps
        TYPE(pdb),                                   INTENT(IN)             :: pdb_2
        CHARACTER(LEN=*),                            INTENT(IN)             :: mode
        INTEGER,                                     INTENT(IN)             :: maxit
        REAL(KIND=wp),                               INTENT(OUT)            :: trace
        REAL(KIND=wp), DIMENSION(2),                 INTENT(INOUT)          :: bounds
        LOGICAL,                                     INTENT(IN)             :: hutchinson
        REAL(KIND=wp), DIMENSION(:),                 INTENT(IN)             :: sk
        TYPE(csr_matrix),                            INTENT(IN)             :: QINV
        REAL(KIND=wp),                               INTENT(IN),   OPTIONAL :: xp
        INTEGER,                                     INTENT(IN),   OPTIONAL :: xpower
!       Local arrays:
        REAL(KIND=wp), DIMENSION(n)                                         :: z
!       This depends on MAXIT:
        REAL(KIND=wp), DIMENSION(maxit)                                     :: L
        REAL(KIND=wp), DIMENSION(maxit)                                     :: U
        REAL(KIND=wp), DIMENSION(2)                                         :: extreme_lambdas
        REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE                          :: H
!       Local variables:
        INTEGER                                                             :: power
        REAL(KIND=wp)                                                       :: p
        REAL(KIND=wp)                                                       :: Lmin
        REAL(KIND=wp)                                                       :: Umax
        REAL(KIND=wp)                                                       :: Lp 
        REAL(KIND=wp)                                                       :: Up
        REAL(KIND=wp)                                                       :: nu
        REAL(KIND=wp)                                                       :: Ij
        REAL(KIND=wp)                                                       :: dummy 
!       Counters:
        INTEGER                                                             :: i
        INTEGER                                                             :: itn
        INTEGER                                                             :: j
        INTEGER                                                             :: positive
        INTEGER                                                             :: total_sign
        CHARACTER(LEN=32)                                                   :: srname='implicit_trace_estimation'
!       CPU:
        REAL(KIND=sp)                                                       :: time0
        REAL(KIND=sp)                                                       :: time1

        CALL CPU_TIME(time0)
        CALL messag(' ', srname)

!       Set default probability level:
        p = 0.95_wp
        IF ( PRESENT ( xp ) ) THEN
            p = xp
            IF ( p < 0.75_wp ) THEN
                CALL warn('Too low probablity level.', srname)
            ENDIF
        ENDIF
        
        power = -1
        IF ( PRESENT ( xpower ) ) THEN
            power = xpower
        ENDIF

!       Start with smallest lambda:
        CALL implicit_shift_invert(n, 1, 8, 'LM', mtz_1, maps, pdb_2, mode, extreme_lambdas(1), &
                                   dummy, sk = sk)

!       Largest lambda:
        CALL implicit_dsdrv1(n, 1, 10, 'LM', mtz_1, maps, pdb_2, mode, dummy, &
                             extreme_lambdas(2), sk = sk)
         
        WRITE(*,"(' IMPLICIT_TRACE_ESTIMATION> ', 'Extreme lambdas (ARPACK est.)=', 2ES9.2)") &
        extreme_lambdas(1:2)

!       Would be nice to calculate Kantorovich bound here... But it's not obvious how to
!       get correct lambdas for the original unconditioned matrix...

        IF ( Hutchinson ) THEN
            CALL messag('Using Hutchinson method.', srname)
            CALL init_random_number()
        ELSE
            CALL messag('Using Hadamard method.', srname)
            i = 1
            DO j = 1, 30
                i = 2 * i
                IF ( i > n ) EXIT
            ENDDO
            CALL allocate_array ( H, maxit, i)
            CALL hadamard(maxit, i, H)
        ENDIF

        CALL messag(' ', srname)
        CALL messag('                               -Confidence bounds- ', srname)
        CALL messag('      itns      trace(est.)    Lower         Upper       CPU    itns to',  srname)
        CALL messag('                               bound         bound       time   converge', srname)
        CALL messag(' ', srname)

!       Random number quality check:
        total_sign = 0
        positive   = 0

        DO j = 1, maxit

            IF ( HUTCHINSON ) THEN

                CALL RANDOM_NUMBER(z)
                DO i = 1, n
                    IF ( z(i) < 0.5_wp ) THEN
                        z(i) = -1.0_wp 
                    ELSE
                        z(i) = 1.0_wp
                    ENDIF
                ENDDO
            ELSE

!               Use Hadamard row vectors:
                z(1:n) = H(j,1:n)

                IF ( debug > 40 ) THEN
                    IF ( n < 20 ) THEN
                        WRITE(*,"(' z=', 100F4.1)") z(1:n)
                    ENDIF
                ENDIF
            ENDIF

            total_sign = total_sign + n
            positive = positive + COUNT ( z(1:n) > 0 )

!           Estimate z'Tr(A**-1)z:
!            IF ( PRESENT ( sk ) ) THEN
                z = sk * z
                CALL golub_implicit(n, mtz_1, maps, pdb_2, mode, extreme_lambdas, z, bounds, 0, sk,QINV,&
                                    itn=itn)
!            ELSE
!                CALL golub_implicit(n, mtz_1, maps, pdb_2, mode, extreme_lambdas, z, bounds, 0, &
!                                    xpower=power, itn=itn)
!            ENDIF

            L(j) = bounds(1)
            U(j) = bounds(2)

            IF ( L(j) <= 0.0_wp .OR. U(j) <= 0.0_wp ) THEN
                WRITE(*,*) L(j), U(j)
                CALL die('Matrix is not positive definite.', srname)
            ENDIF

            Ij = ( 0.5_wp / j ) * SUM ( L(1:j) + U(1:j) )

            Lmin = MIN ( Lmin, L(j) )
            Umax = MAX ( Umax, U(j) )
                
            nu = SQRT ( -0.5_wp * j * ( Umax - Lmin ) ** 2 * LOG ( ( 1.0_wp - p ) / 2.0_wp ) )

            Lp = ( 1.0_wp / j ) * SUM ( L(1:j) ) - nu / j
            Up = ( 1.0_wp / j ) * SUM ( U(1:j) ) + nu / j

            CALL CPU_TIME(time1)
            WRITE(*,"(' IMPLICIT_TRACE_ESTIMATION> ', 5X, I5, 5X, ES11.4, 2X, 2(2X, ES9.2), F10.1 's', I10)") &
            j, Ij, Lp, Up, time1 - time0, itn
        ENDDO
 
        trace = Ij

        bounds(1) = Lp
        bounds(2) = Up

!       Report quality of random number generation:
        CALL messag(' ',srname)

        WRITE(*,"(' IMPLICIT_TRACE_ESTIMATION> ', 'estimated trace=', ES11.4, /,                                    &
       &          ' IMPLICIT_TRACE_ESTIMATION> ' , I2, '%-probability lower and upper bounds for trace=', 2ES9.2)") &
        trace, NINT ( 100 * p ), bounds
        WRITE(*,"(' IMPLICIT_TRACE_ESTIMATION> ', 'RG quality: +/- sign=', 2I8)") &
        positive, total_sign - positive

        IF ( ALLOCATED ( H ) ) CALL deallocate_array(H)
        
        CALL messag(' ',srname)

    END SUBROUTINE implicit_trace_estimation

    SUBROUTINE general_tridiagonal_solver(n, dl, d, du, b, sigma_or_bound, lambda)
!   
!       Supposed to be used for symmetric but indefinite linear systems.
!

        INTEGER,                     INTENT(IN)            :: n
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: dl
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: d
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: du
        REAL(KIND=wp), DIMENSION(:), INTENT(IN)            :: b
        REAL(KIND=wp),               INTENT(OUT)           :: sigma_or_bound
        REAL(KIND=wp),               INTENT(IN),  OPTIONAL :: lambda
!       Local automatic arrays:
        REAL(KIND=wp), DIMENSION(n-1)                      :: xdl
        REAL(KIND=wp), DIMENSION(n)                        :: xd
        REAL(KIND=wp), DIMENSION(n-1)                      :: xdu
        REAL(KIND=wp), DIMENSION(n-2)                      :: du2
        REAL(KIND=wp), DIMENSION(n)                        :: xb
        INTEGER,       DIMENSION(n)                        :: ipiv
        INTEGER                                            :: info
        CHARACTER(LEN=32),                            SAVE :: srname = 'general_tridiagonal_solver'

        IF ( debug > 1000 ) THEN
            IF ( PRESENT ( lambda ) ) THEN
                WRITE(*,*) ' lambda present'
            ELSE
                WRITE(*,*) ' lambda NOT present'
            ENDIF
        ENDIF
!       Move the data:
        xdl = dl(1:n-1)
        IF ( PRESENT ( lambda ) ) THEN
            xd  = d(1:n) - lambda
        ELSE
            xd = d(1:n)
        ENDIF

        xdu = du(1:n-1)
        xb  = b(1:n)
        CALL gttrf(xdl, xd, xdu, du2, ipiv, info=info)

!       Check output:
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('Problems obtaining LU factorization for tridiagonal symmetric system.', srname)       
        ENDIF

        CALL gttrs(xdl, xd, xdu, du2, xb, ipiv, TRANS='N', info=info)
        IF ( info /= 0 ) THEN
            WRITE(*,*) ' info=', info
            CALL die('Problems solving tridiagonal system based on LU factorization.', srname)
        ENDIF

        IF ( PRESENT ( lambda ) ) THEN
            sigma_or_bound = xb(n)
        ELSE

            IF ( b(1) /= 1.0_wp .OR. ANY ( b(2:n) /= 0.0_wp ) ) THEN
                WRITE(*,*) b(1)
                CALL die('Programming error. Incorrect vector has been supplied.', srname) 
            ENDIF

            sigma_or_bound = xb(1)

        ENDIF

    END SUBROUTINE general_tridiagonal_solver
 
END MODULE estimators
