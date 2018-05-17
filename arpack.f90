MODULE arpack
USE implicit_matrix_inversion
USE second_derivatives, op => Ax
USE sparse_basic
!USE util
IMPLICIT NONE
CONTAINS
    SUBROUTINE dsdrv1 ( n, nev, ncv, which, A, lambda_1, lambda_n, resid_in )
        INTEGER,                               INTENT(IN)           :: n
        INTEGER,                               INTENT(IN)           :: nev
        INTEGER,                               INTENT(IN)           :: ncv
        CHARACTER(LEN=2),                      INTENT(IN)           :: which
        TYPE(csr_matrix),                      INTENT(IN)           :: A
        REAL(KIND=wp),                         INTENT(OUT)          :: lambda_1   ! min lambda
        REAL(KIND=wp),                         INTENT(OUT)          :: lambda_n   ! max lambda
        REAL(KIND=wp),                         INTENT(IN), OPTIONAL :: resid_in
!       Local arrays and variables:
        REAL(KIND=wp), DIMENSION(n,ncv)                             :: v
        REAL(KIND=wp), DIMENSION(ncv*(ncv+8))                       :: workl
        REAL(KIND=wp), DIMENSION(3*n)                               :: workd
        REAL(KIND=wp), DIMENSION(ncv,2)                             :: d
        REAL(KIND=wp), DIMENSION(n)                                 :: resid
        REAL(KIND=wp), DIMENSION(n)                                 :: ax
        LOGICAL,       DIMENSION(ncv)                               :: select
        INTEGER,       DIMENSION(11)                                :: iparam
        INTEGER,       DIMENSION(11)                                :: ipntr
!       Local variables:
        CHARACTER(LEN=1)                                            :: bmat
        INTEGER                                                     :: ido
        INTEGER                                                     :: ldv
        INTEGER                                                     :: lworkl
        INTEGER                                                     :: info
        INTEGER                                                     :: ierr
        INTEGER                                                     :: j
        INTEGER                                                     :: ishfts
        INTEGER                                                     :: maxitr
        INTEGER                                                     :: mode
        INTEGER                                                     :: nconv
        LOGICAL                                                     :: rvec
        REAL(KIND=wp)                                               :: tol
        REAL(KIND=wp)                                               :: sigma
        REAL(KIND=wp), PARAMETER                                    :: zero = 0.0_wp

!       BLAS LAPACK routines
        REAL(KIND=wp)                                               :: dnrm2
        EXTERNAL                                                    :: dnrm2
        EXTERNAL                                                    :: daxpy

!       Intrinsic
        INTRINSIC                                                   :: ABS
        INCLUDE 'arpack_debug.h'

!       arpack_debug.h
        ndigit = -3
        logfil = 6
        msgets = 0
        msaitr = 1   ! default = 0 
        msapps = 0
        msaupd = 1   ! default = 0
        msaup2 = 0
        mseigt = 0
        mseupd = 0   ! default = 0 (was 2)

!       Conventional (simple) eigenvalue problem:
        bmat = 'I'
!       which = 'BE'
        ldv = n 

!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero
      info = 0
      ido = 0
 
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 3000
      mode   = 1
!      
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid,     &
                       ncv, v, ldv, iparam, ipntr, workd, workl, &
                       lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!

!            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
            call matvec (n, workd(ipntr(2):ipntr(2)+n-1), workd(ipntr(1):ipntr(1)+n-1), A)
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
                                                                  
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,        &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!

!                call av(nx, v(1,j), ax)
                call matvec ( n, ax, v(1:n,j), A)

                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!

!            MAXNCV replaced with ncv:
             call dmout(6, nconv, 2, d, ncv, -6,      &
                 'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue

!     Report extreme eigenvalues: 
      lambda_1 = MINVAL ( d(1:nconv,1) )
      lambda_n = MAXVAL ( d(1:nconv,1) )
      WRITE(*,"(' DSDRV1> ', 'lambda min=', ES9.2, ' lambda max=', ES9.2)") &
      lambda_1, lambda_n
    END SUBROUTINE dsdrv1
    
    SUBROUTINE implicit_dsdrv1 ( n, nev, ncv, which, mtz_1, maps, pdb_2, refmode, &
                                 lambda_1, lambda_n, sk, Q, resid_in )
        INTEGER,                                  INTENT(IN)           :: n
        INTEGER,                                  INTENT(IN)           :: nev
        INTEGER,                                  INTENT(IN)           :: ncv
        CHARACTER(LEN=2),                         INTENT(IN)           :: which
        TYPE(mtz),                                INTENT(INOUT)        :: mtz_1
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                         INTENT(IN)           :: refmode
        REAL(KIND=wp),                            INTENT(OUT)          :: lambda_1   ! min lambda
        REAL(KIND=wp),                            INTENT(OUT)          :: lambda_n   ! max lambda
        REAL(KIND=wp), DIMENSION(:),              INTENT(IN), OPTIONAL :: sk
        TYPE(csr_matrix),                         INTENT(IN), OPTIONAL :: Q
        REAL(KIND=wp),                            INTENT(IN), OPTIONAL :: resid_in
!       Local arrays and variables:
        REAL(KIND=wp), DIMENSION(n,ncv)                                :: v
        REAL(KIND=wp), DIMENSION(ncv*(ncv+8))                          :: workl
        REAL(KIND=wp), DIMENSION(3*n)                                  :: workd
        REAL(KIND=wp), DIMENSION(ncv,2)                                :: d
        REAL(KIND=wp), DIMENSION(n)                                    :: resid
        REAL(KIND=wp), DIMENSION(n)                                    :: ax
        LOGICAL,       DIMENSION(ncv)                                  :: select
        INTEGER,       DIMENSION(11)                                   :: iparam
        INTEGER,       DIMENSION(11)                                   :: ipntr
!       Local variables:
        CHARACTER(LEN=1)                                               :: bmat
        INTEGER                                                        :: ido
        INTEGER                                                        :: ldv
        INTEGER                                                        :: lworkl
        INTEGER                                                        :: info
        INTEGER                                                        :: ierr
        INTEGER                                                        :: j
        INTEGER                                                        :: ishfts
        INTEGER                                                        :: maxitr
        INTEGER                                                        :: mode
        INTEGER                                                        :: nconv
        LOGICAL                                                        :: rvec
        REAL(KIND=wp)                                                  :: tol
        REAL(KIND=wp)                                                  :: sigma
        REAL(KIND=wp), PARAMETER                                       :: zero = 0.0_wp

!       BLAS LAPACK routines
        REAL(KIND=wp)                                                  :: dnrm2
        EXTERNAL                                                       :: dnrm2
        EXTERNAL                                                       :: daxpy

!       Intrinsic
        INTRINSIC                                                      :: ABS
!       CPU:
        REAL(KIND=sp)                                                  :: time0
        REAL(KIND=sp)                                                  :: time1
!       Counter:
        INTEGER                                                        :: mult
!       SUB:
        CHARACTER(LEN=32)                                              :: srname='implicit_dsdrv1'
        INCLUDE 'arpack_debug.h'

!       arpack_debug.h
        ndigit = -3
        logfil = 6
        msgets = 0
        msaitr = 1   ! default = 0 
        msapps = 0
        msaupd = 1   ! default = 0
        msaup2 = 0
        mseigt = 0
        mseupd = 0   ! default = 0 (was 2)

!       Conventional (simple) eigenvalue problem:
        bmat = 'I'
!       which = 'BE'
        ldv = n 

!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
!      tol = zero
      tol = 1.0E-5
      info = 0
      ido = 0
 
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 30000
      mode   = 1
!      
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      CALL CPU_TIME(time0)
      mult = 0
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid,     &
                       ncv, v, ldv, iparam, ipntr, workd, workl, &
                       lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!

!            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
!            call matvec (n, workd(ipntr(2):ipntr(2)+n-1), workd(ipntr(1):ipntr(1)+n-1), A)
            IF ( PRESENT ( Q ) ) THEN
                CALL cond_Ax(n, workd(ipntr(2):ipntr(2)+n-1), workd(ipntr(1):ipntr(1)+n-1), &
                             mtz_1, maps, pdb_2, refmode, sk, Q)
            ELSE
                CALL die('Program not ready for that option...', srname)
!                CALL op(n,  workd(ipntr(2):ipntr(2)+n-1), workd(ipntr(1):ipntr(1)+n-1), &
!                        mtz_1, maps, pdb_2, refmode)
            ENDIF

            mult = mult + 1

            CALL CPU_TIME(time1)
            IF ( mult < 10 .OR. MOD ( mult, 10 ) == 0 ) THEN
                IF ( mult == 1 ) THEN
                    WRITE(*,"(' IMPLICIT_DSDRV1> ', I8, ' matrix-vector product has been calculated.', &
           &                  '   CPU time=', F10.1, ' s')") mult, time1 - time0
                ELSE
                    WRITE(*,"(' IMPLICIT_DSDRV1> ', I8, ' matrix-vector products have been calculated.', &
           &                  ' CPU time=', F10.1, ' s')") mult, time1 - time0
                ENDIF
                IF ( MOD ( mult, 10 ) == 0 ) CALL messag(' ', srname)
            ENDIF
!         
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
                                                                  
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,        &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!

!               call av(nx, v(1,j), ax)
!               call matvec ( n, ax, v(1:n,j), A)

!                IF ( PRESENT ( sk ) ) THEN
                    CALL cond_Ax(n,  ax, v(1:n,j), mtz_1, maps, pdb_2, refmode, sk, Q)
!                ELSE
!                    CALL op(n,  ax, v(1:n,j), mtz_1, maps, pdb_2, refmode)
!                ENDIF

                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!

!            MAXNCV replaced with ncv:
             call dmout(6, nconv, 2, d, ncv, -6,      &
                 'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue

!     Report extreme eigenvalues: 
      lambda_1 = MINVAL ( d(1:nconv,1) )
      lambda_n = MAXVAL ( d(1:nconv,1) )
      WRITE(*,"(' IMPLICIT_DSDRV1> ', 'lambda min=', ES9.2, ' lambda max=', ES9.2)") &
      lambda_1, lambda_n
      CALL messag(' ', srname)
    END SUBROUTINE implicit_dsdrv1

    SUBROUTINE implicit_shift_invert( n, nev, ncv, which, mtz_1, maps, pdb_2, refmode, &
                                      lambda_1, lambda_n, sk, ACSR, resid_in )
        INTEGER,                                  INTENT(IN)           :: n
        INTEGER,                                  INTENT(IN)           :: nev
        INTEGER,                                  INTENT(IN)           :: ncv
        CHARACTER(LEN=2),                         INTENT(IN)           :: which
        TYPE(mtz),                                INTENT(INOUT)        :: mtz_1
        TYPE(map),     DIMENSION(:), ALLOCATABLE, INTENT(INOUT)        :: maps
        TYPE(pdb),                                INTENT(IN)           :: pdb_2
        CHARACTER(LEN=*),                         INTENT(IN)           :: refmode
        REAL(KIND=wp),                            INTENT(OUT)          :: lambda_1   ! min lambda
        REAL(KIND=wp),                            INTENT(OUT)          :: lambda_n   ! max lambda
        REAL(KIND=wp), DIMENSION(:),              INTENT(IN)           :: sk
        TYPE(csr_matrix),                                     OPTIONAL :: ACSR
        REAL(KIND=wp),                            INTENT(IN), OPTIONAL :: resid_in
!       Local arrays and variables:
        REAL(KIND=wp), DIMENSION(n,ncv)                                :: v
        REAL(KIND=wp), DIMENSION(ncv*(ncv+8))                          :: workl
        REAL(KIND=wp), DIMENSION(3*n)                                  :: workd
        REAL(KIND=wp), DIMENSION(ncv,2)                                :: d
        REAL(KIND=wp), DIMENSION(n)                                    :: resid
        REAL(KIND=wp), DIMENSION(n)                                    :: ax
        LOGICAL,       DIMENSION(ncv)                                  :: select
        INTEGER,       DIMENSION(11)                                   :: iparam
        INTEGER,       DIMENSION(11)                                   :: ipntr
!       Local variables:
        CHARACTER(LEN=1)                                               :: bmat
        INTEGER                                                        :: ido
        INTEGER                                                        :: ldv
        INTEGER                                                        :: lworkl
        INTEGER                                                        :: info
        INTEGER                                                        :: ierr
        INTEGER                                                        :: j
        INTEGER                                                        :: ishfts
        INTEGER                                                        :: maxitr
        INTEGER                                                        :: mode
        INTEGER                                                        :: nconv
        LOGICAL                                                        :: rvec
        REAL(KIND=wp)                                                  :: tol
        REAL(KIND=wp)                                                  :: sigma
        REAL(KIND=wp), PARAMETER                                       :: zero = 0.0_wp

!       BLAS LAPACK routines
        REAL(KIND=wp)                                                  :: dnrm2
        EXTERNAL                                                       :: dnrm2
        EXTERNAL                                                       :: daxpy
!       SYMMLQ:
        LOGICAL                                                        :: checkA
        LOGICAL                                                        :: precon
        INTEGER                                                        :: itnlim
        INTEGER                                                        :: nout
        REAL(KIND=wp)                                                  :: rtol
        INTEGER                                                        :: istop
        INTEGER                                                        :: itn
        REAL(KIND=wp)                                                  :: Anorm
        REAL(KIND=wp)                                                  :: Acond
        REAL(KIND=wp)                                                  :: rnorm
        REAL(KIND=wp)                                                  :: ynorm
!       Intrinsic
        INTRINSIC                                                      :: ABS
!       CPU:
        REAL(KIND=sp)                                                  :: time0
        REAL(KIND=sp)                                                  :: time1
!       Counter:
        INTEGER                                                        :: mult
!       SUB:
        CHARACTER(LEN=32)                                              :: srname='implicit_shift_invert'
        INCLUDE 'arpack_debug.h'

!       arpack_debug.h
        ndigit = -3
        logfil = 6
        msgets = 0
        msaitr = 1   ! default = 0 
        msapps = 0
        msaupd = 1   ! default = 0
        msaup2 = 0
        mseigt = 0
        mseupd = 0   ! default = 0 (was 2)

!       Conventional (simple) eigenvalue problem:
        bmat = 'I'
!       which = 'BE'
        sigma = zero
        ldv = n 

!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
!      tol = zero
      tol = 1.0E-5
      info = 0
      ido = 0
 
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 3 
!      
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      CALL CPU_TIME(time0)
      mult = 0

!     Initialise SYMMLQ parameters:
      checkA = .FALSE.
      precon = .FALSE.
      itnlim =  n / 10
      nout = 6 
      rtol = 0.1_wp ** 5
      
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid,     &
                       ncv, v, ldv, iparam, ipntr, workd, workl, &
                       lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!

!            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
!            call matvec (n, workd(ipntr(2):ipntr(2)+n-1), workd(ipntr(1):ipntr(1)+n-1), A)
            IF ( .NOT. PRESENT ( ACSR )  ) THEN

                CALL symmlq(n, OP, Diagonal_Msolve, workd(ipntr(1):ipntr(1)+n-1),    &
                            sigma, checkA, precon, workd(ipntr(2):ipntr(2)+n-1), &
                            itnlim, nout, rtol, istop, itn,                      &
                            Anorm, Acond, rnorm, ynorm,                          &
                            mtz_1, maps, pdb_2, refmode, sk)

            ELSE

!               ******
!               We need to modify symmlq and Ax=OP in case precon=.FALSE. and PRESENT ( ACSR )
!               The following product has to be computed in Aprod=OP=AX:
!               ACSR^-1/2 * A * (ACSR^-1/2 * x):
!
!               a) in call to Aprod when precon=.FALSE., Acsr must participate.
!               b) Ax has to be modified. Add Acsr
!               c) Need to prepare block diagonal matirx ACSR^-1/2 (this will keep matrix symmetry)
!               *******
!               Overall we need VERY GOOD preconditioner if want to estimate trace and calculate
!               standard deviations quickly.
!
!               True result can be recovered then by multiplying ACSR^1/2 * x before calling
!               corresponding routines, i.e. u'f(A)u = u' ACSR^1/2 * [ACSR^-1/2 * A * ACSR^-1/2] * ACSR^1/2 u
!
!               Matrix in squared brackets should be as close to unit matrix as possible then everything
!               will be easy. 
!
!               resolution:               block:
!               1.2 A                     atom (5x5)
!               1.5 A                     residue
!               2.0 A                     2 residues
!               2.2 -2.5 A                4 residues blocks this should amount to 40-50 * millions elements ~ 1 - GB
!                                         for two matrices. So machines should have 4-8 GB memory at least.
!
!               lower                     better not to think about it... but restraints should help enormously 
!                                         for all resolution ranges
!
!               This clearly shows we need to keep TWO sparse matrices in memory.
!               For larger proteins we need at least 8 - 16 GB of memory and 10 - 100 processors
! FIXME
                CALL messag('Calculating smallest eigenvalue of Q**(-1/2)*H*Q**(-1/2)...', srname)


!               diagonal_Msolve will not be used since PRECON=.FALSE.:                
                CALL symmlq(n, cond_Ax, diagonal_Msolve, workd(ipntr(1):ipntr(1)+n-1),    &
                            sigma, checkA, precon, workd(ipntr(2):ipntr(2)+n-1), &
                            itnlim, nout, rtol, istop, itn,                      &
                            Anorm, Acond, rnorm, ynorm,                          &
                            mtz_1, maps, pdb_2, refmode, sk, ACSR)
                
            ENDIF

            mult = mult + 1

            CALL CPU_TIME(time1)
            IF ( mult == 1 ) THEN
                WRITE(*,"(' IMPLICIT_SHIFT_INVERT> ', I8, ' linear system has been solved. (', A,' itns)', &
           &              '   CPU time=', F10.1, ' s')") mult, TRIM ( int_to_c ( itn ) ), time1 - time0
             ELSE
                    WRITE(*,"(' IMPLICIT_SHIFT_INVERT> ', I8, ' linear systems have been solved. (' A,' itns)', &
           &                  ' CPU time=', F10.1, ' s')") mult, TRIM ( int_to_c ( itn ) ), time1 - time0
             ENDIF
             IF ( MOD ( mult, 10 ) == 0 ) CALL messag(' ', srname)
!         
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
                                                                  
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,        &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!

!               call av(nx, v(1,j), ax)
!               call matvec ( n, ax, v(1:n,j), A)

                IF ( .NOT. PRESENT ( ACSR ) ) THEN
                    CALL op(n,  ax, v(1:n,j), mtz_1, maps, pdb_2, refmode, sk)
                    
                ELSE
!                   CALL cond_Ax(n,  ax, v(1:n,j), mtz_1, maps, pdb_2, refmode, sk, Q)
                    CALL cond_Ax(n,  ax, v(1:n,j), mtz_1, maps, pdb_2, refmode, sk, ACSR)
                ENDIF

                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!

!            MAXNCV replaced with ncv:
             call dmout(6, nconv, 2, d, ncv, -6,      &
                 'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue

!     Report extreme eigenvalues: 
      lambda_1 = MINVAL ( d(1:nconv,1) )
      lambda_n = MAXVAL ( d(1:nconv,1) )
      WRITE(*,"(' IMPLICIT_SHIFT_INVERT> ', 'lambda min=', ES9.2, ' lambda max=', ES9.2)") &
      lambda_1, lambda_n
      CALL messag(' ', srname)
    END SUBROUTINE implicit_shift_invert


END MODULE arpack

