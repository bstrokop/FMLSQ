MODULE implicit_matrix_inversion
USE second_derivatives
USE sparse_basic
IMPLICIT NONE
PUBLIC   :: MINRES
PUBLIC   :: SYMMLQ
CONTAINS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minresModule.f90
!
! MINRES solves symmetric systems Ax = b or min ||Ax - b||_2,
! where the matrix A may be indefinite and/or singular.
!
! Contributors:
!     Chris Paige <chris@cs.mcgill.ca>
!     Sou-Cheng Choi <scchoi@stanford.edu>
!
!     Michael Saunders <saunders@stanford.edu>
!     Systems Optimization Laboratory (SOL)
!     Stanford University
!     Stanford, CA 94305-4026, USA
!     (650)723-1875
!
! 09 Oct 2007: F90 version constructed from the F77 version.
!              Initially used compiler option -r8, but this is nonstandard.
! 15 Oct 2007: Test on Arnorm = ||Ar|| added to recognize singular systems.
! 15 Oct 2007: Temporarily used real(8) everywhere.
! 16 Oct 2007: Use minresDataModule to define dp = selected_real_kind(15).
!              We need "use minresDataModule"
!              at the beginning of modules AND inside interfaces.
!
!              g95 compiles successfully with the following options:
!   g95 -c -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full minresModule.f90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE minres( n, Aprod, Msolve, b, shift, checkA, precon,     &
                     x, itnlim, nout, rtol,                          &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm, &
                     sk, ACSR )

    INTEGER,                     INTENT(IN)  :: n
    INTEGER,                     INTENT(IN)  :: itnlim
    INTEGER,                     INTENT(IN)  :: nout
    LOGICAL,                     INTENT(IN)  :: checkA
    LOGICAL,                     INTENT(IN)  :: precon
    REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: b
    REAL(KIND=wp),               INTENT(IN)  :: shift
    REAL(KIND=wp),               INTENT(IN)  :: rtol
    REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: x
    INTEGER,                     INTENT(OUT) :: istop
    INTEGER,                     INTENT(OUT) :: itn
    REAL(KIND=wp),               INTENT(OUT) :: Anorm
    REAL(KIND=wp),               INTENT(OUT) :: Acond
    REAL(KIND=wp),               INTENT(OUT) :: rnorm
    REAL(KIND=wp),               INTENT(OUT) :: Arnorm
    REAL(KIND=wp),               INTENT(OUT) :: ynorm
    REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: sk
    TYPE(csr_matrix),            INTENT(IN)  :: ACSR

    INTERFACE

       SUBROUTINE Aprod (n, y, x, ACSR)                       ! y := A*x
         USE select_kinds
         USE sparse_basic
         INTEGER,  INTENT(IN)                                 :: n
         REAL(KIND=wp), DIMENSION(:), INTENT(IN)              :: x
         REAL(KIND=wp), DIMENSION(n), INTENT(OUT)             :: y
         TYPE(csr_matrix),            INTENT(IN)              :: ACSR
       END SUBROUTINE Aprod

       SUBROUTINE Msolve(n,y,x,diag,A)                     ! Solve M*y = x
         USE select_kinds 
         USE sparse_basic
         INTEGER,         INTENT(IN)                         :: n
         REAL(KIND=wp),   DIMENSION(n), INTENT(OUT)          :: y
         REAL(KIND=wp),   DIMENSION(n), INTENT(IN)           :: x
         REAL(KIND=wp),   DIMENSION(n), INTENT(IN)           :: diag
         TYPE(csr_matrix),              INTENT(IN), OPTIONAL :: A
       END SUBROUTINE Msolve

    END INTERFACE

    !-------------------------------------------------------------------
    !
    ! MINRES  is designed to solve the system of linear equations
    !
    !    Ax = b
    !
    ! or the least-squares problem
    !
    !    min ||Ax - b||_2,
    !
    ! where A is an n by n symmetric matrix and b is a given vector.
    ! The matrix A may be indefinite and/or singular.
    !
    ! 1. If A is known to be positive definite, the Conjugate Gradient
    ! Method might be preferred, since it requires the same number
    ! of iterations as MINRES but less work per iteration.
    !
    ! 2. If A is indefinite but Ax = b is known to have a solution
    ! (e.g. if A is nonsingular), SYMMLQ might be preferred,
    ! since it requires the same number of iterations as MINRES
    ! but slightly less work per iteration.
    !
    ! The matrix A is intended to be large and sparse.  It is accessed
    ! by means of a subroutine call of the form
    ! SYMMLQ development:
    !
    !    call Aprod ( n, x, y )
    !
    ! which must return the product y = Ax for any given vector x.
    !
    !
    ! More generally, MINRES is designed to solve the system
    !
    !    (A - shift*I) x = b
    ! or
    !    min ||(A - shift*I) x - b||_2,
    !
    ! where  shift  is a specified scalar value.  Again, the matrix
    ! (A - shift*I) may be indefinite and/or singular.
    ! The work per iteration is very slightly less if  shift = 0.
    !
    ! Note: If  shift  is an approximate eigenvalue of  A
    ! and  b  is an approximate eigenvector,  x  might prove to be
    ! a better approximate eigenvector, as in the methods of
    ! inverse iteration and/or Rayleigh-quotient iteration.
    ! However, we're not yet sure on that -- it may be better to use SYMMLQ.
    !
    ! A further option is that of preconditioning, which may reduce
    ! the number of iterations required.  If M = C C' is a positive
    ! definite matrix that is known to approximate  (A - shift*I)
    ! in some sense, and if systems of the form  My = x  can be
    ! solved efficiently, the parameters precon and Msolve may be
    ! used (see below).  When  precon = .true., MINRES will
    ! implicitly solve the system of equations
    !
    !    P (A - shift*I) P' xbar  =  P b,
    !
    ! i.e.             Abar xbar  =  bbar
    ! where                    P  =  C**(-1),
    !                       Abar  =  P (A - shift*I) P',
    !                       bbar  =  P b,
    !
    ! and return the solution       x  =  P' xbar.
    ! The associated residual is rbar  =  bbar - Abar xbar
    !                                  =  P (b - (A - shift*I)x)
    !                                  =  P r.
    !
    ! In the discussion below, eps refers to the machine precision.
    !
    ! Parameters
    ! ----------
    !
    ! n       input      The dimension of the matrix A.
    ! b(n)    input      The rhs vector b.
    ! x(n)    output     Returns the computed solution x.
    !
    ! Aprod   external   A subroutine defining the matrix A.
    !                       call Aprod ( n, x, y )
    !                    must return the product y = Ax
    !                    without altering the vector x.
    !
    ! Msolve  external   An optional subroutine defining a
    !                    preconditioning matrix M, which should
    !                    approximate (A - shift*I) in some sense.
    !                    M must be positive definite.
    !
    !                       call Msolve( n, x, y )
    !
    !                    must solve the linear system My = x
    !                    without altering the vector x.
    !
    !                    In general, M should be chosen so that Abar has
    !                    clustered eigenvalues.  For example,
    !                    if A is positive definite, Abar would ideally
    !                    be close to a multiple of I.
    !                    If A or A - shift*I is indefinite, Abar might
    !                    be close to a multiple of diag( I  -I ).
    !
    ! checkA  input      If checkA = .true., an extra call of Aprod will
    !                    be used to check if A is symmetric.  Also,
    !                    if precon = .true., an extra call of Msolve
    !                    will be used to check if M is symmetric.
    !
    ! precon  input      If precon = .true., preconditioning will
    !                    be invoked.  Otherwise, subroutine Msolve
    !                    will not be referenced; in this case the
    !                    actual parameter corresponding to Msolve may
    !                    be the same as that corresponding to Aprod.
    !
    ! shift   input      Should be zero if the system Ax = b is to be
    !                    solved.  Otherwise, it could be an
    !                    approximation to an eigenvalue of A, such as
    !                    the Rayleigh quotient b'Ab / (b'b)
    !                    corresponding to the vector b.
    !                    If b is sufficiently like an eigenvector
    !                    corresponding to an eigenvalue near shift,
    !                    then the computed x may have very large
    !                    components.  When normalized, x may be
    !                    closer to an eigenvector than b.
    !
    ! nout    input      A file number.
    !                    If nout > 0, a summary of the iterations
    !                    will be printed on unit nout.
    !
    ! itnlim  input      An upper limit on the number of iterations.
    !
    ! rtol    input      A user-specified tolerance.  MINRES terminates
    !                    if it appears that norm(rbar) is smaller than
    !                       rtol * norm(Abar) * norm(xbar),
    !                    where rbar is the transformed residual vector,
    !                       rbar = bbar - Abar xbar.
    !
    !                    If shift = 0 and precon = .false., MINRES
    !                    terminates if norm(b - A*x) is smaller than
    !                       rtol * norm(A) * norm(x).
    !
    ! istop   output     An integer giving the reason for termination...
    !
    !          -1        beta2 = 0 in the Lanczos iteration; i.e. the
    !                    second Lanczos vector is zero.  This means the
    !                    rhs is very special.
    !                    If there is no preconditioner, b is an
    !                    eigenvector of A.
    !                    Otherwise (if precon is true), let My = b.
    !                    If shift is zero, y is a solution of the
    !                    generalized eigenvalue problem Ay = lambda My,
    !                    with lambda = alpha1 from the Lanczos vectors.
    !
    !                    In general, (A - shift*I)x = b
    !                    has the solution         x = (1/alpha1) y
    !                    where My = b.
    !
    !           0        b = 0, so the exact solution is x = 0.
    !                    No iterations were performed.
    !
    !           1        Norm(rbar) appears to be less than
    !                    the value  rtol * norm(Abar) * norm(xbar).
    !                    The solution in  x  should be acceptable.
    !
    !           2        Norm(rbar) appears to be less than
    !                    the value  eps * norm(Abar) * norm(xbar).
    !                    This means that the residual is as small as
    !                    seems reasonable on this machine.
    !
    !           3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
    !                    which should indicate that x has essentially
    !                    converged to an eigenvector of A
    !                    corresponding to the eigenvalue shift.
    !
    !           4        Acond (see below) has exceeded 0.1/eps, so
    !                    the matrix Abar must be very ill-conditioned.
    !                    x may not contain an acceptable solution.
    !
    !           5        The iteration limit was reached before any of
    !                    the previous criteria were satisfied.
    !
    !           6        The matrix defined by Aprod does not appear
    !                    to be symmetric.
    !                    For certain vectors y = Av and r = Ay, the
    !                    products y'y and r'v differ significantly.
    !
    !           7        The matrix defined by Msolve does not appear
    !                    to be symmetric.
    !                    For vectors satisfying My = v and Mr = y, the
    !                    products y'y and r'v differ significantly.
    !
    !           8        An inner product of the form  x' M**(-1) x
    !                    was not positive, so the preconditioning matrix
    !                    M does not appear to be positive definite.
    !
    !                    If istop >= 5, the final x may not be an
    !                    acceptable solution.
    !
    ! itn     output     The number of iterations performed.
    !
    ! Anorm   output     An estimate of the norm of the matrix operator
    !                    Abar = P (A - shift*I) P',   where P = C**(-1).
    !
    ! Acond   output     An estimate of the condition of Abar above.
    !                    This will usually be a substantial
    !                    under-estimate of the true condition.
    !
    ! rnorm   output     An estimate of the norm of the final
    !                    transformed residual vector,
    !                       P (b  -  (A - shift*I) x).
    !
    ! ynorm   output     An estimate of the norm of xbar.
    !                    This is sqrt( x'Mx ).  If precon is false,
    !                    ynorm is an estimate of norm(x).
    !-------------------------------------------------------------------
    ! MINRES is an implementation of the algorithm described in
    ! the following reference:
    !
    ! C. C. Paige and M. A. Saunders (1975),
    ! Solution of sparse indefinite systems of linear equations,
    ! SIAM J. Numer. Anal. 12(4), pp. 617-629.
    !-------------------------------------------------------------------
    !
    !
    ! MINRES development:
    !    1972: First version, similar to original SYMMLQ.
    !          Later lost @#%*!
    !    Oct 1995: Tried to reconstruct MINRES from
    !              1995 version of SYMMLQ.
    ! 30 May 1999: Need to make it more like LSQR.
    !              In middle of major overhaul.
    ! 19 Jul 2003: Next attempt to reconstruct MINRES.
    !              Seems to need two vectors more than SYMMLQ.  (w1, w2)
    !              Lanczos is now at the top of the loop,
    !              so the operator Aprod is called in just one place
    !              (not counting the initial check for symmetry).
    ! 22 Jul 2003: Success at last.  Preconditioning also works.
    !              minres.f added to http://www.stanford.edu/group/SOL/.
    !
    ! 16 Oct 2007: Added a stopping rule for singular systems,
    !              as derived in Sou-Cheng Choi's PhD thesis.
    !              Note that ||Ar|| small => r is a null vector for A.
    !              Subroutine minrestest2 in minresTestModule.f90
    !              tests this option.  (NB: Not yet working.)
    !-------------------------------------------------------------------

!     Local arrays and variables
      REAL(KIND=wp)  :: r1(n), r2(n), v(n), w(n), w1(n), w2(n), y(n)
      real(KIND=wp)  :: alfa  , beta  , beta1 , cs    ,          &
                        dbar  , delta , denom , diag  ,          &
                        eps   , epsa  , epsln , epsr  , epsx  ,  &
                        gamma , gbar  , gmax  , gmin  ,          &
                        oldb  , oldeps, qrnorm, phi   , phibar,  &
                        rhs1  , rhs2  , rnorml, rootl ,          &
                        Arnorml,        relArnorml,              &
                        s     , sn    , t     , tnorm2, ynorm2, z
      INTEGER        :: i
      LOGICAL        :: debug, prnt

    ! Local constants
      real(KIND=wp),         PARAMETER :: zero =  0.0_wp,  one = 1.0_wp
      real(KIND=wp),         PARAMETER :: ten  = 10.0_wp
      CHARACTER(LEN=*), PARAMETER :: enter = ' Enter MINRES.  '
      CHARACTER(LEN=*), PARAMETER :: exitt = ' Exit  MINRES.  '
      CHARACTER(LEN=*), PARAMETER :: msg(-1:8) =                  &
        (/ 'beta2 = 0.  If M = I, b and x are eigenvectors of A', & ! -1
           'beta1 = 0.  The exact solution is  x = 0           ', & !  0
           'Requested accuracy achieved, as determined by rtol ', & !  1
           'Reasonable accuracy achieved, given eps            ', & !  2
           'x has converged to an eigenvector                  ', & !  3
           'Acond has exceeded 0.1/eps                         ', & !  4
           'The iteration limit was reached                    ', & !  5
           'Aprod  does not define a symmetric matrix          ', & !  6
           'Msolve does not define a symmetric matrix          ', & !  7
           'Msolve does not define a pos-def preconditioner    ' /) !  8
    !-------------------------------------------------------------------

    INTRINSIC       :: ABS, DOT_PRODUCT, EPSILON, MIN, MAX, SQRT

    ! Print heading and initialize.

    debug = .FALSE.
    eps   = EPSILON ( 1.0_wp )

    IF (nout > 0) THEN
        WRITE(nout, 1000) enter, n, checkA, precon, itnlim, rtol, shift
    ENDIF

    istop  = 0
    itn    = 0
    Anorm  = zero
    Acond  = zero
    rnorm  = zero
    ynorm  = zero
    x(1:n) = zero
    !-------------------------------------------------------------------
    ! Set up y and v for the first Lanczos vector v1.
    ! y = beta1 P' v1, where P = C**(-1).
    ! v is really P' v1.
    !-------------------------------------------------------------------
    y      = b
    r1     = b

!   if ( precon ) call Msolve( n, b, y )
    IF ( precon ) CALL Msolve( n, y, b, sk )

    beta1  = DOT_PRODUCT ( b, y )

    if (beta1 < zero) then     ! M must be indefinite.
       istop = 8
       go to 900
    end if

    if (beta1 == zero) then    ! b = 0 exactly.  Stop with x = 0.
       istop = 0
       go to 900
    end if

    beta1  = SQRT ( beta1 )     ! Normalize y to get v1 later.

    !-------------------------------------------------------------------
    ! See if Msolve is symmetric.
    !-------------------------------------------------------------------
    IF ( checkA  .AND.  precon ) THEN
!       call Msolve( n, y, r2 )
       CALL Msolve( n, r2, y, sk )
       s      = DOT_PRODUCT ( y ,y )
       t      = DOT_PRODUCT ( r1, r2 )
       z      = ABS ( s - t )
       epsa   = (s + eps) * eps**0.33333
       if (z > epsa) then
          istop = 7
          GOTO 900
        ENDIF
    ENDIF

    !-------------------------------------------------------------------
    ! See if Aprod  is symmetric.  Initialize Arnorm.
    !-------------------------------------------------------------------
    IF ( checkA ) THEN
!       call Aprod ( n, y, w )
       CALL Aprod ( n, w, y, ACSR )
!      call Aprod ( n, w, r2 )
       CALL Aprod ( n, r2, w, ACSR )
       s      = DOT_PRODUCT ( w, w )
       t      = DOT_PRODUCT ( y, r2 )
       z      = ABS ( s - t )
       epsa   = (s + eps) * eps**0.33333
       IF ( nout > 0 ) THEN
           WRITE(nout, "(21X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa
       ENDIF

       IF ( z > epsa ) THEN
           istop = 6
           GOTO 900
       ENDIF
       Arnorml = SQRT ( s ) ;

    ELSE

!      call Aprod ( n, y, w )
       CALL Aprod ( n, w, y, ACSR  )
       Arnorml = SQRT ( DOT_PRODUCT ( w, w ) )

    ENDIF

    !-------------------------------------------------------------------
    ! Initialize other quantities.
    !-------------------------------------------------------------------
    oldb   = zero
    beta   = beta1
    dbar   = zero
    epsln  = zero
    qrnorm = beta1
    phibar = beta1
    rhs1   = beta1
    rhs2   = zero
    tnorm2 = zero
    ynorm2 = zero
    cs     = - one
    sn     = zero
    w(1:n) = zero
    w2(1:n)= zero
    r2(1:n)= r1

    if (debug) then
       write(*,*) ' '
       write(*,*) 'b    ', b
       write(*,*) 'beta ', beta
       write(*,*) ' '
    end if

    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
       itn = itn + 1               ! k = itn = 1 first time through

       !----------------------------------------------------------------
       ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
       ! The general iteration is similar to the case k = 1 with v0 = 0:
       !
       !   p1      = Operator * v1  -  beta1 * v0,
       !   alpha1  = v1'p1,
       !   q2      = p2  -  alpha1 * v1,
       !   beta2^2 = q2'q2,
       !   v2      = (1/beta2) q2.
       !
       ! Again, y = betak P vk,  where  P = C**(-1).
       ! .... more description needed.
       !----------------------------------------------------------------
       s      = one / beta            ! Normalize previous vector (in y).
       v      = s*y(1:n)              ! v = vk if P = I

!      call Aprod ( n, v, y )
       CALL Aprod ( n, y, v, ACSR )
       y      = y - shift*v           ! call daxpy ( n, (- shift), v, 1, y, 1 )
       if (itn >= 2) then
          y   = y - (beta/oldb)*r1    ! call daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
       end if

       alfa   = DOT_PRODUCT ( v,y )      ! alphak
       y      = y - (alfa/beta)*r2       ! call daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
       r1     = r2
       r2     = y
!       if ( precon ) call Msolve( n, r2, y )
       IF ( precon ) CALL Msolve( n, y, r2, sk )

       oldb   = beta                      ! oldb = betak
       beta   = DOT_PRODUCT ( r2, y )     ! beta = betak+1^2
       if (beta < zero) then
          istop = 6
          go to 900
       end if

       beta   = SQRT ( beta )          ! beta = betak+1
       tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

       if (itn == 1) then                   ! Initialize a few things.
          if (beta/beta1 <= ten*eps) then   ! beta2 = 0 or ~ 0.
             istop = -1                     ! Terminate later.
          end if
         !tnorm2 = alfa**2
          gmax   = abs( alfa )              ! alpha1
          gmin   = gmax                     ! alpha1
       end if

       ! Apply previous rotation Qk-1 to get
       !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
       !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

       oldeps = epsln
       delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
       gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
       epsln  =               sn * beta ! epsln2 = 0         epslnk+1
       dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

       ! Compute the next plane rotation Qk

       gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
       cs     = gbar / gamma                ! ck
       sn     = beta / gamma                ! sk
       phi    = cs * phibar                 ! phik
       phibar = sn * phibar                 ! phibark+1

       if (debug) then
          write(*,*) ' '
          write(*,*) 'v    ', v
          write(*,*) 'alfa ', alfa
          write(*,*) 'beta ', beta
          write(*,*) 'gamma', gamma
          write(*,*) 'delta', delta
          write(*,*) 'gbar ', gbar
          write(*,*) 'epsln', epsln
          write(*,*) 'dbar ', dbar
          write(*,*) 'phi  ', phi
          write(*,*) 'phiba', phibar
          write(*,*) ' '
       end if

       ! Update  x.

       denom = one/gamma

       do i = 1, n
          w1(i) = w2(i)
          w2(i) = w(i)
          w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
          x(i)  =   x(i) +   phi * w(i)
       end do

       ! Go round again.

       gmax   = max( gmax, gamma )
       gmin   = min( gmin, gamma )
       z      = rhs1 / gamma
       ynorm2 = z**2  +  ynorm2
       rhs1   = rhs2  -  delta * z
       rhs2   =       -  epsln * z

       ! Estimate various norms and test for convergence.

       Anorm  = sqrt( tnorm2 )
       ynorm  = sqrt( ynorm2 )
       epsa   = Anorm * eps
       epsx   = Anorm * ynorm * eps
       epsr   = Anorm * ynorm * rtol
       diag   = gbar
       if (diag == zero) diag = epsa

       qrnorm = phibar
       rnorml = rnorm
       rnorm  = qrnorm
       rootl       = sqrt( gbar**2 +dbar**2  )  ! norm([gbar; dbar]);
       Arnorml     = rnorml*rootl               ! ||A r_{k-1} ||
       relArnorml  = rootl  /  Anorm;           ! ||Ar|| / (||A|| ||r||)     
       !relArnorml = Arnorml / Anorm;           ! ||Ar|| / ||A|| 

       ! Estimate  cond(A).
       ! In this version we look at the diagonals of  R  in the
       ! factorization of the lower Hessenberg matrix,  Q * H = R,
       ! where H is the tridiagonal matrix from Lanczos with one
       ! extra row, beta(k+1) e_k^T.

       Acond  = gmax / gmin

       ! See if any of the stopping criteria are satisfied.
       ! In rare cases, istop is already -1 from above (Abar = const*I).

       if (istop == 0) then
          if (itn    >= itnlim    ) istop = 5
          if (Acond  >= 0.1d+0/eps) istop = 4
          if (epsx   >= beta1     ) istop = 3
          if (qrnorm <= epsx  .or.  relArnorml <= epsx) istop = 2
          if (qrnorm <= epsr  .or.  relArnorml <= epsr) istop = 1
          IF ( debug ) THEN
              IF ( istop == 1 ) THEN
                  WRITE(*,*) ' qrnorm=', qrnorm 
                  WRITE(*,*) ' relArnorml=', relArnorml
                  WRITE(*,*) ' rootl=', rootl
                  WRITE(*,*) ' gbar=', gbar, ' dbar=', dbar 
                  WRITE(*,*) ' epsr=', epsr
              ENDIF
          ENDIF
       ENDIF


       ! See if it is time to print something.

       if (nout > 0) then
          prnt   = .false.
          if (n      <= 40         ) prnt = .true.
          if (itn    <= 10         ) prnt = .true.
          if (itn    >= itnlim - 10) prnt = .true.
          if (mod(itn,10)  ==     0) prnt = .true.
          if (qrnorm <=  ten * epsx) prnt = .true.
          if (qrnorm <=  ten * epsr) prnt = .true.
          if (relArnorml<= ten*epsx) prnt = .true.
          if (relArnorml<= ten*epsr) prnt = .true.
          if (Acond  >= 1.0d-2/eps ) prnt = .true.
          if (istop  /=  0         ) prnt = .true.

          if ( prnt ) then
             if (    itn     == 1) write(nout, 1200)
             write(nout, 1300) itn, x(1), qrnorm, Anorm, Acond
             if (mod(itn,10) == 0) write(nout, 1500)
          end if
       end if
       if (istop /= 0) exit

    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

    ! Display final status.

900 Arnorm = Arnorml
    if (nout  > 0) then
       write(nout, 2000) exitt, istop, itn,   &
                         exitt, Anorm, Acond, &
                         exitt, rnorm, ynorm, Arnorm
       write(nout, 3000) exitt, msg(istop)
    end if

    return

 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'    &
              / ' n      =', i7, 5x, 'checkA =', l4, 12x,         &
                 'precon =', l4                                   &
              / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,       &
                 'shift  =', e23.14)
 1200 format(// 5x, 'itn', 8x, 'x(1)', 10x,                       &
                'norm(r)', 3x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, 3e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 5x, 'istop =', i3,   14x, 'itn   =', i8     &
             /     a, 5x, 'Anorm =', e12.4, 5x, 'Acond =', e12.4  &
             /     a, 5x, 'rnorm =', e12.4, 5x, 'ynorm =', e12.4, 5x, 'Arnorml =', e12.4)
 3000 format(      a, 5x, a )

  END SUBROUTINE MINRES

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File symmlq.f
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE symmlq( n, Aprod, Msolve, b, shift, checkA, precon,                    &
                       x, itnlim, nout, rtol, istop, itn, Anorm, Acond, rnorm, ynorm, &
                       mtz_1, map_1, pdb_2, mode, sk, Acsr )

      INTEGER,                               INTENT(IN)            :: n
      INTEGER,                               INTENT(IN)            :: itnlim
      INTEGER,                               INTENT(IN)            :: nout 
      LOGICAL,                               INTENT(IN)            :: checkA
      LOGICAL,                               INTENT(IN)            :: precon
      REAL(KIND=wp),                         INTENT(IN)            :: b(n)
      REAL(KIND=wp),                         INTENT(IN)            :: shift
      REAL(KIND=wp),                         INTENT(IN)            :: rtol
      REAL(KIND=wp),                         INTENT(INOUT)         :: x(n)
      INTEGER,                               INTENT(OUT)           :: istop
      INTEGER,                               INTENT(OUT)           :: itn
      REAL(KIND=wp),                         INTENT(OUT)           :: Anorm
      REAL(KIND=wp),                         INTENT(OUT)           :: Acond
      REAL(KIND=wp),                         INTENT(OUT)           :: rnorm
      REAL(KIND=wp),                         INTENT(OUT)           :: ynorm
!     X-ray data:
      TYPE(mtz),                             INTENT(INOUT)         :: mtz_1
      TYPE(map),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT)         :: map_1
      TYPE(pdb),                             INTENT(IN)            :: pdb_2
      CHARACTER(LEN=*),                      INTENT(IN)            :: mode
      REAL(KIND=wp), DIMENSION(:),           INTENT(IN)            :: sk            ! optional??
      TYPE(csr_matrix),                      INTENT(IN),  OPTIONAL :: Acsr
      INTERFACE
      SUBROUTINE Aprod (n, y, x, mtz_1, map_1, pdb_2, mode, sk, Acsr)                   ! y := A*x
         USE basic_pdb_manip,  ONLY : pdb 
         USE mtz_io
         USE select_kinds
         USE map_fft
         USE sparse_basic
         INTEGER,                               INTENT(IN)             :: n
         REAL(KIND=wp), DIMENSION(n),           INTENT(OUT)            :: y
         REAL(KIND=wp), DIMENSION(n),           INTENT(IN)             :: x
         TYPE(mtz),                             INTENT(INOUT)          :: mtz_1
         TYPE(map),  DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: map_1
         TYPE(pdb),                             INTENT(IN)             :: pdb_2
         CHARACTER(LEN=*),                      INTENT(IN)             :: mode
         REAL(KIND=wp), DIMENSION(:),           INTENT(IN),   OPTIONAL :: sk
         TYPE(csr_matrix),                      INTENT(IN),   OPTIONAL :: Acsr 
       END SUBROUTINE Aprod

       SUBROUTINE Msolve(n, y, b, sk, A)                   ! Solve M*y = x
         USE select_kinds
         USE sparse_basic
         INTEGER,                        INTENT(IN)           :: n
         REAL(KIND=wp),    DIMENSION(n), INTENT(OUT)          :: y
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: b
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: sk
         TYPE(csr_matrix),               INTENT(IN), OPTIONAL :: A
       END SUBROUTINE Msolve

    END INTERFACE

!     ------------------------------------------------------------------
!
!     SYMMLQ  is designed to solve the system of linear equations
!
!                Ax = b
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A is not required to be positive definite.
!     (If A is known to be definite, the method of conjugate gradients
!     might be preferred, since it will require about the same number of
!     iterations as SYMMLQ but slightly less work per iteration.)
!
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, SYMMLQ is designed to solve the system
!
!                (A - shift*I) x = b
!
!     where  shift  is a specified scalar value.  If  shift  and  b
!     are suitably chosen, the computed vector x may approximate an
!     (unnormalized) eigenvector of A, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     Again, the matrix (A - shift*I) need not be positive definite.
!     The work per iteration is very slightly less if  shift = 0.
!
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., SYMMLQ will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
!     for IBM mainframes and IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling SYMMLQ must declare
!                        Aprod and Msolve to be external.
!
!     checkA  input      If checkA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     goodb   input      Usually, goodb should be .false.
!                        If x is expected to contain a large multiple of
!                        b (as in Rayleigh-quotient iteration),
!                        better precision may result if goodb = .true.
!                        See Lewis (1977) below.
!                        When goodb = .true., an extra call to Msolve
!                        is required.
!
!     precon  input      If precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  SYMMLQ terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and precon = .false., SYMMLQ
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     istop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!                        
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        Acond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If istop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     Anorm   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     Acond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If precon is false,
!                        ynorm is an estimate of norm(x).
!
!
!
!     To change precision
!     -------------------
!
!     Alter the words
!            double precision,
!            daxpy, dcopy, ddot, dnrm2
!     to their single or double equivalents.
!     ------------------------------------------------------------------
!
!
!     This routine is an implementation of the algorithm described in
!     the following references:
!
!     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
!          Systems of Linear Equations,
!          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
!
!     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
!          Report STAN-CS-77-595, Computer Science Department,
!          Stanford University, Stanford, California, March 1977.
!
!     Applications of SYMMLQ and the theory of preconditioning
!     are described in the following references:
!
!     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
!          Type Methods to Eigenvalue Calculations,
!          in R. Vichnevetsky and R.S. Steplman (editors),
!          Advances in Computer Methods for Partial Differential
!          Equations -- I/II, IMACS, 1979, 167-173.
!
!     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
!          Generalized Eigenvalue Calculations,
!          Ph. D. dissertation, Department of Mathematics,
!          New York University, New York, October 1983.
!
!     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
!          Preconditioners for indefinite systems arising in
!          optimization, SIMAX 13, 1, 292--311, January 1992.
!          (SIAM J. on Matrix Analysis and Applications)
!     ------------------------------------------------------------------
!
!
!     SYMMLQ development:
!            1972: First version.
!            1975: John Lewis recommended modifications to help with
!                  inverse iteration:
!                  1. Reorthogonalize v1 and v2.
!                  2. Regard the solution as x = x1  +  bstep * b,
!                     with x1 and bstep accumulated separately
!                     and bstep * b added at the end.
!                     (In inverse iteration, b might be close to the
!                     required x already, so x1 may be a lot smaller
!                     than the multiple of b.)
!            1978: Daniel Szyld and Olof Widlund implemented the first
!                  form of preconditioning.
!                  This required both a solve and a multiply with M.
!            1979: Implemented present method for preconditioning.
!                  This requires only a solve with M.
!            1984: Sven Hammarling noted corrections to tnorm and x1lq.
!                  SYMMLQ added to NAG Fortran Library.
!     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
!     16 Feb 1989: First F77 version.
!
!     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
!                  if Abar = const*I.  istop = -1 added for this case.
!
!     01 Mar 1989: Hans Mittelmann observed premature termination on
!                  ( 1  1  1 )     (   )                   ( 1  1    )
!                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
!                  ( 1     1 )     (   )                   (    1  1 )
!                  T2 is exactly singular, so estimating cond(A) from
!                  the diagonals of Lbar is unsafe.  We now use
!                  L       or  Lbar         depending on whether
!                  lqnorm  or  cgnorm       is least.
!
!     03 Mar 1989: eps computed internally instead of coming in as a
!                  parameter.
!     07 Jun 1989: ncheck added as a parameter to say if A and M
!                  should be checked for symmetry.
!                  Later changed to checkA (see below).
!     20 Nov 1990: goodb added as a parameter to make Lewis's changes
!                  an option.  Usually b is NOT much like x.  Setting
!                  goodb = .false. saves a call to Msolve at the end.
!     20 Nov 1990: Residual not computed exactly at end, to save time
!                  when only one or two iterations are required
!                  (e.g. if the preconditioner is very good).
!                  Beware, if precon is true, rnorm estimates the
!                  residual of the preconditioned system, not Ax = b.
!     04 Sep 1991: Parameter list changed and reordered.
!                  integer ncheck is now logical checkA.
!     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
!                  showed that beta2 = 0 (istop = -1) means that
!                  b is an eigenvector when M = I.
!                  More complicated if there is a preconditioner;
!                  not clear yet how to describe it.
!     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
!                  Need to estimate Anorm from column of Tk.
!
!     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
!     Department of EESOR                  mike@sol-michael.stanford.edu
!     Stanford University
!     Stanford, CA 94305-4023                             (650) 723-1875
!     ------------------------------------------------------------------
!
!
!     Subroutines and functions
!
!     USER       Aprod, Msolve
!     BLAS       daxpy, dcopy, ddot , dnrm2
!
!
!     Intrinsics and local variables
      LOGICAL        :: goodb = .FALSE.
      REAL(KIND=wp)  :: r1(n)
      REAL(KIND=wp)  :: r2(n)
      REAL(KIND=wp)  :: v(n)
      REAL(KIND=wp)  :: w(n)
      REAL(KIND=wp)  :: y(n)

      INTRINSIC      :: abs
      INTRINSIC      :: max
      INTRINSIC      :: min
      INTRINSIC      :: mod
      INTRINSIC      :: sqrt
      REAL(KIND=wp)  ::  ddot, dnrm2
      REAL(KIND=wp)  ::  alfa, b1, beta, beta1, bstep, cs, &
                         cgnorm, dbar, delta, denom, diag, &
                         eps, epsa, epsln, epsr, epsx, &
                         gamma, gbar, gmax, gmin, gpert, &
                         lqnorm, oldb, qrnorm, rhs1, rhs2, &
                         s, sn, snprod, t, tnorm, &
                         x1cg, x1lq, ynorm2, zbar, z
      integer            i

      REAL(KIND=wp), PARAMETER :: zero = 0.0_wp,  one = 1.0_wp,  two = 2.0_wp

      CHARACTER(LEN=16) :: enter, exit
      CHARACTER(LEN=52) :: msg(-1:8)

      DATA               enter /' Enter SYMMLQ.  '/, &
                         exit  /' Exit  SYMMLQ.  '/

      DATA               msg                                        &
       / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',     &
         'beta1 = 0.  The exact solution is  x = 0',                &
         'Requested accuracy achieved, as determined by rtol',      &
         'Reasonable accuracy achieved, given eps',                 &
         'x has converged to an eigenvector',                       &
         'Acond has exceeded 0.1/eps',                              &
         'The iteration limit was reached',                         &
         'Aprod  does not define a symmetric matrix',               &
         'Msolve does not define a symmetric matrix',               &
         'Msolve does not define a pos-def preconditioner' /
!     ------------------------------------------------------------------

!     Compute eps, the machine precision.  The call to daxpy is
!     intended to fool compilers that use extra-length registers.

      eps    = one / 16.0_wp

   10 eps    = eps / two
      x(1)   = eps
      y(1)   = one
      call daxpy ( 1, one, x, 1, y, 1 )
      if (y(1) .gt. one) go to 10

      eps    = eps * two

!     Print heading and initialize.

      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, goodb, precon,&
                           itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero

      do 50 i = 1, n
         x(i) = zero
   50 continue

!     Set up y for the first Lanczos vector v1.
!     y is really beta1 * P * v1  where  P = C**(-1).
!     y and beta1 will be zero if b = 0.

      CALL dcopy ( n, b, 1, y , 1 )
      CALL dcopy ( n, b, 1, r1, 1 )
!      if ( precon ) call Msolve( n, r1, y )

      IF ( precon ) THEN
          IF ( PRESENT ( Acsr ) ) THEN
!             Case when sparse matrix preconditioner is available:
              CALL Msolve( n, y, r1, sk, Acsr )
          ELSE
!             Just use main diag as preconditioner:
              CALL Msolve( n, y, r1, sk)
          ENDIF
      ENDIF

      IF ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      END IF
      beta1  = ddot  ( n, r1, 1, y, 1 )

!     See if Msolve is symmetric.

      if (checkA  .and.  precon) then

!         call Msolve( n, y, r2 )
         IF ( PRESENT ( Acsr ) ) THEN
             CALL Msolve ( n, r2, y, sk, Acsr )
         ELSE
!            Use main diagonal only for preconditioning:
             CALL Msolve( n, r2, y, sk )
         ENDIF

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = ABS ( s - t )
         epsa   = (s + eps) * eps ** 0.333333_wp

         WRITE(*,"(' Msolve check:', 7X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa

         if ( z > epsa ) then
             istop = 7
             go to 900
         end if
      end if

!     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         istop = 8
         go to 900
      end if

!     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

!     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1
      do 100 i = 1, n
         v(i)  = s * y(i)
  100 continue

!     See if Aprod  is symmetric.

!     call Aprod (  A, v, y )
      IF ( .NOT. precon ) THEN
!         Example of future modification:
!         IF ( PRESENT ( ACSR ) ) CALL Aprod (..., mode, ACSR )
          CALL Aprod( n, y, v, mtz_1, map_1, pdb_2, mode, sk, ACSR )
      ELSE
          CALL Aprod( n, y, v, mtz_1, map_1, pdb_2, mode )
      ENDIF

      if ( checkA ) then
!         call Aprod ( A, y, r2 )

         IF ( .NOT. precon ) THEN
             CALL Aprod( n, r2, y, mtz_1, map_1, pdb_2, mode, sk, ACSR )
         ELSE
             CALL Aprod( n, r2, y, mtz_1, map_1, pdb_2, mode)
         ENDIF

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333

         WRITE(*,"(21X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa

         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if

!     Set up y for the second Lanczos vector.
!     Again, y is beta * P * v2  where  P = C**(-1).
!     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

!     Make sure  r2  will be orthogonal to the first  v.

      z      = ddot  ( n, v, 1, y, 1 )
      s      = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
!      if ( precon ) call Msolve( n, r2, y )
      if ( precon ) THEN
          IF ( PRESENT ( ACSR ) ) THEN
              CALL Msolve( n, y, r2, sk, Acsr )  
          ELSE
              CALL Msolve( n, y, r2, sk )
          ENDIF
      ENDIF

      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 8
         WRITE(*,*) ' beta=', beta, ' istop=', istop
         go to 900
      end if

!     Cause termination (later) if beta is essentially zero.

      beta   = sqrt( beta )
      if (beta .le. eps) then
         istop = -1
      end if

!     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

!     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs ( alfa ) + eps
      gmin   = gmax

      if ( goodb ) then
         do 200 i = 1, n
            w(i)  = zero
  200    continue
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------

!     Estimate various norms and test for convergence.

  300 Anorm  = sqrt( tnorm  )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa

      lqnorm = sqrt( rhs1**2 + rhs2**2 )
      qrnorm = snprod * beta1
      cgnorm = qrnorm * beta / abs( diag )

!     Estimate  cond(A).
!     In this version we look at the diagonals of  L  in the
!     factorization of the tridiagonal matrix,  T = L*Q.
!     Sometimes, T(k) can be misleadingly ill-conditioned when
!     T(k+1) is not, so we must be careful not to overestimate Acond.

      if (lqnorm .le. cgnorm) then
         Acond  = gmax / gmin
      else
         denom  = min( gmin, abs( diag ) )
         Acond  = gmax / denom
      end if

!     See if any of the stopping criteria are satisfied.
!     In rare cases, istop is already -1 from above (Abar = const * I).

      if (istop .eq. 0) then
         if (itn    .ge. itnlim ) istop = 5
         if (Acond  .ge. 0.1/eps) istop = 4
         if (epsx   .ge. beta1  ) istop = 3
         if (cgnorm .le. epsx   ) istop = 2
         if (cgnorm .le. epsr   ) istop = 1
      end if
!     ==================================================================

!     See if it is time to print something.

      if (nout .le.  0)          go to 600
      if (n    .le. 40)          go to 400
      if (itn  .le. 10)          go to 400
      if (itn  .ge. itnlim - 10) go to 400
      if (mod(itn,10)  .eq.   0) go to 400
      if (cgnorm .le. 10.0*epsx) go to 400
      if (cgnorm .le. 10.0*epsr) go to 400
      if (Acond  .ge. 0.01/eps ) go to 400
      if (istop  .ne. 0)         go to 400
      go to 600

!     Print a line for this iteration.

  400 zbar   = rhs1 / diag
      z      = (snprod * zbar  +  bstep) / beta1
      x1lq   = x(1)  +  b1 * bstep / beta1
      x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

      if (    itn     .eq. 0) write(nout, 1200)
!      IF ( PRESENT ( sk ) ) THEN
!          write(nout, 1300) itn, x1cg*sk(1), cgnorm, bstep/beta1, Anorm, Acond
!      ELSE
          write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1, Anorm, Acond
!      ENDIF

      if (mod(itn,10) .eq. 0) write(nout, 1500)
!     ==================================================================


!     Obtain the current Lanczos vector  v = (1 / beta)*y
!     and set up  y  for the next iteration.

  600 if (istop /= 0) go to 800
      s  = one / beta

      do 620 i = 1, n
         v(i)  = s * y(i)
  620 continue

!     call Aprod ( a, v, y )
      IF ( .NOT. precon ) THEN
          CALL Aprod(n, y, v, mtz_1, map_1, pdb_2, mode, sk, ACSR)
      ELSE
          CALL Aprod(n, y, v, mtz_1, map_1, pdb_2, mode)
      ENDIF

      call daxpy ( n, (- shift), v, 1, y, 1 )
      call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
      alfa   = ddot( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
      call dcopy ( n, r2, 1, r1, 1 )
      call dcopy ( n, y, 1, r2, 1 )

!      if ( precon ) call Msolve( n, r2, y )
      IF ( precon ) THEN
          IF ( PRESENT ( ACSR ) ) THEN
              CALL Msolve( n, y, r2, sk, ACSR )
          ELSE
              CALL Msolve( n, y, r2, sk )
          ENDIF
      ENDIF

      oldb   = beta
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 6
         WRITE(*,*) ' beta=' , beta, ' istop=', istop, ' should be 8'
         go to 800
      end if

      beta   = sqrt( beta )
      tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

!     Compute the next plane rotation for  Q.

      gamma  = sqrt( gbar**2 + oldb**2 )
      cs     = gbar / gamma
      sn     = oldb / gamma
      delta  = cs * dbar  +  sn * alfa
      gbar   = sn * dbar  -  cs * alfa
      epsln  = sn * beta
      dbar   =            -  cs * beta

!     Update  x.

      z      = rhs1 / gamma
      s      = z * cs
      t      = z * sn

      do 700 i = 1, n
         x(i)  = (w(i) * s   +   v(i) * t)  +  x(i)
         w(i)  =  w(i) * sn  -   v(i) * cs
  700 continue

!     Accumulate the step along the direction  b,
!     and go round again.

      bstep  = snprod * cs * z  +  bstep
      snprod = snprod * sn
      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z
      itn    = itn   +  1
      go to 300

!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

!     Move to the CG point if it seems better.
!     In this version of SYMMLQ, the convergence tests involve
!     only cgnorm, so we're unlikely to stop at an LQ point,
!     EXCEPT if the iteration limit interferes.

  800 if (cgnorm .le. lqnorm) then
         zbar   = rhs1 / diag
         bstep  = snprod * zbar  +  bstep
         ynorm  = sqrt( ynorm2  +  zbar**2 )
         rnorm  = cgnorm
         call daxpy ( n, zbar, w, 1, x, 1 )
      else
         rnorm  = lqnorm
      end if

      if ( goodb ) then

!        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
!         if ( precon ) call Msolve( n, b, y )
         IF ( precon ) THEN
             IF ( PRESENT ( ACSR ) ) THEN
                 CALL Msolve( n, y, b, sk, ACSR )
             ELSE
                 CALL Msolve( n, y, b, sk )
             ENDIF
         ENDIF
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

!     ==================================================================
!     Display final status.
!     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,  &
                           exit, Anorm, Acond,&
                           exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return

!     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b' &
     &       / ' n      =', i7, 5x, 'checkA =', l4, 12x,       &
     &          'goodb  =', l4, 7x, 'precon =', l4             &
     &       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,     &
     &          'shift  =', e23.14)                          
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2   &
     &       / ' (v1,v2) before and after ', e14.2             &
     &       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,                  &
     &         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8     &
     &       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4  &
     &       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
!     ------------------------------------------------------------------
!     end of SYMMLQ
    END SUBROUTINE symmlq

    SUBROUTINE symmlq2( n, Aprod, Msolve, b, shift, checkA, precon,                   &
                       x, itnlim, nout, rtol, &
                       istop, itn, Anorm, Acond, rnorm, ynorm, &
                       sk, Acsr )
!
!     Purpose:
!     =======
!           
!     Version of symmlq for sparse matrices.
!
!
!
!     Slightly modified on March 12 2008 by B.Strokopytov
!
!
      INTEGER,                               INTENT(IN)            :: n
      INTEGER,                               INTENT(IN)            :: itnlim
      INTEGER,                               INTENT(IN)            :: nout 
      LOGICAL,                               INTENT(IN)            :: checkA
      LOGICAL,                               INTENT(IN)            :: precon
      REAL(KIND=wp),                         INTENT(IN)            :: b(n)
      REAL(KIND=wp),                         INTENT(IN)            :: shift
      REAL(KIND=wp),                         INTENT(IN)            :: rtol
      REAL(KIND=wp),                         INTENT(OUT)           :: x(n)
      INTEGER,                               INTENT(OUT)           :: istop
      INTEGER,                               INTENT(OUT)           :: itn
      REAL(KIND=wp),                         INTENT(OUT)           :: Anorm
      REAL(KIND=wp),                         INTENT(OUT)           :: Acond
      REAL(KIND=wp),                         INTENT(OUT)           :: rnorm
      REAL(KIND=wp),                         INTENT(OUT)           :: ynorm
      REAL(KIND=wp), DIMENSION(:),           INTENT(IN)            :: sk
      TYPE(csr_matrix),                      INTENT(IN)            :: Acsr

      INTERFACE

       SUBROUTINE Aprod (n, y, x, ACSR)                   ! y := A*x
         USE select_kinds
         USE sparse_basic
         INTEGER,  INTENT(IN)                                          :: n
         REAL(KIND=wp), DIMENSION(n),           INTENT(OUT)            :: y
         REAL(KIND=wp), DIMENSION(:),           INTENT(IN)             :: x
         TYPE(csr_matrix),                      INTENT(IN)             :: ACSR
       END SUBROUTINE Aprod

       SUBROUTINE Msolve(n, y, b, sk, A)                   ! Solve M*y = x
         USE select_kinds
         USE sparse_basic
         INTEGER,                        INTENT(IN)           :: n
         REAL(KIND=wp),    DIMENSION(n), INTENT(OUT)          :: y
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: b
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: sk
         TYPE(csr_matrix),               INTENT(IN), OPTIONAL :: A
       END SUBROUTINE Msolve

    END INTERFACE

!     ------------------------------------------------------------------
!
!     SYMMLQ  is designed to solve the system of linear equations
!
!                Ax = b
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A is not required to be positive definite.
!     (If A is known to be definite, the method of conjugate gradients
!     might be preferred, since it will require about the same number of
!     iterations as SYMMLQ but slightly less work per iteration.)
!
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, SYMMLQ is designed to solve the system
!
!                (A - shift*I) x = b
!
!     where  shift  is a specified scalar value.  If  shift  and  b
!     are suitably chosen, the computed vector x may approximate an
!     (unnormalized) eigenvector of A, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     Again, the matrix (A - shift*I) need not be positive definite.
!     The work per iteration is very slightly less if  shift = 0.
!
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., SYMMLQ will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
!     for IBM mainframes and IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling SYMMLQ must declare
!                        Aprod and Msolve to be external.
!
!     checkA  input      If checkA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     goodb   input      Usually, goodb should be .false.
!                        If x is expected to contain a large multiple of
!                        b (as in Rayleigh-quotient iteration),
!                        better precision may result if goodb = .true.
!                        See Lewis (1977) below.
!                        When goodb = .true., an extra call to Msolve
!                        is required.
!
!     precon  input      If precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  SYMMLQ terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and precon = .false., SYMMLQ
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     istop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!                        
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        Acond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If istop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     Anorm   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     Acond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If precon is false,
!                        ynorm is an estimate of norm(x).
!
!
!
!     To change precision
!     -------------------
!
!     Alter the words
!            double precision,
!            daxpy, dcopy, ddot, dnrm2
!     to their single or double equivalents.
!     ------------------------------------------------------------------
!
!
!     This routine is an implementation of the algorithm described in
!     the following references:
!
!     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
!          Systems of Linear Equations,
!          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
!
!     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
!          Report STAN-CS-77-595, Computer Science Department,
!          Stanford University, Stanford, California, March 1977.
!
!     Applications of SYMMLQ and the theory of preconditioning
!     are described in the following references:
!
!     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
!          Type Methods to Eigenvalue Calculations,
!          in R. Vichnevetsky and R.S. Steplman (editors),
!          Advances in Computer Methods for Partial Differential
!          Equations -- I/II, IMACS, 1979, 167-173.
!
!     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
!          Generalized Eigenvalue Calculations,
!          Ph. D. dissertation, Department of Mathematics,
!          New York University, New York, October 1983.
!
!     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
!          Preconditioners for indefinite systems arising in
!          optimization, SIMAX 13, 1, 292--311, January 1992.
!          (SIAM J. on Matrix Analysis and Applications)
!     ------------------------------------------------------------------
!
!
!     SYMMLQ development:
!            1972: First version.
!            1975: John Lewis recommended modifications to help with
!                  inverse iteration:
!                  1. Reorthogonalize v1 and v2.
!                  2. Regard the solution as x = x1  +  bstep * b,
!                     with x1 and bstep accumulated separately
!                     and bstep * b added at the end.
!                     (In inverse iteration, b might be close to the
!                     required x already, so x1 may be a lot smaller
!                     than the multiple of b.)
!            1978: Daniel Szyld and Olof Widlund implemented the first
!                  form of preconditioning.
!                  This required both a solve and a multiply with M.
!            1979: Implemented present method for preconditioning.
!                  This requires only a solve with M.
!            1984: Sven Hammarling noted corrections to tnorm and x1lq.
!                  SYMMLQ added to NAG Fortran Library.
!     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
!     16 Feb 1989: First F77 version.
!
!     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
!                  if Abar = const*I.  istop = -1 added for this case.
!
!     01 Mar 1989: Hans Mittelmann observed premature termination on
!                  ( 1  1  1 )     (   )                   ( 1  1    )
!                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
!                  ( 1     1 )     (   )                   (    1  1 )
!                  T2 is exactly singular, so estimating cond(A) from
!                  the diagonals of Lbar is unsafe.  We now use
!                  L       or  Lbar         depending on whether
!                  lqnorm  or  cgnorm       is least.
!
!     03 Mar 1989: eps computed internally instead of coming in as a
!                  parameter.
!     07 Jun 1989: ncheck added as a parameter to say if A and M
!                  should be checked for symmetry.
!                  Later changed to checkA (see below).
!     20 Nov 1990: goodb added as a parameter to make Lewis's changes
!                  an option.  Usually b is NOT much like x.  Setting
!                  goodb = .false. saves a call to Msolve at the end.
!     20 Nov 1990: Residual not computed exactly at end, to save time
!                  when only one or two iterations are required
!                  (e.g. if the preconditioner is very good).
!                  Beware, if precon is true, rnorm estimates the
!                  residual of the preconditioned system, not Ax = b.
!     04 Sep 1991: Parameter list changed and reordered.
!                  integer ncheck is now logical checkA.
!     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
!                  showed that beta2 = 0 (istop = -1) means that
!                  b is an eigenvector when M = I.
!                  More complicated if there is a preconditioner;
!                  not clear yet how to describe it.
!     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
!                  Need to estimate Anorm from column of Tk.
!
!     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
!     Department of EESOR                  mike@sol-michael.stanford.edu
!     Stanford University
!     Stanford, CA 94305-4023                             (650) 723-1875
!     ------------------------------------------------------------------
!
!
!     Subroutines and functions
!
!     USER       Aprod, Msolve
!     BLAS       daxpy, dcopy, ddot , dnrm2
!
!
!     Intrinsics and local variables
      LOGICAL        :: goodb = .FALSE.
      REAL(KIND=wp)  :: r1(n)
      REAL(KIND=wp)  :: r2(n)
      REAL(KIND=wp)  :: v(n)
      REAL(KIND=wp)  :: w(n)
      REAL(KIND=wp)  :: y(n)

      INTRINSIC      :: abs
      INTRINSIC      :: max
      INTRINSIC      :: min
      INTRINSIC      :: mod
      INTRINSIC      :: sqrt
      REAL(KIND=wp)  ::  ddot, dnrm2
      REAL(KIND=wp)  ::  alfa, b1, beta, beta1, bstep, cs, &
                         cgnorm, dbar, delta, denom, diag, &
                         eps, epsa, epsln, epsr, epsx, &
                         gamma, gbar, gmax, gmin, gpert, &
                         lqnorm, oldb, qrnorm, rhs1, rhs2, &
                         s, sn, snprod, t, tnorm, &
                         x1cg, x1lq, ynorm2, zbar, z
      integer            i

      REAL(KIND=wp), PARAMETER :: zero = 0.0_wp,  one = 1.0_wp,  two = 2.0_wp

      CHARACTER(LEN=16) :: enter, exit
      CHARACTER(LEN=52) :: msg(-1:8)

      DATA               enter /' Enter SYMMLQ.  '/, &
                         exit  /' Exit  SYMMLQ.  '/

      DATA               msg                                        &
       / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',     &
         'beta1 = 0.  The exact solution is  x = 0',                &
         'Requested accuracy achieved, as determined by rtol',      &
         'Reasonable accuracy achieved, given eps',                 &
         'x has converged to an eigenvector',                       &
         'Acond has exceeded 0.1/eps',                              &
         'The iteration limit was reached',                         &
         'Aprod  does not define a symmetric matrix',               &
         'Msolve does not define a symmetric matrix',               &
         'Msolve does not define a pos-def preconditioner' /
!     ------------------------------------------------------------------

!     Compute eps, the machine precision.  The call to daxpy is
!     intended to fool compilers that use extra-length registers.

      eps    = one / 16.0_wp

   10 eps    = eps / two
      x(1)   = eps
      y(1)   = one
      call daxpy ( 1, one, x, 1, y, 1 )
      if (y(1) .gt. one) go to 10

      eps    = eps * two

!     Print heading and initialize.

      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, goodb, precon,&
                           itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero

      do 50 i = 1, n
         x(i) = zero
   50 continue

!     Set up y for the first Lanczos vector v1.
!     y is really beta1 * P * v1  where  P = C**(-1).
!     y and beta1 will be zero if b = 0.

      CALL dcopy ( n, b, 1, y , 1 )
      CALL dcopy ( n, b, 1, r1, 1 )

!      if ( precon ) call Msolve( n, r1, y )

      IF ( precon ) THEN
          CALL Msolve( n, y, r1, sk)
      ENDIF

      IF ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      END IF

      beta1  = ddot  ( n, r1, 1, y, 1 )

!     See if Msolve is symmetric.

      if (checkA  .and.  precon) then

!         call Msolve( n, y, r2 )
         CALL Msolve( n, r2, y, sk )

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = ABS ( s - t )
         epsa   = (s + eps) * eps ** 0.333333_wp

         WRITE(*,"(' Msolve check:', 7X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa

         if ( z > epsa ) then
             istop = 7
             go to 900
         end if
      end if

!     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         istop = 8
         go to 900
      end if

!     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

!     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1
      do 100 i = 1, n
         v(i)  = s * y(i)
  100 continue

!     See if Aprod  is symmetric.

!     call Aprod (  A, v, y )
      CALL Aprod( n, y, v, ACSR )

      if ( checkA ) then
!         call Aprod ( A, y, r2 )

         CALL Aprod( n, r2, y, ACSR)

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333

         IF ( nout > 0 ) THEN
             WRITE(nout, "(21X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa
         ENDIF

         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if

!     Set up y for the second Lanczos vector.
!     Again, y is beta * P * v2  where  P = C**(-1).
!     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

!     Make sure  r2  will be orthogonal to the first  v.

      z  = ddot  ( n, v, 1, y, 1 )
      s  = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
!      if ( precon ) call Msolve( n, r2, y )
      IF ( precon ) THEN
          CALL Msolve( n, y, r2, sk )
      ENDIF

      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 8
         WRITE(*,*) ' beta=', beta, ' istop=', istop
         go to 900
      end if

!     Cause termination (later) if beta is essentially zero.

      beta   = SQRT ( beta )
      if (beta .le. eps) then
         istop = -1
      end if

!     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

!     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs ( alfa ) + eps
      gmin   = gmax

      if ( goodb ) then
         do 200 i = 1, n
            w(i)  = zero
  200    continue
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------

!     Estimate various norms and test for convergence.

  300 Anorm  = sqrt( tnorm  )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa

      lqnorm = sqrt( rhs1**2 + rhs2**2 )
      qrnorm = snprod * beta1
      cgnorm = qrnorm * beta / abs( diag )

!     Estimate  cond(A).
!     In this version we look at the diagonals of  L  in the
!     factorization of the tridiagonal matrix,  T = L*Q.
!     Sometimes, T(k) can be misleadingly ill-conditioned when
!     T(k+1) is not, so we must be careful not to overestimate Acond.

      if (lqnorm .le. cgnorm) then
         Acond  = gmax / gmin
      else
         denom  = min( gmin, abs( diag ) )
         Acond  = gmax / denom
      end if

!     See if any of the stopping criteria are satisfied.
!     In rare cases, istop is already -1 from above (Abar = const * I).

      if (istop .eq. 0) then
         if (itn    .ge. itnlim ) istop = 5
         if (Acond  .ge. 0.1/eps) istop = 4
         if (epsx   .ge. beta1  ) istop = 3
         if (cgnorm .le. epsx   ) istop = 2
         if (cgnorm .le. epsr   ) istop = 1
      end if
!     ==================================================================

!     See if it is time to print something.

      if (nout .le.  0)          go to 600
      if (n    .le. 40)          go to 400
      if (itn  .le. 10)          go to 400
      if (itn  .ge. itnlim - 10) go to 400
      if (mod(itn,10)  .eq.   0) go to 400
      if (cgnorm .le. 10.0*epsx) go to 400
      if (cgnorm .le. 10.0*epsr) go to 400
      if (Acond  .ge. 0.01/eps ) go to 400
      if (istop  .ne. 0)         go to 400
      go to 600

!     Print a line for this iteration.

  400 zbar   = rhs1 / diag
      z      = (snprod * zbar  +  bstep) / beta1
      x1lq   = x(1)  +  b1 * bstep / beta1
      x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

      if (    itn     .eq. 0) write(nout, 1200)
!      IF ( PRESENT ( sk ) ) THEN
!          write(nout, 1300) itn, x1cg*sk(1), cgnorm, bstep/beta1, Anorm, Acond
!      ELSE
          write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1, Anorm, Acond
!      ENDIF

      if (mod(itn,10) .eq. 0) write(nout, 1500)
!     ==================================================================


!     Obtain the current Lanczos vector  v = (1 / beta)*y
!     and set up  y  for the next iteration.

  600 if (istop /= 0) go to 800
      s  = one / beta

      do 620 i = 1, n
         v(i)  = s * y(i)
  620 continue

!     call Aprod ( a, v, y )
      CALL Aprod(n, y, v, ACSR)

      call daxpy ( n, (- shift), v, 1, y, 1 )
      call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
      alfa   = ddot( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
      call dcopy ( n, r2, 1, r1, 1 )
      call dcopy ( n, y, 1, r2, 1 )

!      if ( precon ) call Msolve( n, r2, y )
      IF ( precon ) THEN
          CALL Msolve( n, y, r2, sk )
      ENDIF

      oldb   = beta
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 6
         WRITE(*,*) ' beta=' , beta, ' istop=', istop, ' should be 8'
         go to 800
      end if

      beta   = sqrt( beta )
      tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

!     Compute the next plane rotation for  Q.

      gamma  = sqrt( gbar**2 + oldb**2 )
      cs     = gbar / gamma
      sn     = oldb / gamma
      delta  = cs * dbar  +  sn * alfa
      gbar   = sn * dbar  -  cs * alfa
      epsln  = sn * beta
      dbar   =            -  cs * beta

!     Update  x.

      z      = rhs1 / gamma
      s      = z * cs
      t      = z * sn

      do 700 i = 1, n
         x(i)  = (w(i) * s   +   v(i) * t)  +  x(i)
         w(i)  =  w(i) * sn  -   v(i) * cs
  700 continue

!     Accumulate the step along the direction  b,
!     and go round again.

      bstep  = snprod * cs * z  +  bstep
      snprod = snprod * sn
      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z
      itn    = itn   +  1
      go to 300

!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

!     Move to the CG point if it seems better.
!     In this version of SYMMLQ, the convergence tests involve
!     only cgnorm, so we're unlikely to stop at an LQ point,
!     EXCEPT if the iteration limit interferes.

  800 if (cgnorm .le. lqnorm) then
         zbar   = rhs1 / diag
         bstep  = snprod * zbar  +  bstep
         ynorm  = sqrt( ynorm2  +  zbar**2 )
         rnorm  = cgnorm
         call daxpy ( n, zbar, w, 1, x, 1 )
      else
         rnorm  = lqnorm
      end if

      if ( goodb ) then

!        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
!         if ( precon ) call Msolve( n, b, y )
         IF ( precon ) THEN
             CALL Msolve( n, y, b, sk )
         ENDIF
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

!     ==================================================================
!     Display final status.
!     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,  &
                           exit, Anorm, Acond,&
                           exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return

!     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b' &
     &       / ' n      =', i7, 5x, 'checkA =', l4, 12x,       &
     &          'goodb  =', l4, 7x, 'precon =', l4             &
     &       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,     &
     &          'shift  =', e23.14)                          
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2   &
     &       / ' (v1,v2) before and after ', e14.2             &
     &       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,                  &
     &         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8     &
     &       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4  &
     &       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
!     ------------------------------------------------------------------
!     end of SYMMLQ
    END SUBROUTINE symmlq2

    SUBROUTINE symmlq_dense (n, Aprod, Msolve, b, shift, checkA, precon, &
                             x, itnlim, nout, rtol,                      &
                             istop, itn, Anorm, Acond, rnorm, ynorm,     &
                             sk, A )
!
!     Purpose:
!     =======
!           
!     Version of symmlq for dense matrices.
!
!
!
!     Slightly modified on March 17 2008 by B.Strokopytov
!
!
      INTEGER,                               INTENT(IN)            :: n
      INTEGER,                               INTENT(IN)            :: itnlim
      INTEGER,                               INTENT(IN)            :: nout 
      LOGICAL,                               INTENT(IN)            :: checkA
      LOGICAL,                               INTENT(IN)            :: precon
      REAL(KIND=wp),                         INTENT(IN)            :: b(n)
      REAL(KIND=wp),                         INTENT(IN)            :: shift
      REAL(KIND=wp),                         INTENT(IN)            :: rtol
      REAL(KIND=wp),                         INTENT(OUT)           :: x(n)
      INTEGER,                               INTENT(OUT)           :: istop
      INTEGER,                               INTENT(OUT)           :: itn
      REAL(KIND=wp),                         INTENT(OUT)           :: Anorm
      REAL(KIND=wp),                         INTENT(OUT)           :: Acond
      REAL(KIND=wp),                         INTENT(OUT)           :: rnorm
      REAL(KIND=wp),                         INTENT(OUT)           :: ynorm
      REAL(KIND=wp), DIMENSION(:),           INTENT(IN)            :: sk
      REAL(KIND=wp), DIMENSION(:,:),         INTENT(IN)            :: A

      INTERFACE

       SUBROUTINE Aprod (n, y, x, A)                      ! y := A*x
         USE select_kinds
         USE sparse_basic
         INTEGER,                               INTENT(IN)             :: n
         REAL(KIND=wp), DIMENSION(n),           INTENT(OUT)            :: y
         REAL(KIND=wp), DIMENSION(:),           INTENT(IN)             :: x
         REAL(KIND=wp), DIMENSION(:,:),         INTENT(IN)             :: A
       END SUBROUTINE Aprod

       SUBROUTINE Msolve(n, y, b, sk,A)                      ! Solve SK * y = x
         USE select_kinds
         USE sparse_basic
         INTEGER,                        INTENT(IN)           :: n
         REAL(KIND=wp),    DIMENSION(n), INTENT(OUT)          :: y
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: b
         REAL(KIND=wp),    DIMENSION(n), INTENT(IN)           :: sk
         TYPE(csr_matrix),                      OPTIONAL      :: A
       END SUBROUTINE Msolve

    END INTERFACE

!     ------------------------------------------------------------------
!
!     SYMMLQ  is designed to solve the system of linear equations
!
!                Ax = b
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A is not required to be positive definite.
!     (If A is known to be definite, the method of conjugate gradients
!     might be preferred, since it will require about the same number of
!     iterations as SYMMLQ but slightly less work per iteration.)
!
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, SYMMLQ is designed to solve the system
!
!                (A - shift*I) x = b
!
!     where  shift  is a specified scalar value.  If  shift  and  b
!     are suitably chosen, the computed vector x may approximate an
!     (unnormalized) eigenvector of A, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     Again, the matrix (A - shift*I) need not be positive definite.
!     The work per iteration is very slightly less if  shift = 0.
!
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., SYMMLQ will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
!     for IBM mainframes and IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling SYMMLQ must declare
!                        Aprod and Msolve to be external.
!
!     checkA  input      If checkA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     goodb   input      Usually, goodb should be .false.
!                        If x is expected to contain a large multiple of
!                        b (as in Rayleigh-quotient iteration),
!                        better precision may result if goodb = .true.
!                        See Lewis (1977) below.
!                        When goodb = .true., an extra call to Msolve
!                        is required.
!
!     precon  input      If precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  SYMMLQ terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and precon = .false., SYMMLQ
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     istop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!                        
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        Acond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If istop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     Anorm   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     Acond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If precon is false,
!                        ynorm is an estimate of norm(x).
!
!
!
!     To change precision
!     -------------------
!
!     Alter the words
!            double precision,
!            daxpy, dcopy, ddot, dnrm2
!     to their single or double equivalents.
!     ------------------------------------------------------------------
!
!
!     This routine is an implementation of the algorithm described in
!     the following references:
!
!     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
!          Systems of Linear Equations,
!          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
!
!     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
!          Report STAN-CS-77-595, Computer Science Department,
!          Stanford University, Stanford, California, March 1977.
!
!     Applications of SYMMLQ and the theory of preconditioning
!     are described in the following references:
!
!     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
!          Type Methods to Eigenvalue Calculations,
!          in R. Vichnevetsky and R.S. Steplman (editors),
!          Advances in Computer Methods for Partial Differential
!          Equations -- I/II, IMACS, 1979, 167-173.
!
!     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
!          Generalized Eigenvalue Calculations,
!          Ph. D. dissertation, Department of Mathematics,
!          New York University, New York, October 1983.
!
!     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
!          Preconditioners for indefinite systems arising in
!          optimization, SIMAX 13, 1, 292--311, January 1992.
!          (SIAM J. on Matrix Analysis and Applications)
!     ------------------------------------------------------------------
!
!
!     SYMMLQ development:
!            1972: First version.
!            1975: John Lewis recommended modifications to help with
!                  inverse iteration:
!                  1. Reorthogonalize v1 and v2.
!                  2. Regard the solution as x = x1  +  bstep * b,
!                     with x1 and bstep accumulated separately
!                     and bstep * b added at the end.
!                     (In inverse iteration, b might be close to the
!                     required x already, so x1 may be a lot smaller
!                     than the multiple of b.)
!            1978: Daniel Szyld and Olof Widlund implemented the first
!                  form of preconditioning.
!                  This required both a solve and a multiply with M.
!            1979: Implemented present method for preconditioning.
!                  This requires only a solve with M.
!            1984: Sven Hammarling noted corrections to tnorm and x1lq.
!                  SYMMLQ added to NAG Fortran Library.
!     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
!     16 Feb 1989: First F77 version.
!
!     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
!                  if Abar = const*I.  istop = -1 added for this case.
!
!     01 Mar 1989: Hans Mittelmann observed premature termination on
!                  ( 1  1  1 )     (   )                   ( 1  1    )
!                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
!                  ( 1     1 )     (   )                   (    1  1 )
!                  T2 is exactly singular, so estimating cond(A) from
!                  the diagonals of Lbar is unsafe.  We now use
!                  L       or  Lbar         depending on whether
!                  lqnorm  or  cgnorm       is least.
!
!     03 Mar 1989: eps computed internally instead of coming in as a
!                  parameter.
!     07 Jun 1989: ncheck added as a parameter to say if A and M
!                  should be checked for symmetry.
!                  Later changed to checkA (see below).
!     20 Nov 1990: goodb added as a parameter to make Lewis's changes
!                  an option.  Usually b is NOT much like x.  Setting
!                  goodb = .false. saves a call to Msolve at the end.
!     20 Nov 1990: Residual not computed exactly at end, to save time
!                  when only one or two iterations are required
!                  (e.g. if the preconditioner is very good).
!                  Beware, if precon is true, rnorm estimates the
!                  residual of the preconditioned system, not Ax = b.
!     04 Sep 1991: Parameter list changed and reordered.
!                  integer ncheck is now logical checkA.
!     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
!                  showed that beta2 = 0 (istop = -1) means that
!                  b is an eigenvector when M = I.
!                  More complicated if there is a preconditioner;
!                  not clear yet how to describe it.
!     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
!                  Need to estimate Anorm from column of Tk.
!
!     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
!     Department of EESOR                  mike@sol-michael.stanford.edu
!     Stanford University
!     Stanford, CA 94305-4023                             (650) 723-1875
!     ------------------------------------------------------------------
!
!
!     Subroutines and functions
!
!     USER       Aprod, Msolve
!     BLAS       daxpy, dcopy, ddot , dnrm2
!
!
!     Intrinsics and local variables
      LOGICAL        :: goodb 
      REAL(KIND=wp)  :: r1(n)
      REAL(KIND=wp)  :: r2(n)
      REAL(KIND=wp)  :: v(n)
      REAL(KIND=wp)  :: w(n)
      REAL(KIND=wp)  :: y(n)

      INTRINSIC      :: abs
      INTRINSIC      :: max
      INTRINSIC      :: min
      INTRINSIC      :: mod
      INTRINSIC      :: sqrt
      REAL(KIND=wp)  ::  ddot, dnrm2
      REAL(KIND=wp)  ::  alfa, b1, beta, beta1, bstep, cs, &
                         cgnorm, dbar, delta, denom, diag, &
                         eps, epsa, epsln, epsr, epsx, &
                         gamma, gbar, gmax, gmin, gpert, &
                         lqnorm, oldb, qrnorm, rhs1, rhs2, &
                         s, sn, snprod, t, tnorm, &
                         x1cg, x1lq, ynorm2, zbar, z
      integer            i

      REAL(KIND=wp), PARAMETER :: zero = 0.0_wp,  one = 1.0_wp,  two = 2.0_wp

      CHARACTER(LEN=16) :: enter, exit
      CHARACTER(LEN=52) :: msg(-1:8)

      DATA               enter /' Enter SYMMLQ.  '/, &
                         exit  /' Exit  SYMMLQ.  '/

      DATA               msg                                        &
       / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',     &
         'beta1 = 0.  The exact solution is  x = 0',                &
         'Requested accuracy achieved, as determined by rtol',      &
         'Reasonable accuracy achieved, given eps',                 &
         'x has converged to an eigenvector',                       &
         'Acond has exceeded 0.1/eps',                              &
         'The iteration limit was reached',                         &
         'Aprod  does not define a symmetric matrix',               &
         'Msolve does not define a symmetric matrix',               &
         'Msolve does not define a pos-def preconditioner' /
!     ------------------------------------------------------------------

      goodb = n == 1

!     Compute eps, the machine precision.  The call to daxpy is
!     intended to fool compilers that use extra-length registers.

      eps    = one / 16.0_wp

   10 eps    = eps / two
      x(1)   = eps
      y(1)   = one
      call daxpy ( 1, one, x, 1, y, 1 )
      if (y(1) .gt. one) go to 10

      eps    = eps * two

!     Print heading and initialize.

      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, goodb, precon,&
                           itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero

      do 50 i = 1, n
         x(i) = zero
   50 continue

!     Set up y for the first Lanczos vector v1.
!     y is really beta1 * P * v1  where  P = C**(-1).
!     y and beta1 will be zero if b = 0.

      CALL dcopy ( n, b, 1, y , 1 )
      CALL dcopy ( n, b, 1, r1, 1 )

!      if ( precon ) call Msolve( n, r1, y )

      IF ( precon ) THEN
          CALL Msolve( n, y, r1, sk)
      ENDIF

      IF ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      END IF

      beta1  = ddot  ( n, r1, 1, y, 1 )

!     See if Msolve is symmetric.

      if (checkA  .and.  precon) then

!         call Msolve( n, y, r2 )
         CALL Msolve( n, r2, y, sk )

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = ABS ( s - t )
         epsa   = (s + eps) * eps ** 0.333333_wp

         IF ( nout > 0 ) THEN
             WRITE(nout, "(' Msolve check:', 7X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") &
             z, epsa
         ENDIF

         if ( z > epsa ) then
             istop = 7
             go to 900
         end if
      end if

!     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         istop = 8
         go to 900
      end if

!     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

!     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1
      do 100 i = 1, n
         v(i)  = s * y(i)
  100 continue

!     See if Aprod  is symmetric.

!     call Aprod (  A, v, y )
      CALL Aprod( n, y, v, A)

      if ( checkA ) then
!         call Aprod ( A, y, r2 )

         CALL Aprod( n, r2, y, A)

         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333

         IF ( nout > 0 ) THEN
             WRITE(nout, "(21X, 'z      =', ES11.2, 5X, 'epsa   =', ES11.2)") z, epsa
         ENDIF

         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if

!     Set up y for the second Lanczos vector.
!     Again, y is beta * P * v2  where  P = C**(-1).
!     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

!     Make sure  r2  will be orthogonal to the first  v.

      z  = ddot  ( n, v, 1, y, 1 )
      s  = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
!      if ( precon ) call Msolve( n, r2, y )
      IF ( precon ) THEN
          CALL Msolve( n, y, r2, sk )
      ENDIF

      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 8
         WRITE(*,*) ' beta=', beta, ' istop=', istop
         go to 900
      end if

!     Cause termination (later) if beta is essentially zero.

      beta   = SQRT ( beta )
      if (beta .le. eps) then
         istop = -1
      end if

!     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

!     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs ( alfa ) + eps
      gmin   = gmax

      if ( goodb ) then
         do 200 i = 1, n
            w(i)  = zero
  200    continue
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------

!     Estimate various norms and test for convergence.

  300 Anorm  = sqrt( tnorm  )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa

      lqnorm = sqrt( rhs1**2 + rhs2**2 )
      qrnorm = snprod * beta1
      cgnorm = qrnorm * beta / abs( diag )

!     Estimate  cond(A).
!     In this version we look at the diagonals of  L  in the
!     factorization of the tridiagonal matrix,  T = L*Q.
!     Sometimes, T(k) can be misleadingly ill-conditioned when
!     T(k+1) is not, so we must be careful not to overestimate Acond.

      if (lqnorm .le. cgnorm) then
         Acond  = gmax / gmin
      else
         denom  = min( gmin, abs( diag ) )
         Acond  = gmax / denom
      end if

!     See if any of the stopping criteria are satisfied.
!     In rare cases, istop is already -1 from above (Abar = const * I).

      if (istop .eq. 0) then
         if (itn    .ge. itnlim ) istop = 5
         if (Acond  .ge. 0.1/eps) istop = 4
         if (epsx   .ge. beta1  ) istop = 3
         if (cgnorm .le. epsx   ) istop = 2
         if (cgnorm .le. epsr   ) istop = 1
      end if
!     ==================================================================

!     See if it is time to print something.

      if (nout .le.  0)          go to 600
      if (n    .le. 40)          go to 400
      if (itn  .le. 10)          go to 400
      if (itn  .ge. itnlim - 10) go to 400
      if (mod(itn,10)  .eq.   0) go to 400
      if (cgnorm .le. 10.0*epsx) go to 400
      if (cgnorm .le. 10.0*epsr) go to 400
      if (Acond  .ge. 0.01/eps ) go to 400
      if (istop  .ne. 0)         go to 400
      go to 600

!     Print a line for this iteration.

  400 zbar   = rhs1 / diag
      z      = (snprod * zbar  +  bstep) / beta1
      x1lq   = x(1)  +  b1 * bstep / beta1
      x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

      if (    itn     .eq. 0) write(nout, 1200)
!      IF ( PRESENT ( sk ) ) THEN
!          write(nout, 1300) itn, x1cg*sk(1), cgnorm, bstep/beta1, Anorm, Acond
!      ELSE
          write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1, Anorm, Acond
!      ENDIF

      if (mod(itn,10) .eq. 0) write(nout, 1500)
!     ==================================================================


!     Obtain the current Lanczos vector  v = (1 / beta)*y
!     and set up  y  for the next iteration.

  600 if (istop /= 0) go to 800
      s  = one / beta

      do 620 i = 1, n
         v(i)  = s * y(i)
  620 continue

!     call Aprod ( a, v, y )
      CALL Aprod(n, y, v, A)

      call daxpy ( n, (- shift), v, 1, y, 1 )
      call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
      alfa   = ddot( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
      call dcopy ( n, r2, 1, r1, 1 )
      call dcopy ( n, y, 1, r2, 1 )

!      if ( precon ) call Msolve( n, r2, y )
      IF ( precon ) THEN
          CALL Msolve( n, y, r2, sk )
      ENDIF

      oldb   = beta
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         istop = 6
         WRITE(*,*) ' beta=' , beta, ' istop=', istop, ' should be 8'
         go to 800
      end if

      beta   = sqrt( beta )
      tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

!     Compute the next plane rotation for  Q.

      gamma  = sqrt( gbar**2 + oldb**2 )
      cs     = gbar / gamma
      sn     = oldb / gamma
      delta  = cs * dbar  +  sn * alfa
      gbar   = sn * dbar  -  cs * alfa
      epsln  = sn * beta
      dbar   =            -  cs * beta

!     Update  x.

      z      = rhs1 / gamma
      s      = z * cs
      t      = z * sn

      do 700 i = 1, n
         x(i)  = (w(i) * s   +   v(i) * t)  +  x(i)
         w(i)  =  w(i) * sn  -   v(i) * cs
  700 continue

!     Accumulate the step along the direction  b,
!     and go round again.

      bstep  = snprod * cs * z  +  bstep
      snprod = snprod * sn
      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z
      itn    = itn   +  1
      go to 300

!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

!     Move to the CG point if it seems better.
!     In this version of SYMMLQ, the convergence tests involve
!     only cgnorm, so we're unlikely to stop at an LQ point,
!     EXCEPT if the iteration limit interferes.

  800 if (cgnorm .le. lqnorm) then
         zbar   = rhs1 / diag
         bstep  = snprod * zbar  +  bstep
         ynorm  = sqrt( ynorm2  +  zbar**2 )
         rnorm  = cgnorm
         call daxpy ( n, zbar, w, 1, x, 1 )
      else
         rnorm  = lqnorm
      end if

      if ( goodb ) then

!        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
!         if ( precon ) call Msolve( n, b, y )
         IF ( precon ) THEN
             CALL Msolve( n, y, b, sk )
         ENDIF
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

!     ==================================================================
!     Display final status.
!     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,  &
                           exit, Anorm, Acond,&
                           exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return

!     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b' &
     &       / ' n      =', i7, 5x, 'checkA =', l4, 12x,       &
     &          'goodb  =', l4, 7x, 'precon =', l4             &
     &       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,     &
     &          'shift  =', e23.14)                          
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2   &
     &       / ' (v1,v2) before and after ', e14.2             &
     &       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,                  &
     &         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8     &
     &       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4  &
     &       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
!     ------------------------------------------------------------------
!     end of SYMMLQ
    END SUBROUTINE symmlq_dense

    SUBROUTINE Aprod_dense (n, y, x, A)                      ! y := A*x
!
!        Purpose:
!        =======
!        Shoud be used in conjuction with SYMMLQ_DENSE
!         
!
         USE select_kinds
         USE sparse_basic
         INTEGER,                       INTENT(IN)  :: n
         REAL(KIND=wp), DIMENSION(n),   INTENT(OUT) :: y
         REAL(KIND=wp), DIMENSION(:),   INTENT(IN)  :: x
         REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)  :: A

         y(1:n) = MATMUL(A(1:n,1:n),x(1:n))

    END SUBROUTINE Aprod_dense

    SUBROUTINE resmin_minres (n, y, x, sk, A )
!
!       Purpose:
!       =======
!       Solving My = x
!
!       Note:       
!       ====
!       MINRES seems unstable until now (MAR 2008).
!       Sometimes it stops just after one iteration.
!       Waiting for better version. As a shortcut
!       we use SYMMLQ2 which is modified version of 
!       original SYMMLQ.       
!
        INTEGER,                     INTENT(IN)  :: n
        REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: y
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: x
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: sk
        TYPE(csr_matrix),            INTENT(IN)  :: A
!       Local variables:
        REAL(KIND=wp)                            :: shift
        LOGICAL                                  :: precon
        LOGICAL                                  :: checkA
        INTEGER                                  :: itnlim
        INTEGER                                  :: nout
        REAL(KIND=wp)                            :: rtol
        INTEGER                                  :: istop
        INTEGER                                  :: itn
        REAL(KIND=wp)                            :: Anorm
        REAL(KIND=wp)                            :: Acond
        REAL(KIND=wp)                            :: rnorm
        REAL(KIND=wp)                            :: Arnorm
        REAL(KIND=wp)                            :: ynorm
        REAL(KIND=wp), DIMENSION(n)              :: b
!       Test:
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: r1
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: r2
!       LAPACK:
        REAL(KIND=wp), EXTERNAL                  :: dnrm2
        REAL(KIND=wp)                            :: r1norm
!       Counters:
        INTEGER                                  :: i
        CHARACTER(LEN=32)                        :: srname = 'resmin'
        
!       Use silent mode (no printing):
        nout = 0 

        shift = 0.0_wp
        IF ( lambda_shift > 0.0_wp ) THEN
            shift = -lambda_shift
            IF ( nout > 0 ) THEN
                WRITE(nout,*) ' Applying lambda shift=', shift
            ENDIF
        ENDIF

!        checkA = .TRUE.
        checkA = .FALSE.
        precon = .FALSE.

!       Well conditioned sparse matrix won't need too many itns, but just for precaution:
        itnlim = 5 * n


!       Get as accurate solution as possible:
        rtol = EPSILON ( 1.0_wp )

!       Scale RHS vector since matrix has been conditioned already:
        b = sk * x

!       Solve My = b:
        CALL minres (n, MATVEC, Diagonal_Msolve, b, shift, checkA, precon, &
                     y, itnlim, nout, rtol,                                &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm,       & ! ynorm is actually xnorm
                     sk, A)

        IF ( istop /= 1 .OR. itn < 2 .OR. itn > n) THEN
            WRITE(*,"(' RESMIN_MINRES> ', 'istop=', I2, ' itn=', I5)") istop, itn
            CALL warn('Ooops... Something went wrong in subroutine MINRES', srname)

!           Got unusual condition --> no silent mode anymore:
            nout = 6
        ENDIF

!       Check whether we need checking and printing:
        IF ( nout > 0 ) THEN

            CALL allocate_array ( r1, n )
            CALL allocate_array ( r2, n )

!           Check solution:
            CALL matvec ( n, r2, y, A)

!           b - Ax:
            r1 = b - r2

            r1norm = dnrm2 ( n, r1, 1)
            WRITE(*,*) ' r1norm(old)=', r1norm

!           Euclidean norm of a vector:
            r1norm = SQRT ( MAX ( DOT_PRODUCT ( r1, r1 ), 0.0_wp ) )

            DO i = 1, MIN ( 10, n)
                WRITE(nout, "(' b=', ES12.5, ' y=', ES12.5, ' diff= ', ES9.2)") &
                b(i), r2(i), r1(i)
            ENDDO

!           Print residual:
            WRITE(nout, 2000) r1norm
            2000 FORMAT(/ ' Final residual =', 1p, e8.1)

            CALL deallocate_array ( r1 )
            CALL deallocate_array ( r2 )

        ENDIF

!       Rescale the solution since ACSR has been explicitly conditioned:
        y = sk * y

    END SUBROUTINE resmin_minres

    SUBROUTINE resmin_symmlq (n, y, x, sk, A )
!       Solving My = x (Ay=x)
!
!       Note:
!       ====
!       At the moment we use preconditioned matrix A.
!       Therefore final solution MUST BE rescaled.
!
        INTEGER,                     INTENT(IN)  :: n
        REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: y
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: x
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: sk
        TYPE(csr_matrix),            INTENT(IN), OPTIONAL  :: A
!       Local variables:
        REAL(KIND=wp)                            :: shift
        LOGICAL                                  :: precon
        LOGICAL                                  :: checkA
        INTEGER                                  :: itnlim
        INTEGER                                  :: nout
        REAL(KIND=wp)                            :: rtol
        INTEGER                                  :: istop
        INTEGER                                  :: itn
        REAL(KIND=wp)                            :: Anorm
        REAL(KIND=wp)                            :: Acond
        REAL(KIND=wp)                            :: rnorm
        REAL(KIND=wp)                            :: ynorm
        REAL(KIND=wp), DIMENSION(n)              :: b
!       Test:
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: r1
        REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: r2
!       LAPACK:
        REAL(KIND=wp), EXTERNAL                  :: dnrm2
        REAL(KIND=wp)                            :: r1norm
!       Counters:
        INTEGER                                  :: i
        CHARACTER(LEN=32)                        :: srname = 'resmin_symmlq'
        
!       Use silent mode (no printing):
        nout = 0 

        shift = 0.0_wp
        IF ( lambda_shift > 0.0_wp ) THEN
            shift = -lambda_shift
            IF ( nout > 0 ) THEN
                WRITE(nout,*) ' Applying lambda shift=', shift
            ENDIF
        ENDIF

        checkA = .TRUE.
!        checkA = .FALSE.
        precon = .FALSE.

!       Well conditioned sparse matrix won't need too many itns, but just for precaution:
        itnlim = 5 * n

!       Get as accurate solution as possible:
        rtol = EPSILON ( 1.0_wp )

!       Scale RHS vector since matrix has been conditioned already:
        b = sk * x

!       Solve My = b:
        CALL symmlq2(n, MATVEC, Diagonal_Msolve, b, shift, checkA, precon, &
                     y, itnlim, nout, rtol,                                &
                     istop, itn, Anorm, Acond, rnorm, ynorm,               & ! ynorm is actually xnorm
                     sk, A)

        IF ( istop /= 1 .OR. itn < 2 .OR. itn > n) THEN
            WRITE(*,"(' RESMIN_SYMMLQ> ', 'istop=', I2, ' itn=', I5)") istop, itn
            CALL warn('Ooops... Something went wrong in subroutine SYMMLQ2', srname)

!           Got unusual condition --> no silent mode anymore:
            nout = 6
        ENDIF

!       Check whether we need checking and printing:
        IF ( nout > 0 ) THEN

            CALL allocate_array ( r1, n )
            CALL allocate_array ( r2, n )

!           Check solution:
            CALL matvec ( n, r2, y, A)

!           b - Ax:
            r1 = b - r2

            r1norm = dnrm2 ( n, r1, 1)
            WRITE(*,*) ' r1norm(old)=', r1norm

!           Euclidean norm of a vector:
            r1norm = SQRT ( MAX ( DOT_PRODUCT ( r1, r1 ), 0.0_wp ) )

            DO i = 1, MIN ( 10, n)
                WRITE(nout, "(' b=', ES12.5, ' y=', ES12.5, ' diff= ', ES9.2)") &
                b(i), r2(i), r1(i)
            ENDDO

!           Print residual:
            WRITE(nout, 2000) r1norm
            2000 FORMAT(/ ' Final residual =', 1p, e8.1)

            CALL deallocate_array ( r1 )
            CALL deallocate_array ( r2 )

        ENDIF

!       Rescale the solution since ACSR has been explicitly conditioned:
        y = sk * y

    END SUBROUTINE resmin_symmlq

    SUBROUTINE resmin_inv (n, y, x, sk, A )
!       Solving My = x (Ay=x)
!
!       Note:
!       ====
!       At the moment we use preconditioned matrix A.
!       Therefore final solution MUST BE rescaled.
!
        INTEGER,                     INTENT(IN)  :: n
        REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: y
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: x
! Unnecessary FIXME:
        REAL(KIND=wp), DIMENSION(n), INTENT(IN)  :: sk
        TYPE(csr_matrix),            INTENT(IN), OPTIONAL  :: A

        CALL matvec ( n, y, x, A)

    END SUBROUTINE resmin_inv

    SUBROUTINE Diagonal_Msolve (n, y, x, just_diag, A)
!
!        Purpose:
!        =======
!        Matrix conditioning by:
!
!        Solving My = x
         INTEGER,       INTENT(IN)    :: n
         REAL(KIND=wp), INTENT(IN)    :: x(n)
         REAL(KIND=wp), INTENT(OUT)   :: y(n)
         REAL(KIND=wp), INTENT(IN)    :: just_diag(n)
         TYPE(csr_matrix), INTENT(IN), OPTIONAL :: A
!        Vector X cannot and must not be altered:
         y = x / just_diag
    END SUBROUTINE Diagonal_Msolve

END MODULE implicit_matrix_inversion
