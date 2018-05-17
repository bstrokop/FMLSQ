MODULE sparse
USE fail
USE select_kinds
IMPLICIT NONE
CONTAINS
    SUBROUTINE coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )
!
!   Purpose:
!   =======
!   COOCSR converts COO to CSR.
!
!   Discussion:
!   ==========
!   This routine converts a matrix that is stored in COO coordinate format
!   a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!   Date:               Programmer:          History of changes:
!   ====                ==========           ==================
!   07 January 2004     Youcef Saad          Original code
!   03 March   2008     Boris Strokopytov    Adapted for use as a SPARSE module routine
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NNZ, the number of nonzero elements in the matrix.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on RETURN:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
        INTEGER,                          INTENT(IN)    :: nrow
        INTEGER(KIND=eb),                 INTENT(IN)    :: nnz
        REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
        INTEGER,       DIMENSION(:),      INTENT(IN)    :: ir
        INTEGER,       DIMENSION(:),      INTENT(IN)    :: jc
        REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: ao
        INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: iao
        INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jao
!       Local variables:
        INTEGER                                         :: i
        INTEGER                                         :: iad
        INTEGER                                         :: j
        INTEGER                                         :: k
        INTEGER                                         :: k0
        REAL(KIND=wp)                                   :: x

        iao(1:nrow+1) = 0
!
!       Determine the row lengths.
!
        DO k = 1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        ENDDO

!       The starting position of each row:
        k = 1
        DO j = 1, nrow+1
            k0 = iao(j)
            iao(j) = k
            k = k + k0
        ENDDO

!       Go through the structure once more.  Fill in output matrix:
        DO k = 1, nnz
            i = ir(k)
            j = jc(k)
            x = a(k)
            iad = iao(i)
            ao(iad) = x
            jao(iad) = j
            iao(i) = iad + 1
        ENDDO

!       Shift back IAO:
        DO j = nrow, 1, -1
            iao(j+1) = iao(j)
        ENDDO
        iao(1) = 1

    END SUBROUTINE coocsr

! FIXME: needs beautifying:
SUBROUTINE amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
                  iw, ierr )

!*****************************************************************************80
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  ModIFied:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix.
!
!    Input, INTEGER :: JOB, job indicator.  When JOB = 0, only the structure
!    is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, INTEGER :: NZMAX, the length of the arrays c and jc.
!    The routine will stop IF the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on RETURN:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = INTEGER ::. serving as error message.
!         ierr = 0 means normal RETURN,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = INTEGER :: work array of length equal to the number of
!         columns in A.
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: ncol
  INTEGER,                          INTENT(IN)    :: nzmax

  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: b
  REAL(KIND=wp), DIMENSION(nzmax),  INTENT(INOUT) :: c
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: jb
  INTEGER,       DIMENSION(ncol+1), INTENT(IN)    :: ib
  INTEGER,       DIMENSION(ncol+1), INTENT(INOUT) :: ic
  INTEGER,       DIMENSION(nzmax),  INTENT(INOUT) :: jc
! Work automatic array:
  INTEGER,       DIMENSION(ncol),   INTENT(INOUT) :: iw
  INTEGER,                          INTENT(OUT)   :: ierr
! Local variables and arrays:
  INTEGER                                         :: ii
  INTEGER                                         :: jcol
  INTEGER                                         :: jj
  INTEGER                                         :: job
  INTEGER                                         :: jpos
  INTEGER                                         :: k
  INTEGER                                         :: ka
  INTEGER                                         :: kb
  INTEGER                                         :: len
  REAL(KIND=wp)                                   :: scal
  LOGICAL                                         :: values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  DO ii = 1, nrow
!
!  Row I.
!
    DO ka = ia(ii), ia(ii+1)-1

      IF ( values ) THEN
        scal = a(ka)
      ENDIF

      jj = ja(ka)

      DO kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           IF ( jpos == 0 ) THEN
              len = len + 1
              IF ( nzmax < len ) THEN
                 ierr = ii
                 RETURN
              ENDIF
              jc(len) = jcol
              iw(jcol)= len
              IF ( values ) THEN
                c(len) = scal * b(kb)
              ENDIF
           else
              IF ( values ) THEN
                c(jpos) = c(jpos) + scal * b(kb)
              ENDIF
           ENDIF

         ENDDO

    ENDDO

    DO k = ic(ii), len
      iw(jc(k)) = 0
    ENDDO

    ic(ii+1) = len + 1

  ENDDO

END SUBROUTINE amub

SUBROUTINE amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the DOt product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  INTEGER, INTENT(IN) :: n

  REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: a
  REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: x
  REAL(KIND=wp), DIMENSION(:), INTENT(OUT) :: y(n)
  INTEGER,       DIMENSION(:), INTENT(IN)  :: ia
  INTEGER,       DIMENSION(:), INTENT(IN)  :: ja
! Local variables:
  REAL(KIND=wp)                            :: t
! Counters:
  INTEGER                                  :: i
  INTEGER                                  :: k

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(n, ia, ja, a, y, x)
  DO i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0_wp
    DO k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    ENDDO

    y(i) = t

  ENDDO

  RETURN
END SUBROUTINE amux

SUBROUTINE aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
                  iw, ierr )

!*****************************************************************************80
!
!! APLB performs the CSR matrix sum C = A + B.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of A and B.
!
!    Input, INTEGER :: NCOL, the column dimension of A and B.
!
!    Input, INTEGER :: JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = INTEGER ::. The  length of the arrays c and jc.
!         amub will stop IF the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on RETURN:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = INTEGER ::. serving as error message.
!         ierr = 0 means normal RETURN,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = INTEGER :: work array of length equal to the number of
!         columns in A.
!
  INTEGER,                          INTENT(IN)    :: ncol
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: job
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: b
  REAL(KIND=wp), DIMENSION(:),      INTENT(OUT)   :: c
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ib
  INTEGER,       DIMENSION(nrow+1), INTENT(OUT)   :: ic
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: jb
  INTEGER,       DIMENSION(:),      INTENT(OUT)   :: jc
  INTEGER,       DIMENSION(ncol),   INTENT(INOUT) :: iw
  INTEGER,                          INTENT(IN)    :: nzmax
  INTEGER,                          INTENT(OUT)   :: ierr
! Local variables:
  INTEGER                                         :: ii
  INTEGER                                         :: jcol
  INTEGER                                         :: jpos
  INTEGER                                         :: k
  INTEGER                                         :: ka
  INTEGER                                         :: kb
  INTEGER                                         :: len
  LOGICAL                                         :: values

  values = ( job /= 0 )
  ierr = 0
  len = 0
  ic(1) = 1
  iw(1:ncol) = 0

  DO ii = 1, nrow
!
!  Row I.
!
     DO ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        IF ( nzmax < len ) THEN
          ierr = ii
          RETURN
        ENDIF

        jc(len) = jcol
        IF ( values ) THEN
          c(len) = a(ka)
        ENDIF
        iw(jcol) = len
     ENDDO

     DO kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)

        IF ( jpos == 0 ) THEN

           len = len + 1

           IF ( nzmax < len ) THEN
             ierr = ii
             RETURN
           ENDIF

           jc(len) = jcol
           IF ( values ) THEN
             c(len) = b(kb)
           ENDIF
           iw(jcol)= len
        else
           IF ( values ) THEN
             c(jpos) = c(jpos) + b(kb)
           ENDIF
        ENDIF

     ENDDO

     DO k = ic(ii), len
       iw(jc(k)) = 0
     ENDDO

     ic(ii+1) = len+1
  ENDDO

  RETURN
END SUBROUTINE aplb

SUBROUTINE apldia ( nrow, job, a, ja, ia, diag, b, jb, ib, iw )

!*****************************************************************************80
!
!! APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in place (b, jb, ib, can be the same as
!    a, ja, ia, on entry). See comments for parameter job.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: JOB, job indicator. Job=0 means get array b only
!    (i.e. assume that a has already been copied into array b,
!    or that algorithm is used in place. ) For all practical
!    puposes enter job=0 for an in-place call and job=1 otherwise.
!    In case there are missing diagonal elements in A,
!    THEN the option job =0 will be ignored, since the algorithm
!    must modIFy the data structure (i.e. jb, ib) in this
!    situation.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), a diagonal matrix.
!
! on RETURN:
!
! b,
! jb,
! ib      = resulting matrix B in compressed sparse row sparse format.
!
!
! iw    = INTEGER :: work array of length n. On RETURN iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!

  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: job
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: b
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(IN)    :: diag(nrow)
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia(nrow+1)
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ib(nrow+1)
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jb
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: iw
  INTEGER                                         :: icount
  INTEGER                                         :: ii
  INTEGER                                         :: j
! Local variables:
  INTEGER                                         :: k
  INTEGER                                         :: k1
  INTEGER                                         :: k2
  INTEGER                                         :: ko
  INTEGER                                         :: nnz
  LOGICAL                                         :: test
!
!  Copy INTEGER :: arrays into B's data structure IF required.
!
  IF ( job /= 0 ) THEN
    nnz = ia(nrow+1)-1
    jb(1:nnz) = ja(1:nnz)
    ib(1:nrow+1) = ia(1:nrow+1)
  ENDIF
!
!  Get positions of diagonal elements in data structure.
!
  call diapos ( nrow, ja, ia, iw )
!
!  Count number of holes in diagonal and add DIAG elements to
!  valid diagonal entries.
!
  icount = 0

  DO j = 1, nrow

     IF ( iw(j) == 0 ) THEN
        icount = icount + 1
     else
        b(iw(j)) = a(iw(j)) + diag(j)
     ENDIF

  ENDDO
!
!  If no diagonal elements to insert, RETURN.
!
  IF ( icount == 0 ) THEN
    RETURN
  ENDIF
!
!  Shift the nonzero elements IF needed, to allow for created
!  diagonal elements.
!
  ko = ib(nrow+1) + icount
!
!  Copy rows backward.
!
  DO ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ib(ii)
     k2 = ib(ii+1) - 1
     ib(ii+1) = ko
     test = ( iw(ii) == 0 )

     DO k = k2, k1, -1

        j = jb(k)

         IF ( test .and. j < ii ) THEN
           test = .false.
           ko = ko - 1
           b(ko) = diag(ii)
           jb(ko) = ii
           iw(ii) = ko
        ENDIF

        ko = ko - 1
        b(ko) = a(k)
        jb(ko) = j

      ENDDO
!
!  The diagonal element has not been added yet.
!
     IF ( test ) THEN
        ko = ko - 1
        b(ko) = diag(ii)
        jb(ko) = ii
        iw(ii) = ko
     ENDIF

  ENDDO

  ib(1) = ko

END SUBROUTINE apldia

SUBROUTINE aplsb ( nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, nzmax, &
                   iw, ierr )

!*****************************************************************************80
!
!! APLSB performs the matrix linear combination C = A + s * B.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix B.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real S, scalar factor for B.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = INTEGER ::. The  length of the arrays c and jc.
!         amub will stop IF the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on RETURN:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = INTEGER ::. serving as error message.
!         ierr = 0 means normal RETURN,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = INTEGER :: work array of length equal to the number of
!         columns in A.
!
  INTEGER,                          INTENT(IN)    :: ncol
  INTEGER,                          INTENT(IN)    :: nrow
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia
  REAL(KIND=wp),                    INTENT(IN)    :: s
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: b
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: jb
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ib
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: c
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jc
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ic
  INTEGER,                          INTENT(IN)    :: nzmax
  INTEGER,       DIMENSION(ncol),   INTENT(INOUT) :: iw
  INTEGER,                          INTENT(OUT)   :: ierr
! Local variables:
  INTEGER                                         :: ii
  INTEGER                                         :: jcol
  INTEGER                                         :: jpos
  INTEGER                                         :: k
  INTEGER                                         :: ka
  INTEGER                                         :: kb
  INTEGER                                         :: len

  ierr = 0
  len = 0
  ic(1) = 1

  iw(1:ncol) = 0

  DO ii = 1, nrow
!
!  Row I.
!
     DO ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        IF ( nzmax < len ) THEN
          ierr = ii
          RETURN
        ENDIF

        jc(len) = jcol
        c(len) = a(ka)
        iw(jcol)= len

     ENDDO

     DO kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)
        IF ( jpos == 0 ) THEN
           len = len + 1
           IF ( nzmax < len ) THEN
             ierr = ii
             RETURN
           ENDIF
           jc(len) = jcol
           c(len) = s * b(kb)
           iw(jcol)= len
        else
           c(jpos) = c(jpos) + s * b(kb)
        ENDIF

     ENDDO

     DO k = ic(ii), len
        iw(jc(k)) = 0
     ENDDO

     ic(ii+1) = len + 1
 
  ENDDO

END SUBROUTINE aplsb

SUBROUTINE aplsca ( nrow, a, ja, ia, scal, iw )

!*****************************************************************************80
!
!! APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    important: the matrix A may be expanded slightly to allow for
!    additions of nonzero elements to previously nonexisting diagonals.
!    There is no checking as to whether there is enough space appENDed
!    to the arrays a and ja. IF not sure allow for n additional
!    elements.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real SCAL, a scalar to be added to the diagonal entries.
!
! on RETURN:
!
!
! a,
! ja,
! ia      = matrix A with diagonal elements shifted (or created).
!
! iw    = INTEGER :: work array of length n. On RETURN iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  INTEGER,                          INTENT(IN)    :: nrow
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia(nrow+1)
  REAL(KIND=wp),                    INTENT(IN)    :: scal
  INTEGER,       DIMENSION(nrow),   INTENT(INOUT) :: iw
! Local variables:
  INTEGER                                         :: icount
  INTEGER                                         :: ii
  INTEGER                                         :: j
  INTEGER                                         :: k
  INTEGER                                         :: k1
  INTEGER                                         :: k2
  INTEGER                                         :: ko
  LOGICAL                                         :: test

  call diapos ( nrow, ja, ia, iw )
  icount = 0

  DO j = 1, nrow

     IF ( iw(j) == 0 ) THEN
        icount = icount + 1
     else
        a(iw(j)) = a(iw(j)) + scal
     ENDIF

  ENDDO
!
!  If no diagonal elements to insert in data structure, RETURN.
!
  IF ( icount == 0 ) THEN
    RETURN
  ENDIF
!
!  Shift the nonzero elements IF needed, to allow for created
!  diagonal elements.
!
  ko = ia(nrow+1) + icount
!
!  Copy rows backward.
!
  DO ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     ia(ii+1) = ko
     test = ( iw(ii) == 0 )

     DO k = k2, k1, -1

        j = ja(k)

        IF ( test .AND. j < ii ) THEN
           test = .FALSE.
           ko = ko - 1
           a(ko) = scal
           ja(ko) = ii
           iw(ii) = ko
        ENDIF

        ko = ko - 1
        a(ko) = a(k)
        ja(ko) = j

    ENDDO
!
!  The diagonal element has not been added yet.
!
     IF ( test ) THEN
        ko = ko - 1
        a(ko) = scal
        ja(ko) = ii
        iw(ii) = ko
     ENDIF

  ENDDO

  ia(1) = ko

END SUBROUTINE aplsca

SUBROUTINE atmux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! ATMUX computes A' * x for a CSR matrix A.
!
!  Discussion:
!
!    This routine multiplies the transpose of a matrix by a vector when the
!    original matrix is stored in compressed sparse row storage. Can also be
!    viewed as the product of a matrix by a vector when the original
!    matrix is stored in the compressed sparse column format.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the row dimension of the matrix.
!
!    Input, real X(*), an array whose length is equal to the
!    column dimension of A.
!
!    Output, real Y(N), the product A' * X.
!
!    Input, real A(*), INTEGER :: JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER,                       INTENT(IN)  :: n
  REAL(KIND=wp), DIMENSION(:),   INTENT(IN)  :: x
  REAL(KIND=wp), DIMENSION(n),   INTENT(OUT) :: y
  REAL(KIND=wp), DIMENSION(:),   INTENT(IN)  :: a
  INTEGER,       DIMENSION(:),   INTENT(IN)  :: ja
  INTEGER,       DIMENSION(n+1), INTENT(IN)  :: ia
! Counters:
  INTEGER                                    :: i
  INTEGER                                    :: k

  y(1:n) = 0.0_wp

  DO i = 1, n
    DO k = ia(i), ia(i+1)-1
      y(ja(k)) = y(ja(k)) + x(i) * a(k)
    ENDDO
  ENDDO

END SUBROUTINE atmux

SUBROUTINE cnrms ( nrow, nrm, a, ja, ia, diag )

! Purpose:
! =======
! CNRMS gets the norms of each column of A.
!
!  Discussion:
!
!    There is a choice of three norms.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NRM, choosed the norm:
!    1, means 1-norm,
!    2, means the 2-nrm,
!    0, means max norm
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!   Output, REAL(KIND=wp) :: DIAG(NROW), the row norms.
!
  INTEGER,                          INTENT(IN)  :: nrow
  INTEGER,                          INTENT(IN)  :: nrm
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)  :: a
  INTEGER,       DIMENSION(:),      INTENT(IN)  :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)  :: ia
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(OUT) :: diag
! Counters:
  INTEGER                                       :: ii
  INTEGER                                       :: j
  INTEGER                                       :: k
  INTEGER                                       :: k1
  INTEGER                                       :: k2

  diag(1:nrow) = 0.0_wp

  DO ii = 1, nrow

     k1 = ia(ii)
     k2 = ia(ii+1) - 1

     DO k = k1, k2

        j = ja(k)
!
!  Update the norm of each column.
!
        IF ( nrm == 0 ) THEN
           diag(j) = MAX ( diag(j), ABS ( a(k) ) )
        else IF ( nrm == 1 ) THEN
           diag(j) = diag(j) + ABS ( a(k) )
        else
           diag(j) = diag(j) + a(k)**2
        ENDIF

    ENDDO

  ENDDO

  IF ( nrm /= 2 ) THEN
    RETURN
  ENDIF

  DO k = 1, nrow
      diag(k) = SQRT ( diag(k) )
  ENDDO

END SUBROUTINE cnrms

SUBROUTINE coocsr_inplace ( n, nnz, job, a, ja, ia, iwk )
!
!  Purpose:
!  =======
!  COOCSR_INPLACE converts COO to CSR in place.
!
!  Discussion:
!
!    This routine converts a matrix stored in coordinate format into
!    the CSR format.  The conversion is DOne in place in that the arrays
!    a,ja,ia of the result are overwritten onto the original arrays.
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order) use COOCSR
!    IF you want them sorted.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the row dimension of the matrix.
!
!    Input, INTEGER :: NNZ, the number of nonzero elements in A.
!
!    Input, INTEGER :: JOB.  When JOB = 1, the real values in A are
!    filled.  Otherwise A is not touched and the structure of the
!    array only (i.e. JA, IA)  is obtained.
!
!    Input/output, real A(NNZ).  On input, the matrix numeric values,
!    stored in the COO format.  On output, the numeric values, stored
!    in CSR format.
!
! ja      = INTEGER :: array of length nnz containing the column positions
!         of the corresponding elements in a.
!
! ia      = INTEGER :: array of length nnz containing the row positions
!         of the corresponding elements in a.
!
! iwk      = INTEGER :: work array of length n.
!
! on RETURN:
!
!    Output, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER, INTENT(IN)                          :: n
  INTEGER, INTENT(IN)                          :: nnz
  INTEGER, INTENT(IN)                          :: job

  REAL(KIND=wp), DIMENSION(:),   INTENT(INOUT) :: a
  INTEGER,       DIMENSION(nnz), INTENT(INOUT) :: ja(nnz)
  INTEGER,       DIMENSION(nnz), INTENT(INOUT) :: ia(nnz)
  INTEGER,       DIMENSION(n),   INTENT(INOUT) :: iwk
! Counters:
  INTEGER                                      :: i
  INTEGER                                      :: inext
  INTEGER                                      :: init
  INTEGER                                      :: ipos
  INTEGER                                      :: j
  INTEGER                                      :: jnext
  INTEGER                                      :: k
  REAL(KIND=wp)                                :: t
  REAL(KIND=wp)                                :: tnext
  LOGICAL                                      :: values

  values = (job == 1)
!
!  Find pointer array for resulting matrix.
!
  iwk(1:n+1) = 0

  DO k = 1, nnz
    i = ia(k)
    iwk(i+1) = iwk(i+1) + 1
  ENDDO

  iwk(1) = 1
  DO i = 2, n
    iwk(i) = iwk(i-1) + iwk(i)
  ENDDO
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    CONTINUE

  IF ( values ) THEN
    t = a(init)
  ENDIF

  i = ia(init)
  j = ja(init)
  ia(init) = -1

 6 CONTINUE
   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  ipos = iwk(i)
!
!  Save the chased element.
!
  IF ( values ) THEN
    tnext = a(ipos)
  ENDIF

  inext = ia(ipos)
  jnext = ja(ipos)
!
!  Then occupy its location.
!
  IF ( values ) THEN
    a(ipos) = t
  ENDIF

  ja(ipos) = j
!
!  Update pointer information for next element to come in row I.
!
  iwk(i) = ipos + 1
!
!  Determine the next element to be chased.
!
  IF ( ia(ipos) < 0 ) THEN
    GOTO 65
  ENDIF

  t = tnext
  i = inext
  j = jnext
  ia(ipos) = -1

  IF ( k < nnz ) THEN
    GOTO 6
  ENDIF

  GOTO 70

 65 CONTINUE

  init = init + 1

  IF ( nnz < init ) THEN
    GOTO 70
  ENDIF

  IF ( ia(init) < 0 ) THEN
    GOTO 65
  ENDIF
!
!  Restart chasing.
!
  GOTO 5

 70   CONTINUE

  ia(1) = 1
  ia(2:n+1) = iwk(1:n)

END SUBROUTINE coocsr_inplace 

SUBROUTINE copmat ( nrow, a, ja, ia, ao, jao, iao, ipos )

!*****************************************************************************80
!
!! COPMAT copies the matrix a, ja, ia, into the matrix ao, jao, iao.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, INTEGER :: IPOS, indicates the position in the array AO, JAO
!    where the first element should be copied.  Thus IAO(1) = IPOS on RETURN.
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the copied matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER, INTENT(IN)                        :: nrow

  REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: a
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: ao
  INTEGER, DIMENSION(:),       INTENT(IN)    :: ia
  INTEGER, DIMENSION(nrow+1),  INTENT(INOUT) :: iao
  INTEGER, DIMENSION(:),       INTENT(IN)    :: ja
  INTEGER, DIMENSION(:),       INTENT(INOUT) :: jao
  INTEGER,                     INTENT(IN)    :: ipos
! Local variables: 
  INTEGER                                    :: i
  INTEGER                                    :: k
  INTEGER                                    :: kst

  kst = ipos - ia(1)

  DO i = 1, nrow+1
    iao(i) = ia(i) + kst
  ENDDO

  DO k = ia(1), ia(nrow+1)-1
    ao(kst+k) = a(k)
    jao(kst+k)= ja(k)
  ENDDO

END SUBROUTINE copmat

SUBROUTINE cscal ( nrow, job, nrm, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! CSCAL scales the columns of A such that their norms are one.
!
!  Discussion:
!
!    The result matrix is written on B, or overwritten on A.
!
!    3 choices of norms: 1-norm, 2-norm, max-norm. in place.
!
!    The column dimension of A is not needed.
!
!    The algorithm in place (B can take the place of A).
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
! job   = INTEGER ::. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the INTEGER :: arrays ib, jb.
!
! nrm   = INTEGER ::. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on RETURN:
!
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the columns have been scaled, i.e., on RETURN
!        we have B = A * Diag
!
!    Output, real B(*), INTEGER :: JB(*), IB(NROW+1), the scaled matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: job
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: b
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(INOUT) :: diag(nrow)
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia(nrow+1)
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ib(nrow+1)
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jb

  INTEGER,                          INTENT(IN)    :: nrm

  call cnrms ( nrow, nrm, a, ja, ia, diag )

  diag(1:nrow) = 1.0D+00 / diag(1:nrow)

  call amudia ( nrow, job, a, ja, ia, diag, b, jb, ib )

END SUBROUTINE cscal

SUBROUTINE csrcoo ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )
!
! Purpose:
! =======
! CSRCOO converts Compressed Sparse Row to Coordinate format.
!
!  Discussion:
!
!   This routine converts a matrix that is stored in row general sparse 
!   A, JA, IA format into coordinate format AO, IR, JC. 
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
! job   = INTEGER :: serving as a job indicator.
!         IF job = 1 fill in only the array ir, ignore jc, and ao.
!         IF job = 2 fill in ir, and jc but not ao
!         IF job = 3 fill in everything.
!         The reason why these options are provided is that on RETURN
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         RETURNed. With job=1 only the array ir is RETURNed. Moreover,
!         the algorithm is in place:
!           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!         (Important: note the order in the output arrays a, ja, ia. )
!         i.e., ao can be the same as a, ir can be the same as ia
!         and jc can be the same as ja.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly IF the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on RETURN:
!-
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
!
! ierr       = INTEGER :: error indicator.
!         ierr == 0 means normal retur
!         ierr == 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
  INTEGER,                          INTENT(IN)    :: nrow
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia
  INTEGER,                          INTENT(OUT)   :: ierr
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: ao
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ir
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jc
  INTEGER,                          INTENT(IN)    :: job
  INTEGER,                          INTENT(OUT)   :: nnz
  INTEGER,                          INTENT(IN)    :: nzmax
! Local variables:
  INTEGER                                         :: i
  INTEGER                                         :: k
  INTEGER                                         :: k1
  INTEGER                                         :: k2

  ierr = 0
  nnz = ia(nrow+1) - 1

  IF ( nzmax < nnz ) THEN
    ierr = 1
    RETURN
  ENDIF

  IF ( 3 <= job ) THEN
    ao(1:nnz) = a(1:nnz)
  ENDIF

  IF ( 2 <= job ) THEN
    jc(1:nnz) = ja(1:nnz)
  ENDIF
!
!  Copy backward.
!
  DO i = nrow, 1, -1
    k1 = ia(i+1) - 1
    k2 = ia(i)
    DO k = k1, k2, -1
      ir(k) = i
    ENDDO
  ENDDO

END SUBROUTINE csrcoo

SUBROUTINE csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, INTEGER :: NDNS, the dimension of the DNS array.
!
!    Output, INTEGER :: IERR, error indicator.
!    0, means normal RETURN
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  INTEGER,                             INTENT(IN)    :: nrow 
  INTEGER,                             INTENT(IN)    :: ncol
  INTEGER,                             INTENT(IN)    :: ndns

  REAL(KIND=wp), DIMENSION(:),         INTENT(IN)    :: a
  INTEGER,       DIMENSION(:),         INTENT(IN)    :: ja
  INTEGER,       DIMENSION(:),         INTENT(IN)    :: ia
  REAL(KIND=wp), DIMENSION(ndns,ncol), INTENT(INOUT) :: dns
  INTEGER,                             INTENT(OUT)   :: ierr
! Counters:
  INTEGER                                            :: i
  INTEGER                                            :: j
  INTEGER                                            :: k
  
  ierr = 0
  dns(1:nrow,1:ncol) = 0.0_wp

  DO i = 1, nrow
    DO k = ia(i), ia(i+1)-1
      j = ja(k)
      IF ( ncol < j ) THEN
        ierr = i
        RETURN
      ENDIF
      dns(i,j) = a(k)
    ENDDO
  ENDDO

END SUBROUTINE csrdns

SUBROUTINE diamua ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
! Purpose:
! =======
! DIAMUA performs the matrix by matrix product B = Diag * A.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place; that is, B can take the place of A.
!    in this case use job=0.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: JOB, indicates the job to be DOne.
!    0, means get array B only;
!    1, means get B, and the INTEGER :: arrays IB and JB.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(N), a diagonal matrix stored as a vector.
!
!    Output, real B(*), INTEGER :: JB(*), INTEGER :: IB(NROW+1), the resulting 
!    matrix B in compressed sparse row sparse format.
!
  INTEGER, INTENT(IN) :: nrow

  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: b
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(IN)    :: diag
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ib
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jb
  INTEGER,                          INTENT(IN)    :: job
! Local variables:
  REAL(KIND=wp)                                   :: scal
! Counters:
  INTEGER                                         :: ii
  INTEGER                                         :: k
  INTEGER                                         :: k1
  INTEGER                                         :: k2

  DO ii = 1, nrow
!
!  Normalize each row.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     scal = diag(ii)
     b(k1:k2) = a(k1:k2) * scal

  ENDDO

  IF ( job == 0 ) THEN
    RETURN
  ENDIF

  ib(1) = ia(1)

  DO ii = 1, nrow
    ib(ii) = ia(ii)
    DO k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    ENDDO
  ENDDO

END SUBROUTINE diamua

SUBROUTINE diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
! Purpose:
! =======
!
! DIAPOS RETURNs the positions of the diagonal elements of a sparse matrix.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the row dimension of the matrix.
!
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, INTEGER :: IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  INTEGER,                 INTENT(IN)    :: n
  INTEGER, DIMENSION(:),   INTENT(IN)    :: ja
  INTEGER, DIMENSION(n+1), INTENT(IN)    :: ia
  INTEGER, DIMENSION(n),   INTENT(INOUT) :: idiag
! Counters:
  INTEGER                                :: i
  INTEGER                                :: k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  DO i = 1, n
    DO k = ia(i), ia(i+1) -1
      IF ( ja(k) == i ) THEN
        idiag(i) = k
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE diapos

SUBROUTINE dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
!
!
! Purpose:
! =======
! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine DOes not check whether an element is small.  It considers 
!    that A(I,J) is zero only IF it is exactly equal to zero.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix.
!
!    Input, INTEGER :: NZMAX, the maximum number of nonzero elements 
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, INTEGER :: NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, INTEGER :: IERR, error indicator.
!    0 means normal RETURN;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
  INTEGER,                             INTENT(IN)    :: ncol
  INTEGER,                             INTENT(IN)    :: ndns
  INTEGER,                             INTENT(IN)    :: nrow
  REAL(KIND=wp), DIMENSION(:),         INTENT(INOUT) :: a
  INTEGER,       DIMENSION(:),         INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(nrow+1),    INTENT(INOUT) :: ia
  REAL(KIND=wp), DIMENSION(ndns,ncol), INTENT(INOUT) :: dns
  INTEGER,                             INTENT(IN)    :: nzmax
  INTEGER,                             INTENT(OUT)   :: ierr
! Counters:
  INTEGER                                            :: i
  INTEGER                                            :: j
  INTEGER                                            :: next

  ierr = 0
  next = 1
  ia(1) = 1

  DO i = 1, nrow

    DO j = 1, ncol

      IF ( dns(i,j) /= 0.0_wp ) THEN

        IF ( nzmax < next ) THEN
          ierr = i
          RETURN
        ENDIF

        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1

      ENDIF

    ENDDO

    ia(i+1) = next

  ENDDO

END SUBROUTINE dnscsr

SUBROUTINE dscaldg ( n, a, ja, ia, diag, job )
!
! Purpose:
! =======
! DSCALDG scales rows by a diagonal factor.
!
!  Discussion:
!
!    This routine scales rows of a matrix by a diagonal factor DIAG.
!    DIAG is either given or to be computed.
!
!    If job = 1, we scale row I by by  +/- max |a(i,j) | and put the
!    inverse of the scaling factor in DIAG(i), where +/- is the sign of a(i,i).
!
!    If job = 2, we scale by the 2-norm of each row.
!
!    If DIAG(I) = 0, THEN DIAG(I) is replaced by 1.0.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the order of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, INTEGER :: JOB, describes the task to be performed.
!

  INTEGER,                       INTENT(IN)    :: n
  REAL(KIND=wp), DIMENSION(:),   INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(:),   INTENT(INOUT) :: diag
  INTEGER,       DIMENSION(:),   INTENT(IN)    :: ja
  INTEGER,       DIMENSION(n+1), INTENT(IN)    :: ia
  INTEGER,                       INTENT(IN)    :: job
! Local variables:
  REAL(KIND=wp)                                :: t
! Counters:
  INTEGER                                      :: i
  INTEGER                                      :: j
  INTEGER                                      :: k
  INTEGER                                      :: k1
  INTEGER                                      :: k2

  IF ( job == 2 ) THEN

    DO j = 1, n
      k1 = ia(j)
      k2 = ia(j+1) - 1
      t = 0.0_wp
      DO k = k1, k2
        t = t + a(k) * a(k)
      ENDDO
      diag(j) = SQRT ( t )
    ENDDO

  ELSEIF ( job == 1 ) THEN

    CALL retmx ( n, a, ja, ia, diag )

  ENDIF

   DO j = 1, n

     IF ( diag(j) /= 0.0_wp ) THEN
        diag(j) = 1.0_wp / diag(j)
     else
        diag(j) = 1.0_wp
     ENDIF

  ENDDO

  DO i = 1, n
    t = diag(i)
    DO k = ia(i), ia(i+1) -1
      a(k) = a(k) * t
    ENDDO
  ENDDO

END SUBROUTINE dscaldg

SUBROUTINE dump ( n, a, ja, ia, iout )

!*****************************************************************************80
!
!! DUMP writes the matrix to a file.
!
!  Discussion:
!
!    This routine writes the matrix to a file, one row at a time in a nice 
!    readable format.  This is a simple routine which is useful for debugging.
!
!    The output unit iout will have written in it the matrix in
!    one of two possible formats (depENDing on the max number of
!    elements per row. the values are output with only two digits
!    of accuracy (D9.2).
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the order of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, INTEGER :: IOUT, the FORTRAN output unit number.
!
  INTEGER, INTENT(IN) :: n

  REAL(KIND=wp), DIMENSION(:),   INTENT(IN) :: a
  INTEGER :: i
  INTEGER,       DIMENSION(:),   INTENT(IN) :: ja(*)
  INTEGER,       DIMENSION(n+1), INTENT(IN) :: ia
  INTEGER,                       INTENT(IN) :: iout
! Counters:
  INTEGER                                   :: k
  INTEGER                                   :: k1
  INTEGER                                   :: k2
  INTEGER                                   :: maxr
!
!  Select mode horizontal or vertical.
!
  maxr = 0
  DO i = 1, n
    maxr = MAX ( maxr, ia(i+1) - ia(i) )
  ENDDO

  IF ( maxr <= 8 ) THEN
!
!  Able to print one row across line.
!
    DO i = 1, n
      WRITE(iout,100) i
      k1 = ia(i)
      k2 = ia(i+1) - 1
      WRITE (iout,101) ja(k1:k2)
      WRITE (iout,102) a(k1:k2)
    ENDDO

  ELSE
!
!  Unable to print one row acros line.  Do three items at a time acros line.
!
     DO i = 1, n
        WRITE(iout,200) i
        k1 = ia(i)
        k2 = ia(i+1) - 1
        WRITE (iout,201) (ja(k),a(k), k = k1, k2)
    ENDDO

  ENDIF

 100  FORMAT(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  FORMAT(' col:',8(i5,6h     :))
 102  FORMAT(' val:',8(E9.2,2h :) )
 200  FORMAT(1h ,31(1h-),' row',i3,1x,31(1h-),/ &
         3('  columns :   values   *') )
 201  FORMAT(3(1h ,i5,6h    : ,D9.2,3h  *) )
END SUBROUTINE dump

SUBROUTINE getdia ( nrow, ncol, job, a, ja, ia, len, diag, idiag, ioff )
!
!*****************************************************************************80
!
!  Purpose:
!  =======
!  GETDIA extracts a given diagonal from a matrix stored in CSR format. 
!
!  Discussion:
!
!    The output matrix may be transformed with the diagonal removed
!    from it IF desired (as indicated by job.)
!
!    Our definition of a diagonal of matrix is a vector of length nrow
!    (always) which contains the elements in rows 1 to nrow of
!    the matrix that are contained in the diagonal offset by ioff
!    with respect to the main diagonal. If the diagonal element
!    falls outside the matrix THEN it is defined as a zero entry.
!    Thus the proper definition of diag(*) with offset ioff is
!
!    diag(k) = a(k,ioff+k) k = 1,2,...,nrow
!    with elements falling outside the matrix being defined as zero.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix.
!
! job   = INTEGER ::. Job indicator.  If job = 0 THEN
!         the matrix a, ja, ia, is not altered on RETURN.
!         IF job/=1  THEN getdia will remove the entries
!         collected in diag from the original matrix.
!         This is DOne in place.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! ioff  = INTEGER ::,containing the offset of the wanted diagonal
!        the diagonal extracted is the one corresponding to the
!        entries a(i,j) with j-i = ioff.
!        thus ioff = 0 means the main diagonal
!
! on RETURN:
!
! len   = number of nonzero elements found in diag.
!         (len <= min ( nrow, ncol-ioff ) - max ( 1, 1-ioff) + 1 )
!
! diag  = real array of length nrow containing the wanted diagonal.
!        diag contains the diagonal (a(i,j),j-i = ioff ) as defined
!         above.
!
! idiag = INTEGER :: array of  length len, containing the positions
!         in the original arrays a and ja of the diagonal elements
!         collected in diag. A zero entry in idiag(i) means that
!         there was no entry found in row i belonging to the diagonal.
!
! a, ja,
!    ia = IF job /= 0 the matrix is unchanged. otherwise the nonzero
!         diagonal entries collected in diag are removed from the
!         matrix. the structure is modIFied since the diagonal elements
!        are removed from a,ja,ia. Thus, the  RETURNed matrix will
!         have len fewer elements IF the diagonal is full.
!
  INTEGER,                     INTENT(IN)    :: ncol
  INTEGER,                     INTENT(IN)    :: nrow
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: a
  INTEGER,       DIMENSION(:), INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(:), INTENT(INOUT) :: ia
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: diag
  INTEGER :: i
  INTEGER,       DIMENSION(:), INTENT(INOUT) :: idiag
  INTEGER,                     INTENT(IN)    :: ioff
  INTEGER                                    :: istart
  INTEGER,                     INTENT(IN)    :: job
  INTEGER,                     INTENT(OUT)   :: len
! Counters:
  INTEGER                                    :: iEND
  INTEGER                                    :: k
  INTEGER                                    :: kdiag
  INTEGER                                    :: ko
  INTEGER                                    :: kold

  istart = MAX ( 0, -ioff )
  iEND   = MIN ( nrow, ncol-ioff )
  LEN    = 0
  idiag(1:nrow) = 0
  diag(1:nrow)  = 0.0_wp
!
!  Extract the diagonal elements.
!
  DO i = istart+1, iEND

     DO k = ia(i), ia(i+1) -1
        IF ( ja(k) - i == ioff ) THEN
           diag(i) = a(k)
           idiag(i) = k
           len = len + 1
           EXIT
        ENDIF
     ENDDO

  ENDDO

  IF ( job == 0 .or. len == 0 ) THEN
    RETURN
  ENDIF
!
!  Rewind the structure.
!
  ko = 0

  DO i = istart+1, iEND

    kold = ko
    kdiag = idiag(i)

    IF ( kdiag /= 0 ) THEN

      DO k = ia(i), ia(i+1)-1
        IF ( ja(k) /= kdiag ) THEN
          ko = ko + 1
          a(ko) = a(k)
          ja(ko) = ja(k)
        ENDIF
      ENDDO
      ia(i) = kold + 1
    ENDIF

  ENDDO
!
!  Redefine IA(NROW+1).
!
  ia(nrow+1) = ko + 1

END SUBROUTINE getdia

FUNCTION getelm ( i, j, a, ja, ia, iadd, sorted )

!*****************************************************************************80
!
! Purpose:
! =======
! GETELM RETURNs the element A(I,J) of a CSR matrix A. 
!
!  Discussion:
!
!    The matrix is assumed to be stored in Compressed Sparse Row (CSR) format.
!    This routine performs a binary search in the case where it is known 
!    that the elements are sorted, so that the column indices are in 
!    increasing order.  It also RETURNs, in IADD, the address of the 
!    element A(I,J) in arrays A and JA when the search is successsful.
!    IADD is 0 IF the element could not be found.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Noel Nachtigal, MIT
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: I, J, the row and column indices of the element.
!
!    Input, real A(*), INTEGER :: JA(*), IA(?+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, INTEGER :: IADD, the address of element A(I,J) in arrays A, JA 
!    IF found, zero IF not found.
!
!    Input, LOGICAL :: SORTED, is true IF the matrix is known to have its 
!    column indices sorted in increasing order.
!
!    Output, real GETELM, the value of A(I,J).
!
  REAL(KIND=wp)                            :: getelm
  REAL(KIND=wp), DIMENSION(:), INTENT(IN)  :: a
  INTEGER                                  :: i
  INTEGER,       DIMENSION(:), INTENT(IN)  :: ja
  INTEGER,       DIMENSION(:), INTENT(IN)  :: ia
  INTEGER,                     INTENT(OUT) :: iadd
! Local variables:
  LOGICAL :: sorted
! Counters:
  INTEGER :: ibeg
  INTEGER :: iEND
  INTEGER :: imid
  INTEGER :: j
  INTEGER :: k
!
!  Initialization.
!
  iadd = 0
  getelm = 0.0_wp
  ibeg = ia(i)
  iEND = ia(i+1) - 1
!
!  Case where the matrix is not necessarily sorted.
!
  IF ( .not. sorted ) THEN
!
!  Scan the row, and exit as soon as A(I,J) is found.
!
    DO k = ibeg, iEND
      IF ( ja(k) == j ) THEN
        iadd = k
        GOTO 20
      ENDIF
    ENDDO

  else
!
!  Begin binary search.  Compute the middle index.
!
 10 CONTINUE

   imid = ( ibeg + iEND ) / 2
!
!  Test IF found.
!
     IF ( ja(imid) == j ) THEN
        iadd = imid
        GOTO 20
     ENDIF

     IF ( iEND <= ibeg ) THEN
       GOTO 20
     ENDIF
!
!  else update the interval bounds.
!
     IF ( j < ja(imid) ) THEN
        iEND = imid - 1
     else
        ibeg = imid + 1
     ENDIF

     GOTO 10

  ENDIF

 20 CONTINUE

  IF ( iadd /= 0 ) THEN
    getelm = a(iadd)
  ENDIF

END FUNCTION getelm

SUBROUTINE ope ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!  Purpose:
!  =======
!  OPE sparse matrix * vector multiplication
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the order of the matrix.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real A(*), INTEGER :: JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER,                       INTENT(IN)    :: n
  REAL(KIND=wp), DIMENSION(n),   INTENT(IN)    :: x
  REAL(KIND=wp), DIMENSION(n),   INTENT(INOUT) :: y
  REAL(KIND=wp), DIMENSION(:),   INTENT(IN)    :: a
  INTEGER,       DIMENSION(n+1), INTENT(IN)    :: ia
  INTEGER,       DIMENSION(:),   INTENT(IN)    :: ja
! Counters:
  INTEGER                                      :: i
  INTEGER                                      :: k
! Pointers:
  INTEGER                                      :: k1
  INTEGER                                      :: k2

  DO i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) -1
    y(i) = 0.0_wp
    DO k = k1, k2
      y(i) = y(i) + a(k) * x(ja(k))
    ENDDO
  ENDDO

END SUBROUTINE ope

SUBROUTINE rnrms ( nrow, nrm, a, ja, ia, diag )

!*****************************************************************************80
!
!  Purpose:
!  =======
!  RNRMS gets the norms of each row of A.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
! nrm   = INTEGER ::. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), INTEGER ::, JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DIAG(NROW), the row norms.
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: nrm
  REAL(KIND=wp), DIMENSION(:),      INTENT(IN)    :: a
  INTEGER,       DIMENSION(:),      INTENT(IN)    :: ja
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(INOUT) :: diag
  INTEGER,       DIMENSION(nrow+1), INTENT(IN)    :: ia(nrow+1)
! Local variables:
  REAL(KIND=wp)                                   :: scal
! Counters:
  INTEGER                                         :: ii
  INTEGER                                         :: k
! Pointers:
  INTEGER                                         :: k1
  INTEGER                                         :: k2
!
!  Compute the norm of each element.
!
  DO ii = 1, nrow

    scal = 0.0_wp
    k1 = ia(ii)
    k2 = ia(ii+1) - 1

    IF ( nrm == 0 ) THEN
      DO k = k1, k2
        scal = MAX ( scal, ABS ( a(k) ) )
      ENDDO
    ELSE IF ( nrm == 1 ) THEN
      DO k = k1, k2
        scal = scal + ABS ( a(k) )
      ENDDO
    ELSE
      DO k = k1, k2
        scal = scal + a(k)**2
      ENDDO
    ENDIF

    IF ( nrm == 2 ) THEN
      scal = SQRT ( scal )
    ENDIF

    diag(ii) = scal

  ENDDO

END SUBROUTINE rnrms

SUBROUTINE rscal ( nrow, job, nrm, a, ja, ia, diag, b, jb, ib )
!
!  Purpose:
!  =======
!  RSCAL normalizes the rows of A.
!
!  Discussion:
!
!    There are three choices for the norm to use, the 1-norm, 2-norm 
!    or max-norm.
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place, so the A and B information can share
!    the same memory.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
! job   = INTEGER ::. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the INTEGER :: arrays ib, jb.
!
! nrm   = INTEGER ::. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on RETURN:
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the rows have been scaled, i.e., on RETURN
!        we have B = Diag*A.
!
!    Output, real B(*), INTEGER :: JB(*), IB(NROW+1), the scaled matrix in CSR
!    Compressed Sparse Row format.
!
  INTEGER, INTENT(IN) :: nrow

  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: b
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(INOUT) :: diag
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ib
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jb
  INTEGER,                          INTENT(IN)    :: job
  INTEGER,                          INTENT(IN)    :: nrm

  CALL rnrms ( nrow, nrm, a, ja, ia, diag )

  diag(1:nrow) = 1.0_wp / diag(1:nrow)

  CALL diamua ( nrow, job, a, ja, ia, diag, b, jb, ib )

END SUBROUTINE rscal

SUBROUTINE ssrcsr ( nrow, a, ja, ia, nzmax, ao, jao, iao, indu, ierr )

!*****************************************************************************80
!
!! SSRCSR converts Symmetric Sparse Row to (regular) Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a symmetric  matrix in which only the lower
!    part is  stored in compressed sparse row format, i.e.,
!    a matrix stored in symmetric sparse format, into a fully stored matrix
!    i.e., a matrix where both the lower and upper parts are stored in
!    compressed sparse row format. the algorithm is in place (i.e. result
!    may be overwritten onto the input matrix a, ja, ia ----- ).
!
!    the output matrix delivered by ssrcsr is such that each row starts with
!    the elements of the lower part followed by those of the upper part.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
! a,
! ia,
! ja    = matrix in compressed sparse row format. This is assumed to be
!         a lower triangular matrix.
!
! nzmax      = size of arrays ao and jao. ssrcsr will abort IF the storage
!         provided in a, ja is not sufficient to store A. See ierr.
!
! on RETURN:
!
! ao, iao,
!   jao = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A**T - D, IF
!         A is the original lower triangular matrix. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!
! indu  = INTEGER :: array of length nrow+1. If the input matrix is such
!         that the last element in each row is its diagonal element THEN
!         on RETURN, indu will contain the pointers to the diagonal
!         element in each row of the output matrix. Otherwise used as
!         work array.
! ierr  = INTEGER ::. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr RETURNs the minimum value
!         needed for nzmax. otherwise ierr=0 (normal RETURN).
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: nzmax
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(nzmax),  INTENT(INOUT) :: ao
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: iao
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: indu
  INTEGER,                          INTENT(OUT)   :: ierr
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(nzmax),  INTENT(INOUT) :: jao
! Counters:
  INTEGER                                         :: i
  INTEGER                                         :: ipos
  INTEGER                                         :: j
  INTEGER                                         :: k
  INTEGER                                         :: ko
  INTEGER                                         :: kfirst
  INTEGER                                         :: klast
  INTEGER                                         :: kosav
  INTEGER                                         :: lenrow
  INTEGER                                         :: nnz

  ierr = 0
  indu(1:nrow+1) = 0
!
!  Compute number of elements in each row of strict upper part.
!  Put result in INDU(I+1) for row I.
!
  DO i = 1, nrow
    DO k = ia(i), ia(i+1)-1
       j = ja(k)
       IF ( j < i ) THEN
         indu(j+1) = indu(j+1) + 1
       ENDIF
    ENDDO
  ENDDO
!
!  Find addresses of first elements of ouput matrix.  Result in INDU.
!
  indu(1) = 1
  DO i = 1, nrow
    lenrow = ia(i+1)-ia(i)
    indu(i+1) = indu(i) + indu(i+1) + lenrow
  ENDDO
!
!  Enough storage in A, JA?
!
  nnz = indu(nrow+1) - 1

  IF ( nzmax < nnz ) THEN
    ierr = nnz
    RETURN
  ENDIF
!
!  Now copy lower part (backwards).
!
  kosav = indu(nrow+1)

  DO i = nrow, 1, -1

     klast = ia(i+1) - 1
     kfirst = ia(i)
     iao(i+1) = kosav
     ko = indu(i)
     kosav = ko

     DO k = kfirst, klast
       ao(ko) = a(k)
       jao(ko) = ja(k)
       ko = ko+1
     ENDDO

     indu(i) = ko
  ENDDO

  iao(1) = 1
!
!  Copy upper part.  Go through the structure of AO, JAO, IAO
!  that has already been copied (lower part).  INDU(I) is the address
!  of the next free location in row I for AO, JAO.
!
!  I-th row is now in AO, JAO, IAO structure, lower half part.
!
  DO i = 1, nrow

    DO k = iao(i), iao(i+1)-1

      j = jao(k)

      IF ( i <= j ) THEN
        EXIT
      ENDIF

      ipos = indu(j)
      ao(ipos) = ao(k)
      jao(ipos) = i
      indu(j) = indu(j) + 1

    ENDDO

  ENDDO

END SUBROUTINE ssrcsr

SUBROUTINE submat ( n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao )

!*****************************************************************************80
!
!! SUBMAT extracts the submatrix A(i1:i2,j1:j2).
!
!  Discussion:
!
!    This routine extracts a submatrix and puts the result in
!    matrix ao,iao,jao.  It is an "in place" routine, so ao,jao,iao may be 
!    the same as a,ja,ia.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: N, the row dimension of the matrix.
!
! i1,i2 = two INTEGER ::s with i2 >= i1 indicating the range of rows to be
!          extracted.
!
! j1,j2 = two INTEGER ::s with j2 >= j1 indicating the range of columns
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! job      = job indicator: IF job /= 1 THEN the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!
! on output
!
! nr      = number of rows of submatrix
! nc      = number of columns of submatrix
!        * IF either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!      the column indices,and iao being the pointer to the beginning
!      of the row,in arrays a,ja.
!
  INTEGER,                     INTENT(IN)    :: n
  INTEGER,                     INTENT(IN)    :: job
  INTEGER,                     INTENT(IN)    :: i1
  INTEGER,                     INTENT(IN)    :: i2
  INTEGER,                     INTENT(IN)    :: j1
  INTEGER,                     INTENT(IN)    :: j2
  REAL(KIND=wp), DIMENSION(:), INTENT(IN)    :: a
  INTEGER,       DIMENSION(:), INTENT(IN)    :: ja
  INTEGER,       DIMENSION(:), INTENT(IN)    :: ia
! Output:
  INTEGER,                     INTENT(OUT)   :: nr
  INTEGER,                     INTENT(OUT)   :: nc
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: ao
  INTEGER,       DIMENSION(:), INTENT(INOUT) :: iao
  INTEGER,       DIMENSION(:), INTENT(INOUT) :: jao
! Counters:
  INTEGER                                    :: i
  INTEGER                                    :: ii
  INTEGER                                    :: j
  INTEGER                                    :: k
  INTEGER                                    :: k1
  INTEGER                                    :: k2
  INTEGER                                    :: klen

  nr = i2 - i1 + 1
  nc = j2 - j1 + 1

  IF ( nr <= 0 .or. nc <= 0 ) THEN
    RETURN
  ENDIF

  klen = 0
!
!  Simple procedure that proceeds row-wise.
!
  DO i = 1, nr

    ii = i1 + i - 1
    k1 = ia(ii)
    k2 = ia(ii+1) - 1
    iao(i) = klen + 1

    DO k = k1, k2
      j = ja(k)
      IF ( j1 <= j .and. j <= j2 ) THEN
        klen = klen + 1
        IF ( job == 1 ) THEN
          ao(klen) = a(k)
        ENDIF
        jao(klen) =j
      ENDIF
    ENDDO

  ENDDO

  iao(nr+1) = klen + 1

END SUBROUTINE submat

SUBROUTINE timestamp ( )

!*****************************************************************************80
!
! Purpose:
! =======
! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  ModIFied:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!

! Local variables:
  CHARACTER(LEN=40) :: string

  CALL timestring ( string )

  WRITE(*,'(a)' ) TRIM ( string )

END SUBROUTINE timestamp

SUBROUTINE timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  ModIFied:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  CHARACTER(LEN=*)      :: string
! Local variables:
  CHARACTER(LEN=8)      :: ampm
  INTEGER               :: d
  CHARACTER(LEN=8)      :: date
  INTEGER               :: h
  INTEGER               :: m
  INTEGER               :: mm
  INTEGER               :: n
  INTEGER               :: s
  CHARACTER(LEN=10)     :: time
  INTEGER, DIMENSION(8) :: values
  INTEGER               :: y
  CHARACTER(LEN=5)      :: zone
! Character constants:
  CHARACTER(LEN=9), PARAMETER, DIMENSION(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

  CALL DATE_AND_TIME( date, time, zone, values)

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  IF ( h < 12 ) THEN
    ampm = 'AM'
  else IF ( h == 12 ) THEN
    IF ( n == 0 .and. s == 0 ) THEN
      ampm = 'Noon'
    else
      ampm = 'PM'
    ENDIF
  else
    h = h - 12
    IF ( h < 12 ) THEN
      ampm = 'PM'
    else IF ( h == 12 ) THEN
      IF ( n == 0 .and. s == 0 ) THEN
        ampm = 'Midnight'
      else
        ampm = 'AM'
      ENDIF
    ENDIF
  ENDIF

  WRITE ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    TRIM ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, TRIM ( ampm )

  RETURN
END SUBROUTINE timestring

SUBROUTINE transp ( nrow, ncol, a, ja, ia, iwk, ierr )

!*****************************************************************************80
!
!  Purpose:
!  =======
!  TRANSP carries out in-place transposition routine.
!
!  Discussion:
!
!    This routine transposes a matrix stored in compressed sparse row
!    format.  The transposition is DOne in place in that the arrays 
!    A, JA, and IA of the transpose are overwritten onto the original arrays.
!
!    If you DO not need the transposition to be DOne in place
!    it is preferrable to use the conversion routine csrcsc
!    (see conversion routines in formats).
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order).  Use CSRCSC
!    IF you want them sorted.
!
!  ModIFied:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, INTEGER :: NROW, the row dimension of the matrix.
!
!    Input, INTEGER :: NCOL, the column dimension of the matrix.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Workspace, INTEGER :: IWK(*), of the same length as JA.
!
! on RETURN:
!
!
! ncol      = actual row dimension of the transpose of the input matrix.
!         Note that this may be <= the input value for ncol, in
!         case some of the last columns of the input matrix are zero
!         columns. In the case where the actual number of rows found
!         in transp(A) exceeds the input value of ncol, transp will
!         RETURN without completing the transposition. see ierr.
!
!    Input, real A(*), INTEGER :: JA(*), IA(NCOL+1), the transposed matrix
!    in CSR Compressed Sparse Row format.
!
! ierr      = INTEGER ::. error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!         THEN ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal RETURN).
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(INOUT) :: ncol
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: iwk
  INTEGER,                          INTENT(OUT)   :: ierr
! Local variables:
  INTEGER                                         :: i
  INTEGER                                         :: inext
  INTEGER                                         :: init
  INTEGER                                         :: j
  INTEGER                                         :: jcol
  INTEGER                                         :: k
  INTEGER                                         :: l
  INTEGER                                         :: nnz
  REAL(KIND=wp)                                   :: t
  REAL(KIND=wp)                                   :: t1

  ierr = 0
  nnz = ia(nrow+1) - 1
!
!  Determine the column dimension.
!
  jcol = 0
  DO k = 1, nnz
    jcol = max ( jcol, ja(k) )
  ENDDO

  IF ( ncol < jcol ) THEN
     ierr = jcol
     RETURN
  ENDIF
!
!  Convert to coordinate format.  Use IWK for row indices.
!
  ncol = jcol

  DO i = 1, nrow
    DO k = ia(i), ia(i+1)-1
      iwk(k) = i
    ENDDO
  ENDDO
!
!  Find pointer array for transpose.
!
  ia(1:ncol+1) = 0

  DO k = 1, nnz
    i = ja(k)
    ia(i+1) = ia(i+1) + 1
  ENDDO
  ia(1) = 1

  DO i = 1, ncol
    ia(i+1) = ia(i) + ia(i+1)
  ENDDO
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    CONTINUE

  t = a(init)
  i = ja(init)
  j = iwk(init)
  iwk(init) = -1

 6 CONTINUE

   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  l = ia(i)
!
!  Save the chased element.
!
  t1 = a(l)
  inext = ja(l)
!
!  Then occupy its location.
!
  a(l) = t
  ja(l) = j
!
!  Update pointer information for next element to be put in row i.
!
  ia(i) = l + 1
!
!  Determine next element to be chased.
!
  IF ( iwk(l) < 0 ) THEN
    GOTO 65
  ENDIF

  t = t1
  i = inext
  j = iwk(l)
  iwk(l) = -1

  IF ( k < nnz ) THEN
    GOTO 6
  ENDIF

  DO i = ncol, 1, -1
    ia(i+1) = ia(i)
  ENDDO

  ia(1) = 1

  RETURN

 65   CONTINUE

  init = init + 1

  IF ( nnz < init ) THEN

    DO i = ncol, 1, -1
      ia(i+1) = ia(i)
    ENDDO

    ia(1) = 1

    RETURN
  ENDIF

  IF ( iwk(init) < 0 ) THEN
    GOTO 65
  ENDIF
!
!  Restart chasing.
!
  GOTO 5
END SUBROUTINE transp

SUBROUTINE amudia ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! AMUDIA performs the matrix by matrix product B = A * Diag  (in place)
!
!  Discussion:
!
!    The column dimension of A is not needed.
!    The algorithm is "in place", so B can take the place of A.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer JOB, job indicator. Job=0 means get array b only
!    job = 1 means get b, and the integer arrays ib, jb.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), the diagonal matrix stored as a vector.
!
!    Output, B(*), JB(*), IB(NROW+1), the resulting matrix B in 
!    compressed sparse row sparse format.
!
  INTEGER,                          INTENT(IN)    :: nrow
  INTEGER,                          INTENT(IN)    :: job
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: a
  REAL(KIND=wp), DIMENSION(:),      INTENT(INOUT) :: b
  REAL(KIND=wp), DIMENSION(nrow),   INTENT(IN)    :: diag
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ia
  INTEGER,       DIMENSION(nrow+1), INTENT(INOUT) :: ib
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: ja
  INTEGER,       DIMENSION(:),      INTENT(INOUT) :: jb
! Counters:
  INTEGER ii
  INTEGER k
  INTEGER k1
  INTEGER k2

  DO ii = 1, nrow
!
!  Scale each element.
!
    k1 = ia(ii)
    k2 = ia(ii+1) - 1
    DO k = k1, k2
      b(k) = a(k) * diag(ja(k))
    ENDDO

  END DO

  IF ( job == 0 ) THEN
    RETURN
  END IF

  ib(1) = ia(1)
  DO ii = 1, nrow
    ib(ii) = ia(ii)
    DO k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    END DO
  END DO

END SUBROUTINE amudia

SUBROUTINE retmx ( n, a, ja, ia, dd )

!*****************************************************************************80
!
!  Purpose:
!  =======
!  RETMX returns in dd(*) the max absolute value of elements in row *.
!
!  Discussion:
!
!    This routine is used for scaling.  It has been superseded by RNRMS.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DD(N), the element of each row that has the largest absolute
!    value.  The sign of DD is modified such that it is the same as that 
!    of the diagonal element in its row.
!
  INTEGER,                       INTENT(IN)    :: n
  REAL(KIND=wp), DIMENSION(:),   INTENT(IN)    :: a
  REAL(KIND=wp), DIMENSION(:),   INTENT(INOUT) :: dd
  INTEGER,       DIMENSION(n+1), INTENT(IN)    :: ia
  INTEGER,       DIMENSION(:),   INTENT(IN)    :: ja
! Counters:
  INTEGER                                      :: i
  INTEGER                                      :: k
  INTEGER                                      :: k1
  INTEGER                                      :: k2
! Local variables:
  REAL(KIND=wp)                                :: t
  REAL(KIND=wp)                                :: t1
  REAL(kind=wp)                                :: t2
!
!  Initialize.
!
  k2 = 1

  DO i = 1, n

     k1 = k2
     k2 = ia(i+1) - 1
     t = 0.0_wp

     DO k = k1, k2

        t1 = ABS ( a(k) )
        IF ( t < t1 ) THEN
          t = t1
        END IF

        IF ( ja(k) == i ) THEN

          IF ( a(k) < 0.0D+00 ) THEN
            t2 = -1.0D+00
          ELSE IF ( a(k) == 0.0D+00 ) THEN
            t2 = 0.0D+00
          ELSE
            t2 = 1.0D+00
          ENDIF

        ENDIF

     ENDDO

     dd(i) = t2 * t
!
!  We do not invert the diagonal entries here.
!
  ENDDO

END SUBROUTINE retmx




END MODULE sparse
