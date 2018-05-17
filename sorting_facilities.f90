MODULE sorting_facilities
USE select_kinds
USE fail
USE util
IMPLICIT NONE
INTERFACE sortag
    MODULE PROCEDURE isort
    MODULE PROCEDURE isort0
    MODULE PROCEDURE sort1
END INTERFACE sortag

CONTAINS
      SUBROUTINE isort ( a )
          INTEGER, DIMENSION(:), INTENT(INOUT) :: a
!         Local variables:
          INTEGER                              :: n
          n = SIZE(a)
          CALL isort0 (a, n)
      END SUBROUTINE isort

      SUBROUTINE isort0 ( a, n )
!
!     ISORT0 - sorts array A (integer*4)
!     N = size of A
! ******
      INTEGER,               INTENT(IN)    :: n
      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER                              :: t
      INTEGER                              :: tt
      INTEGER                              :: i
      INTEGER                              :: j
      INTEGER                              :: k
      INTEGER                              :: l
      INTEGER                              :: m
      INTEGER                              :: ij
      INTEGER, DIMENSION(24)               :: iu
      INTEGER, DIMENSION(24)               :: il
!===========================================================
      m = 1
      i = 1
      j = n
    5 IF ( i .GE. j ) GO TO 70
   10 k = i
      ij = (j + i)/2
      t = a(ij)
      IF ( a(i) .LE. t) GO TO 20
      a(ij) = a(i)
      a(i)  = t
      t     = a(ij)
   20 l = j
      IF ( a(j) .GE. t ) GO TO 40
      a(ij) = a(j)
      a(j)  = t
      t     = a(ij)
      IF ( a(i) .LE. t ) GO TO 40
      a(ij) = a(i)
      a(i)  = t
      t     = a(ij)
      GO TO 40
   30 a(l) = a(k)
      a(k) = tt
   40 l = l - 1
      IF ( a(l) .GT. t ) GO TO 40
      tt = a(l)
   50 k  = k + 1
      IF ( a(k) .LT. t ) GO TO 50
      IF ( k .LE. l) GO TO 30
      IF((l-i) .LE. (j-k)) GO TO 60
      il(m) = i
      iu(m) = l
      i = k
      m = m+1
      GO TO 80
   60 il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80
   70 m = m - 1
      IF( m .EQ. 0) RETURN
      i = il(m)
      j = iu(m)
   80 IF( (j-i) .GE. 1 ) GO TO 10
      IF( i .EQ. 1 ) GO TO 5
      i = i-1
   90 i = i+1
      IF ( i.EQ.j ) GO TO 70
      t = a(i+1)
      IF ( a(i) .LE. t ) GO TO 90
      k = i
  100 a(k+1) = a(k)
      k = k-1
      IF ( t .LT. a(k) ) GO TO 100
      a(k+1) = t
      GO TO 90
    END SUBROUTINE isort0

    SUBROUTINE sort1(A,N,TAG)
!
! -P- SORT1 - sort array A (real) with aray TAG
!                                (integer*4)
!     N = number of elements in A
      REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: a
      INTEGER,       DIMENSION(:), INTENT(INOUT) :: tag
      REAL(KIND=wp)                              :: t
      REAL(KIND=wp)                              :: tt
      INTEGER                                    :: tg
!     Local variables:
      INTEGER                                    :: n
      INTEGER, DIMENSION(24)                     :: iu
      INTEGER, DIMENSION(24)                     :: il
!     Counters:
      INTEGER                                    :: i
      INTEGER                                    :: j
      INTEGER                                    :: k
      INTEGER                                    :: l
      INTEGER                                    :: m
      INTEGER                                    :: ij
!===========================================================
      IF ( n < 1 ) THEN 
          CALL warn ( 'Nothing to sort.', 'sort1' )
          RETURN
      ENDIF
      m = 1
      i = 1
      j = n
    5 IF ( i .GE. j ) GO TO 70
   10 k = i
      ij = (j + i)/2
      t = a(ij)
      IF ( a(i) .LE. t ) GO TO 20
      a(ij) = a(i)
      a(i)  = t
      t     = a(ij)
      tg    = tag(ij)
      tag(ij) = tag(i)
      tag(i)  = tg
   20 l = j
      IF ( a(j) .GE. t ) GO TO 40
      a(ij) = a(j)
      a(j)  = t
      t     = a(ij)
      tg    = tag(ij)
      tag(ij) = tag(j)
      tag(j)  = tg
      IF ( a(i) .LE. t ) GO TO 40
      a(ij) = a(i)
      a(i)  = t
      t     = a(ij)
      tg    = tag(ij)
      tag(ij) = tag(i)
      tag(i)  = tg
      GO TO 40
   30 a(l) = a(k)
      a(k) = tt
      tg   = tag(l)
      tag(l) = tag(k)
      tag(k) = tg
   40 l = l - 1
      IF ( a(l) .GT. t ) GO TO 40
      tt = a(l)
   50 k  = k + 1
      IF ( a(k) .LT. t ) GO TO 50
      IF ( k .LE. l ) GO TO 30
      IF( (l-i) .LE. (j-k) ) GO TO 60
      il(m) = i
      iu(m) = l
      i     = k
      m     = m + 1
      GO TO 80
   60 il(m) = k
      iu(m) = j
      j     = l
      m     = m + 1
      GO TO 80
   70 m = m - 1
      IF ( m .EQ. 0 ) RETURN
      i = il(m)
      j = iu(m)
   80 IF( (j-i) .GE. 1 ) GO TO 10
      IF( i .EQ. 1 ) GO TO 5
      i = i - 1
   90 i = i + 1
      IF ( i .EQ. j ) GO TO 70
      t = a(i+1)
      IF ( a(i) .LE. t ) GO TO 90
      tg = tag(i+1)
      k  = i
  100 a(k+1)   = a(k)
      tag(k+1) = tag(k)
      k        = k-1
      IF ( t .LT. a(k) ) GO TO 100
      a(k+1)   = t
      tag(k+1) = tg
      GO TO 90
      END SUBROUTINE sort1

END MODULE sorting_facilities
