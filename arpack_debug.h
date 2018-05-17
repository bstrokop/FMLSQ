!
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-05  Time: 17:44:43
 
!\SCCS Information: @(#)
! FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2

!     %---------------------------------%
!     | See debug.doc for documentation |
!     %---------------------------------%
INTEGER :: logfil, ndigit, mgetv0,  &
    msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,  &
    mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,  &
    mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
COMMON /debug/ logfil, ndigit, mgetv0,  &
    msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,  &
    mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,  &
    mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
