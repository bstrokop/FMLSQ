MODULE resolution_bins
USE constants
USE fail
USE mtz_io
USE util
IMPLICIT NONE
CONTAINS
    SUBROUTINE resbins ( mtz_1, nrange, maxranges )
!
!   Purpose:
!   =======
!   Setups resolution bins subject to a minimum number of reflections
!   and minimum shell thickness
!   
!   nrange - approximate number of refls in a lowest range of resolution
!
    TYPE(mtz),                                INTENT(INOUT) :: mtz_1
    INTEGER,                                  INTENT(IN)    :: nrange
    INTEGER,                                  INTENT(IN)    :: maxranges
!   Local vars:
    INTEGER                                                 :: nres
    INTEGER,       PARAMETER                                :: nhist = 1000
    INTEGER,       DIMENSION(:), ALLOCATABLE                :: hist
    REAL(KIND=wp)                                           :: ds
    REAL(KIND=wp)                                           :: sbin
    REAL(KIND=wp)                                           :: dsmin
!   Temporary array for very likely reduction of rres:
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE                :: temp
!   Counters:
    INTEGER                                                 :: i
    INTEGER                                                 :: ires
    INTEGER                                                 :: ibin
    INTEGER                                                 :: nbin    
    INTEGER                                                 :: total_refl

!   Set up histogram:
    ds = mtz_1%s_max_squared / ( nhist - 1 )
    dsmin = mtz_1%s_max_squared / maxranges

!   Initialize hist array:
    CALL allocate_array ( hist, nhist )
    hist = 0

!   Accumulate histogram:
    DO i = 1, SIZE ( mtz_1%hkl )
        ires = INT ( mtz_1%s(i) / ds ) + 1
        hist(ires) = hist(ires) + 1
    ENDDO

    IF ( .NOT. ALLOCATED ( mtz_1%rres ) ) CALL allocate_array ( mtz_1%rres, maxranges )

!   Now set resolution limits to accomodate minimum s and no. of refls:
    ibin = 0
    nbin = 0
    sbin = 0.0_wp
    total_refl = 0
    DO i = 1, nhist

!       FIXME: add condition whether we want this particular reflection

        nbin = nbin + hist(i)
        sbin = sbin + ds
        IF ( nbin > nrange .AND. sbin > dsmin ) THEN
            ibin = ibin + 1
            mtz_1%rres(ibin) = 1.0_wp / SQRT ( i * ds )
!            WRITE(*,"(' bin=', I4, ' resolution=', F6.2, ' refls num=', I6)") ibin, mtz_1%rres(ibin), nbin
            total_refl = total_refl + nbin
            nbin = 0
            sbin = 0.0_wp
        ENDIF
    ENDDO

!   Check if we could make at least one bin:
    IF ( ibin == 0 ) THEN
        CALL die ( ' Not enough reflections... Try other sub...', 'resbin' )
    ENDIF

!   It is tempting to add one more resolution bin (e.g. ibin=ibin+1), but it's more safe to merge it with
!   previous one:
    nres = ibin 

!   With almost 100% probabilty merges with previous bin:
    mtz_1%rres(nres) = 1.0_wp / SQRT ( mtz_1%s_max_squared )
    IF ( debug > 10000 ) THEN
        WRITE(*,"(' bin=', I4, ' resolution=', F6.2, ' refls num=', I6)") ibin, mtz_1%rres(ibin), nbin
    ENDIF

!   Free memory:
    CALL deallocate_array ( hist )

!   Reduce rres array:
    IF ( nres < maxranges ) THEN
        IF ( debug > 10000 ) THEN
            CALL messag ( 'Reducing RRES array...', 'resbin' )
        ENDIF
        CALL allocate_array ( temp, nres )
        temp = mtz_1%rres(1:nres)    
        CALL allocate_array  ( mtz_1%rres, nres )
        mtz_1%rres = temp
        CALL deallocate_array ( temp )
    ENDIF  

!   FINAL CHECK:
    total_refl = total_refl + nbin
    IF ( total_refl /= SIZE ( mtz_1%hkl ) ) THEN
        WRITE(*,*) total_refl, SIZE ( mtz_1%hkl )
        WRITE(*,*) SIZE ( mtz_1%hkl ) - total_refl, ' reflections disappeared'
        CALL die ( 'Programming error...', 'resbins' )
    ENDIF

    CALL allocate_array ( mtz_1%bin, SIZE ( mtz_1%s ) )
    DO i = 1, SIZE ( mtz_1%s )
        mtz_1%bin(i) = iresbin ( mtz_1%s(i), mtz_1%rres )
    ENDDO

    END SUBROUTINE resbins

    FUNCTION iresbin ( s, rres )
! 
!       Purpose:
!       =======
!       Quickly finds the resolution bin for this reflection
!       
!       Date:         Programmer:        Description of changes:
!       ====          ==========         ======================
!       Mar 2007      B.Strokopytov      Based on K.Cowtan's IRESBIN
!
    INTEGER                                              :: iresbin
    REAL(KIND=wp),                            INTENT(IN) :: s
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: rres
!   Counters:
    INTEGER                                              :: i

    DO i = 1, SIZE ( rres ) - 1
        IF ( s < 1.0_wp / rres(i) ** 2 ) EXIT
    ENDDO

!   Implicitly use the fact that if we fail to EXIT from the loop 
!   the counter 'i' will produce last resolution bin, i.e. i = i + 1
    iresbin = i

    END FUNCTION iresbin

END MODULE resolution_bins
