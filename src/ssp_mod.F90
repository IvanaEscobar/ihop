#include "IHOP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

!BOP
! !INTERFACE:
MODULE ssp_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>

  ! Holds SSP input by user and associated variables

  ! This module is very similar to the one used by the other programs in the
  ! Acoustics Toolbox. However, it returns the SSP *and* its derivatives

  ! Also, a greater premium has been placed on returning this info quickly,
  ! since IHOP calls it at every step so more information is pre-computed

  USE ihop_mod,     only: PRTFile

  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#ifdef ALLOW_USE_MPI
# include "EESUPPORT.h"
#endif
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public initSSP, setSSP, evalSSP, SSP, alphaR, betaR, &
           alphaI, betaI, rhoR, iSegz, iSegr

!=======================================================================

! INPUT/OUTPUT PARAMETERS:

! FUNCTIONS:
  _RL CHEN_MILLERO
  EXTERNAL CHEN_MILLERO

! LOCAL VARIABLES
! == Local Variables ==

! LEGACY VARIABLES
! == Legacy Local Variables ==
  INTEGER, PARAMETER     :: MaxSSP = 201
  INTEGER                :: iSegr = 1, iSegz = 1
  LOGICAL                :: foundr, foundz
#ifdef IHOP_THREED
  INTEGER                :: iSegx = 1, iSegy = 1
  LOGICAL                :: foundx, foundy
#endif /* IHOP_THREED */
  INTEGER                :: iostat, iallocstat
  INTEGER,           PRIVATE :: iz
  REAL (KIND=_RL90), PRIVATE :: Depth, W
  ! DEFAULT values, IHOP only modifies alphaR
  REAL (KIND=_RL90)      :: alphaR = 1500, betaR = 0, alphaI = 0, &
                            betaI = 0, rhoR = 1
  ! SSP interpolation parameters, only used in ssp_mod
  COMPLEX (KIND=_RL90) :: n2(MaxSSP), n2z(MaxSSP)
  COMPLEX (KIND=_RL90) :: cSpln( 4, MaxSSP ), cCoef( 4, MaxSSP )

! TYPE STRUCTURES
! == Type Structures ==
  TYPE rxyz_vector
    REAL (KIND=_RL90), ALLOCATABLE :: r(:)
#ifdef IHOP_THREED
    REAL (KIND=_RL90), ALLOCATABLE :: x(:), y(:), z(:)
#endif
  END TYPE rxyz_vector

  ! SSP
  TYPE SSPStructure
    INTEGER                 :: nPts, nR, nX, nY, nZ
    REAL    (KIND=_RL90)    :: z( MaxSSP ), rho( MaxSSP )
    COMPLEX (KIND=_RL90)    :: c( MaxSSP ), cz( MaxSSP )
    REAL    (KIND=_RL90), ALLOCATABLE   :: cMat( :, : ), czMat( :, : )
#ifdef IHOP_THREED
    REAL    (KIND=_RL90), ALLOCATABLE   :: cMat3( :, :, : ), czMat3( :, :, : )
#endif /* IHOP_THREED */
    TYPE ( rxyz_vector ) :: Seg
    CHARACTER (LEN=1)    :: Type
    CHARACTER (LEN=2)    :: AttenUnit
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP
!EOP

CONTAINS
!**********************************************************************!
  SUBROUTINE initSSP( myThid )

  ! Call the particular profile routine indicated by the SSP%Type and
  ! perform initialize SSP structures

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid
    ! Check if SSPFile exists
    IF (useSSPFile) THEN
      CALL ReadSSP( myThid )
    ELSE
      CALL init_fixed_SSP( myThid )
    END IF

  RETURN
  END !SUBROUTINE initSSP

!**********************************************************************!
  SUBROUTINE setSSP( myThid )
    USE pchip_mod,  only: PCHIP
    USE splinec_mod,only: cspline

    ! Set SSP structures for particular SSP%Type routine

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    INTEGER :: ir, iz
#ifdef ALLOW_USE_MPI
    INTEGER :: mpiRC
    mpiRC  = 0
#endif

!$TAF init setssp1 = 'ssp_mod_setssp'

! IESCO24: Write derived type with allocatable memory by type: SSP from ssp_mod
! Scalar components
! Fixed arrays
!$TAF store ssp%c = setssp1
! Allocatable arrays
!$TAF store ssp%cmat,ssp%czmat,ssp%seg%r = setssp1

    ! init defaults for ssp_mod scoped arrays
    n2    = (-1.,-1.)
    n2z   = (-1.,-1.)
    cSpln = (-1.,-1.)
    cCoef = (-1.,-1.)

    ! Extract gcm SSP field
    IF ( .NOT.useSSPFile ) CALL gcmSSP( myThid )

    ! Write to PRTFile
    IF ( .NOT.usingMPI ) THEN
      CALL writeSSP( myThid )
#ifdef ALLOW_USE_MPI
    ELSE ! using MPI
      CALL MPI_COMM_RANK( MPI_COMM_MODEL, mpiMyId, mpiRC )
      myProcId = mpiMyId

      ! Hard coded write on single proc
      IF ( myProcId.EQ.0 ) CALL writeSSP( myThid )
#endif /* ALOOW_USE_MPI */
    ENDIF ! IF ( .NOT.usingMPI )

    ! Populate rest of SSP derived type based on SSP interpolation scheme
    SELECT CASE ( SSP%Type )
    CASE ( 'N' )  !  N2-linear profile option
      n2( 1:SSP%nPts ) = 1.0 / SSP%c( 1:SSP%nPts )**2
      !IEsco23 Test this: n2(  1:SSP%nZ ) = 1.0 / SSP%c( 1:SSP%nZ )**2

      ! compute gradient, n2z
      DO iz = 2, SSP%nPts
        n2z( iz-1 ) = (  n2(   iz ) - n2(    iz-1 ) ) / &
                      ( SSP%Z( iz ) - SSP%Z( iz-1 ) )
      END DO

    CASE ( 'C' )  !  C-linear profile option
    CASE ( 'P' )  !  monotone PCHIP ACS profile option
      !                                                               2      3
      ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
      !
      CALL PCHIP( SSP%Z, SSP%c, SSP%nPts, cCoef, cSpln )
!IEsco23 Test this:
!     CALL PCHIP( SSP%Z, SSP%c, SSP%nZ, cCoef, cSpln )

    CASE ( 'S' )  !  Cubic spline profile option
      cSpln( 1, 1:SSP%nPts ) = SSP%c( 1:SSP%nPts )
!IEsco23 Test this:
!     cSpln( 1, 1 : SSP%nZ ) = SSP%c( 1 : SSP%nZ )

      ! Compute spline coefs
      CALL cSpline( SSP%Z, cSpln( 1, 1 ), SSP%nPts, 0, 0, SSP%nPts )
!IEsco23 Test this:
!     CALL CSpline( SSP%Z, cSpln( 1,1 ), SSP%nZ,iBCBeg, iBCEnd, SSP%nZ )

    CASE ( 'Q' )
      ! calculate cz
      DO ir = 1, SSP%nR
        DO iz = 2, SSP%nZ
          ! delta_z = ( SSP%Z( iz2 ) - SSP%Z( iz2-1 ) )
          SSP%czMat( iz-1, ir ) = ( SSP%cMat( iz,   ir ) - &
                                    SSP%cMat( iz-1, ir ) ) / &
                                  ( SSP%Z( iz ) - SSP%Z( iz-1 ) )
        END DO
      END DO

    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(A)') 'SSPMOD setSSP: Invalid SSP profile option'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R setSSP'
    END SELECT

  RETURN
  END !SUBROUTINE setSSP

!**********************************************************************!
  SUBROUTINE evalSSP( x, c, cimag, gradc, crr, crz, czz, rho, myThid )

    ! Call the particular profile routine indicated by the SSP%Type and
    ! tabulate cp, cs, rhoT

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x(2)  ! r-z SSP evaluation point
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho

    SELECT CASE ( SSP%Type )
    CASE ( 'N' )  !  N2-linear profile option
       CALL n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    CASE ( 'C' )  !  C-linear profile option
       CALL cLinear(  x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    CASE ( 'P' )  !  monotone PCHIP ACS profile option
       CALL cPCHIP(   x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    CASE ( 'S' )  !  Cubic spline profile option
       CALL cCubic(   x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    CASE ( 'Q' )
       CALL Quad(     x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'Profile option: ', SSP%Type
        ! In adjoint mode we do not write output besides on the first run
        IF (IHOP_dumpfreq.GE.0) &
            CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 'SSPMOD evalSSP: Invalid SSP profile option'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R evalSSP'
        c = 0.         !RG
        cimag = 0.     !RG
        gradc = 0.     !RG
        crr   = 0.     !RG
        crz   = 0.     !RG
        czz   = 0.     !RG
        rho   = 0.     !RG
    END SELECT

  RETURN
  END !SUBROUTINE evalSSP

!**********************************************************************!
  SUBROUTINE n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, myThid )

    ! N2-linear interpolation of SSP data
    ! Return SSP info

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, &
                                        rho ! sound speed and its derivatives


    iSegz = 1                   !RG
    IF ( x( 2 ) < SSP%Z( iSegz ) .OR. x( 2 ) > SSP%Z( iSegz+1 ) ) THEN
      foundz=.false.
!IEsco23 Test this:
!     DO iz = 2, SSP%nZ   ! Search for bracketting Depths
      DO iz = 2, SSP%nPts   ! Search for bracketting Depths
        IF ( x( 2 ) < SSP%Z( iz ) .and. .not. foundz ) THEN
          iSegz  = iz - 1
          foundz = .true.
        END IF
      END DO
    END IF

    W = ( x( 2 ) - SSP%Z( iSegz ) ) / ( SSP%Z( iSegz+1 ) - SSP%Z( iSegz ) )

    c     = REAL(  1.0D0 / SQRT( ( 1.0D0-W ) * n2( iSegz ) &
            + W * n2( iSegz+1 ) ) )
    cimag = AIMAG( 1.0D0 / SQRT( ( 1.0D0-W ) * n2( iSegz ) &
            + W * n2( iSegz+1 ) ) )

    gradc = [ 0.0D0, -0.5D0 * c * c * c * REAL( n2z( iSegz ) ) ]
    crr   = 0.0d0
    crz   = 0.0d0
    czz   = 3.0d0 * gradc( 2 ) * gradc( 2 ) / c

    rho   = ( 1.0D0-W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz+1 )
  RETURN
  END !SUBROUTINE n2Linear

!**********************************************************************!
  SUBROUTINE cLinear( x, c, cimag, gradc, crr, crz, czz, rho, myThid  )

  ! c-linear interpolation of SSP data

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x(2)  ! r-z SSP evaluation point
    ! sound speed and its derivatives
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, rho

    iSegz = 1                   !RG
    IF ( x(2) < SSP%Z( iSegz ) .OR. x(2) > SSP%Z( iSegz+1 ) ) THEN
       foundz=.false.
!IEsco23 Test this:
!      DO iz = 2, SSP%nZ   ! Search for bracketting Depths
       DO iz = 2, SSP%nPts   ! Search for bracketting Depths
          IF ( x(2) < SSP%Z( iz ) .and. .not. foundz ) THEN
             iSegz  = iz - 1
             foundz = .true.
          END IF
       END DO
    END IF

    c     = REAL(  SSP%c( iSegz ) + ( x(2) - SSP%Z( iSegz ) ) * SSP%cz( iSegz ) )
    cimag = AIMAG( SSP%c( iSegz ) + ( x(2) - SSP%Z( iSegz ) ) * SSP%cz( iSegz ) )
    gradc = [ 0.0D0, REAL( SSP%cz( iSegz ) ) ]
    crr   = 0.0d0
    crz   = 0.0d0
    czz   = 0.0d0

    W     = ( x(2) - SSP%Z( iSegz ) ) / ( SSP%Z( iSegz+1 ) - SSP%Z( iSegz ) )
    rho   = ( 1.0D0-W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cLinear

!**********************************************************************!
  SUBROUTINE cPCHIP( x, c, cimag, gradc, crr, crz, czz, rho, myThid )

    ! This implements the monotone piecewise cubic Hermite interpolating
    ! polynomial (PCHIP) algorithm for the interpolation of the sound speed c.


  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x(2)  ! r-z SSP evaluation point
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    REAL    (KIND=_RL90) :: xt
    COMPLEX (KIND=_RL90) :: c_cmplx


    iSegz = 1                   !RG
    IF ( x(2) < SSP%Z( iSegz ) .OR. x(2) > SSP%Z( iSegz+1 ) ) THEN
       foundz=.false.
!IEsco23 Test this:
!      DO iz = 2, SSP%nZ   ! Search for bracketting Depths
       DO iz = 2, SSP%nPts   ! Search for bracketting Depths
          IF ( x(2) < SSP%Z( iz ) .and. .not. foundz ) THEN
             iSegz  = iz - 1
             foundz = .true.
          END IF
       END DO
    END IF

    xt = x( 2 ) - SSP%Z( iSegz )
    c_cmplx = cCoef( 1, iSegz ) &
          + ( cCoef( 2, iSegz ) &
          + ( cCoef( 3, iSegz ) &
          +   cCoef( 4, iSegz ) * xt ) * xt ) * xt

    c     = REAL(  c_cmplx )
    cimag = AIMAG( c_cmplx )

    gradc = [ 0.0D0, &
              REAL( cCoef( 2, iSegz ) + ( 2.0D0 * cCoef( 3, iSegz ) &
                    + 3.0D0 * cCoef( 4, iSegz ) * xt ) * xt ) ]

    crr   = 0.0D0
    crz   = 0.0D0
    czz   = REAL( 2.0D0 * cCoef( 3, iSegz ) + &
                  6.0D0 * cCoef( 4, iSegz ) * xt )   ! dgradc(2)/dxt

    W     = ( x(2) - SSP%Z( iSegz ) ) / ( SSP%Z( iSegz+1 ) - SSP%Z( iSegz ) )
    ! linear interp of density
    rho   = ( 1.0D0-W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cPCHIP

!**********************************************************************!
  SUBROUTINE cCubic( x, c, cimag, gradc, crr, crz, czz, rho, myThid  )

    ! Cubic spline interpolation

    USE splinec_mod,  only: splineall

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid

  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x(2)  ! r-z SSP evaluation point
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    REAL     (KIND=_RL90)   :: hSpline
    COMPLEX  (KIND=_RL90)   :: c_cmplx, cz_cmplx, czz_cmplx

    ! *** Section to return SSP info ***

    iSegz = 1                   !RG
     IF ( x(2) < SSP%Z( iSegz ) .OR. x(2) > SSP%Z( iSegz+1 ) ) THEN
        foundz=.false.
!IEsco23 Test this:
!       DO iz = 2, SSP%nZ   ! Search for bracketting Depths
        DO iz = 2, SSP%nPts   ! Search for bracketting Depths
           IF ( x(2) < SSP%Z( iz ) .and. .not. foundz ) THEN
              iSegz  = iz - 1
              foundz = .true.
           END IF
        END DO
     END IF

    hSpline = x(2) - SSP%Z( iSegz )

    CALL SplineALL( cSpln( 1, iSegz ), hSpline, c_cmplx, cz_cmplx, czz_cmplx )

    c     = DBLE(  c_cmplx )
    cimag = AIMAG( c_cmplx )
    gradc = [ 0.0D0, DBLE( cz_cmplx ) ]
    czz   = DBLE( czz_cmplx )
    crr   = 0.0d0
    crz   = 0.0d0

    ! linear interpolation for density
    W   = ( x(2) - SSP%Z( iSegz ) ) / ( SSP%Z( iSegz+1 ) - SSP%Z( iSegz ) )
    rho = ( 1.0D0-W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz+1 )

  RETURN
  END !SUBROUTINE cCubic

!**********************************************************************!
  SUBROUTINE Quad( x, c, cimag, gradc, crr, crz, czz, rho, myThid )
    ! Bilinear quadrilatteral interpolation of SSP data in 2D, SSP%Type = 'Q'

    !     == Routine Arguments ==
    !     myThid :: Thread number. Unused by IESCO
    !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

    !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x(2)  ! r-z SSP evaluation point
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc(2), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    INTEGER             :: irT, iz2
    INTEGER             :: isegzold
    REAL (KIND=_RL90)   :: c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z

    ! *** Section to return SSP info ***

    ! IESCO22: iSegz is the depth index containing x depth
    ! find depth-layer where x(2) in ( SSP%Z( iSegz ), SSP%Z( iSegz+1 ) )
    iSegz = 1                   !RG
    IF ( x(2) < SSP%Z( iSegz ) .OR. x(2) > SSP%Z( iSegz+1 ) ) THEN
       foundz=.false.
       DO iz = 2, SSP%nZ   ! Search for bracketting Depths
          IF ( x(2) < SSP%Z( iz ) .and. .not. foundz ) THEN
             iSegz  = iz - 1
             foundz = .true.
          END IF
       END DO
    END IF

    ! Check that x is inside the box where the sound speed is defined
    IF ( x(1) < SSP%Seg%R( 1 ) .OR. x(1) > SSP%Seg%R( SSP%nR ) ) THEN
#ifdef IHOP_WRITE_OUT
     ! In adjoint mode we do not write output besides on the first run
     IF (IHOP_dumpfreq.GE.0) THEN
      WRITE(msgBuf,'(2A)') 'ray is outside the box where ocean ',&
        'soundspeed is defined'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,2F13.4)') ' x = ( r, z ) = ', x
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
        'ray is outside the box where the soundspeed is defined'
      CALL PRINT_ERROR( msgBuf,myThid )
     ENDIF
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R Quad'
    END IF

    ! find range-segment where x(1) in [ SSP%Seg%R( iSegr ), SSP%Seg%R( iSegr+1 ) )
    iSegr = 1           !RG
    IF ( x(1) < SSP%Seg%R( iSegr ) .OR. x(1) >= SSP%Seg%R( iSegr+1 ) ) THEN
      foundr=.false.
      DO irT = 2, SSP%nR   ! Search for bracketting segment ranges
        IF ( x(1) < SSP%Seg%R( irT ) .and. .not. foundr ) THEN
          iSegr = irT - 1
          foundr=.true.
        END IF
      END DO
    END IF

    ! for depth, x(2), get the sound speed at both ends of range segment
    cz1 = SSP%czMat( iSegz, iSegr   )
    cz2 = SSP%czMat( iSegz, iSegr+1 )

    !IESCO22: s2 is distance btwn field point, x(2), and ssp depth @ iSegz
    s2      = x(2)             - SSP%Z( iSegz )
    delta_z = SSP%Z( iSegz+1 ) - SSP%Z( iSegz )
    IF (delta_z <= 0 .OR. s2 > delta_z) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf, *) delta_z, s2, iSegz, SSP%Z(iSegz)
      CALL PRINT_ERROR( msgBuf,myThid )
      WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
        'depth is not monotonically increasing in SSP%Z'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R Quad'
    END IF

    c1 = SSP%cMat( iSegz, iSegr   ) + s2*cz1
    c2 = SSP%cMat( iSegz, iSegr+1 ) + s2*cz2

    ! s1 = proportional distance of x(1) in range
    delta_r = SSP%Seg%R( iSegr+1 ) - SSP%Seg%R( iSegr )
    s1 = ( x(1) - SSP%Seg%R( iSegr ) ) / delta_r
    ! piecewise constant extrapolation for ranges outside SSP box
    s1 = MIN( s1, 1.0D0 )
    s1 = MAX( s1, 0.0D0 )

    c = ( 1.0D0-s1 )*c1 + s1*c2 ! c @ x

    ! interpolate the attenuation !!!! SSP in ENVFile needs to match first column of SSPFile
    s2    = s2 / delta_z   ! normalize depth layer
    ! volume attenuation is taken from the single c(z) profile
    cimag = AIMAG( ( 1.0D0-s2 )*SSP%c( iSegz ) + s2*SSP%c( iSegz+1 ) )

    cz  = ( 1.0D0-s1 )*cz1 + s1*cz2 ! cz @ x

    cr  = ( c2  - c1  ) / delta_r ! SSP grid cr
    crz = ( cz2 - cz1 ) / delta_r ! SSP grid crz

    gradc = [ cr, cz ]
    crr   = 0.0
    czz   = 0.0

    ! linear interpolation for density
    W   = ( x(2) - SSP%Z( iSegz ) ) / ( SSP%Z( iSegz+1 ) - SSP%Z( iSegz ) )
    rho = ( 1.0D0-W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz+1 )

    !IESCO22: for thesis, czz=crr=0, and rho=1 at all times
  RETURN
  END !SUBROUTINE Quad

!**********************************************************************!
  SUBROUTINE ReadSSP( myThid )
    USE atten_mod, only: CRCI
    USE ihop_mod,  only: SSPFile
    USE bdry_mod, only: Bdry
    ! reads SSP in m/s from .ssp file and convert to AttenUnit (ie. Nepers/m)
    ! Populates SSPStructure: SSP

  ! IESCO24
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  !USE initenvihop, only: fT

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
    REAL (KIND=_RL90) :: Depth
    REAL (KIND=_RL90) :: bPower, fT

    ! IESCO24 fT init
    bPower = 1.0
    fT = 1000.0
    Depth = Bdry%Bot%HS%Depth

    ! I/O on main thread only
    _BEGIN_MASTER(myThid)
    ! OPEN SSPFile to read
    OPEN ( FILE = TRIM(IHOP_fileroot) // '.ssp', UNIT = SSPFile, &
        FORM = 'FORMATTED', STATUS = 'OLD', IOSTAT = iostat )
    IF ( IOSTAT /= 0 ) THEN   ! successful open?
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'SSPMOD ReadSSP: Unable to open the SSP file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSSP'
    END IF

    READ( SSPFile,  * ) SSP%nR, SSP%nZ

    ALLOCATE( SSP%cMat( SSP%nZ, SSP%nR ), &
              SSP%czMat( SSP%nZ-1, SSP%nR ), &
              SSP%Seg%R( SSP%nR ), &
              STAT = iallocstat )
    IF ( iallocstat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
                            'Insufficient memory to store SSP'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSSP'
    END IF

    ! Initiate to nonsense
    SSP%cMat  = -99.0 _d 0
    SSP%czMat = -99.0 _d 0

    ! Set SSP contents
    READ( SSPFile,  * ) SSP%Seg%R( 1:SSP%nR )
    SSP%Seg%R = 1000.0 * SSP%Seg%R   ! convert km to m

    READ( SSPFile,  * ) SSP%Z( 1:SSP%nZ )

    DO iz = 1, SSP%nZ
       READ(  SSPFile, * ) SSP%cMat( iz, : )
    END DO
    CLOSE( SSPFile )

    SSP%nPts = 1
    DO iz = 1, MaxSSP
       alphaR = SSP%cMat( iz, 1 )

       SSP%c(iz) = CRCI( SSP%Z(iz), alphaR, alphaI, SSP%AttenUnit, bPower, fT, &
                           myThid )
       SSP%rho(iz) = rhoR !IEsco22: set to a default value of 1

       ! verify depths are monotone increasing
       IF ( iz > 1 ) THEN
          IF ( SSP%Z( iz ) .LE. SSP%Z( iz-1 ) ) THEN
#ifdef IHOP_WRITE_OUT
            WRITE(msgBuf,'(2A,F10.2)') 'SSPMOD ReadSSP: ', &
                  'The depths in the SSP must be monotone increasing', SSP%Z(iz)
            CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
            STOP 'ABNORMAL END: S/R ReadSSP'
          END IF
       END IF

       ! compute gradient, cz
       IF ( iz > 1 ) SSP%cz( iz-1 )  = ( SSP%c( iz ) - SSP%c( iz-1 ) ) / &
                                       ( SSP%Z( iz ) - SSP%Z( iz-1 ) )

       ! Did we read the last point?
       IF ( ABS( SSP%Z( iz ) - Depth ) < 100. * EPSILON( 1.0e0 ) ) THEN
          IF ( SSP%nPts == 1 ) THEN
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
                                     'The SSP must have at least 2 points'
                CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
                STOP 'ABNORMAL END: S/R ReadSSP'
          END IF

          RETURN
       ENDIF

       SSP%nPts = SSP%nPts + 1
    END DO

    ! Fall through means too many points in the profile
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
                         'Number of SSP points exceeds limit'
    CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
    STOP 'ABNORMAL END: S/R ReadSSP'

    ! I/O on main thread only
    _END_MASTER(myThid)
  RETURN
  END !SUBROUTINE ReadSSP

!**********************************************************************!
SUBROUTINE init_fixed_ssp( myThid )
  ! Initiate parameters that don't change within a time series
  ! Sets SSP%nR,nZ,Seg%R, and ihop_sumweights
  USE bdry_mod, only: Bdry

  ! == Routine Arguments ==
  ! myThid :: Thread number. Unused by IESCO
  ! msgBuf :: Used to build messages for printing.
  INTEGER, INTENT(IN)     :: myThid
  CHARACTER*(MAX_LEN_MBUF):: msgBuf

  ! == Local Variables ==
  LOGICAL :: skip_range
  INTEGER :: iallocstat
  INTEGER :: bi,bj, i,j,k, ii,jj
  INTEGER :: nii(IHOP_NPTS_RANGE), njj(IHOP_NPTS_RANGE)
  REAL (KIND=_RL90) :: tolerance

  ! init local vars
  skip_range =.false.
  nii(:)    = 0
  njj(:)    = 0
  tolerance  = 5 _d -5

  ! init default SSP values (only fixed memory vars)
  SSP%nPts = -1
  SSP%Z    = -999.0
  SSP%rho  = -999.0
  SSP%c    = (-999.0, 0.0)
  SSP%cz   = (-999.0, 0.0)

  ! set ihop SSP grid size
  SSP%nZ = Nr+2 ! add z=0 z=Depth layers to GCM Nr
  SSP%nR = IHOP_NPTS_RANGE
  SSP%nPts = SSP%nZ

  ! set SSP%Z from rC, rkSign=-1 used bc ihop uses +ive depths
  SSP%Z( 1 )            = 0.0 _d 0
  SSP%Z( 2:(SSP%nZ-1) ) = rkSign*rC( 1:Nr )
  SSP%Z( SSP%nZ )       = Bdry%Bot%HS%Depth ! rkSign*rF(Nr+1)*1.05

  ! set SSP%Seg%R from data.ihop -> ihop_ranges
  !IF (ALLOCATED(SSP%Seg%R)) DEALLOCATE(SSP%Seg%R)
  ALLOCATE( SSP%Seg%R( SSP%nR ), STAT = iallocstat )
  IF ( iallocstat /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD init_fixed_SSP: ', &
      'Insufficient memory to store SSP%Seg%R'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R init_fixed_SSP'
  END IF

  SSP%Seg%R( 1:SSP%nR ) = ihop_ranges( 1:SSP%nR )
  ! Modify from [m] to [km]
  SSP%Seg%R = 1000.0 * SSP%Seg%R

  ! ONLY ALLOCATE cmat and czmat, to be filled per ihop run
  ALLOCATE( SSP%cMat( SSP%nZ, SSP%nR ), &
            SSP%czMat( SSP%nZ-1, SSP%nR ), &
            STAT = iallocstat )
  IF ( iallocstat /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD init_fixed_SSP: ', &
      'Insufficient memory to store SSP%cmat, SSP%czmat'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R init_fixed_SSP'
  END IF
  ! Initiate to nonsense
  SSP%cMat  = -99.0 _d 0
  SSP%czMat = -99.0 _d 0

  !==================================================
  ! IDW Interpolate: COMPARING with LAT LONs (xC, yC)
  !==================================================
  ! Sum IDW weights
  DO ii = 1, SSP%nR
    ihop_sumweights(ii,:) = sum(ihop_idw_weights(ii,:))
  END DO

  ! Adapt IDW interpolation by bathymetry
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      DO j=1,sNy
        DO i=1,sNx
          DO ii=1, SSP%nR
            skip_range = .FALSE.

            DO jj=1,IHOP_npts_idw
              IF (ABS(xC(i, j, bi, bj) - ihop_xc(ii, jj)) .LE. tolerance .AND. &
                  ABS(yC(i, j, bi, bj) - ihop_yc(ii, jj)) .LE. tolerance) THEN
                DO k=1,Nr
    ! No IDW interpolation on xc, yc centered ranges
    IF (nii(ii) .EQ. 1 .AND. k .GT. njj(ii)) THEN
      skip_range = .TRUE.
    END IF

    IF (.NOT. skip_range) THEN
      IF (hFacC(i, j, k, bi, bj) .EQ. 0.0) THEN
        ihop_sumweights(ii, k) = &
            ihop_sumweights(ii, k) - ihop_idw_weights(ii, jj)

        ! No interpolation on xc, yc centered ranges
        IF (ihop_idw_weights(ii, jj) .EQ. 0.0) THEN
          ihop_sumweights(ii, k:Nr) = 0.0
          nii(ii) = 1
          njj(ii) = k
        END IF
      END IF

      ! Set TINY and negative values to 0.0
      IF (ihop_sumweights(ii, k) .LT. 1D-13) ihop_sumweights(ii, k) = 0.0
    END IF
                END DO !k
              END IF
            END DO !jj
          END DO !ii
        END DO !i
      END DO !j
    END DO !bi
  END DO !bj

END !SUBROUTINE init_fixed_ssp

!**********************************************************************!
SUBROUTINE gcmSSP( myThid )
  ! Extracts SSP from MITgcm grid points

  USE atten_mod, only: CRCI
  USE bdry_mod, only: Bdry

#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

  ! == Routine Arguments ==
  ! myThid :: Thread number. Unused by IESCO
  ! msgBuf :: Used to build messages for printing.
  INTEGER, INTENT(IN)     :: myThid
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER*(80)          :: fmtstr

  ! == Local Variables ==
  ! interp_finished : is adaptive IDW interp finished at this range point?
  ! njj             : number interpolation points per depth level
  LOGICAL :: interp_finished
  INTEGER :: iallocstat
  INTEGER :: bi,bj, i,j,k, ii,jj
  INTEGER :: njj(IHOP_NPTS_RANGE)
  REAL (KIND=_RL90)             :: dcdz, tolerance
  REAL (KIND=_RL90), ALLOCATABLE:: tileSSP(:,:,:,:), tmpSSP(:,:,:), globSSP(:)
  ! IESCO24
  ! fT = 1000 ONLY for acousto-elastic halfspaces, I will have to pass this
  ! parameter in a different way after ssp_mod is split btwn fixed and varia
  REAL (KIND=_RL90)             :: bPower, fT
#ifdef ALLOW_AUTODIFF_TAMC
  INTEGER tkey, ijkey, hkey, lockey

!$TAF init loctape_ihop_gcmssp_bibj_ij_iijj_k = STATIC, nSx*nSy*sNx*sNy*IHOP_MAX_RANGE*IHOP_MAX_NC_SIZE*(Nr + 2)
#endif

  ! IESCO24 fT init
  bPower = 1.0
  fT = 1000.0

  ! init local vars
  interp_finished =.false.
  njj       = 0
  dcdz      = 0.0 _d 0
  tolerance = 5 _d -5

  IF(ALLOCATED(tileSSP)) DEALLOCATE(tileSSP)
  IF(ALLOCATED(globSSP)) DEALLOCATE(globSSP)
!  ALLOCATE( tileSSP(SSP%nZ,SSP%nR,nSx,nSy), STAT = iallocstat )
  ALLOCATE( tileSSP(SSP%nZ,SSP%nR,nSx,nSy), tmpSSP(nSx,nSy,SSP%nZ*SSP%nR), &
            globSSP(SSP%nZ*SSP%nR), STAT = iallocstat )
  IF ( iallocstat /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD gcmSSP: ', &
      'Insufficient memory to store tileSSP and/or globSSP'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R gcmSSP'
  END IF
  ! Initiate to ceros
  tileSSP    = 0.0 _d 0
  tmpSSP     = 0.0 _d 0
  globSSP    = 0.0 _d 0

  ! interpolate SSP with adaptive IDW from gcm grid to ihop grid
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
#ifdef ALLOW_AUTODIFF_TAMC
      tkey = bi + (bj-1)*nSx + (ikey_dynamics-1)*nSx*nSy
#endif

! IESCO24: don't worry about gcm overlaps right now
      DO j=1,sNy
        DO i=1,sNx
#ifdef ALLOW_AUTODIFF_TAMC
          ijkey = i + (j-1)*sNx + (tkey-1)*sNx*sNy
!$TAF store njj(ii) = comlev1_bibj_ij_ihop, key=ijkey, kind=isbyte
#endif
          DO ii=1,IHOP_npts_range
            interp_finished = .FALSE.
            DO jj=1,IHOP_npts_idw
#ifdef ALLOW_AUTODIFF_TAMC
!$TAF STORE interp_finished = comlev1_bibj_ij_ihop, key=ijkey, kind=isbyte
#endif

              ! Interpolate from GCM grid cell centers
              IF (ABS(xC(i, j, bi, bj) - ihop_xc(ii, jj)) .LE. tolerance .AND. &
                  ABS(yC(i, j, bi, bj) - ihop_yc(ii, jj)) .LE. tolerance .AND. &
                  .NOT. interp_finished) THEN
                njj(ii) = njj(ii) + 1

#ifdef ALLOW_AUTODIFF_TAMC
#endif
                DO iz = 1, SSP%nZ - 1
#ifdef ALLOW_AUTODIFF_TAMC
                  hkey = jj + (ii-1)*IHOP_npts_idw + (ijkey-1)*sNy*sNx*nSy*nSx
                  lockey = iz + (hkey-1)*(SSP%nZ-1)*IHOP_npts_idw*IHOP_npts_range*sNx*sNy*nSx*nSy
!                  + ((jj-1) + ((ii-1) + (ijkey-1)*IHOP_npts_range)*IHOP_npts_idw)*(SSP%nZ-1)
!$TAF store njj(ii) = loctape_ihop_gcmssp_bibj_ij_iijj_k, key=lockey, kind=isbyte
#endif

    IF (iz .EQ. 1) THEN
      ! Top vlevel zero depth
      tileSSP(1, ii, bi, bj) = tileSSP(1, ii, bi, bj) + &
        CHEN_MILLERO(i, j, 0, bi, bj, myThid) * &
        ihop_idw_weights(ii, jj) / ihop_sumweights(ii, iz)
    ELSE ! 2:(SSP%nZ-1)
      ! Middle depth layers, only when not already underground
      IF (ihop_sumweights(ii, iz-1) .GT. 0.0) THEN
        ! isolate njj increments for TAF 
        IF (ihop_idw_weights(ii, jj) .EQ. 0.0) THEN
          njj(ii) = IHOP_npts_idw + 1

        END IF

        ! Exactly on a cell center, ignore interpolation
        IF (ihop_idw_weights(ii, jj) .EQ. 0.0) THEN
          tileSSP(iz, ii, bi, bj) = ihop_ssp(i, j, iz-1, bi, bj)

        ! Apply IDW interpolation
        ELSE IF (njj(ii) .LE. IHOP_npts_idw) THEN
          tileSSP(iz, ii, bi, bj) = tileSSP(iz, ii, bi, bj) + &
            ihop_ssp(i, j, iz-1, bi, bj) * &
            ihop_idw_weights(ii, jj) / ihop_sumweights(ii, iz-1)

        ELSE
          ! do nothing
          tileSSP(iz, ii, bi, bj) = tileSSP(iz, ii, bi, bj)

        END IF
      END IF

      ! Extrapolate through bathymetry; don't interpolate
      IF (iz .EQ. SSP%nZ-1 .OR. ihop_sumweights(ii, iz-1) .EQ. 0.0) THEN
        k = iz

        IF (njj(ii) .GE. IHOP_npts_idw) THEN
          ! Determine if you are at the last vlevel
          IF (iz .EQ. SSP%nZ-1 .AND. ihop_sumweights(ii, iz-1) .NE. 0.0) &
              k = k + 1

          ! Calc depth gradient
          dcdz = (tileSSP(k-1, ii, bi, bj) - tileSSP(k-2, ii, bi, bj)) / &
                 (SSP%Z(k-1) - SSP%Z(k-2))
          ! Extrapolate
          tileSSP(k:SSP%nZ, ii, bi, bj) = &
            tileSSP(k-1, ii, bi, bj) + dcdz * SSP%Z(k:SSP%nZ)
          ! Move to next range point, ii
          interp_finished = .TRUE.

        ELSE
          interp_finished = .FALSE.

        END IF
      END IF
    END IF
                END DO !iz
              END IF
            END DO !jj
          END DO !ii
        END DO !i
      END DO !j
    END DO !bi
  END DO !bj

! IESCO24: MITgcm checkpoint69a uses a new global sum subroutine...
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      k = 1
      DO jj = 1,SSP%nR
        DO ii = 1,SSP%nZ
          tmpSSP(bi,bj,k) = tileSSP(ii,jj,bi,bj)
          k = k + 1
        END DO
      END DO
    END DO
  END DO

  CALL GLOBAL_SUM_VECTOR_RL(SSP%nZ*SSP%nR,tmpSSP,globSSP,myThid)

  ! reshape ssp to matrix
  k = 1
  DO jj = 1,SSP%nR
    DO ii = 1,SSP%nZ
      SSP%cMat(ii,jj) = globSSP(k)
      k = k + 1
    END DO
  END DO

  IF(ALLOCATED(tileSSP)) DEALLOCATE(tileSSP)
  IF(ALLOCATED(tmpSSP))  DEALLOCATE(tmpSSP)
  IF(ALLOCATED(globSSP)) DEALLOCATE(globSSP)

  !==================================================
  ! END IDW Interpolate
  !==================================================

  ! set vector structured c, rho, and cz for first range point
  IF (.not. useSSPFile ) THEN ! if usessp, these have already been set
    DO iz = 1,SSP%nZ
      alphaR = SSP%cMat( iz, 1 )

      SSP%c(iz) = CRCI( SSP%Z(iz), alphaR, alphaI, SSP%AttenUnit, bPower, fT, &
                          myThid )
      SSP%rho(iz) = rhoR

      IF ( iz > 1 ) THEN
        IF ( SSP%Z( iz ) .LE. SSP%Z( iz-1 ) ) THEN
#ifdef IHOP_WRITE_OUT
          WRITE( msgBuf,'(2A)' ) 'SSPMOD gcmSSP: ', &
            'The depths in the SSP must be monotone increasing'
          CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
          STOP 'ABNORMAL END: S/R gcmSSP'
        END IF
      END IF

      ! Compute gradient, cz
      IF ( iz>1 ) SSP%cz( iz-1 ) = ( SSP%c( iz ) - SSP%c( iz-1 ) ) / &
                                   ( SSP%Z( iz ) - SSP%Z( iz-1 ) )
    END DO
  END IF

RETURN
END !SUBROUTINE gcmSSP

!**********************************************************************!
SUBROUTINE writeSSP( myThid )
  ! Extracts SSP from MITgcm grid points
  USE ihop_mod, only: PRTFile

  ! == Routine Arguments ==
  ! myThid :: Thread number. Unused by IESCO
  ! msgBuf :: Used to build messages for printing.
  INTEGER, INTENT( IN )   :: myThid
  CHARACTER*(MAX_LEN_MBUF):: msgBuf

  ! == Local Variables ==
  INTEGER           :: iz
  CHARACTER*(80)    :: fmtstr
  REAL (KIND=_RL90) :: ssptmp(SSP%nR)

  ! In adjoint mode we do not write output besides on the first run
  IF ( IHOP_dumpfreq.LT.0) RETURN

  ! init local vars
  ssptmp = 0.0

  ! I/O on main thread only
  _BEGIN_MASTER(myThid)

#ifdef IHOP_WRITE_OUT
  ! Write relevant diagnostics
  WRITE(msgBuf,'(2A)')'________________________________________________', &
    '___________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') "Sound Speed Field"
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(2A)') 'Profile option: ', SSP%Type
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  IF (SSP%nR .GT. 1) THEN
    WRITE(msgBuf,'(A)') 'Using range-dependent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF
  IF (SSP%nR .EQ. 1) THEN
    WRITE(msgBuf,'(A)') 'Using range-independent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF

  WRITE(msgBuf,'(A,I10)') 'Number of SSP ranges = ', SSP%nR
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,I10)') 'Number of SSP depths = ', SSP%nZ
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Profile ranges [km]:'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(fmtStr,'(A,I10,A)') '(T11,',SSP%nR, 'F10.2)'
  ssptmp = SSP%Seg%R( 1:SSP%nR ) / 1000.0
  WRITE(msgBuf,fmtStr) ssptmp
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Sound speed matrix:'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') ' Depth [m ]     Soundspeed [m/s]'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  ssptmp = 0.0
  WRITE(fmtStr,'(A,I10,A)') '(',SSP%nR+1, 'F10.2)'
  DO iz = 1, SSP%nZ
    ssptmp = ssp%cMat( iz,: )
    WRITE(msgBuf,fmtStr) SSP%Z( iz ), ssptmp
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END DO

  IF (useSSPFile) THEN
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Sound speed profile:'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'      z         alphaR      betaR     rho      ', &
                        '  alphaI     betaI'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'     [m]         [m/s]      [m/s]   [g/cm^3]   ', &
                        '   [m/s]     [m/s]'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    WRITE(msgBuf,'(2A)')'_______________________________________________', &
                        '____________'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

    DO iz = 1, SSP%nPts
       WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4)') &
           SSP%Z( iz ), SSP%cMat(iz,1), betaR, rhoR, alphaI, betaI
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END DO
  END IF

  WRITE(msgBuf,'(2A)')'_____________________________________________', &
                      '______________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)')
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER

RETURN
END !SUBROUTINE writeSSP

END MODULE ssp_mod
