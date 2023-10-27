#include "IHOP_OPTIONS.h"
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
  ! since BELLHOP calls it at every step so more information is pre-computed

  USE splinec_mod,  only: cspline, splineall
  USE ihop_mod,     only: PRTFile, SSPFile

  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public EvaluateSSP, HSInfo, Bdry, SSP, zTemp, alphaR, betaR, alphaI, &
           betaI, rhoR, betaPowerLaw, fT, iSegz, iSegr

!=======================================================================

! INPUT/OUTPUT PARAMETERS:

! FUNCTIONS:
  _RL CHEN_MILLERO
  EXTERNAL CHEN_MILLERO

! LOCAL VARIABLES
! == Local Variables ==
  INTEGER bi,bj
  INTEGER i,j

! LEGACY VARIABLES
! == Legacy Local Variables ==
  INTEGER, PARAMETER     :: MaxSSP = 201
  INTEGER                :: iSegr = 1, iSegz = 1
#ifdef IHOP_THREED
  INTEGER                :: iSegx = 1, iSegy = 1
#endif /* IHOP_THREED */
  INTEGER                :: iostat, iallocstat
  INTEGER,           PRIVATE :: iz
  REAL (KIND=_RL90), PRIVATE :: Depth, W
  REAL (KIND=_RL90)          :: zTemp, betaPowerLaw = 1, fT = 1D20
  ! DEFAULT values, BELLHOP only uses alphaR
  REAL (KIND=_RL90)          :: alphaR = 1500, betaR = 0, alphaI = 0, &
                                betaI = 0, rhoR = 1
                            
! TYPE STRUCTURES
! == Type Structures ==
  TYPE rxyz_vector
    REAL (KIND=_RL90), ALLOCATABLE :: r(:), x(:), y(:), z(:)
  END TYPE rxyz_vector

  ! SSP
  TYPE SSPStructure
    INTEGER                 :: NPts, Nr, Nx, Ny, Nz
    REAL    (KIND=_RL90)    :: z( MaxSSP ), rho( MaxSSP )
    COMPLEX (KIND=_RL90)    :: c( MaxSSP ), cz( MaxSSP ), n2( MaxSSP ), &
                               n2z( MaxSSP ), cSpline( 4, MaxSSP )
    REAL    (KIND=_RL90), ALLOCATABLE   :: cMat( :, : ),     czMat( :, : ), &
                                           cMat3( :, :, : ), czMat3( :, :, : )
    TYPE ( rxyz_vector ) :: Seg
    CHARACTER (LEN=1)    :: Type
    CHARACTER (LEN=2)    :: AttenUnit
    ! for PCHIP coefs.
    COMPLEX (KIND=_RL90)    :: cCoef( 4, MaxSSP ), CSWork( 4, MaxSSP )   
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP

  ! *** Halfspace properties structure ***

  TYPE HSInfo
     ! compressional and shear wave speeds/attenuations in user units
     REAL   (KIND=_RL90)    :: alphaR, alphaI, betaR, betaI    
     REAL   (KIND=_RL90)    :: rho, Depth        ! density, depth
     COMPLEX(KIND=_RL90)    :: cP, cS            ! P-wave, S-wave speeds
     CHARACTER (LEN=1)      :: BC                ! Boundary condition type
     CHARACTER (LEN=6)      :: Opt
  END TYPE HSInfo

  TYPE BdryPt
     TYPE( HSInfo )   :: HS
  END TYPE

  TYPE BdryType
     TYPE( BdryPt )   :: Top, Bot
  END TYPE BdryType

  TYPE(BdryType) :: Bdry
!EOP

CONTAINS
  SUBROUTINE EvaluateSSP( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, &
      myThid )

    ! Call the particular profile routine indicated by the SSP%Type and 
    ! perform Task
    !   Task = 'TAB' to tabulate cp, cs, rhoT 
    !   Task = 'INI' to initialize

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: freq
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    CHARACTER( LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
#ifdef IHOP_THREED
    REAL (KIND=_RL90)                :: gradc_3d( 3 ), cxx, cyy, cxy, cxz, cyz
    REAL (KIND=_RL90)                :: x3( 3 )
#endif /* IHOP_THREED */

    SELECT CASE ( SSP%Type )
    CASE ( 'N' )  !  N2-linear profile option
       CALL n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
    CASE ( 'C' )  !  C-linear profile option
       CALL cLinear(  x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
    CASE ( 'P' )  !  monotone PCHIP ACS profile option
       CALL cPCHIP(   x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
    CASE ( 'S' )  !  Cubic spline profile option
       CALL cCubic(   x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
    CASE ( 'Q' )
       CALL Quad(     x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
    CASE ( 'A' )  !  Analytic profile option
       CALL Analytic( x, c, cimag, gradc, crr, crz, czz, rho, Task, myThid )
    CASE DEFAULT
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'Profile option: ', SSP%Type
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') 'SSPMOD EvaluateSSP: ', 'Invalid SSP profile option'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R EvaluateSSP'
    END SELECT

  RETURN
  END !SUBROUTINE EvaluateSSP
  
!**********************************************************************!

  SUBROUTINE n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )

    ! N2-linear interpolation of SSP data

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: freq
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    
    IF ( Task == 'INI' ) THEN   ! read in SSP data
       ! *** Task 'INI' for initialization ***

       Depth = x( 2 )
       CALL ReadSSP( Depth, freq, myThid )
              
       SSP%n2(  1 : SSP%NPts ) = 1.0 / SSP%c( 1 : SSP%NPts )**2
       !IEsco23 Test this: SSP%n2(  1 : SSP%Nz ) = 1.0 / SSP%c( 1 : SSP%Nz )**2

       ! compute gradient, n2z
       DO iz = 2, SSP%Npts
          SSP%n2z( iz - 1 ) = ( SSP%n2(   iz ) - SSP%n2(   iz - 1 ) ) / &
                              ( SSP%z(    iz ) - SSP%z(    iz - 1 ) )
       END DO
    ELSE                         ! return SSP info

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketting Depths
!IEsco23 Test this: 
!          DO iz = 2, SSP%Nz   ! Search for bracketting Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       W = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )

       c     = REAL(  1.0D0 / SQRT( ( 1.0D0 - W ) * SSP%n2( iSegz ) &
               + W * SSP%n2( iSegz + 1 ) ) )
       cimag = AIMAG( 1.0D0 / SQRT( ( 1.0D0 - W ) * SSP%n2( iSegz ) &
               + W * SSP%n2( iSegz + 1 ) ) )

       gradc = [ 0.0D0, -0.5D0 * c * c * c * REAL( SSP%n2z( iSegz ) ) ]
       crr   = 0.0d0
       crz   = 0.0d0
       czz   = 3.0d0 * gradc( 2 ) * gradc( 2 ) / c

       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  RETURN
  END !SUBROUTINE n2Linear

  !**********************************************************************!

  SUBROUTINE cLinear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid  )

    ! c-linear interpolation of SSP data

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: freq
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    
    IF ( Task == 'INI' ) THEN   ! read in SSP data
       ! *** Task 'INI' for initialization ***

       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq, myThid )
    ELSE                        ! return SSP info

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketting Depths
!IEsco23 Test this: 
!          DO iz = 2, SSP%Nz   ! Search for bracketting Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       c     = REAL(  SSP%c( iSegz ) + ( x( 2 ) - SSP%z( iSegz ) ) * SSP%cz( iSegz ) )
       cimag = AIMAG( SSP%c( iSegz ) + ( x( 2 ) - SSP%z( iSegz ) ) * SSP%cz( iSegz ) )
       gradc = [ 0.0D0, REAL( SSP%cz( iSegz ) ) ]
       crr   = 0.0d0
       crz   = 0.0d0
       czz   = 0.0d0

       W     = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  RETURN
  END !SUBROUTINE cLinear

  !**********************************************************************!

  SUBROUTINE cPCHIP( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )

    ! This implements the monotone piecewise cubic Hermite interpolating
    ! polynomial (PCHIP) algorithm for the interpolation of the sound speed c.

    USE pchip_mod,  only: PCHIP

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: freq
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    REAL    (KIND=_RL90) :: xt
    COMPLEX (KIND=_RL90) :: c_cmplx

    IF ( Task == 'INI' ) THEN   ! read in SSP data
       ! *** Task 'INI' for initialization ***

       Depth = x( 2 )
       CALL ReadSSP( Depth, freq, myThid )

       !                                                               2      3
       ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
       !
       CALL PCHIP( SSP%z, SSP%c, SSP%NPts, SSP%cCoef, SSP%CSWork )
!IEsco23 Test this: 
!       CALL PCHIP( SSP%z, SSP%c, SSP%Nz, SSP%cCoef, SSP%CSWork )

    ELSE    ! return SSP info, recall iSegz is initiated to 1

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketting Depths
!IEsco23 Test this: 
!          DO iz = 2, SSP%Nz   ! Search for bracketting Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       xt = x( 2 ) - SSP%z( iSegz )
       c_cmplx = SSP%cCoef( 1, iSegz ) &
             + ( SSP%cCoef( 2, iSegz ) &
             + ( SSP%cCoef( 3, iSegz ) &
             +   SSP%cCoef( 4, iSegz ) * xt ) * xt ) * xt

       c     = REAL(  c_cmplx )
       cimag = AIMAG( c_cmplx )

       gradc = [ 0.0D0, &
                 REAL( SSP%cCoef( 2, iSegz ) + ( 2.0D0 * SSP%cCoef( 3, iSegz ) &
                       + 3.0D0 * SSP%cCoef( 4, iSegz ) * xt ) * xt ) ]

       crr   = 0.0D0
       crz   = 0.0D0
       czz   = REAL( 2.0D0 * SSP%cCoef( 3, iSegz ) + &
                     6.0D0 * SSP%cCoef( 4, iSegz ) * xt )   ! dgradc(2)/dxt

       W     = ( x( 2 ) - SSP%z( iSegz ) ) / &
               ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       ! linear interp of density
       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )  

    END IF

  RETURN
  END !SUBROUTINE cPCHIP

  !**********************************************************************!

  SUBROUTINE cCubic( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid  )

    ! Cubic spline interpolation

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
    INTEGER, INTENT( IN )   :: myThid
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN )  :: freq
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, &
                                        rho ! sound speed and its derivatives
    INTEGER                 :: iBCBeg, iBCEnd
    REAL     (KIND=_RL90)   :: hSpline
    COMPLEX  (KIND=_RL90)   :: c_cmplx, cz_cmplx, czz_cmplx
    
    IF ( Task == 'INI' ) THEN
       ! *** Task 'INI' for initialization ***
       
       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq, myThid )

       SSP%cSpline( 1, 1 : SSP%NPts ) = SSP%c( 1 : SSP%NPts )
!IEsco23 Test this: 
!       SSP%cSpline( 1, 1 : SSP%Nz ) = SSP%c( 1 : SSP%Nz )
       
       ! Compute spline coefs
       iBCBeg = 0
       iBCEnd = 0
       CALL CSpline( SSP%z, SSP%cSpline( 1, 1 ), SSP%NPts, iBCBeg, iBCEnd, SSP%NPts )
!IEsco23 Test this: 
!       CALL CSpline( SSP%z, SSP%cSpline( 1, 1 ), SSP%Nz, iBCBeg, iBCEnd, SSP%Nz )
    ELSE

       ! *** Section to return SSP info ***

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketting Depths
!IEsco23 Test this: 
!          DO iz = 2, SSP%Nz   ! Search for bracketting Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO

       END IF

       hSpline = x( 2 ) - SSP%z( iSegz )

       ! c   = Spline(   SSP%cSpline( 1, iSegz ), hSpline )
       ! cz  = SplineX(  SSP%cSpline( 1, iSegz ), hSpline )
       ! czz = SplineXX( SSP%cSpline( 1, iSegz ), hSpline )

       CALL SplineALL( SSP%cSpline( 1, iSegz ), hSpline, c_cmplx, cz_cmplx, czz_cmplx )

       c     = DBLE(  c_cmplx )
       cimag = AIMAG( c_cmplx )
       gradc = [ 0.0D0, DBLE( cz_cmplx ) ]
       czz   = DBLE( czz_cmplx )
       crr   = 0.0d0
       crz   = 0.0d0

       ! linear interpolation for density
       W   = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  RETURN
  END !SUBROUTINE cCubic

!**********************************************************************!

SUBROUTINE Quad( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task, myThid )
  ! Bilinear quadrilatteral interpolation of SSP data in 2D, SSP%Type = 'Q'

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
  INTEGER, INTENT( IN )   :: myThid
  CHARACTER*(MAX_LEN_MBUF):: msgBuf

  !     == Local Variables ==
  REAL (KIND=_RL90), INTENT( IN  ) :: freq
  REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )  ! r-z SSP evaluation point
  CHARACTER (LEN=3), INTENT( IN  ) :: Task
  REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, &
                                      rho ! sound speed and its derivatives
  INTEGER             :: irT, iz2
  REAL (KIND=_RL90)   :: c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z
  
  IF ( Task == 'INI' ) THEN
    ! *** Task 'INI' for initialization ***
      
    Depth = x( 2 )
    ! Defensive check if SSPFile exists
    IF (useSSPFile .EQV. .TRUE.) THEN ! eqv for logical operands
      CALL ReadSSP( Depth, freq, myThid )
    ELSE
      CALL ExtractSSP(Depth, freq, myThid )
    END IF

    ! calculate cz
    DO irT = 1, SSP%Nr
      DO iz2 = 2, SSP%Nz
        delta_z = ( SSP%z( iz2 ) - SSP%z( iz2-1 ) )
        SSP%czMat( iz2-1, irT ) = ( SSP%cMat( iz2  , irT ) - &
                                    SSP%cMat( iz2-1, irT ) ) / delta_z
      END DO
    END DO

    RETURN

  ELSE ! Task == 'TAB'
    ! *** Section to return SSP info ***

    ! IESCO22: iSegz is the depth index containing x depth
    ! find depth-layer where x(2) in ( SSP%z( iSegz ), SSP%z( iSegz+1 ) )
    IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
      DO iz = 2, SSP%Nz   ! Search for bracketting Depths
        IF ( x( 2 ) < SSP%z( iz ) ) THEN
          iSegz = iz - 1
          EXIT
        END IF
      END DO
    END IF

    ! Check that x is inside the box where the sound speed is defined
    IF ( x( 1 ) < SSP%Seg%r( 1 ) .OR. x( 1 ) > SSP%Seg%r( SSP%Nr ) ) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf,'(2A)') 'ray is outside the box where ocean ',&
        'soundspeed is defined'
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,2F13.4)') ' x = ( r, z ) = ', x
      CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
        'ray is outside the box where the soundspeed is defined'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R Quad'
    END IF

    ! find range-segment where x(1) in [ SSP%Seg%r( iSegr ), SSP%Seg%r( iSegr+1 ) )
    IF ( x( 1 ) < SSP%Seg%r( iSegr ) .OR. x( 1 ) >= SSP%Seg%r( iSegr + 1 ) ) THEN
      DO irT = 2, SSP%Nr   ! Search for bracketting segment ranges
        IF ( x( 1 ) < SSP%Seg%r( irT ) ) THEN
          iSegr = irT - 1
          EXIT
        END IF
      END DO
    END IF

    ! for depth, x(2), get the sound speed at both ends of range segment
    cz1 = SSP%czMat( iSegz, iSegr   )
    cz2 = SSP%czMat( iSegz, iSegr+1 )

    ! IESCO22: s2 is distance btwn field point, x(2), and ssp depth @ iSegz
    s2      = x( 2 )           - SSP%z( iSegz )            
    delta_z = SSP%z( iSegz+1 ) - SSP%z( iSegz )
    IF (delta_z <= 0 .OR. s2 > delta_z) THEN
#ifdef IHOP_WRITE_OUT
      WRITE(msgBuf, *) delta_z, s2, iSegz, SSP%z(iSegz)
      CALL PRINT_ERROR( msgBuf,myThid )
      WRITE(msgBuf,'(2A)') 'SSPMOD Quad: ', &
        'depth is not monotonically increasing in SSP%z'
      CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R Quad'
    END IF
    
    c1 = SSP%cMat( iSegz, iSegr   ) + s2*cz1
    c2 = SSP%cMat( iSegz, iSegr+1 ) + s2*cz2

    ! s1 = proportional distance of x(1) in range
    delta_r = SSP%Seg%r( iSegr+1 ) - SSP%Seg%r( iSegr )
    s1 = ( x( 1 ) - SSP%Seg%r( iSegr ) ) / delta_r
    s1 = MIN( s1, 1.0D0 )   ! piecewise constant extrapolation for ranges outside SSPFile box
    s1 = MAX( s1, 0.0D0 )   ! "

    c = ( 1.0D0-s1 )*c1 + s1*c2 ! c @ x

    ! interpolate the attenuation !!!! SSP in ENVFile needs to match first column of SSPFile
    s2    = s2 / delta_z   ! normalize depth layer
    cimag = AIMAG( ( 1.0D0-s2 )*SSP%c( Isegz ) + s2*SSP%c( Isegz+1 ) )   ! volume attenuation is taken from the single c(z) profile

    cz  = ( 1.0D0-s1 )*cz1 + s1*cz2 ! cz @ x

    cr  = ( c2  - c1  ) / delta_r ! SSPFile grid cr
    crz = ( cz2 - cz1 ) / delta_r ! SSPFile grid crz

    gradc = [ cr, cz ]
    crr   = 0.0D0
    czz   = 0.0D0

    ! linear interpolation for density
    W   = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
    rho = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
  END IF

  !IESCO22: for thesis, czz=crr=0, and rho=1 at all times
RETURN
END !SUBROUTINE Quad

!**********************************************************************!

  SUBROUTINE Analytic( x, c, cimag, gradc, crr, crz, czz, rho, Task, myThid )

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT( IN  ) :: x( 2 )
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL (KIND=_RL90), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
    REAL (KIND=_RL90)                :: c0, cr, cz, DxtDz, xt

    IF ( Task == 'INI' ) THEN
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(A)') 'Analytic SSP option'
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
       SSP%NPts = 2
       SSP%z(1) = 0.0
       SSP%z(2) = Bdry%Bot%HS%Depth
    END IF

    iSegz = 1
    c0    = 1500.0
    rho   = 1.0

    ! homogeneous halfspace was removed since BELLHOP needs to get gradc just a little below the boundaries, on ray reflection

!!$  IF ( x( 2 ) < 5000.0 ) THEN
    xt    = 2.0 * ( x( 2 ) - 1300.0 ) / 1300.0
    DxtDz = 2.0 / 1300.0
    c     = C0 * ( 1.0 + 0.00737*( xt - 1.0 + EXP( -xt ) ) )
    cimag = 0.
    cz    = C0 * 0.00737 * ( 1.0 - EXP( -xt ) ) * DxtDz
    czz   = C0 * 0.00737 * EXP( -xt ) * DxtDz ** 2
!!$  ELSE
!!$     ! Homogeneous half-space
!!$     xt   = 2.0 * ( 5000.0 - 1300.0 ) / 1300.0
!!$     c    = C0 * ( 1.0 + 0.00737 * ( xt - 1.0 + EXP( -xt ) ) )
!!$     cz   = 0.0
!!$     czz  = 0.0
!!$  END IF

    cr = 0.0
    gradc = [ cr, cz ]
    crz = 0.0
    crr = 0.0

  RETURN
  END !SUBROUTINE Analytic

!**********************************************************************!

  SUBROUTINE ReadSSP( Depth, freq, myThid )
    ! reads SSP in m/s from .ssp file and convert to AttenUnit (ie. Nepers/m)
    ! Populates SSPStructure: SSP

    USE atten_mod, only: CRCI

  !     == Routine Arguments ==
  !     myThid :: Thread number. Unused by IESCO
  !     msgBuf :: Used to build messages for printing.
    INTEGER, INTENT( IN )   :: myThid
    CHARACTER*(MAX_LEN_MBUF):: msgBuf
  
  !     == Local Variables ==
    REAL (KIND=_RL90), INTENT(IN) :: Depth, freq
    INTEGER :: iz2

    ! I/O on main thread only
    _BEGIN_MASTER(myThid)
    ! OPEN SSPFile to read
    OPEN ( FILE = TRIM( IHOP_fileroot ) // '.ssp', UNIT = SSPFile, &
        FORM = 'FORMATTED', STATUS = 'OLD', IOSTAT = iostat )
    IF ( IOSTAT /= 0 ) THEN   ! successful open?
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'SSPFile = ', TRIM( IHOP_fileroot ) // '.ssp'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') 'SSPMOD ReadSSP: Unable to open the SSP file'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSSP'
    END IF

    ! Write relevant diagnostics
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') "Sound Speed Field" 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(2A)')'____________________________________________________',&
                        '_______'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    READ( SSPFile,  * ) SSP%Nr, SSP%Nz
#ifdef IHOP_WRITE_OUT
    IF (SSP%Nr .GT. 1) THEN
        WRITE(msgBuf,'(A)') 'Using range-dependent sound speed'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF
    IF (SSP%Nr .EQ. 1) THEN
        WRITE(msgBuf,'(A)') 'Using range-independent sound speed'
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    END IF

    WRITE(msgBuf,'(A,I10)') 'Number of SSP ranges = ', SSP%Nr
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A,I10)') 'Number of SSP depths = ', SSP%Nz
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

    ALLOCATE( SSP%cMat( SSP%Nz, SSP%Nr ), &
              SSP%czMat( SSP%Nz-1, SSP%Nr ), &
              SSP%Seg%r( SSP%Nr ), &
              STAT = iallocstat )
    IF ( iallocstat /= 0 ) THEN
#ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
                            'Insufficient memory to store SSP'
        CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ReadSSP'
    END IF

    READ( SSPFile,  * ) SSP%Seg%r( 1 : SSP%Nr )
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Profile ranges (km):'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(F10.2)') SSP%Seg%r( 1 : SSP%Nr )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    SSP%Seg%r = 1000.0 * SSP%Seg%r   ! convert km to m

    READ( SSPFile,  * ) SSP%z( 1 : SSP%Nz )
!#ifdef IHOP_DEBUG
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Profile depths (m):'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(F10.2)') SSP%z( 1 : SSP%Nz )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
!#endif

    ! IEsco23: contain read of ssp in this subroutine only 
    ! IEsco23: change to allocatable memory since we should know Nz
#ifdef IHOP_DEBUG
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Sound speed matrix:'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') ' Depth (m )     Soundspeed (m/s)'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
#endif /* IHOP_DEBUG */
    DO iz2 = 1, SSP%Nz
       READ(  SSPFile, * ) SSP%cMat( iz2, : )
#ifdef IHOP_DEBUG
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'(12F10.2)') SSP%z( iz2 ), SSP%cMat( iz2, : )
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
#endif /* IHOP_DEBUG */
    END DO
    CLOSE( SSPFile )

#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 'Sound speed profile:'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')'      z         alphaR      betaR     rho        alphaI     betaI'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)')'     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    
    WRITE(msgBuf,'(2A)')'____________________________________________________',&
                        '_______'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
    WRITE(msgBuf,'(A)') 
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */
    SSP%NPts = 1
    DO iz = 1, MaxSSP 
       alphaR = SSP%cMat( iz, 1 )
#ifdef IHOP_WRITE_OUT
       WRITE(msgBuf,'( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4)') &
           SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI
       CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
#endif /* IHOP_WRITE_OUT */

       SSP%c(   iz ) = CRCI( SSP%z( iz ), alphaR, alphaI, freq, freq, &
                             SSP%AttenUnit, betaPowerLaw, fT, myThid )
       SSP%rho( iz ) = rhoR !IEsco22: set to a default value of 1

       ! verify depths are monotone increasing
       IF ( iz > 1 ) THEN
          IF ( SSP%z( iz ) .LE. SSP%z( iz - 1 ) ) THEN
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(A,F10.2)') 'Bad depth in SSP: ', SSP%z( iz )
                CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
                WRITE(msgBuf,'(2A)') 'SSPMOD ReadSSP: ', &
                            'The depths in the SSP must be monotone increasing'
                CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
                STOP 'ABNORMAL END: S/R ReadSSP'
          END IF
       END IF

       ! compute gradient, cz
       IF ( iz > 1 ) SSP%cz( iz - 1 )  = ( SSP%c( iz ) - SSP%c( iz - 1 ) ) / &
                                         ( SSP%z( iz ) - SSP%z( iz - 1 ) )

       ! Did we read the last point?
       IF ( ABS( SSP%z( iz ) - Depth ) < 100. * EPSILON( 1.0e0 ) ) THEN
          IF ( SSP%NPts == 1 ) THEN
#ifdef IHOP_WRITE_OUT
                WRITE(msgBuf,'(A,I10)') '#SSP points: ', SSP%NPts
                CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
                WRITE(msgBuf,'(2A)')  'SSPMOD ReadSSP: ', &
                                    'The SSP must have at least 2 points'
                CALL PRINT_ERROR( msgBuf,myThid )
#endif /* IHOP_WRITE_OUT */
                STOP 'ABNORMAL END: S/R ReadSSP'
          END IF

          RETURN
       ENDIF

       SSP%NPts = SSP%NPts + 1
    END DO
 
    ! Fall through means too many points in the profile
#ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(A,I10)') 'Max. #SSP points: ', MaxSSP
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
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

SUBROUTINE ExtractSSP( Depth, freq, myThid )
  ! Extracts SSP from MITgcm grid points

  use atten_mod, only: CRCI

  ! == Routine Arguments ==
  ! myThid :: Thread number. Unused by IESCO
  ! msgBuf :: Used to build messages for printing.
  INTEGER, INTENT(IN)     :: myThid
  CHARACTER*(MAX_LEN_MBUF):: msgBuf
  CHARACTER*(80)          :: fmtstr

  ! == Local Variables ==
  INTEGER                       :: ii, jj, njj(IHOP_NPTS_RANGE), k
  REAL (KIND=_RL90), INTENT(IN) :: Depth, freq
  REAL (KIND=_RL90)             :: sumweights(IHOP_NPTS_RANGE, Nr), &
                                   dcdz
  REAL (KIND=_RL90), ALLOCATABLE:: tmpSSP(:,:,:,:)

  SSP%Nz = Nr+2 ! add z=0 z=Depth layers 
  SSP%Nr = IHOP_NPTS_RANGE

  ALLOCATE( SSP%cMat( SSP%Nz, SSP%Nr ), &
            SSP%czMat( SSP%Nz-1, SSP%Nr ), &
            SSP%Seg%r( SSP%Nr ), tmpSSP(SSP%Nz,SSP%Nr,nSx,nSy),&
            STAT = iallocstat )
  IF ( iallocstat /= 0 ) THEN
# ifdef IHOP_WRITE_OUT
    WRITE(msgBuf,'(2A)') 'SSPMOD ExtractSSP: ', &
      'Insufficient memory to store SSP'
    CALL PRINT_ERROR( msgBuf,myThid )
# endif /* IHOP_WRITE_OUT */
      STOP 'ABNORMAL END: S/R ExtractSSP'
  END IF

  ! Initiate SSP%cMat to ceros
  SSP%cMat  = 0.0 _d 0
  tmpSSP    = 0.0 _d 0
  njj       = 0
  dcdz      = 0.0 _d 0

  ! set SSP%Seg%r from data.ihop -> ihop_ranges
  SSP%Seg%r( 1:SSP%Nr ) = ihop_ranges( 1:SSP%Nr )

  ! set SSP%z from rC, rkSign=-1 used bc ihop uses +ive depths
  SSP%z( 1 )            = 0.0 _d 0
  SSP%z( 2:(SSP%Nz-1) ) = rkSign*rC( 1:Nr )
  SSP%z( SSP%Nz )       = Bdry%Bot%HS%Depth ! rkSign*rF(Nr+1)*1.05

  !==================================================
  ! IDW Interpolate: COMPARING with LAT LONs (xC, yC)
  !==================================================
  ! Sum IDW weights
  DO ii = 1,IHOP_npts_range
    sumweights(ii,:) = sum(ihop_idw_weights(ii,:))
  END DO 

  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      DO j=1,sNy
        DO i=1,sNx
          DO ii=1,IHOP_npts_range
            ! IDW Interpolation weight sum
            DO jj=1,IHOP_npts_idw
              IF (xC(i,j,bi,bj) .eq. ihop_xc(ii,jj) .and. &
                  yC(i,j,bi,bj) .eq. ihop_yc(ii,jj)) THEN 
                DO k=1,Nr
                  IF ( hFacC(i,j,k,bi,bj).eq.0.0 ) THEN
                    sumweights(ii,k) = sumweights(ii,k) - ihop_idw_weights(ii,jj)
                  ENDIF
                  ! Set TINY values to 0.0
                  IF ( ABS(sumweights(ii,k)).LT.1D-24 ) sumweights(ii,k) = 0.0
                ENDDO
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! from ocean grid to acoustic grid with IDW
  DO bj=myByLo(myThid),myByHi(myThid)
    DO bi=myBxLo(myThid),myBxHi(myThid)
      DO j=1,sNy
        DO i=1,sNx
          DO ii=1,IHOP_npts_range
            ! IDW Interpolate SSP at second order
            interp: DO jj=1,IHOP_npts_idw
            IF (xC(i,j,bi,bj) .eq. ihop_xc(ii,jj) .and. &
                yC(i,j,bi,bj) .eq. ihop_yc(ii,jj)) THEN
              njj(ii) = njj(ii) + 1

              DO iz=1,SSP%Nz-1
                IF (iz.eq.1) THEN
                  ! Top vlevel zero depth
                  tmpSSP(1,ii,bi,bj) = tmpSSP(1,ii,bi,bj) + &
                    CHEN_MILLERO(i,j,0,bi,bj,myThid)* &
                    ihop_idw_weights(ii,jj)/sumweights(ii,iz)
                ELSE ! 2:(SSP%Nz-1)
                  ! Middle depth layers
                  IF (sumweights(ii,iz-1).gt.0.0) THEN
                    IF (ihop_idw_weights(ii,jj).eq.0.0) THEN
                      ! Exactly on a cell center, ignore interpolation
                      tmpSSP(iz,ii,bi,bj) = ihop_ssp(i,j,iz-1,bi,bj)
                      EXIT interp
                    ELSE
                      tmpSSP(iz,ii,bi,bj) = tmpSSP(iz,ii,bi,bj) + &
                        ihop_ssp(i,j,iz-1,bi,bj)* &
                        ihop_idw_weights(ii,jj)/sumweights(ii,iz-1)
                    END IF
                  END IF 
                  
                  IF ( iz.eq.SSP%Nz-1 .or. sumweights(ii,iz-1).eq.0.0 ) THEN 
                    ! Extrapolate through bathymetry
                    k=iz
                    
                    IF ( njj(ii).eq.IHOP_npts_idw ) THEN 
                      ! Last gcm vlevel extrapolation
                      IF ( iz.eq.SSP%Nz-1 .and. sumweights(ii,iz-1).ne.0.0 ) k = k+1
                       
                      ! Calc depth gradient
                      dcdz =  ( tmpSSP(k-1,ii,bi,bj)-tmpSSP(k-2,ii,bi,bj) ) / &
                              ( SSP%z(k-1)-SSP%z(k-2) )
                      ! Extrapolate
                      tmpSSP(k:SSP%Nz,ii,bi,bj) = &
                          tmpSSP(k-1,ii,bi,bj) + dcdz*SSP%z(k:SSP%Nz)
                      ! Move to the next range point, ii
                      EXIT interp
                    END IF 
                  END IF

                END IF
              END DO  
            
            END IF
            END DO interp
          END DO 
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  CALL GLOBAL_VEC_SUM_R8(SSP%Nz*SSP%Nr,SSP%Nz*SSP%Nr,tmpSSP,myThid)
  SSP%cMAT = tmpSSP(:,:,1,1)
  !==================================================
  ! END IDW Interpolate
  !==================================================

  ! I/O on main thread only
  _BEGIN_MASTER(myThid)
  ! set vector structured c, rho, and cz for first range point
  DO iz = 1,SSP%Nz
    alphaR = SSP%cMat( iz, 1 )

    SSP%c(   iz ) = CRCI( SSP%z( iz ), alphaR, alphaI, freq, freq, &
                          SSP%AttenUnit, betaPowerLaw, fT, myThid )
    SSP%rho( iz ) = rhoR

    IF ( iz > 1 ) THEN
      IF ( SSP%z( iz ) .LE. SSP%z( iz-1 ) ) THEN
#  ifdef IHOP_WRITE_OUT
        WRITE(msgBuf,'(A)') 'Bad depth in SSP: ', SSP%z(iz)
        CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
        WRITE( msgBuf,'(2A)' ) 'SSPMOD ExtractSSP: ', &
          'The depths in the SSP must be monotone increasing'
        CALL PRINT_ERROR( msgBuf,myThid )
#  endif /* IHOP_WRITE_OUT */
        STOP 'ABNORMAL END: S/R ExtractSSP'
      END IF
    END IF

    ! Compute gradient, cz
    IF ( iz>1 ) SSP%cz( iz-1 ) = ( SSP%c( iz ) - SSP%c( iz-1 ) ) / &
                                  ( SSP%z( iz ) - SSP%z( iz-1 ) )
  END DO

  ! Write relevant diagnostics
#ifdef IHOP_WRITE_OUT
  WRITE(msgBuf,'(2A)')'________________________________________________', &
    '___________'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') "Sound Speed Field" 
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  IF (SSP%Nr.GT.1) THEN
    WRITE(msgBuf,'(A)') 'Using range-dependent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF
  IF (SSP%Nr.EQ.1) THEN
    WRITE(msgBuf,'(A)') 'Using range-independent sound speed'
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END IF

  WRITE(msgBuf,'(A,I10)') 'Number of SSP ranges = ', SSP%Nr
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A,I10)') 'Number of SSP depths = ', SSP%Nz
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )

  WRITE(msgBuf,'(A)') 
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Profile ranges (km):'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(fmtStr,'(A,I10,A)') '(T11,',SSP%Nr, 'F10.2)'
  WRITE(msgBuf,fmtStr) SSP%Seg%r( 1:SSP%Nr )
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') 'Sound speed matrix:'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  WRITE(msgBuf,'(A)') ' Depth (m)     Soundspeed (m/s)'
  CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  DO iz = 1, SSP%Nz
    WRITE(msgBuf,'(12F10.2)'  ) SSP%z( iz ), SSP%cMat( iz, : )
    CALL PRINT_MESSAGE( msgbuf, PRTFile, SQUEEZE_RIGHT, myThid )
  END DO
#endif /* IHOP_WRITE_OUT */
  ! I/O on main thread only
  _END_MASTER(myThid)
  _BARRIER

  SSP%Seg%r = 1000.0 * SSP%Seg%r   ! convert km to m

RETURN
END !SUBROUTINE ExtractSSP

END MODULE ssp_mod
