#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: splinec_mod
MODULE SPLINEC_MOD
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module contains subroutines for cubic spline interpolation.

! !USES:
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC cspline, splineall
!=======================================================================

! == Module variables == None
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R CSPLINE
! S/R SPLINEALL
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: CSPLINE
! !INTERFACE:
  SUBROUTINE CSPLINE (TAU, C, N, IBCBEG, IBCEND, NDIM)
! !DESCRIPTION:
!  TAKEN FROM "A PRACTICAL GUIDE TO SPLINES", BY CARL DE BOOR. 1978.
!  SPRINGER-VERLAG.  THE INPUT PARAMETER "NDIM" HAS BEEN ADDED TO
!  ALLOW FOR MULTIPLE CALLS WITH DIFFERENT VALUES OF N. - DENNIS DUNDORE

!  SUBSTANTIAL MODIFICATIONS MADE BY STEVE WALES, APRIL 1983,
!  PRINCIPALLY TO HANDLE COMPLEX NUMBERS (C) & UPDATE THE FORTRAN.

! *****************************  I N P U T  ****************************

!  N = NUMBER OF DATA POINTS.  ASSUMED TO BE .GE. 2.

!  (TAU(I), C(1,I),I=1,...,N) = ABSCISSAE AND ORDINATES OF THE DATA
!      POINTS.  TAU IS ASSUMED TO BE STRICTLY INCREASING.

!  IBCBEG, IBCEND = BOUNDARY CONDITION INDICATORS, AND
!  C(2,1), C(2,N) = BOUNDARY CONDITION INFORMATION.  SPECIFICALLY,
!      IBCBEG = 0 IMPLIES NO BOUNDARY CONDITION AT TAU(1) IS GIVEN.
!            IN THIS CASE, THE "NOT-A-KNOT" CONDITION IS USED, IE THE
!            JUMP IN THE 3-RD DERIVATIVE ACROSS TAU(2) IS FORCED TO 0.,
!            THUS THE 1-ST AND 2-ND POLYNOMIAL PIECES ARE MADE TO
!            COINCIDE.
!      IBCBEG = 1 IMPLIES THAT THE SLOPE AT TAU(1) IS SET EQUAL TO C(2,1)
!            INPUT BY USER.
!      IBCBEG = 2 IMPLIES THAT THE 2-ND DERIVATIVE AT TAU(1) IS SET EQUAL
!            TO C(2,1), SUPPLIED BY INPUT.
!      IBCEND = 0, 1, OR 2 HAS ANALOGOUS MEANING CONCERNING THE BOUNDARY
!            CONDITION AT TAU(N), WITH INFORMATION SUPPLIED BY C(2,N).

!  NDIM = ROW DIMENSION OF ! MATRIX:  C(4,NDIM)

! **************************  O U T P U T  ****************************

!  C(J,I), J=1,...,4;  I=1,...,L=N-1  =  THE POLY COEFFS OF THE CUBI!
!      SPLINE WITH INTERIOR KNOTS TAU(2),...,TAU(N-1).  PRECISELY, IN THE
!      INTERVAL (TAU(I), TAU(I+1)), THE SPLINE F IS GIVEN BY

!       F(X) = C(1,I) + H*(C(2,I) + H*(C(3,I)/2. + H*C(4,I)/6.))

!      WHERE H = X - TAU(I).

!     THE COEFFICIENTS CALCULATED ARE, 1) THE VALUE, 2) THE SLOPE, AND
!     3) THE CURVATURE AT EACH OF THE KNOTS 1 TO N-1, AND 4) THE RATE OF
!     CHANGE OF THE CURVATURE OVER THE FOLLOWING INTERVAL.  IN ADDITION,
!     WE HAVE THE VALUE AND THE SLOPE AT THE LAST POINT. THE LAST TWO
!     POSTIONS AT THE LAST POINT ARE THEN SET TO THE CURVATURE AT THAT
!     POINT (IN C(3,N)) AND THE MEAN VALUE OVER THE ENTIRE INTERVAL,
!     CALCULATED AS THE INTEGRAL OVER THE INTERVAL DIVIDED BY THE LENGTH
!     OF THE INTERVAL (IN C(4,N)).

! !USES:
  IMPLICIT INTEGER            (I-M)
  IMPLICIT REAL (KIND=_RL90)  (A-H,O-Z)

! !INPUT PARAMETERS:
! TAU :: vector of abscissae, length N
! C   :: matrix of coefficients, 4 rows, NDIM columns
! N   :: number of data points, assumed to be >= 2
! IBCBEG :: boundary condition indicator at the beginning of the interval
! IBCEND :: boundary condition indicator at the end of the interval
! NDIM :: row dimension of the matrix C
  INTEGER,              INTENT(IN) :: N, NDIM, IBCBEG, IBCEND
  REAL    (KIND=_RL90), INTENT(IN) :: TAU(N)
  COMPLEX (KIND=_RL90), INTENT(INOUT) :: C(4,NDIM)
! !OUTPUT PARAMETERS: C

! !LOCAL VARIABLES:
! G :: temporary complex variable
! DTAU :: difference in abscissae
! DIVDF1, DIVDF3 :: temporary complex variables for derivatives
  COMPLEX (KIND=_RL90) :: G, DTAU, DIVDF1, DIVDF3
!EOP

  L = N - 1

!$TAF init cspline1 = static, 4*ndim*L

  DO M = 2,N
    C(3,M) = TAU(M) - TAU(M-1)
    C(4,M) = (C(1,M) - C(1,M-1)) / C(3,M)
  ENDDO

  !   * BEGINNING BOUNDARY CONDITION SECTION *

!$TAF store C = cspline1
  IF ( IBCBEG.EQ.0 ) THEN
!$TAF store C = cspline1
    IF ( N.GT.2 ) THEN
      C(4,1) = C(3,3)
      C(3,1) = C(3,2) + C(3,3)
      C(2,1) = ((C(3,2) + 2.0*C(3,1))*C(4,2)*C(3,3) + &
            C(3,2)**2 * C(4,3)) / C(3,1)
    ELSE ! N = 2
      C(4,1) = (1.0,0.0)
      C(3,1) = (1.0,0.0)
      C(2,1) = 2.0 * C(4,2)
    ENDIF

  ELSEIF ( IBCBEG.EQ.1 ) THEN
    C(4,1) = (1.0,0.0)
    C(3,1) = (0.0,0.0)
    C(2,1) = C(2,1)

  ELSEIF ( IBCBEG.EQ.2 ) THEN
    C(4,1) = (2.0,0.0)
    C(3,1) = (1.0,0.0)
    C(2,1) = 3.0*C(4,2) - C(2,1)*C(3,2)/2.0

  ELSE
    ! do nothing
    C(4,1) = C(4,1)
    C(3,1) = C(3,1)
    C(2,1) = C(2,1)

  ENDIF ! IF (IBCBEG.EQ.0)

  !   * RUNNING CALCULATIONS TO N-1 - LOOP IS NOT EXECUTED IF N = 2 *

  DO M = 2,L
!$TAF store C = cspline1
    G = -C(3,M+1) / C(4,M-1)
    C(2,M) = G*C(2,M-1) + 3.0*(C(3,M)*C(4,M+1) + C(3,M+1)*C(4,M))
    C(4,M) = G*C(3,M-1) + 2.0*(C(3,M) + C(3,M+1))
  ENDDO

  !   * ENDING BOUNDARY CONDITION SECTION *

  IF ( IBCEND.NE.1 ) THEN
!$TAF store C = cspline1
    IF ( IBCEND.EQ.0 ) THEN
      IF ( N.EQ.2 .AND. IBCBEG.EQ.0 ) THEN
        C(2,N) = C(4,N)
      ELSEIF ( (N.EQ.3 .AND. IBCBEG.EQ.0) .OR. N.EQ.2 ) THEN
        C(2,N) = 2.0 * C(4,N)
        C(4,N) = (1.,0.)
        G = -1.0 / C(4,N-1)
      ELSE
!$TAF store C = cspline1
        G = C(3,N-1) + C(3,N)
        C(2,N) = ((C(3,N) + 2.0*G) * C(4,N)*C(3,N-1) + &
               C(3,N)**2 * (C(1,N-1)-C(1,N-2)) / C(3,N-1)) / G
        G = -G / C(4,N-1)
        C(4,N) = C(3,N-1)
      ENDIF

    ELSEIF (IBCEND.EQ.2) THEN
!$TAF store C = cspline1
      C(2,N) = 3.0 * C(4,N) + C(2,N)*C(3,N)/2.0
      C(4,N) = (2.0,0.0)
      G = -1.0 / C(4,N-1)
    ELSE
      ! do nothing
      C(2,N) = C(2,N) 
      C(4,N) = C(4,N) 
      G = G
    ENDIF ! IF (IBCEND.EQ.0)

    IF ( IBCBEG.GT.0 .OR. N.GT.2 ) THEN
!$TAF store C = cspline1
      C(4,N) = G*C(3,N-1) + C(4,N)
      C(2,N) = (G*C(2,N-1) + C(2,N)) / C(4,N)
    ENDIF

  ENDIF ! IF (IBCEND.NE.1)

  !   * RUN THE ENDING BOUNDARY EFFECT BACK THROUGH THE EQUATIONS *

  DO J = L,1,-1
    C(2,J) = (C(2,J) - C(3,J)*C(2,J+1)) / C(4,J)
  ENDDO

  !   * FINAL CALCULATIONS *

  DO I = 2,N
    DTAU = C(3,I)
    DIVDF1 = (C(1,I)-C(1,I-1)) / DTAU
    DIVDF3 = C(2,I-1) + C(2,I) - 2.0*DIVDF1
    C(3,I-1) = 2.0 * (DIVDF1-C(2,I-1)-DIVDF3) / DTAU
    C(4,I-1) = (DIVDF3/DTAU) * (6.0/DTAU)
  ENDDO

  !     * ADD THE CURVATURE AT THE LAST POINT IN THE THIRD POSITION OF THE
  !       LAST NODE *

  C(3,N) = C(3,L) + (TAU(N)-TAU(L)) * C(4,L)


  !     * ADD THE MEAN VALUE OF THE ENTIRE INTERVAL IN THE FOURTH POSITION OF
  !       THE LAST NODE.  MEAN VALUE IS CALCULATED AS THE INTEGRAL OVER THE
  !       INTERVAL DIVIDED BY THE LENGTH OF THE INTERVAL. *

  C(4,N) = (0.0,0.0)
  DO I = 1,L                            ! INTEGRATE OVER THE INTERVAL
!$TAF store C = cspline1
    DTAU = TAU(I+1) - TAU(I)
    C(4,N) = C(4,N) + DTAU*(C(1,I) + DTAU*(C(2,I)/2.0 + &
          DTAU*(C(3,I)/6.0 + DTAU*C(4,I)/24.0)))
  ENDDO

  C(4,N) = C(4,N) / (TAU(N)-TAU(1))        ! DIVIDE BY LENGTH OF INTERVAL

  RETURN
  END !SUBROUTINE CSPLINE

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: SPLINEALL
! !INTERFACE:
  SUBROUTINE SPLINEALL ( C, H, F, FX, FXX )
! !DESCRIPTION:
!  THIS ROUTINE EVALUATES THE SPLINE, SPLINE DERIVATIVE, AND
!  SPLINE 2ND DERIVATIVE AT THE POINT H.

! !USES:
  IMPLICIT REAL (KIND=_RL90) ( A-G, O-Z )

! !INPUT PARAMETERS:
! C :: matrix of coefficients, 4 rows, 4 columns
! H :: point at which to evaluate the spline and its derivatives
! F, FX, FXX :: output variables for the spline value, first derivative, and second derivative
  COMPLEX (KIND=_RL90), INTENT(IN)  :: C(4)
  REAL (KIND=_RL90),    INTENT(IN)  :: H
  COMPLEX (KIND=_RL90), INTENT(OUT) :: F, FX, FXX
! !OUTPUT PARAMETERS: F, FX, FXX

! !LOCAL VARIABLES:
! HALF, SIXTH :: constants for spline evaluation
  REAL (KIND=_RL90), PARAMETER ::  HALF = 0.5, SIXTH = 1.0 / 6.0
!EOP

  F   = C(1) + H * ( C(2) + H * ( HALF * C(3) + SIXTH * H * C(4) ) )
  FX  = C(2) + H * ( C(3) + H * HALF * C(4) )
  FXX = C(3) + H * C(4)

  RETURN
  END !SUBROUTINE SPLINEALL

END MODULE splinec_mod
