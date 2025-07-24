#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: pchip_mod
MODULE pchip_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! This module contains the PCHIP subroutine and related functions.
! It implements the monotone piecewise cubic Hermite interpolating polynomial (PCHIP)
! algorithm, which is a variant of monotone PCHIP.

! !USES:
  USE splinec_mod, only: cspline
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC PCHIP
!=======================================================================

! == Module variables ==
  REAL (KIND=_RL90), PRIVATE :: h, fprime_r, fprime_i
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R PCHIP
! S/R h_del
! FXN fprime_interior_Cmplx
! FXN fprime_left_end_Cmplx
! FXN fprime_right_end_Cmplx
! FXN fprime_interior
! FXN fprime_left_end
! FXN fprime_right_end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: PCHIP
! !INTERFACE:
  SUBROUTINE PCHIP( x, y, N, PolyCoef, csWork )
! !DESCRIPTION:
! This implements the monotone piecewise cubic Hermite interpolating
! polynomial (PCHIP) algorithm. This is a new variant of monotone PCHIP
! (paper submitted to JASA). Also see;
!
! F. N. Fritsch and J. Butland, "A Method for Constructing Local Monotone
! Piecewise Cubic Interpolants", SIAM Journal on Scientific and Statistical
! Computing, 5(2):300-304, (1984) https://doi.org/10.1137/0905021
!
! F. N. Fritsch and R. E. Carlson. "Monotone Piecewise Cubic Interpolation",
! SIAM Journal on Numerical Analysis, 17(2):238-246, (1980)
! https://doi.org/10.1137/0717021

! !USES: None

! !INPUT PARAMETERS:
! x :: Vector of abscissa values (length N)
! y :: Vector of ordinate values (length N)
! N :: Number of nodes (length of x and y)
! PolyCoef :: Coefficients of the standard polynomial (4 x N)
! csWork :: Temporary workspace for the cubic spline (4 x N)
  INTEGER,              INTENT( IN  )   :: N
  REAL    (KIND=_RL90), INTENT( IN  )   :: x( * )
  COMPLEX (KIND=_RL90), INTENT( IN  )   :: y( * )
  COMPLEX (KIND=_RL90), INTENT( INOUT ) :: PolyCoef( 4, * ), csWork( 4, * )
! !OUTPUT PARAMETERS: PolyCoef, csWork

! !LOCAL VARIABLES:
! ix :: Index variable
! iBCBeg, iBCEnd :: Indices for the boundary conditions
! h1, h2 :: Step sizes for the left and right neighbors
! del1, del2 :: Differences in the ordinate values
! f1, f2 :: Function values at the nodes
! f1prime, f2prime :: Derivatives at the nodes
! fprimeT :: Temporary variable for the derivative at the endpoints
  INTEGER               :: ix, iBCBeg, iBCEnd
  REAL     (KIND=_RL90) :: h1, h2
  COMPLEX  (KIND=_RL90) :: del1, del2, f1, f2, f1prime, f2prime, fprimeT
!EOP


  !  Precompute estimates of the derivatives at the nodes
  !  The vector PolyCoef(1,*) holds the ordinate values at the nodes
  !  The vector PolyCoef(2,*) holds the ordinate derivatives at the nodes

  IF ( N.EQ.2 ) THEN
    ! handle special case of two data points seperately (linear interpolation)
    PolyCoef( 1, 1 ) = y( 1 )
    PolyCoef( 2, 1 ) = ( y( 2 ) - y( 1 ) ) / ( x( 2 ) - x( 1 ) )
    PolyCoef( 3, 1 ) = 0.0D0
    PolyCoef( 4, 1 ) = 0.0D0

  ELSE
    ! general case of more than two data points
    PolyCoef( 1, 1:N ) = y( 1:N )

    ! left endpoint (non-centered 3-point difference formula)
    CALL h_del( x, y, 2, h1, h2, del1, del2 )
    fprimeT = ( ( 2.0D0 * h1 + h2 ) * del1 - h1 * del2 ) / ( h1 + h2 )
    PolyCoef( 2, 1 ) = fprime_left_end_Cmplx( del1, del2, fprimeT )

    ! right endpoint (non-centered 3-point difference formula)

    CALL h_del( x, y, N-1, h1, h2, del1, del2 )
    fprimeT = ( -h2 * del1 + ( h1 + 2.0D0 * h2 ) * del2 ) / ( h1 + h2 )
    PolyCoef( 2, N ) = fprime_right_end_Cmplx( del1, del2, fprimeT )

    ! compute coefficients of the cubic spline interpolating polynomial

    iBCBeg = 1   ! specified derivatives at the end points
    iBCEnd = 1
    csWork( 1, 1:N ) = PolyCoef( 1, 1:N )
    csWork( 2, 1 ) = PolyCoef( 2, 1 )
    csWork( 2, N ) = PolyCoef( 2, N )
    CALL CSpline( x, csWork( 1, 1 ), N, iBCBeg, iBCEnd, N )

    ! interior nodes (use derivatives from the cubic spline as initial estimate)
    DO ix = 2, N - 1
      CALL h_del( x, y, ix, h1, h2, del1, del2 )
      ! check if the derivative from the cubic spline satisfies monotonicity
      PolyCoef( 2, ix ) = fprime_interior_Cmplx( del1, del2, csWork( 2, ix ) )
    ENDDO

  !                                                               2      3
  ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
  !
    DO ix = 1, N - 1
      h  = x( ix+1 ) - x( ix )

      f1 = PolyCoef( 1, ix )
      f2 = PolyCoef( 1, ix+1 )

      f1prime = PolyCoef( 2, ix )
      f2prime = PolyCoef( 2, ix+1 )

      PolyCoef( 3, ix ) = ( 3.0D0 * ( f2 - f1 )  &
                        - h * ( 2.0D0 * f1prime + f2prime ) ) / h**2
      PolyCoef( 4, ix ) = ( h * ( f1prime + f2prime ) - 2.0D0 * ( f2 - f1 ) ) &
                        / h**3
    ENDDO

  ENDIF ! IF ( N.EQ.2 )

  RETURN
  END !SUBROUTINE PCHIP

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: h_del
! !INTERFACE:
  SUBROUTINE h_del( x, y, ix, h1, h2, del1, del2 )
! !DESCRIPTION:
! This subroutine computes the step sizes and differences in the ordinate values
! for the left and right neighbors of a given index in the input vectors x and y.

! !USES: None

! !INPUT PARAMETERS:
! x :: Vector of abscissa values (length N)
! y :: Vector of ordinate values (length N)
! ix :: Index of the center point (1-based index)
! h1, h2 :: Step sizes for the left and right neighbors
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
  REAL (KIND=_RL90),    INTENT( IN  ) :: x( * )
  COMPLEX (KIND=_RL90), INTENT( IN  ) :: y( * )
  INTEGER,              INTENT( IN  ) :: ix
  REAL (KIND=_RL90),    INTENT( OUT ) :: h1, h2
  COMPLEX (KIND=_RL90), INTENT( OUT ) :: del1, del2
! !OUTPUT PARAMETERS: h1, h2, del1, del2

! !LOCAL VARIABLES: None
!EOP

  h1   =   x( ix     ) - x( ix - 1 )
  h2   =   x( ix + 1 ) - x( ix     )

  del1 = ( y( ix     ) - y( ix - 1 ) ) / h1
  del2 = ( y( ix + 1 ) - y( ix     ) ) / h2

  RETURN
  END !SUBROUTINE h_del

 !---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_interior_Cmplx
! !INTERFACE:
  FUNCTION fprime_interior_Cmplx( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the interior derivative of the PCHIP interpolant
! for complex input values. It checks if the derivative is within the trust region
! and projects it into the trust region if not.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_interior_Cmplx

! !LOCAL VARIABLES:
! a, b, c :: Real parts of the input complex numbers
  REAL (KIND=_RL90)    :: a,b,c
  COMPLEX (KIND=_RL90) :: fprime_interior_Cmplx
!EOP

  a=REAL(del1)
  b=REAL(del2)
  c=REAL(fprime)
  fprime_r = fprime_interior( a,b,c )
  a=AIMAG(del1)
  b=AIMAG(del2)
  c=AIMAG(fprime)
  fprime_i = fprime_interior( a,b,c )

  fprime_interior_Cmplx = CMPLX( fprime_r, fprime_i, KIND=_RL90 )

  RETURN
  END !FUNCTION fprime_interior_Cmplx

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_left_end_Cmplx
! !INTERFACE:
  FUNCTION fprime_left_end_Cmplx( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the left endpoint derivative of the PCHIP interpolant
! for complex input values. It checks if the derivative is within the trust region
! and projects it into the trust region if not.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_left_end_Cmplx

! !LOCAL VARIABLES:
! a, b, c :: Real parts of the input complex numbers
  REAL (KIND=_RL90)    :: a,b,c
  COMPLEX (KIND=_RL90) :: fprime_left_end_Cmplx
!EOP

  a=REAL(del1)
  b=REAL(del2)
  c=REAL(fprime)
  fprime_r = fprime_left_end( a,b,c )
  a=AIMAG(del1)
  b=AIMAG(del2)
  c=AIMAG(fprime)
  fprime_i = fprime_left_end( a,b,c )

  fprime_left_end_Cmplx = CMPLX( fprime_r, fprime_i, KIND=_RL90 )

  RETURN
  END !FUNCTION fprime_left_end_Cmplx

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_right_end_Cmplx
! !INTERFACE:
  FUNCTION fprime_right_end_Cmplx( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the right endpoint derivative of the PCHIP interpolant
! for complex input values. It checks if the derivative is within the trust region
! and projects it into the trust region if not.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_right_end_Cmplx

! !LOCAL VARIABLES:
! a, b, c :: Real parts of the input complex numbers
  REAL (KIND=_RL90)    :: a,b,c
  COMPLEX (KIND=_RL90) :: fprime_right_end_Cmplx
!EOP

  a=REAL(del1)
  b=REAL(del2)
  c=REAL(fprime)
  fprime_r = fprime_right_end( a,b,c )
  a=AIMAG(del1)
  b=AIMAG(del2)
  c=AIMAG(fprime)
  fprime_i = fprime_right_end( a,b,c )

  fprime_right_end_Cmplx = CMPLX( fprime_r, fprime_i, KIND=_RL90 )

  RETURN
  END !FUNCTION fprime_right_end_Cmplx

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_interior
! !INTERFACE:
  FUNCTION fprime_interior( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the interior derivative of the PCHIP interpolant.
! It checks if the derivative is within the trust region and projects it into the trust region if not.
! This is a real-valued version of the function.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_interior

! !LOCAL VARIABLES:
  REAL (KIND=_RL90) :: fprime_interior
!EOP

  ! check if derivative is within the trust region, project into it if not

  IF ( del1*del2.GT.0.0 ) THEN
    ! adjacent secant slopes have the same sign, enforce monotonicity
    IF ( del1.GT.0.0 ) THEN
      fprime_interior = MIN( MAX(fprime, 0.0D0), 3.0D0 * MIN(del1, del2) )
    ELSE
      fprime_interior = MAX( MIN(fprime, 0.0D0), 3.0D0 * MAX(del1, del2) )
    ENDIF
  ELSE
    ! force the interpolant to have an extrema here
    fprime_interior = 0.0D0;
  ENDIF

  RETURN
  END !FUNCTION fprime_interior

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_left_end
! !INTERFACE:
  FUNCTION fprime_left_end( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the left endpoint derivative of the PCHIP interpolant.
! It checks if the derivative is within the trust region and projects it into the trust region if not.
! This is a real-valued version of the function.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_left_end

! !LOCAL VARIABLES:
  REAL (KIND=_RL90) :: fprime_left_end
!EOP

  fprime_left_end = fprime

  IF ( del1*fprime.LE.0.0D0 ) THEN
      ! set derivative to zero if the sign differs from sign of secant slope
      fprime_left_end = 0.0;
  ELSEIF ( ( del1*del2.LE. 0.0D0 ) .AND. &
          ( ABS( fprime ).GT.ABS( 3.0D0 * del1 ) ) ) THEN
    ! adjust derivative value to enforce monotonicity
    fprime_left_end = 3.0D0 * del1;
  ELSE
    ! do nothing
  ENDIF

  RETURN
  END !FUNCTION fprime_left_end

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: fprime_right_end
! !INTERFACE:
  FUNCTION fprime_right_end( del1, del2, fprime )
! !DESCRIPTION:
! This function computes the right endpoint derivative of the PCHIP interpolant.
! It checks if the derivative is within the trust region and projects it into the trust region if not.
! This is a real-valued version of the function.

! !USES: None

! !INPUT PARAMETERS:
! del1, del2 :: Differences in the ordinate values for the left and right neighbors
! fprime :: Derivative at the center point
  REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
! !OUTPUT PARAMETERS: fprime_right_end

! !LOCAL VARIABLES:
  REAL (KIND=_RL90) :: fprime_right_end
!EOP

  fprime_right_end = fprime

  IF ( del2*fprime.LE.0.0D0 ) THEN
    ! set derivative to zero if the sign differs from sign of secant slope
    fprime_right_end = 0.0;
  ELSEIF ( ( del1*del2.LE.0.0D0 ) .AND. &
            ( ABS( fprime ).GT.ABS( 3.0D0 * del2 ) ) ) THEN
    ! adjust derivative value to enforce monotonicity
    fprime_right_end = 3.0D0 * del2;
  ELSE
    ! do nothing
  ENDIF

  RETURN
  END !FUNCTION fprime_right_end

END !MODULE pchip_mod