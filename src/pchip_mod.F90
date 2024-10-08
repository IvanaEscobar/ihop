#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE pchip_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! subroutines and functions related to the calculation of the
  ! Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)

  USE splinec_mod, only: cspline

  IMPLICIT NONE
  PRIVATE

! public interfaces
!=======================================================================

  public PCHIP

!=======================================================================

  REAL (KIND=_RL90), PRIVATE    :: h, fprime_r, fprime_i

CONTAINS
  SUBROUTINE PCHIP( x, y, N, PolyCoef, csWork )

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
    !
    ! N is the number of nodes
    ! x is a vector of the abscissa values
    ! y is a vector of the associated ordinate values
    ! PolyCoef are the coefficients of the standard polynomial
    ! csWork is a temporary work space for the cubic spline

    INTEGER,              INTENT( IN  )   :: N
    REAL    (KIND=_RL90), INTENT( IN  )   :: x( * )
    COMPLEX (KIND=_RL90), INTENT( IN  )   :: y( * )
    COMPLEX (KIND=_RL90), INTENT( INOUT ) :: PolyCoef( 4, * ), csWork( 4, * )

    INTEGER               :: ix, iBCBeg, iBCEnd
    REAL     (KIND=_RL90) :: h1, h2
    COMPLEX  (KIND=_RL90) :: del1, del2, f1, f2, f1prime, f2prime, fprimeT

    !  Precompute estimates of the derivatives at the nodes
    !  The vector PolyCoef(1,*) holds the ordinate values at the nodes
    !  The vector PolyCoef(2,*) holds the ordinate derivatives at the nodes

    IF ( N .EQ. 2 ) THEN

    ! handle special case of two data points seperately (linear interpolation)

    PolyCoef( 1, 1 ) = y( 1 )
    PolyCoef( 2, 1 ) = ( y( 2 ) - y( 1 ) ) / ( x( 2 ) - x( 1 ) )
    PolyCoef( 3, 1 ) = 0.0D0
    PolyCoef( 4, 1 ) = 0.0D0

    ELSE

    ! general case of more than two data points

    PolyCoef( 1, 1 : N ) = y( 1 : N )

    ! left endpoint (non-centered 3-point difference formula)

    CALL h_del( x, y, 2, h1, h2, del1, del2 )
    fprimeT = ( ( 2.0D0 * h1 + h2 ) * del1 - h1 * del2 ) / ( h1 + h2 )
    PolyCoef( 2, 1 ) = fprime_left_end_Cmplx( del1, del2, fprimeT )

    ! right endpoint (non-centered 3-point difference formula)

    CALL h_del( x, y, N - 1, h1, h2, del1, del2 )
    fprimeT = ( -h2 * del1 + ( h1 + 2.0D0 * h2 ) * del2 ) / ( h1 + h2 )
    PolyCoef( 2, N ) = fprime_right_end_Cmplx( del1, del2, fprimeT )

    ! compute coefficients of the cubic spline interpolating polynomial

    iBCBeg = 1   ! specified derivatives at the end points
    iBCEnd = 1
    csWork( 1, 1 : N ) = PolyCoef( 1, 1 : N )
    csWork( 2,     1 ) = PolyCoef( 2, 1 )
    csWork( 2,     N ) = PolyCoef( 2, N )
    CALL CSpline( x, csWork( 1, 1 ), N, iBCBeg, iBCEnd, N )

    ! interior nodes (use derivatives from the cubic spline as initial estimate)

    DO ix = 2, N - 1
       CALL h_del( x, y, ix, h1, h2, del1, del2 )
       ! check if the derivative from the cubic spline satisfies monotonicity
       PolyCoef( 2, ix ) = fprime_interior_Cmplx( del1, del2, csWork( 2, ix ) )
    END DO

    !                                                               2      3
    ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
    !

    DO ix = 1, N - 1
       h  = x( ix + 1 ) - x( ix )

       f1 = PolyCoef( 1, ix )
       f2 = PolyCoef( 1, ix + 1 )

       f1prime = PolyCoef( 2, ix )
       f2prime = PolyCoef( 2, ix + 1 )

       PolyCoef( 3, ix ) = ( 3.0D0 * ( f2 - f1 )  &
                         - h * ( 2.0D0 * f1prime + f2prime ) ) / h**2
       PolyCoef( 4, ix ) = ( h * ( f1prime + f2prime ) - 2.0D0 * ( f2 - f1 ) ) &
                         / h**3
    END DO

    END IF

  RETURN
  END !SUBROUTINE PCHIP

  !**********************************************************************!

  SUBROUTINE h_del( x, y, ix, h1, h2, del1, del2 )

    INTEGER,              INTENT( IN  ) :: ix   ! index of the center point
    REAL (KIND=_RL90),    INTENT( IN  ) :: x( * )
    COMPLEX (KIND=_RL90), INTENT( IN  ) :: y( * )
    REAL (KIND=_RL90),    INTENT( OUT ) :: h1, h2
    COMPLEX (KIND=_RL90), INTENT( OUT ) :: del1, del2

    h1   =   x( ix     ) - x( ix - 1 )
    h2   =   x( ix + 1 ) - x( ix     )

    del1 = ( y( ix     ) - y( ix - 1 ) ) / h1
    del2 = ( y( ix + 1 ) - y( ix     ) ) / h2

  RETURN
  END !SUBROUTINE h_del

  !**********************************************************************!

  FUNCTION fprime_interior_Cmplx( del1, del2, fprime )

    COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)                  :: a,b,c
    COMPLEX (KIND=_RL90)               :: fprime_interior_Cmplx

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

  !**********************************************************************!

  FUNCTION fprime_left_end_Cmplx( del1, del2, fprime )

    COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)                  :: a,b,c
    COMPLEX (KIND=_RL90)               :: fprime_left_end_Cmplx

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

  !**********************************************************************!

  FUNCTION fprime_right_end_Cmplx( del1, del2, fprime )

    COMPLEX (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)                  :: a,b,c
    COMPLEX (KIND=_RL90)               :: fprime_right_end_Cmplx

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

 !**********************************************************************!

  FUNCTION fprime_interior( del1, del2, fprime )

    REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)               :: fprime_interior

    ! check if derivative is within the trust region, project into it if not

    IF ( del1 * del2 > 0.0 ) THEN
      ! adjacent secant slopes have the same sign, enforce monotonicity
      IF ( del1 > 0.0 ) THEN
        fprime_interior = MIN( MAX(fprime, 0.0D0), 3.0D0 * MIN(del1, del2) )
      ELSE
        fprime_interior = MAX( MIN(fprime, 0.0D0), 3.0D0 * MAX(del1, del2) )
      END IF
    ELSE
      ! force the interpolant to have an extrema here
      fprime_interior = 0.0D0;
    END IF

  RETURN
  END !FUNCTION fprime_interior

  !**********************************************************************!

  FUNCTION fprime_left_end( del1, del2, fprime )

    REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)               :: fprime_left_end

    fprime_left_end = fprime

    IF ( del1 * fprime <= 0.0D0 ) THEN
       ! set derivative to zero if the sign differs from sign of secant slope
       fprime_left_end = 0.0;
    ELSE IF ( ( del1 * del2 <= 0.0D0 ) &
              .AND. ( ABS( fprime ) > ABS( 3.0D0 * del1 ) ) ) THEN
       ! adjust derivative value to enforce monotonicity
       fprime_left_end = 3.0D0 * del1;
    ELSE
       ! do nothing
    END IF

  RETURN
  END !FUNCTION fprime_left_end

  !**********************************************************************!

  FUNCTION fprime_right_end( del1, del2, fprime )

    ! This is essentially the same as fprime_left_end( del2, del1, fprime )
    ! Written separately for clarity

    REAL (KIND=_RL90), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=_RL90)               :: fprime_right_end

    fprime_right_end = fprime

    IF ( del2 * fprime <= 0.0D0 ) THEN
       ! set derivative to zero if the sign differs from sign of secant slope
       fprime_right_end = 0.0;
    ELSE IF ( ( del1 * del2 <= 0.0D0 ) &
              .AND. ( ABS( fprime ) > ABS( 3.0D0 * del2 ) ) ) THEN
       ! adjust derivative value to enforce monotonicity
       fprime_right_end = 3.0D0 * del2;
    ELSE
       ! do nothing
    END IF

  RETURN
  END !FUNCTION fprime_right_end

END !MODULE pchip_mod
