#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: poly_mod
MODULE poly_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module contains functions for polynomial interpolation. - mbp 7/2015

! !USES: 
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC Poly
!=======================================================================

! == Module variables ==
  INTEGER, PRIVATE :: i, j

! == Derived types ==
  INTERFACE Poly
     MODULE PROCEDURE PolyR, PolyC, PolyZ
  END INTERFACE Poly
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! FXN PolyR
! FXN PolyC
! FXN PolyZ
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: PolyR
! !INTERFACE:
  FUNCTION PolyR( x0, x, f, N )
! !DESCRIPTION:
!   This function computes the value of a polynomial at a given point using
!   Newton's divided difference method. The polynomial is defined by its
!   values at a set of points (x, f).

! !USES: None

! !INPUT PARAMETERS:
!   x0 :: The point at which the polynomial is evaluated.
!   x  :: An array of x values defining the polynomial.
!   f  :: An array of function values corresponding to the x values.
!   N  :: The order of the polynomial (number of points - 1).
  INTEGER, INTENT( IN ) :: N
  REAL,    INTENT( IN ) :: x0, x( N ), f( N )
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! ft :: Array to hold the divided differences.
! h  :: Array to hold the differences between x values and x0.
! PolyR :: The result of the polynomial evaluation at x0.
  REAL :: ft( N ), h( N ), PolyR
!EOP

  ! Initialize arrays

  h  = x - x0
  ft = f

  ! Recursion for solution
  IF ( N.GE.2 ) THEN
    DO i = 1, N - 1
      DO j = 1, N - i
        ft( j ) = ( h( j + i ) * ft( i ) - h( i ) * ft( i + 1 ) ) / &
                  ( h( j + i ) - h( j ) )
      ENDDO
    ENDDO
  ENDIF

  PolyR = ft( 1 )

  RETURN
  END FUNCTION PolyR

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: PolyC
! !INTERFACE:
  COMPLEX FUNCTION PolyC( x0, x, f, N )
! !DESCRIPTION:
!   This function computes the value of a polynomial at a given point using
!   Newton's divided difference method. The polynomial is defined by its
!   values at a set of points (x, f).

! !USES: None

! !INPUT PARAMETERS:
!   x0 :: The point at which the polynomial is evaluated.
!   x  :: An array of x values defining the polynomial.
!   f  :: An array of function values corresponding to the x values.
!   N  :: The order of the polynomial (number of points - 1).
  INTEGER, INTENT( IN ) :: N
  COMPLEX, INTENT( IN ) :: x0, x( N ), f( N )
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! ft :: Array to hold the divided differences.
! h  :: Array to hold the differences between x values and x0.
  COMPLEX :: ft( N ), h( N )
!EOP

  ! Initialize arrays
  h  = x - x0
  ft = f

  ! Recursion for solution
  IF ( N.GE.2 ) THEN
    DO i = 1, N-1
      DO j = 1, N-I
        !         ft( J ) = ( h( J+I ) * ft( J ) - h( J ) * ft( J+1 ) ) / &
        !                                   ( h( J+I ) - h( J ) )
        ft( J ) = ft( J ) + h( J ) * ( ft( J )  - ft( J+1 ) ) / &
                                     ( h( J+I ) - h(  J   ) )
      ENDDO
    ENDDO
  ENDIF

  PolyC = ft( 1 )

  RETURN
  END FUNCTION PolyC

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: PolyZ
! !INTERFACE:
  FUNCTION PolyZ( x0, x, f, N )
! !DESCRIPTION:
!   This function computes the value of a polynomial at a given point using
!   Newton's divided difference method. The polynomial is defined by its
!   values at a set of points (x, f).

! !USES: None

! !INPUT PARAMETERS:
!   x0 :: The point at which the polynomial is evaluated.
!   x  :: An array of x values defining the polynomial.
!   f  :: An array of function values corresponding to the x values.
!   N  :: The order of the polynomial (number of points - 1).
  INTEGER, INTENT( IN ) :: N
  COMPLEX (KIND=_RL90), INTENT( IN ) :: x0, x( N ), f( N )
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! ft :: Array to hold the divided differences.
! h  :: Array to hold the differences between x values and x0.
  COMPLEX (KIND=_RL90) :: ft( N ), h( N )
  COMPLEX (KIND=_RL90) :: PolyZ
!EOP

  ! Initialize arrays
  h  = x - x0
  ft = f

  ! Recursion for solution
  IF ( N.GE.2 ) THEN
    DO i = 1, N - 1
      DO j = 1, N - i
        ft( j ) = ( h( j + i ) * ft( j ) - h( j ) * ft( j + 1 ) ) &
                  / ( h( j + i ) - h( j ) )
      ENDDO
    ENDDO
  ENDIF

  PolyZ = ft( 1 )

  RETURN
  END FUNCTION PolyZ

END !MODULE poly_mod