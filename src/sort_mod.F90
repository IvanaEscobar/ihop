#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: sort_mod
MODULE sort_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   Module for sorting vectors of real numbers using insertion sort.
! At the Ith step, the first I-1 positions contain a sorted
! vector.  We shall insert the Ith value into its place in that
! vector shifting up to produce a new vector of length I. - mbp 1/2015

! !USES:
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC Sort
!=======================================================================

! == Module variables ==
  INTEGER, PRIVATE :: ILeft, IMiddle, IRight, I

! == Derived types ==
  INTERFACE Sort
    MODULE PROCEDURE Sort_sngl, Sort_dble, Sort_cmplx
  END INTERFACE Sort
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R Sort_sngl, Sort_dble, Sort_cmplx
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Sort_sngl
! !INTERFACE:
  SUBROUTINE Sort_sngl( x, N )
! !DESCRIPTION:
!   Sorts a vector of single precision real numbers in ascending order.

! !USES: None
! !INPUT PARAMETERS:
! x - vector to be sorted
! N - number of elements in the vector
  INTEGER, INTENT( IN ) :: N
  REAL    :: x( * ), xTemp
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES: None
!EOP

  IF ( N.EQ.1 ) RETURN

  DO I = 2, N
    xTemp = x( I )

    IF ( xTemp.LT.x( 1 ) ) THEN
      x( 2 : I ) = x( 1 : I - 1 )
      x( 1 )     = xTemp  ! goes in the first position

    ELSEIF ( xTemp.LT.x( I-1 ) ) THEN ! Binary search for its place
      IRight = I - 1
      ILeft  = 1

      DO WHILE ( IRight.GT.ILeft+1 )
        IMiddle = ( ILeft + IRight ) / 2
        IF ( xTemp.LT.x( IMiddle ) ) THEN
          IRight = IMiddle
        ELSE
          ILeft  = IMiddle
        ENDIF

      ENDDO

      ! Shift and insert
      x( IRight+1 : I ) = x( IRight : I-1 )
      x( IRight ) = xTemp

    ELSE
      ! do nothing
    ENDIF

  ENDDO

  RETURN
  END !SUBROUTINE Sort_sngl

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Sort_dble
! !INTERFACE:
  SUBROUTINE Sort_dble( x, N )
! !DESCRIPTION:
!   Sorts a vector of double precision real numbers in ascending order.

! !USES: None
! !INPUT PARAMETERS:
! x - vector to be sorted
! N - number of elements in the vector
  INTEGER, INTENT( IN ) :: N
  REAL (KIND=_RL90) :: x( * ), xTemp
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES: None
!EOP

  IF ( N.EQ.1 ) RETURN

  DO I = 2, N
    xTemp = x( I )

    IF ( xTemp.LT.x( 1 ) ) THEN
      x( 2:I ) = x( 1:I-1 )
      x( 1 )     = xTemp  ! goes in the first position

    ELSEIF ( xTemp.LT.x( I-1 ) ) THEN ! Binary search for its place
      IRight = I - 1
      ILeft  = 1

      DO WHILE ( IRight.GT.ILeft+1 )
        IMiddle = ( ILeft + IRight ) / 2
        IF ( xTemp.LT.x( IMiddle ) ) THEN
          IRight = IMiddle
        ELSE
          ILeft  = IMiddle
        ENDIF

      ENDDO

      ! Shift and insert
      x( IRight+1:I ) = x( IRight:I-1 )
      x( IRight ) = xTemp

    ELSE
     ! do nothing
    ENDIF

  ENDDO

  RETURN
  END !SUBROUTINE Sort_dble

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: Sort_cmplx
! !INTERFACE:
  SUBROUTINE Sort_cmplx( x, N )
! !DESCRIPTION:
!   Sorts a vector of complex numbers in ascending order.

! !USES: None
! !INPUT PARAMETERS:
! x - vector to be sorted
! N - number of elements in the vector
  INTEGER, INTENT( IN ) :: N
  COMPLEX (KIND=_RL90) :: x( N ), xTemp
! !OUTPUT PARAMETERS: x

! !LOCAL VARIABLES: None
!EOP

  IF ( N.EQ.1 ) RETURN

  DO I = 2, N
    xTemp = x( I )

    IF ( REAL( xTemp ).GT.REAL( x( 1 ) ) ) THEN
      x( 2:I ) = x( 1:I-1 )
      x( 1 )     = xTemp  ! goes in the first position
    !else Binary search for its place

    ELSEIF ( REAL( xTemp ).GT.REAL( x( I-1 ) ) ) THEN
      IRight = I - 1
      ILeft  = 1

      DO WHILE ( IRight.GT.ILeft+1 )
        IMiddle = ( ILeft + IRight ) / 2
        IF ( REAL( xTemp ).GT.REAL( x( IMiddle ) ) ) THEN
          IRight = IMiddle
        ELSE
          ILeft  = IMiddle
        ENDIF

      ENDDO

      x( IRight+1:I ) = x( IRight:I-1 )
      x( IRight ) = xTemp

    ELSE
      ! do nothing
    ENDIF

  ENDDO

  RETURN
  END !SUBROUTINE Sort_cmplx

END !MODULE sort_mod