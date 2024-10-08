#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE sort_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! mbp 1/2015 incorporating subroutines from decades past
  ! Does an insertion sort on a vector of real numbers

  ! At the Ith step, the first I-1 positions contain a sorted
  ! vector.  We shall insert the Ith value into its place in that
  ! vector shifting up to produce a new vector of length I.

  IMPLICIT NONE
  PRIVATE

! public interfaces
!=======================================================================

    public Sort

!=======================================================================

  INTEGER, PRIVATE :: ILeft, IMiddle, IRight, I

  INTERFACE Sort
     MODULE PROCEDURE Sort_sngl, Sort_dble, Sort_cmplx
  END INTERFACE Sort

CONTAINS
  SUBROUTINE Sort_sngl( x, N )

    INTEGER, INTENT( IN ) :: N
    REAL    :: x( * ), xTemp

    IF ( N == 1 ) RETURN

    DO I = 2, N

       xTemp = x( I )

       IF ( xTemp < x( 1 ) ) THEN
          x( 2 : I ) = x( 1 : I - 1 )
          x( 1 )     = xTemp  ! goes in the first position
       ELSE IF ( xTemp < x( I - 1 ) ) THEN ! Binary search for its place

          IRight = I - 1
          ILeft  = 1

          DO WHILE ( IRight > ILeft + 1 )
             IMiddle = ( ILeft + IRight ) / 2
             IF ( xTemp < x( IMiddle ) ) THEN
                IRight = IMiddle
             ELSE
                ILeft  = IMiddle
             ENDIF
          END DO

          ! Shift and insert
          x( IRight + 1 : I ) = x( IRight : I - 1 )
          x( IRight ) = xTemp

       ELSE
        ! do nothing
       ENDIF

    END DO

  RETURN
  END SUBROUTINE Sort_sngl

  ! ________________________________________________________________________

  SUBROUTINE Sort_dble( x, N )

    INTEGER, INTENT( IN ) :: N
    REAL (KIND=_RL90) :: x( * ), xTemp
    IF ( N == 1 ) RETURN

    DO I = 2, N
       xTemp = x( I )

       IF ( xTemp < x( 1 ) ) THEN
          x( 2 : I ) = x( 1 : I - 1 )
          x( 1 )     = xTemp  ! goes in the first position
       ELSE IF ( xTemp < x( I - 1 ) ) THEN ! Binary search for its place

          IRight = I - 1
          ILeft  = 1

          DO WHILE ( IRight > ILeft + 1 )
             IMiddle = ( ILeft + IRight ) / 2
             IF ( xTemp < x( IMiddle ) ) THEN
                IRight = IMiddle
             ELSE
                ILeft  = IMiddle
             ENDIF
          END DO

          ! Shift and insert
          x( IRight + 1 : I ) = x( IRight : I - 1 )
          x( IRight ) = xTemp

       ELSE
        ! do nothing
       ENDIF

    END DO

  RETURN
  END SUBROUTINE Sort_dble

  ! ________________________________________________________________________

  SUBROUTINE Sort_cmplx( x, N )

    ! Based on order of decreasing real part

    INTEGER, INTENT( IN ) :: N
    COMPLEX (KIND=_RL90)      :: x( N ), xTemp

    IF ( N == 1 ) RETURN

    DO I = 2, N

       xTemp = x( I )

       IF ( REAL( xTemp ) > REAL( x( 1 ) ) ) THEN
          x( 2 : I ) = x( 1 : I - 1 )
          x( 1 )     = xTemp  ! goes in the first position
       !else Binary search for its place
       ELSE IF ( REAL( xTemp ) > REAL( x( I - 1 ) ) ) THEN

          IRight = I - 1
          ILeft  = 1

          DO WHILE ( IRight > ILeft + 1 )
             IMiddle = ( ILeft + IRight ) / 2

             IF ( REAL( xTemp ) > REAL( x( IMiddle ) ) ) THEN
                IRight = IMiddle
             ELSE
                ILeft  = IMiddle
             END IF
          END DO

          x( IRight + 1 : I ) = x( IRight : I - 1 )
          x( IRight ) = xTemp

       ELSE
        ! do nothing
       END IF

    END DO

  RETURN
  END SUBROUTINE Sort_cmplx

END MODULE sort_mod
