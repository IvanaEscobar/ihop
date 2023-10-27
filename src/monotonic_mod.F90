#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE monotonic_mod
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

! tests whether an input vector is strictly monotonically increasing
!
! mbp February 2014

  IMPLICIT NONE
  PRIVATE
  
! public interfaces
!=======================================================================

  public monotonic

!=======================================================================
  
  INTERFACE monotonic
     MODULE PROCEDURE monotonic_sngl, monotonic_dble
  END INTERFACE monotonic

CONTAINS
  FUNCTION monotonic_sngl( x, N )
    LOGICAL               :: monotonic_sngl
    INTEGER, INTENT( IN ) :: N
    REAL (KIND=4), DIMENSION( N ), INTENT( IN ) :: x

    monotonic_sngl = .TRUE.
    IF ( N == 1 ) RETURN
    IF ( ANY( x( 2 : N ) <= x( 1 : N - 1 ) ) ) monotonic_sngl = .FALSE.   
  RETURN
  END !FUNCTION monotonic_sngl

  FUNCTION monotonic_dble( x, N )
    LOGICAL               :: monotonic_dble
    INTEGER, INTENT( IN ) :: N
    REAL (KIND=_RL90), DIMENSION( N ), INTENT( IN ) :: x

    monotonic_dble = .TRUE.
    IF ( N == 1 ) RETURN
    IF ( ANY( x( 2 : N ) <= x( 1 : N - 1 ) ) ) monotonic_dble = .FALSE.   
  RETURN
  END !FUNCTION monotonic_dble
  
END !MODULE monotonic_mod
