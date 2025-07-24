#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: monotonic_mod
MODULE monotonic_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
!   This module contains functions for testing the monotonicity of arrays.

! !USES:
  IMPLICIT NONE

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC monotonic
!=======================================================================

! == Module variables ==
  INTERFACE monotonic
    MODULE PROCEDURE monotonic_sngl, monotonic_dble
  END INTERFACE monotonic
!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! FXN monotonic_sngl
! FXN monotonic_dble
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: monotonic_sngl
! !INTERFACE:
  FUNCTION monotonic_sngl( x, N )
! !DESCRIPTION:
!   This function checks if the array x is monotonically increasing.

! !USES: None

! !INPUT PARAMETERS:
! x :: array to be checked
! N :: size of the array x
  INTEGER, INTENT( IN ) :: N
  REAL (KIND=4), DIMENSION( N ), INTENT( IN ) :: x
  ! !OUTPUT PARAMETERS:
  ! monotonic_sngl :: returns .TRUE. if x is monotonically increasing,
  !                   .FALSE. otherwise
  LOGICAL :: monotonic_sngl

! !LOCAL VARIABLES: None
!EOP

  monotonic_sngl = .TRUE.
  IF ( N.EQ.1 ) RETURN
  IF ( ANY( x( 2:N ).LE.x( 1:N-1 ) ) ) monotonic_sngl = .FALSE.

  RETURN
  END !FUNCTION monotonic_sngl

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: monotonic_dble
! !INTERFACE:
  FUNCTION monotonic_dble( x, N )
! !DESCRIPTION:
!   This function checks if the array x is monotonically increasing.

! !USES: None

! !INPUT PARAMETERS:
! x :: array to be checked
! N :: size of the array x
  INTEGER, INTENT( IN ) :: N
  REAL (KIND=_RL90), DIMENSION( N ), INTENT( IN ) :: x
  ! !OUTPUT PARAMETERS:
  ! monotonic_dble :: returns .TRUE. if x is monotonically increasing,
  !                   .FALSE. otherwise
  LOGICAL :: monotonic_dble

! !LOCAL VARIABLES: None
!EOP

  monotonic_dble = .TRUE.
  IF ( N.EQ.1 ) RETURN
  IF ( ANY( x( 2:N ).LE.x( 1:N-1 ) ) ) monotonic_dble = .FALSE.
  
  RETURN
  END !FUNCTION monotonic_dble

END !MODULE monotonic_mod