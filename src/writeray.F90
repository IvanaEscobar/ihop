#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: writeRay
MODULE writeRay
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! Compress the ray data keeping every iSkip point, points near surface or
! bottom, and last point.
! Write to RAYFile.

! During an eigenray calculation, subsets of the full ray may be passed
! These have lengths Nsteps1 vs. Nsteps for the entire ray

! !USES: 
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC WriteRayOutput
!=======================================================================

! == Module variables ==
! the maximum length of the ray vector that is written out
  INTEGER, PRIVATE :: MaxPoints = 50000   
  INTEGER, PRIVATE :: is, N2, iSkip

!EOP

CONTAINS
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! S/R WriteRayOutput
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
! !ROUTINE: WriteRayOutput
! !INTERFACE:
  SUBROUTINE WriteRayOutput( funit, nSteps, col1In, col2In, tBnc, bBnc )
  ! !DESCRIPTION:
  ! Compress and write ray data keeping every iSkip point, points near
  ! surface or bottom, and last point.

  ! Assuming col2 is always depth

  ! !USES:
  USE ihop_mod, only: rad2deg, SrcDeclAngle
  USE bdry_mod, only: Bdry

! !INPUT PARAMETERS:
! funit  :: file unit to write to
! nSteps :: number of steps in the ray
! col1In :: array of ray points (e.g. range)
! col2In :: array of ray points (e.g. depth)
! tBnc   :: top bounces
! bBnc   :: bottom bounces
  INTEGER,           INTENT( IN ) :: funit, nSteps, tBnc, bBnc
  REAL (KIND=_RL90), INTENT( IN ) :: col1In(nSteps), col2In(nSteps)
! !OUTPUT PARAMETERS: None

! !LOCAL VARIABLES:
! col1, col2 :: arrays of ray points to write out
  REAL (KIND=_RL90) :: col1(nSteps+1), col2(nSteps+1)
!EOP

  ! In adjoint mode we do not write output besides on the first run
  IF (IHOP_dumpfreq.LT.0) RETURN

  ! compression: don't print reflection points
  N2    = 1
  iSkip = MAX( nSteps / MaxPoints, 1 )
  col1  = 0.0
  col2  = 0.0

  col1(1) = col1In(1)
  col2(1) = col2In(1)
  Stepping: DO is = 2, nSteps
    ! ensure that we always write ray points near bdry reflections (works
    ! only for flat bdry)
    IF ( MIN( Bdry%Bot%HS%Depth - col2In( is ),  &
              col2In( is ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
          MOD( is, iSkip ) == 0 .OR. is==nSteps ) THEN
      N2 = N2 + 1
      col1(N2) = col1In(is)
      col2(N2) = col2In(is)

!      ELSE ! IESCO22 debugging
!        col1(is) = col1In(is)
!        col2(is) = col2In(is)

    END IF
  END DO Stepping

  ! write to output file
#ifdef IHOP_WRITE_OUT
  WRITE( fUnit, '(F12.6)' ) SNGL( SrcDeclAngle )
  WRITE( fUnit, '(3I10)' ) N2, tBnc, bBnc

  DO is = 1, N2
    WRITE( fUnit, '(2G16.10)' ) col1(is), col2(is) 
  END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE ! WriteRayOutput

END MODULE ! writeRay