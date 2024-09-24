#include "IHOP_OPTIONS.h"
!BOP
! !INTERFACE:
MODULE writeRay
    ! <CONTACT EMAIL="ivana@utexas.edu">
    !   Ivana Escobar
    ! </CONTACT>

  ! Compress the ray data keeping every iSkip point, points near surface or
  ! bottom, and last point.
  ! Write to RAYFile.

  ! During an eigenray calculation, subsets of the full ray may be passed
  ! These have lengths Nsteps1 vs. Nsteps for the entire ray


  IMPLICIT NONE
!   == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

  PRIVATE

! public interfaces
!=======================================================================

    public WriteRay2D, WriteDel2D, WriteRayOutput

!=======================================================================

! the maximum length of the ray vector that is written out
  INTEGER, PRIVATE :: MaxPoints = 50000   
  INTEGER, PRIVATE :: is, N2, iSkip

CONTAINS
  SUBROUTINE WriteRay2D( alpha0, Nsteps1 )
    ! The 2D version is for ray traces in (r,z) coordinates

    USE ihop_mod, only: RAYFile, ray2D
    USE bdry_mod, only: Bdry

    ! == Local ==
    INTEGER,           INTENT( IN ) :: Nsteps1
    REAL (KIND=_RL90), INTENT( IN ) :: alpha0   ! take-off angle of this ray
    REAL (KIND=_RL90) :: x(Nsteps1), z(Nsteps1)

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

    ! compression: don't print reflection points
    N2    = 1
    iSkip = MAX( Nsteps1 / MaxPoints, 1 )
    x = 0.0
    z = 0.0

    x(1) = ray2D(1)%x(1)
    z(1) = ray2D(1)%x(2)
    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near bdry reflections (works
       ! only for flat bdry)
       IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  &
                 ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          x( N2 ) = ray2D( is )%x(1)
          z( N2 ) = ray2D( is )%x(2)
      ELSE
          x( is ) = ray2D( is )%x(1)
          z( is ) = ray2D( is )%x(2)
       END IF
    END DO Stepping

    ! write to RAYFile
#ifdef IHOP_WRITE_OUT
    WRITE( RAYFile, '(G16.10)') alpha0
    WRITE( RAYFile, '(3I10)' ) N2, ray2D( Nsteps1 )%NumTopBnc, &
                        ray2D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( RAYFile, '(2G16.10)' ) x( is ), z( is )
    END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE WriteRay2D

! **************************************************************************** !
  SUBROUTINE WriteDel2D( alpha0, Nsteps1 )
    ! The 2D version is for ray traces in (r,z) coordinates

    USE ihop_mod, only: DELFile, ray2D
    USE bdry_mod,  only: Bdry

    INTEGER,           INTENT( IN ) :: Nsteps1
    REAL (KIND=_RL90), INTENT( IN ) :: alpha0   ! take-off angle of this ray
    REAL (KIND=_RL90) :: x(Nsteps1), z(Nsteps1)

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

    ! compression: don't print reflection points
    N2    = 1
    iSkip = MAX( Nsteps1 / MaxPoints, 1 )
    x = 0.0
    z = 0.0

    x(1) = ray2D(1)%tau
    z(1) = ray2D(1)%x(2)
    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near bdry reflections (works
       ! only for flat bdry)
       IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  &
                 ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          x( N2 ) = ray2D( is )%tau
          z( N2 ) = ray2D( is )%x(2)
      ELSE
          x( is ) = ray2D( is )%tau
          z( is ) = ray2D( is )%x(2)
       END IF
    END DO Stepping

    ! write to delay file
    x = REAL( x )
#ifdef IHOP_WRITE_OUT
    WRITE( DELFile, '(G16.10)' ) alpha0
    WRITE( DELFile, '(3I10)' ) N2, ray2D( Nsteps1 )%NumTopBnc, &
                        ray2D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( DELFile, '(2G16.10)' ) x( is ), z( is ) 
    END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE WriteDel2D

! **************************************************************************** !
  SUBROUTINE WriteRayOutput( funit, nSteps, col1In, col2In, tBnc, bBnc )
    ! The 2D version is for ray traces in (r,z) coordinates
    ! Assuming col2 is always depth
    USE ihop_mod, only: rad2deg, SrcDeclAngle
    USE bdry_mod, only: Bdry

    INTEGER,           INTENT( IN ) :: funit, nSteps, tBnc, bBnc
    REAL (KIND=_RL90), INTENT( IN ) :: col1In(:), col2In(:)
    REAL (KIND=_RL90) :: col1(nSteps), col2(nSteps)
    REAL (KIND=_RL90) :: alpha   ! take-off angle of this ray

    ! In adjoint mode we do not write output besides on the first run
    IF (IHOP_dumpfreq.LT.0) RETURN

    ! compression: don't print reflection points
    N2    = 1
    iSkip = MAX( nSteps / MaxPoints, 1 )
    col1 = 0.0
    col2 = 0.0

    col1(1) = col1In(1)
    col2(1) = col2In(1)
    Stepping: DO is = 2, nSteps
      ! ensure that we always write ray points near bdry reflections (works
      ! only for flat bdry)
      IF ( MIN( Bdry%Bot%HS%Depth - col2In( is ),  &
                col2In( is ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
           MOD( is, iSkip ) == 0 .OR. is == nSteps ) THEN
        N2 = N2 + 1
        col1(N2) = col1In(is)
        col2(N2) = col2In(is)

      ELSE
        col1(is) = col1In(is)
        col2(is) = col2In(is)

      END IF
    END DO Stepping

    ! write to output file
#ifdef IHOP_WRITE_OUT
    WRITE( fUnit, '(G16.10)' ) SrcDeclAngle 
    WRITE( fUnit, '(3I10)' ) N2, tBnc, bBnc

    DO is = 1, N2
      WRITE( fUnit, '(2G16.10)' ) col1(is), col2(is) 
    END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE WriteRayOutput


END MODULE writeRay
