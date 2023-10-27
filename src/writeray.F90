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

  USE ihop_mod, only: RAYFile, DELFile, ray2D
  USE ssp_mod,  only: Bdry

  IMPLICIT NONE
  PRIVATE

! public interfaces
!=======================================================================

    public WriteRay2D, WriteDel2D

!=======================================================================

  INTEGER, PRIVATE :: MaxNRayPoints = 500000   ! this is the maximum length of 
  ! the ray vector that is written out
  INTEGER, PRIVATE :: is, N2, iSkip

CONTAINS
  SUBROUTINE WriteRay2D( alpha0, Nsteps1 )

    ! The 2D version is for ray traces in (r,z) coordinates

    INTEGER,           INTENT( IN ) :: Nsteps1
    REAL (KIND=_RL90), INTENT( IN ) :: alpha0   ! take-off angle of this ray

    ! compression

    N2    = 1
    iSkip = MAX( Nsteps1 / MaxNRayPoints, 1 )

    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near bdry reflections (works 
       ! only for flat bdry)
       IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  &
                 ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          ray2D( N2 )%x = ray2D( is )%x
       END IF
    END DO Stepping

    ! write to ray file

#ifdef IHOP_WRITE_OUT
    WRITE( RAYFile, * ) alpha0
    WRITE( RAYFile, * ) N2, ray2D( Nsteps1 )%NumTopBnc, &
                        ray2D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( RAYFile, * ) ray2D( is )%x
    END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE WriteRay2D

  SUBROUTINE WriteDel2D( alpha0, Nsteps1 )

    ! The 2D version is for ray traces in (r,z) coordinates

    INTEGER,           INTENT( IN ) :: Nsteps1
    REAL (KIND=_RL90), INTENT( IN ) :: alpha0   ! take-off angle of this ray

    ! compression

    N2    = 1
    iSkip = MAX( Nsteps1 / MaxNRayPoints, 1 )

    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near bdry reflections (works 
       ! only for flat bdry)
       IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  &
                 ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          ray2D( N2 )%x = ray2D( is )%x
       END IF
    END DO Stepping

    ! write to ray file

#ifdef IHOP_WRITE_OUT
    WRITE( DELFile, * ) alpha0
    WRITE( DELFile, * ) N2, ray2D( Nsteps1 )%NumTopBnc, &
                        ray2D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( DELFile, * ) REAL(ray2D( is )%tau), ray2D( is )%x(2)
    END DO
#endif /* IHOP_WRITE_OUT */

  RETURN
  END SUBROUTINE WriteDel2D

  ! **********************************************************************!

END MODULE writeRay
