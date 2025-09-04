#include "IHOP_OPTIONS.h"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!BOP
!MODULE: ihop_mod
MODULE ihop_mod
! <CONTACT EMAIL="ivana@utexas.edu">
!   Ivana Escobar
! </CONTACT>
! !DESCRIPTION:
! Defines modules used in iHOP

! !USES:
  IMPLICIT NONE
! == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

! !SCOPE: 
  PRIVATE
!=======================================================================
  PUBLIC  rad2deg, oneCMPLX, &
          PRTFile, RAYFile, DELFile, SHDFile, ARRFile, SSPFile, &
          ATIFile, BTYFile, BRCFile, TRCFile, IRCFile, SBPFile, &
          MaxN, nRz_per_range, iStep, afreq, SrcDeclAngle,      &
          Title, Beam, ray2D, ray2DPt, iSmallStepCtr, rxyz
!=======================================================================

! == Module variables ==
  ! *** fixed parameters useful for ONLY ihop ***
  REAL(KIND=_RL90),     PARAMETER :: rad2deg = 180.D0 / PI
  COMPLEX (KIND=_RL90), PARAMETER :: oneCMPLX = ( 0.0D0, 1.0D0 )

  INTEGER, PARAMETER :: PRTFile = 61, &    ! standard output file
                        RAYFile = 21, &    ! ray paths file
                        DELFile = 22, &    ! ray paths file
                        SHDFile = 25, &    ! TL calc output file
                        ARRFile = 37, &    ! Arrivals calc output file
                        SSPFile = 40, &    ! optional 2D/3D SSP file
                        ATIFile = 41, &    ! optional 2D/3D altimetry
                        BTYFile = 42, &    ! optional 2D/3D bathymetry
                        BRCFile = 38, TRCFile = 39, IRCFile = 16, &
                        SBPFile = 50, &
                        MaxN = 50000

  ! *** varying parameters for ihop ***
  INTEGER            :: iSmallStepCtr = 0
  INTEGER            :: nRz_per_range, iStep
  REAL (KIND=_RL90)  :: afreq, SrcDeclAngle, SrcAzimAngle
  CHARACTER*(80)     :: Title

! == Derived types ==
  ! *** Beam structure ***
  TYPE rxyz
    REAL (KIND=_RL90) :: r, x, y, z
  END TYPE rxyz

  TYPE BeamStructure
    INTEGER           :: NBeams, Nimage, Nsteps, iBeamWindow
    REAL (KIND=_RL90) :: deltas, epsMultiplier = 1, rLoop
    CHARACTER*(1)     :: Component ! Pressure or displacement
    CHARACTER*(4)     :: Type = 'G S '
    CHARACTER*(7)     :: RunType
    TYPE( rxyz )      :: Box
  END TYPE BeamStructure

  TYPE( BeamStructure ) :: Beam

  ! *** ray structure ***
  TYPE ray2DPt
    INTEGER                :: NumTopBnc, NumBotBnc, NumTurnPt
    REAL    (KIND=_RL90)   :: x( 2 ), t( 2 ), p( 2 ), q( 2 ), c, Amp, Phase
    COMPLEX (KIND=_RL90)   :: tau
  END TYPE ray2DPt

  TYPE( ray2DPt )      :: ray2D( MaxN )
!EOP

END !MODULE ihop_mod
