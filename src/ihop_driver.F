#include "IHOP_OPTIONS.h"

CBOP
C !ROUTINE: IHOP_DRIVER

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_DRIVER( myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate ray travel times, tau

C !FORTRAN90 USE MODULES: ===============================================
      USE bellhop

      IMPLICIT NONE
C  == Global variables ==
C#include "SIZE.h"
C#include "GRID.h"
C#include "EEPARAMS.h"
C#include "PARAMS.h"
C#include "IHOP_SIZE.h"
C#include "IHOP.h"

C !INPUT PARAMETERS: ===================================================
C  myTime   :: Current time in simulation
C  myIter   :: Current time-step number
C  myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================

C !LOCAL VARIABLES: ====================================================
C  bi,bj    :: Tile indices
C  lat2d    :: latitude of grid-cell center          [rad]
C pHalf3d   :: pressure at interface between 2 levels [Pa]
C pFull3d   :: pressure at level center               [Pa]
C zHalf3d   :: height of interface between 2 levels   [m]
C zFull3d   :: height of level center                 [m]
C  t3d      :: absolute temperature                   [K]
C  q3d      :: specific humidity                    [kg/kg]
C  u3d      :: wind speed, 1rst component (X-dir)    [m/s]
C  v3d      :: wind speed, 2nd  component (Y-dir)    [m/s]
      INTEGER bi, bj
      _RL lat2d   (sNx,sNy)
      _RL pHalf3d (sNx,sNy,Nr+1)
      _RL pFull3d (sNx,sNy,Nr)
      _RL zHalf3d (sNx,sNy,Nr+1)
      _RL zFull3d (sNx,sNy,Nr)
      _RL t3d     (sNx,sNy,Nr)
      _RL q3d     (sNx,sNy,Nr)
      _RL u3d     (sNx,sNy,Nr)
      _RL v3d     (sNx,sNy,Nr)
      _RL tdt3d   (sNx,sNy,Nr)
      _RL qdt3d   (sNx,sNy,Nr)
      _RL udt3d   (sNx,sNy,Nr)
      _RL vdt3d   (sNx,sNy,Nr)
CEOP

#ifdef ALLOW_IHOP
      IF ( ihop_iter .EQ. myIter ) THEN
       CALL IHOP_INIT( myThid )
      ENDIF
#endif /* ALLOW_IHOP */

      RETURN
      END
