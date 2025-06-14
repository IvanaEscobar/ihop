#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !INTERFACE:
      SUBROUTINE IHOP_DIAGNOSTICS_STATE(
     I             myTime,myIter, myThid )

C     !DESCRIPTION:
C     Provide the information for iHOP diagnostics

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef IHOP_3D_STATE
C sound speed [m/s]
       CALL DIAGNOSTICS_FILL( ihop_ssp,'SVEL    ',
     &                        0, Nr, 0, 1, 1, myThid )
#endif /* IHOP_3D_STATE */
C sonic layer depth [m]
       CALL DIAGNOSTICS_FILL( ihop_sld,'SLD     ',
     &                        0, 1, 1, 1, 1, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_IHOP */

      RETURN
      END
