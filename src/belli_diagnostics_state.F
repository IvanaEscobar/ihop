#include "BELLI_OPTIONS.h"

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
#ifdef ALLOW_BELLI
# include "BELLI_SIZE.h"
# include "BELLI.h"
#endif

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

#ifdef ALLOW_BELLI
C     !LOCAL VARIABLES:
CEOP

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef BELLI_3D_STATE
C sound speed [m/s]
       CALL DIAGNOSTICS_FILL( belli_ssp,'SVEL    ',
     &                        0, Nr, 0, 1, 1, myThid )
#endif /* BELLI_3D_STATE */
C sonic layer depth [m]
       CALL DIAGNOSTICS_FILL( belli_sld,'SLD     ',
     &                        0, 1, 1, 1, 1, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_BELLI */

      RETURN
      END