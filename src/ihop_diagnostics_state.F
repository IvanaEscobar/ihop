#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !INTERFACE:
      SUBROUTINE IHOP_DIAGNOSTICS_STATE(
     I             myTime,myIter, myThid )

C     !DESCRIPTION:
C     Provide the information for the diagnostics

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
#include "IHOP_SIZE.h"
#include "IHOP.h"
#endif

C     !INPUT PARAMETERS:
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_IHOP
C     !LOCAL VARIABLES:

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef IHOP_3D_STATE
       CALL DIAGNOSTICS_FILL( ihop_ssp,'ihop_ssp',
     &                        0,Nr, 0, 1, 1, myThid )
#endif /* IHOP_3D_STATE */
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_IHOP */

      RETURN
      END
