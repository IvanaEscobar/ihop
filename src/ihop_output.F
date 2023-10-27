#include "IHOP_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_OUTPUT( myTime, myIter, myThid )

C !DESCRIPTION:
C calls subroutine that calculate diagnostic specific to IHOP
C EMPTY FOR NOW

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
cccccccccccccccccccccc
c diagnostics        c
cccccccccccccccccccccc


#endif /* ALLOW_IHOP */

      RETURN
      END
