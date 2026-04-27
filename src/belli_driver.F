#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_DRIVER

C !INTERFACE:
      SUBROUTINE BELLI_DRIVER( myTime, myIter, myThid )

C     !DESCRIPTION:
C     Calculate ray travel times, tau

C     !USES:
      USE belli, only: belli_main
      IMPLICIT NONE
C     == Global variables ===
#ifdef ALLOW_BELLI
# include "SIZE.h"
# include "GRID.h"
# include "EEPARAMS.h"
# include "EESUPPORT.h"
# include "PARAMS.h"
# include "BELLI_SIZE.h"
# include "BELLI.h"
#endif /* ALLOW_BELLI */

C     !INPUT PARAMETERS:
C     myTime   :: Current time in simulation
C     myIter   :: Current time-step number
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
C     !OUTPUT PARAMETERS:

C     !LOCAL VARIABLES:
      INTEGER t
CEOP

#ifdef ALLOW_BELLI
      DO t=1,nts
        IF ( IHOP_iter(t).GE.0 ) THEN
          IF ( IHOP_iter(t).EQ.myIter ) THEN
            CALL BELLI_MAIN( myTime, myIter, myThid )
          ENDIF
        ENDIF
      ENDDO
#endif /* ALLOW_BELLI */

      RETURN
      END
