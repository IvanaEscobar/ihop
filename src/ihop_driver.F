#include "IHOP_OPTIONS.h"

CBOP
C !ROUTINE: IHOP_DRIVER

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_DRIVER( myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate ray travel times, tau

C !FORTRAN90 USE MODULES: ===============================================
      USE ihop, only: ihop_main

      IMPLICIT NONE
C !GLOBAL VARIABLES: ===============================================
#ifdef ALLOW_IHOP
!   == Global variables ==
# include "SIZE.h"
# include "GRID.h"
# include "EEPARAMS.h"
# include "EESUPPORT.h"
# include "PARAMS.h"
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif /* ALLOW_IHOP */

C !INPUT PARAMETERS: ===================================================
C  myTime   :: Current time in simulation
C  myIter   :: Current time-step number
C  myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================

C !LOCAL VARIABLES: ====================================================
      INTEGER t
      LOGICAL do_ihop_main

CEOP

#ifdef ALLOW_IHOP
      do_ihop_main=.FALSE.
      DO t=1,nts
        IF ( IHOP_iter(t).GE.0 .AND. IHOP_iter(t).EQ.myIter ) THEN
          do_ihop_main=.TRUE.
        ENDIF
      ENDDO
      IF ( do_ihop_main ) CALL IHOP_MAIN( myTime, myIter, myThid )
#endif /* ALLOW_IHOP */

      RETURN
      END
