#include "IHOP_OPTIONS.h"

CBOP
C !ROUTINE: IHOP_INI_IO

C !INTERFACE: 
      SUBROUTINE IHOP_INI_IO( myThid )

C !DESCRIPTION:
C Create directory for ihop output if needed

C !USES: 
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#ifdef ALLOW_COST
# include "IHOP_COST.h"
#endif


C !INPUT PARAMETERS: 
C myThid    :: my Thread Id number
      INTEGER myThid

C !FUNCTIONS
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK
CEOP

#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
C !LOCAL VARIABLES: 
#ifndef HAVE_SYSTEM
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#endif
      CHARACTER*(MAX_LEN_FNAM) namBuf
      INTEGER iL, pIL

      _BEGIN_MASTER( myThid )

      IF ( ihopObsDir .NE. ' ' ) THEN
       iL = ILNBLNK( ihopObsDir )
C      append "/", if necessary
       IF ( iL.LT.MAX_LEN_FNAM .AND. ihopObsDir(iL:iL).NE.'/' ) THEN
        namBuf(1:iL) = ihopObsDir(1:iL)
        WRITE(ihopObsDir(1:iL+1),'(2A)') namBuf(1:iL),'/'
       ENDIF
#ifdef HAVE_SYSTEM
C      create directory
       iL = ILNBLNK( ihopObsDir ) -1
       WRITE(namBuf,'(2A)') ' mkdir -p ', ihopObsDir(1:iL)
       pIL = 1 + ILNBLNK( namBuf )
       WRITE(standardMessageUnit,'(3A)')
     & '==> SYSTEM CALL (from IHOP_INI_IO): >',namBuf(1:pIL),'<'
       CALL SYSTEM( namBuf(1:pIL) )
#else
       WRITE(msgBuf,'(2A)') '** WARNING ** IHOP_INI_IO: cannot ',
     &       'call mkdir -> please create ihopObsDir manually'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
#endif
      ENDIF

      _END_MASTER( myThid )

#endif /* ALLOW_IHOP && ALLOW_COST */

      RETURN
      END
