#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_INI_IO

C !INTERFACE:
      SUBROUTINE BELLI_INI_IO( myThid )

C     !DESCRIPTION:
C     Create directory for belli and belli_cost output if set at runtime

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "BELLI_SIZE.h"
#include "BELLI.h"
#ifdef ALLOW_COST
# include "BELLI_COST.h"
#endif


C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:

C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#if (defined ALLOW_BELLI) && (defined ALLOW_COST)
# ifndef HAVE_SYSTEM
      CHARACTER*(MAX_LEN_MBUF) msgBuf
# endif
      CHARACTER*(MAX_LEN_FNAM) namBuf
      INTEGER iL, pIL
CEOP

      _BEGIN_MASTER( myThid )

      IF (B_Obs_Dir.NE.' ') THEN
        iL = ILNBLNK( B_Obs_Dir )
C Append "/", if necessary
        IF (iL.LT.MAX_LEN_FNAM .AND. B_Obs_Dir(iL:iL).NE.'/') THEN
          namBuf(1:iL) = B_Obs_Dir(1:iL)
          WRITE(B_Obs_Dir(1:iL+1),'(2A)') namBuf(1:iL),'/'
        ENDIF
# ifdef HAVE_SYSTEM
C Create directory
        iL = ILNBLNK( B_Obs_Dir ) - 1
        WRITE(namBuf,'(2A)') ' mkdir -p ', B_Obs_Dir(1:iL)
        pIL = 1 + ILNBLNK( namBuf )
        WRITE(standardMessageUnit,'(3A)')
     &    '==> SYSTEM CALL (from BELLI_INI_IO): >',namBuf(1:pIL),'<'
        CALL SYSTEM( namBuf(1:pIL) )
# else
        WRITE(msgBuf,'(2A)') '** WARNING ** BELLI_INI_IO: cannot ',
     &    'call mkdir -> please create B_Obs_Dir manually'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
# endif
      ENDIF

      _END_MASTER( myThid )

#endif /* ALLOW_BELLI && ALLOW_COST */

      RETURN
      END
