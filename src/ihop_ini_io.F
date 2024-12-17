#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_INI_IO

C !INTERFACE:
      SUBROUTINE IHOP_INI_IO( myThid )

C     !DESCRIPTION:
C     Create directory for ihop and ihop_cost output if set at runtime

C     !USES:
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


C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#if (defined ALLOW_IHOP) && (defined ALLOW_COST)
# ifndef HAVE_SYSTEM
      CHARACTER*(MAX_LEN_MBUF) msgBuf
# endif
      CHARACTER*(MAX_LEN_FNAM) namBuf
      INTEGER iL, pIL

      _BEGIN_MASTER( myThid )

      IF (IHOPObs_Dir.NE.' ') THEN
        iL = ILNBLNK( IHOPObs_Dir )
C Append "/", if necessary
        IF (iL.LT.MAX_LEN_FNAM .AND. IHOPObs_Dir(iL:iL).NE.'/') THEN
          namBuf(1:iL) = IHOPObs_Dir(1:iL)
          WRITE(IHOPObs_Dir(1:iL+1),'(2A)') namBuf(1:iL),'/'
        ENDIF
# ifdef HAVE_SYSTEM
C Create directory
        iL = ILNBLNK( IHOPObs_Dir ) - 1
        WRITE(namBuf,'(2A)') ' mkdir -p ', IHOPObs_Dir(1:iL)
        pIL = 1 + ILNBLNK( namBuf )
        WRITE(standardMessageUnit,'(3A)')
     &  '==> SYSTEM CALL (from IHOP_INI_IO): >',namBuf(1:pIL),'<'
        CALL SYSTEM( namBuf(1:pIL) )
# else
        WRITE(msgBuf,'(2A)') '** WARNING ** IHOP_INI_IO: cannot ',
     &        'call mkdir -> please create IHOPObs_Dir manually'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
# endif
      ENDIF

      _END_MASTER( myThid )

#endif /* ALLOW_IHOP && ALLOW_COST */

      RETURN
      END
