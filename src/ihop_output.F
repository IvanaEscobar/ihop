#include "IHOP_OPTIONS.h"

C !INTERFACE:
      SUBROUTINE IHOP_OUTPUT( myTime, myIter, myThid )

C !DESCRIPTION:
C calls subroutine that calculate diagnostic specific to IHOP

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
#endif
#ifdef ALLOW_MNC
#endif /* ALLOW_MNC */

C !INPUT PARAMETERS:
C  myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !FUNCTIONS:
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C Local Variables
      CHARACTER*(10) suff
#ifdef ALLOW_MNC
      CHARACTER*(1) pf
#endif /* ALLOW_MNC */
CEOP

      IF ( DIFFERENT_MULTIPLE(IHOP_dumpFreq,myTime,deltaTClock)
     &  .OR. dumpInitAndLast.AND.( myTime.EQ.endTime .OR.
     &                             myTime.EQ.startTime  )
     &     ) THEN
#ifdef ALLOW_IHOP
cccccccccccccccccccccc
c diagnostics        c
cccccccccccccccccccccc
# ifdef ALLOW_MNC
       IF ( IHOP_MNC ) THEN
        IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
          pf(1:1) = 'D'
        ELSE
          pf(1:1) = 'R'
        ENDIF
        CALL MNC_CW_SET_UDIM('ihop', -1, myThid)
        CALL MNC_CW_RL_W_S('D','ihop',0,0,'T', myTime, myThid)
        CALL MNC_CW_SET_UDIM('ihop', 0, myThid)
        CALL MNC_CW_I_W_S('I','ihop',0,0,'iter', myIter, myThid)
        CALL MNC_CW_RL_W_S('D','ihop',0,0,'model_time',
     &       myTime,myThid)

        CALL MNC_CW_RL_W( pf, 'ihop', 0, 0,
     &       'ihop_soundSpeed', ihop_ssp, myThid )

       ENDIF ! IF ( IHOP_MNC ) THEN
# endif /* ALLOW_MNC */
       IF ( IHOP_MDSIO ) THEN 
        IF ( rwSuffixType.EQ.0 ) THEN
          WRITE(suff,'(I10.10)') myIter
        ELSE
          CALL RW_GET_SUFFIX( suff, myTime, myIter, myThid )
        ENDIF

        CALL WRITE_FLD_XYZ_RL('ihop_ssp.',suff,ihop_ssp,myIter,myThid)

       ENDIF ! IF ( IHOP_MDSIO ) THEN

#endif /* ALLOW_IHOP */
      ENDIF ! DIFFERENT_MULTIPLE

      RETURN
      END
