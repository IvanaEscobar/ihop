#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_OUTPUT

C !INTERFACE:
      SUBROUTINE BELLI_OUTPUT( myTime, myIter, myThid )

C     !DESCRIPTION:
C     IESCO25: unused
C     calls subroutine that calculate diagnostic specific to belli

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_BELLI
# include "BELLI_SIZE.h"
# include "BELLI.h"
#endif
#ifdef ALLOW_MNC
#endif /* ALLOW_MNC */

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !FUNCTIONS:
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
      CHARACTER*(10) suff
#ifdef ALLOW_MNC
      CHARACTER*(1) pf
#endif /* ALLOW_MNC */
CEOP

      IF ( DIFFERENT_MULTIPLE(BELLI_dumpFreq,myTime,deltaTClock)
     &  .OR. dumpInitAndLast.AND.( myTime.EQ.endTime .OR.
     &                             myTime.EQ.startTime  )
     &     ) THEN
#ifdef ALLOW_BELLI
cccccccccccccccccccccc
c diagnostics        c
cccccccccccccccccccccc
# ifdef ALLOW_MNC
       IF ( BELLI_MNC ) THEN
        IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
          pf(1:1) = 'D'
        ELSE
          pf(1:1) = 'R'
        ENDIF
C        CALL MNC_CW_SET_UDIM('belli', -1, myThid)
C        CALL MNC_CW_RL_W_S('D','belli',0,0,'T', myTime, myThid)
C        CALL MNC_CW_SET_UDIM('belli', 0, myThid)
C        CALL MNC_CW_I_W_S('I','belli',0,0,'iter', myIter, myThid)
C        CALL MNC_CW_RL_W_S('D','belli',0,0,'model_time',
C     &       myTime,myThid)
C
C        CALL MNC_CW_RL_W( pf, 'belli', 0, 0,
C     &       'belli_soundSpeed', belli_ssp, myThid )

       ENDIF ! IF ( BELLI_MNC ) THEN
# endif /* ALLOW_MNC */

       IF ( BELLI_MDSIO ) THEN 
        IF ( rwSuffixType.EQ.0 ) THEN
          WRITE(suff,'(I10.10)') myIter
        ELSE
          CALL RW_GET_SUFFIX( suff, myTime, myIter, myThid )
        ENDIF

        CALL WRITE_FLD_XYZ_RL('SVEL.',suff,belli_ssp,myIter,myThid)

       ENDIF ! IF ( BELLI_MDSIO ) THEN

#endif /* ALLOW_BELLI */
      ENDIF ! DIFFERENT_MULTIPLE

      RETURN
      END
