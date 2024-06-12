#include "IHOP_OPTIONS.h"

      subroutine ihop_cost_init_fixed( myThid )

C     ==================================================================
C     SUBROUTINE ihop_cost_init_fixed
C     ==================================================================
C
C     o Set contributions to the cost function and the cost function
C       itself to zero. The cost function and the individual contribu-
C       tions are defined in the header file "IHOP_COST.h".
C
C     ==================================================================
C     SUBROUTINE ihop_cost_init_fixed
C     ==================================================================

      implicit none

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#ifdef ALLOW_COST
# include "IHOP_COST.h"
#endif
#ifdef ALLOW_CAL
# include "cal.h"
#endif

C     == routine arguments ==
      integer myThid

#ifdef ALLOW_COST
C     == local variables ==
      integer tempDate1(4)
      integer tempDate2(4)
C     == end of interface ==

#ifdef ALLOW_CAL

      _BEGIN_MASTER( myThid )

C      IF ( useCAL ) THEN
CC--   Ihop cost function start and stop times in model seconds.
C        call cal_FullDate  ( costIhopStart1 , costIhopStart2,
C     &                       tempDate1     , myThid )
C        call cal_TimePassed( modelStartDate, tempDate1    ,
C     &                       tempDate2     , myThid )
C        call cal_ToSeconds ( tempDate2 , costIhopStart , myThid )
C        costIhopStart = costIhopStart + startTime
C
C        if ( costIhopEnd1.EQ.0 ) then
CC-      use "modelEndDate" as default value for costIhopEnd
C          costIhopEnd1 = modelEndDate(1)
C          costIhopEnd2 = modelEndDate(2)
C        endif
C        call cal_FullDate  ( costIhopEnd1   , costIhopEnd2  ,
C     &                       tempDate1     , myThid )
C        call cal_TimePassed( modelStartDate, tempDate1    ,
C     &                       tempDate2     , myThid )
C        call cal_ToSeconds ( tempDate2 , costIhopEnd , myThid )
C        costIhopEnd = costIhopEnd + startTime
CC--
C      ENDIF
      _END_MASTER( myThid )

      _BARRIER

#endif /* ALLOW_CAL */
#endif /* ALLOW_COST */

      RETURN
      END
