#include "IHOP_OPTIONS.h"
#ifdef ALLOW_COST
#include "COST_OPTIONS.h"
#endif
#ifdef ALLOW_CTRL
#include "CTRL_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: IHOP_COST_FINAL
C     !INTERFACE:
      SUBROUTINE IHOP_COST_FINAL( myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE IHOP_COST_FINAL
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"
#include "IHOP_COST.h"
#ifdef ALLOW_COST
# include "cost.h"
#endif
#ifdef ALLOW_CTRL
# include "ctrl.h"
! c68s
# include "optim.h"
! >=c68u
!# include "optimcycle.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
      INTEGER myThid

#ifdef ALLOW_COST
C     ! FUNCTIONS:
      LOGICAL  MASTER_CPU_THREAD
      EXTERNAL MASTER_CPU_THREAD

C     !LOCAL VARIABLES:
      INTEGER bi, bj
      INTEGER ifc
      _RL f_ihop
      _RL no_ihop
      CHARACTER*23 cfname
c     CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C--   open a local file to write fc
      ifc = 30

      f_ihop      = 0. _d 0
      no_ihop     = 0. _d 0

C--   Sum up all contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

          tile_fc(bi,bj) = tile_fc(bi,bj)
     &          + mult_ihop(bi,bj) * objf_ihop(bi,bj)

       ENDDO
      ENDDO

C--   Note: global summation (tile_fc --> fc) is done only in cost_final

C--   Do global summation for each part of the cost function

      CALL GLOBAL_SUM_TILE_RL( objf_ihop,     f_ihop,     myThid )
      CALL GLOBAL_SUM_TILE_RL( num_ihop,     no_ihop,     myThid )

      WRITE(standardMessageUnit,'(A,D22.15)')
     &     ' --> f_ihop     =',f_ihop

C--   Each process has calculated the global part for itself.
      IF ( MASTER_CPU_THREAD(myThid) ) THEN

#ifdef ALLOW_CTRL
        WRITE(cfname,'(A,i4.4)') 'costfunction_ihop', optimcycle
#endif
        OPEN(unit=ifc,file=cfname)

        WRITE(ifc,*) 'fc =', fc
        WRITE(ifc,*) 'f_ihop   =', f_ihop, no_ihop

        CLOSE(ifc)

      ENDIF

#endif /* ALLOW_COST */

      RETURN
      END
