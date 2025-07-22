#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_TURNOFF_IO

C !INTERFACE:
      SUBROUTINE IHOP_TURNOFF_IO( seqFlag, myThid )

C     !DESCRIPTION:
C     o Turn off some of the ptracers output flags
C     *==========================================================*
C     | Used in adjoint simulation (and called after the first
C     |  forward sweep) to avoid writing output multiple times (if
C     |  recomputations and/or grdchk) with the same iter number.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
c#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C     !INPUT PARAMETERS:
C     seqFlag :: flag that indicates where this S/R is called from:
C             :: =0 called at the end of S/R COST_FINAL
C             :: =1 called at initialisation when using DIVA
C     myThid  :: My Thread Id number
      INTEGER seqFlag
      INTEGER myThid
C     !OUTPUT PARAMETERS: None

C     !LOCAL VARIABLES:
C     CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C--   only master-thread resets shared flags (in common block)
      _BEGIN_MASTER( myThid )

C--   Set output freq. to zero to avoid re-write of
C     averaged fields in reverse checkpointing loops
      IHOP_dumpFreq    = -1.

      _END_MASTER( myThid )

      RETURN
      END