#include "AUTODIFF_OPTIONS.h"

C--  File global_vec_sum_ad.F: Perform adjoint of global sum on a 
C--                            tiled-array of vectors.
C--   Contents
C--   o GLOBAL_VEC_ADSUM_R4
C--   o GLOBAL_VEC_ADSUM_R8
C--   File global_sum_vector_ad.F: Perform adjoint of global sum
C                               on a tile array, vector field
C--    Contents
C--    o GLB_SUM_ADVEC
C--    o GLOBAL_SUM_ADVECTOR_RL
C--    o GLOBAL_SUM_ADVECTOR_RS
C--    o GLOBAL_SUM_ADVEC_ALT_RL
C--    o GLOBAL_SUM_ADVEC_ALT_RS

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLB_SUM_ADVEC

C     !INTERFACE:
      SUBROUTINE GLB_SUM_ADVEC(
     I                        nVec, adtiledVec,
     O                        adsumVec,
     I                        myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLB\_SUM\_ADVEC
C     | o Handle global sum from master thread
C     *==========================================================*
C     | Apply sum on a shared array of tiled vector
C     |  and operate over all tiles & all the processes.
C     | Warning: Only works if argument array "sumVec" is shared 
C     |          by all threads.
C     *==========================================================*
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     sumPhi :: Result of sum.
C     myThid :: My thread id.
      INTEGER nVec
      Real*8  adtiledVec(nSx,nSy,nVec)
      Real*8  adsumVec(nVec)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:

      CALL GLB_SUM_VEC( nVec, adtiledVec, adsumVec, myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_SUM_ADVECTOR_RL

C     !INTERFACE:
      SUBROUTINE GLOBAL_SUM_ADVECTOR_RL(
     I                       nVec, adtiledVec,
     O                       adsumVec,
     I                       myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL\_SUM\_ADVECTOR\_RL
C     | o Handle sum for _RL vector.
C     *==========================================================*
C     | Apply sum on an array of tiled vector (input: nSx,nSy,nVec)
C     |  and operate over all tiles & all the processes.
C     | Oversize vector is cut and processed in "nCuts" pieces
C     |  that fit in buffer
C     | Warning: Only works if argument array "sumVec" is shared 
C     |          by all threads.
C     *==========================================================*
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"

C     !INPUT/OUTPUT PARAMETERS:
C     nVec     :: length of input/output vector
C     tiledVec :: Input tiled vector (dims: nSx,nSy,nVec)
C     sumVec   :: Result of sum over all tiles, threads & procs
C     myThid   :: My thread id.
      INTEGER nVec
      _RL     adtiledVec(nSx,nSy,nVec)
      _RL     adsumVec(nVec)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==

      CALL GLOBAL_SUM_VECTOR_RL( nVec, adtiledVec, adsumVec, myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_SUM_ADVECTOR_RS

C     !INTERFACE:
      SUBROUTINE GLOBAL_SUM_ADVECTOR_RS(
     I                       nVec, adtiledVec,
     O                       adsumVec,
     I                       myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL\_SUM\_ADVECTOR\_RS
C     | o Handle sum for _RS vector.
C     *==========================================================*
C     | Apply sum on an array of tiled vector (input: nSx,nSy,nVec)
C     |  and operate over all tiles & all the processes.
C     | Oversize vector is cut and processed in "nCuts" pieces
C     |  that fit in buffer
C     | Warning: Only works if argument array "sumVec" is shared 
C     |          by all threads.
C     *==========================================================*
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"

C     !INPUT/OUTPUT PARAMETERS:
C     nVec     :: length of input/output vector
C     tiledVec :: Input tiled vector (dims: nSx,nSy,nVec)
C     sumVec   :: Result of sum over all tiles, threads & procs
C     myThid   :: My thread id.
      INTEGER nVec
      _RL     adtiledVec(nSx,nSy,nVec)
      _RL     adsumVec(nVec)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==

      CALL GLOBAL_SUM_VECTOR_RS( nVec, adtiledVec, adsumVec, myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_SUM_ADVEC_ALT_RL

C     !INTERFACE:
      SUBROUTINE GLOBAL_SUM_ADVEC_ALT_RL(
     I                       nDim, nVec, advecTile,
     O                       advecSum,
     I                       myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL\_SUM\_ADVEC_ALT\_RL
C     | o Handle sum for _RL vector.
C     *==========================================================*
C     | Apply sum on an array of tiled vector (input: nSx,nSy,nVec)
C     |  and operate over all tiles & all the processes.
C     | Oversize vector is cut and processed in "nCuts" pieces
C     |  that fit in buffer
C     | Warning: Only works if argument array "sumVec" is shared 
C     |          by all threads.
C     *==========================================================*
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"

C     !INPUT/OUTPUT PARAMETERS:
C     nDim     :: dimension of input vector
C     nVec     :: length of input/output vector to process
C     vecTile  :: Input tiled vector (dims: nDim,nSx,nSy)
C     vecSum   :: Result of sum over all tiles, threads & procs
C     myThid   :: My thread id.
      INTEGER nDim, nVec
      _RL     advecTile(nDim,nSx,nSy)
      _RL     advecSum(nVec)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==

      CALL GLOBAL_SUM_VEC_ALT_RL( nDim, nVec, advecTile, advecSum, 
     &                          myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_SUM_ADVEC_ALT_RS

C     !INTERFACE:
      SUBROUTINE GLOBAL_SUM_ADVEC_ALT_RS(
     I                       nDim, nVec, advecTile,
     O                       advecSum,
     I                       myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL\_SUM\_ADVEC_ALT\_RS
C     | o Handle sum for _RS vector.
C     *==========================================================*
C     | Apply sum on an array of tiled vector (input: nSx,nSy,nVec)
C     |  and operate over all tiles & all the processes.
C     | Oversize vector is cut and processed in "nCuts" pieces
C     |  that fit in buffer
C     | Warning: Only works if argument array "sumVec" is shared 
C     |          by all threads.
C     *==========================================================*
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"

C     !INPUT/OUTPUT PARAMETERS:
C     nDim     :: dimension of input vector
C     nVec     :: length of input/output vector to process
C     vecTile  :: Input tiled vector (dims: nDim,nSx,nSy)
C     vecSum   :: Result of sum over all tiles, threads & procs
C     myThid   :: My thread id.
      INTEGER nDim, nVec
      _RL     advecTile(nDim,nSx,nSy)
      _RL     advecSum(nVec)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==

      CALL GLOBAL_SUM_VEC_ALT_RS( nDim, nVec, advecTile, advecSum, 
     &                          myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
