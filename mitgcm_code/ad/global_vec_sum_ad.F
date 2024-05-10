#include "AUTODIFF_OPTIONS.h"

C--  File global_vec_sum_ad.F: Perform adjoint of global sum on a 
C--                            tiled-array of vectors.
C--   Contents
C--   o GLOBAL_VEC_ADSUM_R4
C--   o GLOBAL_VEC_ADSUM_R8

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_VEC_ADSUM_R4

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE GLOBAL_VEC_ADSUM_R4(
     I                            myThid,
     U                            adsumPhi,
     I                            nval, ndim
     &                          )
#else
      SUBROUTINE GLOBAL_VEC_ADSUM_R4(
     I                            ndim, nval,
     U                            adsumPhi,
     I                            myThid
     &                          )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL_VEC_ADSUM_R4
C     | o Handle vector sum for real*4 data.
C     | Warning: Only works if argument array "sumPhi" is shared 
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
      Real*4 adsumPhi
      INTEGER ndim, nval, myThid
CEOP

C     !LOCAL VARIABLES:

      CALL GLOBAL_VEC_SUM_R4( ndim, nval, adsumPhi, myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_VEC_ADSUM_R8

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE GLOBAL_VEC_ADSUM_R8(
     I                            myThid,
     U                            adsumPhi,
     I                            nval, ndim
     &                          )
#else
      SUBROUTINE GLOBAL_VEC_ADSUM_R8(
     I                            ndim, nval,
     U                            adsumPhi,
     I                            myThid
     &                          )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL_VEC_ADSUM_R8
C     | o Handle vector sum for real*8 data.
C     | Warning: Only works if argument array "sumPhi" is shared 
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
      Real*8 adsumPhi
      INTEGER ndim, nval, myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==

      CALL GLOBAL_VEC_SUM_R8( ndim, nval, adsumPhi, myThid )

      RETURN
      END
