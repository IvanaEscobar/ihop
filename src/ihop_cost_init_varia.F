#include "IHOP_OPTIONS.h"

      SUBROUTINE IHOP_COST_INIT_VARIA( myThid )

C     ==================================================================
C     SUBROUTINE IHOP_COST_INIT_VARIA
C     ==================================================================
C     o Initialise the variable cost function part.
C     added ihop term: ivana@oden.utexas.edu May-2024
C
C     ==================================================================

      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#ifdef ALLOW_COST
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif

C     == routine arguments ==
      INTEGER myThid

#ifdef ALLOW_COST
C     == local variables ==
      INTEGER bi,bj
C     == end of interface ==

C--   Initialize the tiled cost function contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          objf_ihop(bi,bj)        = 0. _d 0
          num_ihop(bi,bj)         = 0. _d 0
        ENDDO
      ENDDO

      _BARRIER
#endif /* ALLOW_COST */

      RETURN
      END
