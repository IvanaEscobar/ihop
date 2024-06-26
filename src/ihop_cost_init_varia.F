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
      INTEGER l,num_file
C     == end of interface ==

C--   Initialize the tiled cost function contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO num_file=1,NFILESMAX_IHOP
          objf_ihop(num_file,bi,bj)  =0. _d 0
           ihop_dummy(num_file,bi,bj)=0. _d 0
         ENDDO
         num_ihop(bi,bj) = 0. _d 0
        ENDDO
      ENDDO

C Initialize global variables for ihop_read_obs.F
      ihop_curfile_buff=0
      DO l=1,1000
       ihop_data_buff(l)=0. _d 0
       ihop_uncert_buff(l)=0. _d 0
      ENDDO


      _BARRIER
#endif /* ALLOW_COST */

      RETURN
      END
