#include "IHOP_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_INIT_VARIA

C !INTERFACE:
      SUBROUTINE IHOP_COST_INIT_VARIA( myThid )

C    !DESCRIPTION:
C    Initialise variable cost function part

C    !USES:
      IMPLICIT NONE
C == Global Variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#ifdef ALLOW_COST
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
#endif
#ifdef ALLOW_GRDCHK
# include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my thread ID
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

C     !FUNCTIONS:
#ifdef ALLOW_GRDCHK
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK
#endif /* ALLOW_GRDCHK */

C     !LOCAL VARIABLES:
#ifdef ALLOW_COST
      INTEGER bi, bj
      INTEGER l, num_file
# ifdef ALLOW_GRDCHK
      INTEGER err
# endif /* ALLOW_GRDCHK */

C Initialize tiled cost fxn contributions.
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO num_file = 1, NFILESMAX_IHOP
            ihop_dummy(num_file,bi,bj) = 0. _d 0
# ifdef ALLOW_GRDCHK
            IF ( (obsNo_tiled(num_file,bi,bj).GT.0) .AND.
     &           (ihopDoNcOutput) ) THEN
C Sync tiled .equi file so data in buffer is not lost
              err = NF_SYNC( ncidFWD(num_file,bi,bj) )
            ENDIF
# endif /* ALLOW_GRDCHK */

          ENDDO !DO num_file
        ENDDO !DO bi
      ENDDO !DO bj

C Initialize global cost fxn contributions.
      DO num_file = 1, NFILESMAX_IHOP
        num_ihop(num_file)  = 0. _d 0
        objf_ihop(num_file) = 0. _d 0
        ihop_globaldummy(num_file) = 0. _d 0

# ifdef ALLOW_GRDCHK
         IF (obsNo(num_file).GT.0) THEN
C Sync global .equi file so data in buffer is not lost
          err = NF_SYNC( ncidGLOB(num_file) )
         ENDIF

# endif /* ALLOW_GRDCHK */
      ENDDO !DO num_file

C Initialize global variables for ihop_read_obs.F
      ihop_curfile_buff = 0
      DO l = 1, 1000
        ihop_data_buff(l) = 0. _d 0
        ihop_uncert_buff(l) = 0. _d 0
      ENDDO

      _BARRIER
#endif /* ALLOW_COST */

      RETURN
      END
