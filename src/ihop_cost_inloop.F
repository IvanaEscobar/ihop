#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_INLOOP

C !INTERFACE:
      SUBROUTINE IHOP_COST_INLOOP( myTime, myThid )

C     !DESCRIPTION:
C     Extract ihop model counterparts for cost function data misfit

C     !USES:
      IMPLICIT NONE
C == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PARAMS.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_CTRL
# include "OPTIMCYCLE.h"
#endif
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

C     !INPUT PARAMETERS:
C     myTime :: Current time in simulation
C     myThid :: my thread ID
      _RL     myTime
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

C    !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
      INTEGER m,k,bi,bj,num_file,varid1,err
      INTEGER maxindex(Nrr, Nrd)
      _RL     ihop_modval
# ifndef ALLOW_CTRL
      INTEGER optimcycle
# endif
# ifndef ALLOW_ECCO
      INTEGER i,j
# endif
# ifdef ALLOW_AUTODIFF_TAMC
      INTEGER act1, act2, act3, act4
      INTEGER max1, max2, max3
      INTEGER ikey, itdkey
# endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf

# ifndef ALLOW_CTRL
      optimcycle = 0
# endif
# ifdef ALLOW_AUTODIFF_TAMC
      itdkey = bi + (bj - 1)*nSx + (ikey_dynamics - 1)*nSx*nSy
# endif

# ifdef ALLOW_COST
      _BEGIN_MASTER( myThid )

      DO bj = 1, nSy
        DO bi = 1, nSx
          DO num_file = 1, NFILESMAX_IHOP
            DO m = 1, NOBSMAX_IHOP
              IF (m.LE.ObsNo_tiled(num_file,bi,bj)) THEN
                IF ( (ihopObs_time(num_file,m,bi,bj).GE.myTime) .AND.
     &               (ihopObs_time(num_file,m,bi,bj).LT.
     &                (myTime+deltaTclock)) ) THEN
#  ifdef ALLOW_AUTODIFF_TAMC
!$TAF STORE ihop_modval  = comlev1_bibj, key=itdkey, byte=isbyte
#  endif
                  ihop_modval = 0

                  CALL IHOP_COST_MODVAL( ihop_modval, num_file,
     &                 1, 1, myTime, myThid )

                  CALL ACTIVE_WRITE_IHOP_TILE( num_file,
     &                 ihop_modval, m, optimcycle, bi, bj, myThid,
     &                 ihop_dummy(num_file,bi,bj) )

                ENDIF !IF ( (ihopObs_time
              ENDIF !IF (ObsNo_tiled(num_file,bi,bj).NE.0)
            ENDDO !DO m
          ENDDO !DO num_file
        ENDDO !DO bi
      ENDDO !DO bj

      _END_MASTER( myThid )

# endif /* ALLOW_COST */
#endif /* ALLOW_IHOP */

      RETURN
      END
