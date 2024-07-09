#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: IHOP_COST_INLOOP

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_COST_INLOOP( myTime, myThid )

C !DESCRIPTION:
C Extract ihop model counterparts for cost function data misfit

C !USES: ===============================================================
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PARAMS.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_CTRL
#include "optim.h"
!# include "OPTIMCYCLE.h"
#endif
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
# include "netcdf.inc"
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

C !INPUT PARAMETERS: 
C myTime    :: Current time in simulation
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES: 
      INTEGER m,k,bi,bj,num_file,varid1,err
      INTEGER maxindex(Nrr, Nrd)
      _RL     ihop_modval
#ifndef ALLOW_CTRL
      INTEGER optimcycle
#endif
#ifndef ALLOW_ECCO
      INTEGER i,j
#endif
#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER act1, act2, act3, act4
      INTEGER max1, max2, max3
      INTEGER ikey, itdkey
#endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifndef ALLOW_CTRL
      optimcycle = 0
#endif
#ifdef ALLOW_AUTODIFF_TAMC
      itdkey = bi + (bj - 1)*nSx + (ikey_dynamics - 1)*nSx*nSy
#endif

        _BEGIN_MASTER( myThid )

      DO bj=1,nSy
       DO bi=1,nSx

        DO num_file=1,NFILESMAX_IHOP
         DO m=1,NOBSMAX_IHOP
          IF (m.LE.ObsNo_tiled(num_file,bi,bj)) THEN

           IF ((ihopObs_time(num_file,m,bi,bj).GE.myTime)
     &      .AND.(ihopObs_time(num_file,m,bi,bj).LT. 
     &      (myTime+deltaTclock))) THEN

!ADJ STORE ihop_modval  = comlev1_bibj, key=itdkey, kind=isbyte

            ihop_modval=0

            CALL ihop_cost_modval(ihop_modval, num_file,
     &          1, 1, myTime, myThid)
            
!            err = NF_INQ_VARID(fidfwd_obs(num_file,bi,bj),ihop_nameequi, 
!     &              varid1)
!            err = NF_PUT_VARA_DOUBLE(fidfwd_obs(num_file,bi,bj), varid1,
!     &              m, 1, ihop_modval)
            CALL active_write_ihop_tile( num_file,
     &         ihop_modval, m, optimcycle,bi,bj,myThid,
     &         ihop_dummy(num_file,bi,bj) )

           ENDIF !if ((ihopObs_time...
          ENDIF !if (ObsNo_tiled(num_file,bi,bj).NE.0) then
         ENDDO !do m...
        ENDDO !do num_file=1,NFILESMAX_OBS
      
       ENDDO
      ENDDO

        _END_MASTER( myThid )

#endif

      RETURN
      END
