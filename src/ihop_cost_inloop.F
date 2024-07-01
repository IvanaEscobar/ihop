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
# include "OPTIMCYCLE.h"
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
      INTEGER m,k,bi,bj,num_file
      INTEGER maxindex(Nrr, Nrd)
      _RL     ihop_modval, saved_modval
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
!
!            DO k=1,NUM_INTERP_PTS_OBS
!             sample_i(k)= sample_interp_i(num_file,m,k,bi,bj)
!             sample_j(k)= sample_interp_j(num_file,m,k,bi,bj)
!             sample_k(k)= sample_interp_k(num_file,m,k,bi,bj)
!             sample_w(k)= sample_interp_w(num_file,m,k,bi,bj)
!            ENDDO

CADJ STORE ihop_modval  = comlev1_bibj, key=itdkey, kind=isbyte
CADJ STORE saved_modval  = comlev1_bibj, key=itdkey, kind=isbyte
CADJ STORE ihopObs_modmask_tiled  = comlev1_bibj, key=itdkey, kind=isbyte

            ihop_modval=0
            saved_modval=0

            print *, 'escobar: ', m
            CALL ihop_cost_modval(ihop_modval, num_file,
     &          1, 1, myTime, myThid)

!            CALL obsfit_interp(ihop_modval,
!     &         sample_i,sample_j,sample_k,sample_w,
!     &         sample_type(num_file,m,bi,bj),
!     &         num_file,myTime,bi,bj,myThid)
!
            ihop_modval=ihop_modval*3
!     &                   obsfit_facmod(num_file)
!
!            IF (obsfitOperation(num_file).gt.0) THEN
!               ihop_modval = saved_modval+ihop_modval*deltaTclock
!            ENDIF
!
!            ihopObs_modmask_tiled(bi,bj)=ihopObs_modmask_tiled(bi,bj)*
!     &                   sample_weight(num_file,m,bi,bj)
!
!            CALL active_write_obs_tile(num_file,
!     &         ihop_modval,
!     &         m,optimcycle,bi,bj,myThid,
!     &         obsfit_dummy(num_file,bi,bj))
!
           ENDIF !if ((sample_time...
          ENDIF !if (ObsNo_tiled(num_file,bi,bj).NE.0) then
         ENDDO !do m...
        ENDDO !do num_file=1,NFILESMAX_OBS
      
       ENDDO
      ENDDO

        _END_MASTER( myThid )

#endif

      RETURN
      END
