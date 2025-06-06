#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_MODVAL

C !INTERFACE:
      SUBROUTINE IHOP_COST_MODVAL(
     O                         modVal,
     I                         num_file, ri, rj, myTime, myThid )

C     !DESCRIPTION:
C     Computes ihop predicted datum

C     !USES:
      USE arr_mod,  only: Arr, nArr, maxNArr
      USE srPos_mod,only: Pos
      USE ihop_mod, only: rad2deg
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
      INTEGER num_file
      INTEGER ri,rj
      _RL     myTime
      INTEGER myThid
C     !OUTPUT PARAMETERS:
      _RL     modVal
CEOP

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER ii
      INTEGER maxID
      _RL     maxarrVal, arrVal(maxNArr)
      COMPLEX*8 arrAmplitude(maxNArr)

! ===================================================================
! Apply parameter-to-observable map on predicted datum
! Want to find arrival amplitude, A:
!       A = 1/sqrt(r) * Amp * exp(oneCMPLX*rad2deg*phase)
! then want to save the index of maximum magnitude:
!       ia = argmax(abs(A))
! Finally, will save the travel time of max magnitude
!       tau_predicted = tau(ia)

! Added uses of modules above, have access to Arr and Pos
! ===================================================================

! Are ri and rj contained in nRR AND nRD
      IF ( (ri.GT.nRR).AND.(rj.gt.nRD) ) THEN
        WRITE(msgBuf,'(3A)')
     &   'IHOP_COST_MODVAL:',
     &   ' receiver range and depth ID larger than nRR and nRD.'
        CALL PRINT_ERROR( msgBuf, myThid )
      ENDIF

! Calc arrival amplitude with cylindrical spreading to receiver range
      DO ii = 1, nArr(ri,rj)
        arrAmplitude(ii) = Arr(ii,ri,rj)%A / SQRT( Pos%RR(ri) )
     &   * EXP( CMPLX(0,1)*rad2deg*Arr(ii,ri,rj)%Phase )

        arrVal(ii) = ABS( arrAmplitude(ii) )

      ENDDO !DO ii

!$TAF init ihop_modVal = 'ihop_cost_modVal'
! Find index of max amplitude for each ri,rj arrival, custom maxloc 
      maxarrVal = arrVal(1)
      maxID = 1

      DO ii = 1, nArr(ri,rj)
!$TAF store maxarrVal = ihop_modVal
        IF (arrVal(ii)>maxarrVal) THEN
          maxarrVal = arrVal(ii)
          maxID = ii
        ENDIF
      ENDDO !DO ii

! Save predicted travel time according to max index
      modVal = REAL( Arr( maxID,ri,rj )%delay )

! ray storage of top arrival...

#endif /* ALLOW_IHOP */

      RETURN
      END
