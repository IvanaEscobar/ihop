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
     O                         modval,
     I                         num_file, ri, rj, myTime, myThid )

C     !DESCRIPTION:
C     Computes ihop predicted datum

C     !USES:
      USE arr_mod,  only: Arr, Narr, MaxNArr
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
      _RL     modval
CEOP

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER ii,j
      INTEGER k
      INTEGER maxindex
      _RL     maxarrval, arrval(MaxNArr)
      COMPLEX*8 arrival_amplitude(MaxNArr)

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

! Are ri and rj contained in Nrr AND Nrd
      IF ( (ri.GT.Nrr).AND.(rj.gt.Nrd) ) THEN
        WRITE(msgBuf,'(3A)')
     &   'IHOP_COST_MODVAL:',
     &   ' receiver range and depth ID larger than Nrr and Nrd.'
        CALL PRINT_ERROR( msgBuf, myThid )
      ENDIF

! Calc arrival amplitude with cylindrical spreading to receiver range
      DO ii = 1, Narr(ri,rj)
        arrival_amplitude(ii) = Arr(ii,ri,rj)%A / SQRT( Pos%Rr(ri) )
     &   * EXP( CMPLX(0,1)*rad2deg*Arr(ii,ri,rj)%Phase )
        arrval(ii) = ABS( arrival_amplitude(ii) )

      ENDDO !DO ii

!$TAF init ihop_modval = 'ihop_cost_modval'
! Find index of max amplitude for each ri,rj arrival, custom maxloc 
      maxarrval = arrval(1)
      maxindex  = 1

      DO ii = 1, Narr(ri,rj)
!$TAF store maxarrval = ihop_modval
        IF (arrval(ii)>maxarrval) THEN
          maxarrval = arrval(ii)
          maxindex  = ii
        ENDIF
      ENDDO !DO ii

! Save predicted travel time according to max index
      modval = REAL( Arr( maxindex,ri,rj )%delay )
! ray storage of top arrival...
#endif /* ALLOW_IHOP */

      RETURN
      END
