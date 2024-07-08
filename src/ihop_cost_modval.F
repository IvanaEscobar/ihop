#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: IHOP_COST_MODVAL

C !INTERFACE: ==========================================================
      SUBROUTINE IHOP_COST_MODVAL(
     O                         modval, 
     I                         num_file, ri, rj, myTime, myThid )

C !DESCRIPTION:
C Computes ihop predicted datum 

C !FORTRAN90 USE MODULES: ===============================================
      USE arr_mod,  only: Arr, Narr, MaxNArr
      USE srPos_mod,only: Pos
      USE ihop_mod, only: rad2deg, i, ray2D

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
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES: 
      INTEGER ri,rj
      INTEGER ii,j
      INTEGER k,num_file
      INTEGER maxindex
      _RL     maxarrval, arrval(MaxNArr)
      _RL     modval
      COMPLEX*8 arrival_amplitude(MaxNArr)

! ======================================================================
! apply parameter-to-observable map on predicted datum
! want to find arrival amplitude, A:
!       A = 1/sqrt(r) * Amp * exp(i*rad2deg*phase)
! then want to save the index of maximum magnitude:
!       ia = argmax(abs(A))
! Finally, will save the travel time of max magnitude
!       tau_predicted = tau(ia)

! added uses of modules above , have access to Arr, Pos, and ray2D
! ======================================================================

! are ri and rj contained in Nrr AND Nrd
      IF ( (ri.GT.Nrr).and.(rj.gt.Nrd) ) THEN
            WRITE(msgBuf,'(3A)')
     &       'IHOP_COST_MODVAL:', 
     &       ' receiver range and depth ID larger than Nrr and Nrd.'
            CALL PRINT_ERROR( msgBuf, myThid)
      ENDIF

! Calc Arrival Amplitude with cylindrical spreading to receiver range
      DO ii=1,Narr(ri,rj)
       arrival_amplitude(ii) = Arr(ii,ri,rj)%A / sqrt(pos%rr(ri)) 
     &        * exp(cmplx(0,1)*rad2deg*Arr(ii,ri,rj)%Phase)
       arrval(ii) = ABS( arrival_amplitude(ii) )

      ENDDO !ii

! Find index of maximum amplitude for each ri,rj arrival, a hand-writted
! maxloc
      maxarrval = arrval(1)
      maxindex  = 1

      DO ii=1,Narr(ri,rj)
          IF (arrval(ii) > maxarrval) THEN
              maxarrval = arrval(ii)
              maxindex  = ii
          ENDIF
      ENDDO !ii

! save predicted travel time according to max index
      modval = real( Arr( maxindex,ri,rj )%delay )

#endif

      RETURN
      END
