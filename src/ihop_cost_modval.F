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
     I                         num_file, ri, rj, myIter, myThid )

C     !DESCRIPTION:
C     Computes ihop predicted datum

C     !USES:
#ifdef ALLOW_IHOP
      USE arr_mod,   only: Arr, nArrival, nMaxArr
      USE srPos_mod, only: Pos
      USE ihop_mod,  only: rad2deg
#endif
#ifdef TEST_IHOP_COST_INLOOP
C      USE ssp_mod,  only: SSP, Grid
C      USE ihop_mod, only: ray2D, nMax, Beam
#endif
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
# include "IHOP.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
# endif
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

C     !INPUT PARAMETERS:
C     myIter :: Current time in simulation
C     myThid :: my thread ID
      INTEGER num_file
      INTEGER ri,rj
      _RL     myIter
      INTEGER myThid
C     !OUTPUT PARAMETERS:
      _RL     modVal

C     !LOCAL VARIABLES:
# ifdef TEST_IHOP_COST_INLOOP
      INTEGER r,k
      INTEGER bi,bj, i,j
      INTEGER itlo,ithi
      INTEGER jtlo,jthi
      _RL locQoI, locVol
      _RL locTHETA
# endif
#ifdef ALLOW_IHOP
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER ii
      INTEGER maxID
      _RL     maxarrVal, arrVal(nMaxArr)
      COMPLEX(KIND=_RL90) arrAmplitude(nMaxArr)
CEOP

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
# ifdef TEST_IHOP_COST_INLOOP
C-- OPTION 0)
C-- global summation of sound speed field within a prescribed box of
C   bcg.hd
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C -- calc modVal for timesteps with signal transmissions
      locVol=0.
      locQoI = 0.
      locTHETA = 0.

      DO bj = jtlo,jthi
      DO bi = itlo,ithi
C-- top 85 layers
C        DO k=1,85 
        DO k=1,5 
        DO j=1,sNy
        DO i=1,sNx
C-- in bounding box [12-16degE] [44-52degN]
          IF (yC(i,j,bi,bj).GE.44 .AND. yC(i,j,bi,bj).LE.52 .AND.
     &        xC(i,j,bi,bj).GE.12 .AND. xC(i,j,bi,bj).LE.16 ) THEN
      WRITE(standardMessageUnit,'(A,5I4,2F6.2)')
     & "ESCOBAR modval hit ", i,j,k, bi,bj,
     & xC(i,j, bi,bj), yC(i,j,bi,bj)

C-- sum for total volume of 3D cell 
            locVol = rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
      
C-- sum tile-wise QoI
            locQoI   = locQoI + ihop_ssp(i,j,k,bi,bj)
            locTHETA = locTHETA + theta(i,j,k,bi,bj)

            locQoI = locQoI / locVol
            locTHETA = locTHETA / locVol
          ENDIF

        ENDDO ! DO i
        ENDDO ! DO j
        ENDDO ! DO k

      ENDDO ! DO bi
      ENDDO ! DO bj

C-- report values
      modVal   = locQoI / deltaT
      locTHETA = locTHETA / deltaT
      WRITE(standardMessageUnit,'(A,E22.17)')
     & "ESCOBAR THETA/vdt ", locTHETA
C-- Sound Speed QoI END


CC-- OPTION 2b)
CC-- global summation: F90 Sound Speed if cMat is NOT default -99
C       IF ( SSP%cMat(1,3).NE.-99. _d 0 ) THEN
C       DO r = 1,Grid%nR
C       DO k = 1,Grid%nZ
CC- square of the sums
C         modVal = modVal
C     &    + SSP%cMat(k,r)/100.
C
C       ENDDO !DO j
C       ENDDO !DO k
C       ENDIF
CC-- global summation: F90 Sound Speed END

CC-- OPTION 3)
CC-- global ray2D%tau: F90 travel times if tau is NOT default 0
C      r = Beam%nSteps
C      modval = modval + real(ray2d(r)%tau)
CC-- global ray2D%tau: F90 travel times END 

CC-- OPTION 4)
CC-- global ray2D%q(1): F90 dynmaic ray parameter
C      r = Beam%nSteps
C      modval = modval + ray2d(r)%q(1)
C      PRINT *, "ESCOBAR modval", modval
CC      PRINT *, "ESCOBAR q1", ray2d(r)%q(1)
CC-- global ray2D%q(1): F90 dynamic ray parameter END 

CC-- OPTION 5)
CC-- global geninfluence: F90 eigenray travel times
C      DO ii = 1,35604
C        modval = modval + geninfluence(ii)
C        modval = modval + REAL(Arr(2,1,1)%delay)
C      ENDDO
CC-- global geninfluence: F90 eigenray travel times END

CC-- OPTION 6)
CC-- global Arr%delay: F90 eigenray travel time
C      ii = 4
C      modval = modval + REAL(Arr(ii,ri,rj)%delay)
CC-- global Arr%delay: F90 eigenray travel time END

# else /* TEST_IHOP_COST_INLOOP */

! Are ri and rj contained in nRR AND nRD
      IF ( (ri.GT.nRR).AND.(rj.gt.nRD) ) THEN
        WRITE(msgBuf,'(2A)')
     &    'IHOP_COST_MODVAL:',
     &    ' receiver range and depth ID larger than nRR and nRD.'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP "ABNORMAL END: S/R IHOP_COST_MODVAL"
      ENDIF

! Calc arrival amplitude with cylindrical spreading to receiver range
      DO ii = 1, nArrival(ri,rj)
        IF ( Pos%RR(ri).GT.0) THEN
          arrAmplitude(ii) = Arr(ii,ri,rj)%A / SQRT( Pos%RR(ri) )
     &     * EXP( CMPLX(0,1)*rad2deg*Arr(ii,ri,rj)%Phase )
        ENDIF

        arrVal(ii) = ABS( arrAmplitude(ii) )

      ENDDO !DO ii

!$TAF init ihop_modVal = 'ihop_cost_modVal'
! Find index of max amplitude for each ri,rj arrival, custom maxloc 
      maxarrVal = arrVal(1)
      maxID = 1

      DO ii = 1, nArrival(ri,rj)
!$TAF store maxarrVal = ihop_modVal
        IF (arrVal(ii)>maxarrVal) THEN
          maxarrVal = arrVal(ii)
          maxID = ii
        ENDIF
      ENDDO !DO ii

! Save predicted travel time according to max index
      modVal = REAL( Arr( maxID,ri,rj )%delay )

! ray storage of top arrival...
# endif /* TEST_IHOP_COST_INLOOP */
      WRITE(standardMessageUnit,'(A,E22.17)')
     & "ESCOBAR modVal ", modVal

#endif /* ALLOW_IHOP */

      RETURN
      END
