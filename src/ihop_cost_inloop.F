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
      SUBROUTINE IHOP_COST_INLOOP( myTime, myIter, myThid )

C !DESCRIPTION:
C Computes and writes model counterparts for obsfit data  

C !FORTRAN90 USE MODULES: ===============================================
      USE arr_mod,  only: Arr, Narr, MaxNArr
      USE srPos_mod,only: Pos
      USE ihop_mod, only: rad2deg, ray2D

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
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_IHOP
C !LOCAL VARIABLES: 
      LOGICAL inTime
      INTEGER k,bi,bj,num_file
      INTEGER maxindex(Nrr, Nrd)
      _RL     ihop_modval, saved_modval
      _RL     taumodval(Nrr,Nrd)
      _RL     maxmodval(Nrr,Nrd)
      _RL     arrival_amplitude(Nrr,Nrd,MaxNArr)
!      INTEGER sample_i(NUM_INTERP_PTS_OBS)
!      INTEGER sample_j(NUM_INTERP_PTS_OBS)
!      INTEGER sample_k(NUM_INTERP_PTS_OBS)
!      _RL     sample_w(NUM_INTERP_PTS_OBS)
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
  
!! TEMPORARY: check if in iteration
!      inTime = .false.
!      DO i=1,NTS
!        IF ( (IHOP_iter(i).GE.0) .AND. (IHOP_iter(i).EQ.myIter) ) THEN
!            PRINT *, myIter
!            inTime = .true.
!        ENDIF
!      ENDDO
!
!      IF (inTime) THEN
!      DO k=1,Nrd
!        DO j=1,Nrr
!! Calc Arrival Amplitude with cylindrical spreading to receiver range
!! I know the looping is suboptimal... that was inherited... fix tbd
!            DO i=1,Narr(j,k)
!                arrival_amplitude(i,j,k) = Arr(i,j,k)%A / sqrt(pos%rr(j)) 
!     &              * exp(cmplx(rad2deg*Arr(i,j,k)%Phase))
!            ENDDO !i
!! Find index of maximum amplitude for each j,k arrival, a hand-writted
!! maxloc
!            maxmodval(j,k) = arrival_amplitude(1,j,k)
!            maxindex(j,k)  = 1
!
!            DO i=1,Narr(j,k)
!                IF (arrival_amplitude(i,j,k) > maxmodval(j,k)) THEN
!                    maxmodval(j,k) = arrival_amplitude(i,j,k)
!                    maxindex(j,k)  = i
!                ENDIF
!            ENDDO !i
!
!! save predicted travel time according to max index
!            taumodval(j,k) = real( Arr( maxindex(j,k),j,k )%delay )
!
!        ENDDO !j
!      ENDDO !k
!
!
!! TEMPORARY: check if in iteration
!      ENDIF ! inTime
!
!!      DO bj=1,nSy
!!       DO bi=1,nSx
!!
!!        DO num_file=1,NFILESMAX_IHOP
!!
!!         DO sample_num=1,NSAMPLESMAX
!!          IF (sample_num.LE.sampleNo(num_file,bi,bj)) THEN
!!
!!           IF ((sample_timeE(num_file,sample_num,bi,bj).GE.myTime)
!!     &      .AND.(sample_timeS(num_file,sample_num,bi,bj).LT. 
!!     &      (myTime+deltaTclock))) THEN
!!
!!            DO k=1,NUM_INTERP_PTS_OBS
!!             sample_i(k)= sample_interp_i(num_file,sample_num,k,bi,bj)
!!             sample_j(k)= sample_interp_j(num_file,sample_num,k,bi,bj)
!!             sample_k(k)= sample_interp_k(num_file,sample_num,k,bi,bj)
!!             sample_w(k)= sample_interp_w(num_file,sample_num,k,bi,bj)
!!            ENDDO
!!
!!CADJ STORE sample_modval  = comlev1_bibj, key=itdkey, kind=isbyte
!!CADJ STORE saved_modval  = comlev1_bibj, key=itdkey, kind=isbyte
!!CADJ STORE sample_modmask  = comlev1_bibj, key=itdkey, kind=isbyte
!!
!!            sample_modval=0
!!            saved_modval=0
!!
!!C If averaging or integrating in time, read previously 
!!C recorded value and add to it
!!            IF (obsfitOperation(num_file).gt.0) THEN
!!             CALL active_read_obs_tile(num_file,
!!     &          saved_modval,
!!     &          sample_num,.false.,optimcycle,bi,bj,myThid,
!!     &          obsfit_dummy(num_file,bi,bj))
!!            ENDIF
!!
!!            CALL obsfit_interp(sample_modval,
!!     &         sample_i,sample_j,sample_k,sample_w,
!!     &         sample_type(num_file,sample_num,bi,bj),
!!     &         num_file,myTime,bi,bj,myThid)
!!
!!            sample_modval=sample_modval*
!!     &                   obsfit_facmod(num_file)
!!
!!            IF (obsfitOperation(num_file).gt.0) THEN
!!               sample_modval = saved_modval+sample_modval*deltaTclock
!!            ENDIF
!!
!!            sample_modmask(bi,bj)=sample_modmask(bi,bj)*
!!     &                   sample_weight(num_file,sample_num,bi,bj)
!!
!!            CALL active_write_obs_tile(num_file,
!!     &         sample_modval,
!!     &         sample_num,optimcycle,bi,bj,myThid,
!!     &         obsfit_dummy(num_file,bi,bj))
!!
!!           ENDIF !if ((sample_time...
!!          ENDIF !if (sampleNo(num_file,bi,bj).NE.0) then
!!         ENDDO !do sample_num...
!!        ENDDO !do num_file=1,NFILESMAX_OBS
!!      
!!       ENDDO
!!      ENDDO

        _END_MASTER( myThid )

#endif

      RETURN
      END
