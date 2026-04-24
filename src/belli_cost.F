#include "BELLI_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: BELLI_COST

C !INTERFACE:
      SUBROUTINE BELLI_COST( myTime, myIter, myThid )

C    !DESCRIPTION:
C    Computes belli data misfit for cost function

C    !USES:
#if (defined ALLOW_BELLI) && (defined TEST_BELLI_COST)
C     For Grid%nR,nZ, SSP%cMat
      USE ssp_mod, only: Grid, SSP

#endif
      IMPLICIT NONE
C == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef ALLOW_BELLI
# include "BELLI_SIZE.h"
# include "BELLI.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
#  include "BELLI_COST.h"
# endif
#endif
#ifdef ALLOW_CTRL
! >=c68u
# include "OPTIMCYCLE.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
#endif

C     !INPUT PARAMETERS:
C     myTime :: Current time in simulation
C     myIter :: Current time-step number
C     myThid :: my thread ID
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
C     !OUTPUT PARAMETERS:
CEOP

#ifdef ALLOW_BELLI
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
C     bi,bj :: Tile indices
      INTEGER bi,bj
      INTEGER num_file, iObs
      _RL belli_modval, belli_modvaltmp, myssp
      _RL belli_data, belli_uncert, belli_weight
# ifndef ALLOW_CTRL
      INTEGER optimcycle
# endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_FNAM) B_Obs_File
      INTEGER IL, err
      INTEGER irec
      LOGICAL exst
      _RL  objf_belli_glo
      _RL  num_belli_glo
      _RL  tmpgs
# ifdef ALLOW_COST
      _RL  B_Obs_buff(NOBSMAX_BELLI)
      _RL  B_Obs_modval_glob(NOBSMAX_BELLI)
# endif
# ifdef TEST_BELLI_COST
      INTEGER ig, jg
      INTEGER itlo,ithi
      INTEGER jtlo,jthi
      INTEGER i,j,k
      _RL objf_belli_tloc(nSx,nSy)
      _RL objf_belli_gloc
# endif

# ifndef ALLOW_CTRL
      optimcycle = 0
# endif

      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') '== belli_cost: begin =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

C Initialise local storage
#ifdef ALLOW_AUTODIFF_TAMC
CADJ INIT tapelev_belli_cost = COMMON, NFILESMAX_BELLI
#endif

# ifdef ALLOW_COST
#  ifdef TEST_BELLI_COST
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      DO num_file=1,NFILESMAX_BELLI
C-- Reset local variables
       DO bj = jtlo,jthi
       DO bi = itlo,ithi
         objf_belli_tloc(bi,bj) = 0.

       ENDDO !DO bi
       ENDDO !DO bj

       objf_belli_gloc = 0.

       DO bj = jtlo,jthi
       DO bi = itlo,ithi
CC-- OPTION 1)
CC-- tile-wise summation: F77 Sound Speed
C        DO k=1,5
C        DO j=1,sNy
C        jg = myYGlobalLo-1+(bj-1)*sNy+j
C        DO i=1,sNx
C         ig = myXGlobalLo-1+(bi-1)*sNx+i
C
C         IF (jg.GE.35 .AND. jg.LE.40 .AND. 
C     &       ig.GE.20 .AND. ig.LE.25 ) THEN
CC         IF (jg.GE.25 .AND. jg.LE.45 .AND. 
CC     &       ig.GE.10 .AND. ig.LE.30 ) THEN
C           objf_belli_tloc(bi,bj) = objf_belli_tloc(bi,bj)
C     &      + belli_ssp(i,j,k,bi,bj)*belli_ssp(i,j,k,bi,bj)
C         ENDIF
C
C        ENDDO !DO i
C        ENDDO !DO j
C        ENDDO !DO k
CC-- tile-wise summation: F77 Sound Speed END

        WRITE(msgBuf,'(A,E22.17)') 'QoI sum(C) tiled: ',
     &                      objf_belli_tloc(bi,bj)
        CALL PRINT_MESSAGE( msgBuf,
     &      standardMessageUnit, SQUEEZE_RIGHT, myThid )

       ENDDO !DO bi
       ENDDO !DO bj

C-- OPTION 2a)
C-- global summation: F90 Sound Speed
       DO j = 1,Grid%nR
       DO k = 1,Grid%nZ
C- sum of the squares
         objf_belli_gloc = objf_belli_gloc 
     &    + SSP%cMat(k,j)*SSP%cMat(k,j)

       ENDDO !DO j
       ENDDO !DO k
C-- global summation: F90 Sound Speed END

Cc--   Do global summation, if using tile-wise summation
C       IF ( objf_belli_tloc(1,1).NE.0 ) THEN
C         CALL GLOBAL_SUM_TILE_RL( objf_belli_tloc,objf_belli_gloc,myThid )
C       ENDIF
C       CALL GLOBAL_SUM_TILE_RL( objf_belli_tloc,objf_belli_gloc,myThid )

       _BEGIN_MASTER( myThid )
       objf_belli(num_file) = objf_belli(num_file) + objf_belli_gloc
       _END_MASTER( myThid )

       WRITE(msgBuf,'(A,E22.17)') 'QoI sum(C2) global: ', 
     &                      objf_belli(num_file)
       CALL PRINT_MESSAGE( msgBuf,
     &      standardMessageUnit, SQUEEZE_RIGHT, myThid )

      ENDDO !DO num_file

#  else /* TEST_BELLI_COST */
      _BEGIN_MASTER( myThid )

      DO num_file = 1, NFILESMAX_BELLI

C File maintenance
        DO bj = 1, nSy
          DO bi = 1, nSx
            IF ( (ObsNo_tiled(num_file,bi,bj).GT.0) .AND.
     &           (belliDoNcOutput) ) THEN
C Sync tiled file so data is not lost in buffer
              err = NF_SYNC( ncidFWD(num_file,bi,bj) )
              CALL BELLI_NF_ERROR('cost 1',err,bi,bj,
     & myThid)

            ENDIF

          ENDDO !DO bi
        ENDDO !DO bj

        DO iObs = 1, NOBSMAX_BELLI
          B_Obs_buff(iObs) = zeroRL
        ENDDO

        DO bj = 1, nSy
          DO bi = 1, nSx
C Open tiled files and read to buffer
            DO iObs = 1, NOBSMAX_BELLI
              IF (iObs.LE.ObsNo_tiled(num_file,bi,bj)) THEN
                belli_modvaltmp = zeroRL
C IESCO25: dummy +/- of belli_ssp to get TAF TLM to work
                myssp = belli_ssp(1,1,1,bi,bj)
                belli_modvaltmp = belli_modvaltmp + myssp
                belli_modvaltmp = belli_modvaltmp - myssp

                CALL ACTIVE_READ_BELLI_TILE( 
     I              num_file,
     O              belli_modvaltmp, 
     I              iObs,.false., optimcycle, bi, bj, myThid,
     I              belli_dummy(num_file,bi,bj) )

C Save tiled model equi from buffer, see PROFILES_READVECTOR
C Combine all threads to global
                irec = B_Obs_ind_glob_tiled(num_file,iObs,bi,bj)
C IESCO25: ANOTHER dummy +/- of belli_ssp to get TAF TLM to work
                B_Obs_buff(irec)=B_Obs_buff(irec)+belli_modvaltmp
     &            + myssp
                B_Obs_buff(irec)=B_Obs_buff(irec)
     &            - myssp
              ENDIF !IF (iObs.LE.ObsNo_tiled(num_file,bi,bj))
            ENDDO !DO iObs

          ENDDO !DO bj
        ENDDO !DO bi

C Combine MPI processes
        DO iObs = 1, NOBSMAX_BELLI
          tmpgs = B_Obs_buff(iObs)
          _GLOBAL_SUM_RL(tmpgs, myThid)
          B_Obs_modval_glob(iObs) = tmpgs

        ENDDO !DO iObs

#ifdef ALLOW_AUTODIFF_TAMC
C     Not really necessary, but cheap (for small NFILESMAX_BELLI)
C     and avoids some hidden recomputations.
CADJ STORE B_Obs_modval_glob
CADJ &     = tapelev_belli_cost, key = num_file, kind = isbyte
#endif

        IF (myProcId.EQ.0) THEN
C Loop over global obs
          DO iObs = 1, NOBSMAX_BELLI
CC TAF active var init
C            belli_globaldummy(num_file) = 0
C
            IF (iObs.LE.ObsNo(num_file)) THEN
              belli_modval = zeroRL
C ESCO Future dev: add time and space averaging here ... see pkg/obsfit

C Extract model equi
              belli_modval = belli_modval+B_Obs_modval_glob(iObs)

C Write to global netcdf file
              CALL ACTIVE_WRITE_BELLI_GLOB( 
     I             num_file, belli_modval,
     I             iObs, optimcycle, myThid,
     I             belli_globaldummy(num_file) )

            ENDIF !IF (iObs.LE.ObsNo(num_file))

          ENDDO !DO iObs

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE belli_globaldummy(num_file)
CADJ &     = tapelev_belli_cost, key = num_file, kind = isbyte
#endif

C Sync global obs
          err = NF_SYNC( ncidGLOB(num_file) )
          CALL BELLI_NF_ERROR(
     &         'COST: NF_SYNC ncidGLOB',err,0,0,myThid )

          DO iObs = 1, NOBSMAX_BELLI
            IF (iObs.LE.ObsNo(num_file)) THEN
              belli_modval = zeroRL
              belli_data   = zeroRL
              belli_uncert = zeroRL

C Read model equivalent from global file
              CALL ACTIVE_READ_BELLI_GLOB( 
     I             num_file,
     O             belli_modval, 
     I             iObs,.false., optimcycle, myThid,
     I             belli_globaldummy(num_file) )

#ifdef TEST_BELLI_COST_INLOOP
C Set uncertainty to one
              belli_uncert = 1. _d 0

#else /* TEST_BELLI_COST_INLOOP */
C Read observation and uncertainty
              CALL BELLI_COST_READ_OBS(
     I             num_file,  1, B_Obs_ind_glob(num_file,iObs),
     O             belli_data,
     I             myThid )

              CALL BELLI_COST_READ_OBS(
     I             num_file, -1, B_Obs_ind_glob(num_file,iObs),
     O             belli_uncert,
     I             myThid )
#endif /* TEST_BELLI_COST_INLOOP */

              IF (belli_data.EQ.-9999) belli_uncert = zeroRL

C Calc misfit or QoI
              IF (belli_uncert.GT.zeroRL) THEN
                belli_weight = 1. _d 0 / (belli_uncert*belli_uncert)

                objf_belli(num_file) = objf_belli(num_file) +
     &                                belli_weight
     &                                *(belli_modval-belli_data)
#ifndef TEST_BELLI_COST_INLOOP
     &                                *(belli_modval-belli_data)
#endif /* TEST_BELLI_COST_INLOOP */

      write(unit=standardmessageunit,fmt='(A,1x,E30.17,1x,A,I0)') 
     &'ESCOBAR J=|m-d|^2(t)= ', objf_belli(num_file), 
     & 'Obs#', num_belli(num_file)
C write(unit=standardmessageunit,fmt='(a,1x,E30.17,1x,a,i0)') 'ESCOBAR dJ_TL(t)=', objf_belli_tl(num_file), 'Obs#', num_belli(num_file)


                num_belli(num_file) = num_belli(num_file) + 1

              ENDIF

            ENDIF !IF (iObs.LE.ObsNo(num_file))
          ENDDO !DO iObs

        ENDIF !IF (myProcId.EQ.0)

      ENDDO !DO num_file

      _END_MASTER( myThid )
#  endif /* TEST_BELLI_COST */
# endif /* ALLOW_COST */

# ifdef ALLOW_COST
C Print cost function values
      DO num_file = 1, NFILESMAX_BELLI
        objf_belli_glo = objf_belli(num_file)

        num_belli_glo = num_belli(num_file)

        WRITE(msgBuf,'(A,I2,A,2D12.5)')
     &    ' belli_cost(', num_file, ' )= ',
     &    objf_belli_glo, num_belli_glo

        IF (num_belli_glo.GT.zeroRL) 
     &    CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

      ENDDO !DO num_file
# endif /* ALLOW_COST */

      WRITE(msgBuf,'(A)') '== belli_cost: end   =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
#endif /* ALLOW_BELLI */

      RETURN
      END
