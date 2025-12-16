#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST

C !INTERFACE:
      SUBROUTINE IHOP_COST( myTime, myIter, myThid )

C    !DESCRIPTION:
C    Computes iHOP data misfit for cost function

C    !USES:
#if (defined ALLOW_IHOP) && (defined TEST_IHOP_COST)
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
#ifdef ALLOW_IHOP
# include "IHOP_SIZE.h"
# include "IHOP.h"
# include "netcdf.inc"
# ifdef ALLOW_COST
#  include "IHOP_COST.h"
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

#ifdef ALLOW_IHOP
C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
C     bi,bj :: Tile indices
      INTEGER bi,bj
      INTEGER num_file, iObs
      _RL ihop_modval, ihop_modvaltmp
      _RL ihop_data, ihop_uncert, ihop_weight
# ifndef ALLOW_CTRL
      INTEGER optimcycle
# endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_FNAM) IHOPObs_File
      INTEGER IL, err
      INTEGER irec
      LOGICAL exst
      _RL  objf_ihop_glo
      _RL  num_ihop_glo
      _RL  tmpgs
# ifdef ALLOW_COST
      _RL  ihopObs_buff(NOBSMAX_IHOP)
      _RL  ihopObs_modval_glob(NOBSMAX_IHOP)
# endif
# ifdef TEST_IHOP_COST
      INTEGER ig, jg
      INTEGER itlo,ithi
      INTEGER jtlo,jthi
      INTEGER i,j,k
      _RL objf_ihop_tloc(nSx,nSy)
      _RL objf_ihop_gloc
# endif

# ifndef ALLOW_CTRL
      optimcycle = 0
# endif

      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') '== ihop_cost: begin =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

C Initialise local storage
#ifdef ALLOW_AUTODIFF_TAMC
CADJ INIT tapelev_ihop_cost = COMMON, NFILESMAX_IHOP
#endif

# ifdef ALLOW_COST
#  ifdef TEST_IHOP_COST
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      DO num_file=1,NFILESMAX_IHOP
C-- Reset local variables
       DO bj = jtlo,jthi
       DO bi = itlo,ithi
         objf_ihop_tloc(bi,bj) = 0.

       ENDDO !DO bi
       ENDDO !DO bj

       objf_ihop_gloc = 0.

       DO bj = jtlo,jthi
       DO bi = itlo,ithi
CC-- OPTION 1: tile-wise summation: F77 Sound Speed
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
C           objf_ihop_tloc(bi,bj) = objf_ihop_tloc(bi,bj)
C     &      + ihop_ssp(i,j,k,bi,bj)*ihop_ssp(i,j,k,bi,bj)
C         ENDIF
C
C        ENDDO !DO i
C        ENDDO !DO j
C        ENDDO !DO k
CC-- tile-wise summation: F77 Sound Speed END

        WRITE(msgBuf,'(A,E22.17)') 'QoI sum(C2) tiled: ', 
     &                      objf_ihop_tloc(bi,bj)
        CALL PRINT_MESSAGE( msgBuf,
     &      standardMessageUnit, SQUEEZE_RIGHT, myThid )

       ENDDO !DO bi
       ENDDO !DO bj

C-- OPTION 2: global summation: F90 Sound Speed
       DO j = 1,Grid%nR
       DO k = 1,Grid%nZ
C- sum of the squares
         objf_ihop_gloc = objf_ihop_gloc 
     &    + SSP%cMat(k,j)*SSP%cMat(k,j)

       ENDDO !DO j
       ENDDO !DO k
C-- global summation: F90 Sound Speed END

Cc--   Do global summation, if using tile-wise summation
C       IF ( objf_ihop_tloc(1,1).NE.0 ) THEN
C         CALL GLOBAL_SUM_TILE_RL( objf_ihop_tloc,objf_ihop_gloc,myThid )
C       ENDIF
C       CALL GLOBAL_SUM_TILE_RL( objf_ihop_tloc,objf_ihop_gloc,myThid )

       _BEGIN_MASTER( myThid )
       objf_ihop(num_file) = objf_ihop(num_file) + objf_ihop_gloc
       _END_MASTER( myThid )

       WRITE(msgBuf,'(A,E22.17)') 'QoI sum(C2) global: ', 
     &                      objf_ihop(num_file)
       CALL PRINT_MESSAGE( msgBuf,
     &      standardMessageUnit, SQUEEZE_RIGHT, myThid )

      ENDDO !DO num_file

#  else /* TEST_IHOP_COST */
      _BEGIN_MASTER( myThid )

      DO num_file = 1, NFILESMAX_IHOP

C File maintenance
        DO bj = 1, nSy
          DO bi = 1, nSx
            IF ( (ObsNo_tiled(num_file,bi,bj).GT.0) .AND.
     &           (ihopDoNcOutput) ) THEN
C Sync tiled file so data is not lost in buffer
              err = NF_SYNC( ncidFWD(num_file,bi,bj) )
              CALL IHOP_COST_NF_ERROR('ihop_cost 1',err,bi,bj,myThid)

            ENDIF

          ENDDO !DO bi
        ENDDO !DO bj

        DO iObs = 1, NOBSMAX_IHOP
          ihopObs_buff(iObs) = zeroRL
        ENDDO

        DO bj = 1, nSy
          DO bi = 1, nSx
C Open tiled files and read to buffer
            DO iObs = 1, NOBSMAX_IHOP
              IF (iObs.LE.ObsNo_tiled(num_file,bi,bj)) THEN
                ihop_modvaltmp = zeroRL

                CALL ACTIVE_READ_IHOP_TILE( 
     I              num_file,
     O              ihop_modvaltmp, 
     I              iObs,.false., optimcycle, bi, bj, myThid,
     I              ihop_dummy(num_file,bi,bj) )

C Save tiled model equi from buffer, see PROFILES_READVECTOR
C Combine all threads to global
                irec = ihopObs_ind_glob_tiled(num_file,iObs,bi,bj)
C IESCO25: dummy +/- of ihop_ssp to get TAF TLM to work
                ihopObs_buff(irec)=ihopObs_buff(irec)+ihop_modvaltmp
     &                             + ihop_ssp(1,1,1,bi,bj)
                ihopObs_buff(irec)=ihopObs_buff(irec)
     &                             - ihop_ssp(1,1,1,bi,bj)
              ENDIF !IF (iObs.LE.ObsNo_tiled(num_file,bi,bj))
            ENDDO !DO iObs

          ENDDO !DO bj
        ENDDO !DO bi

C Combine MPI processes
        DO iObs = 1, NOBSMAX_IHOP
          tmpgs = ihopObs_buff(iObs)
          _GLOBAL_SUM_RL(tmpgs, myThid)
          ihopObs_modval_glob(iObs) = tmpgs

        ENDDO !DO iObs

#ifdef ALLOW_AUTODIFF_TAMC
C     Not really necessary, but cheap (for small NFILESMAX_IHOP)
C     and avoids some hidden recomputations.
CADJ STORE ihopObs_modval_glob
CADJ &     = tapelev_ihop_cost, key = num_file, kind = isbyte
#endif

        IF (myProcId.EQ.0) THEN
C Loop over global obs
          DO iObs = 1, NOBSMAX_IHOP
CC TAF active var init
C            ihop_globaldummy(num_file) = 0
C
            IF (iObs.LE.ObsNo(num_file)) THEN
              ihop_modval = zeroRL
C ESCO Future dev: add time and space averaging here ... see pkg/obsfit

C Extract model equi
              ihop_modval = ihop_modval+ihopObs_modval_glob(iObs)

C Write to global netcdf file
              CALL ACTIVE_WRITE_IHOP_GLOB( 
     I             num_file, ihop_modval,
     I             iObs, optimcycle, myThid,
     I             ihop_globaldummy(num_file) )

            ENDIF !IF (iObs.LE.ObsNo(num_file))

          ENDDO !DO iObs

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE ihop_globaldummy(num_file)
CADJ &     = tapelev_ihop_cost, key = num_file, kind = isbyte
#endif

C Sync global obs
          err = NF_SYNC( ncidGLOB(num_file) )
          CALL IHOP_COST_NF_ERROR(
     &         'IHOP_COST: NF_SYNC ncidGLOB',err,0,0,myThid )

          DO iObs = 1, NOBSMAX_IHOP
            IF (iObs.LE.ObsNo(num_file)) THEN
              ihop_modval = zeroRL
              ihop_data   = zeroRL
              ihop_uncert = zeroRL

C Read model equivalent from global file
              CALL ACTIVE_READ_IHOP_GLOB( 
     I             num_file,
     O             ihop_modval, 
     I             iObs,.false., optimcycle, myThid,
     I             ihop_globaldummy(num_file) )

#ifdef TEST_IHOP_COST_INLOOP
C Set uncertainty to one
              ihop_uncert = 1. _d 0

#else /* TEST_IHOP_COST_INLOOP */
C Read observation and uncertainty
              CALL IHOP_COST_READ_OBS( 
     I             num_file,  1, ihopObs_ind_glob(num_file,iObs),
     O             ihop_data, 
     I             myThid )

              CALL IHOP_COST_READ_OBS(
     I             num_file, -1, ihopObs_ind_glob(num_file,iObs),
     O             ihop_uncert,
     I             myThid )
#endif /* TEST_IHOP_COST_INLOOP */

              IF (ihop_data.EQ.-9999) ihop_uncert = zeroRL

C Calc misfit
              IF (ihop_uncert.GT.zeroRL) THEN
                ihop_weight = 1. _d 0 / (ihop_uncert*ihop_uncert)

                objf_ihop(num_file) = objf_ihop(num_file) + 
     &                                ihop_weight
     &                                *(ihop_modval-ihop_data)
     &                                *(ihop_modval-ihop_data)

                num_ihop(num_file) = num_ihop(num_file) + 1

              ENDIF

            ENDIF !IF (iObs.LE.ObsNo(num_file))
          ENDDO !DO iObs

        ENDIF !IF (myProcId.EQ.0)

      ENDDO !DO num_file

      _END_MASTER( myThid )
#  endif /* TEST_IHOP_COST */
# endif /* ALLOW_COST */

# ifdef ALLOW_COST
C Print cost function values
      DO num_file = 1, NFILESMAX_IHOP
        objf_ihop_glo = objf_ihop(num_file)

        num_ihop_glo = num_ihop(num_file)

        WRITE(msgBuf,'(A,I2,A,2D12.5)')
     &    ' ihop_cost(', num_file, ' )= ',
     &    objf_ihop_glo, num_ihop_glo

        IF (num_ihop_glo.GT.zeroRL) 
     &    CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

      ENDDO !DO num_file
# endif /* ALLOW_COST */

      WRITE(msgBuf,'(A)') '== ihop_cost: end   =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
#endif /* ALLOW_IHOP */

      RETURN
      END
