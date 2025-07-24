#include "IHOP_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: COST_IHOP

C !INTERFACE:
      SUBROUTINE COST_IHOP( myIter, myTime, myThid )

C    !DESCRIPTION:
C    Computes iHOP data misfit for cost function

C    !USES:
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
      INTEGER myIter, myThid
C     !OUTPUT PARAMETERS:
CEOP

C     !FUNCTIONS:
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
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

# ifndef ALLOW_CTRL
      optimcycle = 0
# endif

      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') '== cost_ihop: begin =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

        _BEGIN_MASTER( myThid )

# ifdef ALLOW_COST
      DO num_file = 1, NFILESMAX_IHOP
        DO bj = 1, nSy
          DO bi = 1, nSx
C File maintenance
            IF ( (ObsNo_tiled(num_file,bi,bj).GT.0) .AND.
     &           (ihopDoNcOutput) ) THEN
C Sync tiled file so data is not lost in buffer
              err = NF_SYNC( ncidFWD(num_file,bi,bj) )
              CALL IHOP_COST_NF_ERROR('cost_ihop 1',err,bi,bj,myThid)

            ENDIF

          ENDDO !DO bi
        ENDDO !DO bj

        DO iObs = 1, NOBSMAX_IHOP
          ihopObs_buff(iObs) = 0.0
        ENDDO

        DO bj = 1, nSy
          DO bi = 1, nSx
C Open tiled files and read to buffer
            DO iObs = 1, NOBSMAX_IHOP
              IF (iObs.LE.ObsNo_tiled(num_file,bi,bj)) THEN
                ihop_modvaltmp = 0.

                CALL ACTIVE_READ_IHOP_TILE( num_file,
     &              ihop_modvaltmp, iObs,
     &              .false., optimcycle, bi, bj, myThid,
     &              ihop_dummy(num_file,bi,bj) )

C Save tiled model equi from buffer, see PROFILES_READVECTOR
C Combine all threads to global
                irec = ihopObs_ind_glob_tiled(num_file,iObs,bi,bj)
                ihopObs_buff(irec) =ihopObs_buff(irec)+ihop_modvaltmp

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

        IF (myProcId.EQ.0) THEN
C Loop over global obs
          DO iObs = 1, NOBSMAX_IHOP
C TAF active var init
            ihop_globaldummy(num_file) = 0

            IF (iObs.LE.ObsNo(num_file)) THEN
              ihop_modval = 0.
C ESCO Future dev: add time and space averaging here ...

C Extract model equi
              ihop_modval = ihop_modval+ihopObs_modval_glob(iObs)

C Write to global netcdf file
              CALL ACTIVE_WRITE_IHOP_GLOB( num_file,
     &             ihop_modval,
     &             iObs, optimcycle, myThid,
     &             ihop_globaldummy(num_file) )

            ENDIF !IF (iObs.LE.ObsNo(num_file))

          ENDDO !DO iObs

C Sync global obs
          err = NF_SYNC( ncidGLOB(num_file) )
          CALL IHOP_COST_NF_ERROR('cost_ihop 2',err,bi,bj,myThid)

          DO iObs = 1, NOBSMAX_IHOP
            IF (iObs.LE.ObsNo(num_file)) THEN
#  ifdef ALLOW_AUTODIFF_TAMC
!$TAF STORE ihop_globaldummy = comlev1, key=ikey_dynamics, byte=isbyte
#  endif
              ihop_modval = 0.
              ihop_data = 0.
              ihop_uncert = 0.

C Read model equivalent from global file
              CALL ACTIVE_READ_IHOP_GLOB( num_file,
     &             ihop_modval, iObs,
     &             .false., optimcycle, myThid,
     &             ihop_globaldummy(num_file) )

C Read observation and uncertainty
              CALL IHOP_COST_READ_OBS( num_file,  1,
     &             ihopObs_ind_glob(num_file,iObs),
     &             ihop_data, myThid )

              CALL IHOP_COST_READ_OBS( num_file, -1,
     &             ihopObs_ind_glob(num_file,iObs),
     &             ihop_uncert, myThid )

              IF (ihop_data.EQ.-9999) ihop_uncert = 0.

              IF (ihop_uncert.GT.0.) THEN
                ihop_weight = 1. / (ihop_uncert*ihop_uncert)

                objf_ihop(num_file) = objf_ihop(num_file)
     &           +ihop_weight
     &           *(ihop_modval-ihop_data)
     &           *(ihop_modval-ihop_data)

                num_ihop(num_file) = num_ihop(num_file)+1

              ENDIF

            ENDIF !IF (iObs.LE.ObsNo(num_file))
          ENDDO !DO iObs

        ENDIF !IF (myProcId.EQ.0)

      ENDDO !DO num_file

# endif /* ALLOW_COST */

      _END_MASTER( myThid )

# ifdef ALLOW_COST
C Print cost function values
      DO num_file = 1, NFILESMAX_IHOP
        objf_ihop_glo = objf_ihop(num_file)

        num_ihop_glo = num_ihop(num_file)

        WRITE(msgBuf,'(A,I2,A,2D12.5)')
     &    ' cost_ihop(', num_file, ' )= ',
     &    objf_ihop_glo, num_ihop_glo

        IF (num_ihop_glo.GT.0.) 
     &    CALL PRINT_MESSAGE( msgBuf,
     &         standardMessageUnit, SQUEEZE_RIGHT, myThid )

      ENDDO !DO num_file
# endif /* ALLOW_COST */

      WRITE(msgBuf,'(A)') '== cost_ihop: end   =='
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf,
     &     standardMessageUnit, SQUEEZE_RIGHT, myThid )

#endif /* ALLOW_IHOP */

      RETURN
      END