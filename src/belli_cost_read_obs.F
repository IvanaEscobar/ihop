#include "BELLI_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_READ_OBS

C !INTERFACE:
      SUBROUTINE IHOP_COST_READ_OBS(
     I                      fNb,
     I                      vNb,
     I                      irec,
     O                      vec_loc,
     I                      myThid )

C     !DESCRIPTION:
C     Reads an observation and its uncertainty from a NetCDF ihop
C     input file

C     !USES:
      IMPLICIT NONE
C == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_BELLI
# include "netcdf.inc"
# include "BELLI_SIZE.h"
# include "BELLI.h"
# ifdef ALLOW_COST
#  include "BELLI_COST.h"
# endif
#endif

C     !INPUT PARAMETERS:
C     fNb, vNb  ::
C     irec      ::
C     myThid    :: my thread ID
      INTEGER fNb
      INTEGER vNb
      INTEGER irec
      INTEGER myThid
C     !OUTPUT PARAMETERS:
C     vec_loc   ::
      _RL     vec_loc

C     !LOCAL VARIABLES:
#if (defined ALLOW_BELLI) && (defined ALLOW_COST)
C tmpbObsNo :: number of obs in the file
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER v1, err, varId1, tmpbObsNo
      INTEGER vec_start, vec_count
      _RL     vec_tmp1(1000), vec_tmp2(1000)
CEOP

      IF ( (irec.LT.belli_minind_buff) .OR.
     &     (irec.GT.belli_maxind_buff) .OR.
     &     (belli_curfile_buff.NE.fNb) ) THEN
        err = NF_INQ_DIMID( ncidData(fNb), 'iOBS', varId1 )
        CALL IHOP_COST_NF_ERROR('read_obs 1',err,1,1,myThid)
        err = NF_INQ_DIMLEN( ncidData(fNb), varId1, tmpbObsNo )
        CALL IHOP_COST_NF_ERROR('read_obs 2',err,1,1,myThid)

        IF (belli_curfile_buff.NE.fNb) THEN
C No asumption on whether a forward or a backward loop is calling
          belli_minind_buff = MAX(1,irec-500+1)
          belli_maxind_buff = MIN(tmpbObsNo,irec+500)

        ELSEIF (irec.LT.belli_minind_buff) THEN
C Implies that a backward loop is calling
          belli_minind_buff = MAX(1,irec-999)
          belli_maxind_buff = irec

        ELSE
C Implies that a forward loop is calling
          belli_minind_buff = irec
          belli_maxind_buff = MIN(tmpbObsNo,irec+999)

       ENDIF !IF (belli_curfile_buff.NE.fNb)

#ifdef ALLOW_DEBUG
        IF ( debugLevel.GE.debLevD ) THEN
          WRITE(msgBuf,'(A,5I9)')
     &      'buffer readobsfile ',
     &      belli_minind_buff, belli_maxind_buff,
     &      irec, bObsNo(fNb), tmpbObsNo
          CALL PRINT_MESSAGE( msgBuf,
     &      standardmessageunit, SQUEEZE_RIGHT, myThid )
        ENDIF
#endif

        vec_start = belli_minind_buff
        vec_count = belli_maxind_buff-belli_minind_buff+1

        err = NF_INQ_VARID( ncidData(fNb), belli_nameval,
     &        varId1 )
        CALL IHOP_COST_NF_ERROR('read_obs 3',err,1,1,myThid)
        err = NF_GET_VARA_DOUBLE(ncidData(fNb), varId1,
     &        vec_start, vec_count, vec_tmp1 )
        CALL IHOP_COST_NF_ERROR('read_obs 4',err,1,1,myThid)

        err = NF_INQ_VARID( ncidData(fNb), belli_nameuncert,
     &        varId1 )
        CALL IHOP_COST_NF_ERROR('read_obs 5',err,1,1,myThid)
        err = NF_GET_VARA_DOUBLE(ncidData(fNb), varId1,
     &        vec_start, vec_count, vec_tmp2 )
        CALL IHOP_COST_NF_ERROR('read_obs 6',err,1,1,myThid)

        IF (err.NE.NF_NOERR) THEN
          WRITE(errorMessageUnit,'(A)')
     &      'WARNING in ihop_read_obs: record not found!!'
        ENDIF

        DO v1 = 1, vec_count
          ihop_data_buff(v1) = vec_tmp1(v1)
          ihop_uncert_buff(v1) = vec_tmp2(v1)
        ENDDO

        belli_curfile_buff = fNb
      ENDIF !IF ( (irec.LT.belli_minind_buff)

C Get vec_loc from the buffer
      IF (vNb.LT.0) THEN
        vec_loc = ihop_uncert_buff(irec-belli_minind_buff+1)

      ELSE
        vec_loc = ihop_data_buff(irec-belli_minind_buff+1)

      ENDIF

#endif /* ALLOW_BELLI */

      RETURN
      END
