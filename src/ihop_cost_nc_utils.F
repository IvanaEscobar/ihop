#include "IHOP_OPTIONS.h"
#include "AD_CONFIG.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_NF_ERROR

C !INTERFACE:
      SUBROUTINE IHOP_COST_NF_ERROR( message, STATUS, bi, bj, myThid )

C     !DESCRIPTION:
C     Print NetCDF error message

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#ifdef ALLOW_IHOP
# include "SIZE.h"
# include "EEPARAMS.h"
# include "PARAMS.h"
# include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     message   :: optional message
C     STATUS    :: NetCDF error status
C     bi,bj     :: Tile indices
C     myThid    :: my thread ID
      CHARACTER*(*) message
      INTEGER STATUS, bi, bj, myThid

C     !OUTPUT PARAMETERS:

C     !FUNCTIONS:
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
      INTEGER IL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      IF (debugLevel .GE. debLevA .AND. STATUS .NE. NF_NOERR) THEN
        IL = ILNBLNK(message)
        IF ( IL .GT. 0 ) THEN
          WRITE(msgBuf,'(2A,2I3,1X,A)')
     &      'NF: IHOP_COST_', message(1:IL),
     &      bi, bj, NF_STRERROR(STATUS)
        ELSE
          WRITE(msgBuf,'(A,2I3,1X,A)') 'NF: IHOP_COST_',
     &      bi, bj, NF_STRERROR(STATUS)
        ENDIF
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )
        CALL PRINT_ERROR( msgBuf, myThid )
      ENDIF
! Print if there are no errors
!      IF (STATUS.EQ.NF_NOERR .AND. ILNBLNK(message).GT.0 ) THEN
!        IL = ILNBLNK(message)
!        WRITE(msgBuf,'(2A)') "NF IHOP NOERR: ", message(1:IL)
!        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
!     &       SQUEEZE_RIGHT, myThid )
!      ENDIF
#endif /* ALLOW_IHOP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: IHOP_COST_NC_CLOSE

C !INTERFACE:
      SUBROUTINE IHOP_COST_NC_CLOSE( myThid )

C     !DESCRIPTION:
C     Close NetCDF files

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#ifdef ALLOW_IHOP
# include "SIZE.h"
# include "EEPARAMS.h"
# include "IHOP_SIZE.h"
# include "IHOP_COST.h"
# include "netcdf.inc"
#endif

C     !INPUT PARAMETERS:
C     myThid :: my Thread Id number
      INTEGER myThid

C     !OUTPUT PARAMETERS:

C     !FUNCTIONS:
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
#ifdef ALLOW_IHOP
C     bi,bj :: Tile indices
      INTEGER bi, bj
      INTEGER IL, numFile
      INTEGER STATUS
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      IF (ihopDoNcOutput) THEN
! Close tiled files
        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            DO numFile = 1, NFILESMAX_IHOP
              IL = ILNBLNK( ihopObs_Files(numFile) )

              IF (IL.NE.0) THEN
                WRITE(msgBuf,'(A,A,5(1X,I8))')
     &            'S/R IHOP_COST_NC_CLOSE:',
     &            ' Closing '//ihopObs_Files(numFile)(1:IL), numFile,
     &           ncidFWD(numFile,bi,bj),
     &           ncidAD(numFile,bi,bj), ncidTL(numFile,bi,bj)
                CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &               SQUEEZE_RIGHT, myThid )

C Forward equivalent (.equi) file
                STATUS = NF_CLOSE( ncidFWD(numFile,bi,bj) )
                WRITE(msgBuf,'(A,I3,I8)')
     &            'NC_CLOSE: NF_CLOSE fwd unit',
     &            numFile, ncidFWD(numFile,bi,bj)
                CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &               bi, bj, myThid )
# ifdef ALLOW_ADJOINT_RUN
C Adjoint equivalent (.equi) file
                STATUS = NF_CLOSE( ncidAD(numFile,bi,bj) )
                WRITE(msgBuf,'(A,I3,I8)')
     &            'NC_CLOSE: NF_CLOSE adj unit',
     &            numFile, ncidAD(numFile,bi,bj)
                CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &               bi, bj, myThid )
# endif
# ifdef ALLOW_TANGENTLINEAR_RUN
C Tangent linear equivalent (.equi) file
                STATUS = NF_CLOSE( ncidTL(numFile,bi,bj) )
                WRITE(msgBuf,'(A,I3,I8)')
     &            'NC_CLOSE: NF_CLOSE tlm unit',
     &            numFile, ncidTL(numFile,bi,bj)
                CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &               bi, bj, myThid )
# endif
              ENDIF !(IL.NE.0)

            ENDDO !DO numFile
          ENDDO !DO bi 
        ENDDO !DO bj

! Close global files
        DO numFile = 1, NFILESMAX_IHOP
          IL = ILNBLNK( ihopObs_Files(numFile) )

          IF (IL.NE.0) THEN
            WRITE(msgBuf,'(A,A,5(1X,I8))')
     &        'S/R IHOP_COST_NC_CLOSE:',
     &        ' Closing '//ihopObs_Files(numFile)(1:IL), numFile,
     &       ncidData(numFile), ncidGLOB(numFile),
     &       ncidADGLOB(numFile), ncidTLGLOB(numFile)
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &           SQUEEZE_RIGHT, myThid )
C Data file
            STATUS = NF_CLOSE( ncidData(numFile) )
            WRITE(msgBuf,'(A,I3,I8)')
     &        'NC_CLOSE: NF_CLOSE data unit',
     &        numFile, ncidData(numFile)
            CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &           1, 1, myThid )

C Forwaruivalent (.equi) file
            STATUS = NF_CLOSE( ncidGLOB(numFile) )
            WRITE(msgBuf,'(A,I3,I8)')
     &        'NC_CLOSE: NF_CLOSE fwd unit',
     &        numFile, ncidGLOB(numFile)
            CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &           1, 1, myThid )
# ifdef W_ADJOINT_RUN
C Adjoint equivalent (.equi) file
            STATUS = NF_CLOSE( ncidADGLOB(numFile) )
            WRITE(msgBuf,'(A,I3,I8)')
     &        'NC_CLOSE: NF_CLOSE adj unit',
     &        numFile, ncidADGLOB(numFile)
            CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &           1, 1, myThid )
# endif
# ifdef W_TANGENTLINEAR_RUN
C Tangent linear equivalent (.equi) file
            STATUS = NF_CLOSE( ncidTLGLOB(numFile) )
            WRITE(msgBuf,'(A,I3,I8)')
     &        'NC_CLOSE: NF_CLOSE tlm unit',
     &        numFile, ncidTLGLOB(numFile)
            CALL IHOP_COST_NF_ERROR( msgBuf, STATUS,
     &           1, 1, myThid )
# endif
          ENDIF !(IL.NE.0)

        ENDDO !DO numFile

      ENDIF !IF (ihopDoNcOutput)
#endif /* ALLOW_IHOP */

      RETURN
      END