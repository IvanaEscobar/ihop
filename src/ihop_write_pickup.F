#include "IHOP_OPTIONS.h"

CBOP
C !ROUTINE: IHOP_WRITE_PICKUP

C !INTERFACE:
      SUBROUTINE IHOP_WRITE_PICKUP( permPickup,
     &                    suff, myTime, myIter, myThid )

C !DESCRIPTION:
C     Writes current state of passive tracers to a pickup file

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "IHOP_SIZE.h"
#include "IHOP.h"

C !INPUT PARAMETERS:
C     permPickup      :: write a permanent pickup
C     suff            :: suffix for pickup file (eg. ckptA)
C     myTime          :: model time
C     myIter          :: time-step number
C     myThid          :: thread number
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !OUTPUT PARAMETERS:
C  none

#ifdef ALLOW_IHOP
#if (defined IHOP_3D_STATE) || (defined IHOP_2D_STATE)

C     === Functions ====
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES:
C     j           :: loop index / field number
C     nj          :: record number
C     fp          :: pickup-file precision
C     glf         :: local flag for "globalFiles"
C     fn          :: character buffer for creating filename
C     nWrFlds     :: number of fields being written
C     listDim     :: dimension of "wrFldList" local array
C     wrFldList   :: list of written fields
C     msgBuf      :: Informational/error message buffer
      INTEGER j, nj, fp, lChar
      LOGICAL glf
      _RL     timList(1)
      CHARACTER*(MAX_LEN_FNAM) fn
      INTEGER listDim, nWrFlds
      PARAMETER( listDim = 12 )
      CHARACTER*(8) wrFldList(listDim)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C        lChar = ILNBLNK(suff)
C        IF ( lChar.EQ.0 ) THEN
C          WRITE(fn,'(2A)') 'pickup_ihop'
C        ELSE
C          WRITE(fn,'(2A)') 'pickup_ihop.',suff(1:lChar)
C        ENDIF
C        fp = precFloat64
C        j  = 0
C
CC       Firstly, write 3-D fields as consecutive records,
C
C#ifdef IHOP_3D_STATE
CC     record number < 0 : a hack not to write meta files now:
C        j = j + 1
C        CALL WRITE_REC_3D_RL( fn, fp, Nr,
C     &                        ihop_ssp, -j, myIter, myThid )
C        IF (j.LE.listDim) wrFldList(j) = 'ihop_ssp'
C
C#endif /* IHOP_3D_STATE */
C
CC-    switch to 2-D fields:
C        nj = -j*Nr
C
C#ifdef IHOP_2D_STATE
C#endif /* IHOP_2D_STATE */
C
CC--------------------------
C        nWrFlds = j
C        IF ( nWrFlds.GT.listDim ) THEN
C          WRITE(msgBuf,'(2A,I5,A)') 'IHOP_WRITE_PICKUP: ',
C     &     'trying to write ',nWrFlds,' fields'
C          CALL PRINT_ERROR( msgBuf, myThid )
C          WRITE(msgBuf,'(2A,I5,A)') 'IHOP_WRITE_PICKUP: ',
C     &     'field-list dimension (listDim=',listDim,') too small'
C          CALL PRINT_ERROR( msgBuf, myThid )
C          CALL ALL_PROC_DIE( myThid )
C          STOP 'ABNORMAL END: S/R IHOP_WRITE_PICKUP (list-size Pb)'
C        ENDIF
C#ifdef ALLOW_MDSIO
CC Uses this specific S/R to write only meta files
C        j  = 1
C        nj = ABS(nj)
C        IF ( nWrFlds*Nr .EQ. nj ) THEN
C          j  = Nr
C          nj = nWrFlds
C        ENDIF
C        glf  = globalFiles
C        timList(1) = myTime
C        CALL MDS_WR_METAFILES( fn, fp, glf, .FALSE.,
C     &                         0, 0, j, ' ',
C     &                         nWrFlds, wrFldList,
C     &                         1, timList, oneRL,
C     &                         nj, myIter, myThid )
C#endif /* ALLOW_MDSIO */
CC--------------------------
C
C#endif /* IHOP_3D_STATE or IHOP_2D_STATE */
#endif /* ALLOW_IHOP */

      RETURN
      END
